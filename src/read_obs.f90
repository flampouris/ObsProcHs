!> @file read_obs.f90
!! @brief The read_obs module imports the observations (InSitu and from Satellite).
!! It includes the following subroutines functions:
!! read_obs       : Main \n
!! read_prepbufr  : Imports insitu obs from files with prepbufr format \n
!! read_satmar    : Imports satellite obs from files with prepbufr format \n
!! lldistm        : Calculates Pythagorean and Haversine distance (m) from LatLon \n
!! datesec        : Calculates time in minutes from 1.1.1978 (the seconds are decimals) \n
!! counter        : Counting number of messages in a prepbufr file (total number of Obs) \n
!! file_existence : Check if file exist \n
!! count_kwrd     : Counts words in a string \n
!!
!! @author stelios flampouris
!! @version 0.3
!> @copyright GNU GPL
!<
!! @date 05-Mar-2016 - Creation
!! @date 26-Oct-2016 - Last Update
!! @date 18-Nov-2016 by stelios the satmar output is sorted in time
!<

module read_obs_mod
   use variables_obsop_hs, only : sat_obs, ins_obs, dodebug1,dofilterS
   use kinds, only : i_kind,r_kind
   use constants_hs, only : i0, i1, i2, i25
   use general_func, only : writeMatrixdp
   use sorting
   use despiking
!
   implicit none
   private
   public read_obs
   contains
!> 
!! @brief 
!! !>
!! @brief The read_obs reads all the available observations
!! param[in] obsSatinfile     - filename with the satelite data
!! param[in] obsInSinfile     - filename with the insitu data
!! param[in] gstime           - time of prediction, 
!! param[in] winlen           - temporal window for keeping data
!! param[in] mbr_domain       - array with the coordinates of the domain.
!! param[out] ins_obs         - see the variables_hs module for details 
!<
!> @details
!! 1. This routine reads all observations end sends them to the main 
!<
   subroutine read_obs (nsatflnm, obsSatinfile, obsInSinfile, gstime         &
                      , winlen, mbr_domain, do_superobs, do_despiking        &
                      , despikeStdFac,do_ma,MaLngth, SprObsDim, n_obs        )
      use variables_obsop_hs, only: allocate_all, all_obs, destroy_obs, sat_obs_dum
!
      integer(i_kind)            , intent(in)  :: nsatflnm
      character(len=*), dimension(nsatflnm), intent(in    ) :: obsSatinfile
      character(len=*)           , intent(in ) :: obsInSinfile
      real(r_kind)               , intent(in ) :: gstime
      real(r_kind)               , intent(in ) :: winlen
      real(r_kind)               , intent(in ) :: despikeStdFac
      real(r_kind), dimension(*) , intent(in ) :: mbr_domain
      logical                    , intent(in ) :: do_superobs
      logical                    , intent(in ) :: do_despiking
      logical                    , intent(in ) :: do_ma      
      real(r_kind), dimension(*) , intent(in ) :: SprObsDim
      integer(i_kind)            , intent(in ) :: MaLngth
      integer(i_kind)            , intent(out) :: n_obs
!
      character(i25), parameter                 :: myname='read_obs'
      integer(i_kind)                           :: n_sat=i0,n_ins=i0
      integer(i_kind)                           :: k
      logical                                   :: file_exists
!
      if (.not.allocated(sat_obs_dum)) allocate (sat_obs_dum(500000_i_kind))
      do k=i1,nsatflnm,i1 
         print*, trim(myname),' Reading Satellite file: ', trim(obsSatinfile(k))
         call read_satmar  (obsSatinfile(k), gstime, winlen, mbr_domain &
                           , do_superobs, do_despiking, despikeStdFac   &
                           , do_ma, MaLngth, SprObsDim)
         sat_obs_dum(n_sat+i1:n_sat+size(sat_obs)) = sat_obs(1:n_sat)
         n_sat = n_sat+size(sat_obs)
         print*,myname, n_sat
      end do 
!
      inquire(file=trim(obsInSinfile), exist=file_exists)
      if (file_exists) then
         print*, trim(myname),' Reading Insitu Observations: ', trim(obsInSinfile)
         call read_prepbufr(obsInSinfile, gstime, winlen, mbr_domain)
         n_ins = size(ins_obs)
      end if 
!
      n_obs = n_sat+n_ins
      
      if (n_obs > i0) then 
         call allocate_all(n_obs)
      else 
         print*, trim(myname), ': There are not Observations. Termination'
         stop
      end if
!
      if (n_sat>0) then
         all_obs(1:n_sat)=sat_obs_dum(1:n_sat)
      else
         print*, trim(myname), ': There are not Satellite Observations'
      end if
      if (n_ins>0) then
         all_obs(n_sat+1:n_sat+n_ins)=ins_obs
      else
         print*, trim(myname), ': There are not insitu Observations'
      end if
!
      call destroy_obs
      if (allocated(sat_obs_dum)) deallocate (sat_obs_dum)
   end subroutine read_obs
!
!----
!> @brief The subroutine creates superobs from satellite observations
!! param[in] sat_obs     - The kept satellite data
!! param[in] SprObsDim   - The delta_x and _y defining the domain for the each
!! param[in] mbr_domain  - array with the coordinates of the domain.
!! mbr_domain = [minval(latg), maxval(latg), minval(long),
!! superob (in degrees), e.g. /0.5, 0.5/ 
!<
!> @details
!! This routine creates superobs in the domains with dimensions defined by SprObsDim
!! The position of each superobs is the mean position of the observations for
!! the specific domain 
!! The value of  superobs is the mean value of the observations in the specific
!! domain
!<
   subroutine super_obs (SprObsDim,mbr_domain, nobs, lonobs, latobs, valobs     &
                        , said, nSobs)
      use kinds, only : i_kind, r_kind
      use constants_hs, only : i1,i2,i3,i4
!
      real(r_kind), dimension(*) , intent(in) :: SprObsDim
      real(r_kind), dimension(*) , intent(in) :: mbr_domain
      integer (i_kind), dimension(:)  , intent(in) :: nobs 
      real(r_kind), dimension(:), intent(inout) :: lonobs, latobs

      real(r_kind), dimension(:,:), intent(inout) :: valobs
      character(i25) , dimension(nobs(i2)), optional   :: said
      integer(i_kind)                  , intent(out)   :: nSobs
! local
      character(i25), parameter                 :: myname='super_obs'

      integer(i_kind)                           :: kend,lend,kk,ll,lll
      real(r_kind), allocatable, dimension(:,:) :: gridlon,gridlat
      logical, dimension(nobs(i2)) :: maskcoor
      real(r_kind), dimension(size(lonobs,1)) :: lonobsdum, latobsdum
      character(i25) , dimension(nobs(i2))    :: saidd

      lonobsdum = 0.
      latobsdum = 0.
      saidd   = ''
!
      kend = ceiling((mbr_domain(i2)-mbr_domain(i1))/SprObsDim(i2))+1 
      lend = ceiling((mbr_domain(i4)-mbr_domain(i3))/SprObsDim(i1))+1
!
      if(.not.allocated(gridlon))allocate(gridlon(lend,kend))
      if(.not.allocated(gridlat))allocate(gridlat(lend,kend))
      do ll = 0,lend-1,i1
         do kk = 0,kend-1,i1
            gridlon(ll+1,kk+1) = mbr_domain(i3)+ll*SprObsDim(i1)
            gridlat(ll+1,kk+1) = mbr_domain(i1)+kk*SprObsDim(i2)
         end do
      end do
!
      nSobs = 1
      do ll = i1,lend-1,i1
         do kk = i1,kend-1,i1
            maskcoor = .false.
            maskcoor=lonobs>gridlon(ll,kk).and.lonobs<gridlon(ll+1,kk+1).and.  &
                     latobs>gridlat(ll,kk).and.latobs<gridlat(ll+1,kk+1)
            if (any(maskcoor)) then
               lonobsdum(nSobs) =  sum(lonobs, maskcoor)/(max(1,count(maskcoor)))
               latobsdum(nSobs) =  sum(latobs, maskcoor)/(max(1,count(maskcoor)))
               if (present(said))then
                  do lll=i1,nobs(i2)
                     if (maskcoor(lll)) then
                        saidd(nSobs) = said(lll) 
                     end if
                  end do
               end if
               do lll = i1,nobs(i1),i1
                  valobs(lll,nSobs) =  sum(valobs(lll,:), maskcoor)/(max(1,count(maskcoor)))
               end do
               nSobs = nSobs+1
            end if
         end do
      end do
! 
      nSobs = nSobs-i1
      lonobs = lonobsdum
      latobs = latobsdum
      said   = saidd
!
      if(allocated(gridlon))deallocate(gridlon)
      if(allocated(gridlat))deallocate(gridlat)
!
   end subroutine super_obs
!
!-----
! 
!> @brief The @a reads insitu observations
!! param[in] infile      - filename with the satelite data
!! param[in] gstime      - time of prediction, 
!! param[in] winlen      - temporal window for keeping data
!! param[in] mbr_domain  - array with the coordinates of the domain.
!! param[out] ins_obs    - see the variables_hs module for details 
!!     
!<
!> @details
!! 1. This routine reads insitu observations 
!!
!! For reading the data of interest, the "headers" (hdr_ variables) have to be
!! modified accordingly; in this case, Significant Wave Height (howv or hs) data are
!! imported.
!!
!! @todo The value of the dflt_error varies upon the platform, gulev_2003
!<
!
   subroutine read_prepbufr(infile, gstime, winlen, mbr_domain)
      use variables_obsop_hs, only : allocate_ins, ins_obs
      use kinds, only: r_kind,i_kind
      use constants_hs, only : uid_rd                                   &
                              ,i0, i1, i6, i10, i256                    &
                              ,r0, r1, r90, r360                        &
                              ,tiny_r_kind                              &
                              ,r10inv, r60inv
!
      character(len=*)           , intent(in    ) :: infile
      real(r_kind)               , intent(in    ) :: gstime
      real(r_kind)               , intent(in    ) :: winlen
      real(r_kind), dimension(*) , intent(in    ) :: mbr_domain
!
! local variables
      character (i25), parameter :: myname='read_prepbufr'
	   real(r_kind), allocatable, dimension(:) :: howv
	   real(r_kind), allocatable, dimension(:) :: hdr, loc_1d, tm_1d
      integer (i_kind), allocatable, dimension(:) :: tm_1dMN
!
      character(i25) :: howvstr='HOWV'
      character(i256):: hdstr='SID XOB YOB DHR TYP ELV SAID T29'
      character(i256), parameter:: hdr_tm       = 'YEAR MNTH DAYS HOUR MINU SECU'
      character(i256), parameter:: hdr_loc      = 'CLAT CLON'
      
      integer(i_kind),parameter :: howvMax = 40_i_kind

      integer(i_kind) :: levs
      real ::a
!
! shared variables among module's subroutines
      logical           :: bool_open
      integer(i_kind)   :: cnt,ll
      integer(i_kind), parameter                   :: nreal = i6  !23
      integer(i_kind), allocatable, dimension(:)   :: isort, iloc
      character(i25) , allocatable, dimension(:)   :: said
!
      integer(i_kind) :: nrec, irec, nmind
      integer(i_kind) :: ireadmg,ireadsb,idate
      real(r_kind)   , allocatable, dimension(:,:) :: data_all
      integer(i_kind) :: n_howv, n_hdr, n_lc, n_tm
      character(len=i25) :: subset
      real(r_kind) :: rminobs
      integer(i_kind) :: ndata, iout
      real(r_kind)   ,parameter :: dflt_err = 0.35_r_kind
      real(r_kind) :: said_1d
!
      n_hdr    = count_kwrd(hdstr)
      n_howv   = count_kwrd(howvstr)
      n_lc     = count_kwrd(hdr_loc)
      n_tm     = count_kwrd(hdr_tm)
!
      if(.not.allocated(howv     )) allocate(howv     (n_howv   ))
      if(.not.allocated(hdr      )) allocate(hdr      (n_hdr    ))
      if(.not.allocated(loc_1d   )) allocate(loc_1d   (n_lc     ))
      if(.not.allocated(tm_1d    )) allocate(tm_1d    (n_tm     ))
      if(.not.allocated(tm_1dMN  )) allocate(tm_1dMN  (n_tm-i1  ))
!
      irec     = i0
      nrec     = i0
      ndata    = i0
      said_1d  = r0

      inquire( unit=uid_rd, opened=bool_open )
      if(bool_open) call closbf(uid_rd)
!
      call file_existence(infile, myname)
!
      call counter(uid_rd, infile, cnt)
!
      if(.not.allocated(data_all ))   allocate(data_all  (nreal, cnt ))
      if(.not.allocated(isort    ))   allocate(isort     (cnt     ))
      if(.not.allocated(said     ))   allocate(said      (cnt     ))
!
      isort = i0
!
      open(uid_rd,file=trim(infile),action='read',form='unformatted') 
!
      call openbf(uid_rd,'IN',uid_rd)
!
      call datelen(i10)
!
      read_msg: do while(ireadmg(uid_rd,subset,idate) == i0)
         read_loop: do while (ireadsb(uid_rd) == i0)
            nrec = nrec + i1
            howv=howvMax             
            call ufbint(uid_rd, howv, n_howv, 1, irec, trim(howvstr))
            if ( howv(n_howv) < howvmax) then   
               irec = irec + 1 
            else 
               cycle read_loop
            end if
!
            if (howv(n_howv)<=r10inv) cycle read_loop
!
            call ufbint(uid_rd, loc_1d, n_lc, 1, irec, trim(hdr_loc))
            if(abs(loc_1d(1))>r90 .or. abs(loc_1d(2))>r360) cycle read_loop
!
            if(loc_1d(2)== r360)loc_1d(2)=loc_1d(2)-r360
 	         if(loc_1d(2) < r0)loc_1d(2)=loc_1d(2)+r360
!
            if (loc_1d(1)<mbr_domain(1).or.loc_1d(1)>mbr_domain(2)) cycle  read_loop!xcld
            if (loc_1d(2)<mbr_domain(3).or.loc_1d(2)>mbr_domain(4)) cycle  read_loop!xcld
!
            call ufbint(uid_rd, tm_1d, n_tm, 1, irec, trim(hdr_tm))
            tm_1d(6)=r0
            tm_1dMN = int(tm_1d(1:5))
            call w3fs21(tm_1dMN, nmind)
            rminobs=real(nmind,r_kind)+(real(tm_1d(6),r_kind)*r60inv)
            if (abs(gstime-rminobs)*r60inv <= winlen) then 
!               space for further processing in temporal space, if required 
            else 
               cycle read_loop                                              !qc
            end if
!
!  Slayer  : Just in Case 
!
            ndata=ndata+i1
            isort(nrec)=ndata
!
            data_all(1,ndata) = loc_1d(2)
            data_all(2,ndata) = loc_1d(1)
            data_all(3,ndata) = howv(1)
            data_all(4,ndata) = dflt_err
            data_all(5,ndata) = rminobs
            data_all(6,ndata) = r1
            data_all(7,ndata) = said_1d
            said(ndata) = subset
         enddo read_loop
      enddo read_msg
!      
      call closbf(uid_rd)
!
      if(.not.allocated(iloc)) allocate(iloc(ndata))

      cnt = i0
      do ll = i1,size(data_all,2)
         if(isort(ll) > i0) then
            cnt = cnt + i1
            iloc(cnt) = isort(ll)
         end if
      end do
!
      if (.not.allocated(ins_obs)) call allocate_ins(cnt)
      do ll = i1,ndata
         iout=iloc(ll)
            ins_obs(ll)%x_grd(1) =  data_all(1,iout)
            ins_obs(ll)%x_grd(2) =  data_all(2,iout)
            ins_obs(ll)%value    =  data_all(3,iout)
            ins_obs(ll)%oerr     =  data_all(4,iout)
            ins_obs(ll)%hour     =  data_all(5,iout)
            ins_obs(ll)%qkey     =  data_all(6,iout)
            ins_obs(ll)%typ      =  said(iout)
            ins_obs(ll)%kept     = .true.
      end do
!
      if(allocated(howv     )) deallocate(howv    )
      if(allocated(hdr      )) deallocate(hdr     )
      if(allocated(loc_1d   )) deallocate(loc_1d  )
      if(allocated(tm_1d    )) deallocate(tm_1d   )
      if(allocated(tm_1dMN  )) deallocate(tm_1dMN )
      if(allocated(data_all )) deallocate(data_all)
      if(allocated(isort    )) deallocate(isort   )
      if(allocated(said     )) deallocate(said    )
      if(allocated(iloc     )) deallocate(iloc    )
!
   end subroutine read_prepbufr         
!
!>
!! @brief The satmar reads altimeter data
!! param[in] infile      - filename with the satelite data
!! param[in] gstime      - time of prediction, 
!! param[in] winlen      - temporal window for keeping data
!! param[in] mbr_domain  - array with the coordinates of the domain.
!! param[out] sat_obs    - see the variables_hs module for details 
!!     
!<
!> @details
!! 1. This routine reads data from five different altimeters (Jason-2,3,
!!    Cryosat-2, Saral/Altika and Sentinel3a
!! 2. It uses the provided flags for QC. 
!! 3. Observations only within the domain of interest are retained.
!!
!! For reading the data of interest, the "headers" (hdr_ variables) have to be
!! modified accordingly; in this case, Significant Wave Height (howv or hs) data are
!! imported. 
!<
subroutine read_satmar ( infile, gstime, winlen, mbr_domain, do_superobs   &
                       , do_despiking, despikeStdFac, do_ma, MaLngth       &
                       , SprObsDim                                         )
      use variables_obsop_hs, only : allocate_alt, sat_obs,dodebug1,dofilters
      use constants_hs, only : uid_rd                                         &
                              ,i0, i1, i2, i4, i3, i5, i6, i8, i10, i256      &
                              ,r0, r1, r60, r90, r360                         &
                              ,tiny_r_kind                                    &
                              ,r60inv, deg2rad                                &
                              ,grav                                           &
                              ,namesat
!
      character(*)               , intent(in    ) :: infile
      real(r_kind)               , intent(in    ) :: gstime
      real(r_kind)               , intent(in    ) :: winlen
      real(r_kind)               , intent(in    ) :: despikeStdFac      
      real(r_kind), dimension(*) , intent(in    ) :: mbr_domain
      integer(i_kind)            , intent(in    ) :: MaLngth      
      logical                    , intent(in    ) :: do_superobs
      logical                    , intent(in    ) :: do_despiking
      logical                    , intent(in    ) :: do_ma
      real(r_kind), dimension(*) , intent(in    ) :: SprObsDim
!
!  local variable
      character (i25), parameter :: myname='read_satmar'
!
      integer(i_kind) :: cnt,ll,tot 
      integer(i_kind) :: ndata, iout
      
      integer(i_kind) :: ireadmg,ireadsb,idate
      integer(i_kind) :: nrec, irec, nmind
      character(len=i25) :: subset
!
      integer(i_kind), parameter:: nreal        = i8  !23
      character(i256), parameter:: hdr_fltJS2   = 'RSST AETP ASFL ADQF ALRF IPIN ODLE SAID'
      character(i256), parameter:: hdr_time     = 'YEAR MNTH DAYS HOUR MINU SECW'
      character(i256), parameter:: hdr_loc      = 'CLATH CLONH'
      character(i256), parameter:: hdr_howvJS2  = 'KBSW RKSW WSPA' !CBSW RCSW'
      character(i256), parameter:: hdr_biasJS2  = 'KNCS SBCK' !KBIC'
!      
      character(i256), parameter:: hdr_fltSAL   = 'RSST BSADQF NVPSWH IPIN ODLE SAID'
      character(i256), parameter:: hdr_howvSAL  = 'SBSWH RMSSWH WSPA'
      character(i256), parameter:: hdr_biasSAL  = 'SBCSB NISWH' !SBNIC' 
!
      character(i256),parameter:: hdr_timeCS2  = 'YEAR MNTH DAYS HOUR MINU SECO'
      character(i256),parameter:: hdr_fltCS2   = 'DSST ODLE L1PQ L1PF L2PF SAID'
      character(i256),parameter:: hdr_howvCS2  = 'SWHS S20SWHS WSPA'
!      character(i256),parameter:: hdr_howvCS2  = 'KBSW RKSW WSPA'
      character(i256),parameter:: hdr_biasCS2  = 'SBCK'
!
      character(i256), parameter:: hdr_howvNO  = 'HOWV SDWH WS10'
      character(i256), parameter:: hdr_locNO   = 'CLAT CLON'
      character(i256),parameter :: hdr_timeNO  = 'YEAR MNTH DAYS HOUR MINU SECO'
      
!      
      integer(i_kind) :: n_tm, n_howv, n_lc, n_fltJS2, n_fltSAL, n_fltCS2     &
                       , n_biasJS2, n_biasSAL, n_biasCS2
!                                                                        
      real(r_kind) :: depth, said_1d                                     
!
      real(r_kind), allocatable, dimension(:,:) :: data_all,data_sobs
      real(r_kind), allocatable, dimension(:)   :: dum_sort
      real(r_kind), allocatable, dimension(:)   :: flt_1dJS2, flt_1dSAL, flt_1dCS2  &
                  , time_1d, howv_1d, howv_1d_m1 , loc_1d, loc_1d_m1                &
                  , biasJS2, biasSAL, biasCS2
      integer(i_kind), allocatable, dimension(:):: time_1dMN, isort, iloc
      character(i25) , allocatable, dimension(:):: said, said_sobs
      character(i256)                           :: hdr_loc_dum
!  
      real(r_kind) :: rminobs
      real(r_kind) :: nsec, nsec_m1, dt_sec
      real(r_kind) :: Hdistm, Pdistm
!
!  Swords
      integer(i_kind),parameter :: howvMax = 15_i_kind
      integer(i_kind),parameter :: howvRatMiuSigma = 1./i3 
      integer(i_kind),parameter :: howvRathowvDpth = i2
      real(r_kind)   ,parameter :: howvDistm = 10000.0_r_kind
      real(r_kind)   ,parameter :: hs_max_diff = 0.2_r_kind
      real(r_kind)   ,parameter :: dflt_err = 0.35_r_kind
!
      logical :: bool_open
!
      n_tm     = count_kwrd(hdr_time)
      n_lc     = count_kwrd(hdr_loc)
      n_howv   = count_kwrd(hdr_howvJS2)
      n_fltJS2 = count_kwrd(hdr_fltJS2) 
      n_fltSAL = count_kwrd(hdr_fltSAL)
      n_fltCS2 = count_kwrd(hdr_fltCS2)
      n_biasJS2= count_kwrd(hdr_biasJS2)
      n_biasSAL= count_kwrd(hdr_biasSAL)
      n_biasCS2= count_kwrd(hdr_biasCS2)
!      
      if(.not.allocated(time_1d))   allocate(time_1d  (n_tm    ))
      if(.not.allocated(time_1dMN)) allocate(time_1dMN(n_tm-i1 ))
      if(.not.allocated(loc_1d))    allocate(loc_1d   (n_lc    ))
      if(.not.allocated(loc_1d_m1)) allocate(loc_1d_m1(n_lc-i1 ))
      if(.not.allocated(howv_1d))   allocate(howv_1d  (n_howv  ))
      if(.not.allocated(howv_1d_m1))allocate(howv_1d_m1(n_howv ))
      if(.not.allocated(flt_1dJS2)) allocate(flt_1dJS2(n_fltJS2))
      if(.not.allocated(flt_1dSAL)) allocate(flt_1dSAL(n_fltSAL))
      if(.not.allocated(flt_1dCS2)) allocate(flt_1dCS2(n_fltCS2))
      if(.not.allocated(biasJS2))   allocate(biasJS2  (n_biasJS2))
      if(.not.allocated(biasSAL))   allocate(biasSAL  (n_biasSAL))
      if(.not.allocated(biasCS2))   allocate(biasCS2  (n_biasCS2))
!
      irec     = i0
      nrec     = i0
      ndata    = i0
      said_1d  = i0
      tot      = i1
!
      if (allocated(sat_obs))  deallocate(sat_obs)
!
      inquire( unit=uid_rd, opened=bool_open )
      if(bool_open) call closbf(uid_rd)
!
      call file_existence(infile, myname) 
!
      call counter(uid_rd, infile, cnt) 
!
      allocate (data_all (nreal, cnt),isort(cnt),said(cnt))
!
      isort = i0
      data_all = r0
!
      open(uid_rd,file=trim(infile),action='read',form='unformatted') 
!
      call openbf(uid_rd,'IN',uid_rd)
!
      call datelen(i10)
!
      read_msg: do while(ireadmg(uid_rd,subset,idate) == 0)
         read_loop: do while (ireadsb(uid_rd) == 0)
            nrec = nrec + i1
            time_1d  = r0 
            howv_1d  = r0
            loc_1d   = r0
!
            if (    (index(trim(subset),trim(namesat(1))) > 0 )       &           !JS2
               .or. (index(trim(subset),trim(namesat(4))) > 0 )       )    then   !JS4
               biasJS2  = r0
               call ufbint(uid_rd,flt_1dJS2,n_fltJS2,1,irec,trim(hdr_fltJS2))
!  Qc flags
               if (  (flt_1dJS2(1).gt.1.)  .or. &
                     (flt_1dJS2(2).ne.0.)  .or. &
                     (flt_1dJS2(3).ne.0.)  .or. &
                     (flt_1dJS2(4).ne.0.)  .or. &
                     (flt_1dJS2(5).ne.0.)  .or. &
                     (flt_1dJS2(6).ne.0.)       )    cycle
               depth = abs(flt_1dJS2(7))
               said_1d = flt_1dJS2(n_fltJS2)
!  Time 
               call ufbint(uid_rd,time_1d,n_tm,1,irec,trim(hdr_time))
!  Howv
               call ufbint(uid_rd,howv_1d,n_howv,1,irec,trim(hdr_howvJS2))
               call ufbint(uid_rd,biasJS2,n_biasJS2,1,irec,trim(hdr_biasJS2))
!
            else if   (index(trim(subset),trim(namesat(2)))>0) then                 !SRLTK
               biasSAL  = r0
               call ufbint(uid_rd,flt_1dSAL,n_fltSAL,1,irec,trim(hdr_fltSAL))
               if ( (flt_1dSAL(1).gt.1. ) .or.   &
                    (flt_1dSAL(2).ne.0. ) .or.   &
                    (flt_1dSAL(3).le.30.) .or.   &
                    (flt_1dSAL(4).ne.0. )        ) cycle  
               depth = abs(flt_1dSAL(5))
               said_1d = flt_1dJS2(n_fltJS2)
!  Time 
               call ufbint(uid_rd,time_1d,n_tm,1,irec,trim(hdr_time))
!  Howv
               call ufbint(uid_rd,howv_1d,n_howv,1,irec,trim(hdr_howvSAL))
               call ufbint(uid_rd,biasSAL,n_biasSAL,1,irec,trim(hdr_biasSAL))
!               print*,'biasSAL',biasSAL
!               howv_1d(1) = howv_1d(1)+sum(biasSAL)
!
            else if   (index(trim(subset),trim(namesat(3)))>0) then                 !CS2
               biasCS2  = r0
               call ufbint(uid_rd,flt_1dCS2,n_fltCS2,1,irec,trim(hdr_fltCS2))
               if ( (flt_1dCS2(1) > 1    ) .or.   &
                    (flt_1dCS2(2) > 0    ) .or.   &
                    (flt_1dCS2(3) < 90.  ) .or.   &
                    (flt_1dCS2(4) /= 0.  ) .or.   &
                    (flt_1dCS2(5) /= 0.  )        ) cycle 
               depth = abs(flt_1dCS2(2))
               said_1d = flt_1dJS2(n_fltJS2)
!  Time 
               call ufbint(uid_rd,time_1d,n_tm,1,irec,trim(hdr_timeCS2))
!  Howv
               call ufbint(uid_rd,howv_1d,n_howv,1,irec,trim(hdr_howvCS2))
               howv_1d(1)=sqrt(howv_1d(1))
!               print*,'howv_1d: ',howv_1d
               call ufbint(uid_rd,biasCS2,n_biasCS2,1,irec,trim(hdr_biasCS2))
               if (sum(biasCS2)>10) biasCS2=0
            else if   (index(trim(subset),trim(namesat(5)))>0) then                 !S3a
!  Time 
               call ufbint(uid_rd,time_1d,n_tm,1,irec,trim(hdr_timeNO))
!  Howv
               call ufbint(uid_rd,howv_1d,n_howv,1,irec,trim(hdr_howvNO))
!               print*,'howv_1d: ',howv_1d
               if (sum(biasCS2)>10) biasCS2=0
               
            end if
! Exclude Observations with Wind>20 according to S3-SP-1322_3.pdf
            if (howv_1d(3) > 30_r_kind) cycle                             ! qc max(wind)  
!
!  Temporal space
            time_1dMN = int(time_1d(1:5))
            call w3fs21(time_1dMN, nmind)
            rminobs=real(nmind,r_kind)+(real(time_1d(6),r_kind)*r60inv)
            if (abs(gstime-rminobs)*r60inv <= winlen) then
!              space for further processing in temporal space, if required 
               if (dodebug1) then
                  print*, myname,int(time_1d(1:6))
                  print*, abs(gstime-rminobs)*r60inv 
                  print*, winlen
               end if 
            else 
               cycle                                                       !qc
            end if
!
!  Physical Space
            if  (index(trim(subset),trim(namesat(5)))>0) then  
               hdr_loc_dum=hdr_locNO
            else
               hdr_loc_dum=hdr_loc   
            end if
            call ufbint(uid_rd,loc_1d,n_lc,1,irec,trim(hdr_loc_dum))
!
            if(abs(loc_1d(1))>r90 .or. abs(loc_1d(2))>r360) cycle                                   !qc
!
            if (loc_1d(2)==r360  ) loc_1d(2)=loc_1d(2)-r360
            if (loc_1d(2)< r0    ) loc_1d(2)=loc_1d(2)+r360
!
            if (loc_1d(1)<mbr_domain(1).or.loc_1d(1)>mbr_domain(2)) cycle  !xcld
            if (loc_1d(2)<mbr_domain(3).or.loc_1d(2)>mbr_domain(4)) cycle  !xcld
!        
!  Slayer  
            if (dofilterS) then
               if (howv_1d(1)/=howv_1d(1)                                  ) cycle  !qc0
!               if (howv_1d(1)<=dflt_err                                    ) cycle  !qc01
               if (howv_1d(1)              > howvMax                       ) cycle  !qc1
               if (howv_1d(2)<=tiny_r_kind) howv_1d(2)=dflt_err
               if ((howv_1d(2)              > r0)            .and.       &
                  ( howv_1d(1) / howv_1d(2) < howvRatMiuSigma)             ) cycle  !qc2
               if ( howv_1d(1) / depth      > howvRathowvDpth              ) cycle  !qc3
!              howv_1d(1) = howv_1d(1)+sum(biasCS2)
!
!               call datesec(time_1d, nsec)
!               if (tot==i1) then
!                  nsec_m1 = nsec
!                  loc_1d_m1 = loc_1d
!                  howv_1d_m1 = howv_1d(1)
!               else if(tot>i1)then
!                  dt_sec = nsec-nsec_m1
!                  call lldistm(loc_1d,loc_1d_m1, Hdistm, Pdistm)
!                  if ( (abs(Hdistm)<howvDistm).and.(dt_sec<i2) ) then
!                     if (abs(howv_1d(1)-howv_1d_m1(1))/howv_1d_m1(1) > hs_max_diff  ) cycle    !qc4
!                     if (abs(howv_1d(1)-howv_1d_m1(1))/(dt_sec)**i2  > grav/i2      ) cycle    !qc5
!                  end if
!                  nsec_m1 = nsec
!                  loc_1d_m1 = loc_1d 
!                  howv_1d_m1 = howv_1d(1)
!               end if
            end if
!
           if (dodebug1) print*,myname, ndata,' KBSW=',howv_1d(1), ' RKSW=',howv_1d(2),' WSPA=',howv_1d(3),' KNCS=',howv_1d(4)
!  
            tot=tot+i1
            ndata=ndata+i1
            isort(nrec)=ndata
!
            data_all(1,ndata) = loc_1d(2)
            data_all(2,ndata) = loc_1d(1)
            data_all(3,ndata) = howv_1d(1)
            data_all(4,ndata) = howv_1d(2)
            data_all(5,ndata) = rminobs
            data_all(6,ndata) = r1
            data_all(7,ndata) = said_1d
            data_all(8,ndata) = howv_1d(3)
            said(ndata) = subset
!
         enddo read_loop
      enddo read_msg
!
      call closbf(uid_rd)
!
      if (dodebug1)then
         call writeMatrixdp('PureObs.txt', data_all)
      end if
!
      if(.not.allocated(iloc)) allocate(iloc(ndata))
      cnt = i0
      do ll = i1,size(data_all,2)
         if(isort(ll) > i0) then
            cnt = cnt + i1
            iloc(cnt) = isort(ll)
         end if
      end do
      if (.not.allocated(data_sobs)) allocate(data_sobs(nreal,ndata), said_sobs(ndata)) 
      do ll = i1,ndata
         iout=iloc(ll)
         data_sobs(1,ll)=data_all(1,iout)
         data_sobs(2,ll)=data_all(2,iout)
         data_sobs(3,ll)=data_all(3,iout)
         data_sobs(4,ll)=data_all(4,iout)
         data_sobs(5,ll)=data_all(5,iout)
         data_sobs(6,ll)=data_all(6,iout)
         data_sobs(7,ll)=data_all(7,iout)
         data_sobs(8,ll)=data_all(8,iout)
         said_sobs(  ll)=said    (  iout)
      end do
      if (dodebug1)then
         call writeMatrixdp('SuperObs.txt', data_sobs)
      end if
!
      if (do_superobs) then
          call super_obs(SprObsDim, mbr_domain                     &
             , shape(data_sobs(i3:i8,:)),data_sobs(i1,:), data_sobs(i2,:), data_sobs(i3:i8,:)   &
             , said_sobs, ndata)
      end if
!
      if(.not.allocated(dum_sort)) allocate(dum_sort(size(data_sobs(5,:))))
      dum_sort = data_sobs(5,:)
      call qsort(dum_sort,data_sobs)
      if (dodebug1)call writeMatrixdp('ObsSort.txt', data_sobs)
      if(allocated(dum_sort)) deallocate(dum_sort)
!
      if ((do_despiking).or. (do_ma)) then
         print*,trim(myname), 'Despiking Starts'
         call despike(data_sobs,ndata, do_despiking, despikeStdFac  &
                     ,do_ma, MaLngth, dodebug1)
         print*,trim(myname), 'Despiking Ends'         
       end if
!
      if (.not.allocated(sat_obs)) call allocate_alt(ndata)
      do ll = i1,ndata
            sat_obs(ll)%x_grd(1) =  data_sobs(1,ll)
            sat_obs(ll)%x_grd(2) =  data_sobs(2,ll)
            sat_obs(ll)%value    =  data_sobs(3,ll)
            sat_obs(ll)%oerr     =  data_sobs(4,ll)
            sat_obs(ll)%hour     =  data_sobs(5,ll)
            sat_obs(ll)%qkey     =  data_sobs(6,ll)
            sat_obs(ll)%wspa     =  data_sobs(8,ll)
            sat_obs(ll)%typ      =  said_sobs(ll)
            sat_obs(ll)%kept     = .true.
      end do
!
      if(allocated(time_1d)   )   deallocate(time_1d   )
      if(allocated(time_1dMN) )   deallocate(time_1dMN )
      if(allocated(loc_1d)    )   deallocate(loc_1d    )
      if(allocated(loc_1d_m1) )   deallocate(loc_1d_m1 )
      if(allocated(howv_1d)   )   deallocate(howv_1d   )
      if(allocated(howv_1d_m1))   deallocate(howv_1d_m1)
      if(allocated(flt_1dJS2) )   deallocate(flt_1dJS2 )
      if(allocated(flt_1dSAL) )   deallocate(flt_1dSAL )
      if(allocated(flt_1dCS2) )   deallocate(flt_1dCS2 )
      if(allocated(data_all)  )   deallocate(data_all  )
      if(allocated(data_sobs) )   deallocate(data_sobs )      
      if(allocated(isort)     )   deallocate(isort     )
      if(allocated(iloc)      )   deallocate(iloc      )
      if(allocated(said)      )   deallocate(said      )
      if(allocated(said_sobs) )   deallocate(said_sobs )
!
   end subroutine read_satmar
!
!-----------   
!> @brief file_existence Subroutine checks if a file exist
!! @param[in] filename name of the file to be cheched
!! @param[in] myname_in is the name of the calling subroutine
!<
   subroutine file_existence(filename, myname_in)
!      integer(i_kind), intent(in)      :: lenflnm
      character(len=*), intent(in)   :: filename
      character(len=*), intent(in)   :: myname_in
      logical :: ex
      character (i25),parameter :: myname='file_existence'
   !
      inquire(file=trim(filename),exist=ex)
      if (.not.ex) then
         print*, trim(myname),' : Input file ', trim(filename)&
               , ' does not exist! Instant Termination! Called by ',trim(myname_in)                  
         stop
       end if
   ! 
   end subroutine file_existence
!
!-----------   
!> @brief The subroutine Counter counts (xoxo) total number of messages into a bufr file
!! @param[in] uid_rd is the open unit
!! @param[in] infile is the name of the file to be read
!! @param[out] cnt is the total number of messages 
!<
   subroutine counter(uid_rd,infile,cnt)
      integer(i_kind), intent(out)   :: cnt
      integer(i_kind), intent(in )   :: uid_rd
      character(*),    intent(in )   :: infile 
! local
      integer*8   :: uid_rd1
      character(len=i25) :: subset1
      integer(i_kind) :: ireadmg, ireadsb, idate, ierr
      character(len=i25) :: myname='counter'     
!
      uid_rd1 = int8(uid_rd)
      open(uid_rd1,file=trim(infile),action='read',form='unformatted', iostat=ierr)
      if (ierr/=0) print*, trim(myname),': Error Opening the ', trim(infile) 
!
      call openbf(uid_rd1,'IN',uid_rd1)
!   
      cnt = 0
      do while(ireadmg(uid_rd1,subset1,idate) == 0)
         do while (ireadsb(uid_rd1) == 0)
            cnt = cnt+1
         end do
      end do
      call closbf(uid_rd1)
!
   end subroutine counter
!>
!! @brief The datesec reads calculates the number of seconds since 00:00:00,
!!     1 January 1978 (RefDate : in days)
!! @author stelios flampouris (stylianos.flampouris@noaa.gov)
!! @date 2016.03.08  
!!
!! param[in] idate       - size(idate) = 6, see below
!! param[out] nse        - number of seconds 
!!     
!! @details [YYYY,MM,DD]-->Convert to Julian Days (JD)-->NDays=JD-Refdate-->Convert to
!! Seconds (NDaysInSec)--> nsec=NDaysInSec+HH*3600+MN*60+SS
!!
!! => Input Arguments    : idate real array with size 6: 
!!                         idate(1)=YYYY   !Year
!!                         idate(2)=MM     !Month
!!                         idate(3)=DD     !Day
!!                         idate(4)=HH     !Hour
!!                         idate(5)=MN     !Minute
!!                         idate(6)=SS     !Second
!<
   subroutine datesec(idate, nsec)
      real(r_kind), intent(in) :: idate(6)
      real(r_kind), intent(out) :: nsec
!
!  local parameters 
      real(r_kind), parameter :: ReFDate = 2443510.
      real(r_kind) :: JD
!  Initializing variables
      nsec = 0
!
!  Convert idate(1:3) to JD
      JD = idate(3) - 32075                                                       &
                    + 1461 * ( idate(1) + 4800 + (idate(2) - 14) / 12) / 4        &
                    +  367 * ( idate(2) - 2    - (idate(2) - 14) / 12 * 12) / 12  &
                    -    3 * ((idate(1) + 4900 + (idate(2) - 14) / 12) / 100) / 4
!
!  Number of days from the reference days
      JD = JD - ReFDate
!
!  Number of seconds
      nsec = (JD* 86400) + (idate(4) * 3600) + (idate(5) * 60) + idate(6)
!
   end subroutine datesec
!
!------
!>
!!@param[in] latlon1 latitude & longtitude of point 1
!!@param[in] latlon2 latitude & longtitude of point 2
!!@param[out] Hdistm Haversine distance between the point1 and point2 in meters(m)
!!@param[out] Pdistm Pythagorian distance between the point1 and point2 in meters(m)
!!@brief Subroutine lldistm : Calculates the distance in m between two points in
!!spherical coordinates (Earth in this case).
!!@details Fortran code for Haversine and Pythagorian distance 
!!@author stelios flampouris (stylianos.flampouris@noaa.gov)
!!
!@date 2016.03.08            : stelios - Prototype 
!@date 2016.08.24            : stelios - Fully compatible with GSI
!@date 2016.10.15            : stelios - Fully compatible with LETKF
!<

   subroutine lldistm(latlon1,latlon2, Hdistm, Pdistm)
      use constants_hs, only : deg2rad,pi,rearth
      implicit none
      real(r_kind), dimension(2), intent(in) :: latlon1, latlon2
      real(r_kind), intent(out):: Hdistm, Pdistm
!  local variables 
      real(r_kind) :: lat1, lat2, lon1, lon2, dLat, dLon, dum1, dum2, x, y
!
      lat1=latlon1(1)*deg2rad
      lat2=latlon2(1)*deg2rad;
      lon1=latlon1(2)*deg2rad;
      lon2=latlon2(2)*deg2rad;
      dLat=lat2-lat1;
      dLon=lon2-lon1;
!  Haversine distance
      dum1=sin((dLat)/2)**2 + cos(lat1)*cos(lat2) * sin(dLon/2)**2;
      dum2=2*atan2(sqrt(dum1),sqrt(1-dum1))
      Hdistm=rearth*dum2
!  Pythagoran distance
      x=dLon*cos((lat1+lat2)/2);
      y=dLat;
      Pdistm=rearth*sqrt(x*x + y*y);
!   
   end subroutine lldistm
!
!-----------   
!> @brief count_kwrd counting total number of words into a string
!! @param[in] text is the string
!! count_kwrd the total number of words 
!<

function count_kwrd(text)
      use constants_hs, only: i0, i1
      character (len=*), intent(in) :: text
      integer(i_kind)               :: count_kwrd
! local
      integer              :: nwords, pos, i
!
      pos = i1
      count_kwrd = i0
      loop: do
         i = verify(text(pos:), ' ')   !-- Find next non-blank.
         if (i == i0) exit loop        !-- No word found.
         count_kwrd = count_kwrd + 1   !-- Found something.
         pos = pos + i - i1            !-- Move to start of the word.
         i = scan(text(pos:), ' ')     !-- Find next blank.
         if (i == i0) exit loop        !-- No blank found.
         pos = pos + i - i1            !-- Move to the blank.
      end do loop
!
   end  function count_kwrd

end module read_obs_mod



!# Graveyard
!

!!! Thinning
!!   !! 
!!   !!  Dragon on Diet
!!   !         cnt = 0
!!   !
!!   !         iuse=icuse(nc)
!!   !         if (ithin > 0 .and. iuse >=0) then
!!   !            ntmp=ndata
!!   !            if (thin4d) then
!!   !               timedif = r0                ! crit1=0.01_r_kind
!!   !            else
!!   !               timedif=abs(t4dv-toff)        ! timedif = 6.0_r_kind*abs(tdiff) & crit1=0.01_r_kind+timedif
!!   !            end if
!!   !            crit1 = timedif/r6+half
!!   !!
!!   !            call map3grids(-1,0,DumForThin,nlevp,dlat_earth,dlon_earth     &
!!   !                           ,i1 ,crit1,ndata,iout,nrec,iiout,luse,.false.,.false.)
!!   !
!!   !               if (.not. luse) cycle 
!!   !               if(iiout > 0) isort(iiout)=0
!!   !               if (ndata > ntmp) then
!!   !                  nodata=nodata+1
!!   !               endif
!!   !               isort(nrec)=iout
!!   !         else  ! - no thinnning
!!   !               ndata=ndata+1
!!   !               nodata=nodata+1
!!   !               iout=ndata
!!   !               isort(nrec)=iout
!!   !         endif 


!!!GSI format
!!   !         usage = r0 !-  Set usage variable :: practically useless

!!   !!         
!!   !         call ufbint(uid_rd,rstation_id,1,1,irec,hdr_station)
!!   !!
!!   !         data_all(1,iout) = howv_1d(2)                 ! significant wave height error (m)
!!   !         data_all(2,iout) = dlon                       ! grid relative longitude
!!   !         data_all(3,iout) = dlat                       ! grid relative latitude
!!   !         data_all(4,iout) = r0                !105.0       ! pressure (in cb)
!!   !         data_all(5,iout) = howv_1d(1)                 ! significant wave height (in m)
!!   !         data_all(6,iout) = rstation_id                ! station id
!!   !         data_all(7,iout) = t4dv                       ! time
!!   !         data_all(8,iout) = nc                         ! type
!!   !         data_all(9,iout) = 0_r_kind                   ! quality mark
!!   !         data_all(10,iout) = 0.2_r_kind                ! original obs error (m)
!!   !         data_all(11,iout) = usage                     ! usage parameter
!!   !         if (lhilbert) thisobtype_usage=11         ! save INDEX of where usage is stored for hilbertcurve cross validation (if requested)
!!   !         data_all(12,iout) = r0                      ! dominate surface type
!!   !         data_all(13,iout) = 295_r_kind                ! skin temperature
!!   !         data_all(14,iout) = 1.0                         ! 10 meter wind factor
!!   !         data_all(15,iout) = bmiss                     ! surface roughness
!!   !         data_all(16,iout) = dlon_earth*rad2deg        ! earth relative longitude (degrees)
!!   !         data_all(17,iout) = dlat_earth*rad2deg        ! earth relative latitude (degrees)
!!   !         data_all(18,iout) = r0                      ! station elevation (m)
!!   !         data_all(19,iout) = r0                      ! observation height (m)
!!   !         data_all(20,iout) = r0  !-depth                    ! terrain height at ob location
!!   !         data_all(21,iout) = 100000000000.000                     ! provider name !r_prvstg(1,1)
!!   !         data_all(22,iout) = 100000000000.000                     ! subprovider name !r_sprvstg(1,1)
!!   !         data_all(23,iout) = 0_r_kind                     ! cat 
!!   !!


!!            dlon_earth=loc_1d(2)*deg2rad
!!            dlat_earth=loc_1d(1)*deg2rad
!!            nread = nread + 1 
!!   !         
!!   !         if(regional)then                                            ! Regional 
!!   !            call tll2xy(dlon_earth,dlat_earth,dlon,dlat,outside)     ! *Convert to rotated coordinate
!!   !            if(outside) cycle
!!   !         else                                                        ! Global case
!!   !            dlon=dlon_earth
!!   !            dlat=dlat_earth
!!   !            call grdcrd1(dlat,rlats,nlat,1)
!!   !            call grdcrd1(dlon,rlons,nlon,1)
!!   !         endif
!

