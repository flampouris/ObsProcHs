!> @file general_functions.f90
!! @brief The module includes subroutines and functions used multiple times throught
!! out the program
!! 
!! @author stelios flampouris
!! @version 0.2
!!
!> @copyright GNU (See the main)
!<
!!@date 26-Oct-2016
!!@date 17-Nov-2016 by stelios. The writeMatrix(dp) subroutines were added
!!@todo Clean up all the modules... Maybe, in the future. 
!<
module general_func
   use kinds, only : i_kind,r_kind
   use constants_hs, only : i25,i256
   use variables_obsop_hs, only : obs_data, out_data, dodebug1
!
   implicit none
   private
   public collocation,write_obs,writeMatrix,writeMatrixdp,write_error_stat,min2date,writeArray 
   contains
!
!>
!! @brief The collocation collocates (Original Name hoho) oservations and model
!! @param[in] nlon        - number of grid points in lon
!! @param[in] nlat        - number of grid points in lat
!! @param[in] nv2d        - number variables read from the model output
!! @param[in] long        - gridded longitude 
!! @param[in] latg        - gridded latitude
!! @param[in] mbr_data    - gridded model data from one grib2 file
!! @param[in] nobs        - number of observations
!! @param[in] all_obs     - all the observations (structure, check variables_obs mod)
!! @param[out] collocated - The collocated information (see the variables_obs for details)
!! @param[out] cnt		  - The number of the collocated data
!!
!! @details 
!! @date 17-Nov-2016 Debugging Added
!! @date 20-Nov-2016 stelios/Jessica M. Exclude the no wave prediction points from the collocation
!! @date 20-Jan-2017 Stelios/Optimized for cartesian grid
!!
!! @todo Unstructured grid : Some code is here
!<
!! 
!
   subroutine collocation(nlon, nlat, nv2d, latg, long, mbr_data           &
                         ,nobs, all_obs, collocated, cnt)
      use constants_hs, only : i0, i1, i2,i3, i4                           &
                              ,r0
! 
      integer(i_kind),                            intent(in) :: nlon
      integer(i_kind),                            intent(in) :: nlat
      integer(i_kind),                            intent(in) :: nv2d
      real(r_kind)   , dimension(nlon,nlat     ), intent(in) :: long
      real(r_kind)   , dimension(nlon,nlat     ), intent(in) :: latg
      real(r_kind)   , dimension(nlon,nlat,nv2d), intent(in) :: mbr_data
      integer(i_kind),                            intent(in) :: nobs
      type (obs_data), dimension(nobs)          , intent(in) :: all_obs 
      type (out_data), dimension(nobs)          , intent(inout) :: collocated

! local
      character(i25), parameter :: myname='collocation'
      integer(i_kind) :: ll,k1,k2,cnt
      integer(i_kind), dimension(i2)        :: ind1, i_Md
      integer(i_kind), dimension(i4, i2)    :: keep_ind
      integer(i_kind), dimension(nobs)      :: isort
      
      real(r_kind)   , dimension(nlon,nlat) :: hdist,latg0180
      real(r_kind)   , dimension(i4,nv2d  ) :: keep_val
      real(r_kind)   , dimension(i4       ) :: keep_lon,keep_lat,keep_dist
      real(r_kind)   , dimension(i2       ) :: obs_ind
      real(r_kind)                          :: dx,dy
      type (out_data), dimension(nobs)      :: dumcollocated
!
    isort = i0
    cnt = i0
    latg0180=latg+90
!   
    i_Md=(/ nint(real(nlon/i2)), nint(real(nlat/i2)) /)
!   
    dx=abs(long(i_Md(i1)+i1,i_Md(i2))-long(i_Md(i1),i_Md(i2)+i1))
    dy=abs(latg(i_Md(i1)+i1,i_Md(i2))-latg(i_Md(i1),i_Md(i2)+i1))

    if (dodebug1)print*,'At ',trim(myname),' dx=',dx,' dy=',dy, ' nobs=', nobs
!
    loop_obs : do ll = i1,nobs,i1
         dumcollocated(ll)%x_grd      =r0 
         dumcollocated(ll)%ind_grd    =r0
         dumcollocated(ll)%obs_value  =r0
         dumcollocated(ll)%obs_oerr   =r0
         dumcollocated(ll)%timestamp  =r0
         dumcollocated(ll)%wspa       =r0
         dumcollocated(ll)%mbr_value  =r0
!
         obs_ind =(/((all_obs(ll)%x_grd(1)      - minval(long    ))/dx)+i1    &
                  , ((all_obs(ll)%x_grd(2) + 90 - minval(latg0180))/dy)+i1    /)
         ind1=nint(obs_ind)
!         
! Start Code for unstructured grid
!         do k1 = i1, nlon, i1
!            do k2 = i1, nlat, i1
!               hdist(k1,k2) = hdistance (all_obs(ll)%x_grd(1), all_obs(ll)%x_grd(2), long(k1,k2),latg(k1,k2))
!            end do
!         end do
!         ind1 = minloc(hdist)
! End Code for unstructured grid
!
         if (mbr_data(ind1(i1), ind1(i2),i1)>100) then 
            if (dodebug1) print*,trim(myname), mbr_data(ind1(i1), ind1(i2),i1)
            cycle loop_obs
         end if
         if (long(ind1(i1),ind1(i2))-all_obs(ll)%x_grd(1)<=0) then 
            keep_ind(i1   ,i1) = ind1(i1)
            keep_ind(i2:i3,i1) = ind1(i1)+i1
            keep_ind(i4   ,i1) = ind1(i1)
         else
            keep_ind(i1   ,i1) = ind1(i1)-i1
            keep_ind(i2:i3,i1) = ind1(i1)
            keep_ind(i4   ,i1) = ind1(i1)-i1
         end if
!
         if (latg(ind1(i1),ind1(i2))-all_obs(ll)%x_grd(2)<=0) then 
            if (all_obs(ll)%x_grd(i2)<=0) then
               keep_ind(i1:i2,i2) = ind1(i2) 
               keep_ind(i3:i4,i2) = ind1(i2)+i1
            else
               keep_ind(i1:i2,i2) = ind1(i2)
               keep_ind(i3:i4,i2) = ind1(i2)+i1
            end if
         else
            if (all_obs(ll)%x_grd(2)>0) then
               keep_ind(i1:i2,i2) = ind1(i2)-i1
               keep_ind(i3:i4,i2) = ind1(i2)
            else
               keep_ind(i1:i2,i2) = ind1(i2)
               keep_ind(i3:i4,i2) = ind1(i2)-i1
            end if
         end if
!
         do k1 = i1,i4,i1
            if ((mbr_data(ind1(i1),ind1(2), 1 )<0.01  )  &
               .and.   (dodebug1)                    )  then
               print*,myname, ind1,long(ind1(1), ind1(2)), latg(ind1(1),ind1(2)) 
               print*,mbr_data(ind1(i1),ind1(2), 1 )
               print*,myname, ind1,long(ind1(2), ind1(1)), latg(ind1(2),ind1(1))
               print*,myname,'shape mbr_data',shape(mbr_data) 
               print*,mbr_data(ind1(2),ind1(1), 1 )
            end if
!
            do k2 = i1,nv2d,i1
               keep_val(k1, k2) = mbr_data(keep_ind(k1,i1), keep_ind(k1,i2), k2 )
               if (keep_val(k1, k2) >40) cycle loop_obs
            end do
            keep_lon(k1  ) = long    (keep_ind(k1,i1), keep_ind(k1,i2)     ) 
            keep_lat(k1  ) = latg    (keep_ind(k1,i1), keep_ind(k1,i2)     )
            keep_dist(k1 ) = hdist   (keep_ind(k1,i1), keep_ind(k1,i2)     )
         end do
!        
! Start Code for unstructured grid
!         call real_indices (size(keep_ind,1),size(keep_ind,2), keep_ind       &
!                          , keep_lon, keep_lat                                &
!                          , size(all_obs(ll)%x_grd), all_obs(ll)%x_grd        &
!                          , obs_ind)
!         collocated(ll)%ind_grd(:) = obs_ind
! End Code for unstructured grid
!
         do k2 = i1,nv2d,i1
            call bilnr_int (size(keep_lon,1), keep_lon, keep_lat, keep_val(:,k2) &
                           , size(all_obs(ll)%x_grd),all_obs(ll)%x_grd           &
                           , collocated(ll)%mbr_value(k2))
         end do
!
         collocated(ll)%x_grd(:)   = all_obs(ll)%x_grd
         collocated(ll)%obs_value  = all_obs(ll)%value
         collocated(ll)%obs_oerr   = all_obs(ll)%oerr
         collocated(ll)%timestamp  = all_obs(ll)%hour
         collocated(ll)%typ        = all_obs(ll)%typ
         collocated(ll)%wspa       = all_obs(ll)%wspa
         collocated(ll)%ind_grd(:) = obs_ind
!
         cnt = cnt+1
         isort(ll) = cnt
      end do loop_obs
!
      cnt = i0
      do ll = i1,size(isort,1)
         if(isort(ll) > i0) then
            cnt = cnt + i1
            dumcollocated(cnt) = collocated(ll)
         end if

      end do
      if(dodebug1) print*, trim(myname),'Number of Collocated Obs:', cnt
! 
      collocated = dumcollocated
   end subroutine collocation
!
!------
!!@brief The bilnr_int subroutine is the simplest implementation of the bilinear interpolation 
!!@param[in] nx number of grid points in x (and y in this case)
!!@param[in] grd_lon longitude of the 4 points around the obs point // dimension(nx)
!!@param[in] grd_lat latitude of the 4 points around the obs point // dimension(nx)
!!@param[in] grd_val values at the 4 points of the grid
!!@param[in] no number of observations  
!!@param[in] obs_coor Coordinates of the observation point // dimension(no)
!!@param[out] obs_val the computed value of wahatever at obs location
!! 
!! @details See numerical recipes p.96
!!
!! @author stelios flampouris
!! @version 0.1
!! @copyright GNU 
!! @date 28-Oct-2016 
!<
   subroutine bilnr_int (nx, grd_lon, grd_lat, grd_val, no, obs_coor, obs_val)
      use constants_hs, only : i1,i2,i3,i4
      integer(i_kind)               , intent(in)  :: nx,no
      real(r_kind)   , dimension(nx), intent(in)  :: grd_lon, grd_lat, grd_val
      real(r_kind)   , dimension(no), intent(in)  :: obs_coor
      real(r_kind)                  , intent(out) ::obs_val
! local
      character(i25), parameter :: myname='bilnr_int'
      real(r_kind) :: t,u
!
      t = (obs_coor(i1)-minval(grd_lon))/(maxval(grd_lon)-minval(grd_lon))
      u = (obs_coor(i2)-minval(grd_lat))/(maxval(grd_lat)-minval(grd_lat))
!
      obs_val = (i1-t)*(i1-u)*grd_val(i1)    &
              + t     *(i1-u)*grd_val(i2)    &
              + t     *u     *grd_val(i3)    &
              + (i1-t)*u     *grd_val(i4)
!
   end subroutine bilnr_int
!
!-----
!> @brief The subroutine real_indices is defining the "indices of the obseravion
!!point"
!!@details The obs point is not located at the grid therefore, its coorinates
!!are real (NOT integer). This is useful for the DA.
!!
!!@param[in] nx, ny dimensions of grd_ind
!!@param[in] grd_ind Indices of the 4 points around the obs point // dimension(nx,ny)
!!@param[in] grd_lon longitude of the 4 points around the obs point // dimension(nx)
!!@param[in] grd_lat latitude of the 4 points around the obs point // dimension(nx)
!!@param[in] no dimensions of obs_coor
!!@param[in] obs_coor Coordinates of the observation point // dimension(no)
!!@param[out] obs_ind Indices of the observation point // dimension(no)
!!
!! @author stelios flampouris
!! @version 0.1
!! @copyright GNU 
!! @date 28-Oct-2016 
!<
   subroutine real_indices (nx,ny, grd_ind, grd_lon, grd_lat, no, obs_coor, obs_ind)
      use constants_hs, only : i1,i2
      integer(i_kind)                  , intent(in) :: nx,ny,no
      integer(i_kind), dimension(nx,ny), intent(in) :: grd_ind
      real(r_kind)   , dimension(nx   ), intent(in) :: grd_lon, grd_lat
      real(r_kind)   , dimension(no   ), intent(in) :: obs_coor
      real(r_kind)   , dimension(no   ), intent(out):: obs_ind
!  local
      character(i25), parameter :: myname='real_indices'
!
      obs_ind(i1) = minval(grd_ind(:,i1))+(obs_coor(i1)-minval(grd_lon))/(maxval(grd_lon)-minval(grd_lon))
      obs_ind(i2) = minval(grd_ind(:,i2))+(obs_coor(i2)-minval(grd_lat))/(maxval(grd_lat)-minval(grd_lat))
   end subroutine 
!
!-----
!> @brief The hdistance function calculates the Haversine distance between the 
!! point1 and point2 in meters(m)
!! @param[in] lon1        - longitude of point1 in deg
!! @param[in] lat1        - latitude of point1 in deg
!! @param[in] lon2        - longitude of point2 in deg
!! @param[in] lat2        - latitude of point2 in deg
!!
!! @author stelios flampouris
!! @version 0.1
!! @copyright GNU 
!! @date 26-Oct-2016 
!<
   function hdistance (lon1, lat1, lon2,lat2)
      use constants_hs, only : deg2rad, rearth
      real(r_kind) :: hdistance
      real(r_kind),intent(in) :: lat1, lat2, lon1, lon2
!  local variables
      real(r_kind) :: dLat, dLon, dum1, dum2, lat1in, lat2in, lon1in, lon2in
!
      hdistance=0.
      lat1in=lat1*deg2rad
      lat2in=lat2*deg2rad
      lon1in=lon1*deg2rad
      lon2in=lon2*deg2rad
      dLat=lat2in-lat1in
      dLon=lon2in-lon1in
!
      dum1=sin((dLat)/2)**2 + cos(lat1in)*cos(lat2in) * sin(dLon/2)**2;
      dum2=2*atan2(sqrt(dum1),sqrt(1-dum1))
      hdistance=rearth*dum2

!       dum1 = dLon*cos((lat1in+lat2in)/2)
!       dum2 = dLat;
!       hdistance=rearth*sqrt(dum1**2 + dum2**2)

   end function hdistance
!
!----
!
!> @brief The subroutine read_obs imports the data from a binary file.
!! param[in] all_obs
!!
!> @details
!! 1. From the structure to real 
!! 2. Write a binary file
!! alt_ID  YYYY MM DD HH mm ss lat lon mod_hs mod_ws altim_hs altim_ws indRealX indReaY
!! 
!<
   subroutine read_obs!(cfile,nn,elem,rlon,rlat,rlev,odat,oerr,ohx,oqc,obhr)
   use constants_hs, only : uid_rd,i256
   use kinds, only:r_single
   integer(i_kind) :: nobs,ndim,i
   real(r_single),dimension(:,:),allocatable::wk
   
   character(i256)::inputflnm
!local
   character(i25),parameter  :: myname='read_obs'
   inputflnm='obsophs_2017012012.20.dat'

   open(unit=uid_rd,FILE=inputflnm,FORM='unformatted',ACCESS='stream')
   read(unit=uid_rd,pos=1) nobs
   read(unit=uid_rd,pos=5) ndim
   
   if (dodebug1) print*,trim(myname), 'nobs:',nobs,'ndim:',ndim 
!
   if (.not.allocated(wk)) allocate(wk(nobs,ndim))
   do i=1,nobs
      read(unit=uid_rd) wk(i,:)
      if (dodebug1)print*,trim(myname),' ',wk(i,:)
   end do
!    elem(n) = REAL(wk(1),r_size)
!    rlon(n) = REAL(wk(2),r_size)
!    rlat(n) = REAL(wk(3),r_size)
!    rlev(n) = REAL(wk(4),r_size)
!    odat(n) = REAL(wk(5),r_size)
!    oerr(n) = REAL(wk(6),r_size)
!    ohx(n)  = REAL(wk(7),r_size)
!    oqc(n)  = NINT(wk(8))

   close(unit=uid_rd)
!
   end subroutine read_obs
   
!
!----
!
!> @brief The subroutine write_obs exports the data at a binary file.
!! param[in] all_obs
!!
!> @details
!! Steps: 
!! 1. From the structure to real 
!! 2. Write a binary file
!! alt_ID  YYYY MM DD HH mm ss lat lon mod_hs mod_ws altim_hs altim_ws indRealX indReaY observation_typical_error
!! alt_ID  : Sensor ID
!! YYYY MM DD HH mm ss : time of data acquisition
!! lat lon 
!! altim_hs : observed hs 
!! mod_hs 
!! altim_ws : 
!! mod_ws : 
!! indRealX : decimal index in X direction for the obs at the model grid 
!! indReaY  : decimal index in Y direction for the obs at the model grid 
!! observation_typical_error  : error variance of observation.
!! 
!<
   subroutine write_obs(collocated,nobs,cfile,dotextout)
      use variables_obsop_hs,only:out_data
      use constants_hs, only : i1,i2,i6,r0, uid_rd
      use kinds, only:r_single
      type (out_data), dimension(:), intent(in) :: collocated
      integer(i_kind),intent(in)::nobs
      character(*), intent(in) :: cfile
      logical, intent(in) :: dotextout
!local
      character(i25),parameter  :: myname='write_obs'
      integer(i_kind) :: kk
      integer(i_kind),dimension(i6)    :: idate
      real(r_single),dimension(16) :: dataout
      character(7) :: dumstr
      integer :: mypos
!
      open(unit=uid_rd,file=trim(cfile)//'.dat',form='unformatted',access='stream')
      if (dotextout) open(unit=uid_rd+1,file=trim(cfile)//'.txt',form='formatted')
!
      write(unit=uid_rd) nobs
      if (dodebug1)then 
         inquire(uid_rd,POS=mypos)
         print*, 'Mypos 1:',trim(myname),' ', mypos
      end if
      write(unit=uid_rd) size(dataout,1) 
      if (dodebug1)then 
         inquire(uid_rd,POS=mypos)
         print*, 'Mypos 2:',trim(myname),' ', mypos
      end if

      do kk=i1, nobs, i1
         dumstr=trim(collocated(kk)%typ(4:i25))
         if (index(trim(collocated(kk)%typ),'NC031')>0) then
            call min2date(collocated(kk)%timestamp,idate)
            read(dumstr(1:5),*)dataout(1) 
            dataout(2:7)=idate
            dataout(8:9)=real(collocated(kk)%x_grd,r_single)
            dataout(10)=real(collocated(kk)%obs_value,r_single)
            dataout(11)=real(collocated(kk)%mbr_value(i1),r_single)
            dataout(12)=real(collocated(kk)%wspa,r_single)
            dataout(13)=real(collocated(kk)%mbr_value(i2),r_single)
            dataout(14:15)=real(collocated(kk)%ind_grd,r_single)
            dataout(16)=real(collocated(kk)%obs_oerr,r_single)
            if (dodebug1) print*, kk, dataout
            write(unit=uid_rd) dataout
            if (dotextout)  write(uid_rd+1,*) dataout
         end if
      end do
      close(unit=uid_rd)
      close(unit=uid_rd+1)
   end subroutine write_obs 
!
!---------
!> @brief The subroutine secdate inverts the way that rnminobs is calculated
!! param[in] nsec from 1978
!! param[out] idate(6) [YYYY,MM,DD,HH, mm, ss]
!> @details Creates a size 6 array with the date as described at the idate 
!<
   subroutine min2date(nsec,idate)
      real(r_kind), intent(in) :: nsec
      integer(i_kind),dimension(*), intent(out) :: idate
! local
      character(i25), parameter :: myname='min2date '
      real(r_kind)              :: dum
      integer(i_kind)           ::  iyear,month,iday,idaywk,idayyr

      dum = mod(nsec,60.)
      idate(5) = dum
      idate(6) = (dum - idate(5))*60 
      dum = (nsec - idate(5))/60
      idate(4) = mod(dum,24.)
      dum = dum-idate(4)
      dum = dum/24. + 2443510.
!
      call w3fs26 (int(dum),idate(1),idate(2),idate(3),idaywk,idayyr )
!      
   end subroutine min2date
!
!---------
!> @brief Writes a 1D array to text file, column by column
!! @param [in] fileName path to the output file
!! @param [in] rda_A 1D array to write, rda_A  is real as defined at the program
!! (e.g. r_kind)
!! @details  Self-Explenatif, no?!?
!<
   subroutine writeArray(fileName, rda_A)
      real(r_kind), dimension(:), intent(in) :: rda_A
      character(len=*)         , intent(in) :: fileName
      integer ib_i, ib_j, il_ios
      integer, parameter :: ip_fid = 41
      open( unit = ip_fid, file = fileName, status = 'replace', form = 'formatted', iostat = il_ios)
      if (il_ios /= 0) print*,'In writeMatrix : Error creating file'//fileName
         do ib_i = 1, size(rda_A,1)
            write(unit=ip_fid, fmt='(E18.8, $)') rda_A(ib_i)!$ to avoid end of line
            write(unit=ip_fid, fmt=*)''! just for adding end of line
         end do
      write(unit=ip_fid, fmt=*)''
      close(ip_fid )
    end subroutine writeArray
!
!---------
!> @brief Writes a 2D array to text file, column by column
!! @param [in] fileName path to the output file
!! @param [in] rda_A 2D array to write, rda_A  is real as defined at the program
!! (e.g. r_kind)
!! @details  Self-Explenatif, no?!?
!<
   subroutine writeMatrixdp(fileName, rda_A)
      real(r_kind), dimension(:, :), intent(in) :: rda_A
      character(len=*)         , intent(in) :: fileName
      integer ib_i, ib_j, il_ios
      integer, parameter :: ip_fid = 41
      open( unit = ip_fid, file = fileName, status = 'replace', form = 'formatted', iostat = il_ios)
      if (il_ios /= 0) print*,'In writeMatrix : Error creating file'//fileName
         do ib_j = 1, size(rda_A,2)
            do ib_i = 1, size(rda_A,1)
               write(unit=ip_fid, fmt='(E18.8, $)') rda_A(ib_i,ib_j)!$ to avoid end of line
            end do
            write(unit=ip_fid, fmt=*)''! just for adding end of line
         end do
      close(ip_fid )
    end subroutine writeMatrixdp
!
!---------
!!> @brief Writes a 2D array to text file, column by column
!! @param [in] fileName path to the output file,rda_A  is real as defined during
!! compilation
!! @param [in] rda_A 2D array to write
!! @details Self-Explenatif, no?!?
!<
   subroutine writeMatrix(fileName, rda_A)
      real, dimension(:, :), intent(in) :: rda_A
      character(len=*)         , intent(in) :: fileName
      integer ib_i, ib_j, il_ios
      integer, parameter :: ip_fid = 41
!      logical :: iAmMaster
!
!        iAmMaster = parallel_amIMaster()
        !
!        if(iAmMaster)then
      open( unit = ip_fid, file = fileName, status = 'replace', form = 'formatted', iostat = il_ios)
      if (il_ios /= 0) print*,'In writeMatrix : Error creating file'//fileName
         do ib_j = 1, size(rda_A,2)
            do ib_i = 1, size(rda_A,1)
               write(unit=ip_fid, fmt='(E18.8, $)') rda_A(ib_i,ib_j)!$ to avoid end of line
            end do
            write(unit=ip_fid, fmt=*)''! just for adding end of line
         end do
      close(ip_fid )
   end subroutine writeMatrix
!---------
!> @brief Export as text file the Error Statistics for each DAQuantity.
!! @param [in] DAquantity
!! @param [in] guesfile
!! @param [in] f_number
!! @param [in] winlen
! @param [in]  nobs
!! @param [in] bias
!! @param [in] rmse
!! @param [in] si
!! @param [in] prcntErr array, shape(prcntErr)=(/3,1/)
!! 
!! @details  Format: GuessFileName, F number, Window Length, nobs, Bias, RMSE, SI, Percentile 90, P95, P99
!<
   subroutine write_error_stat(DAquantity,guesfile,f_number,winlen,nobs, bias, rmse, si, prcntErr)
      character(*)   , intent(in) :: DAquantity
      character(*)   , intent(in) :: guesfile
      integer(i_kind), intent(in) :: f_number, nobs
      real(r_kind)   , intent(in) :: winlen,bias,rmse,si
      real(r_kind)   , dimension(:), intent(in) :: prcntErr 
!local
      character(i25), parameter  :: myname='write_error_stat'
      character(i256)            ::fileName
      integer ib_i, ib_j, il_ios
      logical                    :: exst
      integer, parameter :: ip_fid = 41
!
      fileName='ErrorStats_'//trim(DAquantity)//'.txt'      
!
      inquire(file=fileName,exist=exst)
      if (exst) then
         open( unit = ip_fid, file = fileName, status = 'old', position = 'append', action = 'write',form = 'formatted', iostat = il_ios)
      else
         open( unit = ip_fid, file = fileName, status = 'new', action = 'write',form = 'formatted', iostat = il_ios)
      end if
      if (il_ios /= 0) print*,'In',trim(myname),' : Error creating file '//trim(fileName)
!
      write(unit=ip_fid, fmt='(a, I6, f5.1, I6, f7.2, f7.2, f7.2, f7.2, f7.2, f7.2)')  &
      trim(guesfile),f_number, winlen, nobs, bias, rmse, si, prcntErr

      close(ip_fid)
   end subroutine write_error_stat
!
end module general_func
!
