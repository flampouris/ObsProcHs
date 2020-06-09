!> @file variables_obsop_hs.f90
!! @brief Declaration of global variables for obsop_hs
!! @details Defining the global variables 
!! @author stelios  
!! @date October 2016
!! @date 18-Nov-2016  by stelios. The default DAquantity is HTSGW 
!! @bug N/A (yet)
!! @version 0.1
!! @copyright GNU
!! @todo update according to the needs of main.
!! @todo collocated%mbr_value allocatable. Now it is i2
!<
module variables_obsop_hs
   use kinds,only : r_kind, r_single, i_kind
   use constants_hs, only : i256, i25, i1, i0, i4, i2, i3, r1

   implicit none
   public :: process_command_line, allocate_member,allocate_alt      &
           , destroy_arrays, allocate_ins, destroy_obs, allocate_all &
           , allocate_collo, obs_data, out_data
   private
   !
!-------------------------------------------------------------------------
!  Variables from command line
!-------------------------------------------------------------------------
!
   character(i256), public :: guesfile='gues.grb2'           !IN (default) 
   character(i256), dimension(:), allocatable, public :: DAquantity !IN 2D Quantities to be assimilated
   integer(i_kind), public :: nv2d = i1 !how many 2D quantities to be assimilated (note : nv2d==nDAq) 
   integer(i_kind), public :: f_number = i0
   real(r_kind)   , public :: winlen = r1      !window length for syncronization in hours
   integer(i_kind), public ::  nsatflnm = i1! number of files with satellite data   
   character(i256), dimension(:), allocatable, public :: obssatinfile!='obssatin.prepbufr'      !IN (default satellite)
   character(i256), public :: obsInSinfile='obsinsin.prepbufr'      !IN (default insitu)
   character(i256), public :: obsoutfile='obsout'   !OUT(default)
   logical        , public :: do_superobs = .false.
   logical        , public :: dotextout = .false.
   logical        , public :: dofilterS = .false.
   logical        , public :: do_despiking = .false.
   logical        , public :: do_ma = .false.
   real(r_kind)   , dimension(i2) , public ::SprObsDim = 0.5_r_kind !delta_x, delta_y for superobs
   real(r_kind)   , public :: despikeStdFac = 2.0_r_kind
   integer(i_kind), public :: MaLngth = 5_i_kind
   real(r_kind)   , public :: obserr_scaling=1_r_kind      !STEVE: use this to scale the input observations
   real(r_kind)   , public :: obs_randselect=1_r_kind      !STEVE: use this to select random input observations
   real(r_kind)   , public :: min_oerr = 1_r_kind
   integer(i_kind), public :: min_quality_level=4  !STEVE: (default) for AVHRR
   logical        , public :: DO_REMOVE_65N = .false.
   logical        , public :: dodebug1 = .false.
!
!-------------------------------------------------------------------------!  
!  Prediction Members
!-------------------------------------------------------------------------!
   integer(i_kind), public :: nlat, nlon, nobs
   real(r_kind) , dimension(:,:),    allocatable, public :: latg, long     !  latg and long of model grid
   real(r_kind) , dimension(i4)                 , public :: mbr_domain     !  [minval(latg), maxval(latg), minval(long), maxval(long)]
   real(r_kind) , dimension(:,:,:),  allocatable, public :: mbr_data       !  data in lat and lon for the quantities of DA
   real(r_kind)                                 , public :: gstime         !  analysis time in minutes from reference date
!  
!--------------------------------------------------------------------------
!  Satellite Data
!--------------------------------------------------------------------------
   !> @brief LETKF has some kind of internal format for the obesrvations defined
   !! by Miyoshi.
   !<
   !> @note
   !! The format :  
   !!1 = obelm
   !!2 = lon
   !!3 = lat
   !!4 = lev
   !!5 = value
   !!6 = oberr
   !! 
   !! e.g. for the type wk 
   !! do i=1,nobs 
   !!    wk(1) = obs_data(i)%typ
   !!    wk(2) = obs_data(i)%x_grd(1)
   !!    wk(3) = obs_data(i)%x_grd(2)
   !!    wk(4) = obs_data(i)%x_grd(3)
   !!    wk(5) = obs_data(i)%value
   !!    wk(6) = obs_data(i)%oerr
   !!    write(fid) wk
   !! end do
   !!
   !!@todo The inrnal GSI format includes 23 variables, probably is a better
   !!choice for including secondary useful information, e.g. depth.  
   !<
   type :: obs_data
      real(r_kind)   :: x_grd(2)   ! longitude, latitude
      real(r_kind)   :: value      ! actual physical value of the parameter measured at this grid point
      real(r_kind)   :: oerr       ! observation standard error
      real(r_kind)   :: hour       ! time of observation
      integer(i_kind):: qkey       ! Quality key
      character(i25) :: typ        !
      real(r_kind)   :: wspa 
      logical :: kept            ! tells letkf whether this obs is kept for assimilation
   end type obs_data
   type (obs_data), allocatable, dimension(:), public :: all_obs, sat_obs, ins_obs, sat_obs_dum
!
   type :: out_data
      real(r_kind)   :: x_grd(i2)   ! longitude, latitude
      real(r_kind)   :: ind_grd(i2)
      real(r_kind)   :: mbr_value(i2)      ! actual physical value of the model at this grid point
      real(r_kind)   :: obs_value      ! actual physical value of the obs at this grid point
      real(r_kind)   :: obs_oerr       ! observation standard error
      real(r_kind)   :: wspa      ! actual physical value of the obs at this grid point
      real(r_kind)   :: timestamp      ! seconds from 0000 1978
      character(i25) :: typ            !
   end type out_data
   type (out_data), allocatable, dimension(:),public :: collocated
!
!-------------------------------------------------------------------------!  
!  Error Statistics
!-------------------------------------------------------------------------!
   real(r_kind)               , public :: bias, rmse, si
   real(r_kind), dimension(3) , public :: prcntErr
   real(r_kind), dimension(3) , public :: prcntl=(/90, 95, 99/)

!
   contains
!
!> @brief allocate the arrays for all the collocated data  
!! To use: call allocate_collo(n)
!!
!! @author stelios 
!! @date october 2016
!! @version 0.1
!! @todo 
!<
   subroutine allocate_collo(n)
!
      integer(i_kind),intent(in) :: n
! local
      character(i25) :: myname =' allocate_collo '
!if (.not.allocated())allocate( ))
      if (allocated(collocated))then 
         deallocate(collocated)
      end if
      allocate(collocated(n))
      if (.not.allocated(collocated))allocate(collocated(n))
!
   end subroutine allocate_collo
!
!> @brief allocate the arrays for importing all the data  
!! To use: call allocate_all(n)
!!
!! @author stelios 
!! @date october 2016
!! @version 0.1
!! @todo 
!<
   subroutine allocate_all(n)
!
      integer(i_kind),intent(in) :: n
! local
      character(i25) :: myname =' allocate_all '
!if (.not.allocated())allocate( ))   
      if (.not.allocated(all_obs))allocate(all_obs(n))
!
   end subroutine allocate_all

!> @brief allocate the arrays for importing the in_situ data  
!! To use: call allocate_ins(n)
!!
!! @author stelios 
!! @date october 2016
!! @version 0.1
!! @todo 
!<
   subroutine allocate_ins(n)
!
      integer(i_kind),intent(in) :: n
! local
      character(i25) :: myname =' allocate_ins '
!if (.not.allocated())allocate( ))   
      if (.not.allocated(ins_obs))allocate(ins_obs(n))
!
   end subroutine allocate_ins

!> @brief allocate the arrays for importing the altimeter data  
!! To use: call allocate_alt(n)
!!
!! @author stelios 
!! @date october 2016
!! @version 0.1
!! @todo 
!<
   subroutine allocate_alt(n)
!
      integer(i_kind),intent(in) :: n
! local
      character(i25) :: myname =' allocate_alt '
!if (.not.allocated())allocate( ))   
      if (.not.allocated(sat_obs))allocate(sat_obs(n))
!
   end subroutine allocate_alt
!
!> @brief allocate the arrays for importing the ensemble members  
!! To use: call allocate_member(nlon,nlat,nv2d)
!!
!! @author stelios 
!! @date october 2016
!! @version 0.1
!! @todo 
!<

   subroutine allocate_member(nlon_loc,nlat_loc,nv2d_loc)
!
      integer(i_kind),intent(in) :: nlon_loc,nlat_loc,nv2d_loc
! local
      character(i25) :: myname =' allocate_memberobsop_hs '
!
!if (.not.allocated())allocate( ))   
      if (.not.allocated(mbr_data)) allocate(mbr_data(nlon_loc,nlat_loc,nv2d_loc))
      if (.not.allocated(long))allocate(long(nlon_loc,nlat_loc))
      if (.not.allocated(latg))allocate(latg(nlon_loc,nlat_loc))
!
   end subroutine allocate_member
!---
!> @brief deallocate the global allocatable arrays for obs
!! To use: call destroy_obs
!!
!! @author stelios 
!! @date october 2016
!! @version 0.1
!! @todo 
!<
!
   subroutine destroy_obs 
!
!   if (allocated()) deallocate()
      if (allocated(sat_obs)     )  deallocate(sat_obs)
      if (allocated(sat_obs_dum) )  deallocate(sat_obs_dum)      
      if (allocated(ins_obs)     )  deallocate(ins_obs)
!    
   end subroutine destroy_obs
!> @brief deallocate the global allocatable arrays  
!! To use: call destroy_arrays
!!
!! @author stelios 
!! @date october 2016
!! @version 0.1
!! @todo 
!<
   subroutine destroy_arrays 
!
!   if (allocated()) deallocate()
      if (allocated(latg)        )  deallocate(latg   )
      if (allocated(long)        )  deallocate(long   )
      if (allocated(mbr_data)    )  deallocate(mbr_data)
      if (allocated(DAquantity)  )  deallocate(DAquantity)
      if (allocated(obssatinfile))  deallocate(obssatinfile)
      if (allocated(sat_obs)     )  deallocate(sat_obs)
      if (allocated(ins_obs)     )  deallocate(ins_obs)
      if (allocated(collocated)  )  deallocate(collocated)
   end subroutine destroy_arrays
!
!> @brief Loads parameters from the command line 
!! To use: inputs are in the format "-x xxx"
!!
!! @authors Steve Penny, Takemasa Miyoshi 
!! @date Earlier than Sep 2016
!! @author stelios 
!! @date october 2016
!! @date 18-Nov 2016 by stelios. Delete forgoten debugging code at Satellite
!! list 
!! @todo transfer it to namelist, all the parameters for the guess files are
!! hard written at the original code. The LETKF-Hs needs more flexibility
!! than that. 
!<
   subroutine process_command_line ()
! 
   character(4*i256) :: command, arg1,arg2
   character(i25), parameter :: myname='process_command_line'
   integer :: i,j, ierr,a
!
   print*, trim(myname), ' starts! '
   call get_command(command,status=ierr)
   call check_err(ierr,myname)
   if (index(command,'debug .t')>0) then 
       dodebug1=.true.
   end if
!
   do i=1,command_argument_count()
      call get_command_argument(i,arg1)
!
      select case (arg1)
      case('-gues')
            call get_command_argument(i+1,arg2,status=ierr)
            call check_err(ierr,myname) 
            guesfile = arg2
            if (dodebug1) call print_arg(i, myname, arg1, arg2)
      case('-davar')
            nv2d = 0
            call get_command_argument(i+nv2d+1,arg2,status=ierr)
            do while (index(trim(arg2),'-')==0)
               nv2d = nv2d + 1                
               call get_command_argument(i+nv2d,arg2,status=ierr)
               if (len_trim(arg2)==0) exit
            end do
            nv2d = nv2d-1
            if (.not.allocated(DAquantity)) allocate (DAquantity(nv2d))
            do j=i+1,i+nv2d,1
               call get_command_argument(j,arg2,status=ierr)
               DAquantity(j-i) = arg2
               if (dodebug1) call print_arg(j, myname, arg1, arg2)
            end do
       case('-fhours') 
            call get_command_argument(i+1,arg2,status=ierr)
            call check_err(ierr,myname)
            read(arg2,*,iostat=ierr) f_number
            if (dodebug1) call print_arg(i, myname, arg1, arg2)
       case('-obssatin')
            nsatflnm = 0
            call get_command_argument(i+nsatflnm+1,arg2,status=ierr)
            do while (index(trim(arg2),'-')==0)
               nsatflnm = nsatflnm + 1
               call get_command_argument(i+nsatflnm,arg2,status=ierr)
               if (len_trim(arg2)==0) exit
            end do
            nsatflnm = nsatflnm-1
            if (.not.allocated(obssatinfile)) allocate (obssatinfile(nsatflnm))
            do j=i+1,i+nsatflnm,1
               call get_command_argument(j,arg2,status=ierr)
               obssatinfile(j-i) = arg2
               if (dodebug1) call print_arg(j, myname, arg1, arg2) 
            end do
            call get_command_argument(i+1,arg2,status=ierr)
            call check_err(ierr,myname) 
            if (dodebug1) call print_arg(i, myname, arg1, arg2)
       case('-obsinsin')
            call get_command_argument(i+1,arg2,status=ierr)
            call check_err(ierr,myname) 
            obsInSinfile = arg2
            if (dodebug1) call print_arg(i, myname, arg1, arg2)
       case('-obsout')
            call get_command_argument(i+1,arg2,status=ierr)
            call check_err(ierr,myname) 
            obsoutfile = arg2
            if (dodebug1) call print_arg(i, myname, arg1, arg2)
       case('-timewnd')
            call get_command_argument(i+1,arg2,status=ierr)
            call check_err(ierr,myname) 
            read(arg2,*,iostat=ierr) winlen 
            if (dodebug1) call print_arg(i, myname, arg1, arg2)
       case('-rm65N')
            call get_command_argument(i+1,arg2,status=ierr)
            call check_err(ierr,myname) 
            read (arg2,*) DO_REMOVE_65N
            if (dodebug1) call print_arg(i, myname, arg1, arg2)
       case('-superob')
            a = 0
            call get_command_argument(i+a+1,arg2,status=ierr)
            do while (index(trim(arg2),'-')==0)
               a = a + 1                
               call get_command_argument(i+a,arg2,status=ierr)
               if (len_trim(arg2)==0) exit
            end do
            a = a-1
            do j=i+1,i+a,1
               call get_command_argument(j,arg2,status=ierr)
               call check_err(ierr,myname) 
               if (j-i==i1) then
                  read (arg2,*) do_superobs
               else if  (j-i==i2) then
                  read(arg2,*,iostat=ierr) SprObsDim(i1)
                  SprObsDim(i2) = SprObsDim(i1)
               else if  (j-i==i3) then
                   read(arg2,*,iostat=ierr) SprObsDim(i2)
               end if
               if (dodebug1) call print_arg(j, myname, arg1, arg2) 
            end do
       case('-filtersat')
            call get_command_argument(i+1,arg2,status=ierr)
            call check_err(ierr,myname) 
            read (arg2,*)dofilterS
            if (dodebug1) call print_arg(i, myname, arg1, arg2)
       case('-thin')
            call get_command_argument(i+1,arg2,status=ierr)
            call check_err(ierr,myname) 
            read (arg2,*) obs_randselect
            if (dodebug1) call print_arg(i, myname, arg1, arg2)
       case('-scale')
            call get_command_argument(i+1,arg2,status=ierr)
            call check_err(ierr,myname) 
            read (arg2,*) obserr_scaling
            if (dodebug1) call print_arg(i, myname, arg1, arg2)
       case('-minerr')
            call get_command_argument(i+1,arg2,status=ierr)
            call check_err(ierr,myname) 
            read (arg2,*) min_oerr
            if (dodebug1) call print_arg(i, myname, arg1, arg2)
       case('-minqc')
            call get_command_argument(i+1,arg2,status=ierr)
            call check_err(ierr,myname) 
            read (arg2,*) min_quality_level
            if (dodebug1) call print_arg(i, myname, arg1, arg2)
       case('-debug')
            call get_command_argument(i+1,arg2,status=ierr)
            call check_err(ierr,myname) 
            read (arg2,*) dodebug1
            if (dodebug1) call print_arg(i, myname, arg1, arg2)
       case('-textout')
            dotextout=.true.
            call get_command_argument(i+1,arg2,status=ierr)
            call check_err(ierr,myname) 
            read (arg2,*) dotextout
            if (dodebug1) call print_arg(i, myname, arg1, arg2)
       case('-despike')
            call get_command_argument(i+1,arg2,status=ierr)
            call check_err(ierr,myname) 
            read (arg2,*) do_despiking
            if (dodebug1) call print_arg(i, myname, arg1, arg2)
            call get_command_argument(i+2,arg2,status=ierr)
            call check_err(ierr,myname) 
            read (arg2,*)despikeStdFac
            if (dodebug1) call print_arg(i+2, myname, arg1, arg2)
       case('-movavg')
            call get_command_argument(i+1,arg2,status=ierr)
            call check_err(ierr,myname) 
            read (arg2,*) do_ma
            if (dodebug1) call print_arg(i, myname, arg1, arg2)
            call get_command_argument(i+2,arg2,status=ierr)
            call check_err(ierr,myname) 
            read (arg2,*)MaLngth
            if (dodebug1) call print_arg(i+2, myname, arg1, arg2)
!       case default
!            print*,trim(myname)
!            print*, "ERROR: option is not supported: ", trim(arg1)
!            print*, "(with value : ", trim(arg2), " )"
!            stop 1
      end select
   enddo
!
   if (.not.allocated(DAquantity)) then
      nv2d=i1
      allocate (DAquantity(nv2d))
      DAquantity(nv2d)='HTSGW'
   end if  
   print*, trim(myname), ' ends! '
!
end subroutine process_command_line
!----
!> @brief check the errors of the "get_command_argument" family
!! @param [in] ierr integer of the status 
!! @param [in] myname, character with varying length of the calling
!! subroutine name
!! @author stelios 
!! @date october 2016
!!
!! @todo 
!<
!----
   subroutine check_err(ierr, myname)
      integer, intent(in) :: ierr
      character(*),intent(in) :: myname
      if (ierr==-1)then
         print*, myname, ': There is an error with the input arguments'
         stop 2 
      end if
   end subroutine check_err
!----
!> @brief Print of the arguments from the "get_command_argument" family,
!! in case of debugging only
!! @param [in] k number of the input argument 
!! @param [in] myname, character with varying length of the calling
!! @param [in] arg1 the name of the argument format: '-xxx' 
!! @param [in] arg2 the value of the arg1 
!! @author stelios 
!! @date october 2016
!!
!<
   subroutine print_arg(k, myname, arg1, arg2)
      integer, intent(in) :: k
      character(*),intent(in) :: myname, arg1, arg2
      print*, myname,'Argument at position ',k, trim(arg1),' with value '  &
             , trim(arg2), ' at position', k+1
   end subroutine print_arg
!----
end module variables_obsop_hs
