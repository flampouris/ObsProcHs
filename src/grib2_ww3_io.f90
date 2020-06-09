!> @file grib2_ww3_io.f90 
!! @brief Module grib2_ww3_io contains all the subroutines for I/O wgrib2 files out of WW3.
!! @details The grib2_ww3_io module is based on the grib2 api by Wesley Ebisuzaki
!! It includes the following subroutines functions:
!! read_time_anl : Imports the time of (whatever quantity) analysis
!! setup_read_member : setups all the variables needed for reading the data
!! read_member : reads the file (Main) 
!! filter_data : includes basic filtering criteria 
!! create_filename : 
!! file_existence Check if file exists \n
!! check_iret : Check for errors returning from wrgibapi \n
!! find_msg : Determine the descriptive string of data record to be read \n
!! define_domain : Define the min and max LonLat of domain \n
!! @author stelios flampouris
!! @version 0.42
!! @copyright The same as the wgrib2 api if none, then GNU GPL
!!
!! @date ??-Sep-2016 - Creation
!! @date 02-Nov-2016 - Last Update
!! @date 17-Nov-2016 - Gave access to the subroutines writeMatrix,
!! @date 29-Nov-2016 - Feature/Bug for this application, with the inv file
!! writeMatrixdp. They write 2d arrays, the first with default precision the
!! second with the user defined precision. 
!! @date 17-Nov-2016 - Gave access to the variable dodebug1 from variables_obsop_hs, 
!! for debugging 
!! @todo Document all the parameters In/Out
!! @bug The inv file is not closed DONE 
!<
module grib2_ww3_io
   use kinds,     only: i_kind, r_single, r_kind
   use constants_hs, only: i200, i25
   use general_func, only: writeMatrix, writeMatrixdp
   use variables_obsop_hs, only:dodebug1
!   
   implicit none
   character(i200), parameter  :: inv = 'inv.in' !'@mem:0'
!
   private 
   public inv, setup_read_member, read_member, define_domain
!
   contains
!-----------------------------------------------------------------------
! read_time_anl
!-----------------------------------------------------------------------
!! @brief subroutine read_time_anl
!! @details The read_time_anl determines the analysis time in minutes from reference date
!! @date 2016 October
!! @author stelios flampouris
!! @version 0.1
!! @copyright GNU
!! @todo DONE :) Whatever it was...
!<
   subroutine read_time_anl(lenmsg, msg, f_number, gstime)
      use constants_hs, only: i0, i2, i5, i10, r60, HrDay
      integer(i_kind), intent(in)   :: lenmsg
      character(lenmsg), intent(in) :: msg
      integer(i_kind), intent(in)   :: f_number
      real(r_kind), intent(out)     :: gstime
!  local
      character(i25)       :: myname ='read_time_anl'
      character(i10)       :: date_str
      integer(i_kind)      :: ind_1
      integer(i_kind)      :: time_1dMN(i5), nmind
      
      real(r_kind)         :: tm_fn
!
      time_1dMN = i0
      ind_1=index(msg,'d=')
      date_str=msg(ind_1+i2:ind_1+i2+i10)
!
      time_1dMN(1)= str2int(date_str(1:4 ))
      time_1dMN(2)= str2int(date_str(5:6 ))
      time_1dMN(3)= str2int(date_str(7:8 ))
      time_1dMN(4)= str2int(date_str(9:10))
!
      call w3fs21(time_1dMN, nmind)
      gstime=real(nmind,r_kind)+(real(f_number,r_kind)*r60)
!
   end subroutine read_time_anl
! 
!-----------------------------------------------------------------------
! Setup read member
!-----------------------------------------------------------------------
!! @brief subroutine setup_read_member
!! @details provides, nlat, nlon
!! @date 2016 October
!! @author stelios flampouris
!! @version 0.1
!! @copyright The same as the wgrib2 api if none, then GNU 
!! @todo export for all the inputs the dimensions not only the analysis
!<
   subroutine setup_read_member(filename,lenDAQ, DAQuantity_1,f_number, nlon,nlat)
      use constants_hs, only : i0, i1
      use wgrib2api
      implicit none
      character(*),intent(in) :: filename
      integer, intent(in):: lenDAQ
      character(lenDAQ), intent(in) :: DAQuantity_1
      integer(i_kind), intent(in) :: f_number
      integer(i_kind), intent(out) :: nlon,nlat
! local
      character(i200)           :: msg
      character(i25) :: myname ='setup_read_member'
!
      call file_existence  (len(filename), filename)
!  
      call check_iret (grb2_mk_inv(filename, inv)+i1)
      call system('ls -lh inv.in')
!      pause
!
      
      call find_msg (lenDAQ, DAQuantity_1, f_number, len(msg), msg)
!
      call check_iret(grb2_inq(filename, inv, msg                       &
                        , nx=nlon, ny=nlat                              & 
                        , close=inv))
!
   end subroutine setup_read_member
!!-----------------------------------------------------------------------
!! Read ensemble data
!!-----------------------------------------------------------------------
!> @brief Subroutine read_member reads the data
!! @details The subroutine read_member imports the fields requested by the user
!! @param[in] filename The name of the grib2 file
!! @param[in] nv2d The number of 2d fields to be read in
!! @param[in] lenDAQ The length of the string of DAQuantity
!! @param[in] DAQuantity The string with the names of quantities to be read
!! @param[in] f_number Time (hours) of prediction, it has to be integer (see variables_obsop_hs)
!! @param[in] nlon Number of grid points in longitude
!! @param[in] nlat Number of points in latitude
!! @param[in] latg The latitude in grid
!! @param[in] long The longitude in grid
!! @param[inout] data_mbr The multidimensional array with the data
!! shape(data_mbr) = [nlon,nlat,nv2d]
!! @param[out] gstime The time from 1.1.1978-00:00
!!
!! @date 2016 October
!! @date 20-Nov-2016 - stelios : dodebug1 is added 
!! @author stelios flampouris
!! @version 0.1
!! @copyright The same as the wgrib2 api if none, then GNU 
!! @todo export for all the inputs the dimensions not only the analysis
!<
   subroutine read_member(filename, nv2d, lenDAQ, DAQuantity               &
                         ,f_number, nlon, nlat                             &
                         ,latg, long, data_mbr, gstime)
!
      use kinds, only : i_kind, r_single 
      use constants_hs, only : i1, uid_rd
      use wgrib2api
!
      character(*),intent(in) :: filename
      character(lenDAQ),dimension(nv2d), intent(in) :: DAQuantity
      integer(i_kind), intent(in) :: lenDAQ, nv2d, f_number, nlon, nlat
      real(r_kind), dimension(nlon,nlat,nv2d),intent(out)  :: data_mbr
      real(r_kind), dimension(nlon,nlat),intent(out) :: latg, long
      real(r_kind), intent(out) :: gstime
!  local variables
      integer(i_kind) :: ll
      real(r_single) , dimension(:,:),  allocatable :: lcl_hsg, lcl_latg,lcl_long
      character(i200)           :: msg
      character (i25),parameter :: myname='read_member'
!
      do ll=1,nv2d,1
!
         call check_iret (grb2_mk_inv(filename, inv)+i1)
!         
         call find_msg (len(DAQuantity), DAQuantity(ll), f_number, len(msg),msg)
!
         call check_iret(grb2_inq(filename, inv, msg                      &
                        , data2=lcl_hsg, lon=lcl_long, lat=lcl_latg       &
                        , close=inv))
!
         call filter_data(lcl_hsg,nlon,nlat,len(DAQuantity),DAQuantity(ll))
!
         data_mbr(:,:,ll)=real(lcl_hsg(:,:),r_kind)
!         
         if (dodebug1) then 
            call writeMatrixdp('model_out_'//trim(DAQuantity(ll))//'.txt',data_mbr(:,:,ll)) 
         end if
      end do 
!
      call read_time_anl(len(msg),msg,f_number, gstime) 
!
      long=real(lcl_long)
      latg=real(lcl_latg)
      if (dodebug1)then
         call writeMatrixdp('model_out_long.txt',long)
         call writeMatrixdp('model_out_latg.txt',latg)
      end if
!
      if (allocated(lcl_hsg )) deallocate(lcl_hsg )
      if (allocated(lcl_latg)) deallocate(lcl_latg)
      if (allocated(lcl_long)) deallocate(lcl_long)
!
      close (uid_rd)
!
   end subroutine read_member
!----
!> @brief Subroutine filter_data is used to filter data
!! @param[inout] arr2d Two Dimensional Array with the data
!! @param[in   ] nx, ny dimensions of arr2d
!! @param[in   ] lenkwrd is the length of the keyword
!! @param[in   ] kwrd The DA quantity
!! @todo extended it as needed for the other assimilation quantities
!<
   subroutine filter_data(arr2d,nx,ny,lenkwrd,kwrd)
      use constants_hs, only: r0
!
      integer(i_kind)                  , intent(in   ) :: nx,ny,lenkwrd
      real(r_single) , dimension(nx,ny), intent(inout) :: arr2d
      character(lenkwrd)               , intent(in   ) :: kwrd
!  local variables
      real(r_single),parameter :: hs_nan = 9999 !Out of grib2 doc
      character(i25),parameter :: myname = 'filter_data'
! 1.
      where (arr2d/=arr2d) arr2d=r0
! 2.
      if (index(kwrd, 'HTSGW')>0.or.index(kwrd,'WIND')>0) then
         where(arr2d>hs_nan) arr2d=hs_nan
      end if
   end subroutine filter_data
!----
!> @brief Subroutine creates name
!! @param[in] ens_no the number of the ensemble member 
!! @param[in] cyc Prediction cycle
!! @param[in] lenflnm The number of characters of filename
!! @param[out] filename of the file to read in
!<
   subroutine create_filename(ens_no, cyc, lenflnm, filename)
      use constants_hs, only : guesflnm,i25
!    
      integer(i_kind)   , intent(in ) :: ens_no, cyc,lenflnm
      character(lenflnm), intent(out) :: filename
!  local variable
      integer :: ind_1, ind_2, ind_3
!
      ind_1 = index(guesflnm,'.' )
      ind_2 = index(guesflnm,'.t')
      ind_3 = index(guesflnm,'z.')
      
     filename = guesflnm(1    :ind_1-3 )//trim(strlen2(ens_no) )&
              //guesflnm(ind_1:ind_2   )//trim(strlen2(cyc   ) )&
              //guesflnm(ind_3:i25     ) 
!  
   end subroutine create_filename
! ---- 
!> @brief Subroutine checking if a file exist
!! @param[in] lenflnm length of the "filename's" string
!! @param[in] filename name of the file to be cheched
!<
   subroutine file_existence(lenflnm, filename)
      integer(i_kind), intent(in)      :: lenflnm
      character(lenflnm), intent(in)   :: filename
      logical :: ex
      character (i25),parameter :: myname='file_existence'

      inquire(file=trim(filename),exist=ex)
      if (.not.ex) then
         print*, trim(myname),' : Input file ', trim(filename)&
               , 'does not exist! Instant Termination!'
         stop
      end if  
!   
   end subroutine file_existence
! ---- 
!> @brief Subroutine checking the return status from the grb2_mk_inv
!! and grb2_inq.
!! @param[in] iret return status
!<
   subroutine check_iret (iret)
      integer, intent(in) :: iret
      character (i25),parameter ::myname='check_iret'
      if (iret/=1) then
         if (iret<=0) then 
            print*,trim(myname),': No match! Check input to grb2'
         else if (iret>1)then
            print*,trim(myname),'There are ', iret,' messages that matched'
         end if
         print*, 'Instant Termination!'
         stop
      end if
   end subroutine check_iret
! ---- 
!
!> @brief Finds the message to be read, given specific keyword and f-value
!! @param [in] lenkwrd is length of the keyword
!! @param [in] keyword is the keyword of interest, e.g. WVHGT for
!! significant wave height
!! @param [in] f_number timestamp in hours for the prediction for the
!! specific cycle. 0 corresponds to the anaysis.  
!! @param [in] lenmsg is the lenght of the message (msg)
!! @param [out] msg is the option string for reading the grib file. 
!! @details To be written
!! @todo use the function ReadLine
!<
   subroutine find_msg (lenkwrd, keyword,f_number,lenmsg,msg)
      use constants_hs, only : i0, i200, uid_rd
      integer(i_kind)      ,  intent(in)  :: lenkwrd, lenmsg, f_number
      character(lenkwrd)   ,  intent(in)  :: keyword
      character(lenmsg),  intent(out)  :: msg
      integer(i_kind) ::ierr,cnt,ind_1,ind_2,i
      character(i200) :: line
      character(i25)  :: f_str
      character(LEN=:), allocatable :: InLine
      character (i25),parameter ::myname='find_msg'
!       
      if (f_number<i0 ) then 
            print*, 'Negative f is not possible! Instant Termination' 
            stop
      else if (f_number==i0) then 
            f_str = 'anl'
      else if (f_number>i0) then
         f_str = trim(str(f_number))//' hour'
      end if
!
      
      ind_1=i0 
      ind_2=i0
!      call sleep (40)

     call system ('ls -lh inv.in')
!      call sleep (10)
      call file_existence  (len(inv), inv)
      print*,inv
      open(unit=uid_rd,file=inv,status='old',action='read',iostat=ierr)
      if(dodebug1) print*,trim(myname), ' ierr = ', ierr 
!For graveyard START
      !do while (ReadLine(uid_rd, InLine))
      !print*,'line',trim(InLine)
      !end do
!END
      do
         line=''
         read(unit=uid_rd,FMT='(A)',iostat=ierr) line
         if(dodebug1) print*,trim(myname), 'line',trim(line), 'ierr', ierr
         if(ierr/=0)exit
         ind_1 = index(line,trim(keyword))
         if (ind_1>i0) then
            ind_2 = index(line(ind_1+1:i200),trim(f_str))
            if (ind_2>i0) then 
               msg = trim(line)
               ind_1=i0
               ind_2=i0
               exit
            end if
         end if
      end do
      close (unit=uid_rd)
   end subroutine find_msg
! ---- 
!
!> @brief defines the coordinates of the model domain (considet it as
!! regional vs global)
!! @param [in] nlon number of grid points in longitude 
!! @param [in] nlat number of grid points in latitude
!! @param [in] latg gridded latitudes
!! @param [in] long gridded longitudes 
!! @param [in] n size of mbr_domain
!! @param [out] mbr_domain [minval(latg), maxval(latg), minval(long), maxval(long)]
!! @details To be written
!! @todo In case of unstructural grid, the @a will be expanded
!<
   subroutine define_domain(nlon, nlat, latg, long, n, mbr_domain)
      integer(i_kind)                      , intent(in) :: nlon, nlat, n
      real   (r_kind), dimension(nlon,nlat), intent(in) :: latg, long  
      real   (r_kind), dimension(n)        , intent(out):: mbr_domain
!
      mbr_domain = [minval(latg), maxval(latg), minval(long), maxval(long)]
   end subroutine define_domain
!----
!> @brief num2str
!! @param [in] k is a number
!! @param [out] str is the k as string! 
!! @details Makes the code beatifull :) 
!<
   character(len=i25) function str(k)
      integer, intent(in) :: k
      write (str, *) k
      str = adjustl(str)
   end function str
!----
!> @brief num2str with len(str)=2
!! @param [in] k is a number
!! @param [out] strlen2 is the k as string with length 2! 
!! @details Makes the code beatifull :) 
!<
   character(len=i25) function strlen2(k)
   integer, intent(in) :: k
!  local variables
   character (i25),parameter :: myname='strlen2'
   character (i25) :: format_string
   if (k < 10) then
      format_string = "(A5,I1)"
   else if (k<=10.and.k<100) then
      format_string = "(A5,I2)"
   else if (k<=100.and.k<1000) then
      format_string = "(A5,I3)"
   else 
      print*,trim(myname), ': The input is grater than 999'
      stop
   endif
   write (strlen2,format_string) k
!
   end function strlen2
! ---- 
!> @brief str2int
!! @param [in] str is a number in characters
!! @param [out] str2int is the str as integer
!! @details Makes the code beatifull :) 
!<
   integer(i_kind) function str2int(str)
   character(len=*), intent(in) :: str
   read(str,*) str2int
   end function str2int
!
end module grib2_ww3_io
!
! Graveyard
!
!function ReadLine(aunit, InLine, trimmed) result(OK)
!integer, intent(IN) :: aunit
!character(LEN=:), allocatable, optional :: InLine
!logical, intent(in), optional :: trimmed
!integer, parameter :: line_buf_len= 1024*4
!character(LEN=line_buf_len) :: InS
!logical :: OK, set
!integer status, size
!
!OK = .false.
!set = .true.
!print*,'aunit',aunit
!do
!    read (aunit,'(a)',advance='NO',iostat=status, size=size) InS
!    print*,InS
!    OK = .not. IS_IOSTAT_END(status)
!    if (.not. OK) return
!    if (present(InLine)) then
!        if (set) then
!            InLine = InS(1:size)
!            set=.false.
!        else
!            InLine = InLine // InS(1:size)
!        end if
!    end if
!    if (IS_IOSTAT_EOR(status)) exit
!end do
!if (present(trimmed) .and. present(InLine)) then
!    if (trimmed) InLine = trim(adjustl(InLine))
!end if
!
!end function ReadLine

