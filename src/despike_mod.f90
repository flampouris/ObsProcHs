!> @file despike_mod.f90
!! @brief The module is an implementation of the despiking algorithm 
!! presented at  Mori, N., T. Suzuki and S. Kakuno (2007) Noise of acoustic
!! Doppler velocimeter data in bubbly flow, Journal of Engineering 
!! Mechanics, American Society of Civil Engineers, Volume 133, Issue 1, pp.1
!! 22-125.  
!! 
!! @author stelios flampouris
!! @version X.X
!!
!! @copyright GNU (See the main)
!! @date 01-Nov-2017 
!! @todo documentation 
!<
module despiking
   use kinds, only : i_kind,r_kind
   use constants_hs, only: i1
   use sorting
   use general_func, only : writeMatrixdp 
   implicit none
   real(r_kind), parameter :: nan666=666.0

   private
   public despike
   contains
!>@brief 

!!@param[inout] 
!!@param[in] 
!!@param[inout]   
!<

   subroutine despike(data_in, ndata, do_despiking,despikeStdFac &
                     , do_ma, MaLngth, dodebug1                   )
      use constants_hs, only :i25
       
   use general_func, only : min2date
      real(r_kind)   , dimension(:,:), intent(inout)  :: data_in
      integer(i_kind)                , intent(out  )  :: ndata
      logical                        , intent(in   )  :: dodebug1
      logical                        , intent(in   )  :: do_despiking
      logical                        , intent(in   )  :: do_ma
      integer(i_kind)                , intent(in   )  :: MaLngth      
      real(r_kind)                   , intent(in   )  :: despikeStdFac
      
!  local variable
      real(r_kind),dimension(size(data_in,2)-i1)   :: difft
      real(r_kind),allocatable, dimension(:)       :: dum
      real(r_kind),dimension(size(data_in,1) ,size(data_in,2))   :: data_out
      integer(i_kind), dimension(size(difft))      :: ind_chunks
      integer(i_kind)                              :: n_chunks,n
      integer(i_kind)                              :: k,start_ind,end_ind
      integer(i_kind)                              :: i,first_ind_chunk
      real(r_kind)                                 :: sdev,mn,minv,maxv
      character (i25), parameter :: myname='despike'
      integer(i_kind),dimension(6) :: dt
!
      data_out=0.
!!
      if (dodebug1) then       
         print*,trim(myname), ' Shape of data_in',shape(data_in)
         print*,trim(myname), ' Shape of difft',shape(difft)
         call writeMatrixdp('Despike_DataIn.txt', data_in)
      end if 
!
      difft=diff1d(data_in(5,:))
      call rounding(difft,4)
      call chunking(difft,ind_chunks,n_chunks)
!
      if (do_despiking) then
         ndata=0
         n=0
   !
         do k=i1,n_chunks
            if (k==1) then 
               start_ind=1
               end_ind  =ind_chunks(k)
            else if (k==n_chunks) then 
               start_ind=ind_chunks(k-1)
               end_ind  = size(data_in,2)
            else
               start_ind=ind_chunks(k-1)+1
               end_ind  =ind_chunks(k  )-1
            end if   
   !       
            sdev=stddev(data_in(3, start_ind:end_ind))
            mn=meanpos(data_in(3, start_ind:end_ind))
 
            if (dodebug1) then 
               print*,ind_chunks(k), difft(ind_chunks(k))
               print*,'#', trim(myname), '#' 
               print*,'#########'
               print*,'Theoulis1',minval(data_in(3, start_ind:end_ind))
   !
               do i=start_ind,end_ind
                  call min2date(data_in(5,i), dt)
                  print*,dt, data_in(3,i)
               end do
               print*, trim(myname), 'The sdev=',sdev, 'for chunk=',k
               print*,'mean=',mn
            end if
   !        
            first_ind_chunk=ndata+1
   !             
            do i=start_ind,end_ind
               minv=mn-despikeStdFac*sdev
               maxv=mn+despikeStdFac*sdev
   
               if (minv<0.0)minv=0.01
   !            print*,i,'min:',minv,'max:',maxv,'val=', data_in(3,i)
               if (((data_in(3,i)>=minv))               &
                     .and.((data_in(3,i)<=maxv))) then
                  ndata=ndata+1
                  data_out(:,ndata)=data_in(:,i)
                  
                  else
                     n=n+1
               end if
            end do
            if (do_ma) then
               if (.not.allocated(dum)) allocate (dum(ndata-first_ind_chunk+1))
               dum=data_out(3,first_ind_chunk:ndata)
               call MovAvg(dum,MaLngth)
               data_out(3,first_ind_chunk:ndata)=dum
               dum=data_out(8,first_ind_chunk:ndata)
               call MovAvg(dum,MaLngth)
               data_out(8,first_ind_chunk:ndata)=dum

               deallocate(dum)
            end if 
   ! EXPERIMENTAL, to use comment out the sdev approach:  
   !     call despike_phase(data_in(3, start_ind:end_ind))
            
            if (dodebug1) then 
               print*,'#', trim(myname), '#' 
               print*,'#########'
               print*, data_in(3, start_ind:end_ind)
               print*,'#########'
   !          print*,'Theoulis2',minval(data_in(3, start_ind:end_ind))
   !          print*,'666:',count(data_in(3, start_ind:end_ind)==nan666)
   !          do i=start_ind,end_ind
   !            call min2date(data_out(5,i), dt)
   !            print*,dt, data_out(3,i)
   !         end do
            end if 
         end do
         data_in=data_out
         if (dodebug1)call writeMatrixdp('Despike_DataOut.txt', data_out)
      end if
!
!
   end subroutine despike
!
! ###
   subroutine  MovAvg(x,ma)
      real(r_kind), dimension(:), intent(inout) :: x
      integer(i_kind)           , intent(in   ) :: ma      
! Local
      integer(i_kind) :: ma0p5,i
      real(r_kind), dimension(size(x,1))      :: y
      integer(i_kind) :: sindx,eindx,winhs


      ma0p5=INT(floor(real(ma)/2))
      
      y=0.0
 
      do i=1,size( x )
         sindx=i-ma0p5
         if (sindx<=0)then
            sindx=1
            eindx=i + abs(i-sindx)
!            print*,'A i=',i,'sindx=',sindx,'eindx=',eindx!,'ma0p5=',ma0p5
         else if ((i+ma0p5)>size(x))then
            eindx=size(x)
            sindx=i-abs(i-eindx)
!            print*,'B i= ',i, 'sindx=',sindx,'eindx=',eindx
         else
            sindx=i-ma0p5
            eindx=i+ma0p5
!            print*,'C i=',i,'sindx=',sindx,'eindx=',eindx!,'ma0p5=',ma0p5
         end if
         if (eindx > size(x))eindx=size(x)
         y(i)=meanpos(x(sindx:eindx))
!         print*,i,x(i),y(i) 
      end do
         !y(1)=meanpos(x(1:1+ma0p5))
         !y(size(x,1))=meanpos(x(size(x,1)-ma0p5:size(x,1)))
         x=y
   end subroutine MovAvg 
! ###   
   function diff1d (arr1d)
      real(r_kind), dimension(:), intent(in) :: arr1d
      real(r_kind), dimension(size(arr1d)-i1):: diff1d
!
      integer  :: i
! 
      do i = i1,size(arr1d)-i1,i1
         diff1d(i)=arr1d(i+i1)-arr1d(i)
      end do
   end function diff1d
!
! ###
   subroutine rounding(arr1d,prcs)
      real(r_kind), dimension(:), intent(inout) :: arr1d
      integer(i_kind), intent(in)               :: prcs
! 
      integer(i_kind), dimension(size(arr1d))   :: arr1d_int
      real(r_kind), parameter :: ten=10.0
!      
      arr1d_int = ANINT(arr1d*ten**prcs)
      arr1d     = arr1d_int / ten**prcs
   end subroutine rounding
!
! ###
   subroutine chunking (arr1d, ind1d, n_chunks)
      real(r_kind)   , dimension(:), intent(in ) :: arr1d
      integer(i_kind), dimension(:), intent(out) :: ind1d
      integer(i_kind),               intent(out) :: n_chunks
      !     
      real(i_kind), dimension(size(arr1d))    :: dum1d
           real(r_kind)    :: mode_arr1d
      integer(i_kind) :: k,dum
!
      n_chunks=0
!
      mode_arr1d=stat_mode(arr1d)
      
      do k=i1,size(arr1d,1)
         if (ANINT(arr1d(k)/mode_arr1d) < 6.0_r_kind ) cycle
         n_chunks=n_chunks+1
         ind1d(n_chunks)=k !Adjusting the +1 because of the diff
      end do

   end subroutine chunking
! 
   real function stat_mode(a, others, otherslen, ok)
      real(r_kind), dimension(:), intent(in) :: a
      logical, optional, intent(out)    :: ok
      integer, dimension(size(a,1)), optional, intent(out) :: others
      integer, optional, intent(out)    :: otherslen
!
    ! ta is a copy of a, we sort ta modifying it, freq
    ! holds the frequencies and idx the index (for ta) so that
    ! the value appearing freq(i)-time is ta(idx(i))
      real(r_kind),dimension(size(a, 1))  :: ta
      integer, dimension(size(a, 1))      :: freq, idx
      real(r_kind)                        :: tm
      integer                             :: rs, i, ml, tf

      if ( present(ok) ) ok = .false.

      select case ( size(a, 1) )
      case (0) 
         return
      case (1)
         if ( present(ok) ) ok = .true.
         stat_mode = a(1)
         return
      case default
         if ( present(ok) ) ok = .true.
         ta = a
         call qsort(ta)
         freq = 1
         idx = 0
         rs = 1 
!       
         do i = 2, size(ta, 1)
            if ( ta(i-1) == ta(i) ) then
               freq(rs) = freq(rs) + 1
            else
               idx(rs) = i-1
               rs = rs + 1
            end if
         end do
         idx(rs) = i-1
!         
         ml = maxloc(freq(1:rs), 1)  ! index of the max value of freq
         tf = freq(ml)               ! the max frequency
         tm = ta(idx(ml))            ! the value with that freq
!
         if ( present(others) ) then
            i = 1
            others(1) = tm
            do
               freq(ml) = 0
               ml = maxloc(freq(1:rs), 1)
               if ( tf == freq(ml) ) then ! the same freq
                  i = i + 1               ! as the max one
                  others(i) = ta(idx(ml))
               else
                  exit
               end if
            end do
                
            if ( present(otherslen) ) then
               otherslen = i
            end if

         end if
      stat_mode = tm
      end select

  end function stat_mode
!
   subroutine despike_phase(fi)
      real(r_kind), dimension(:), intent(inout ) :: fi
!      real(r_kind), dimension(:), intent(out) :: fo
      integer(i_kind), dimension(size(fi,1)) :: ip
      real(r_kind), dimension(size(fi,1)) :: fi_dum,gfi,ggfi
! Local
      integer (i_kind) :: n_iter, n_out, n_loop, sz_fi,n_ip,k,cnt_h1,cnt_h2
      real(r_kind) :: lambda , fi_mean, theta
! Initialize
      n_iter = 10
      n_out  = 999
      sz_fi  = size(fi,1)
      fi_dum = fi
      lambda = sqrt(2*log(real(sz_fi,r_kind)))
!
! Minimization
      n_loop = 1
      fi_mean= 0
!
      print*, 'size(fi)=',size(fi)
      do while (n_loop <= n_iter)
! Step 1 detrend
         fi_mean = fi_mean + meanpos(fi_dum)
         print*,n_loop,'fi_mean=',fi_mean
         fi_dum = fi_dum - meanpos(fi_dum)
! Step 2 gradient
!         gfi(1:size(fi_dum,1))=diff1d(fi_dum)
!         gfi(size(fi_dum,1))=gfi(size(fi_dum,1)-1)
         
         gfi=gradient(fi_dum)
!print*,gfi
!         ggfi(1:size(fi_dum,1))=diff1d(gfi)
!         ggfi(size(fi_dum,1))=ggfi(size(fi_dum,1)-1)

         ggfi=gradient(gfi)
! Step 3 theta
         if (n_loop==1) then
            theta = atan2( sum(fi_dum*ggfi), sum(fi_dum**2) )
         end if
! Step 4 core - ellipsoid
         ip=0
         call despike_core(fi_dum,gfi,ggfi,theta,ip,n_ip)
         print*, n_loop,'n_ip=',n_ip
         print*, 'ip=',ip
         pause
         if (n_ip==0) exit
! Step 5 exclude values
!         print*, 'A1',maxval(fi_dum)
         cnt_h1=count(fi_dum==nan666)

         do k=i1,n_ip
               fi(ip(k))=nan666
               fi_dum(ip(k))=nan666
            !print*,k, fi_dum(k), fi(k)
         end do
         cnt_h2=count(fi_dum==nan666)
!
         if(cnt_h1==cnt_h2) exit

         print*, maxval(fi_dum)
         print*, maxval(fi)
         
         n_loop = n_loop + 1

      end do


   end subroutine 

   subroutine despike_core (xi,yi,zi,theta,ip,m)
      real(r_kind)   , dimension(:), intent(in ) :: xi,yi,zi
      real(r_kind)   ,               intent(in ) :: theta
      integer(i_kind), dimension(:), intent(out) :: ip
      integer(i_kind) :: m
!
      real(r_kind), dimension(size(xi,1)) :: X,Y,Z
      real(r_kind), dimension(3,3) :: R
      real(r_kind), dimension(3) :: coef
      real(r_kind), parameter :: zero_r=0.0
      real(r_kind) :: lambda,a,b,c
      integer(i_kind) :: n,i
      real(r_kind) :: x1, y1, z1, x2, y2, z2, zt,dis
      
      n = size(xi,1)
      lambda = sqrt(2*log(REAL(n)))
      if (theta == zero_r) then 
         X=xi
         Y=yi
         Z=zi
      else
!         R = reshape((/cos(theta), zero_r    , sin(theta)   &
!                     , zero_r    , 1.0_r_kind, zero_r       &
!                     , -sin(theta),zero_r    , cos(theta)/) &
!                     ,(/3,3/)                               )

        R = reshape((/cos(theta), zero_r    , -sin(theta)   &
                     , zero_r    , 1.0_r_kind, zero_r       &
                     , sin(theta),zero_r    , cos(theta)/) &
                     ,(/3,3/)                               )
       
!
         X = xi*R(1,1) + yi*R(1,2) + zi*R(1,3)
         Y = xi*R(2,1) + yi*R(2,2) + zi*R(2,3)
         Z = xi*R(3,1) + yi*R(3,2) + zi*R(3,3)
      end if
!
      a = lambda*stddev(X)
      b = lambda*stddev(Y)
      c = lambda*stddev(Z)
! Main
!PRINT*,'HELLO',a,b,c
!pause
      m = 0
      do i=i1,n
         x1 = X(i)
         y1 = Y(i)
         z1 = Z(i)
!         pause
!  point on the ellipsoid
         x2 = a*b*c*x1/sqrt((a*c*y1)**2+b**2*(c**2*x1**2+a**2*z1**2))
         y2 = a*b*c*y1/sqrt((a*c*y1)**2+b**2*(c**2*x1**2+a**2*z1**2))
         zt = c**2* ( 1 - (x2/a)**2 - (y2/b)**2 )
  
         if (z1 < zero_r) then
            z2 = -sqrt(zt)
         else if (z1 > zero_r) then
            z2 = sqrt(zt)
         else
            z2 = zero_r
         end if
!         print*, x2,x1, y2,y1, z2,z1
!         pause

!  % check outlier from ellipsoid
!print*,i1,'dis1=',(x1**2+y1**2+z1**2), 'dis2=',(x2**2+y2**2+z2**2), 'dis=', (x2**2+y2**2+z2**2) - (x1**2+y1**2+z1**2)
         dis = (x2**2+y2**2+z2**2) - (x1**2+y1**2+z1**2)
!         print*, 'dis=',dis
         if (dis < zero_r) then  
            m = m + 1
            ip(m) = i
!            xp(m) = xi(i)
!            yp(m) = yi(i)
!            zp(m) = zi(i)
         end if
      end do

      coef = (/a,b,c/)
   end subroutine despike_core
   

   function stddev(vals)
      real(r_kind)                           :: stddev
      real(r_kind), dimension(:), intent(in) :: vals
      real(r_kind) :: mean
      integer(i_kind) :: n
 
      n = max(1,size(vals))
      mean = sum(vals)/n
      stddev = sqrt(sum((vals - mean)**2)/n)
  end function stddev

!   function stddev1(array)
!      real(r_kind)                           :: stddev1
!      real(r_kind), dimension(:), intent(in) :: array
!      real(r_kind), dimension(size(array))   :: dum
!      real(r_kind) :: mean
!      integer(i_kind) :: n,i
!      dum=nan666
!      n = max(1,count(array/=nan666))
!      mean = sum(array,array/=nan666)/n
!      do i=i1,n
!         if (array(i)==nan666)cycle
!         dum(i) = array(i) - mean
!      end do
!      stddev = sqrt(sum(dum**2,dum<300_r_kind )/n)
!   print*,'stddev',stddev
!  end function stddev1

   function meanpos (array)
      real(r_kind) :: meanpos
      real(r_kind), intent(in), dimension(:) :: array
      meanpos = sum(array,array/=nan666) / max(1,count(array/=nan666))
   end function meanpos

   function gradient (array)
!      
      real(r_kind), intent(in), dimension(:) :: array
      real(r_kind), dimension(size(array)) :: gradient
!      real(r_kind), parameter :: nan666=666.0
      integer(i_kind) :: lnga, k
      lnga = size(array,1)
!
      gradient       = 0.0
      
      if (((array(i1) < 40. )) .and.((array(2)<40.))) then 
         !print*, array(i1), nan666,'helo'
         gradient(i1) = array(2) - array(i1)
      end if
      
      if (((array(lnga)< 40. )) .and.((array(lnga-i1)<40.))) then       
         gradient(lnga) =  array(lnga) - array(lnga-i1)
      end if

      do k=2,lnga-i1
!         print*, k, array(k+i1), array(k-i1)
         if (((array(k)<40.)) .and. ((array(k-i1)<40.))) then
            gradient(k)=0.5*(array(k+i1)- array(k-i1))
         end if
!          gradient(k)=array(k+i1)- array(k)
      end do
  !pause
   end function gradient

end module despiking
