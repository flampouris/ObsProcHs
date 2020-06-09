!> @file qsort_mod.f90
!! @brief The module is an implementation of the quick sort algorithm, taken
!! from Numerical recipes 
!! 
!! @author of the modified version stelios flampouris
!! @version 0.9
!!
!!@copyright GNU (See the main)
!!@date 17-Nov-2016 
!!@todo documentation 
!<
module sorting
   use kinds, only : i_kind,r_kind
   use constants_hs, only :i1
   implicit none
   private
   public qsort
   contains
!
!>@brief Sorts an array arr into ascending order using Quicksort, while making the corresponding
!!rearrangement of the same-size array slave. The sorting and rearrangement
!!are performed by means of an index array.
!!@param[inout] arr - One dimensional array to be sorted and used to sort the
!!slave mattrix
!!@param[in] dir - Character use it ("ascending" or anything with "a") to sort the data in ascending
!!@param[inout] slave - Multidimensional array with the data to be sorted
!!according to the arr. But size(slave,1)=size(arr)
!<
   subroutine qsort(arr,slave)
      real(r_kind), dimension(:), intent(inout) :: arr
      real(r_kind), dimension(:,:), intent(inout), optional:: slave
!  local
      integer(i_kind), dimension(size(arr)) :: ind
      integer(i_kind), dimension(2)         :: sz_slave
      integer(i_kind)                       :: ll
      real(r_kind), dimension(size(arr))    :: dum_arr
!
      call indexx(arr,ind)
      arr=arr(ind)
      if(present(slave))then
      !
      sz_slave=shape(slave)
         do ll=1,sz_slave(1),1
            dum_arr = slave(ll,:) 
            dum_arr = dum_arr(ind)
            slave(ll,:) = dum_arr
         end do
      end if

   end subroutine qsort
!-----

subroutine indexx(arr,index)
   real(r_kind), dimension(:), intent(in) :: arr
   integer(i_kind), dimension(:), intent(out) :: index
!local
   integer(i_kind), parameter :: NN=15, NSTACK=1000
   real(r_kind) :: a
   integer(i_kind) :: n,k,i,j,indext,jstack,l,r
   integer(i_kind), dimension(nstack) :: istack
   real(r_kind), parameter :: r1=1.0
!
   n=size(index)
!   
   index=arth(r1,r1,n)
   jstack=0
   l=1
   r=n
   do
      if (r-l < NN) then
         do j=l+1,r
            indext=index(j)
            a=arr(indext)
            do i=j-1,l,-1
               if (arr(index(i)) <= a) exit
               index(i+1)=index(i)
            end do
            index(i+1)=indext
         end do
         if (jstack == 0) return
         r=istack(jstack)
         l=istack(jstack-1)
         jstack=jstack-2
      else
         k=(l+r)/2
         call swap(index(k),index(l+1))
         call icomp_xchg(index(l),index(r))
         call icomp_xchg(index(l+1),index(r))
         call icomp_xchg(index(l),index(l+1))
         i=l+1
         j=r
         indext=index(l+1)
         a=arr(indext)
         do
            do
               i=i+1
               if (arr(index(i)) >= a) exit
            end do
            do
               j=j-1
               if (arr(index(j)) <= a) exit
            end do
            if (j < i) exit
            call swap(index(i),index(j))
         end do
         index(l+1)=index(j)
         index(j)=indext
         jstack=jstack+2
         if (jstack > NSTACK) stop
         if (r-i+1 >= j-l) then
            istack(jstack)=r
            istack(jstack-1)=i
            r=j-1
         else
            istack(jstack)=j-1
            istack(jstack-1)=l
            l=i
         end if
      end if
   end do
!
contains
   subroutine icomp_xchg(i,j)
      integer(i_kind), intent(inout) :: i,j
!  local
      integer(i_kind) :: swp
!
      if (arr(j) < arr(i)) then
         swp=i
         i=j
         j=swp
      end if
   end subroutine icomp_xchg

end subroutine indexx
!!!!!

!@brief  Array function returning an arithmetic progression.
   function arth(first,increment,n)
      real(r_kind), intent(in) :: first,increment
      integer(i_kind), intent(in) :: n
!  local
      real(r_kind), dimension(n) :: arth
      integer(i_kind) :: k,k2
      real(r_kind) :: temp
      real(r_kind), parameter :: NPAR_ARTH=8
      real(r_kind), parameter :: NPAR2_ARTH=16
!
      if (n > 0) arth(1)=first
      if (n <= NPAR_ARTH) then
         do k=2,n
            arth(k)=arth(k-1)+increment
         end do
      else
         do k=2,NPAR2_ARTH
            arth(k)=arth(k-1)+increment
         end do
         temp=increment*NPAR2_ARTH
         k=NPAR2_ARTH
         do
            if (k >= n) exit
            k2=k+k
            arth(k+1:min(k2,n))=temp+arth(1:min(k,n-k))
            temp=temp+temp
            k=k2
         end do
      end if
   end function arth
!----
!@brief swaps a and b
   subroutinE swap(a,b)
     integer(i_kind), intent(inout) :: a,b
!  local
      integer(i_kind) :: dum
      dum=a
      a=b
      b=dum
   end subroutine swap
!
end module sorting
