!> @file error_statistics.f90
!! @brief The error_statistics calculates the error statistics of one quatity
!! (e.g.) Model Hs vs Obs Hs
!! Quantities to be calculated
!!-diff = model_data - altim_data
!!- bias = diff.mean()
!!- rmse = (diff**2).mean()**0.5
!!- scatter_index = 100.0 *(((diff**2).mean())**0.5-bias**2)/altim_data.mean()
!!- highest percentiles {90, 95, 99} for:
!! * alt_hs, alt_ws, mod_hs, mod_ws
!!- number of observations collocated, used in all of the above.

!!
!! It includes the following subroutines functions:
!! AmIWrong       : Main \n
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
!! @date 16-Nov-2016 - Creation
!<


module error_statistics
   use kinds, only : i_kind, r_kind
   use constants_hs, only : i25
   use sorting, only : qsort
   use variables_obsop_hs, only : dodebug1

   implicit none
   private
   public AmIWrong
   contains
!
!! @brief The AmIWrong is the main for all the calculations of Error
!! param[in] x - 1D array with the values of quantity x
!! param[in] y - 1D array with the values of quantity y
!! param[out] 
!<
!> @details
!! 1. This routine reads all observations end sends them to the main 
!<

   subroutine AmIWrong(nx,x,y,bias, rmse, si,prcntErr,prcntl )
!
      integer(i_kind)            , intent(in ) :: nx 
      real(r_kind), dimension(nx), intent(in ) :: x,y
      real(r_kind)               , intent(out) :: bias, rmse, si
      real(r_kind), dimension(:) , intent(out) :: prcntErr
      real(r_kind), dimension(:),  intent(in)  :: prcntl!=(/90, 95, 99/)
!local 
      character(i25), parameter :: myname='AmIWrong'
      real(r_kind), dimension(nx) :: diffxy,dum
      real(r_kind)               :: mean_obs
!
      integer(i_kind) :: k,rank
!
      dum=0.0
!
      mean_obs = sum(y)/max(1,nx)
      diffxy = x - y
      bias = sum(diffxy)/max(1,nx)
      rmse = sqrt(sum(diffxy**2)/max(1,nx))
      si = 100.0*sqrt(sum((diffxy-bias)**2)/(max(1,nx)-1))/mean_obs
!
      call qsort(diffxy)
!
      if (diffxy(1)> diffxy(nx)) then
         dum=diffxy(nx:-1:1)
         diffxy=dum
      end if
!
      do k=1,size(prcntErr)
         rank=ceiling(max(1,nx)*prcntl(k)/100)
         prcntErr(k)= diffxy(rank)
      end do
!
      if (dodebug1) then
         print*,trim(myname), 'bias=', bias, ' rmse=', rmse, ' si =', si
         do k=1,size(prcntErr)
            print*,prcntl(k),prcntErr(k)
         end do
      end if
!
   end subroutine AmIWrong
end module error_statistics
