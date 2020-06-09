!> @file constants.f90
!! @brief Declaring constants for LETKF-WW3 (grib2_ww3_io)
!! @details Defining the constants used through out the program
!! @author Stelios  
!! @date September 2016
!!  
!! @bug N/A (yet)
!! @version 0.1
!! @copyright GNU
!>@todo update according to the needs of main.
!<
module constants_hs
   use kinds, only: i_kind, r_kind
   implicit none
!
   !> @brief Constants for Characters 
   !! 
   !! @param[in]  i_kind   : generic specification kind for default integer
   !! @param[out] zero     : integer equal to 0
   !! @param[out] one      : integer equal to 1
   !! @param[out] i25      : integer equal to 25
   !! @param[out] i200     : integer equal to 200   
   !! @param[out] uid_rd   : defauld unit for reading (20)
   !! @param[out] guesflnm : 
   !<
!  set default as private
   private
! Specific integers
   integer, parameter, public :: i0    =   0_i_kind
   integer, parameter, public :: i1    =   1_i_kind
   integer, parameter, public :: i2    =   2_i_kind
   integer, parameter, public :: i3    =   3_i_kind
   integer, parameter, public :: i4    =   4_i_kind
   integer, parameter, public :: i5    =   5_i_kind   
   integer, parameter, public :: i6    =   6_i_kind
   integer, parameter, public :: i7    =   7_i_kind
   integer, parameter, public :: i8    =   8_i_kind
   integer, parameter, public :: i10   =  10_i_kind
   integer, parameter, public :: i25   =  25_i_kind
   integer, parameter, public :: i200  = 200_i_kind
   integer, parameter, public :: i256  = 256_i_kind
!
   real   , parameter, public :: r0    =   0_r_kind
   real   , parameter, public :: r1    =   1_r_kind
   real   , parameter, public :: r60   =  60_r_kind
   real   , parameter, public :: r90   =  90_r_kind
   real   , parameter, public :: r360  =  360_r_kind
!
   real(r_kind), parameter, public :: tiny_r_kind = tiny(r0)
   real(r_kind), parameter, public :: huge_r_kind = huge(r0)
!
! Reading unit
   integer(8), parameter, public :: uid_rd  = 20  
!
! Declaring Filenames and other characters
  character(i25), public :: guesflnm='gwes00.glo_30m.t00z.grib2'
  
! Declare derived constants
   real(r_kind)   , parameter, public :: r60inv   =  1./r60
   real(r_kind)   , parameter, public :: r10inv   =  1./i10
! Trigonometric Constants
   real(r_kind)  , parameter, public :: pi = acos(-1.)
   real(r_kind)  , parameter, public :: deg2rad = pi/(i2*r90)
! Global Constants
   real(r_kind)  , parameter, public :: rearth = 6.3712e+6_r_kind
   real(r_kind)  , parameter, public :: grav = 9.80665e+0_r_kind
   real(r_kind)  , parameter, public :: HrDay = 24_r_kind
! Available Altimeters
   integer(i_kind), parameter, public :: n_sat=i5
   character(len=i25),dimension(n_sat),parameter, public :: namesat=(/'NC031115','NC031122','NC031123','NC031124','NC031130'/)
!
end module constants_hs

