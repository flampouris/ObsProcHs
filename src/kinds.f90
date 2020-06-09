 module kinds
!> @file kinds.f90
!! @brief Definition of kinds for variable definition used at LETKF-WAVES (grib2_ww3_io)
!! @details Defining the kinds of variables
!! @author Stelios  
!! @date September 2016
!! 
!! @param[out] i_kind : generic specification kind for default integer
!! @param[out] r_single : single precision for real 
!! @param[out] r_double : double precision for real
!! @param[out] r_kind : default precision for real
!! 
!! @bug N/A (yet)
!! @version 0.1
!! @copyright GNU
!! @todo update according to the needs of main.
!<
   implicit none
   private 
!
! Default values 
!> @brief Integer types
   integer, parameter, private :: i_byte  = selected_int_kind(1)     !> byte  integer
   integer, parameter, private :: i_short = selected_int_kind(4)     !> short integer
   integer, parameter, private :: i_long  = selected_int_kind(8)     !> long  integer
   integer, parameter, private :: i_llong = selected_int_kind(16)    !> llong integer
!
   integer, parameter, private :: num_i_kinds = 4
!
   integer, parameter, dimension( num_i_kinds ), private :: integer_types = (/   &
                                         i_byte, i_short, i_long,  i_llong  /)
   integer, parameter, private :: num_integer = 3                    !> Change this value for selecting integer kind
!  
   integer, parameter, public  :: i_kind = integer_types( num_integer )
!
!> @brief Real types
  integer, parameter, public  :: r_single = selected_real_kind(6)  ! single precision
  integer, parameter, public  :: r_double = selected_real_kind(15) ! double precision
  integer, parameter, private :: quad_t   = selected_real_kind(20) ! quad precision
  integer, parameter, private :: r_quad   = max( quad_t, r_double )

! Expected 8-bit byte sizes of the real kinds
  integer, parameter, private :: num_bytes_for_r_single = 4
  integer, parameter, private :: num_bytes_for_r_double = 8
  integer, parameter, private :: num_bytes_for_r_quad   = 16

! Define arrays for default definition
  integer, parameter, private :: num_r_kinds = 3
  integer, parameter, dimension( num_r_kinds ), private :: real_kinds = (/ &
       r_single, r_double, r_quad    /)
  integer, parameter, dimension( num_r_kinds ), private :: real_byte_sizes = (/ &
       num_bytes_for_r_single, num_bytes_for_r_double, &
       num_bytes_for_r_quad    /)
! Default values
	  integer, parameter, private :: default_real = 2  ! 2=double,  !> Change this value for selecting real kind
	  integer, parameter, public  :: r_kind = real_kinds( default_real )

	end module kinds
