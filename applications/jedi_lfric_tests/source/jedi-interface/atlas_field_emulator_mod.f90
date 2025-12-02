!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a Atlas field emulator class.
!>
!> @details This module defines an Atlas field emulator class that is included
!>          to provide a simple field container for external field data. In
!           JEDI we are using Atlas fields. The data is stored continuous
!>          vertically and unstructured horizontally. The data is a 2D array
!>          where the first index stores columns and the second the horizontal.
!>
module atlas_field_emulator_mod

  use, intrinsic :: iso_fortran_env, only : real64

  use constants_mod,                 only : r_def, i_def, str_def
  use lfric_mpi_mod,                 only : lfric_mpi_type

  implicit none

  private

  type, public :: atlas_field_emulator_type
    private

    !> The 64-bit floating point values of the field
    real( kind=real64 ), allocatable :: data(:,:)
    !> The name of the field
    character( len=str_def )         :: field_name
    !> Number of vertical points in the external data
    integer( kind=i_def )            :: n_levels
    !> Number of horizontal points in the external data
    integer( kind=i_def )            :: n_horizontal
    !> Communicator to use with the fields
    type( lfric_mpi_type )           :: lfric_mpi_obj

  contains

    !> Field initialiser.
    procedure, public :: initialise

    !> "=" assignment operator overload
    procedure, private :: field_copy
    generic,   public  :: assignment(=) => field_copy

    !> Get a pointer to the field data
    procedure, public :: get_data

    !> Set field data to zero
    procedure, public :: zero

    !> Set field data to random values
    procedure, public :: random

    !> Multiply field by some scalar
    procedure, public :: multiply_by

    !> Compute dot_product with a supplied input field
    procedure, public :: dot_product_with

    !> Compute the sum of the squares
    procedure, public :: sum_of_squares

    !> Get the number of points in the field
    procedure, public :: number_of_points

    !> Compute the root mean square across all field values
    procedure, public :: root_mean_square

    !> Compute the maximum across all field values
    procedure, public :: maximum

    !> Compute the minimum across all field values
    procedure, public :: minimum

    !> Get the name of the field
    procedure, public :: get_field_name

    !> Get the mpi_obj for use internally only
    procedure, private :: get_mpi_obj

    !> Finalizer
    final             :: atlas_field_emulator_destructor

  end type atlas_field_emulator_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> @brief Initialiser for atlas_field_emulator_type
!>
!> @param [in] n_levels      The number of levels
!> @param [in] n_horizontal  The number of horizontal points
!> @param [in] field_name    The name of the field
subroutine initialise( self, n_levels, n_horizontal, lfric_mpi_obj, field_name )

  implicit none

  class( atlas_field_emulator_type ), intent(inout) :: self
  integer( kind=i_def ), intent(in)                 :: n_levels
  integer( kind=i_def ), intent(in)                 :: n_horizontal
  type( lfric_mpi_type ), intent(in)                :: lfric_mpi_obj
  character( len=* ), intent(in)                    :: field_name

  allocate( self%data(n_levels, n_horizontal) )
  self % field_name = field_name
  self % n_levels = n_levels
  self % n_horizontal = n_horizontal

  ! Store lfric_mpi_obj
  call self % lfric_mpi_obj % initialise(lfric_mpi_obj % get_comm())

end subroutine initialise

!> @brief Set field values by copying
!>
subroutine field_copy(self, rhs)

  implicit none

  class( atlas_field_emulator_type ), intent(inout) :: self
  class( atlas_field_emulator_type ),    intent(in) :: rhs

  self%data = rhs%data

end subroutine field_copy

!> @brief Get pointer to field data array
!>
!> @return  data_ptr A pointer to the 2D field data array
function get_data(self) result(data_ptr)

  implicit none

  class( atlas_field_emulator_type ), target, intent(inout) :: self
  real( real64 ), pointer                                   :: data_ptr(:,:)

  data_ptr => self%data

end function get_data

!> @brief Set field values to zero
!>
subroutine zero(self)

  implicit none

  class( atlas_field_emulator_type ), intent(inout) :: self

  self%data = 0.0_r_def

end subroutine zero

!> @brief Set field data to random values
!>
subroutine random(self)

  implicit none

  class( atlas_field_emulator_type ), intent(inout) :: self

  call random_number( self%data )

end subroutine random

!> @brief Multiply field (self%data) by some scalar
!>
!> @param [in] scalar Scalar to multiply field by
subroutine multiply_by(self, scalar)

  implicit none

  class( atlas_field_emulator_type ), intent(inout) :: self
  real( kind=real64 ),                intent(in)    :: scalar

  self%data = self%data * scalar

end subroutine multiply_by

!> @brief Compute dot_product with a supplied input field
!>
!> @param [in] rhs     An Atlas emulator field
!> @return dot_product The dot product of self with rhs
function dot_product_with( self, rhs ) result( dot_product )

  implicit none

  class( atlas_field_emulator_type ), intent(in) :: self
  class( atlas_field_emulator_type ), intent(in) :: rhs
  real( kind=real64 )                            :: dot_product

  ! Local
  real( kind=real64 )             :: l_dot_product
  type( lfric_mpi_type ), pointer :: lfric_mpi_obj

  ! Get local sum
  l_dot_product = sum( self%data*rhs%data )

  ! Compute global sum
  lfric_mpi_obj => self%get_mpi_obj()
  call lfric_mpi_obj%global_sum(l_dot_product, dot_product)

end function dot_product_with

!> @brief Compute the sum of the squares
!>
!> @return sum_squares The sum of (field values)**2
function sum_of_squares( self ) result( sum_squares )

  implicit none

  class( atlas_field_emulator_type ), intent(in) :: self
  real( kind=real64 )                            :: sum_squares

  ! Local
  real( kind=real64 )             :: l_sum_squares
  type( lfric_mpi_type ), pointer :: lfric_mpi_obj

  ! Compute local sum here
  l_sum_squares = sum( self%data*self%data )

  ! Compute global sum
  lfric_mpi_obj => self%get_mpi_obj()
  call lfric_mpi_obj%global_sum(l_sum_squares, sum_squares)

end function sum_of_squares

!> @brief Get the number of points in the field
!>
!> @return g_number_of_points The number of points in the field
function number_of_points( self ) result(g_number_of_points)

  implicit none

  class( atlas_field_emulator_type ), intent(in) :: self
  integer( kind=i_def )                          :: g_number_of_points

  ! Local
  integer( kind=i_def )           :: l_number_of_points
  type( lfric_mpi_type ), pointer :: lfric_mpi_obj

  ! Compute local number of points
  l_number_of_points = self%n_levels*self%n_horizontal

  ! Compute global sum
  lfric_mpi_obj => self%get_mpi_obj()
  call lfric_mpi_obj%global_sum(l_number_of_points, g_number_of_points)

end function number_of_points

!> @brief Compute the root mean square
!>
!> @return g_rms The root mean square of points in the field
function root_mean_square( self ) result( g_rms )

  implicit none

  class( atlas_field_emulator_type ), intent(in) :: self
  real( kind=real64 )                            :: g_rms

  g_rms = sqrt( self%sum_of_squares() / real( self%number_of_points(), kind=real64 ) )

end function root_mean_square

!> @brief Compute the field max
!>
!> @return g_max The maximum value of all points in the field
function maximum( self ) result( g_max )

  implicit none

  class( atlas_field_emulator_type ), intent(in) :: self
  real( kind=real64 )                            :: g_max

  ! Local
  real( kind=real64 )             :: l_max
  type( lfric_mpi_type ), pointer :: lfric_mpi_obj

  ! Compute local max
  l_max = maxval(self%data)

  ! Compute global max
  lfric_mpi_obj => self%get_mpi_obj()
  call lfric_mpi_obj%global_max(l_max, g_max)

end function maximum

!> @brief Retuns the field minimum
!>
!> @return g_min The minimum value of all points in the field
function minimum( self ) result( g_min )

  implicit none

  class( atlas_field_emulator_type ), intent(in) :: self
  real( kind=real64 )                            :: g_min

  ! Local
  real( kind=real64 )             :: l_min
  type( lfric_mpi_type ), pointer :: lfric_mpi_obj

  ! Compute local min
  l_min = minval(self%data)

  ! Compute global min
  lfric_mpi_obj => self%get_mpi_obj()
  call lfric_mpi_obj%global_min(l_min, g_min)

end function minimum

!> @brief Returns the name of the field
!>
!> @return field_name The name of the field
function get_field_name( self ) result( field_name )

  implicit none

  class( atlas_field_emulator_type ), intent(in) :: self
  character( len=str_def )                       :: field_name

  field_name = self%field_name

end function get_field_name

!> @brief Get pointer to the internal mpi object
!>
!> @return  lfric_mpi_obj A pointer to the internal lfric_mpi_obj
function get_mpi_obj(self) result(lfric_mpi_obj)

  implicit none

  class( atlas_field_emulator_type ), target, intent(in) :: self
  type( lfric_mpi_type ),                        pointer :: lfric_mpi_obj

  lfric_mpi_obj => self%lfric_mpi_obj

end function get_mpi_obj

!> @brief Finaliser for atlas_field_emulator_type
!>
subroutine atlas_field_emulator_destructor(self)

  implicit none

  type(atlas_field_emulator_type), intent(inout) :: self

  if ( allocated(self%data) ) deallocate(self%data)
  call self % lfric_mpi_obj % finalise()
  self % field_name = ""

end subroutine atlas_field_emulator_destructor

end module atlas_field_emulator_mod
