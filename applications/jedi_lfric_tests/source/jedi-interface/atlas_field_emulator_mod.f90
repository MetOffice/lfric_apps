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

  implicit none

  private

  type, public :: atlas_field_emulator_type
    private

    !> The 64-bit floating point values of the field
    real( kind=real64 ), allocatable   :: data(:,:)
    !> The name of the field
    character( len=str_def )           :: field_name
    !> Number of vertical points in the external data
    integer( kind=i_def )              :: n_levels
    !> Number of horizontal points in the external data
    integer( kind=i_def )              :: n_horizontal

  contains

    !> Field initialiser.
    procedure, public :: initialise

    !> Get a pointer to the field data
    procedure, public :: get_data

    !> Set field data to zero
    procedure, public :: zero

    !> Get the name of the field
    procedure, public :: get_field_name

    !> Finalizer
    final             :: atlas_field_emulator_destructor

  end type atlas_field_emulator_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> @brief    Initialiser for atlas_field_emulator_type
!>
!> @param [in] n_levels      The number of levels
!> @param [in] n_horizontal  The number of horizontal points
!> @param [in] field_name    The name of the field
subroutine initialise( self, n_levels, n_horizontal, field_name )

  implicit none

  class( atlas_field_emulator_type ), intent(inout) :: self
  integer( kind=i_def ), intent(in)                 :: n_levels
  integer( kind=i_def ), intent(in)                 :: n_horizontal
  character( len=* ), optional, intent(in)          :: field_name

  allocate( self%data(n_levels, n_horizontal) )
  self % field_name = field_name
  self % n_levels = n_levels
  self % n_horizontal = n_horizontal

end subroutine initialise

!> @brief    Get pointer to field data array
!>
!> @return  data_ptr A pointer to the 2D field data array
function get_data(self) result(data_ptr)


  implicit none

  class( atlas_field_emulator_type ), target, intent(inout) :: self
  real( real64 ), pointer                                   :: data_ptr(:,:)

  data_ptr => self % data

end function get_data

!> @brief    set field values to zero
!>
subroutine zero(self)

  implicit none

  class( atlas_field_emulator_type ), intent(inout) :: self

  self % data = 0.0_r_def

end subroutine zero

!> @brief    Returns the name of the field
!>
!> @param [out] field_name The name of the field
function get_field_name( self ) result( field_name )

  implicit none

  class( atlas_field_emulator_type ), intent(in) :: self
  character( len=str_def )                       :: field_name

  field_name = self%field_name

end function get_field_name

!> @brief    Finaliser for atlas_field_emulator_type
!>
subroutine atlas_field_emulator_destructor(self)

  implicit none

  type(atlas_field_emulator_type), intent(inout) :: self

  if ( allocated(self%data) ) deallocate(self%data)
  self % field_name = ""

end subroutine atlas_field_emulator_destructor

end module atlas_field_emulator_mod
