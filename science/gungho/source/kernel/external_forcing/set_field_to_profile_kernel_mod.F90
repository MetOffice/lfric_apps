!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Initialise a 3D field using a specified 1D profile
!> @details Populate a 3D field from a vertical profile: the field is
!>          constant within each layer, with the value in each layer
!>          being given by the input vertical profile.
module set_field_to_profile_kernel_mod

use argument_mod,  only: arg_type,                  &
                         GH_FIELD, GH_REAL,         &
                         GH_SCALAR, GH_INTEGER,     &
                         GH_WRITE, GH_READ,         &
                         ANY_DISCONTINUOUS_SPACE_1, &
                         CELL_COLUMN

use constants_mod, only: r_def, i_def
use kernel_mod,    only: kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: set_field_to_profile_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                         &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                             &
         /)
    integer :: operates_on = CELL_COLUMN
contains
    procedure, nopass :: set_field_to_profile_code
end type set_field_to_profile_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: set_field_to_profile_code

contains

!> @brief Kernel to set a field using a vertical profile
!> @param[in]     nlayers Number of layers in the grid
!> @param[in,out] field Field to be set using input vertical profile
!> @param[in]     input_layer Index number of layer to be set
!> @param[in]     layer_value Value of field in input_layer
!> @param[in]     ndf Number of degrees of freedom per cell
!> @param[in]     undf Number of unique degrees of freedom
!> @param[in]     map Dofmap of the cell at base of column
  subroutine set_field_to_profile_code( nlayers, field,           &
                                        input_layer, layer_value, &
                                        ndf, undf, map )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                     :: nlayers, ndf, undf
  integer(kind=i_def), intent(in)                     :: input_layer
  integer(kind=i_def), intent(in),    dimension(ndf)  :: map
  real(kind=r_def),    intent(in)                     :: layer_value
  real(kind=r_def),    intent(inout), dimension(undf) :: field

  field( map(1) + input_layer ) = layer_value

  end subroutine set_field_to_profile_code

end module set_field_to_profile_kernel_mod
