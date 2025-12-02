!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Extract 2D field from a specified layer of a 3D field
!> @details Kernel to copy one layer of a 3D field into a 2D field.
module extract_2d_field_kernel_mod

use argument_mod,      only: arg_type,                  &
                             GH_FIELD, GH_REAL,         &
                             GH_SCALAR, GH_INTEGER,     &
                             GH_WRITE, GH_READ,         &
                             ANY_DISCONTINUOUS_SPACE_1, &
                             ANY_DISCONTINUOUS_SPACE_2, &
                             CELL_COLUMN

use constants_mod,     only: r_def, i_def
use kernel_mod,        only: kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: extract_2d_field_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                     &
         arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                           &
         /)
    integer :: operates_on = CELL_COLUMN
contains
    procedure, nopass :: extract_2d_field_code
end type extract_2d_field_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: extract_2d_field_code

contains

  !> @brief Extract a single layer from a 3D field
  !> @param[in]     nlayers Number of layers in the 3D field
  !> @param[in,out] twod_field 2D field into which layer of 3D field is copied
  !> @param[in]     field 3D field from which a single layer is copied
  !> @param[in]     extract_layer Index number of layer to extract
  !> @param[in]     ndf_2d Number of degrees of freedom per cell for 2D field
  !> @param[in[     undf_2d Number of unique degrees of freedom for 2D field
  !> @param[in]     map_2d Dofmap for 2D field
  !> @param[in]     ndf_3d Number of degrees of freedom per cell for 3D field
  !> @param[in[     undf_3d Number of unique degrees of freedom for 3D field
  !> @param[in]     map_3d Dofmap for base layer of 3D field
  subroutine extract_2d_field_code( nlayers, twod_field, field, extract_layer, &
                                    ndf_2d, undf_2d, map_2d,                   &
                                    ndf_3d, undf_3d, map_3d )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)                        :: nlayers
  integer(kind=i_def), intent(in)                        :: ndf_3d
  integer(kind=i_def), intent(in)                        :: undf_3d
  integer(kind=i_def), intent(in)                        :: ndf_2d
  integer(kind=i_def), intent(in)                        :: undf_2d
  integer(kind=i_def), intent(in)                        :: extract_layer
  integer(kind=i_def), intent(in),    dimension(ndf_2d)  :: map_2d
  integer(kind=i_def), intent(in),    dimension(ndf_3d)  :: map_3d
  real(kind=r_def),    intent(in),    dimension(undf_3d) :: field
  real(kind=r_def),    intent(inout), dimension(undf_2d) :: twod_field

  twod_field( map_2d(1) ) = field( map_3d(1) + extract_layer )

  end subroutine extract_2d_field_code

end module extract_2d_field_kernel_mod
