!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes weights for convolution for spectral nudging.
!> @details This kernel computes weights for performing a convolution that is
!!          used for nudging the spectrum of a field towards the spectrum of
!!          a reference field.
!!          Only implemented for the lowest-order elements
module convolution_2d_kernel_mod

  use argument_mod,              only: arg_type,                               &
                                       GH_FIELD,                               &
                                       GH_REAL, GH_READ, GH_WRITE,             &
                                       CELL_COLUMN, STENCIL, REGION,           &
                                       ANY_DISCONTINUOUS_SPACE_1,              &
                                       ANY_DISCONTINUOUS_SPACE_9

  use constants_mod,             only: r_def, i_def
  use kernel_mod,                only: kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel.
  !> Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: convolution_2d_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                        &
        arg_type(GH_FIELD, GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),   &
        arg_type(GH_FIELD, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,    &
                                                            STENCIL(REGION)),  &
        arg_type(GH_FIELD, GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_9)    &
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: convolution_2d_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: convolution_2d_code

contains

!> @brief Computes weights for performing convolution for spectral nudging.
!> @param[in]     nlayers         Number of layers
!> @param[in,out] filtered_field  The filtered field to compute
!> @param[in]     field_in        Field to be filtered
!> @param[in]     stencil_size    Size of the stencil used in the convolution
!> @param[in]     stencil_map     Map of the stencil used in the convolution
!> @param[in]     weights         Weights for performing the convolution
!> @param[in]     ndf_3d          Number of DoFs per cell in 3D field
!> @param[in]     undf_3d         Number of DoFs in this partition in 3D field
!> @param[in]     map_3d          DoF map for the 3D field
!> @param[in]     ndf_2d          Num of DoFs per cell in weights field
!> @param[in]     undf_2d         Num of DoFs in this partition in weights field
!> @param[in]     map_2d          DoF map for the weights field
subroutine convolution_2d_code(nlayers,                                        &
                               filtered_field,                                 &
                               field_in,                                       &
                               stencil_size,                                   &
                               stencil_map,                                    &
                               weights,                                        &
                               ndf_3d,                                         &
                               undf_3d,                                        &
                               map_3d,                                         &
                               ndf_2d,                                         &
                               undf_2d,                                        &
                               map_2d)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_2d, ndf_3d
  integer(kind=i_def), intent(in)    :: undf_2d, undf_3d
  integer(kind=i_def), intent(in)    :: stencil_size
  integer(kind=i_def), intent(in)    :: map_2d(ndf_2d)
  integer(kind=i_def), intent(in)    :: map_3d(ndf_3d)
  integer(kind=i_def), intent(in)    :: stencil_map(ndf_3d, stencil_size)
  real(kind=r_def),    intent(in)    :: weights(undf_2d)
  real(kind=r_def),    intent(in)    :: field_in(undf_3d)
  real(kind=r_def),    intent(inout) :: filtered_field(undf_3d)

  ! Local variables
  integer(kind=i_def) :: i, idx_w, idx_s, idx_f, nl

  nl = nlayers - 2 + ndf_3d  ! nlayers for Wtheta, nlayers - 1 for W3
  idx_f = map_3d(1)

  do i = 1, stencil_size
    idx_w = map_2d(1) + i - 1  ! Multidata index
    idx_s = stencil_map(1, i)  ! Stencil index

    filtered_field(idx_f:idx_f+nl) = (                                          &
        filtered_field(idx_f:idx_f+nl)                                          &
        + weights(idx_w) * field_in(idx_s:idx_s+nl)                             &
    )
  end do

end subroutine convolution_2d_code

end module convolution_2d_kernel_mod
