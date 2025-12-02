!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes a mask to indicate halo cells
module halo_mask_kernel_mod

  use argument_mod,          only : arg_type, func_type,            &
                                    GH_FIELD, GH_WRITE, GH_INTEGER, &
                                    ANY_DISCONTINUOUS_SPACE_3,      &
                                    HALO_CELL_COLUMN

  use constants_mod,         only : i_def
  use kernel_mod,            only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: halo_mask_kernel_type
    private
    type(arg_type) :: meta_args(1) = (/                                        &
        arg_type(GH_FIELD, GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3)    &
    /)
    integer :: operates_on = HALO_CELL_COLUMN
  contains
    procedure, nopass :: halo_mask_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: halo_mask_code

contains

!> @brief Sets a mask which indicates halo cells
!> @param[in]     nlayers           Number of layers in mesh
!> @param[in,out] halo_mask         2D field indicating whether this column is
!!                                  a halo column
!> @param[in]     ndf_w3            Num DoFs per cell for W3
!> @param[in]     undf_w3           Num DoFs for this partition for W3
!> @param[in]     map_w3            DoFmap for W3
subroutine halo_mask_code( nlayers,                                            &
                           halo_mask,                                          &
                           ndf_w3, undf_w3, map_w3                             &
                         )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w3
  integer(kind=i_def), intent(in)    :: undf_w3

  integer(kind=i_def), intent(inout) :: halo_mask(undf_w3)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)

  ! Since this kernel is looping over halo cells, simply set the value to 1
  halo_mask(map_w3(1)) = 1

end subroutine halo_mask_code

end module halo_mask_kernel_mod
