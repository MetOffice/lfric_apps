!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Computes the finite-volume difference in the y direction.
!> @details The flux form semi-Lagrangian (FFSL) scheme updates the field in
!!          the x, y and z directions separately. This code calculates the
!!          difference for the y direction. The scheme is a simple finite
!!          difference of the fluxes at opposite cell edges and is designed to
!!          work only with lowest order W2h and W3 spaces.

module fv_difference_y_kernel_mod

  use argument_mod,       only : arg_type,            &
                                 GH_FIELD, GH_REAL,   &
                                 GH_WRITE, GH_READ,   &
                                 CELL_COLUMN
  use constants_mod,      only : r_tran, i_def
  use fs_continuity_mod,  only : W2h, W3
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: fv_difference_y_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                   &
         arg_type(GH_FIELD,    GH_REAL, GH_WRITE, W3),    & ! difference
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2h)    & ! flux
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: fv_difference_y_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: fv_difference_y_code

contains

  !> @brief Computes the finite-volume difference in the y direction.
  !> @param[in]     nlayers           The number of layers
  !> @param[in,out] difference        The difference or difference values in W3 space
  !> @param[in]     mass_flux         The flux values which are calculated
  !> @param[in]     ndf_w3            Number of degrees of freedom for W3 per cell
  !> @param[in]     undf_w3           Number of unique degrees of freedom for W3
  !> @param[in]     map_w3            Dofmap for W3
  !> @param[in]     ndf_w2h           Number of degrees of freedom for W2h per cell
  !> @param[in]     undf_w2h          Number of unique degrees of freedom for W2h
  !> @param[in]     map_w2h           Dof map for W2h
  subroutine fv_difference_y_code( nlayers,    &
                                   difference, &
                                   mass_flux,  &
                                   ndf_w3,     &
                                   undf_w3,    &
                                   map_w3,     &
                                   ndf_w2h,    &
                                   undf_w2h,   &
                                   map_w2h )

    implicit none

    ! Arguments
    integer(kind=i_def),                      intent(in)    :: nlayers
    integer(kind=i_def),                      intent(in)    :: ndf_w3
    integer(kind=i_def),                      intent(in)    :: undf_w3
    integer(kind=i_def),                      intent(in)    :: ndf_w2h
    integer(kind=i_def),                      intent(in)    :: undf_w2h
    integer(kind=i_def), dimension(ndf_w3),   intent(in)    :: map_w3
    integer(kind=i_def), dimension(ndf_w2h),  intent(in)    :: map_w2h
    real(kind=r_tran),   dimension(undf_w3),  intent(inout) :: difference
    real(kind=r_tran),   dimension(undf_w2h), intent(in)    :: mass_flux

    integer(kind=i_def) :: nl, w3_idx, N_idx, S_idx

    ! This is based on the lowest order W2h dof map
    !
    !    ---4---
    !    |     |
    !    1     3  horizontal
    !    |     |
    !    ---2---

    w3_idx = map_w3(1)
    S_idx  = map_w2h(2)
    N_idx  = map_w2h(4)
    nl = nlayers - 1

    difference(w3_idx : w3_idx+nl) = (                                         &
        mass_flux(S_idx : S_idx+nl) - mass_flux(N_idx : N_idx+nl)              &
    )

  end subroutine fv_difference_y_code

end module fv_difference_y_kernel_mod
