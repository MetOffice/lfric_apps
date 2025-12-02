!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Computes the finite-volume divergence in the x direction.
!> @details The flux form semi-Lagrangian (FFSL) scheme updates the field in
!!          the x, y and z directions separately. This code calculates the
!!          divergence for the x direction. The scheme is a simple finite
!!          difference of the fluxes at opposite cell edges and is designed to
!!          work only with lowest order W2 and W3 spaces.

module fv_divergence_x_kernel_mod

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
  type, public, extends(kernel_type) :: fv_divergence_x_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                   &
         arg_type(GH_FIELD,    GH_REAL, GH_WRITE, W3),    & ! difference
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2h),   & ! flux
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W3)     & ! detj_w3
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: fv_divergence_x_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: fv_divergence_x_code

contains

  !> @brief Computes the finite-volume divergence in the x direction.
  !> @param[in]     nlayers           The number of layers
  !> @param[in,out] divergence        The divergence or difference values in W3 space
  !> @param[in]     mass_flux         The flux values which are calculated
  !> @param[in]     detj_w3           det(J) at W3 points, essentially cell volume
  !> @param[in]     ndf_w3            Number of degrees of freedom for W3 per cell
  !> @param[in]     undf_w3           Number of unique degrees of freedom for W3
  !> @param[in]     map_w3            Dofmap for W3
  !> @param[in]     ndf_w2h           Number of degrees of freedom for W2h per cell
  !> @param[in]     undf_w2h          Number of unique degrees of freedom for W2h
  !> @param[in]     map_w2h           Dofmap for W2h
  subroutine fv_divergence_x_code( nlayers,    &
                                   divergence, &
                                   mass_flux,  &
                                   detj_w3,    &
                                   ndf_w3,     &
                                   undf_w3,    &
                                   map_w3,     &
                                   ndf_w2h,    &
                                   undf_w2h,   &
                                   map_w2h )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)                         :: nlayers
    integer(kind=i_def), intent(in)                         :: ndf_w3
    integer(kind=i_def), intent(in)                         :: undf_w3
    integer(kind=i_def), dimension(ndf_w3),   intent(in)    :: map_w3
    integer(kind=i_def), intent(in)                         :: ndf_w2h
    integer(kind=i_def), intent(in)                         :: undf_w2h
    integer(kind=i_def), dimension(ndf_w2h),  intent(in)    :: map_w2h
    real(kind=r_tran),   dimension(undf_w3),  intent(inout) :: divergence
    real(kind=r_tran),   dimension(undf_w2h), intent(in)    :: mass_flux
    real(kind=r_tran),   dimension(undf_w3),  intent(in)    :: detj_w3

    integer(kind=i_def) :: nl, w3_idx, E_idx, W_idx

    ! This is based on the lowest order W2h dof map
    !
    !    ---4---
    !    |     |
    !    1     3  horizontal
    !    |     |
    !    ---2---

    w3_idx = map_w3(1)
    W_idx  = map_w2h(1)
    E_idx  = map_w2h(3)
    nl = nlayers - 1

    divergence(w3_idx : w3_idx+nl) = (                                         &
        mass_flux(E_idx : E_idx+nl) - mass_flux(W_idx : W_idx+nl)              &
    ) / detj_w3(w3_idx : w3_idx+nl)

  end subroutine fv_divergence_x_code

end module fv_divergence_x_kernel_mod
