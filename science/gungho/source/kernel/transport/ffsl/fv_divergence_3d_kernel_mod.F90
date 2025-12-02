!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief Computes the finite-volume divergence in the 3D direction.
!> @details This is a lowest-order 3D divergence calculation, computing the
!!          divergence as the different of opposite fluxes and dividing by the
!!          cell volume.

module fv_divergence_3d_kernel_mod

  use argument_mod,       only : arg_type,            &
                                 GH_FIELD, GH_REAL,   &
                                 GH_WRITE, GH_READ,   &
                                 CELL_COLUMN
  use constants_mod,      only : r_tran, i_def
  use fs_continuity_mod,  only : W2, W3
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: fv_divergence_3d_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                   &
         arg_type(GH_FIELD,    GH_REAL, GH_WRITE, W3),    & ! difference
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2),    & ! flux
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W3)     & ! detj_w3
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: fv_divergence_3d_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: fv_divergence_3d_code

contains

  !> @brief Computes the finite-volume divergence in 3D.
  !> @param[in]     nlayers           The number of layers
  !> @param[in,out] divergence        The divergence or difference values in W3 space
  !> @param[in]     mass_flux         The flux values which are calculated
  !> @param[in]     detj_w3           det(J) at W3 points, essentially cell volume
  !> @param[in]     ndf_w3            Number of degrees of freedom for W3 per cell
  !> @param[in]     undf_w3           Number of unique degrees of freedom for W3
  !> @param[in]     map_w3            Dofmap for W3
  !> @param[in]     ndf_w2            Number of degrees of freedom for W2 per cell
  !> @param[in]     undf_w2           Number of unique degrees of freedom for W2
  !> @param[in]     map_w2            Dofmap for W2
  subroutine fv_divergence_3d_code( nlayers,    &
                                    divergence, &
                                    mass_flux,  &
                                    detj_w3,    &
                                    ndf_w3,     &
                                    undf_w3,    &
                                    map_w3,     &
                                    ndf_w2,     &
                                    undf_w2,    &
                                    map_w2 )

    implicit none

    ! Arguments
    integer(kind=i_def),                      intent(in)    :: nlayers
    integer(kind=i_def),                      intent(in)    :: ndf_w3
    integer(kind=i_def),                      intent(in)    :: undf_w3
    integer(kind=i_def),                      intent(in)    :: ndf_w2
    integer(kind=i_def),                      intent(in)    :: undf_w2
    integer(kind=i_def), dimension(ndf_w3),   intent(in)    :: map_w3
    integer(kind=i_def), dimension(ndf_w2),   intent(in)    :: map_w2
    real(kind=r_tran),   dimension(undf_w3),  intent(inout) :: divergence
    real(kind=r_tran),   dimension(undf_w2),  intent(in)    :: mass_flux
    real(kind=r_tran),   dimension(undf_w3),  intent(in)    :: detj_w3

    integer(kind=i_def) :: nl, w3_idx, E_idx, W_idx, N_idx, S_idx, B_idx

    ! This is based on the lowest order W2 dof map
    !
    !    ---4---
    !    |     |
    !    1     3  horizontal
    !    |     |
    !    ---2---
    !
    !    ---6---
    !    |     |
    !    |     |  vertical
    !    |     |
    !    ---5---

    w3_idx = map_w3(1)
    W_idx  = map_w2(1)
    S_idx  = map_w2(2)
    E_idx  = map_w2(3)
    N_idx  = map_w2(4)
    B_idx  = map_w2(5)
    nl = nlayers - 1

    divergence(w3_idx : w3_idx+nl) = (                                         &
        mass_flux(E_idx : E_idx+nl) - mass_flux(W_idx : W_idx+nl)              &
        + mass_flux(S_idx : S_idx+nl) - mass_flux(N_idx : N_idx+nl)            &
        + mass_flux(B_idx+1 : B_idx+nl+1) - mass_flux(B_idx : B_idx+nl)        &
    ) / detj_w3(w3_idx : w3_idx+nl)

  end subroutine fv_divergence_3d_code

end module fv_divergence_3d_kernel_mod
