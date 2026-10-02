!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-----------------------------------------------------------------------------

!> @brief Damps a wind field wherever its 1D or 3D Lipschitz numbers exceed
!!        given thresholds.
!> @details Takes a wind field and returns a new wind with damped values
!!          wherever its 1D or 3D Lipschitz numbers exceed the given thresholds.
!!          This delivers a targeted damping, only impacting areas with large
!!          Lipschitz numbers. Only the outflow values are modified in a given
!!          cell, which may increase the Lipschitz number for a neighbouring
!!          cell (since the inflow there is reduced).
!!          The thresholds are set to 1 for the 1D Lipschitz number and 0.5
!!          for the 3D Lipschitz number.
module lipschitz_damping_kernel_mod

  use argument_mod,          only : arg_type,                                  &
                                    GH_FIELD, GH_SCALAR,                       &
                                    GH_REAL,                                   &
                                    GH_READ, GH_INC, GH_WRITE,                 &
                                    CELL_COLUMN
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : W2, W3
  use kernel_mod,            only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> @brief The kernel metadata for the Lipschitz damping operator.
  !> @details Contains the metadata needed by the PSy layer for this kernel.
  type, public, extends(kernel_type) :: lipschitz_damping_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                        &
        arg_type(GH_FIELD,  GH_REAL, GH_INC,   W2),                            &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2),                            &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),                            &
        arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W3),                            &
        arg_type(GH_SCALAR, GH_REAL, GH_READ)                                  &
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: lipschitz_damping_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: lipschitz_damping_code

contains

!> @brief Damps a wind field wherever its 1D or 3D Lipschitz numbers exceed
!!        given thresholds.
!> @param[in]     nlayers     Number of layers in the mesh.
!> @param[in,out] u_out       Output wind field, damped
!> @param[in]     u_in        Input wind field
!> @param[in]     detj_at_w3  Cell volume, V, used to form the Lipschitz numbers
!> @param[in,out] breached    1 where a Lipschitz number breached a threshold,
!!                            0 otherwise
!> @param[in]     dt          The model timestep length.
!> @param[in]     ndf_w2      Number of degrees of freedom per cell for W2
!> @param[in]     undf_w2     Total num DoFs in this partition for W2
!> @param[in]     map_w2      Dofmap for W2
!> @param[in]     ndf_w3      Number of degrees of freedom per cell for W3
!> @param[in]     undf_w3     Total num DoFs in this partition for W3
!> @param[in]     map_w3      Dofmap for W3
subroutine lipschitz_damping_code(nlayers,                                     &
                                  u_out, u_in,                                 &
                                  detj_at_w3, breached, dt,                    &
                                  ndf_w2, undf_w2, map_w2,                     &
                                  ndf_w3, undf_w3, map_w3)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in)    :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)

  real(kind=r_def),    intent(inout) :: u_out(undf_w2)
  real(kind=r_def),    intent(in)    :: u_in(undf_w2)
  real(kind=r_def),    intent(in)    :: detj_at_w3(undf_w3)
  real(kind=r_def),    intent(inout) :: breached(undf_w3)
  real(kind=r_def),    intent(in)    :: dt

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

  integer(kind=i_def) :: nl, w3_idx, w_idx, s_idx, e_idx, n_idx, b_idx, t_idx

  real(kind=r_def), dimension(0:nlayers-1) :: vol
  real(kind=r_def), dimension(0:nlayers-1) :: uw, us, ue, un, ub, ut
  real(kind=r_def), dimension(0:nlayers-1) :: Lx, Ly, Lz, L3D
  real(kind=r_def), dimension(0:nlayers-1) :: excess, total
  real(kind=r_def), dimension(0:nlayers-1) :: contrib_e, contrib_w, contrib_s
  real(kind=r_def), dimension(0:nlayers-1) :: contrib_n, contrib_t, contrib_b
  real(kind=r_def), dimension(0:nlayers-1) :: breach

  real(kind=r_def), parameter :: threshold_1d = 1.0_r_def
  real(kind=r_def), parameter :: threshold_3d = 0.5_r_def

  nl = nlayers - 1

  w3_idx = map_w3(1)
  w_idx  = map_w2(1)
  s_idx  = map_w2(2)
  e_idx  = map_w2(3)
  n_idx  = map_w2(4)
  b_idx  = map_w2(5)
  t_idx  = map_w2(6)

  vol(:) = detj_at_w3(w3_idx : w3_idx+nl)

  uw(:) = u_in(w_idx : w_idx+nl)
  us(:) = u_in(s_idx : s_idx+nl)
  ue(:) = u_in(e_idx : e_idx+nl)
  un(:) = u_in(n_idx : n_idx+nl)
  ub(:) = u_in(b_idx : b_idx+nl)
  ut(:) = u_in(t_idx : t_idx+nl)

  ! 1D Lipschitz numbers
  Lx(:) = (ue(:) - uw(:))*dt/vol(:)
  Ly(:) = (us(:) - un(:))*dt/vol(:)
  Lz(:) = (ut(:) - ub(:))*dt/vol(:)

  ! Wherever a 1D Lipschitz number exceeds threshold_1d, scale back the
  ! outflowing component(s) responsible so that it falls back to
  ! threshold_1d, leaving any inflowing component untouched. Where both
  ! components are outflowing, the correction is split between them in
  ! proportion to their outflow magnitude.
  excess(:) = max(Lx(:) - threshold_1d, 0.0_r_def)*vol(:)/dt
  contrib_e(:) = max(ue(:), 0.0_r_def)
  contrib_w(:) = max(-uw(:), 0.0_r_def)
  total(:) = contrib_e(:) + contrib_w(:)
  where (total(:) > 0.0_r_def)
    ue(:) = ue(:) - (contrib_e(:)/total(:))*excess(:)
    uw(:) = uw(:) + (contrib_w(:)/total(:))*excess(:)
  end where

  excess(:) = max(Ly(:) - threshold_1d, 0.0_r_def)*vol(:)/dt
  contrib_s(:) = max(us(:), 0.0_r_def)
  contrib_n(:) = max(-un(:), 0.0_r_def)
  total(:) = contrib_s(:) + contrib_n(:)
  where (total(:) > 0.0_r_def)
    us(:) = us(:) - (contrib_s(:)/total(:))*excess(:)
    un(:) = un(:) + (contrib_n(:)/total(:))*excess(:)
  end where

  excess(:) = max(Lz(:) - threshold_1d, 0.0_r_def)*vol(:)/dt
  contrib_t(:) = max(ut(:), 0.0_r_def)
  contrib_b(:) = max(-ub(:), 0.0_r_def)
  total(:) = contrib_t(:) + contrib_b(:)
  where (total(:) > 0.0_r_def)
    ut(:) = ut(:) - (contrib_t(:)/total(:))*excess(:)
    ub(:) = ub(:) + (contrib_b(:)/total(:))*excess(:)
  end where

  ! 3D Lipschitz number, recomputed from the (possibly 1D-clamped) winds,
  ! and clamped in the same way as the 1D numbers above
  L3D(:) = (ue(:) - uw(:) + us(:) - un(:) + ut(:) - ub(:))*dt/vol(:)

  ! Flag cells where any threshold was breached, before the 3D clamp below
  breach(:) = merge( 1.0_r_def, 0.0_r_def,                                     &
                     Lx(:) > threshold_1d .or. Ly(:) > threshold_1d .or.       &
                     Lz(:) > threshold_1d .or. L3D(:) > threshold_3d )
  breached(w3_idx : w3_idx+nl) = breach(:)

  excess(:) = max(L3D(:) - threshold_3d, 0.0_r_def)*vol(:)/dt
  contrib_e(:) = max(ue(:), 0.0_r_def)
  contrib_w(:) = max(-uw(:), 0.0_r_def)
  contrib_s(:) = max(us(:), 0.0_r_def)
  contrib_n(:) = max(-un(:), 0.0_r_def)
  contrib_t(:) = max(ut(:), 0.0_r_def)
  contrib_b(:) = max(-ub(:), 0.0_r_def)
  total(:) = contrib_e(:) + contrib_w(:) + contrib_s(:) + contrib_n(:) + &
             contrib_t(:) + contrib_b(:)
  where (total(:) > 0.0_r_def)
    ue(:) = ue(:) - (contrib_e(:)/total(:))*excess(:)
    uw(:) = uw(:) + (contrib_w(:)/total(:))*excess(:)
    us(:) = us(:) - (contrib_s(:)/total(:))*excess(:)
    un(:) = un(:) + (contrib_n(:)/total(:))*excess(:)
    ut(:) = ut(:) - (contrib_t(:)/total(:))*excess(:)
    ub(:) = ub(:) + (contrib_b(:)/total(:))*excess(:)
  end where

  ! Increment the output field based on the new wind components
  u_out(w_idx : w_idx+nl) = u_out(w_idx : w_idx+nl) + uw(:)
  u_out(s_idx : s_idx+nl) = u_out(s_idx : s_idx+nl) + us(:)
  u_out(e_idx : e_idx+nl) = u_out(e_idx : e_idx+nl) + ue(:)
  u_out(n_idx : n_idx+nl) = u_out(n_idx : n_idx+nl) + un(:)
  u_out(b_idx : b_idx+nl) = u_out(b_idx : b_idx+nl) + ub(:)
  u_out(t_idx : t_idx+nl) = u_out(t_idx : t_idx+nl) + ut(:)

end subroutine lipschitz_damping_code

end module lipschitz_damping_kernel_mod
