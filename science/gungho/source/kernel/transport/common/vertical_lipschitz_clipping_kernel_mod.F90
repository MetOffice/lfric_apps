!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! Some of the content of this file has been produced with the assistance of
! Met Office Github Copilot Enterprise
!-----------------------------------------------------------------------------

!> @brief Clips the vertical components of a wind field wherever its vertical
!!        Lipschitz number exceeds 1.
!> @details Takes a wind field and returns an increment that clips the W2v
!!          values wherever the 1D vertical Lipschitz number exceeds 1.
!!          Only the outflow values are modified in a given cell, which may
!!          increase the Lipschitz number for a neighbouring cell (since the
!!          inflow there is reduced).
module vertical_lipschitz_clipping_kernel_mod

  use argument_mod,          only : arg_type,                                  &
                                    GH_FIELD, GH_SCALAR,                       &
                                    GH_REAL,                                   &
                                    GH_READ, GH_INC,                           &
                                    CELL_COLUMN
  use constants_mod,         only : r_def, i_def
  use fs_continuity_mod,     only : W2, W3
  use kernel_mod,            only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> @brief The kernel metadata for the vertical Lipschitz clipping operator.
  !> @details Contains the metadata needed by the PSy layer for this kernel.
  type, public, extends(kernel_type) :: vertical_lipschitz_clipping_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                        &
        arg_type(GH_FIELD,  GH_REAL, GH_INC,   W2),                            &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  W2),                            &
        arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),                            &
        arg_type(GH_SCALAR, GH_REAL, GH_READ)                                  &
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: vertical_lipschitz_clipping_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: vertical_lipschitz_clipping_code

contains

!> @brief Clips the vertical components of a wind field wherever its vertical
!!        Lipschitz number exceeds 1.
!> @param[in]     nlayers     Number of layers in the mesh
!> @param[in,out] u_inc       Output increment wind field
!> @param[in]     u_in        Input wind field
!> @param[in]     detj_at_w3  Cell volume, V, used to form the Lipschitz number
!> @param[in]     dt          The model timestep length.
!> @param[in]     ndf_w2      Number of degrees of freedom per cell for W2
!> @param[in]     undf_w2     Total num DoFs in this partition for W2
!> @param[in]     map_w2      Dofmap for W2
!> @param[in]     ndf_w3      Number of degrees of freedom per cell for W3
!> @param[in]     undf_w3     Total num DoFs in this partition for W3
!> @param[in]     map_w3      Dofmap for W3
subroutine vertical_lipschitz_clipping_code(nlayers,                           &
                                            u_inc, u_in,                       &
                                            detj_at_w3, dt,                    &
                                            ndf_w2, undf_w2, map_w2,           &
                                            ndf_w3, undf_w3, map_w3)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)    :: nlayers
  integer(kind=i_def), intent(in)    :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in)    :: ndf_w3, undf_w3
  integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)
  integer(kind=i_def), intent(in)    :: map_w3(ndf_w3)

  real(kind=r_def),    intent(inout) :: u_inc(undf_w2)
  real(kind=r_def),    intent(in)    :: u_in(undf_w2)
  real(kind=r_def),    intent(in)    :: detj_at_w3(undf_w3)
  real(kind=r_def),    intent(in)    :: dt

  ! Lowest order W2 dof map: dofs 5 and 6 are the bottom and top faces
  integer(kind=i_def) :: nl, w3_idx, b_idx, t_idx

  real(kind=r_def), dimension(0:nlayers-1) :: vol
  real(kind=r_def), dimension(0:nlayers-1) :: ub, ut
  real(kind=r_def), dimension(0:nlayers-1) :: Lz
  real(kind=r_def), dimension(0:nlayers-1) :: excess, total
  real(kind=r_def), dimension(0:nlayers-1) :: contrib_t, contrib_b

  real(kind=r_def), parameter :: threshold_1d_z = 1.0_r_def

  nl = nlayers - 1

  w3_idx = map_w3(1)
  b_idx  = map_w2(5)
  t_idx  = map_w2(6)

  vol(:) = detj_at_w3(w3_idx : w3_idx+nl)

  ub(:) = u_in(b_idx : b_idx+nl)
  ut(:) = u_in(t_idx : t_idx+nl)

  ! 1D vertical Lipschitz number
  Lz(:) = (ut(:) - ub(:))*dt/vol(:)

  ! Wherever the Lipschitz number exceeds the threshold, scale back the
  ! outflowing component(s) responsible so that it falls back to the
  ! threshold, leaving any inflowing component untouched. Where both
  ! components are outflowing, the correction is split between them in
  ! proportion to their outflow magnitude.
  excess(:) = max(Lz(:) - threshold_1d_z, 0.0_r_def)*vol(:)/dt
  contrib_t(:) = max(ut(:), 0.0_r_def)
  contrib_b(:) = max(-ub(:), 0.0_r_def)
  total(:) = contrib_t(:) + contrib_b(:)
  where (total(:) > 0.0_r_def)
    ut(:) = ut(:) - (contrib_t(:)/total(:))*excess(:)
    ub(:) = ub(:) + (contrib_b(:)/total(:))*excess(:)
  end where

  ! Increment the output field based on the new wind components
  u_inc(b_idx : b_idx+nl) = u_inc(b_idx : b_idx+nl) + ub(:) - u_in(b_idx : b_idx+nl)
  u_inc(t_idx : t_idx+nl) = u_inc(t_idx : t_idx+nl) + ut(:) - u_in(t_idx : t_idx+nl)

end subroutine vertical_lipschitz_clipping_code

end module vertical_lipschitz_clipping_kernel_mod
