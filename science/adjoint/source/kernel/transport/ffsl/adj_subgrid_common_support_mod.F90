!------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!> @brief   Common routines for horizontal and vertical subgrid reconstructions
!!          for use in FFSL
!------------------------------------------------------------------------------
module adj_subgrid_common_support_mod

use constants_mod,                  only: i_def, r_tran, l_def, EPS_R_TRAN
use transport_enumerated_types_mod, only: monotone_none
use log_mod,                        only: log_event,                           &
                                          log_scratch_space,                   &
                                          LOG_LEVEL_ERROR
implicit none

private

public :: adj_subgrid_quadratic_recon

contains

  !----------------------------------------------------------------------------
  !> @brief  Returns the adjoint of the horizontal subgrid quadratic reconstruction.
  !!
  !> @param[in,out] reconstruction The quadratic reconstruction
  !> @param[in]     dep            The absolute value of the fractional
  !!                               departure distance for the cell
  !> @param[in,out] field          Field value in the cell
  !> @param[in,out] edge_left      Field value at left edge of cell
  !> @param[in,out] edge_right     Field value at right edge of cell
  !> @param[in]     monotone       Monotone option to ensure no over/undershoots
  !> @param[in]     order          Order of the reconstruction
  !> @param[in]     nlayers        Number of layers in the mesh
  !----------------------------------------------------------------------------
  subroutine adj_subgrid_quadratic_recon(reconstruction, dep, field, edge_left, &
                                         edge_right, monotone, order, nlayers)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: order
    real(kind=r_tran),   intent(inout) :: reconstruction(nlayers)
    real(kind=r_tran),   intent(in)    :: dep(nlayers)
    real(kind=r_tran),   intent(inout) :: field(nlayers,2*order+1)
    real(kind=r_tran),   intent(inout) :: edge_left(nlayers)
    real(kind=r_tran),   intent(inout) :: edge_right(nlayers)
    integer(kind=i_def), intent(in)    :: monotone

    real(kind=r_tran) :: cm(nlayers), cc(nlayers), cp(nlayers)

    integer(kind=i_def) :: i

    i = order + 1

    select case (monotone)
    ! No limiter ---------------------------------------------------------------
    case (monotone_none)
      cp = 1.0_r_tran - 2.0_r_tran*dep + dep**2
      cc = 3.0_r_tran*dep - 2.0_r_tran*dep**2
      cm = -dep + dep**2

    case default
      write(log_scratch_space, *) "monotone option not supported: result is not linear!"
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end select

    ! Apply weights to field and field edge values
    edge_left(:) = edge_left(:) + cm(:)*reconstruction(:)
    field(:,i) = field(:,i) + cc(:)*reconstruction(:)
    edge_right(:) = edge_right(:) + cp(:)*reconstruction(:)
    reconstruction(:) = 0.0_r_tran

  end subroutine adj_subgrid_quadratic_recon

end module adj_subgrid_common_support_mod
