!------------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!> @brief   Routines for computing adjoint of vertical Nirvana and PPM
!!          reconstructions of fields, for use in adjoint FFSL
!------------------------------------------------------------------------------
module adj_subgrid_vertical_support_mod

use constants_mod,                  only: i_def, r_tran, l_def, EPS_R_TRAN
use log_mod,                        only: log_event,                           &
                                          log_scratch_space,                   &
                                          LOG_LEVEL_ERROR

implicit none

private

! Edge reconstructions
public :: adj_third_order_vertical_edge


contains

  ! ========================================================================== !
  ! EDGE RECONSTRUCTIONS
  ! ========================================================================== !

  !----------------------------------------------------------------------------
  !> @brief   Calculates the adjoint of the vertical edge values, taking into account
  !!          the height between layers, using a third-order interpolation.
  !! @details Whilst log_space is included to keep the interfaces identical,
  !!          log_space = true will cause an error as this is not a linear operation.
  !!          Since logging may not be allowed to be handled in kernels and kernel like calls
  !!          this might have to be controlled at the algorithm level. For now we keep the error
  !!          throwing here.
  !!
  !> @param[in,out] field      Field values of three cells which have
  !!                           the ordering | 1 | 2 | 3 |
  !> @param[in]     dla_dz     Three sets of coefficients used in calculating
  !!                           the reconstruction at the edge above
  !> @param[in]     dlb_dz     Three sets of coefficients used in calculating
  !!                           the reconstruction at the edge below
  !> @param[in]     dz         Height of each layer, with same index as field
  !> @param[in,out] edge_above Interpolated value at edge above the cell
  !> @param[in,out] edge_below Interpolated value at edge below the cell
  !> @param[in]     log_space  Whether to perform interpolation on log(field)
  !> @param[in]     nlayers    Number of layers in mesh
  !----------------------------------------------------------------------------
  subroutine adj_third_order_vertical_edge(field, dla_dz_1, dla_dz_2, dla_dz_3,    &
                                           dlb_dz_1, dlb_dz_2, dlb_dz_3, dz,       &
                                           edge_above, edge_below,                 &
                                           log_space, nlayers)

    implicit none

    integer(kind=i_def),       intent(in)    :: nlayers
    real(kind=r_tran), target, intent(inout) :: field(nlayers)
    real(kind=r_tran),         intent(in)    :: dz(nlayers)
    real(kind=r_tran),         intent(in)    :: dla_dz_1(nlayers)
    real(kind=r_tran),         intent(in)    :: dla_dz_2(nlayers)
    real(kind=r_tran),         intent(in)    :: dla_dz_3(nlayers)
    real(kind=r_tran),         intent(in)    :: dlb_dz_1(nlayers)
    real(kind=r_tran),         intent(in)    :: dlb_dz_2(nlayers)
    real(kind=r_tran),         intent(in)    :: dlb_dz_3(nlayers)
    real(kind=r_tran),         intent(inout) :: edge_above(nlayers)
    real(kind=r_tran),         intent(inout) :: edge_below(nlayers)
    logical(kind=l_def),       intent(in)    :: log_space

    real(kind=r_tran), dimension(:),   allocatable :: mass
    real(kind=r_tran), dimension(:,:), allocatable :: cmass

    integer(kind=i_def) :: j, k, b_idx, t_idx

    if (log_space) then
      write(log_scratch_space, *) "log_space = T option not supported: result is not linear!"
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

    allocate(mass(1:nlayers), cmass(1:nlayers, 1:3))

    mass = 0.0_r_tran

    ! Compute edge values ------------------------------------------------------
    cmass(:,1) = edge_below(:)*dlb_dz_1(:) + edge_above(:)*dla_dz_1(:)
    cmass(:,2) = edge_below(:)*dlb_dz_2(:) + edge_above(:)*dla_dz_2(:)
    cmass(:,3) = edge_below(:)*dlb_dz_3(:) + edge_above(:)*dla_dz_3(:)
    edge_below(:) = 0.0_r_tran
    edge_above(:) = 0.0_r_tran

    ! Top layer
    k = nlayers
    do j = 3, 2, -1
      cmass(k, j-1) = cmass(k, j-1) + cmass(k, j)
      mass(k+j-3) = mass(k+j-3) + cmass(k, j)
    end do
    mass(k-2) = mass(k-2) + cmass(k, 1)

    ! Internal layers
    b_idx = 2
    t_idx = nlayers - 1
    do j = 3, 2, -1
      cmass(b_idx:t_idx, j-1) = cmass(b_idx:t_idx, j-1) + cmass(b_idx:t_idx, j)
      mass(b_idx+j-2:t_idx+j-2) = mass(b_idx+j-2:t_idx+j-2) + cmass(b_idx:t_idx, j)
    end do
    mass(b_idx-1:t_idx-1) = mass(b_idx-1:t_idx-1) + cmass(b_idx:t_idx, 1)

    ! Bottom layer
    k = 1
    do j = 3, 2, -1
      cmass(k, j-1) = cmass(k, j-1) + cmass(k, j)
      mass(k+j-1) = mass(k+j-1) + cmass(k, j)
    end do
    mass(k) = mass(k) + cmass(k, 1)

    ! Compute an effective mass, which is field scaled by layer depth
    field(:) = field(:) + mass(:) * dz(:)

    deallocate(mass, cmass)

  end subroutine adj_third_order_vertical_edge

end module adj_subgrid_vertical_support_mod
