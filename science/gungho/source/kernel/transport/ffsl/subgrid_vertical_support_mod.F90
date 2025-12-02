!------------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!> @brief   Routines for computing vertical Nirvana and PPM reconstructions
!!          of fields, for use in FFSL
!------------------------------------------------------------------------------
module subgrid_vertical_support_mod

use constants_mod,                  only: i_def, r_tran, l_def, EPS_R_TRAN

implicit none

private

! Edge reconstructions
public :: second_order_vertical_edge
public :: third_order_vertical_edge
public :: fourth_order_vertical_edge

! Pre-computed coefficients
public :: third_precompute_height
public :: fourth_precompute_height

contains

  ! ========================================================================== !
  ! EDGE RECONSTRUCTIONS
  ! ========================================================================== !

  !----------------------------------------------------------------------------
  !> @brief Calculates the vertical edge values, taking into account the height
  !!        between layers, using a second-order interpolation.
  !> @details Uses a second-order interpolation to find the vertical cell edge
  !!          values of the field. The vertical grid spacing is used to compute the
  !!          mass, and a quadratic is fit through the cumulative
  !!          mass points. This polynomial is differentiated and evaluated
  !!          at the height of the cell edge, to give the cell edge field value.
  !!
  !> @param[in]   field      Field values of two cells which have the ordering
  !!                         | 1 | 2 |
  !> @param[in]   dz         Height of each layer, with index the same as field
  !> @param[in]   edge_to_do Tells routine which edge to do based on
  !!                         cells       | 1 | 2 |
  !!                         with edges  0   1   2
  !> @param[out]  edge_value The interpolated field value at edge_to_do
  !----------------------------------------------------------------------------
  subroutine second_order_vertical_edge(field, dz, edge_to_do, edge_value)

    implicit none

    ! Arguments
    real(kind=r_tran),   intent(in)  :: field(2)
    real(kind=r_tran),   intent(in)  :: dz(2)
    integer(kind=i_def), intent(in)  :: edge_to_do
    real(kind=r_tran),   intent(out) :: edge_value

    ! Internal Variables
    real(kind=r_tran) :: z(0:2), edge_height
    real(kind=r_tran) :: cmass(0:2)

    ! Get heights of edges starting at 0 for lowest edge in stencil
    z(0) = 0.0_r_tran
    z(1) = z(0) + dz(1)
    z(2) = z(1) + dz(2)

    ! Get edge height to interpolate field
    edge_height = z(edge_to_do)

    ! Get cumulative mass
    cmass(0) = 0.0_r_tran
    cmass(1) = cmass(0) + dz(1)*field(1)
    cmass(2) = cmass(1) + dz(2)*field(2)

    ! Calculate derivative of the quadratic at z = edge_height
    edge_value =   ( 2.0_r_tran*edge_height - z(2) ) / ( z(1) * ( z(1)-z(2) ) ) * cmass(1) &
                 + ( 2.0_r_tran*edge_height - z(1) ) / ( z(2) * ( z(2)-z(1) ) ) * cmass(2)

  end subroutine second_order_vertical_edge

  !----------------------------------------------------------------------------
  !> @brief Calculates the vertical edge values, taking into account the height
  !!        between layers, using a third-order interpolation.
  !> @details Uses a third-order interpolation to find the vertical cell edge
  !!          values of the field. The vertical grid spacing is used to compute
  !!          the mass, and a high-order polynomial is fit through the
  !!          cumulative mass points. This polynomial is differentiated and
  !!          evaluated at the height of the cell edge, to give the edge value.
  !!          Both the top and bottom edges are computed for the cell, as this
  !!          edge is not continuous across cells.
  !!
  !> @param[in]     field      Field values of three cells which have
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
  subroutine third_order_vertical_edge(field, dla_dz_1, dla_dz_2, dla_dz_3,    &
                                       dlb_dz_1, dlb_dz_2, dlb_dz_3, dz,       &
                                       edge_above, edge_below,                 &
                                       log_space, nlayers)

    implicit none

    integer(kind=i_def),       intent(in)    :: nlayers
    real(kind=r_tran), target, intent(in)    :: field(nlayers)
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

    real(kind=r_tran), pointer :: field_ptr(:)
    real(kind=r_tran), target  :: log_field(nlayers)
    real(kind=r_tran)          :: mass(nlayers)
    real(kind=r_tran)          :: cmass(nlayers, 3)

    integer(kind=i_def) :: j, k, b_idx, t_idx

    if (log_space) then
      log_field = LOG(MAX(ABS(field), EPS_R_TRAN))
      field_ptr => log_field
    else
      field_ptr => field
    end if

    ! Compute an effective mass, which is field scaled by layer depth
    mass(:) = field_ptr(:) * dz(:)

    ! Compute cumulative mass --------------------------------------------------
    ! The cumulative mass calculation is different for each layer
    ! Need to shift at the domain bottoms and tops

    ! Bottom layer
    k = 1
    cmass(k, 1) = mass(k)
    do j = 2, 3
      cmass(k, j) = cmass(k, j-1) + mass(k+j-1)
    end do

    ! Internal layers
    b_idx = 2
    t_idx = nlayers - 1
    cmass(b_idx:t_idx, 1) = mass(b_idx-1:t_idx-1)
    do j = 2, 3
      cmass(b_idx:t_idx, j) = cmass(b_idx:t_idx, j-1) + mass(b_idx+j-2:t_idx+j-2)
    end do

    ! Top layer
    k = nlayers
    cmass(k, 1) = mass(k-2)
    do j = 2, 3
      cmass(k, j) = cmass(k, j-1) + mass(k+j-3)
    end do

    ! Compute edge values ------------------------------------------------------
    edge_above(:) = (                                                          &
      cmass(:,1)*dla_dz_1(:) + cmass(:,2)*dla_dz_2(:) + cmass(:,3)*dla_dz_3(:) &
    )
    edge_below(:) = (                                                          &
      cmass(:,1)*dlb_dz_1(:) + cmass(:,2)*dlb_dz_2(:) + cmass(:,3)*dlb_dz_3(:) &
    )

    if (log_space) then
      edge_above = EXP(edge_above)
      edge_below = EXP(edge_below)
    end if

  end subroutine third_order_vertical_edge

  !-----------------------------------------------------------------------------
  !> @brief Calculates the vertical edge values, taking into account the height
  !!        between layers, using a fourth-order interpolation.
  !> @details Uses a fourth-order interpolation to find the vertical cell edge
  !!          values of the field. The vertical grid spacing is used to compute
  !!          the mass, and a high-order polynomial is fitted through the
  !!          cumulative mass points. This polynomial is differentiated and
  !!          evaluate at the height of the cell edge, to give the edge value.
  !!
  !> @param[in]     field      Field values of four cells which have the
  !!                           ordering | 1 | 2 | 3 | 4 |
  !> @param[in]     dl_dz      Four sets of coefficients used in calculating
  !!                           the reconstruction at the edge below
  !> @param[in]     dz         Height of each layer, with same index as field
  !> @param[in,out] edge_value Interpolated field value at cell edges
  !> @param[in]     log_space  Whether to perform interpolation on log(field)
  !> @param[in]     nlayers    Number of layers in mesh
  !-----------------------------------------------------------------------------
  subroutine fourth_order_vertical_edge(field, dl_dz_1, dl_dz_2,               &
                                        dl_dz_3, dl_dz_4, dz, edge_value,      &
                                        log_space, nlayers)

    implicit none

    integer(kind=i_def),       intent(in)    :: nlayers
    real(kind=r_tran), target, intent(in)    :: field(nlayers)
    real(kind=r_tran),         intent(in)    :: dz(nlayers)
    real(kind=r_tran),         intent(in)    :: dl_dz_1(0:nlayers)
    real(kind=r_tran),         intent(in)    :: dl_dz_2(0:nlayers)
    real(kind=r_tran),         intent(in)    :: dl_dz_3(0:nlayers)
    real(kind=r_tran),         intent(in)    :: dl_dz_4(0:nlayers)
    real(kind=r_tran),         intent(inout) :: edge_value(0:nlayers)
    logical(kind=l_def),       intent(in)    :: log_space

    real(kind=r_tran), pointer :: field_ptr(:)
    real(kind=r_tran), target  :: log_field(nlayers)
    real(kind=r_tran)          :: mass(nlayers)
    real(kind=r_tran)          :: cmass(0:nlayers,4)

    integer(kind=i_def) :: j, k, b_idx, t_idx

    if (log_space) then
      log_field = LOG(MAX(ABS(field), EPS_R_TRAN))
      field_ptr => log_field
    else
      field_ptr => field
    end if

    ! Compute an effective mass, which is field scaled by layer depth
    mass(:) = field_ptr(:) * dz(:)

    ! Compute cumulative mass --------------------------------------------------
    ! The cumulative mass calculation is different for each layer
    ! Need to shift at the domain bottoms and tops

    ! Bottom two layers
    k = 0
    cmass(k, 1) = mass(k+1)
    do j = 2, 4
      cmass(k, j) = cmass(k, j-1) + mass(k+j)
    end do

    k = 1
    cmass(1, k) = mass(k)
    do j = 2, 4
      cmass(k, j) = cmass(k, j-1) + mass(k+j-1)
    end do

    ! Internal layers
    b_idx = 2
    t_idx = nlayers - 2
    cmass(b_idx:t_idx, 1) = mass(b_idx-1:t_idx-1)
    do j = 2, 4
      cmass(b_idx:t_idx, j) = cmass(b_idx:t_idx, j-1) + mass(b_idx+j-2:t_idx+j-2)
    end do

    ! Top two layers
    k = nlayers - 1
    cmass(k, 1) = mass(k-2)
    do j = 2, 4
      cmass(k, j) = cmass(k, j-1) + mass(k+j-3)
    end do

    k = nlayers
    cmass(k, 1) = mass(k-3)
    do j = 2, 4
      cmass(k, j) = cmass(k, j-1) + mass(k+j-4)
    end do

    ! Compute edge values ------------------------------------------------------
    edge_value(:) = (                                                          &
      cmass(:,1)*dl_dz_1(:) + cmass(:,2)*dl_dz_2(:)                            &
      + cmass(:,3)*dl_dz_3(:) + cmass(:,4)*dl_dz_4(:)                          &
    )

    if (log_space) then
      edge_value = EXP(edge_value)
    end if

  end subroutine fourth_order_vertical_edge

  ! ========================================================================== !
  ! EDGE PRECOMPUTATIONS
  ! ========================================================================== !

  !-----------------------------------------------------------------------------
  !> @brief Calculates constants needed for fourth-order vertical edge
  !!        interpolation, depending on the height of different layers
  !> @param[in,out] dla_dz     Three sets of constants for performing third-
  !!                           order vertical edge reconstruction, to obtain the
  !!                           upper ("above") edge values
  !> @param[in,out] dlb_dz     Three sets of constants for performing third-
  !!                           order vertical edge reconstruction, to obtain the
  !!                           lower ("below") edge values
  !> @param[in]     dz         Height of each layer, with index the same as rho
  !> @param[in]     cell_to_do Tells routine which cell to do based on
  !!                           cells       | 1 | 2 | 3 |
  !!                           with edges  0   1   2   3
  !-----------------------------------------------------------------------------
  subroutine third_precompute_height(dla_dz_1, dla_dz_2, dla_dz_3,             &
                                     dlb_dz_1, dlb_dz_2, dlb_dz_3,             &
                                     dz, cell_to_do)

    implicit none

    real(kind=r_tran),   intent(in)    :: dz(3)
    integer(kind=i_def), intent(in)    :: cell_to_do
    real(kind=r_tran),   intent(inout) :: dla_dz_1, dla_dz_2, dla_dz_3
    real(kind=r_tran),   intent(inout) :: dlb_dz_1, dlb_dz_2, dlb_dz_3

    real(kind=r_tran) :: z(0:3), dzs(1:3), edge_height, inv_sum_dz

    integer(kind=i_def) :: i

    ! Get scaled dz
    inv_sum_dz = 1.0_r_tran / sum(dz)
    dzs = dz * inv_sum_dz

    ! Get heights of edges starting at 0 for lowest edge in stencil
    z(0) = 0.0_r_tran
    do i = 1, 3
      z(i) = z(i-1) + dzs(i)
    end do

    ! Get edge height to interpolate field to for edge below the cell
    edge_height = z(cell_to_do-1)

    ! Calculate derivative of numerator of polynomial at edge height
    dlb_dz_1 = inv_sum_dz *                                                    &
        ( 3.0_r_tran*edge_height**2                                            &
          - 2.0_r_tran*(z(2)+z(3))*edge_height + z(2)*z(3) )                   &
        / ( (z(1)) * (z(1) - z(2)) * (z(1) - z(3)) )
    dlb_dz_2 = inv_sum_dz *                                                    &
        ( 3.0_r_tran*edge_height**2                                            &
          - 2.0_r_tran*(z(1)+z(3))*edge_height + z(1)*z(3) )                   &
        / ( (z(2)) * (z(2) - z(1)) * (z(2) - z(3)) )
    dlb_dz_3 = inv_sum_dz *                                                    &
        ( 3.0_r_tran*edge_height**2                                            &
          - 2.0_r_tran*(z(1)+z(2))*edge_height + z(1)*z(2) )                   &
        / ( (z(3)) * (z(3) - z(1)) * (z(3) - z(2)) )

    ! Get edge height to interpolate field to for edge above the cell
    edge_height = z(cell_to_do)

    ! Calculate derivative of numerator of polynomial at edge height
    dla_dz_1 = inv_sum_dz *                                                    &
        ( 3.0_r_tran*edge_height**2                                            &
          - 2.0_r_tran*(z(2)+z(3))*edge_height + z(2)*z(3) )                   &
        / ( (z(1)) * (z(1) - z(2)) * (z(1) - z(3)) )
    dla_dz_2 = inv_sum_dz *                                                    &
        ( 3.0_r_tran*edge_height**2                                            &
        - 2.0_r_tran*(z(1)+z(3))*edge_height + z(1)*z(3) )                     &
        / ( (z(2)) * (z(2) - z(1)) * (z(2) - z(3)) )
    dla_dz_3 = inv_sum_dz *                                                    &
        ( 3.0_r_tran*edge_height**2                                            &
          - 2.0_r_tran*(z(1)+z(2))*edge_height + z(1)*z(2) )                   &
        / ( (z(3)) * (z(3) - z(1)) * (z(3) - z(2)) )

  end subroutine third_precompute_height

  !-----------------------------------------------------------------------------
  !> @brief Calculates constants needed for fourth-order vertical edge
  !!        interpolation, depending on the height of different layers
  !> @param[in,out] dl_dz      Four sets of constants for performing fourth-
  !!                           order vertical edge reconstruction
  !> @param[in]     dz         Height of each layer, with index the same as rho
  !> @param[in]     edge_to_do Tells routine which edge to do based on
  !!                           cells       | 1 | 2 | 3 | 4 |
  !!                           with edges  0   1   2   3   4
  !-----------------------------------------------------------------------------
  subroutine fourth_precompute_height(dl_dz_1, dl_dz_2, dl_dz_3, dl_dz_4,      &
                                      dz, edge_to_do)

    implicit none

    real(kind=r_tran),   intent(in)    :: dz(4)
    integer(kind=i_def), intent(in)    :: edge_to_do
    real(kind=r_tran),   intent(inout) :: dl_dz_1, dl_dz_2, dl_dz_3, dl_dz_4

    real(kind=r_tran) :: z(0:4), dzs(1:4), edge_height, inv_sum_dz

    integer(kind=i_def) :: i

    ! Get scaled dz
    inv_sum_dz = 1.0_r_tran / sum(dz)
    dzs = dz * inv_sum_dz

    ! Get heights of edges starting at 0 for lowest edge in stencil
    z(0) = 0.0_r_tran
    do i = 1, 4
      z(i) = z(i-1) + dzs(i)
    end do

    ! Get edge height to interpolate rho to
    edge_height = z(edge_to_do)

    ! Calculate derivative of numerator of polynomial at edge height
    dl_dz_1 = inv_sum_dz * (                                                   &
        4.0_r_tran*edge_height**3 - 3.0_r_tran*(z(2)+z(3)+z(4))*edge_height**2 &
        + 2.0_r_tran*(z(3)*z(4) + z(2)*z(3) + z(2)*z(4))*edge_height           &
        - z(2)*z(3)*z(4)                                                       &
    ) / ( (z(1)) * (z(1) - z(2)) * (z(1) - z(3)) * (z(1) - z(4)) )
    dl_dz_2 = inv_sum_dz * (                                                   &
        4.0_r_tran*edge_height**3 - 3.0_r_tran*(z(1)+z(3)+z(4))*edge_height**2 &
        + 2.0_r_tran*(z(3)*z(4) + z(1)*z(3) + z(1)*z(4))*edge_height           &
        - z(1)*z(3)*z(4)                                                       &
    ) / ( (z(2)) * (z(2) - z(1)) * (z(2) - z(3)) * (z(2) - z(4)) )
    dl_dz_3 = inv_sum_dz * (                                                   &
        4.0_r_tran*edge_height**3 - 3.0_r_tran*(z(1)+z(2)+z(4))*edge_height**2 &
        + 2.0_r_tran*(z(2)*z(4) + z(1)*z(2) + z(1)*z(4))*edge_height           &
        - z(1)*z(2)*z(4)                                                       &
    ) / ( (z(3)) * (z(3) - z(1)) * (z(3) - z(2)) * (z(3) - z(4)) )
    dl_dz_4 = inv_sum_dz * (                                                   &
        4.0_r_tran*edge_height**3 - 3.0_r_tran*(z(1)+z(2)+z(3))*edge_height**2 &
        + 2.0_r_tran*(z(2)*z(3) + z(1)*z(2) + z(1)*z(3))*edge_height           &
        - z(1)*z(2)*z(3)                                                       &
    ) / ( (z(4)) * (z(4) - z(1)) * (z(4) - z(2)) * (z(4) - z(3)) )

  end subroutine fourth_precompute_height

end module subgrid_vertical_support_mod
