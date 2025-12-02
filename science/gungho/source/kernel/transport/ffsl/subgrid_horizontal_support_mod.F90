!------------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!> @brief   Routines for computing horizontal Nirvana and PPM reconstructions
!!          of fields, for use in FFSL
!------------------------------------------------------------------------------
module subgrid_horizontal_support_mod

use constants_mod,                  only: i_def, r_tran, l_def, EPS_R_TRAN

implicit none

private

! Edge interpolation routines
public :: fourth_order_horizontal_edge
public :: nirvana_horizontal_edge
public :: nirvana_special_edge
public :: linear_special_edge
public :: fourth_order_special_edge
public :: fourth_nirvana_special_edge

contains

  !----------------------------------------------------------------------------
  !> @brief  Calculates the interpolated field at the edges of a cell required for
  !!         using horizontal PPM to estimate the quadratic subgrid representation
  !!         of a field. The function is passed five field values from consecutive cells
  !!         (which all lie in the same direction) with the dofmap
  !!         | 1 | 2 | 3 | 4 | 5 | and returns the estimated field value between
  !!         cells 2 and 3 (edge_left) and cells 3 and 4 (edge_right).
  !!         The cells are assumed to be uniform in spacing.
  !!
  !> @param[in]     field         Has dof map of the form | 1 | 2 | 3 | 4 | 5 |
  !> @param[in,out] edge_left     Field value at left edge of cell
  !> @param[in,out] edge_right    Field value at right edge of cell
  !> @param[in]     nlayers       Number of layers in the mesh
  !----------------------------------------------------------------------------
  subroutine fourth_order_horizontal_edge(field, edge_left, edge_right, nlayers)

    implicit none

    integer(kind=i_def), intent(in)  :: nlayers
    real(kind=r_tran),   intent(in)  :: field(nlayers,5)
    real(kind=r_tran),   intent(out) :: edge_left(nlayers), edge_right(nlayers)

    real(kind=r_tran), parameter :: twelfth = 1.0_r_tran/12.0_r_tran

    ! As the cell widths are assumed to be constant the edge value reduce to
    ! that given in Colella and Woodward, JCP, 54, 1984, equation (1.9)
    edge_left(:)  = twelfth * (                                                &
        7.0_r_tran * (field(:,2) + field(:,3)) - (field(:,1) + field(:,4))     &
    )
    edge_right(:) = twelfth * (                                                &
        7.0_r_tran * (field(:,3) + field(:,4)) - (field(:,2) + field(:,5))     &
    )

  end subroutine fourth_order_horizontal_edge

  !----------------------------------------------------------------------------
  !> @brief  Calculates the interpolated field at the edges of a cell required
  !!         for using horizontal Nirvana to estimate the quadratic subgrid
  !!         representation of a field. The function is passed three field
  !!         values from consecutive cells (which all lie in the same direction)
  !!         with the dofmap | 1 | 2 | 3 | and returns the estimated field value
  !!         between cells 1 and 2 (edge_left) and cells 2 and 3 (edge_right).
  !!         The cells are assumed to be uniform in spacing.
  !!
  !> @param[in]     field            Has dof map of the form | 1 | 2 | 3 |
  !> @param[in,out] edge_left        Field value at left edge of cell
  !> @param[in,out] edge_right       Field value at right edge of cell
  !> @param[in]     nlayers          Number of layers in the mesh
  !----------------------------------------------------------------------------
  subroutine nirvana_horizontal_edge(field, edge_left, edge_right, nlayers)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    real(kind=r_tran),   intent(in)    :: field(nlayers,3)
    real(kind=r_tran),   intent(inout) :: edge_left(nlayers), edge_right(nlayers)

    real(kind=r_tran), parameter :: sixth = 1.0_r_tran/6.0_r_tran

    ! The Nirvana edges can be derived from Leonard et al. 1995.
    edge_left(:)  = (-field(:,3) + 5.0_r_tran*field(:,2) + 2.0_r_tran*field(:,1)) * sixth
    edge_right(:) = (-field(:,1) + 5.0_r_tran*field(:,2) + 2.0_r_tran*field(:,3)) * sixth

  end subroutine nirvana_horizontal_edge

  !----------------------------------------------------------------------------
  !> @brief  Calculates the interpolated field at the edges of a cell required for
  !!         using horizontal Nirvana to estimate the quadratic subgrid representation
  !!         of a field. The special edge treatment is used to shift the stencil left
  !!         or right depending on where the panel edge lies. As the interpolation is
  !!         third order, even with special edges, this requires a stencil of
  !!         size 5 instead of the usual 3 for Nirvana. The field dofmap is
  !!         | 1 | 2 | 3 | 4 | 5 | and this returns the estimated field value at
  !!         the edges of cell 3. The cells are assumed to be uniform in spacing.
  !!
  !> @param[in]     field        Has dof map of the form | 1 | 2 | 3 | 4 | 5 |
  !> @param[in]     spt_case     Special edges case to shift stencil
  !> @param[in,out] edge_left    Field value at left edge of cell
  !> @param[in,out] edge_right   Field value at right edge of cell
  !> @param[in]     nlayers      Number of layers in the mesh
  !----------------------------------------------------------------------------
  subroutine nirvana_special_edge(field, spt_case, edge_left, edge_right, nlayers)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    real(kind=r_tran),   intent(in)    :: field(nlayers,5)
    integer(kind=i_def), intent(in)    :: spt_case(nlayers)
    real(kind=r_tran),   intent(inout) :: edge_left(nlayers)
    real(kind=r_tran),   intent(inout) :: edge_right(nlayers)

    real(kind=r_tran) :: s1(nlayers), s2(nlayers), s3(nlayers)
    real(kind=r_tran), parameter :: sixth = 1.0_r_tran/6.0_r_tran

    ! Switches that are 1 for particular spt_case values, and 0 otherwise
    s1(:) = MERGE(1.0_r_tran, 0.0_r_tran, spt_case(:) == 1)
    s2(:) = MERGE(1.0_r_tran, 0.0_r_tran, spt_case(:) == 2)
    s3(:) = MERGE(1.0_r_tran, 0.0_r_tran, spt_case(:) == 3)

    ! The switches turn on the coefficients for the particular spt_case. If
    ! none are turned on, revert to normal quadratic reconstruction
    edge_left(:) = sixth * (                                                   &
      -s3(:) * field(:,1)                                                      &
      + (2.0_r_tran*s2(:) + 5.0_r_tran*s3(:)) * field(:,2)                     &
      + (11.0_r_tran*s1(:) + 5.0_r_tran*s2(:) + 2.0_r_tran*s3(:)) * field(:,3) &
      + (-7.0_r_tran*s1(:) - s2(:)) * field(:,4)                               &
      + 2.0_r_tran*s1(:) * field(:,5)                                          &
    )
    edge_right(:) = sixth * (                                                  &
      2.0_r_tran*s3(:) * field(:,1)                                            &
      + (-s2(:) -7.0_r_tran*s3(:)) * field(:,2)                                &
      + (2.0_r_tran*s1(:) + 5.0_r_tran*s2(:) + 11.0_r_tran*s3(:)) * field(:,3) &
      + (5.0_r_tran*s1(:) + 2.0_r_tran*s2(:)) * field(:,4)                     &
      - s1(:) * field(:,5)                                                     &
    )

  end subroutine nirvana_special_edge

  !----------------------------------------------------------------------------
  !> @brief  Calculates the interpolated field at the edges of a cell required for
  !!         using horizontal Nirvana to estimate the quadratic subgrid representation
  !!         of a field. The special edge treatment is used to shift the stencil left
  !!         or right depending on where the panel edge lies. At panel edges the
  !!         interpolation reduces to linear (it is quadratic elsewhere) and so
  !!         this requires a stencil of size 3. The field dofmap is
  !!         | 1 | 2 | 3 | and this returns the estimated field value at
  !!         the edges of cell 2. The cells are assumed to be uniform in spacing.
  !!
  !> @param[in]     field       Has dof map of the form | 1 | 2 | 3 |
  !> @param[in]     spt_case    Special edges case to shift stencil
  !> @param[in,out] edge_left   Field value at left edge of cell
  !> @param[in,out] edge_right  Field value at right edge of cell
  !> @param[in]     nlayers     Number of layers in the mesh
  !----------------------------------------------------------------------------
  subroutine linear_special_edge(field, spt_case, edge_left, edge_right, nlayers)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    real(kind=r_tran),   intent(in)    :: field(nlayers,3)
    integer(kind=i_def), intent(in)    :: spt_case(nlayers)
    real(kind=r_tran),   intent(inout) :: edge_left(nlayers)
    real(kind=r_tran),   intent(inout) :: edge_right(nlayers)

    real(kind=r_tran) :: s1(nlayers), s2(nlayers), s3(nlayers)
    real(kind=r_tran), parameter :: sixth = 1.0_r_tran/6.0_r_tran

    ! Switches that are 1 for particular spt_case values, and 0 otherwise
    s1(:) = MERGE(1.0_r_tran, 0.0_r_tran, spt_case(:) == 1)
    s2(:) = MERGE(1.0_r_tran, 0.0_r_tran, spt_case(:) == 2)
    s3(:) = MERGE(1.0_r_tran, 0.0_r_tran, spt_case(:) == 3)

    ! The switches turn on the coefficients for the particular spt_case. If
    ! none are turned on, revert to normal quadratic reconstruction
    edge_left(:) = sixth * (                                                   &
      (2.0_r_tran*s2(:) + 3.0_r_tran*s3(:)) * field(:,1)                       &
      + (9.0_r_tran*s1(:) + 5.0_r_tran*s2(:) + 3.0_r_tran*s3(:)) * field(:,2)  &
      + (-3.0_r_tran*s1(:) - s2(:)) * field(:,3)                               &
    )
    edge_right(:) = sixth * (                                                  &
      (-s2(:) - 3.0_r_tran*s3(:)) * field(:,1)                                 &
      + (3.0_r_tran*s1(:) + 5.0_r_tran*s2(:) + 9.0_r_tran*s3(:)) * field(:,2)  &
      + (3.0_r_tran*s1(:) + 2.0_r_tran*s2(:)) * field(:,3)                     &
    )

  end subroutine linear_special_edge

  !----------------------------------------------------------------------------
  !> @brief  Calculates the interpolated field at the edges of a cell required for
  !!         using horizontal PPM to estimate the quadratic subgrid representation
  !!         of a field. The special edge treatment is used to shift the stencil left
  !!         or right depending on where the panel edge lies. Fourth order interpolation
  !!         is used and so the stencil size increases to 7 instead of 5 for usual fourth
  !!         order interpolation. The field dofmap is | 1 | 2 | 3 | 4 | 5 | 6 | 7 | and
  !!         this returns the estimated field value at the edges of cell 4.
  !!         The cells are assumed to be uniform in spacing.
  !!
  !> @param[in]     field       Has dof map of the form | 1 | 2 | 3 | 4 | 5 | 6 | 7 |
  !> @param[in]     spt_case    Special edges case to shift stencil
  !> @param[in,out] edge_left   Field value at left edge of cell
  !> @param[in,out] edge_right  Field value at right edge of cell
  !> @param[in]     nlayers     Number of layers in the mesh
  !----------------------------------------------------------------------------
  subroutine fourth_order_special_edge(field, spt_case, edge_left, edge_right, nlayers)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    real(kind=r_tran),   intent(in)    :: field(nlayers,7)
    integer(kind=i_def), intent(in)    :: spt_case(nlayers)
    real(kind=r_tran),   intent(inout) :: edge_left(nlayers)
    real(kind=r_tran),   intent(inout) :: edge_right(nlayers)

    real(kind=r_tran) :: s1(nlayers), s2(nlayers), s3(nlayers)
    real(kind=r_tran) :: s4(nlayers), s5(nlayers)
    real(kind=r_tran), parameter :: twelfth = 1.0_r_tran/12.0_r_tran

    ! Switches that are 1 for particular spt_case values, and 0 otherwise
    s1(:) = MERGE(1.0_r_tran, 0.0_r_tran, spt_case(:) == 1)
    s2(:) = MERGE(1.0_r_tran, 0.0_r_tran, spt_case(:) == 2)
    s3(:) = MERGE(1.0_r_tran, 0.0_r_tran, spt_case(:) == 3)
    s4(:) = MERGE(1.0_r_tran, 0.0_r_tran, spt_case(:) == 4)
    s5(:) = MERGE(1.0_r_tran, 0.0_r_tran, spt_case(:) == 5)

    ! The switches turn on the coefficients for the particular spt_case. If
    ! none are turned on, revert to normal quadratic reconstruction
    edge_left(:) = twelfth * (                                                 &
      s4(:) * field(:,1)                                                       &
      + (-s3(:) - 5.0_r_tran*s4(:) - s5(:))*field(:,2)                         &
      + (3.0_r_tran*s1(:) + 7.0_r_tran*(s3(:) + s5(:))                         &
        + 13.0_r_tran*s4(:)) * field(:,3)                                      &
      + (13.0_r_tran*s1(:) + 25.0_r_tran*s2(:)                                 &
        + 7.0_r_tran*s3(:) + 3.0_r_tran*s4(:) + 7.0_r_tran*s5(:)) * field(:,4) &
      + (-5.0_r_tran*s1(:) - 23.0_r_tran*s2(:) - s3(:) - s5(:)) * field(:,5)   &
      + (s1(:) + 13.0_r_tran*s2(:)) * field(:,6)                               &
      - 3.0_r_tran*s2(:) * field(:,7)                                          &
    )
    edge_right(:) = twelfth * (                                                &
      (-3.0_r_tran*s4(:)) * field(:,1)                                         &
      + (13.0_r_tran*s4(:) + s5(:)) * field(:,2)                               &
      + (-s1(:) - s3(:) - 23.0_r_tran*s4(:) - 5.0_r_tran*s5(:)) * field(:,3)   &
      + (7.0_r_tran*s1(:) + 3.0_r_tran*s2(:) + 7.0_r_tran*s3(:)                &
        + 25.0_r_tran*s4(:) + 13.0_r_tran*s5(:)) * field(:,4)                  &
      + (7.0_r_tran*(s1(:) + s3(:)) + 13.0_r_tran*s2(:)                        &
        + 3.0_r_tran*s5(:)) * field(:,5)                                       &
      + (-s1(:) - 5.0_r_tran*s2(:) - s3(:)) * field(:,6)                       &
      + s2(:) * field(:,7)                                                     &
    )

  end subroutine fourth_order_special_edge

  !----------------------------------------------------------------------------
  !> @brief  Calculates the interpolated field at the edges of a cell required for
  !!         using horizontal PPM to estimate the quadratic subgrid representation
  !!         of a field. The special edge treatment is used to shift the stencil left
  !!         or right depending on where the panel edge lies, and reverts from fourth
  !!         order interpolation to Nirvana (third-order) interpolation at panel edges.
  !!         This does not require a larger stencil than usual fourth-order interpolation.
  !!         The field dofmap is | 1 | 2 | 3 | 4 | 5 | and this returns the estimated
  !!         field value at the edges of cell 3. The cells are assumed to be uniform in spacing.
  !!
  !> @param[in]     field       Has dof map of the form | 1 | 2 | 3 | 4 | 5 |
  !> @param[in]     spt_case    Special edges case to shift stencil
  !> @param[in,out] edge_left   Field value at left edge of cell
  !> @param[in,out] edge_right  Field value at right edge of cell
  !> @param[in]     nlayers     Number of layers in the mesh
  !----------------------------------------------------------------------------
  subroutine fourth_nirvana_special_edge(field, spt_case, edge_left, edge_right, nlayers)

    implicit none

    integer(kind=i_def), intent(in)    :: nlayers
    real(kind=r_tran),   intent(in)    :: field(nlayers,5)
    integer(kind=i_def), intent(in)    :: spt_case(nlayers)
    real(kind=r_tran),   intent(inout) :: edge_left(nlayers)
    real(kind=r_tran),   intent(inout) :: edge_right(nlayers)

    real(kind=r_tran) :: s2(nlayers), s3(nlayers), s4(nlayers), sn(nlayers)
    real(kind=r_tran), parameter :: twelfth = 1.0_r_tran/12.0_r_tran

    ! Switches that are 1 for particular spt_case values, and 0 otherwise
    s2(:) = MERGE(1.0_r_tran, 0.0_r_tran, spt_case(:) == 2)
    s3(:) = MERGE(1.0_r_tran, 0.0_r_tran, spt_case(:) == 3)
    s4(:) = MERGE(1.0_r_tran, 0.0_r_tran, spt_case(:) == 4)
    sn(:) = (                                                                  &
      MERGE(1.0_r_tran, 0.0_r_tran, spt_case(:) == 1)                          &
      + MERGE(1.0_r_tran, 0.0_r_tran, spt_case(:) == 5)                        &
    )

    ! The switches turn on the coefficients for the particular spt_case. If
    ! none are turned on, revert to normal quadratic reconstruction
    edge_left(:) = twelfth * (                                                 &
      (-s3(:) - 2.0_r_tran*s4(:)) * field(:,1)                                 &
      + (4.0_r_tran*sn(:) + 7.0_r_tran*s3(:) + 10.0_r_tran*s4(:)) * field(:,2) &
      + (10.0_r_tran*sn(:) + 22.0_r_tran*s2(:)                                 &
        + 7.0_r_tran*s3(:) + 4.0_r_tran*s4(:)) * field(:,3)                    &
      + (-2.0_r_tran*sn(:) - 14.0_r_tran*s2(:) - s3(:)) * field(:,4)           &
      + 4.0_r_tran*s2(:)*field(:,5)                                            &
    )
    edge_right(:) = twelfth * (                                                &
      4.0_r_tran*s4(:) * field(:,1)                                            &
      + (-2.0_r_tran*sn(:) - s3(:) - 14.0_r_tran*s4(:)) * field(:,2)           &
      + (10.0_r_tran*sn(:) + 4.0_r_tran*s2(:)                                  &
        + 7.0_r_tran*s3(:) + 22.0_r_tran*s4(:)) * field(:,3)                   &
      + (4.0_r_tran*sn(:) + 10.0_r_tran*s2(:) + 7.0_r_tran*s3(:)) * field(:,4) &
      + (-2.0_r_tran*s2(:) - s3(:)) * field(:,5)                               &
    )

  end subroutine fourth_nirvana_special_edge

end module subgrid_horizontal_support_mod
