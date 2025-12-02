!-------------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief   Calculates the fields in x and y at time n+1 using linear
!!          semi-Lagrangian transport.
!> @details This kernel using linear interpolation to solve the one-dimensional
!!          advection equation in both x and y.
!!
!> @note This kernel only works when field is a W3/Wtheta field at lowest order.

module horizontal_linear_sl_kernel_mod

  use argument_mod,          only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   CELL_COLUMN, GH_WRITE,     &
                                   GH_READ, GH_SCALAR,        &
                                   ANY_DISCONTINUOUS_SPACE_1, &
                                   STENCIL, CROSS2D, GH_INTEGER
  use constants_mod,         only: i_def, r_tran, l_def
  use fs_continuity_mod,     only: W2H
  use kernel_mod,            only: kernel_type
  use reference_element_mod, only: W, E, S, N

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: horizontal_linear_sl_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                        &
        arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),  & ! field_out_x
        arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),  & ! field_out_y
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,   &
                                                            STENCIL(CROSS2D)), & ! field
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2H)                         & ! dep_pts
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: horizontal_linear_sl_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: horizontal_linear_sl_code

contains

  !> @brief Compute advective transport in x and y directions using 1D
  !!        Semi-Lagrangian schemes, with a linear reconstruction. This is the
  !!        "inner" step of a COSMIC splitting scheme.
  !> @param[in]     nlayers        Number of layers
  !> @param[in,out] field_x        Field at time n+1 in x direction
  !> @param[in,out] field_y        Field at time n+1 in y direction
  !> @param[in]     field          Field to transport
  !> @param[in]     stencil_sizes  Sizes of the branches of the cross stencil
  !> @param[in]     stencil_max    Maximum size of a cross stencil branch
  !> @param[in]     stencil_map    Dofmap for the field stencil
  !> @param[in]     dep_pts        Departure points
  !> @param[in]     ndf_wf         Num of DoFs for field per cell
  !> @param[in]     undf_wf        Num of DoFs in this partition for field
  !> @param[in]     map_wf         Map for field
  !> @param[in]     ndf_w2h        Num of DoFs for W2H per cell
  !> @param[in]     undf_w2h       Num of DoFs in this partition for W2H
  !> @param[in]     map_w2h        Map for W2H
  subroutine horizontal_linear_sl_code( nlayers,        &
                                        field_x,        &
                                        field_y,        &
                                        field,          &
                                        stencil_sizes,  &
                                        stencil_max,    &
                                        stencil_map,    &
                                        dep_pts,        &
                                        ndf_wf,         &
                                        undf_wf,        &
                                        map_wf,         &
                                        ndf_w2h,        &
                                        undf_w2h,       &
                                        map_w2h )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_wf
    integer(kind=i_def), intent(in) :: ndf_wf
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: stencil_max
    integer(kind=i_def), intent(in) :: stencil_sizes(4)

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_wf(ndf_wf)
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: stencil_map(ndf_wf,stencil_max,4)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: field_x(undf_wf)
    real(kind=r_tran),   intent(inout) :: field_y(undf_wf)
    real(kind=r_tran),   intent(in)    :: field(undf_wf)
    real(kind=r_tran),   intent(in)    :: dep_pts(undf_w2h)

    ! Internal arguments
    integer(kind=i_def) :: i
    integer(kind=i_def) :: stencil_extent_xl, stencil_extent_xr
    integer(kind=i_def) :: stencil_extent_yl, stencil_extent_yr
    integer(kind=i_def) :: stencil_map_x_1d(-stencil_max:stencil_max)
    integer(kind=i_def) :: stencil_map_y_1d(-stencil_max:stencil_max)

    ! Form X and Y 1D stencils
    stencil_extent_xl = stencil_sizes(1) - 1
    stencil_extent_xr = stencil_sizes(3) - 1
    stencil_extent_yl = stencil_sizes(2) - 1
    stencil_extent_yr = stencil_sizes(4) - 1

    do i = -stencil_extent_xl, 0
      stencil_map_x_1d(i) = stencil_map(1, 1-i, 1)
    end do
    do i = 1, stencil_extent_xr
      stencil_map_x_1d(i) = stencil_map(1, i+1, 3)
    end do

    do i = -stencil_extent_yl, 0
      stencil_map_y_1d(i) = stencil_map(1, 1-i, 2)
    end do
    do i = 1, stencil_extent_yr
      stencil_map_y_1d(i) = stencil_map(1, i+1, 4)
    end do

    ! X-calculation ------------------------------------------------------------
    call horizontal_linear_sl_1d( nlayers,           &
                                  .true.,            &
                                  field_x,           &
                                  field,             &
                                  stencil_extent_xl, &
                                  stencil_extent_xr, &
                                  stencil_max,       &
                                  stencil_map_x_1d,  &
                                  dep_pts,           &
                                  ndf_wf,            &
                                  undf_wf,           &
                                  map_wf,            &
                                  ndf_w2h,           &
                                  undf_w2h,          &
                                  map_w2h )

    ! Y-calculation ------------------------------------------------------------
    call horizontal_linear_sl_1d( nlayers,           &
                                  .false.,           &
                                  field_y,           &
                                  field,             &
                                  stencil_extent_yl, &
                                  stencil_extent_yr, &
                                  stencil_max,       &
                                  stencil_map_y_1d,  &
                                  dep_pts,           &
                                  ndf_wf,            &
                                  undf_wf,           &
                                  map_wf,            &
                                  ndf_w2h,           &
                                  undf_w2h,          &
                                  map_w2h )

  end subroutine horizontal_linear_sl_code

! ============================================================================ !
! SINGLE UNDERLYING 1D ROUTINE
! ============================================================================ !

  !> @brief General 1D calculation of linear Semi-Lagrangian advected field
  subroutine horizontal_linear_sl_1d( nlayers,          &
                                      x_direction,      &
                                      field_out,        &
                                      field,            &
                                      stencil_extent_l, &
                                      stencil_extent_r, &
                                      stencil_max,      &
                                      stencil_map,      &
                                      dep_pts,          &
                                      ndf_wf,           &
                                      undf_wf,          &
                                      map_wf,           &
                                      ndf_w2h,          &
                                      undf_w2h,         &
                                      map_w2h )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_wf
    integer(kind=i_def), intent(in) :: ndf_wf
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: stencil_extent_l
    integer(kind=i_def), intent(in) :: stencil_extent_r
    integer(kind=i_def), intent(in) :: stencil_max
    logical(kind=l_def), intent(in) :: x_direction

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_wf(ndf_wf)
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: stencil_map(-stencil_max:stencil_max)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: field_out(undf_wf)
    real(kind=r_tran),   intent(in)    :: field(undf_wf)
    real(kind=r_tran),   intent(in)    :: dep_pts(undf_w2h)

    ! Local arrays
    integer(kind=i_def) :: int_disp(nlayers+ndf_wf-1)
    integer(kind=i_def) :: sign_disp(nlayers+ndf_wf-1)
    integer(kind=i_def) :: rel_idx_hi(nlayers+ndf_wf-1)
    integer(kind=i_def) :: rel_idx(nlayers+ndf_wf-1)
    real(kind=r_tran)   :: displacement(nlayers+ndf_wf-1)
    real(kind=r_tran)   :: field_local(nlayers+ndf_wf-1,2)
    real(kind=r_tran)   :: xx(nlayers+ndf_wf-1)

    ! Local scalars
    integer(kind=i_def) :: j, k, nl
    integer(kind=i_def) :: w2h_df_l, w2h_df_r
    real(kind=r_tran)   :: direction

    ! nl = nlayers      for w3
    !    = nlayers+1    for wtheta
    nl = nlayers + ndf_wf - 1

    if (x_direction) then
      w2h_df_l = map_w2h(W)
      w2h_df_r = map_w2h(E)
      direction = 1.0_r_tran
    else
      ! y-direction
      w2h_df_l = map_w2h(S)
      w2h_df_r = map_w2h(N)
      direction = -1.0_r_tran
    end if

    ! ======================================================================== !
    ! Extract departure info
    ! ======================================================================== !

    if (ndf_wf == 1) then
      ! Advecting W3 field: average the dep distances from this cell's faces
      displacement(:) = 0.5_r_tran * direction * (                             &
        dep_pts(w2h_df_l : w2h_df_l+nl-1)                                      &
        + dep_pts(w2h_df_r : w2h_df_r+nl-1)                                    &
      )
    else
      ! Advecting Wtheta field:
      ! In top and bottom layers, take the dep distances for top/bottom layer
      displacement(1) = 0.5_r_tran * direction * (                             &
        dep_pts(w2h_df_l) + dep_pts(w2h_df_r)                                  &
      )
      if (nlayers > 1) then
        ! NB: nl = nlayers + 1
        displacement(2:nl-1) = 0.25_r_tran * direction * (                     &
          dep_pts(w2h_df_l : w2h_df_l+nl-3)                                    &
          + dep_pts(w2h_df_l+1 : w2h_df_l+nl-2)                                &
          + dep_pts(w2h_df_r : w2h_df_r+nl-3)                                  &
          + dep_pts(w2h_df_r+1 : w2h_df_r+nl-2)                                &
        )
      end if
      ! Top layer
      displacement(nl) = 0.5_r_tran * direction * (                            &
        dep_pts(w2h_df_l+nl-2) + dep_pts(w2h_df_r+nl-2)                        &
      )
    end if

    int_disp(:) = INT(displacement(:), i_def)
    xx(:) = ABS(displacement(:) - REAL(int_disp, r_tran))
    sign_disp(:) = INT(SIGN(1.0_r_tran, displacement(:)))

    ! The relative index of the furthest cell to use in the stencil
    rel_idx_hi(:) = - sign_disp(:) - int_disp(:)

    ! ======================================================================== !
    ! Populate local arrays for interpolation
    ! ======================================================================== !

    ! Loop over points to use in reconstruction
    do j = 1, 2
      ! departure cell, between -stencil_extent_l and stencil_extent_r, e.g.
      ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
      rel_idx(:) = MIN(stencil_extent_r, MAX(-stencil_extent_l,                &
          rel_idx_hi(:) + (2 - j)*sign_disp(:)                                 &
      ))

      ! Loop over layers
      do k = 1, nl
        field_local(k,j) = field(stencil_map(rel_idx(k))+k-1)
      end do
    end do

    ! ======================================================================== !
    ! Perform linear interpolation
    ! ======================================================================== !

    ! Linear interpolation in x
    ! interp = (x-x1) / (x0-x1) f(x0) + (x-x0) / (x1-x0) f(x1)
    ! Set x0 = 0, x1 = 1, and 0 <= x <= 1
    ! interp = - (x-1) f(0) + x f(1)
    field_out(map_wf(1) : map_wf(1)+nl-1) = (                                  &
      -(xx(:)-1.0_r_tran) * field_local(:,1) + xx(:) * field_local(:,2)        &
    )

  end subroutine horizontal_linear_sl_1d

end module horizontal_linear_sl_kernel_mod
