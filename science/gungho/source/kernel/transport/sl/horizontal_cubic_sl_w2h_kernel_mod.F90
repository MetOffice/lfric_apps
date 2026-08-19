!-------------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief   Calculates the advective increments in x and y at time n+1 using
!!          cubic semi-Lagrangian transport.
!> @details This kernel using cubic interpolation to solve the one-dimensional
!!          advection equation in both x and y, giving advective increments
!!          in both directions. This is the second part of the COSMIC splitting,
!!          so the x increment works on the field previously advected in the
!!          y-direction(and vice versa).
!!
!> @note This kernel only works when field is a W2h field at lowest order.

module horizontal_cubic_sl_w2h_kernel_mod

  use argument_mod,          only: arg_type,                     &
                                   GH_FIELD, GH_REAL,            &
                                   CELL_COLUMN, GH_WRITE,        &
                                   GH_READ, GH_SCALAR,           &
                                   STENCIL, CROSS2D, GH_INTEGER
  use constants_mod,         only: r_tran, i_def, l_def
  use fs_continuity_mod,     only: W2H
  use kernel_mod,            only: kernel_type
  use reference_element_mod, only: W, E, S, N

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: horizontal_cubic_sl_w2h_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                   &
        arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2H),                   & ! increment_x
        arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, W2H),                   & ! increment_y
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2H, STENCIL(CROSS2D)), & ! field_x
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2H, STENCIL(CROSS2D)), & ! field_y
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2H, STENCIL(CROSS2D)), & ! dep_pts
        arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                          & ! monotone
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: horizontal_cubic_sl_w2h_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: horizontal_cubic_sl_w2h_code
  public :: horizontal_cubic_sl_w2h_1d

contains

  !> @brief Compute advective transport in x and y directions using 1D
  !!        Semi-Lagrangian schemes, with a cubic reconstruction. This is the
  !!        "outer" step of a COSMIC splitting scheme.
  !> @param[in]     nlayers           Number of layers
  !> @param[in,out] increment_x       Advective increment in x direction
  !> @param[in,out] increment_y       Advective increment in y direction
  !> @param[in]     field_x           Field from x direction
  !> @param[in]     stencil_sizes_x   Sizes of the branches of the cross stencil
  !> @param[in]     stencil_max_x     Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_x     Dofmap for the field_x stencil
  !> @param[in]     field_y           Field from y direction
  !> @param[in]     stencil_sizes_y   Sizes of the branches of the cross stencil
  !> @param[in]     stencil_max_y     Maximum size of a cross stencil branch
  !> @param[in]     stencil_map_y     Dofmap for the field_y stencil
  !> @param[in]     dep_pts           Departure points
  !> @param[in]     d_stencil_sizes   Sizes of branches of the departure-point
  !!                                  cross stencil
  !> @param[in]     d_stencil_max     Maximum size of a departure-point stencil
  !!                                  branch
  !> @param[in]     d_stencil_map     Dofmap for the departure-point stencil
  !> @param[in]     monotone          Horizontal monotone option for cubic SL
  !> @param[in]     ndf_w2h           Num of DoFs for W2H per cell
  !> @param[in]     undf_w2h          Num of DoFs for this partition for W2H
  !> @param[in]     map_w2h           Map for W2H
  subroutine horizontal_cubic_sl_w2h_code( nlayers,         &
                                           increment_x,     &
                                           increment_y,     &
                                           field_x,         &
                                           stencil_sizes_x, &
                                           stencil_max_x,   &
                                           stencil_map_x,   &
                                           field_y,         &
                                           stencil_sizes_y, &
                                           stencil_max_y,   &
                                           stencil_map_y,   &
                                           dep_pts,         &
                                           d_stencil_sizes, &
                                           d_stencil_max,   &
                                           d_stencil_map,   &
                                           monotone,        &
                                           ndf_w2h,         &
                                           undf_w2h,        &
                                           map_w2h )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: stencil_max_x
    integer(kind=i_def), intent(in) :: stencil_max_y
    integer(kind=i_def), intent(in) :: stencil_sizes_x(4)
    integer(kind=i_def), intent(in) :: stencil_sizes_y(4)
    integer(kind=i_def), intent(in) :: d_stencil_max
    integer(kind=i_def), intent(in) :: d_stencil_sizes(4)
    integer(kind=i_def), intent(in) :: monotone

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: stencil_map_x(ndf_w2h,stencil_max_x,4)
    integer(kind=i_def), intent(in) :: stencil_map_y(ndf_w2h,stencil_max_y,4)
    integer(kind=i_def), intent(in) :: d_stencil_map(ndf_w2h,d_stencil_max,4)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: increment_x(undf_w2h)
    real(kind=r_tran),   intent(inout) :: increment_y(undf_w2h)
    real(kind=r_tran),   intent(in)    :: field_x(undf_w2h)
    real(kind=r_tran),   intent(in)    :: field_y(undf_w2h)
    real(kind=r_tran),   intent(in)    :: dep_pts(undf_w2h)

    ! Internal arguments
    integer(kind=i_def) :: i, df
    integer(kind=i_def) :: stencil_extent_xl, stencil_extent_xr
    integer(kind=i_def) :: stencil_extent_yl, stencil_extent_yr
    integer(kind=i_def) :: stencil_map_x_1d(-stencil_max_x:stencil_max_x, ndf_w2h)
    integer(kind=i_def) :: stencil_map_y_1d(-stencil_max_y:stencil_max_y, ndf_w2h)

    ! Form X and Y 1D stencils
    stencil_extent_xl = stencil_sizes_x(1) - 1
    stencil_extent_xr = stencil_sizes_x(3) - 1
    stencil_extent_yl = stencil_sizes_y(2) - 1
    stencil_extent_yr = stencil_sizes_y(4) - 1

    do df = 1, ndf_w2h
      do i = -stencil_extent_xl, 0
        stencil_map_x_1d(i,df) = stencil_map_x(df, 1-i, 1)
      end do
      do i = 1, stencil_extent_xr
        stencil_map_x_1d(i,df) = stencil_map_x(df, i+1, 3)
      end do

      do i = -stencil_extent_yl, 0
        stencil_map_y_1d(i,df) = stencil_map_y(df, 1-i, 2)
      end do
      do i = 1, stencil_extent_yr
        stencil_map_y_1d(i,df) = stencil_map_y(df, i+1, 4)
      end do
    end do

    ! X-calculation ------------------------------------------------------------
    call horizontal_cubic_sl_w2h_1d( nlayers,           &
                                     .true.,            &
                                     increment_x,       &
                                     field_y,           &
                                     stencil_extent_xl, &
                                     stencil_extent_xr, &
                                     stencil_max_x,     &
                                     stencil_map_x_1d,  &
                                     dep_pts,           &
                                     monotone,          &
                                     ndf_w2h,           &
                                     undf_w2h,          &
                                     d_stencil_sizes,   &
                                     d_stencil_map,     &
                                     d_stencil_max )

    ! Y-calculation ------------------------------------------------------------
    call horizontal_cubic_sl_w2h_1d( nlayers,           &
                                     .false.,           &
                                     increment_y,       &
                                     field_x,           &
                                     stencil_extent_yl, &
                                     stencil_extent_yr, &
                                     stencil_max_y,     &
                                     stencil_map_y_1d,  &
                                     dep_pts,           &
                                     monotone,          &
                                     ndf_w2h,           &
                                     undf_w2h,          &
                                     d_stencil_sizes,   &
                                     d_stencil_map,     &
                                     d_stencil_max )

  end subroutine horizontal_cubic_sl_w2h_code

! ============================================================================ !
! SINGLE UNDERLYING 1D ROUTINE
! ============================================================================ !

  !> @brief General 1D calculation of cubic Semi-Lagrangian advective increment
  subroutine horizontal_cubic_sl_w2h_1d( nlayers,          &
                                         x_direction,      &
                                         increment,        &
                                         field,            &
                                         stencil_extent_l, &
                                         stencil_extent_r, &
                                         stencil_max,      &
                                         stencil_map,      &
                                         dep_pts,          &
                                         monotone,         &
                                         ndf_w2h,          &
                                         undf_w2h,         &
                                         w2h_sizes,        &
                                         map_w2h,          &
                                         w2h_max )

    use transport_enumerated_types_mod, only: monotone_strict,  &
                                              monotone_relaxed, &
                                              monotone_positive

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: stencil_extent_l
    integer(kind=i_def), intent(in) :: stencil_extent_r
    integer(kind=i_def), intent(in) :: stencil_max
    integer(kind=i_def), intent(in) :: w2h_max
    integer(kind=i_def), intent(in) :: w2h_sizes(4)
    integer(kind=i_def), intent(in) :: monotone
    logical(kind=l_def), intent(in) :: x_direction

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h, w2h_max, 4)
    integer(kind=i_def), intent(in) :: stencil_map(-stencil_max:stencil_max, ndf_w2h)

    ! Arguments: Fields
    real(kind=r_tran),   intent(inout) :: increment(undf_w2h)
    real(kind=r_tran),   intent(in)    :: field(undf_w2h)
    real(kind=r_tran),   intent(in)    :: dep_pts(undf_w2h)

    ! Local arrays
    integer(kind=i_def) :: int_disp(nlayers)
    integer(kind=i_def) :: sign_disp(nlayers)
    integer(kind=i_def) :: rel_idx_hi_p(nlayers)
    integer(kind=i_def) :: rel_idx(nlayers)
    real(kind=r_tran)   :: displacement(nlayers)
    real(kind=r_tran)   :: field_out(nlayers)
    real(kind=r_tran)   :: field_local(nlayers,4)
    real(kind=r_tran)   :: q_max(nlayers)
    real(kind=r_tran)   :: q_min(nlayers)
    real(kind=r_tran)   :: xx(nlayers)

    ! Local scalars
    integer(kind=i_def) :: j, k, nl, df
    integer(kind=i_def) :: w2h_df_l, w2h_df_r
    integer(kind=i_def) :: w2h_df_lm1, w2h_df_rm1
    integer(kind=i_def) :: w2h_df_lp1, w2h_df_rp1
    real(kind=r_tran)   :: direction

    ! Cubic interpolation weights
    real(kind=r_tran), parameter :: x0 = 0.0_r_tran
    real(kind=r_tran), parameter :: x1 = 1.0_r_tran
    real(kind=r_tran), parameter :: x2 = 2.0_r_tran
    real(kind=r_tran), parameter :: x3 = 3.0_r_tran
    real(kind=r_tran), parameter :: den0 = 1.0_r_tran/((x0-x1)*(x0-x2)*(x0-x3))
    real(kind=r_tran), parameter :: den1 = 1.0_r_tran/((x1-x0)*(x1-x2)*(x1-x3))
    real(kind=r_tran), parameter :: den2 = 1.0_r_tran/((x2-x0)*(x2-x1)*(x2-x3))
    real(kind=r_tran), parameter :: den3 = 1.0_r_tran/((x3-x0)*(x3-x1)*(x3-x2))

    nl = nlayers

    if (x_direction) then
      w2h_df_l   = map_w2h(W, 1, 1)
      w2h_df_r   = map_w2h(E, 1, 1)
      if ( w2h_sizes(2) > 1 ) then
        w2h_df_lm1 = map_w2h(W, 2, 2)
        w2h_df_rm1 = map_w2h(E, 2, 2)
      else
        w2h_df_lm1 = map_w2h(W, 1, 1)
        w2h_df_rm1 = map_w2h(E, 1, 1)
      end if
      if ( w2h_sizes(4) > 1 ) then
        w2h_df_lp1 = map_w2h(W, 2, 4)
        w2h_df_rp1 = map_w2h(E, 2, 4)
      else
        w2h_df_lp1 = map_w2h(W, 1, 1)
        w2h_df_rp1 = map_w2h(E, 1, 1)
      end if
      direction = 1.0_r_tran
    else
      ! y-direction
      w2h_df_l   = map_w2h(S, 1, 1)
      w2h_df_r   = map_w2h(N, 1, 1)
      if ( w2h_sizes(1) > 1 ) then
        w2h_df_lm1 = map_w2h(S, 2, 1)
        w2h_df_rm1 = map_w2h(N, 2, 1)
      else
        w2h_df_lm1 = map_w2h(S, 1, 1)
        w2h_df_rm1 = map_w2h(N, 1, 1)
      end if
      if ( w2h_sizes(3) > 1 ) then
        w2h_df_lp1 = map_w2h(S, 2, 3)
        w2h_df_rp1 = map_w2h(N, 2, 3)
      else
        w2h_df_lp1 = map_w2h(S, 1, 1)
        w2h_df_rp1 = map_w2h(N, 1, 1)
      end if
      direction = -1.0_r_tran
    end if

    do df = 1, 4
      ! ===================================================================== !
      ! Extract departure info
      ! ===================================================================== !

      ! Advecting W2 field: average the dep distances to the cell's faces
      if (df == 1 ) then
        ! Transporting u
        if ( x_direction ) then
          displacement(:) = direction * dep_pts(w2h_df_l : w2h_df_l+nl-1)
        else
          displacement(:) = 0.25_r_tran * direction * (           &
                            dep_pts(w2h_df_l   : w2h_df_l+nl-1)   &
                          + dep_pts(w2h_df_r   : w2h_df_r+nl-1)   &
                          + dep_pts(w2h_df_lm1 : w2h_df_lm1+nl-1) &
                          + dep_pts(w2h_df_rm1 : w2h_df_rm1+nl-1) )
        end if
     elseif ( df == 2 ) then
        ! Transporting v
        if ( x_direction ) then
          displacement(:) = 0.25_r_tran * direction * (           &
                            dep_pts(w2h_df_l   : w2h_df_l+nl-1)   &
                          + dep_pts(w2h_df_r   : w2h_df_r+nl-1)   &
                          + dep_pts(w2h_df_lm1 : w2h_df_lm1+nl-1) &
                          + dep_pts(w2h_df_rm1 : w2h_df_rm1+nl-1) )
        else
          displacement(:) = direction * dep_pts(w2h_df_l : w2h_df_l+nl-1)
        end if
      elseif ( df == 3 ) then
        if ( x_direction ) then
          displacement(:) = direction * dep_pts(w2h_df_r : w2h_df_r+nl-1)
        else
          displacement(:) = 0.25_r_tran * direction * (           &
                            dep_pts(w2h_df_l   : w2h_df_l+nl-1)   &
                          + dep_pts(w2h_df_r   : w2h_df_r+nl-1)   &
                          + dep_pts(w2h_df_lp1 : w2h_df_lp1+nl-1) &
                          + dep_pts(w2h_df_rp1 : w2h_df_rp1+nl-1) )
        end if
      else
        if ( x_direction ) then
          displacement(:) = 0.25_r_tran * direction * (           &
                            dep_pts(w2h_df_l   : w2h_df_l+nl-1)   &
                          + dep_pts(w2h_df_r   : w2h_df_r+nl-1)   &
                          + dep_pts(w2h_df_lp1 : w2h_df_lp1+nl-1) &
                          + dep_pts(w2h_df_rp1 : w2h_df_rp1+nl-1) )
        else
          displacement(:) = direction * dep_pts(w2h_df_r : w2h_df_r+nl-1)
        end if
      end if
      int_disp(:) = INT(displacement(:), i_def)
      xx(:) = 1.0_r_tran + ABS(displacement(:) - REAL(int_disp, r_tran))
      sign_disp(:) = INT(SIGN(1.0_r_tran, displacement(:)))

      ! The relative index of the most downwind cell to use in the stencil
      rel_idx_hi_p(:) = - 2*sign_disp(:) - int_disp(:)

      ! ===================================================================== !
      ! Populate local arrays for interpolation
      ! ===================================================================== !

      ! Loop over points to use in reconstruction
      do j = 1, 4
        ! If this column has idx 0, find relative index for the column of the
        ! departure cell, between -stencil_extent_l and stencil_extent_r, e.g.
        ! Relative idx is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
        rel_idx(:) = MIN(stencil_extent_r, MAX(-stencil_extent_l,             &
            rel_idx_hi_p(:) + (4 - j)*sign_disp(:)                            &
        ))

        ! Loop over layers
        do k = 1, nl
          field_local(k,j) = field(stencil_map(rel_idx(k),df)+k-1)
        end do
      end do

      ! ===================================================================== !
      ! Perform cubic interpolation
      ! ===================================================================== !

      field_out = (                                                           &
          (xx(:)-x1) * (xx(:)-x2) * (xx(:)-x3) * den0 * field_local(:,1)      &
        + (xx(:)-x0) * (xx(:)-x2) * (xx(:)-x3) * den1 * field_local(:,2)      &
        + (xx(:)-x0) * (xx(:)-x1) * (xx(:)-x3) * den2 * field_local(:,3)      &
        + (xx(:)-x0) * (xx(:)-x1) * (xx(:)-x2) * den3 * field_local(:,4)      &
      )

      ! ===================================================================== !
      ! Apply monotonicity constraints
      ! ===================================================================== !

      select case (monotone)
      case (monotone_strict)
        ! Bound field by immediately neighbouring values
        q_min(:) = MIN(field_local(:,2), field_local(:,3))
        q_max(:) = MAX(field_local(:,2), field_local(:,3))
        field_out(:) = MIN(q_max(:), MAX(field_out(:), q_min(:)))

      case (monotone_relaxed)
        ! Bound field by all values in the stencil
        q_min(:) = MIN(                                                       &
            field_local(:,1), field_local(:,2),                               &
            field_local(:,3), field_local(:,4)                                &
        )
        q_max(:) = MAX(                                                       &
            field_local(:,1), field_local(:,2),                               &
            field_local(:,3), field_local(:,4)                                &
        )
        field_out(:) = MIN(q_max(:), MAX(field_out(:), q_min(:)))

      case (monotone_positive)
        ! Just make sure field out is positive
        field_out(:) = MAX(field_out(:), 0.0_r_tran)

      end select

      ! ===================================================================== !
      ! Compute increment
      ! ===================================================================== !

      increment(map_w2h(df,1,1) : map_w2h(df,1,1)+nl-1) = (                   &
        field_out(:) - field(map_w2h(df,1,1) : map_w2h(df,1,1)+nl-1)          &
      )
    end do

  end subroutine horizontal_cubic_sl_w2h_1d

end module horizontal_cubic_sl_w2h_kernel_mod
