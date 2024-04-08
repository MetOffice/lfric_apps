!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Calculates the advective increments in x and y at time n+1 using cubic
!!          semi-Lagrangian transport.
!> @details This kernel using cubic interpolation to solve the one-dimensional
!!          advection equation in both x and y, giving advective increments
!!          in both directions. This is the second part of the COSMIC splitting,
!!          so the x increment works on the field previously advected in the y-direction
!!          (and vice versa).
!!
!> @note This kernel only works when field is a W3/Wtheta field at lowest order.

module horizontal_cubic_sl_kernel_mod

  use argument_mod,       only : arg_type,                   &
                                 GH_FIELD, GH_REAL,          &
                                 CELL_COLUMN, GH_WRITE,      &
                                 GH_READ, GH_SCALAR,         &
                                 STENCIL, CROSS, GH_INTEGER, &
                                 ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,      only : r_tran, i_def, r_def
  use transport_enumerated_types_mod, only : horizontal_monotone_strict,   &
                                             horizontal_monotone_relaxed,  &
                                             horizontal_monotone_positive
  use fs_continuity_mod,  only : W2h
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: horizontal_cubic_sl_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                                        &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),                 & ! increment_x
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),                 & ! increment_y
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS)), & ! field_x
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS)), & ! field_y
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2h),                                       & ! dep_pts
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     ),                                        & ! monotone
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     ),                                        & ! extent_size
         arg_type(GH_SCALAR, GH_REAL,    GH_READ     )                                         & ! dt
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: horizontal_cubic_sl_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: horizontal_cubic_sl_code

contains

  !> @brief Compute the advective increment in x using PPM for the advective fluxes.
  !> @param[in]     nlayers           Number of layers
  !> @param[in,out] increment_x       Advective increment in x direction
  !> @param[in,out] increment_y       Advective increment in y direction
  !> @param[in]     field_x           Field from x direction
  !> @param[in]     stencil_size_x    Local length of field_x stencil
  !> @param[in]     stencil_map_x     Dofmap for the field_x stencil
  !> @param[in]     field_y           Field from y direction
  !> @param[in]     stencil_size_y    Local length of field_y stencil
  !> @param[in]     stencil_map_y     Dofmap for the field_y stencil
  !> @param[in]     dep_pts           Departure points
  !> @param[in]     monotone          Horizontal monotone option for cubic SL
  !> @param[in]     extent_size       Stencil extent needed for the LAM edge
  !> @param[in]     dt                Time step
  !> @param[in]     ndf_wf            Number of degrees of freedom for field per cell
  !> @param[in]     undf_wf           Number of unique degrees of freedom for Wf
  !> @param[in]     map_wf            Map for Wf
  !> @param[in]     ndf_w2h           Number of degrees of freedom for W2h per cell
  !> @param[in]     undf_w2h          Number of unique degrees of freedom for W2h
  !> @param[in]     map_w2h           Map for W2h

  subroutine horizontal_cubic_sl_code( nlayers,        &
                                       increment_x,    &
                                       increment_y,    &
                                       field_x,        &
                                       stencil_size_x, &
                                       stencil_map_x,  &
                                       field_y,        &
                                       stencil_size_y, &
                                       stencil_map_y,  &
                                       dep_pts,        &
                                       monotone,       &
                                       extent_size,    &
                                       dt,             &
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
    integer(kind=i_def), intent(in) :: stencil_size_x
    integer(kind=i_def), intent(in) :: stencil_size_y
    integer(kind=i_def), intent(in) :: monotone
    integer(kind=i_def), intent(in) :: extent_size
    real(kind=r_tran),   intent(in) :: dt

    ! Arguments: Maps
    integer(kind=i_def), dimension(ndf_wf),  intent(in) :: map_wf
    integer(kind=i_def), dimension(ndf_w2h), intent(in) :: map_w2h
    integer(kind=i_def), dimension(ndf_wf,stencil_size_x), intent(in) :: stencil_map_x
    integer(kind=i_def), dimension(ndf_wf,stencil_size_y), intent(in) :: stencil_map_y

    ! Arguments: Fields
    real(kind=r_tran),   dimension(undf_wf),  intent(inout) :: increment_x
    real(kind=r_tran),   dimension(undf_wf),  intent(inout) :: increment_y
    real(kind=r_tran),   dimension(undf_wf),  intent(in)    :: field_x
    real(kind=r_tran),   dimension(undf_wf),  intent(in)    :: field_y
    real(kind=r_tran),   dimension(undf_w2h), intent(in)    :: dep_pts

    ! Variables for flux calculation
    real(kind=r_tran) :: departure_dist, departure_dist_w3, departure_dist_wt

    ! Cubic interpolation coefficients
    real(kind=r_tran)   :: frac_d, x0, x1, x2, x3, xx, yy
    real(kind=r_tran)   :: field_out_x, field_out_y
    real(kind=r_tran)   :: qx_max, qx_min, qy_max, qy_min
    real(kind=r_tran)   :: den0, den1, den2, den3

    ! Indices
    integer(kind=i_def) :: k, km1, kp1, k_w2h, int_d, nl
    integer(kind=i_def) :: rel_idx_hi, rel_idx_hi_p, rel_idx_lo, rel_idx_lo_m
    integer(kind=i_def) :: rel_idy_hi, rel_idy_hi_p, rel_idy_lo, rel_idy_lo_m
    integer(kind=i_def) :: sten_idx_hi, sten_idx_hi_p, sten_idx_lo, sten_idx_lo_m
    integer(kind=i_def) :: sten_idy_hi, sten_idy_hi_p, sten_idy_lo, sten_idy_lo_m

    ! Stencils
    integer(kind=i_def) :: stencil_half, stencil_size, lam_edge_size

    ! Stencil has order e.g.
    !                           | 17 |
    !                           | 16 |
    !                           | 15 |
    !                           | 14 |
    !       |  5 |  4 |  3 |  2 |  1 | 10 | 11 | 12 | 13 | for extent 4
    !                           |  6 |
    !                           |  7 |
    !                           |  8 |
    !                           |  9 |
    !
    ! Relative idx is     | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
    ! Stencil x has order |  5 |  4 |  3 |  2 |  1 | 10 | 11 | 12 | 13 |
    ! Stencil y has order |  9 |  8 |  7 |  6 |  1 | 14 | 15 | 16 | 17 |
    ! Advection calculated for centre cell, e.g. cell 1 of stencil

    ! nl = nlayers-1  for w3
    !    = nlayers    for wtheta
    nl = nlayers - 1 + (ndf_wf - 1)

    ! Use stencil_size_y as each stencil size should be equal - as stencil_size_y
    ! is a cross stencil we need 1D stencil size from this
    stencil_size = (stencil_size_y + 1_i_def) / 2_i_def
    stencil_half = (stencil_size + 1_i_def) / 2_i_def

    ! Get size the stencil should be to check if we are at the edge of a LAM domain
    lam_edge_size = 4_i_def*extent_size+1_i_def

    if (lam_edge_size > stencil_size_x) then

      ! At edge of LAM, so set output to zero
      do k = 0, nl
        increment_x( map_wf(1) + k ) = 0.0_r_tran
        increment_y( map_wf(1) + k ) = 0.0_r_tran
      end do

    else

      ! Not at edge of LAM so compute increments

      ! Set up cubic interpolation weights
      x0 = 0.0_r_tran
      x1 = 1.0_r_tran
      x2 = 2.0_r_tran
      x3 = 3.0_r_tran
      den0 = ((x0-x1)*(x0-x2)*(x0-x3))
      den1 = ((x1-x0)*(x1-x2)*(x1-x3))
      den2 = ((x2-x0)*(x2-x1)*(x2-x3))
      den3 = ((x3-x0)*(x3-x1)*(x3-x2))

      ! Loop over k levels
      do k = 0, nl

        k_w2h = min(k,nlayers-1)
        km1 = max(k-1,0)
        kp1 = min(k,nlayers-1)

        ! x direction departure distance at centre
        departure_dist_w3 = ( dep_pts( map_w2h(1) + k_w2h)+dep_pts( map_w2h(3) + k_w2h) )/2.0_r_tran
        departure_dist_wt = ( dep_pts( map_w2h(1) + km1)+dep_pts( map_w2h(3) + km1) + &
                              dep_pts( map_w2h(1) + kp1)+dep_pts( map_w2h(3) + kp1) )/4.0_r_tran
        ! Combine W3 and Wtheta distances so that this works for either space
        departure_dist = ((2 - ndf_wf) * departure_dist_w3                     &
                          + (ndf_wf - 1) * departure_dist_wt)

        ! Calculates number of cells of interest to move
        frac_d = departure_dist - int(departure_dist)
        int_d = int(departure_dist,i_def)

        ! Set up cubic interpolation in correct cell
        ! For extent=4 the indices are:
        ! Relative id  is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
        ! Stencil has order |  5 |  4 |  3 |  2 |  1 | 10 | 11 | 12 | 13 |
        ! Get relative id depending on sign
        if (departure_dist >= 0.0_r_tran) then
          rel_idx_hi   = - int_d
          rel_idx_hi_p = rel_idx_hi + 1
          rel_idx_lo   = rel_idx_hi - 1
          rel_idx_lo_m = rel_idx_lo - 1
          xx = 2.0_r_tran - frac_d
        else
          rel_idx_hi   = 1 - int_d
          rel_idx_hi_p = rel_idx_hi + 1
          rel_idx_lo   = rel_idx_hi - 1
          rel_idx_lo_m = rel_idx_lo - 1
          xx = 1.0_r_tran - frac_d
        end if
        ! Convert relative id into stencil id
        sten_idx_hi_p = 1 + ABS(rel_idx_hi_p) + (stencil_size - 1)*(1 - SIGN(1, -rel_idx_hi_p))/2
        sten_idx_hi   = 1 + ABS(rel_idx_hi) + (stencil_size - 1)*(1 - SIGN(1, -rel_idx_hi))/2
        sten_idx_lo   = 1 + ABS(rel_idx_lo) + (stencil_size - 1)*(1 - SIGN(1, -rel_idx_lo))/2
        sten_idx_lo_m = 1 + ABS(rel_idx_lo_m) + (stencil_size - 1)*(1 - SIGN(1, -rel_idx_lo_m))/2

        ! Cubic interpolation in x
        field_out_x = ((xx-x1)*(xx-x2)*(xx-x3))/(den0)*field_y(stencil_map_y(1,sten_idx_lo_m) + k) + &
                      ((xx-x0)*(xx-x2)*(xx-x3))/(den1)*field_y(stencil_map_y(1,sten_idx_lo) + k) +   &
                      ((xx-x0)*(xx-x1)*(xx-x3))/(den2)*field_y(stencil_map_y(1,sten_idx_hi) + k) +   &
                      ((xx-x0)*(xx-x1)*(xx-x2))/(den3)*field_y(stencil_map_y(1,sten_idx_hi_p) + k)

        ! y direction departure distance at centre
        departure_dist_w3 = ( dep_pts( map_w2h(2) + k_w2h)+dep_pts( map_w2h(4) + k_w2h) )/2.0_r_tran
        departure_dist_wt = ( dep_pts( map_w2h(2) + km1)+dep_pts( map_w2h(4) + km1) + &
                              dep_pts( map_w2h(2) + kp1)+dep_pts( map_w2h(4) + kp1) )/4.0_r_tran
        ! NB: minus sign as y-direction of stencils is defined in the opposite
        ! direction to the y-direction of the wind field
        ! Combine W3 and Wtheta distances so that this works for either space
        departure_dist = -((2 - ndf_wf) * departure_dist_w3                    &
                           + (ndf_wf - 1) * departure_dist_wt)

        ! Calculates number of cells of interest to move
        frac_d = departure_dist - int(departure_dist)
        int_d = int(departure_dist,i_def)

        ! Set up cubic interpolation in correct cell
        ! For extent=4 the indices are:
        ! Relative id  is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
        ! Stencil has order |  9 |  8 |  7 |  6 |  1 | 14 | 15 | 16 | 17 |
        ! Get relative id depending on sign
        if (departure_dist >= 0.0_r_tran) then
          rel_idy_hi   = - int_d
          rel_idy_hi_p = rel_idy_hi + 1
          rel_idy_lo   = rel_idy_hi - 1
          rel_idy_lo_m = rel_idy_lo - 1
          yy = 2.0_r_tran - frac_d
        else
          rel_idy_hi   = 1 - int_d
          rel_idy_hi_p = rel_idy_hi + 1
          rel_idy_lo   = rel_idy_hi - 1
          rel_idy_lo_m = rel_idy_lo - 1
          yy = 1.0_r_tran - frac_d
        end if
        ! Convert relative id into stencil id
        sten_idy_hi_p = stencil_half + ABS(rel_idy_hi_p)                    &
                        + (stencil_size - 1)*(1 - SIGN(1, -rel_idy_hi_p))/2 &
                        - (stencil_half - 1)*(1 - SIGN(1, abs(rel_idy_hi_p)-1))/2
        sten_idy_hi   = stencil_half + ABS(rel_idy_hi)                      &
                        + (stencil_size - 1)*(1 - SIGN(1, -rel_idy_hi))/2   &
                        - (stencil_half - 1)*(1 - SIGN(1, abs(rel_idy_hi)-1))/2
        sten_idy_lo   = stencil_half + ABS(rel_idy_lo)                      &
                        + (stencil_size - 1)*(1 - SIGN(1, -rel_idy_lo))/2   &
                        - (stencil_half - 1)*(1 - SIGN(1, abs(rel_idy_lo)-1))/2
        sten_idy_lo_m = stencil_half + ABS(rel_idy_lo_m)                    &
                        + (stencil_size - 1)*(1 - SIGN(1, -rel_idy_lo_m))/2 &
                        - (stencil_half - 1)*(1 - SIGN(1, abs(rel_idy_lo_m)-1))/2

        ! Cubic interpolation in y
        field_out_y = ((yy-x1)*(yy-x2)*(yy-x3))/(den0)*field_x(stencil_map_x(1,sten_idy_lo_m) + k) + &
                      ((yy-x0)*(yy-x2)*(yy-x3))/(den1)*field_x(stencil_map_x(1,sten_idy_lo) + k) +   &
                      ((yy-x0)*(yy-x1)*(yy-x3))/(den2)*field_x(stencil_map_x(1,sten_idy_hi) + k) +   &
                      ((yy-x0)*(yy-x1)*(yy-x2))/(den3)*field_x(stencil_map_x(1,sten_idy_hi_p) + k)

        ! Monotone
        if (monotone == horizontal_monotone_strict) then
          ! Get neighbouring field bounds
          qx_min = min( field_y(stencil_map_y(1,sten_idx_lo) + k), &
                        field_y(stencil_map_y(1,sten_idx_hi) + k) )
          qx_max = max( field_y(stencil_map_y(1,sten_idx_lo) + k), &
                        field_y(stencil_map_y(1,sten_idx_hi) + k) )
          qy_min = min( field_x(stencil_map_x(1,sten_idy_lo) + k), &
                        field_x(stencil_map_x(1,sten_idy_hi) + k) )
          qy_max = max( field_x(stencil_map_x(1,sten_idy_lo) + k), &
                        field_x(stencil_map_x(1,sten_idy_hi) + k) )
          field_out_x = min( qx_max, max( field_out_x, qx_min ) )
          field_out_y = min( qy_max, max( field_out_y, qy_min ) )
        else if (monotone == horizontal_monotone_relaxed) then
          ! Get stencil bounds
          qx_min = min( field_y(stencil_map_y(1,sten_idx_lo_m) + k), field_y(stencil_map_y(1,sten_idx_lo) + k), &
                        field_y(stencil_map_y(1,sten_idx_hi) + k), field_y(stencil_map_y(1,sten_idx_hi_p) + k) )
          qx_max = max( field_y(stencil_map_y(1,sten_idx_lo_m) + k), field_y(stencil_map_y(1,sten_idx_lo) + k), &
                        field_y(stencil_map_y(1,sten_idx_hi) + k), field_y(stencil_map_y(1,sten_idx_hi_p) + k) )
          qy_min = min( field_x(stencil_map_x(1,sten_idy_lo_m) + k), field_x(stencil_map_x(1,sten_idy_lo) + k), &
                        field_x(stencil_map_x(1,sten_idy_hi) + k), field_x(stencil_map_x(1,sten_idy_hi_p) + k) )
          qy_max = max( field_x(stencil_map_x(1,sten_idy_lo_m) + k), field_x(stencil_map_x(1,sten_idy_lo) + k), &
                        field_x(stencil_map_x(1,sten_idy_hi) + k), field_x(stencil_map_x(1,sten_idy_hi_p) + k) )
          field_out_x = min( qx_max, max( field_out_x, qx_min ) )
          field_out_y = min( qy_max, max( field_out_y, qy_min ) )
        else if (monotone == horizontal_monotone_positive) then
          ! Make sure field out is positive
          field_out_x = max( field_out_x, 0.0_r_tran )
          field_out_y = max( field_out_y, 0.0_r_tran )
        end if

        ! Get increment
        increment_x(map_wf(1)+k) = (field_y(stencil_map_y(1,1) + k) - field_out_x) / dt
        increment_y(map_wf(1)+k) = (field_x(stencil_map_x(1,1) + k) - field_out_y) / dt

      end do ! vertical levels k

    end if

  end subroutine horizontal_cubic_sl_code

end module horizontal_cubic_sl_kernel_mod
