!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Calculates the fields in x and y at time n+1 using linear
!!          semi-Lagrangian transport.
!> @details This kernel using linear interpolation to solve the one-dimensional
!!          advection equation in both x and y.
!!
!> @note This kernel only works when field is a W3/Wtheta field at lowest order.

module horizontal_linear_sl_kernel_mod

  use argument_mod,       only : arg_type,                  &
                                 GH_FIELD, GH_REAL,         &
                                 CELL_COLUMN, GH_WRITE,     &
                                 GH_READ, GH_SCALAR,        &
                                 ANY_DISCONTINUOUS_SPACE_1, &
                                 STENCIL, CROSS, GH_INTEGER
  use constants_mod,      only : i_def, r_tran
  use fs_continuity_mod,  only : W2h
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: horizontal_linear_sl_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                                        &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),                 & ! field_out_x
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),                 & ! field_out_y
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS)), & ! field
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2h),                                       & ! dep_pts
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ     )                                         & ! extent_size
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: horizontal_linear_sl_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: horizontal_linear_sl_code

contains

  !> @brief Compute the advective field in x and y using linear interpolation.
  !> @param[in]     nlayers        Number of layers
  !> @param[in,out] field_x        Field at time n+1 in x direction
  !> @param[in,out] field_y        Field at time n+1 in y direction
  !> @param[in]     field          Field to transport
  !> @param[in]     stencil_size_c Local length of field cross stencil
  !> @param[in]     stencil_map    Dofmap for the field stencil
  !> @param[in]     dep_pts        Departure points
  !> @param[in]     extent_size    Stencil extent needed for the LAM edge
  !> @param[in]     ndf_wf         Number of degrees of freedom for field per cell
  !> @param[in]     undf_wf        Number of unique degrees of freedom for field
  !> @param[in]     map_wf         Map for field
  !> @param[in]     ndf_w2h        Number of degrees of freedom for W2h per cell
  !> @param[in]     undf_w2h       Number of unique degrees of freedom for W2h
  !> @param[in]     map_w2h        Map for W2h

  subroutine horizontal_linear_sl_code( nlayers,        &
                                        field_x,        &
                                        field_y,        &
                                        field,          &
                                        stencil_size_c, &
                                        stencil_map,    &
                                        dep_pts,        &
                                        extent_size,    &
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
    integer(kind=i_def), intent(in) :: stencil_size_c

    ! Arguments: Maps
    integer(kind=i_def), dimension(ndf_wf),  intent(in) :: map_wf
    integer(kind=i_def), dimension(ndf_w2h), intent(in) :: map_w2h
    integer(kind=i_def), dimension(ndf_wf,stencil_size_c), intent(in) :: stencil_map

    ! Arguments: Fields
    real(kind=r_tran), dimension(undf_wf),  intent(inout) :: field_x
    real(kind=r_tran), dimension(undf_wf),  intent(inout) :: field_y
    real(kind=r_tran), dimension(undf_wf),  intent(in)    :: field
    real(kind=r_tran), dimension(undf_w2h), intent(in)    :: dep_pts
    integer(kind=i_def), intent(in)                       :: extent_size

    ! Local fields
    real(kind=r_tran)   :: departure_dist, departure_dist_w3, departure_dist_wt

    ! Interpolation coefficients
    real(kind=r_tran)   :: frac_d, xx, yy

    ! Indices
    integer(kind=i_def) :: nl, k_w2h, k, km1, kp1, int_d
    integer(kind=i_def) :: rel_idx_hi, rel_idx_lo, sten_idx_hi, sten_idx_lo

    ! Stencils
    integer(kind=i_def) :: stencil_size, stencil_half, lam_edge_size

    ! Cross stencil has order e.g.
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

    ! stencil is a cross stencil we need 1D stencil size from this
    stencil_size = (stencil_size_c + 1_i_def) / 2_i_def
    stencil_half = (stencil_size + 1_i_def) / 2_i_def

    ! Get size the stencil should be to check if we are at the edge of a LAM domain
    lam_edge_size = 4_i_def*extent_size+1_i_def

    if ( lam_edge_size > stencil_size_c) then

      ! At edge of LAM, so set output to zero
      do k = 0, nl
         field_x( map_wf(1) + k ) = 0.0_r_tran
         field_y( map_wf(1) + k ) = 0.0_r_tran
      end do

    else

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

        ! Set up linear interpolation in correct cell
        ! For extent=4 the indices are:
        ! Relative id  is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
        ! Stencil has order |  5 |  4 |  3 |  2 |  1 | 10 | 11 | 12 | 13 |
        ! Get relative id depending on sign
        if (departure_dist >= 0.0_r_tran) then
          rel_idx_hi = - int_d
          rel_idx_lo = rel_idx_hi - 1
          xx = 1.0_r_tran - frac_d
        else
          rel_idx_hi = 1 - int_d
          rel_idx_lo = rel_idx_hi - 1
          xx = - frac_d
        end if
        ! Convert relative id into stencil id
        sten_idx_hi = 1 + ABS(rel_idx_hi) + (stencil_size - 1)*(1 - SIGN(1, -rel_idx_hi))/2
        sten_idx_lo = 1 + ABS(rel_idx_lo) + (stencil_size - 1)*(1 - SIGN(1, -rel_idx_lo))/2

        ! Linear interpolation in x
        ! interp = (x-x1) / (x0-x1) f(x0) + (x-x0) / (x1-x0) f(x1)
        ! Set x0 = 0, x1 = 1, and 0 <= x <= 1
        ! interp = - (x-1) f(0) + x f(1)
        field_x(map_wf(1)+k) = -(xx-1.0_r_tran)*field(stencil_map(1,sten_idx_lo) + k) &
                               + xx*field(stencil_map(1,sten_idx_hi) + k)

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

        ! Set up linear interpolation in correct cell
        ! For extent=4 the indices are:
        ! Relative id  is   | -4 | -3 | -2 | -1 |  0 |  1 |  2 |  3 |  4 |
        ! Stencil has order |  9 |  8 |  7 |  6 |  1 | 14 | 15 | 16 | 17 |
        ! Get relative id depending on sign
        if (departure_dist >= 0.0_r_tran) then
          rel_idx_hi = - int_d
          rel_idx_lo = rel_idx_hi - 1
          yy = 1.0_r_tran - frac_d
        else
          rel_idx_hi = 1 - int_d
          rel_idx_lo = rel_idx_hi - 1
          yy = - frac_d
        end if
        ! Convert relative id into stencil id
        sten_idx_hi = stencil_half + ABS(rel_idx_hi)                    &
                      + (stencil_size - 1)*(1 - SIGN(1, -rel_idx_hi))/2 &
                      - (stencil_half - 1)*(1 - SIGN(1, abs(rel_idx_hi)-1))/2
        sten_idx_lo = stencil_half + ABS(rel_idx_lo)                    &
                      + (stencil_size - 1)*(1 - SIGN(1, -rel_idx_lo))/2 &
                      - (stencil_half - 1)*(1 - SIGN(1, abs(rel_idx_lo)-1))/2

        ! Linear interpolation in y
        ! interp = (y-y1) / (y0-y1) f(y0) + (y-y0) / (y1-y0) f(y1)
        ! Set y0 = 0, y1 = 1, and 0 <= y <= 1
        ! interp = - (y-1) f(0) + y f(1)
        field_y(map_wf(1)+k) = -(yy-1.0_r_tran)*field(stencil_map(1,sten_idx_lo) + k) &
                               + yy*field(stencil_map(1,sten_idx_hi) + k)

      end do ! vertical levels k

    end if

  end subroutine horizontal_linear_sl_code

end module horizontal_linear_sl_kernel_mod
