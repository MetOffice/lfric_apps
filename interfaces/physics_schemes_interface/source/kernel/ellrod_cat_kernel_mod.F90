!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some content in this file was generated or refactored with assistance from
! - Claude Code (Claude Opus 5.5), 2026-10-20.
!-----------------------------------------------------------------------------
!> @brief Probability of clear-air turbulence from Ellrod's first index.

module ellrod_cat_kernel_mod

  use argument_mod,  only: arg_type,                  &
                           GH_FIELD, GH_SCALAR,       &
                           GH_READ, GH_WRITE,         &
                           GH_REAL, GH_INTEGER,       &
                           CELL_COLUMN,               &
                           ANY_DISCONTINUOUS_SPACE_1, &
                           ANY_DISCONTINUOUS_SPACE_2, &
                           STENCIL, CROSS2D
  use constants_mod, only: r_def, i_def, rmdi
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: ellrod_cat_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                        &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                   STENCIL(CROSS2D)),          &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,  &
                                                   STENCIL(CROSS2D)),          &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2,  &
                                                   STENCIL(CROSS2D)),          &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2,  &
                                                   STENCIL(CROSS2D)),          &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             &
         arg_type(GH_SCALAR, GH_REAL,    GH_READ)                              &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: ellrod_cat_code
  end type ellrod_cat_kernel_type

  public :: ellrod_cat_code

  !> Slope of the fit of turbulence probability to Ellrod's index (% s2)
  real(kind=r_def), parameter :: prob_slope = 1383149.3_r_def
  !> Offset of the fit of turbulence probability to Ellrod's index (%)
  real(kind=r_def), parameter :: prob_offset = 2.7258455_r_def
  !> Probabilities below this are set to zero (%)
  real(kind=r_def), parameter :: prob_min = 3.5_r_def

contains

  !> @brief   Probability of clear-air turbulence from Ellrod's first index.
  !> @details Ellrod's index is the vertical wind shear times the total
  !>          deformation of the horizontal wind (Ellrod and Knapp, 1992).
  !>          The horizontal wind derivatives come from the two pairs of
  !>          opposite neighbours in the cross stencil, with each neighbour
  !>          placed in the cell's local east-north tangent plane, so they
  !>          hold on any spherical mesh. The index is converted to a
  !>          percentage probability by a linear fit to turbulence reports;
  !>          probabilities below 3.5% are set to zero. Missing data where a
  !>          neighbour is missing (limited-area edges) or the shear is
  !>          missing.
  !> @param[in]     nlayers      Number of layers
  !> @param[in,out] cat_prob     Turbulence probability on pressure
  !>                             levels (%)
  !> @param[in]     u_plev       Eastward wind on pressure levels (m s-1)
  !> @param[in]     smap_u_size  Number of cells in each stencil branch
  !> @param[in]     smap_u_max   Maximum number of cells in a branch
  !> @param[in]     smap_u       Stencil dofmap for the wind fields
  !> @param[in]     v_plev       Northward wind on pressure levels (m s-1)
  !> @param[in]     smap_v_size  Number of cells in each stencil branch
  !> @param[in]     smap_v_max   Maximum number of cells in a branch
  !> @param[in]     smap_v       Stencil dofmap for the wind fields
  !> @param[in]     wind_shear   Vertical wind shear on pressure levels (s-1)
  !> @param[in]     latitude     Latitude of the cell centre (radians)
  !> @param[in]     smap_ll_size Number of cells in each stencil branch
  !> @param[in]     smap_ll_max  Maximum number of cells in a branch
  !> @param[in]     smap_ll      Stencil dofmap for latitude
  !> @param[in]     longitude    Longitude of the cell centre (radians)
  !> @param[in]     smap_lo_size Number of cells in each stencil branch
  !> @param[in]     smap_lo_max  Maximum number of cells in a branch
  !> @param[in]     smap_lo      Stencil dofmap for longitude
  !> @param[in]     nplev        Number of pressure levels
  !> @param[in]     radius       Planet radius (m)
  !> @param[in]     ndf_pl       Number of degrees of freedom per cell for
  !>                             the pressure-level fields
  !> @param[in]     undf_pl      Number of unique degrees of freedom for
  !>                             the pressure-level fields
  !> @param[in]     map_pl       Dofmap for the pressure-level fields
  !> @param[in]     ndf_2d       Number of degrees of freedom per cell for
  !>                             latitude and longitude
  !> @param[in]     undf_2d      Number of unique degrees of freedom for
  !>                             latitude and longitude
  !> @param[in]     map_2d       Dofmap for latitude and longitude
  subroutine ellrod_cat_code(nlayers, cat_prob,                         &
                             u_plev, smap_u_size, smap_u_max, smap_u,   &
                             v_plev, smap_v_size, smap_v_max, smap_v,   &
                             wind_shear,                                &
                             latitude, smap_ll_size, smap_ll_max,       &
                             smap_ll,                                   &
                             longitude, smap_lo_size, smap_lo_max,      &
                             smap_lo,                                   &
                             nplev, radius,                             &
                             ndf_pl, undf_pl, map_pl,                   &
                             ndf_2d, undf_2d, map_2d)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_pl, undf_pl
    integer(kind=i_def), intent(in), dimension(ndf_pl) :: map_pl
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d) :: map_2d
    integer(kind=i_def), intent(in) :: smap_u_max, smap_v_max
    integer(kind=i_def), intent(in) :: smap_ll_max, smap_lo_max
    integer(kind=i_def), intent(in), dimension(4) :: smap_u_size
    integer(kind=i_def), intent(in), dimension(4) :: smap_v_size
    integer(kind=i_def), intent(in), dimension(4) :: smap_ll_size
    integer(kind=i_def), intent(in), dimension(4) :: smap_lo_size
    integer(kind=i_def), intent(in),                                   &
                         dimension(ndf_pl, smap_u_max, 4) :: smap_u
    integer(kind=i_def), intent(in),                                   &
                         dimension(ndf_pl, smap_v_max, 4) :: smap_v
    integer(kind=i_def), intent(in),                                   &
                         dimension(ndf_2d, smap_ll_max, 4) :: smap_ll
    integer(kind=i_def), intent(in),                                   &
                         dimension(ndf_2d, smap_lo_max, 4) :: smap_lo

    ! Arguments passed explicitly from algorithm
    real(kind=r_def), intent(inout), dimension(undf_pl) :: cat_prob
    real(kind=r_def), intent(in),    dimension(undf_pl) :: u_plev
    real(kind=r_def), intent(in),    dimension(undf_pl) :: v_plev
    real(kind=r_def), intent(in),    dimension(undf_pl) :: wind_shear
    real(kind=r_def), intent(in),    dimension(undf_2d) :: latitude
    real(kind=r_def), intent(in),    dimension(undf_2d) :: longitude
    integer(kind=i_def), intent(in) :: nplev
    real(kind=r_def),    intent(in) :: radius

    ! Internal variables
    integer(kind=i_def) :: branch, kp
    integer(kind=i_def) :: nbr_u(4), nbr_v(4)
    real(kind=r_def)    :: lat0, lat_n, dlon
    ! Eastward and northward position of each neighbour (m)
    real(kind=r_def)    :: x_nbr(4), y_nbr(4)
    ! Separations of the two pairs of opposite neighbours (m)
    real(kind=r_def)    :: ax, ay, bx, by, det
    real(kind=r_def)    :: dudx, dudy, dvdx, dvdy
    real(kind=r_def)    :: stretching, shearing, ellrod_index, prob

    ! Missing data at limited-area edges, where a neighbour is missing
    if ( minval(smap_u_size) < 2_i_def ) then
      do kp = 1, nplev
        cat_prob(map_pl(1) + kp - 1) = rmdi
      end do
      return
    end if

    ! Position of each neighbour in the local tangent plane. Branch 1 is
    ! opposite branch 3, and branch 2 opposite branch 4.
    lat0 = latitude(smap_ll(1, 1, 1))
    do branch = 1, 4
      lat_n = latitude(smap_ll(1, 2, branch))
      dlon  = longitude(smap_lo(1, 2, branch)) - longitude(smap_lo(1, 1, 1))
      x_nbr(branch) = radius * cos(lat_n) * sin(dlon)
      y_nbr(branch) = radius * ( cos(lat0) * sin(lat_n)                     &
                               - sin(lat0) * cos(lat_n) * cos(dlon) )
      nbr_u(branch) = smap_u(1, 2, branch)
      nbr_v(branch) = smap_v(1, 2, branch)
    end do

    ax = x_nbr(3) - x_nbr(1)
    ay = y_nbr(3) - y_nbr(1)
    bx = x_nbr(4) - x_nbr(2)
    by = y_nbr(4) - y_nbr(2)
    det = ax * by - ay * bx

    do kp = 1, nplev

      if ( wind_shear(map_pl(1) + kp - 1) == rmdi ) then
        cat_prob(map_pl(1) + kp - 1) = rmdi
        cycle
      end if

      ! Solve for the gradient that reproduces the differences across
      ! both pairs of opposite neighbours
      call gradient(u_plev(nbr_u(3) + kp - 1) - u_plev(nbr_u(1) + kp - 1), &
                    u_plev(nbr_u(4) + kp - 1) - u_plev(nbr_u(2) + kp - 1), &
                    ax, ay, bx, by, det, dudx, dudy)
      call gradient(v_plev(nbr_v(3) + kp - 1) - v_plev(nbr_v(1) + kp - 1), &
                    v_plev(nbr_v(4) + kp - 1) - v_plev(nbr_v(2) + kp - 1), &
                    ax, ay, bx, by, det, dvdx, dvdy)

      stretching = dudx - dvdy
      shearing   = dvdx + dudy
      ellrod_index = wind_shear(map_pl(1) + kp - 1)                        &
                     * sqrt( stretching**2 + shearing**2 )

      prob = prob_slope * ellrod_index + prob_offset
      if ( prob < prob_min ) prob = 0.0_r_def
      cat_prob(map_pl(1) + kp - 1) = prob

    end do

  end subroutine ellrod_cat_code

  !> @brief   Horizontal gradient from differences across two separations.
  !> @param[in]  diff_a Difference of the field across separation a
  !> @param[in]  diff_b Difference of the field across separation b
  !> @param[in]  ax     Eastward component of separation a (m)
  !> @param[in]  ay     Northward component of separation a (m)
  !> @param[in]  bx     Eastward component of separation b (m)
  !> @param[in]  by     Northward component of separation b (m)
  !> @param[in]  det    ax*by - ay*bx
  !> @param[out] dfdx   Eastward derivative
  !> @param[out] dfdy   Northward derivative
  subroutine gradient(diff_a, diff_b, ax, ay, bx, by, det, dfdx, dfdy)

    implicit none

    real(kind=r_def), intent(in)  :: diff_a, diff_b
    real(kind=r_def), intent(in)  :: ax, ay, bx, by, det
    real(kind=r_def), intent(out) :: dfdx, dfdy

    dfdx = ( diff_a * by - diff_b * ay ) / det
    dfdy = ( ax * diff_b - bx * diff_a ) / det

  end subroutine gradient

end module ellrod_cat_kernel_mod
