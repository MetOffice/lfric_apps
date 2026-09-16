!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Dutton clear air turbulence predictor on one pressure level.

module dutton_cat_kernel_mod

  use argument_mod,      only: arg_type,                   &
                               GH_FIELD, GH_SCALAR,        &
                               GH_REAL, GH_INTEGER,        &
                               GH_READ, GH_READWRITE,      &
                               CELL_COLUMN,                &
                               ANY_DISCONTINUOUS_SPACE_1,  &
                               ANY_DISCONTINUOUS_SPACE_2,  &
                               STENCIL, CROSS2D
  use constants_mod,     only: r_def, i_def, rmdi
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: dutton_cat_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                             &
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE,              &
                  ANY_DISCONTINUOUS_SPACE_1),                      & ! cat
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,                   &
                  ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS2D)),    & ! u_plev
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,                   &
                  ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS2D)),    & ! v_plev
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,                   &
                  ANY_DISCONTINUOUS_SPACE_1),                      & ! dwdz
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,                   &
                  ANY_DISCONTINUOUS_SPACE_2, STENCIL(CROSS2D)),    & ! latitude
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,                   &
                  ANY_DISCONTINUOUS_SPACE_2, STENCIL(CROSS2D)),    & ! longitude
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                  & ! kp
         arg_type(GH_SCALAR, GH_REAL,    GH_READ)                   & ! radius
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: dutton_cat_code
  end type dutton_cat_kernel_type

  public :: dutton_cat_code

  ! Empirical constants of Dutton (1980), Probability forecasts of clear
  ! air turbulence based on numerical model output, Met Mag 109, 293-310,
  ! as used by the UM PWS routine DuttonCAT (pws_cat_mod).
  !> Number of interpolation knots
  integer(kind=i_def), parameter :: knots = 3_i_def
  !> Values of Dutton's empirical indicator at the knots
  real(kind=r_def), parameter :: eval(knots) =                                &
       (/ 5.0_r_def, 7.5_r_def, 65.0_r_def /)
  !> Corresponding values of the final predictor
  real(kind=r_def), parameter :: catval(knots) =                              &
       (/ 0.0_r_def, 1.75_r_def, 7.5_r_def /)
  !> Offset and slope of the empirical indicator in the combined shear
  real(kind=r_def), parameter :: empire_offset = 10.5_r_def
  real(kind=r_def), parameter :: empire_slope = 1.5_r_def
  !> Horizontal shear is required in (m/s)/(100 km): SI is (m/s)/m
  real(kind=r_def), parameter :: shearh_scale = 100000.0_r_def
  !> Vertical shear is required in (m/s)/km: SI is (m/s)/m
  real(kind=r_def), parameter :: shearv_scale = 1000.0_r_def

  !> CROSS2D stencil branches (west, south, east, north in mesh orientation)
  integer(kind=i_def), parameter :: n_branch = 4_i_def
  !> Index of the neighbouring cell within a branch of a depth-1 stencil
  integer(kind=i_def), parameter :: neighbour = 2_i_def

contains

  !> @brief   Dutton clear air turbulence predictor on one pressure level.
  !> @details Forms the horizontal derivatives of both wind components on
  !>          the pressure level by centred differences over the four
  !>          face-neighbouring cells, combines them with the wind and the
  !>          vertical wind shear into Dutton's empirical indicator, and
  !>          interpolates that onto a predictor between 0 and 7.5. This is
  !>          a port of the UM PWS routines DiffX, DiffY (diff_mod) and
  !>          DuttonCAT (pws_cat_mod). The UM differences along the rows and
  !>          columns of a regular latitude-longitude grid; here the two
  !>          opposite neighbour pairs give two directional differences
  !>          whose displacements, in the local east-north tangent plane of
  !>          the cell, are solved for the eastward and northward
  !>          derivatives. On such a grid this reduces to the UM
  !>          centred-difference operator up to the chord-to-arc factor
  !>          sin(d)/d in the grid spacing d (5e-5 at one degree). A cell
  !>          lacking any of the four neighbours (the edge of a
  !>          limited-area domain) receives missing data, as the UM sets
  !>          at the edges of its domain, and so does a cell whose vertical
  !>          shear is missing. Only slot kp of the output is written.
  !> @param[in]     nlayers        Number of layers
  !> @param[in,out] cat            Clear air turbulence predictor on pressure
  !>                               levels
  !> @param[in]     u_plev         Zonal wind on pressure levels
  !> @param[in]     smap_u_size    Size of each branch of the u_plev stencil
  !> @param[in]     u_len          Maximum branch length of the u_plev stencil
  !> @param[in]     smap_u         Stencil dofmap for u_plev
  !> @param[in]     v_plev         Meridional wind on pressure levels
  !> @param[in]     smap_v_size    Size of each branch of the v_plev stencil
  !> @param[in]     v_len          Maximum branch length of the v_plev stencil
  !> @param[in]     smap_v         Stencil dofmap for v_plev
  !> @param[in]     dwdz           Vertical derivative of wind speed with
  !>                               respect to height on pressure levels
  !> @param[in]     latitude       Latitude of cell centres (radians)
  !> @param[in]     smap_lat_size  Size of each branch of the latitude stencil
  !> @param[in]     lat_len        Maximum branch length of the latitude
  !>                               stencil
  !> @param[in]     smap_lat       Stencil dofmap for latitude
  !> @param[in]     longitude      Longitude of cell centres (radians)
  !> @param[in]     smap_lon_size  Size of each branch of the longitude
  !>                               stencil
  !> @param[in]     lon_len        Maximum branch length of the longitude
  !>                               stencil
  !> @param[in]     smap_lon       Stencil dofmap for longitude
  !> @param[in]     kp             Slot of the pressure-level fields to use
  !> @param[in]     planet_radius  Planet radius (m)
  !> @param[in]     ndf_pl         Number of degrees of freedom per cell for
  !>                               the pressure-level fields
  !> @param[in]     undf_pl        Number of unique degrees of freedom for
  !>                               the pressure-level fields
  !> @param[in]     map_pl         Dofmap for the cell for the pressure-level
  !>                               fields
  !> @param[in]     ndf_2d         Number of degrees of freedom per cell for
  !>                               the coordinate fields
  !> @param[in]     undf_2d        Number of unique degrees of freedom for
  !>                               the coordinate fields
  !> @param[in]     map_2d         Dofmap for the cell for the coordinate
  !>                               fields
  subroutine dutton_cat_code(nlayers,                                 &
                             cat,                                     &
                             u_plev, smap_u_size, u_len, smap_u,      &
                             v_plev, smap_v_size, v_len, smap_v,      &
                             dwdz,                                    &
                             latitude, smap_lat_size, lat_len,        &
                             smap_lat,                                &
                             longitude, smap_lon_size, lon_len,       &
                             smap_lon,                                &
                             kp, planet_radius,                       &
                             ndf_pl, undf_pl, map_pl,                 &
                             ndf_2d, undf_2d, map_2d)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_pl, undf_pl
    integer(kind=i_def), intent(in), dimension(ndf_pl) :: map_pl
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d) :: map_2d

    integer(kind=i_def), intent(in) :: u_len, v_len, lat_len, lon_len
    integer(kind=i_def), intent(in), dimension(n_branch) :: smap_u_size
    integer(kind=i_def), intent(in), dimension(n_branch) :: smap_v_size
    integer(kind=i_def), intent(in), dimension(n_branch) :: smap_lat_size
    integer(kind=i_def), intent(in), dimension(n_branch) :: smap_lon_size
    integer(kind=i_def), intent(in),                                          &
         dimension(ndf_pl, u_len, n_branch)   :: smap_u
    integer(kind=i_def), intent(in),                                          &
         dimension(ndf_pl, v_len, n_branch)   :: smap_v
    integer(kind=i_def), intent(in),                                          &
         dimension(ndf_2d, lat_len, n_branch) :: smap_lat
    integer(kind=i_def), intent(in),                                          &
         dimension(ndf_2d, lon_len, n_branch) :: smap_lon

    real(kind=r_def), intent(inout), dimension(undf_pl) :: cat
    real(kind=r_def), intent(in),    dimension(undf_pl) :: u_plev
    real(kind=r_def), intent(in),    dimension(undf_pl) :: v_plev
    real(kind=r_def), intent(in),    dimension(undf_pl) :: dwdz
    real(kind=r_def), intent(in),    dimension(undf_2d) :: latitude
    real(kind=r_def), intent(in),    dimension(undf_2d) :: longitude

    integer(kind=i_def), intent(in) :: kp
    real(kind=r_def),    intent(in) :: planet_radius

    ! Internal variables
    integer(kind=i_def) :: slot, n
    real(kind=r_def) :: lat0, lon0, latn, dlon
    real(kind=r_def), dimension(n_branch) :: x_n, y_n, u_n, v_n
    real(kind=r_def) :: dx_a, dy_a, dx_b, dy_b, det
    real(kind=r_def) :: du_a, du_b, dv_a, dv_b
    real(kind=r_def) :: dudx, dudy, dvdx, dvdy
    real(kind=r_def) :: u0, v0, vmag, shearh, shearv, empire
    real(kind=r_def) :: alpha, catprob

    slot = map_pl(1) + kp - 1_i_def

    ! All four stencils are on the same 2D mesh, so their branch sizes
    ! agree and the u_plev stencil speaks for the others
    if ( minval(smap_u_size) < neighbour .or. dwdz(slot) == rmdi ) then
      ! A neighbour is missing: no centred difference, as at the UM
      ! domain edges (DiffX, DiffY), so the predictor is missing too
      cat(slot) = rmdi
    else
      lat0 = latitude(map_2d(1))
      lon0 = longitude(map_2d(1))

      ! Displacement of each neighbouring cell centre from this cell in
      ! the local east-north tangent plane, and the winds there
      do n = 1, n_branch
        latn = latitude(smap_lat(1, neighbour, n))
        dlon = longitude(smap_lon(1, neighbour, n)) - lon0
        x_n(n) = planet_radius * cos(latn) * sin(dlon)
        y_n(n) = planet_radius * ( cos(lat0) * sin(latn)                      &
                                 - sin(lat0) * cos(latn) * cos(dlon) )
        u_n(n) = u_plev(smap_u(1, neighbour, n) + kp - 1_i_def)
        v_n(n) = v_plev(smap_v(1, neighbour, n) + kp - 1_i_def)
      end do

      ! Centred differences over the two opposite pairs: a = branch 3 minus
      ! branch 1, b = branch 4 minus branch 2 (UM DiffX, DiffY lines
      ! 287-289 and 395-396 on a regular latitude-longitude grid, where
      ! pair a is purely zonal and pair b purely meridional)
      dx_a = x_n(3) - x_n(1)
      dy_a = y_n(3) - y_n(1)
      dx_b = x_n(4) - x_n(2)
      dy_b = y_n(4) - y_n(2)
      det  = dx_a * dy_b - dy_a * dx_b

      du_a = u_n(3) - u_n(1)
      du_b = u_n(4) - u_n(2)
      dv_a = v_n(3) - v_n(1)
      dv_b = v_n(4) - v_n(2)

      ! Eastward and northward derivatives from the two directional ones
      dudx = ( du_a * dy_b - du_b * dy_a ) / det
      dudy = ( du_b * dx_a - du_a * dx_b ) / det
      dvdx = ( dv_a * dy_b - dv_b * dy_a ) / det
      dvdy = ( dv_b * dx_a - dv_a * dx_b ) / det

      ! Dutton's empirical predictor (UM DuttonCAT lines 354-404)
      u0 = u_plev(slot)
      v0 = v_plev(slot)

      shearv = shearv_scale * dwdz(slot)
      vmag = u0**2 + v0**2
      if ( vmag > 0.0_r_def ) then
        shearh = ( shearh_scale / vmag ) *                                    &
                 (   u0 * v0 * dudx                                           &
                   - u0 * u0 * dudy                                           &
                   + v0 * v0 * dvdx                                           &
                   - u0 * v0 * dvdy )
      else
        shearh = 0.0_r_def
      end if
      ! The vector form of the horizontal shear changes sign in the
      ! southern hemisphere
      if ( lat0 < 0.0_r_def ) then
        shearh = -shearh
      end if

      empire = empire_offset + empire_slope * ( shearh + shearv )

      ! Piecewise-linear interpolation between the knots; at or above the
      ! last knot the predictor takes its maximum value
      if ( empire < eval(1) ) then
        catprob = catval(1)
      else if ( empire < eval(2) ) then
        alpha   = ( empire - eval(1) ) / ( eval(2) - eval(1) )
        catprob = alpha * catval(2) + ( 1.0_r_def - alpha ) * catval(1)
      else if ( empire < eval(3) ) then
        alpha   = ( empire - eval(2) ) / ( eval(3) - eval(2) )
        catprob = alpha * catval(3) + ( 1.0_r_def - alpha ) * catval(2)
      else
        catprob = catval(knots)
      end if

      cat(slot) = catprob
    end if

  end subroutine dutton_cat_code

end module dutton_cat_kernel_mod
