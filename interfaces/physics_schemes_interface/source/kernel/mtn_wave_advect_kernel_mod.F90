!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some content in this file was generated or refactored with assistance from
! - Met Office Claude Code via AI Foundry (Claude Opus 5.5).

!> @brief Mountain wave turbulence predictor: the column maximum of the
!>        orographic drag stress, carried one cell downwind along the
!>        200 hPa wind.

module mtn_wave_advect_kernel_mod

  use argument_mod,        only: arg_type,                   &
                                 GH_FIELD, GH_SCALAR,        &
                                 GH_REAL, GH_INTEGER,        &
                                 GH_READ, GH_WRITE,          &
                                 CELL_COLUMN,                &
                                 ANY_DISCONTINUOUS_SPACE_1,  &
                                 ANY_DISCONTINUOUS_SPACE_2,  &
                                 STENCIL, REGION
  use constants_mod,       only: r_def, i_def, pi
  use coord_transform_mod, only: ll2xyz, cart2sphere_vector
  use kernel_mod,          only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: mtn_wave_advect_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                        &
         ! mtn_wave_turb
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,             &
                  ANY_DISCONTINUOUS_SPACE_1),                  &
         ! max_stress
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,              &
                  ANY_DISCONTINUOUS_SPACE_1, STENCIL(REGION)), &
         ! u_plev
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,              &
                  ANY_DISCONTINUOUS_SPACE_2, STENCIL(REGION)), &
         ! v_plev
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,              &
                  ANY_DISCONTINUOUS_SPACE_2, STENCIL(REGION)), &
         ! latitude
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,              &
                  ANY_DISCONTINUOUS_SPACE_1, STENCIL(REGION)), &
         ! longitude
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,              &
                  ANY_DISCONTINUOUS_SPACE_1, STENCIL(REGION)), &
         ! wind_level
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ)              &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: mtn_wave_advect_code
  end type mtn_wave_advect_kernel_type

  public :: mtn_wave_advect_code

  ! Empirical constants of Turner (1999), Development of a mountain wave
  ! turbulence prediction scheme for civil aviation, Met Office Forecasting
  ! Research Technical Report 265.
  !> Number of stress contours
  integer(kind=i_def), parameter :: n_cutoff = 3_i_def
  !> Stress contours, in ascending order (Pa)
  real(kind=r_def), parameter :: cutoff(n_cutoff) =                         &
       (/ 0.007_r_def, 0.065_r_def, 0.25_r_def /)
  !> Wind speed the source cell must exceed for its stress to be carried
  !> (m s-1)
  real(kind=r_def), parameter :: min_speed = 0.0_r_def

  !> Full turn (radians)
  real(kind=r_def), parameter :: full_turn = 2.0_r_def * pi
  !> Widest bearing gap between adjacent neighbours that still has a
  !> neighbour on each side (radians). Gaps are 45 degrees on a regular
  !> grid, about 60 at a cube corner and under 90 for a stretched cell; the
  !> edge of a limited-area domain leaves 180 or more.
  real(kind=r_def), parameter :: max_gap = 0.75_r_def * pi

contains

  !> @brief   Mountain wave turbulence predictor.
  !> @details Each cell takes the column-maximum stress of the neighbour the
  !>          200 hPa wind blows from, when a stress contour separates the
  !>          two and the wind at the neighbour is not calm; otherwise it
  !>          keeps its own.
  !>
  !>          The original method takes its neighbour upwind east-west but
  !>          downwind north-south; here it is upwind in both. It also
  !>          repeats each copy for sources faster than 30 m/s, which cannot
  !>          change a copied value, so that pass is omitted. It updates the
  !>          cells in place, one contour at a time, but a source exceeds
  !>          the contour and so is never itself overwritten: computing each
  !>          cell once from its neighbours' original values gives the same
  !>          result.
  !>
  !>          A cell with no neighbour on the upwind side, at the edge of a
  !>          limited-area domain, keeps its own value.
  !> @param[in]     nlayers        Number of layers
  !> @param[in,out] mtn_wave_turb  Mountain wave turbulence predictor (Pa)
  !> @param[in]     max_stress     Column maximum of the orographic drag
  !>                               stress (Pa)
  !> @param[in]     smap_ms_size   Number of cells in the max_stress stencil
  !> @param[in]     smap_ms        Stencil dofmap for max_stress
  !> @param[in]     u_plev         Zonal wind on pressure levels (m s-1)
  !> @param[in]     smap_u_size    Number of cells in the u_plev stencil
  !> @param[in]     smap_u         Stencil dofmap for u_plev
  !> @param[in]     v_plev         Meridional wind on pressure levels (m s-1)
  !> @param[in]     smap_v_size    Number of cells in the v_plev stencil
  !> @param[in]     smap_v         Stencil dofmap for v_plev
  !> @param[in]     latitude       Latitude of cell centres (radians)
  !> @param[in]     smap_lat_size  Number of cells in the latitude stencil
  !> @param[in]     smap_lat       Stencil dofmap for latitude
  !> @param[in]     longitude      Longitude of cell centres (radians)
  !> @param[in]     smap_lon_size  Number of cells in the longitude stencil
  !> @param[in]     smap_lon       Stencil dofmap for longitude
  !> @param[in]     wind_level     Slot of 200 hPa on the pressure-level axis
  !> @param[in]     ndf_2d         Number of degrees of freedom per cell for
  !>                               the 2D fields
  !> @param[in]     undf_2d        Number of unique degrees of freedom for the
  !>                               2D fields
  !> @param[in]     map_2d         Dofmap for the cell for the 2D fields
  !> @param[in]     ndf_pl         Number of degrees of freedom per cell for
  !>                               the pressure-level fields
  !> @param[in]     undf_pl        Number of unique degrees of freedom for the
  !>                               pressure-level fields
  !> @param[in]     map_pl         Dofmap for the cell for the pressure-level
  !>                               fields
  subroutine mtn_wave_advect_code(nlayers,                               &
                                  mtn_wave_turb,                         &
                                  max_stress, smap_ms_size, smap_ms,     &
                                  u_plev, smap_u_size, smap_u,           &
                                  v_plev, smap_v_size, smap_v,           &
                                  latitude, smap_lat_size, smap_lat,     &
                                  longitude, smap_lon_size, smap_lon,    &
                                  wind_level,                            &
                                  ndf_2d, undf_2d, map_2d,               &
                                  ndf_pl, undf_pl, map_pl)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: ndf_pl, undf_pl
    integer(kind=i_def), intent(in), dimension(ndf_2d) :: map_2d
    integer(kind=i_def), intent(in), dimension(ndf_pl) :: map_pl

    integer(kind=i_def), intent(in) :: smap_ms_size, smap_u_size
    integer(kind=i_def), intent(in) :: smap_v_size, smap_lat_size
    integer(kind=i_def), intent(in) :: smap_lon_size
    integer(kind=i_def), intent(in),                                       &
         dimension(ndf_2d, smap_ms_size)  :: smap_ms
    integer(kind=i_def), intent(in),                                       &
         dimension(ndf_pl, smap_u_size)   :: smap_u
    integer(kind=i_def), intent(in),                                       &
         dimension(ndf_pl, smap_v_size)   :: smap_v
    integer(kind=i_def), intent(in),                                       &
         dimension(ndf_2d, smap_lat_size) :: smap_lat
    integer(kind=i_def), intent(in),                                       &
         dimension(ndf_2d, smap_lon_size) :: smap_lon

    real(kind=r_def), intent(inout), dimension(undf_2d) :: mtn_wave_turb
    real(kind=r_def), intent(in),    dimension(undf_2d) :: max_stress
    real(kind=r_def), intent(in),    dimension(undf_pl) :: u_plev
    real(kind=r_def), intent(in),    dimension(undf_pl) :: v_plev
    real(kind=r_def), intent(in),    dimension(undf_2d) :: latitude
    real(kind=r_def), intent(in),    dimension(undf_2d) :: longitude

    integer(kind=i_def), intent(in) :: wind_level

    ! Internal variables
    integer(kind=i_def) :: n, n_nbr, source, slot, c
    real(kind=r_def)    :: x_centre(3), x_nbr(3), offset(3)
    real(kind=r_def)    :: east(smap_lat_size - 1), north(smap_lat_size - 1)
    real(kind=r_def)    :: u_wind, v_wind, upwind, value, source_stress

    ! All the stencils are on the same horizontal mesh, so they list the
    ! same cells in the same order: slot 1 is this cell, the rest are its
    ! neighbours
    n_nbr = smap_lat_size - 1

    ! Offset of each neighbour's centre from this cell's, east and north in
    ! the plane tangent to the unit sphere at this cell
    call ll2xyz(longitude(smap_lon(1, 1)), latitude(smap_lat(1, 1)),       &
                x_centre(1), x_centre(2), x_centre(3))
    do n = 1, n_nbr
      call ll2xyz(longitude(smap_lon(1, n + 1)),                           &
                  latitude(smap_lat(1, n + 1)),                            &
                  x_nbr(1), x_nbr(2), x_nbr(3))
      offset = cart2sphere_vector(x_centre, x_nbr - x_centre)
      east(n)  = offset(1)
      north(n) = offset(2)
    end do

    ! Bearing the wind blows from, clockwise from north. A calm wind counts
    ! as blowing towards the north, as in the original method.
    u_wind = u_plev(smap_u(1, 1) + wind_level - 1)
    v_wind = v_plev(smap_v(1, 1) + wind_level - 1)
    if ( u_wind == 0.0_r_def .and. v_wind == 0.0_r_def ) then
      upwind = pi
    else
      upwind = atan2(-u_wind, -v_wind)
    end if

    value = max_stress(smap_ms(1, 1))

    source = upwind_neighbour(east, north, upwind)
    if ( source > 0_i_def ) then
      slot = source + 1_i_def
      u_wind = u_plev(smap_u(1, slot) + wind_level - 1)
      v_wind = v_plev(smap_v(1, slot) + wind_level - 1)
      if ( sqrt(u_wind**2 + v_wind**2) > min_speed ) then
        source_stress = max_stress(smap_ms(1, slot))
        do c = 1, n_cutoff
          if ( source_stress > cutoff(c) .and. value < cutoff(c) ) then
            value = source_stress
          end if
        end do
      end if
    end if

    mtn_wave_turb(map_2d(1)) = value

  end subroutine mtn_wave_advect_code

  !> @brief   Neighbour the wind blows from.
  !> @details Each neighbour owns the directions closer to it than to the
  !>          adjacent neighbours either side, the boundary being the
  !>          bearing of the point halfway between the two. On a regular
  !>          latitude-longitude grid these are the sector edges of the
  !>          original method, 26.6 and 63.4 degrees from the grid axes for
  !>          square cells. A gap wider than max_gap between adjacent
  !>          neighbours has no neighbour in it.
  !> @param[in] east    Eastward offsets of the neighbours
  !> @param[in] north   Northward offsets of the neighbours
  !> @param[in] upwind  Bearing the wind blows from, clockwise from north
  !>                    (radians)
  !> @return    Index of the neighbour, or 0 if there is none that way
  pure function upwind_neighbour(east, north, upwind) result(source)

    implicit none

    real(kind=r_def), intent(in) :: east(:), north(:)
    real(kind=r_def), intent(in) :: upwind
    integer(kind=i_def)          :: source

    integer(kind=i_def) :: n, before, after
    real(kind=r_def)    :: bearing(size(east)), past(size(east))
    real(kind=r_def)    :: ahead(size(east)), halfway

    source = 0_i_def

    do n = 1, size(east)
      bearing(n) = atan2(east(n), north(n))
    end do

    ! The neighbours either side of the wind's bearing: the last one at or
    ! before it, clockwise, and the next one after that
    past = modulo(upwind - bearing, full_turn)
    before = minloc(past, 1)
    ahead = modulo(bearing - bearing(before), full_turn)
    ahead(before) = full_turn
    after = minloc(ahead, 1)

    if ( ahead(after) > max_gap ) return

    halfway = modulo(atan2(east(before) + east(after),                      &
                           north(before) + north(after)) - bearing(before), &
                     full_turn)
    if ( past(before) < halfway ) then
      source = before
    else
      source = after
    end if

  end function upwind_neighbour

end module mtn_wave_advect_kernel_mod
