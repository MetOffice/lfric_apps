!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some content in this file was generated or refactored with assistance from
! - Met Office Claude Code Enterprise (Claude Opus 5.5).
!> @brief Maximum wind in the column and the pressures of the base and top
!>        of the jet around it.

module max_wind_diags_kernel_mod

  use argument_mod,             only: arg_type,                  &
                                      GH_FIELD, GH_SCALAR,       &
                                      GH_REAL, GH_READ, GH_WRITE, &
                                      CELL_COLUMN,               &
                                      ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,            only: r_def, i_def, l_def, rmdi
  use empty_data_mod,           only: empty_real_data
  use fs_continuity_mod,        only: W3
  use kernel_mod,               only: kernel_type
  use science_conversions_mod,  only: knots_to_metre_per_second

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: max_wind_diags_kernel_type
    private
    ! u_in_w3 comes first so that nlayers is that of W3
    type(arg_type) :: meta_args(10) = (/                                  &
         ! u_in_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3),                        &
         ! max_wind_u
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         ! max_wind_v
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         ! max_wind_pressure
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         ! base_pressure
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         ! top_pressure
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         ! v_in_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3),                        &
         ! exner
         arg_type(GH_FIELD, GH_REAL, GH_READ, W3),                        &
         ! p_zero
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                           &
         ! recip_kappa
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                            &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: max_wind_diags_code
  end type max_wind_diags_kernel_type

  public :: max_wind_diags_code

  ! The splines are fitted through the level of interest and this many
  ! levels either side of it
  integer(kind=i_def), parameter :: half_width = 2
  integer(kind=i_def), parameter :: n_knots = 2 * half_width + 1
  ! Number of evaluation steps between adjacent levels
  integer(kind=i_def), parameter :: n_steps = 16

contains

  !> @brief Finds the maximum wind in the column, and the pressures of the
  !>        base and top of the jet around it.
  !> @details The first guess is the fastest level whose neighbours lie in
  !>          the 100 to 700 hPa window. Cubic splines of u and v against
  !>          pressure refine it between the levels either side. Where the
  !>          refined wind is at least 80 knots, the base and top are where
  !>          the wind falls to 80 knots below and above it. Each output is
  !>          only written if it was requested, i.e. it is not
  !>          empty_real_data. Points without a value are set to rmdi.
  !> @param[in]     nlayers            Number of layers
  !> @param[in]     u_in_w3            Eastward wind (m s-1)
  !> @param[in,out] max_wind_u         Eastward wind at the maximum (m s-1)
  !> @param[in,out] max_wind_v         Northward wind at the maximum (m s-1)
  !> @param[in,out] max_wind_pressure  Pressure at the maximum (Pa)
  !> @param[in,out] base_pressure      Pressure at the jet base (Pa)
  !> @param[in,out] top_pressure       Pressure at the jet top (Pa)
  !> @param[in]     v_in_w3            Northward wind (m s-1)
  !> @param[in]     exner              Exner pressure in W3
  !> @param[in]     p_zero             Reference pressure for Exner (Pa)
  !> @param[in]     recip_kappa        Ratio of cp to R
  !> @param[in]     ndf_w3             Number of DOFs per cell for W3
  !> @param[in]     undf_w3            Number of unique DOFs for W3
  !> @param[in]     map_w3             Dofmap for the W3 column base
  !> @param[in]     ndf_2d             Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d            Number of unique DOFs for 2D fields
  !> @param[in]     map_2d             Dofmap for the 2D fields
  subroutine max_wind_diags_code( nlayers,                    &
                                  u_in_w3,                    &
                                  max_wind_u,                 &
                                  max_wind_v,                 &
                                  max_wind_pressure,          &
                                  base_pressure,              &
                                  top_pressure,               &
                                  v_in_w3,                    &
                                  exner,                      &
                                  p_zero,                     &
                                  recip_kappa,                &
                                  ndf_w3, undf_w3, map_w3,    &
                                  ndf_2d, undf_2d, map_2d )

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in), dimension(ndf_w3) :: map_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d) :: map_2d

    ! Arguments passed explicitly from algorithm
    real(kind=r_def), intent(in), dimension(undf_w3) :: u_in_w3
    real(kind=r_def), intent(inout), pointer :: max_wind_u(:)
    real(kind=r_def), intent(inout), pointer :: max_wind_v(:)
    real(kind=r_def), intent(inout), pointer :: max_wind_pressure(:)
    real(kind=r_def), intent(inout), pointer :: base_pressure(:)
    real(kind=r_def), intent(inout), pointer :: top_pressure(:)
    real(kind=r_def), intent(in), dimension(undf_w3) :: v_in_w3
    real(kind=r_def), intent(in), dimension(undf_w3) :: exner
    real(kind=r_def), intent(in) :: p_zero
    real(kind=r_def), intent(in) :: recip_kappa

    ! The first guess must have the level above at or below 100 hPa and
    ! the level below above 700 hPa (Pa)
    real(kind=r_def), parameter :: p_window_top = 10000.0_r_def
    real(kind=r_def), parameter :: p_window_bottom = 70000.0_r_def
    ! Wind speed at the base and top of the jet (m s-1). The UM converts
    ! knots with 0.5418 m s-1 rather than 0.5144, so its threshold is
    ! 43.3 m s-1 (about 84 knots) where this one is 41.2 m s-1.
    real(kind=r_def), parameter :: jet_speed =                            &
                                   80.0_r_def * knots_to_metre_per_second

    ! Internal variables
    integer(kind=i_def) :: k, k_max, k_below, k_edge, i_step
    logical(kind=l_def) :: l_u, l_v, l_pressure, l_base, l_top
    real(kind=r_def)    :: speed, speed_max, p_max
    real(kind=r_def)    :: u_col(nlayers), v_col(nlayers), p_col(nlayers)
    real(kind=r_def)    :: speed_col(nlayers)
    real(kind=r_def)    :: u_step(2 * n_steps), v_step(2 * n_steps)
    real(kind=r_def)    :: p_step(2 * n_steps)

    l_u        = .not. associated( max_wind_u, empty_real_data )
    l_v        = .not. associated( max_wind_v, empty_real_data )
    l_pressure = .not. associated( max_wind_pressure, empty_real_data )
    l_base     = .not. associated( base_pressure, empty_real_data )
    l_top      = .not. associated( top_pressure, empty_real_data )

    ! Index k of the column arrays is layer k, counting from 1 at the bottom
    do k = 1, nlayers
      u_col(k) = u_in_w3(map_w3(1) + k - 1)
      v_col(k) = v_in_w3(map_w3(1) + k - 1)
      p_col(k) = p_zero * exner(map_w3(1) + k - 1)**recip_kappa
      speed_col(k) = sqrt( u_col(k)**2 + v_col(k)**2 )
    end do

    ! First guess: the fastest level, the highest on a tie. Only levels
    ! with half_width levels either side can carry a spline.
    k_max = 0
    speed_max = 0.0_r_def
    do k = 1 + half_width, nlayers - half_width
      if ( p_col(k + 1) >= p_window_top .and.                             &
           p_col(k - 1) < p_window_bottom .and.                           &
           speed_col(k) >= speed_max ) then
        k_max = k
        speed_max = speed_col(k)
      end if
    end do

    ! The UM also computes a base and top here, from missing data; they
    ! are left missing instead
    if ( k_max == 0 ) then
      if ( l_u ) max_wind_u(map_2d(1)) = rmdi
      if ( l_v ) max_wind_v(map_2d(1)) = rmdi
      if ( l_pressure ) max_wind_pressure(map_2d(1)) = rmdi
      if ( l_base ) base_pressure(map_2d(1)) = rmdi
      if ( l_top ) top_pressure(map_2d(1)) = rmdi
      return
    end if

    ! Refine: the fastest spline step between the levels either side
    call spline_steps( nlayers, k_max, u_col, v_col, p_col,               &
                       u_step, v_step, p_step )
    i_step = maxloc( u_step**2 + v_step**2, dim=1 )

    if ( l_u ) max_wind_u(map_2d(1)) = u_step(i_step)
    if ( l_v ) max_wind_v(map_2d(1)) = v_step(i_step)
    if ( l_pressure ) max_wind_pressure(map_2d(1)) = p_step(i_step)

    if ( .not. ( l_base .or. l_top ) ) return

    p_max = p_step(i_step)
    speed = sqrt( u_step(i_step)**2 + v_step(i_step)**2 )

    if ( speed < jet_speed ) then
      if ( l_base ) base_pressure(map_2d(1)) = rmdi
      if ( l_top ) top_pressure(map_2d(1)) = rmdi
      return
    end if

    ! The highest level below the maximum. The refined maximum lies above
    ! the level below the first guess, so the search always finds one.
    k_below = 0
    do k = nlayers, 1, -1
      if ( p_col(k) > p_max ) then
        k_below = k
        exit
      end if
    end do

    if ( l_base ) then
      k_edge = nearest_jet_level( nlayers, speed_col, jet_speed,          &
                                  k_below, 1 + half_width, -1 )
      base_pressure(map_2d(1)) = edge_pressure( nlayers, k_edge, u_col,   &
                                                v_col, p_col, jet_speed )
    end if

    if ( l_top ) then
      k_edge = nearest_jet_level( nlayers, speed_col, jet_speed,          &
                                  k_below + 1, nlayers - half_width, 1 )
      top_pressure(map_2d(1)) = edge_pressure( nlayers, k_edge, u_col,    &
                                               v_col, p_col, jet_speed )
    end if

  end subroutine max_wind_diags_code

  !> @brief Walks from one level towards another and returns the level whose
  !>        wind speed is nearest the jet speed.
  !> @details The walk ends after the first level slower than the jet
  !>          speed. On a tie the first level visited is kept.
  !> @param[in] nlayers    Number of layers
  !> @param[in] speed_col  Wind speed on each layer (m s-1)
  !> @param[in] jet_speed  Wind speed at the jet edge (m s-1)
  !> @param[in] k_start    First level visited
  !> @param[in] k_end      Last level that may be visited
  !> @param[in] k_step     Direction of the walk, 1 up or -1 down
  !> @return    Level nearest the jet speed, or 0 if no level was visited
  function nearest_jet_level( nlayers, speed_col, jet_speed,              &
                              k_start, k_end, k_step ) result( k_nearest )

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    real(kind=r_def),    intent(in) :: speed_col(nlayers)
    real(kind=r_def),    intent(in) :: jet_speed
    integer(kind=i_def), intent(in) :: k_start, k_end, k_step
    integer(kind=i_def) :: k_nearest

    integer(kind=i_def) :: k
    real(kind=r_def)    :: diff, diff_min

    k_nearest = 0
    diff_min = huge( diff_min )
    do k = k_start, k_end, k_step
      diff = abs( speed_col(k) - jet_speed )
      if ( diff < diff_min ) then
        k_nearest = k
        diff_min = diff
      end if
      if ( speed_col(k) < jet_speed ) exit
    end do

  end function nearest_jet_level

  !> @brief Pressure of the spline step about a level whose wind speed is
  !>        nearest the jet speed.
  !> @param[in] nlayers    Number of layers
  !> @param[in] k_edge     Level the splines are centred on, or 0 for none
  !> @param[in] u_col      Eastward wind on each layer (m s-1)
  !> @param[in] v_col      Northward wind on each layer (m s-1)
  !> @param[in] p_col      Pressure on each layer (Pa)
  !> @param[in] jet_speed  Wind speed at the jet edge (m s-1)
  !> @return    Pressure at the jet edge (Pa), or rmdi if k_edge is 0
  function edge_pressure( nlayers, k_edge, u_col, v_col, p_col,           &
                          jet_speed ) result( p_edge )

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: k_edge
    real(kind=r_def),    intent(in) :: u_col(nlayers), v_col(nlayers)
    real(kind=r_def),    intent(in) :: p_col(nlayers)
    real(kind=r_def),    intent(in) :: jet_speed
    real(kind=r_def) :: p_edge

    real(kind=r_def) :: u_step(2 * n_steps), v_step(2 * n_steps)
    real(kind=r_def) :: p_step(2 * n_steps)

    ! The UM reads outside the column when there is no level to walk;
    ! here the edge is missing
    if ( k_edge == 0 ) then
      p_edge = rmdi
      return
    end if

    call spline_steps( nlayers, k_edge, u_col, v_col, p_col,              &
                       u_step, v_step, p_step )
    p_edge = p_step( minloc( abs( sqrt( u_step**2 + v_step**2 )           &
                                  - jet_speed ), dim=1 ) )

  end function edge_pressure

  !> @brief Evaluates splines of u and v against pressure at evenly spaced
  !>        pressures between the levels either side of a centre level.
  !> @details The first n_steps points run from the level below to the
  !>          centre level and the rest from the centre level to the level
  !>          above, each half excluding its first point and including its
  !>          last.
  !> @param[in]  nlayers   Number of layers
  !> @param[in]  k_centre  Centre level
  !> @param[in]  u_col     Eastward wind on each layer (m s-1)
  !> @param[in]  v_col     Northward wind on each layer (m s-1)
  !> @param[in]  p_col     Pressure on each layer (Pa)
  !> @param[out] u_step    Eastward wind at each step (m s-1)
  !> @param[out] v_step    Northward wind at each step (m s-1)
  !> @param[out] p_step    Pressure at each step (Pa)
  subroutine spline_steps( nlayers, k_centre, u_col, v_col, p_col,        &
                           u_step, v_step, p_step )

    implicit none

    integer(kind=i_def), intent(in)  :: nlayers
    integer(kind=i_def), intent(in)  :: k_centre
    real(kind=r_def),    intent(in)  :: u_col(nlayers), v_col(nlayers)
    real(kind=r_def),    intent(in)  :: p_col(nlayers)
    real(kind=r_def),    intent(out) :: u_step(2 * n_steps)
    real(kind=r_def),    intent(out) :: v_step(2 * n_steps)
    real(kind=r_def),    intent(out) :: p_step(2 * n_steps)

    integer(kind=i_def) :: i, k_lo
    real(kind=r_def)    :: fraction
    real(kind=r_def)    :: p_knot(n_knots), u_knot(n_knots), v_knot(n_knots)
    real(kind=r_def)    :: u_b(n_knots), u_c(n_knots), u_d(n_knots)
    real(kind=r_def)    :: v_b(n_knots), v_c(n_knots), v_d(n_knots)

    ! Knots in order of increasing pressure, i.e. from the top level down
    do i = 1, n_knots
      p_knot(i) = p_col(k_centre + half_width + 1 - i)
      u_knot(i) = u_col(k_centre + half_width + 1 - i)
      v_knot(i) = v_col(k_centre + half_width + 1 - i)
    end do

    call spline_coefficients( p_knot, u_knot, u_b, u_c, u_d )
    call spline_coefficients( p_knot, v_knot, v_b, v_c, v_d )

    do i = 1, 2 * n_steps
      if ( i <= n_steps ) then
        ! Between the level below and the centre level, whose spline
        ! interval starts at the centre knot
        fraction = real( i, r_def ) / real( n_steps, r_def )
        p_step(i) = p_col(k_centre - 1)                                   &
                  + ( p_col(k_centre) - p_col(k_centre - 1) ) * fraction
        k_lo = half_width + 1
      else
        ! Between the centre level and the level above, whose spline
        ! interval starts at the knot above the centre
        fraction = real( i - n_steps, r_def ) / real( n_steps, r_def )
        p_step(i) = p_col(k_centre)                                       &
                  + ( p_col(k_centre + 1) - p_col(k_centre) ) * fraction
        k_lo = half_width
      end if
      u_step(i) = spline_value( p_step(i) - p_knot(k_lo), u_knot(k_lo),   &
                                u_b(k_lo), u_c(k_lo), u_d(k_lo) )
      v_step(i) = spline_value( p_step(i) - p_knot(k_lo), v_knot(k_lo),   &
                                v_b(k_lo), v_c(k_lo), v_d(k_lo) )
    end do

  end subroutine spline_steps

  !> @brief Value of one interval of a cubic spline.
  !> @param[in] dx  Distance from the start of the interval
  !> @param[in] y   Value at the start of the interval
  !> @param[in] b   First-order coefficient
  !> @param[in] c   Second-order coefficient
  !> @param[in] d   Third-order coefficient
  !> @return    Spline value
  function spline_value( dx, y, b, c, d ) result( s )

    implicit none

    real(kind=r_def), intent(in) :: dx, y, b, c, d
    real(kind=r_def) :: s

    s = y + dx * ( b + dx * ( c + dx * d ) )

  end function spline_value

  !> @brief Coefficients of the cubic spline through n_knots points.
  !> @details On interval i the spline is
  !>          y(i) + b(i) dx + c(i) dx**2 + d(i) dx**3, with dx = x - x(i).
  !>          At each end the third derivative equals that of the cubic
  !>          through the four nearest points (Forsythe, Malcolm and Moler,
  !>          1977), so a cubic is reproduced exactly.
  !> @param[in]  x  Abscissas, strictly increasing
  !> @param[in]  y  Ordinates
  !> @param[out] b  First-order coefficients
  !> @param[out] c  Second-order coefficients
  !> @param[out] d  Third-order coefficients
  subroutine spline_coefficients( x, y, b, c, d )

    implicit none

    real(kind=r_def), intent(in)  :: x(n_knots), y(n_knots)
    real(kind=r_def), intent(out) :: b(n_knots), c(n_knots), d(n_knots)

    integer(kind=i_def) :: i
    real(kind=r_def)    :: h(n_knots - 1), slope(n_knots - 1)
    real(kind=r_def)    :: curv2(n_knots - 2), curv3(n_knots - 3)
    real(kind=r_def)    :: diag(n_knots), rhs(n_knots), sigma(n_knots)
    real(kind=r_def)    :: factor

    ! Interval widths and first, second and third divided differences
    do i = 1, n_knots - 1
      h(i) = x(i + 1) - x(i)
      slope(i) = ( y(i + 1) - y(i) ) / h(i)
    end do
    do i = 1, n_knots - 2
      curv2(i) = ( slope(i + 1) - slope(i) ) / ( x(i + 2) - x(i) )
    end do
    do i = 1, n_knots - 3
      curv3(i) = ( curv2(i + 1) - curv2(i) ) / ( x(i + 3) - x(i) )
    end do

    ! Symmetric tridiagonal system for sigma, a sixth of the second
    ! derivative at each knot; the off-diagonal of row i is h(i)
    diag(1) = -h(1)
    rhs(1) = h(1)**2 * curv3(1)
    do i = 2, n_knots - 1
      diag(i) = 2.0_r_def * ( h(i - 1) + h(i) )
      rhs(i) = slope(i) - slope(i - 1)
    end do
    diag(n_knots) = -h(n_knots - 1)
    rhs(n_knots) = -h(n_knots - 1)**2 * curv3(n_knots - 3)

    ! Forward elimination and back substitution
    do i = 2, n_knots
      factor = h(i - 1) / diag(i - 1)
      diag(i) = diag(i) - factor * h(i - 1)
      rhs(i) = rhs(i) - factor * rhs(i - 1)
    end do
    sigma(n_knots) = rhs(n_knots) / diag(n_knots)
    do i = n_knots - 1, 1, -1
      sigma(i) = ( rhs(i) - h(i) * sigma(i + 1) ) / diag(i)
    end do

    do i = 1, n_knots - 1
      b(i) = slope(i) - h(i) * ( sigma(i + 1) + 2.0_r_def * sigma(i) )
      c(i) = 3.0_r_def * sigma(i)
      d(i) = ( sigma(i + 1) - sigma(i) ) / h(i)
    end do
    b(n_knots) = slope(n_knots - 1)                                       &
               + h(n_knots - 1) * ( sigma(n_knots - 1)                    &
                                    + 2.0_r_def * sigma(n_knots) )
    c(n_knots) = 3.0_r_def * sigma(n_knots)
    d(n_knots) = d(n_knots - 1)

  end subroutine spline_coefficients

end module max_wind_diags_kernel_mod
