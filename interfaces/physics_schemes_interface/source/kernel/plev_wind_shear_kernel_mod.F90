!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some content in this file was generated or refactored with assistance from
! - Claude Code (Claude Opus 5.5), 2026-10-20.
!-----------------------------------------------------------------------------
!> @brief Winds and vertical wind shear on one pressure level.

module plev_wind_shear_kernel_mod

  use argument_mod,  only: arg_type,                  &
                           GH_FIELD, GH_SCALAR,       &
                           GH_READ, GH_READWRITE,     &
                           GH_REAL, GH_INTEGER,       &
                           CELL_COLUMN,               &
                           ANY_DISCONTINUOUS_SPACE_1
  use constants_mod, only: r_def, i_def, rmdi
  use fs_continuity_mod, only: W3
  use kernel_mod,    only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: plev_wind_shear_kernel_type
    private
    type(arg_type) :: meta_args(11) = (/                                       &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                    &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                    &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                    &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                    &
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE,                         &
                                         ANY_DISCONTINUOUS_SPACE_1),           &
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE,                         &
                                         ANY_DISCONTINUOUS_SPACE_1),           &
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE,                         &
                                         ANY_DISCONTINUOUS_SPACE_1),           &
         arg_type(GH_SCALAR, GH_REAL,    GH_READ),                             &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                             &
         arg_type(GH_SCALAR, GH_REAL,    GH_READ),                             &
         arg_type(GH_SCALAR, GH_REAL,    GH_READ)                              &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: plev_wind_shear_code
  end type plev_wind_shear_kernel_type

  public :: plev_wind_shear_code

  !> Number of model levels used by the cubic spline
  integer(kind=i_def), parameter :: nknots = 6_i_def

contains

  !> @brief   Winds and vertical wind shear on one pressure level.
  !> @details The horizontal wind components are interpolated to the
  !>          pressure level linearly in Exner pressure, holding the end
  !>          values outside the column (as pres_interp_kernel_mod). The
  !>          vertical shear of the horizontal wind is
  !>          |dV/dp| / |dz/dp|, with each derivative taken from a cubic
  !>          spline through the six model levels around the pressure
  !>          level and evaluated at that pressure. Where the six levels
  !>          do not fit in the column the shear is set to missing data.
  !> @param[in]     nlayers    Number of layers
  !> @param[in]     u_in_w3    Eastward wind on model levels (m s-1)
  !> @param[in]     v_in_w3    Northward wind on model levels (m s-1)
  !> @param[in]     exner_w3   Exner pressure on model levels
  !> @param[in]     height_w3  Height of model levels (m)
  !> @param[in,out] u_plev     Eastward wind on pressure levels (m s-1)
  !> @param[in,out] v_plev     Northward wind on pressure levels (m s-1)
  !> @param[in,out] wind_shear Vertical wind shear on pressure levels (s-1)
  !> @param[in]     press_lev  Pressure of this level (Pa)
  !> @param[in]     kp         Index of this level in the output fields
  !> @param[in]     p_zero     Reference pressure (Pa)
  !> @param[in]     kappa      Ratio of gas constant to specific heat
  !> @param[in]     ndf_w3     Number of degrees of freedom per cell for W3
  !> @param[in]     undf_w3    Number of unique degrees of freedom for W3
  !> @param[in]     map_w3     Dofmap for the cell at the base of the column
  !> @param[in]     ndf_pl     Number of degrees of freedom per cell for the
  !>                           pressure-level fields
  !> @param[in]     undf_pl    Number of unique degrees of freedom for the
  !>                           pressure-level fields
  !> @param[in]     map_pl     Dofmap for the pressure-level fields
  subroutine plev_wind_shear_code(nlayers,                       &
                                  u_in_w3, v_in_w3,              &
                                  exner_w3, height_w3,           &
                                  u_plev, v_plev, wind_shear,    &
                                  press_lev, kp, p_zero, kappa,  &
                                  ndf_w3, undf_w3, map_w3,       &
                                  ndf_pl, undf_pl, map_pl)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in), dimension(ndf_w3) :: map_w3
    integer(kind=i_def), intent(in) :: ndf_pl, undf_pl
    integer(kind=i_def), intent(in), dimension(ndf_pl) :: map_pl

    ! Arguments passed explicitly from algorithm
    real(kind=r_def), intent(in),    dimension(undf_w3) :: u_in_w3
    real(kind=r_def), intent(in),    dimension(undf_w3) :: v_in_w3
    real(kind=r_def), intent(in),    dimension(undf_w3) :: exner_w3
    real(kind=r_def), intent(in),    dimension(undf_w3) :: height_w3
    real(kind=r_def), intent(inout), dimension(undf_pl) :: u_plev
    real(kind=r_def), intent(inout), dimension(undf_pl) :: v_plev
    real(kind=r_def), intent(inout), dimension(undf_pl) :: wind_shear
    real(kind=r_def), intent(in)    :: press_lev
    integer(kind=i_def), intent(in) :: kp
    real(kind=r_def), intent(in)    :: p_zero, kappa

    ! Internal variables
    integer(kind=i_def) :: k, k_lower, slot
    real(kind=r_def)    :: exner_lev, dudp, dvdp, dzdp
    real(kind=r_def)    :: pressure(nlayers)
    real(kind=r_def)    :: p_knots(nknots)

    slot = map_pl(1) + kp - 1
    exner_lev = ( press_lev / p_zero )**kappa

    ! Winds on the pressure level, linear in Exner pressure
    u_plev(slot) = interp_in_exner(nlayers, u_in_w3(map_w3(1):),            &
                                   exner_w3(map_w3(1):), exner_lev)
    v_plev(slot) = interp_in_exner(nlayers, v_in_w3(map_w3(1):),            &
                                   exner_w3(map_w3(1):), exner_lev)

    do k = 1, nlayers
      pressure(k) = p_zero * exner_w3(map_w3(1) + k - 1)**(1.0_r_def / kappa)
    end do

    ! Model level k_lower is the level below the pressure level:
    ! pressure(k_lower) > press_lev >= pressure(k_lower + 1)
    k_lower = 0_i_def
    do k = 1, nlayers - 1
      if ( pressure(k) > press_lev .and. pressure(k + 1) <= press_lev ) then
        k_lower = k
        exit
      end if
    end do

    ! The six spline knots are levels k_lower-2 to k_lower+3, so the
    ! pressure level lies in the middle interval
    if ( k_lower < 3_i_def .or. k_lower + 3_i_def > nlayers ) then
      wind_shear(slot) = rmdi
    else
      ! Knots in order of increasing pressure
      do k = 1, nknots
        p_knots(k) = pressure(k_lower + 4 - k)
      end do

      dzdp = spline_slope(p_knots,                                         &
                          knot_values(height_w3, map_w3(1), k_lower),      &
                          press_lev)
      dudp = spline_slope(p_knots,                                         &
                          knot_values(u_in_w3, map_w3(1), k_lower),        &
                          press_lev)
      dvdp = spline_slope(p_knots,                                         &
                          knot_values(v_in_w3, map_w3(1), k_lower),        &
                          press_lev)

      wind_shear(slot) = sqrt( dudp**2 + dvdp**2 ) / abs(dzdp)
    end if

  end subroutine plev_wind_shear_code

  !> @brief   Interpolate a column linearly in Exner pressure.
  !> @details Values outside the column take the value at the nearest end.
  !> @param[in] nlayers   Number of layers
  !> @param[in] values    Column values from the lowest level upwards
  !> @param[in] exner     Exner pressure of the column levels
  !> @param[in] exner_lev Exner pressure to interpolate to
  !> @return    Interpolated value
  function interp_in_exner(nlayers, values, exner, exner_lev) result(value)

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    real(kind=r_def),    intent(in) :: values(:), exner(:)
    real(kind=r_def),    intent(in) :: exner_lev
    real(kind=r_def) :: value

    integer(kind=i_def) :: k, k_above

    ! First level with Exner pressure below the target
    k_above = 0_i_def
    do k = 1, nlayers
      if ( exner(k) < exner_lev ) then
        k_above = k
        exit
      end if
    end do

    if ( k_above == 0_i_def ) then
      ! Above the model top
      value = values(nlayers)
    else if ( k_above == 1_i_def ) then
      ! Below the lowest level
      value = values(1)
    else
      value = ( (exner_lev - exner(k_above - 1)) * values(k_above)          &
              - (exner_lev - exner(k_above)) * values(k_above - 1) )        &
              / ( exner(k_above) - exner(k_above - 1) )
    end if

  end function interp_in_exner

  !> @brief   Column values at the spline knots, in order of increasing
  !>          pressure.
  !> @param[in] field   Field on model levels
  !> @param[in] base    Dof of the lowest level of the column
  !> @param[in] k_lower Model level just below the pressure level
  !> @return    Values at the six knots
  function knot_values(field, base, k_lower) result(values)

    implicit none

    real(kind=r_def),    intent(in) :: field(:)
    integer(kind=i_def), intent(in) :: base, k_lower
    real(kind=r_def) :: values(nknots)

    integer(kind=i_def) :: k

    do k = 1, nknots
      values(k) = field(base + k_lower + 3 - k)
    end do

  end function knot_values

  !> @brief   Slope of a cubic spline through six knots, in the middle
  !>          interval.
  !> @details Cubic spline with the end conditions of Forsythe, Malcolm and
  !>          Moler (1977): the third derivative at each end matches that of
  !>          the cubic through the four end knots. The tridiagonal system
  !>          for the second derivatives is solved by forward elimination
  !>          and back substitution, and the derivative is evaluated at
  !>          point xp, which lies between knots 3 and 4.
  !> @param[in] x  Knot positions, strictly increasing
  !> @param[in] y  Values at the knots
  !> @param[in] xp Position at which the slope is wanted
  !> @return    Slope dy/dx at xp
  function spline_slope(x, y, xp) result(slope)

    implicit none

    real(kind=r_def), intent(in) :: x(nknots), y(nknots)
    real(kind=r_def), intent(in) :: xp
    real(kind=r_def) :: slope

    integer(kind=i_def), parameter :: mid = nknots / 2_i_def

    ! Knot spacing and secant slope of each interval
    real(kind=r_def) :: h(nknots - 1), secant(nknots - 1)
    ! Diagonal, right-hand side and solution (second derivative / 6)
    real(kind=r_def) :: diag(nknots), rhs(nknots), sigma(nknots)
    real(kind=r_def) :: factor, dx
    integer(kind=i_def) :: i

    do i = 1, nknots - 1
      h(i) = x(i + 1) - x(i)
      secant(i) = ( y(i + 1) - y(i) ) / h(i)
    end do

    ! Interior rows: h(i-1), 2(h(i-1)+h(i)), h(i)
    do i = 2, nknots - 1
      diag(i) = 2.0_r_def * ( h(i - 1) + h(i) )
      rhs(i)  = secant(i) - secant(i - 1)
    end do

    ! End rows from the third divided differences of the end knots
    diag(1) = -h(1)
    rhs(1)  = ( rhs(3) / ( x(4) - x(2) ) - rhs(2) / ( x(3) - x(1) ) )     &
              * h(1)**2 / ( x(4) - x(1) )
    diag(nknots) = -h(nknots - 1)
    rhs(nknots)  = -( rhs(nknots - 1) / ( x(nknots) - x(nknots - 2) )     &
                    - rhs(nknots - 2) / ( x(nknots - 1) - x(nknots - 3) ) ) &
                   * h(nknots - 1)**2 / ( x(nknots) - x(nknots - 3) )

    ! Forward elimination; the off-diagonals are h(i-1) and h(i)
    do i = 2, nknots
      factor  = h(i - 1) / diag(i - 1)
      diag(i) = diag(i) - factor * h(i - 1)
      rhs(i)  = rhs(i) - factor * rhs(i - 1)
    end do

    ! Back substitution
    sigma(nknots) = rhs(nknots) / diag(nknots)
    do i = nknots - 1, mid, -1
      sigma(i) = ( rhs(i) - h(i) * sigma(i + 1) ) / diag(i)
    end do

    ! Derivative of the cubic on the middle interval
    dx = xp - x(mid)
    slope = secant(mid) - h(mid) * ( sigma(mid + 1) + 2.0_r_def * sigma(mid) ) &
            + dx * ( 6.0_r_def * sigma(mid)                                    &
            + 3.0_r_def * dx * ( sigma(mid + 1) - sigma(mid) ) / h(mid) )

  end function spline_slope

end module plev_wind_shear_kernel_mod
