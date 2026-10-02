!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some content in this file was generated or refactored with assistance from
! - Claude Code (Claude Opus 5.5), 2026-09-29.

!> @brief In-cloud turbulence potential on pressure levels.

module in_cloud_turb_pot_kernel_mod

  use argument_mod,               only: arg_type,                  &
                                        GH_FIELD, GH_SCALAR,       &
                                        GH_READ, GH_WRITE,         &
                                        GH_REAL, GH_INTEGER,       &
                                        CELL_COLUMN,               &
                                        ANY_DISCONTINUOUS_SPACE_1, &
                                        ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,              only: r_def, i_def, l_def, rmdi
  use driver_water_constants_mod, only: latent_heat_h2o_condensation, &
                                        gas_constant_h2o
  use kernel_mod,                 only: kernel_type
  use science_conversions_mod,    only: zero_degrees_celsius

  implicit none

  private

  !> Kernel metadata for PSyclone
  !> The first field is a model-level field so that nlayers is taken from it.
  type, public, extends(kernel_type) :: in_cloud_turb_pot_kernel_type
    private
    type(arg_type) :: meta_args(10) = (/                                   &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                         &
!        arg_type(GH_SCALAR_ARRAY, GH_REAL, GH_READ, 1), see PSyclone #1312
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                            &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                            &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                            &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                             &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: in_cloud_turb_pot_code
  end type in_cloud_turb_pot_kernel_type

  public :: in_cloud_turb_pot_code

  ! Number of model levels in the spline stencil around each pressure level
  integer(kind=i_def), parameter :: npts = 6_i_def
  ! Lowest model level allowed as the lower bracketing level
  integer(kind=i_def), parameter :: k_lowest = 3_i_def
  ! Saturation vapour pressure over water at zero degrees Celsius (Pa)
  real(kind=r_def), parameter :: es_zero_c = 611.0_r_def
  ! Latent heat of condensation over the gas constant for water vapour (K)
  real(kind=r_def), parameter :: lc_over_rv = latent_heat_h2o_condensation &
                                              / gas_constant_h2o

contains

  !> @brief   Calculates the in-cloud turbulence potential on pressure levels.
  !> @details The saturated equivalent potential temperature, theta_e, is
  !>          calculated on the six model levels around each pressure level.
  !>          Its vertical gradient at the pressure level is found from
  !>          cubic-spline derivatives of theta_e and height with respect
  !>          to pressure. Where either model level bracketing the pressure
  !>          level has cloud and theta_e decreases with height, the output
  !>          is the magnitude of the gradient; otherwise it is zero. Points
  !>          too close to the surface or the model top, or where the
  !>          saturation vapour pressure reaches the pressure on a stencil
  !>          level, are set to rmdi.
  !> @param[in]     nlayers       Number of layers
  !> @param[in]     theta         Potential temperature (K)
  !> @param[in]     exner         Exner pressure on theta levels
  !> @param[in]     height        Height of theta levels (m)
  !> @param[in]     cloud         Bulk cloud fraction
  !> @param[in]     nplev         Number of pressure levels
  !> @param[in]     plevs         Pressure levels (Pa)
  !> @param[in,out] turb_pot      In-cloud turbulence potential (K m-1)
  !> @param[in]     p_zero        Reference pressure (Pa)
  !> @param[in]     kappa         Rd / cp
  !> @param[in]     cp            Specific heat of dry air (J kg-1 K-1)
  !> @param[in]     recip_epsilon Ratio of gas constants, vapour / dry air
  !> @param[in]     ndf_in        Number of degrees of freedom per cell for
  !>                              the theta-level fields
  !> @param[in]     undf_in       Number of total degrees of freedom for the
  !>                              theta-level fields
  !> @param[in]     map_in        Dofmap for the cell at the base of the
  !>                              column for the theta-level fields
  !> @param[in]     ndf_out       Number of degrees of freedom per cell for
  !>                              the pressure-level field
  !> @param[in]     undf_out      Number of total degrees of freedom for the
  !>                              pressure-level field
  !> @param[in]     map_out       Dofmap for the cell at the base of the
  !>                              column for the pressure-level field
  subroutine in_cloud_turb_pot_code(nlayers,                 &
                                    theta,                   &
                                    exner,                   &
                                    height,                  &
                                    cloud,                   &
                                    nplev,                   &
                                    plevs,                   &
                                    turb_pot,                &
                                    p_zero,                  &
                                    kappa,                   &
                                    cp,                      &
                                    recip_epsilon,           &
                                    ndf_in, undf_in, map_in, &
                                    ndf_out, undf_out, map_out)

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers, nplev
    integer(kind=i_def), intent(in) :: ndf_in, undf_in
    integer(kind=i_def), intent(in), dimension(ndf_in)  :: map_in
    integer(kind=i_def), intent(in) :: ndf_out, undf_out
    integer(kind=i_def), intent(in), dimension(ndf_out) :: map_out

    ! Arguments passed explicitly from algorithm
    real(kind=r_def), intent(in),    dimension(undf_in)  :: theta, exner, &
                                                            height, cloud
    real(kind=r_def), intent(inout), dimension(undf_out) :: turb_pot
    real(kind=r_def), intent(in), dimension(nplev) :: plevs
    real(kind=r_def), intent(in) :: p_zero, kappa, cp, recip_epsilon

    ! Internal variables
    integer(kind=i_def) :: k, kp, k_upper, k_lower, i, df
    real(kind=r_def)    :: exner_plev, temperature, es, ws
    real(kind=r_def)    :: dthetae_dp, dz_dp, dthetae_dz
    real(kind=r_def), dimension(npts) :: p_knot, thetae_knot, z_knot
    logical(kind=l_def) :: valid

    do kp = 1, nplev
      df = map_out(1) + kp - 1
      exner_plev = ( plevs(kp) / p_zero )**kappa

      ! Lowest theta level at or above the pressure level
      k_upper = nlayers + 1
      do k = 1, nlayers
        if ( exner(map_in(1) + k) <= exner_plev ) then
          k_upper = k
          exit
        end if
      end do
      k_lower = k_upper - 1

      ! The stencil needs two levels either side of the bracketing pair
      if ( k_lower < k_lowest .or. k_upper > nlayers - 2 ) then
        turb_pot(df) = rmdi
        cycle
      end if

      ! Pressure, height and saturated equivalent potential temperature on
      ! the stencil levels, ordered by increasing pressure
      valid = .true.
      do i = 1, npts
        k = map_in(1) + k_upper + 3 - i
        p_knot(i) = p_zero * exner(k)**(1.0_r_def / kappa)
        z_knot(i) = height(k)
        temperature = theta(k) * exner(k)
        es = es_zero_c * exp( lc_over_rv * ( 1.0_r_def / zero_degrees_celsius &
                                           - 1.0_r_def / temperature ) )
        ! Saturation mixing ratio is undefined where es exceeds pressure
        if ( es >= p_knot(i) ) then
          valid = .false.
          exit
        end if
        ws = es / ( recip_epsilon * ( p_knot(i) - es ) )
        thetae_knot(i) = theta(k) * exp( latent_heat_h2o_condensation * ws &
                                         / ( cp * temperature ) )
      end do

      if ( .not. valid ) then
        turb_pot(df) = rmdi
        cycle
      end if

      dthetae_dp = spline_derivative( p_knot, thetae_knot, plevs(kp) )
      dz_dp = spline_derivative( p_knot, z_knot, plevs(kp) )
      dthetae_dz = dthetae_dp / dz_dp

      ! Turbulent where there is cloud and theta_e decreases with height
      if ( ( cloud(map_in(1) + k_lower) > 0.0_r_def .or.   &
             cloud(map_in(1) + k_upper) > 0.0_r_def ) .and. &
           dthetae_dz < 0.0_r_def ) then
        turb_pot(df) = -dthetae_dz
      else
        turb_pot(df) = 0.0_r_def
      end if

    end do

  end subroutine in_cloud_turb_pot_code

  !> @brief   Derivative of a cubic spline through npts points, evaluated in
  !>          the central interval.
  !> @details The spline's third derivative in each end interval matches that
  !>          of the cubic through the four end points (Forsythe, Malcolm and
  !>          Moler, 1977), so cubic data are fitted exactly.
  !> @param[in] x     Abscissae in increasing order
  !> @param[in] y     Ordinates
  !> @param[in] x_ref Point at which to evaluate the derivative, between
  !>                  x(npts/2) and x(npts/2 + 1)
  !> @return    Derivative dy/dx at x_ref
  function spline_derivative( x, y, x_ref ) result( dydx )

    implicit none

    real(kind=r_def), intent(in) :: x(npts), y(npts), x_ref
    real(kind=r_def) :: dydx

    integer(kind=i_def), parameter :: mid = npts / 2
    integer(kind=i_def) :: i
    real(kind=r_def) :: h(npts - 1), slope(npts - 1)
    real(kind=r_def) :: lower(npts), diag(npts), upper(npts), rhs(npts)
    real(kind=r_def) :: curv(npts)
    real(kind=r_def) :: third_start, third_end, w, t

    do i = 1, npts - 1
      h(i) = x(i + 1) - x(i)
      slope(i) = ( y(i + 1) - y(i) ) / h(i)
    end do

    ! Third divided differences of the first and last four points
    third_start = ( ( slope(3) - slope(2) ) / ( x(4) - x(2) )             &
                  - ( slope(2) - slope(1) ) / ( x(3) - x(1) ) )           &
                  / ( x(4) - x(1) )
    third_end = ( ( slope(npts - 1) - slope(npts - 2) )                   &
                  / ( x(npts) - x(npts - 2) )                             &
                - ( slope(npts - 2) - slope(npts - 3) )                   &
                  / ( x(npts - 1) - x(npts - 3) ) )                       &
                / ( x(npts) - x(npts - 3) )

    ! Tridiagonal system for the second derivatives, curv
    lower(1) = 0.0_r_def
    diag(1)  = -h(1)
    upper(1) = h(1)
    rhs(1)   = 6.0_r_def * h(1)**2 * third_start
    do i = 2, npts - 1
      lower(i) = h(i - 1)
      diag(i)  = 2.0_r_def * ( h(i - 1) + h(i) )
      upper(i) = h(i)
      rhs(i)   = 6.0_r_def * ( slope(i) - slope(i - 1) )
    end do
    lower(npts) = h(npts - 1)
    diag(npts)  = -h(npts - 1)
    upper(npts) = 0.0_r_def
    rhs(npts)   = -6.0_r_def * h(npts - 1)**2 * third_end

    do i = 2, npts
      w = lower(i) / diag(i - 1)
      diag(i) = diag(i) - w * upper(i - 1)
      rhs(i)  = rhs(i) - w * rhs(i - 1)
    end do
    curv(npts) = rhs(npts) / diag(npts)
    do i = npts - 1, 1, -1
      curv(i) = ( rhs(i) - upper(i) * curv(i + 1) ) / diag(i)
    end do

    t = x_ref - x(mid)
    dydx = slope(mid) - h(mid) * ( 2.0_r_def * curv(mid) + curv(mid + 1) ) &
                        / 6.0_r_def                                       &
         + curv(mid) * t                                                  &
         + ( curv(mid + 1) - curv(mid) ) * t**2 / ( 2.0_r_def * h(mid) )

  end function spline_derivative

end module in_cloud_turb_pot_kernel_mod
