!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Column inputs for the clear air turbulence predictor at one
!>        pressure level: winds on the level and the vertical wind shear.

module cat_column_kernel_mod

  use argument_mod,      only: arg_type,                  &
                               GH_FIELD, GH_SCALAR,       &
                               GH_REAL, GH_INTEGER,       &
                               GH_READ, GH_READWRITE,     &
                               CELL_COLUMN,               &
                               ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,     only: r_def, i_def, rmdi
  use fs_continuity_mod, only: W3
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  !>
  !> The W3 inputs are listed before the outputs, against the usual
  !> outputs-first convention, because PSyclone takes nlayers from the first
  !> field argument and the outputs are single-layer fields on the 2D mesh.
  type, public, extends(kernel_type) :: cat_column_kernel_type
    private
    type(arg_type) :: meta_args(11) = (/                          &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ, W3),            & ! u
         arg_type(GH_FIELD,  GH_REAL,    GH_READ, W3),            & ! v
         arg_type(GH_FIELD,  GH_REAL,    GH_READ, W3),            & ! exner
         arg_type(GH_FIELD,  GH_REAL,    GH_READ, W3),            & ! height
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE,            &
                  ANY_DISCONTINUOUS_SPACE_1),                    & ! u_plev
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE,            &
                  ANY_DISCONTINUOUS_SPACE_1),                    & ! v_plev
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE,            &
                  ANY_DISCONTINUOUS_SPACE_1),                    & ! dwdz
         arg_type(GH_SCALAR, GH_REAL,    GH_READ),                & ! press
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                & ! kp
         arg_type(GH_SCALAR, GH_REAL,    GH_READ),                & ! p_zero
         arg_type(GH_SCALAR, GH_REAL,    GH_READ)                 & ! kappa
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: cat_column_code
  end type cat_column_kernel_type

  public :: cat_column_code

  !> Number of model levels used by the spline derivative (UM DiffP NPts)
  integer(kind=i_def), parameter :: npts = 6_i_def
  !> Half of npts; the spline is evaluated on segment halfn (UM DiffP halfn)
  integer(kind=i_def), parameter :: halfn = npts / 2_i_def

contains

  !> @brief   Winds on a pressure level and the vertical wind shear there.
  !> @details Interpolates the cell-centre winds to the requested pressure
  !>          level (linear in Exner pressure, the same interpolation as
  !>          pres_interp_kernel_mod uses for every pressure-level
  !>          diagnostic) and computes the vertical derivative of wind speed
  !>          with respect to height at that level, using the derivatives of
  !>          height and of both wind components with respect to pressure
  !>          from a six-point cubic spline through the model levels. The
  !>          spline derivative is a direct port of the UM PWS routine
  !>          DiffP (diff_mod), and the shear is that of the UM routine
  !>          pws_cat (pws_cat_mod). Results are written into slot kp of the
  !>          pressure-level output fields; other slots are left unchanged.
  !>          Columns with fewer than six layers cannot support the spline
  !>          and receive missing data.
  !> @param[in]     nlayers    Number of layers
  !> @param[in]     u_in_w3    Zonal wind at cell centres
  !> @param[in]     v_in_w3    Meridional wind at cell centres
  !> @param[in]     exner_w3   Exner pressure at cell centres
  !> @param[in]     height_w3  Height of cell centres above mean sea level
  !> @param[in,out] u_plev     Zonal wind on pressure levels
  !> @param[in,out] v_plev     Meridional wind on pressure levels
  !> @param[in,out] dwdz       Vertical derivative of wind speed with respect
  !>                           to height on pressure levels
  !> @param[in]     press_lev  Pressure of the level to fill (Pa)
  !> @param[in]     kp         Slot of the pressure-level fields to fill
  !> @param[in]     p_zero     Reference surface pressure
  !> @param[in]     kappa      Rd / cp
  !> @param[in]     ndf_w3     Number of degrees of freedom per cell for W3
  !> @param[in]     undf_w3    Number of unique degrees of freedom for W3
  !> @param[in]     map_w3     Dofmap for the cell at the base of the column
  !>                           for W3
  !> @param[in]     ndf_2d     Number of degrees of freedom per cell for the
  !>                           pressure-level fields
  !> @param[in]     undf_2d    Number of unique degrees of freedom for the
  !>                           pressure-level fields
  !> @param[in]     map_2d     Dofmap for the cell for the pressure-level
  !>                           fields
  subroutine cat_column_code(nlayers,                       &
                             u_in_w3, v_in_w3,              &
                             exner_w3, height_w3,           &
                             u_plev, v_plev, dwdz,          &
                             press_lev, kp,                 &
                             p_zero, kappa,                 &
                             ndf_w3, undf_w3, map_w3,       &
                             ndf_2d, undf_2d, map_2d)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in), dimension(ndf_w3) :: map_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d) :: map_2d

    real(kind=r_def), intent(in), dimension(undf_w3) :: u_in_w3
    real(kind=r_def), intent(in), dimension(undf_w3) :: v_in_w3
    real(kind=r_def), intent(in), dimension(undf_w3) :: exner_w3
    real(kind=r_def), intent(in), dimension(undf_w3) :: height_w3
    real(kind=r_def), intent(inout), dimension(undf_2d) :: u_plev
    real(kind=r_def), intent(inout), dimension(undf_2d) :: v_plev
    real(kind=r_def), intent(inout), dimension(undf_2d) :: dwdz

    real(kind=r_def),    intent(in) :: press_lev
    integer(kind=i_def), intent(in) :: kp
    real(kind=r_def),    intent(in) :: p_zero, kappa

    ! Internal variables
    integer(kind=i_def) :: k, slot
    real(kind=r_def) :: desired_ex
    real(kind=r_def), dimension(nlayers) :: p_col, u_col, v_col, z_col
    real(kind=r_def), dimension(nlayers) :: ex_col
    real(kind=r_def) :: dudp, dvdp, dzdp

    slot = map_2d(1) + kp - 1_i_def

    ! Column profiles, level 1 at the bottom, pressure decreasing upwards
    do k = 1, nlayers
      ex_col(k) = exner_w3(map_w3(1) + k - 1)
      p_col(k)  = p_zero * ex_col(k)**(1.0_r_def / kappa)
      u_col(k)  = u_in_w3(map_w3(1) + k - 1)
      v_col(k)  = v_in_w3(map_w3(1) + k - 1)
      z_col(k)  = height_w3(map_w3(1) + k - 1)
    end do

    ! Winds on the pressure level
    desired_ex = ( press_lev / p_zero )**kappa
    u_plev(slot) = interp_to_exner(nlayers, u_col, ex_col, desired_ex)
    v_plev(slot) = interp_to_exner(nlayers, v_col, ex_col, desired_ex)

    if (nlayers < npts) then
      dwdz(slot) = rmdi
    else

      ! Vertical derivatives of height and wind components with respect to
      ! pressure (UM pws_cat lines 186-199)
      call diffp(nlayers, press_lev, z_col, p_col, dzdp)
      call diffp(nlayers, press_lev, u_col, p_col, dudp)
      call diffp(nlayers, press_lev, v_col, p_col, dvdp)

      ! Vertical derivative of wind speed with respect to height
      ! (UM pws_cat line 205)
      dwdz(slot) = sqrt(dudp**2 + dvdp**2) / abs(dzdp)
    end if

  end subroutine cat_column_code

  !> @brief   Interpolate a column profile to a given Exner pressure.
  !> @details Linear interpolation in Exner pressure between the bracketing
  !>          levels, holding the top or bottom level value beyond the
  !>          model top or below the lowest level. This is the same
  !>          interpolation as pres_interp_code (pres_interp_kernel_mod),
  !>          restated here for one level so the kernel stays
  !>          self-contained.
  !> @param[in] numlevs     Number of levels in the profile
  !> @param[in] ffield      Profile to interpolate
  !> @param[in] exner       Exner pressure at each level, decreasing with
  !>                        level
  !> @param[in] desired_ex  Exner pressure of the required level
  !> @return    Interpolated value
  function interp_to_exner(numlevs, ffield, exner, desired_ex) result(fout)

    implicit none

    integer(kind=i_def), intent(in) :: numlevs
    real(kind=r_def),    intent(in), dimension(numlevs) :: ffield
    real(kind=r_def),    intent(in), dimension(numlevs) :: exner
    real(kind=r_def),    intent(in) :: desired_ex
    real(kind=r_def) :: fout

    integer(kind=i_def) :: k, level_above

    ! First level whose pressure is below the desired pressure
    level_above = 0_i_def
    do k = 1, numlevs
      if ( exner(k) < desired_ex ) then
        level_above = k
        exit
      end if
    end do

    if ( level_above == 0_i_def ) then
      ! Desired level is above the model top: use the highest level value
      fout = ffield(numlevs)
    else if ( level_above == 1_i_def ) then
      ! Desired level is below the lowest level: use its value
      fout = ffield(1)
    else
      fout = ( ( desired_ex - exner(level_above - 1) ) * ffield(level_above)   &
             - ( desired_ex - exner(level_above) ) * ffield(level_above - 1) ) &
             / ( exner(level_above) - exner(level_above - 1) )
    end if

  end function interp_to_exner

  !> @brief   Derivative of a column profile with respect to pressure at a
  !>          given pressure, by cubic spline.
  !> @details Direct port of the UM PWS routine DiffP (diff_mod). A bisection
  !>          over levels halfn to numlevs+1-halfn finds the pair of levels
  !>          bracketing pref, so that the six levels jupr-3 to jupr+2 always
  !>          exist; a pref outside that window is extrapolated by the end
  !>          spline segment, as in the UM. The spline coefficients are
  !>          computed with the pressure axis reversed so that the knots
  !>          increase, and only the coefficients of segment halfn, on which
  !>          the derivative is evaluated, are completed.
  !> @param[in]  numlevs  Number of levels in the profile
  !> @param[in]  pref     Pressure at which the derivative is required
  !> @param[in]  ffield   Profile to differentiate
  !> @param[in]  pfield   Pressure at each level, decreasing with level
  !> @param[out] dfdp     Derivative of ffield with respect to pressure
  subroutine diffp(numlevs, pref, ffield, pfield, dfdp)

    implicit none

    integer(kind=i_def), intent(in) :: numlevs
    real(kind=r_def),    intent(in) :: pref
    real(kind=r_def),    intent(in), dimension(numlevs) :: ffield
    real(kind=r_def),    intent(in), dimension(numlevs) :: pfield
    real(kind=r_def),    intent(out) :: dfdp

    integer(kind=i_def) :: jlwr, jupr, jmid, k
    real(kind=r_def) :: dx, t
    real(kind=r_def), dimension(npts)     :: b, c, d
    real(kind=r_def), dimension(npts - 1) :: deltay, deltax, dydx

    ! Divide and conquer to find the level containing pref
    jlwr = halfn
    jupr = numlevs + 1_i_def - halfn
    do while ( (jupr - jlwr) > 1_i_def )
      jmid = (jlwr + jupr) / 2_i_def
      if ( pref >= pfield(jmid) ) then
        jupr = jmid
      end if
      if ( pref < pfield(jmid) ) then
        jlwr = jmid
      end if
    end do

    ! Use fields in reverse order so that pressure is increasing
    dx = pref - pfield(jupr)
    do k = 1, npts - 1
      deltax(k) = pfield(jupr + halfn - k) - pfield(jupr + halfn - 1 - k)
      deltay(k) = ffield(jupr + halfn - k) - ffield(jupr + halfn - 1 - k)
    end do

    ! Cubic interpolating spline
    !   s(u) = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
    ! Tridiagonal system: b = diagonal, d = offdiagonal, c = right-hand side.
    ! Third derivatives at x(1) and x(n) obtained from divided differences.
    dydx(1:npts-1) = deltay(1:npts-1) / deltax(1:npts-1)
    c(2:npts-1) = dydx(2:npts-1) - dydx(1:npts-2)
    b(1)        = -deltax(1)
    b(npts)     = -deltax(npts-1)
    b(2:npts-1) = deltax(1:npts-2) + deltax(2:npts-1)

    c(1)    = c(3) / b(3) - c(2) / b(2)
    c(npts) = c(npts-2) / b(npts-2) - c(npts-1) / b(npts-1)

    c(1)    = c(1) * b(1)**2 / ( b(3) - b(1) )
    c(npts) = c(npts) * b(npts)**2 / ( b(npts-2) - b(npts) )

    b(2:npts-1) = 2.0_r_def * b(2:npts-1)

    ! Forward elimination
    do k = 1, npts - 1
      t = deltax(k) / b(k)
      b(k+1) = b(k+1) - t * deltax(k)
      c(k+1) = c(k+1) - t * c(k)
    end do
    c(npts) = c(npts) / b(npts)

    ! Back substitution, only as far as the segment we need
    do k = npts - 1, halfn, -1
      c(k) = ( c(k) - deltax(k) * c(k+1) ) / b(k)
    end do

    b(halfn) = dydx(halfn)                                                    &
             - deltax(halfn) * ( c(halfn+1) + 2.0_r_def * c(halfn) )
    d(halfn) = ( c(halfn+1) - c(halfn) ) / deltax(halfn)
    c(halfn) = 3.0_r_def * c(halfn)

    ! Derivative of the spline at pref
    dfdp = b(halfn)                                                           &
         + dx * ( 2.0_r_def * c(halfn) + 3.0_r_def * dx * d(halfn) )

  end subroutine diffp

end module cat_column_kernel_mod
