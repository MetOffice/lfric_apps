!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Galewsky test case geopotential and stream function profile.
!> @details Routine to compute the geopotential and stream function for the
!!          Galewsky instability shallow water test.
!!          Geopotential is set such that it is in balance with a given
!!          velocity. This results in an integral form, which is evaluated
!!          using a quadrature. Also returns a stream function, which can
!!          be used to obtain the velocity field.
!!
module galewsky_test_case_mod

  use constants_mod,          only: r_def, i_def, pi
  use planet_config_mod,      only: scaled_radius, scaled_omega, gravity
  use shallow_water_settings_config_mod, &
                              only: ref_gp

  implicit none

contains

  !> @brief Compute geoptential and stream function for the Galewsky test.
  !> @details Taken from John Thuburn's pdfem code. For a given (lat,long)
  !!          compute Galewsky test case streamfunction (-psi) and geopotential.
  !> @param[in]  lat     Latitude
  !> @param[in]  long    Longitude
  !> @param[out] psi_loc Value of minus psi field at coordinates (lat,long)
  !> @param[out] h_loc   Value of geopotential field at coordinates (lat,long)
  subroutine galewsky_profile(lat, long, psi_loc, h_loc)

    implicit none

    real(kind=r_def), intent(in)  :: lat
    real(kind=r_def), intent(in)  :: long
    real(kind=r_def), intent(out) :: psi_loc
    real(kind=r_def), intent(out) :: h_loc

    real(kind=r_def), allocatable :: psigg(:), hgg(:)
    real(kind=r_def)              :: alpha, beta, clat, dygg
    real(kind=r_def)              :: hpert, piby2, pow

    ! Latitude and longitude values
    real(kind=r_def)              :: l1, l2, lat0, lat1, lat2

    ! Volume and area
    real(kind=r_def)              :: totvol,  totarea

    ! Interpolation coefficients
    real(kind=r_def)              :: cc1, cc2, cc3, cc4

    ! Velocity variables
    real(kind=r_def)              :: u00, en, umen, den, &
                                     uu1, uu2, e1, e2

    ! Indices
    integer(kind=i_def)           :: nygg, j, jy

    ! ncells is used to determine nr of quadrature points
    integer(kind=i_def), parameter :: ncell = 6400

    ! Allocate array size
    nygg = 2.0_r_def*floor(sqrt(real(ncell)))
    allocate(hgg(nygg+1), psigg(nygg+1))

    ! ----------------------------------------------------------------
    ! Galewsky et al. (2004) test case

    ! Set constants
    piby2 = 0.5_r_def*pi
    u00   = 80.0_r_def
    lat0  = pi/7.0_r_def
    lat1  = pi/2.0_r_def - lat0
    en    = exp(-4.0_r_def/(lat1 - lat0)**2)
    umen  = u00/en
    dygg  = pi/nygg

    ! Initialise values before integrating
    totvol   = 0.0_r_def
    totarea  = 0.0_r_def
    hgg(1)   = 0.0_r_def
    psigg(1) = 0.0_r_def
    ! Integrate to tabulate h and psi as functions of geographical latitude
    do j = 2, nygg
      ! Loop over quad points setting latitude values l1 and l2
      l1 = (j-2_r_def)*dygg - piby2
      ! Set denominator in eq(2) of Galewsky et al.
      den = (l1 - lat0)*(l1 - lat1)
      ! Set u(lat) as per eq(2) of Galewsky et al.
      if (den .lt. 0.0_r_def) then
        uu1 = umen*exp(1.0_r_def/den)
      else
        uu1 = 0.0_r_def
      endif
      l2 = (j-1)*dygg - piby2
      ! Set denominator in eq(2) of Galewsky et al.
      den = (l2 - lat0)*(l2 - lat1)
      ! Set u(lat) as per eq(2) of Galewsky et al.
      if (den .lt. 0.0_r_def) then
        uu2 = umen*exp(1.0_r_def/den)
      else
        uu2 = 0.0_r_def
      endif
      ! Sum psi over quad points
      psigg(j) = psigg(j-1) - 0.5_r_def*(uu1 + uu2)*dygg
      ! Add metric terms to u
      uu1 = uu1*(2.0_r_def*scaled_omega*sin(l1) + tan(l1)*uu1/scaled_radius)
      uu2 = uu2*(2.0_r_def*scaled_omega*sin(l2) + tan(l2)*uu2/scaled_radius)
      ! Sum h over quad points
      hgg(j) = hgg(j-1) - scaled_radius*0.5_r_def*(uu1 + uu2)*dygg
      ! Compute total area and volume
      totarea = totarea + cos(l2)*dygg
      totvol = totvol + hgg(j)*cos(l2)*dygg
    enddo
    psigg(nygg+1) = psigg(nygg)
    hgg(nygg+1) = hgg(nygg)
    totvol = totvol/totarea
    hgg = hgg + (ref_gp - totvol)

    ! Compute h value at (long,lat) using interpolation from tabulated values

    l1 = lat
    l2 = long

    ! Initialize at computational point...
    l1 = l1 + piby2
    jy = floor(l1/dygg) + 1
    if ( jy > nygg ) jy = nygg
    beta = (l1 - (jy - 1_r_def)*dygg)/dygg
    if (jy == 1 .or. jy == nygg) then
    ! linear interpolation
      cc2 = 1.0_r_def - beta
      cc3 = beta
      h_loc = (cc2*hgg(jy) + cc3*hgg(jy+1))
    else
      ! cubic interpolation
      cc1 = -beta*(beta - 1.0_r_def)*(beta - 2.0_r_def)/6.0_r_def
      cc2 = 0.5_r_def*(beta + 1.0_r_def)*(beta - 1.0_r_def)*(beta - 2.0_r_def)
      cc3 = -0.5_r_def*(beta + 1.0_r_def)*beta*(beta - 2.0_r_def)
      cc4 = (beta + 1.0_r_def)*beta*(beta - 1.0_r_def)/6.0_r_def
      h_loc = (cc1*hgg(jy-1) + cc2*hgg(jy) + cc3*hgg(jy+1) + cc4*hgg(jy+2))
    end if

    ! Compute psi value at (long,lat) using interpolation from tabulated values

    l1 = lat
    l2 = long

    l1 = l1 + piby2
    jy = floor(l1/dygg) + 1
    if ( jy > nygg ) jy = nygg
    beta = (l1 - (jy - 1_r_def)*dygg)/dygg
    if (jy == 1 .or. jy == nygg) then
      ! linear interpolation
      cc2 = 1.0_r_def - beta
      cc3 = beta
      psi_loc = cc2*psigg(jy) + cc3*psigg(jy+1)
    else
      ! cubic interpolation
      cc1 = -beta*(beta - 1.0_r_def)*(beta - 2.0_r_def)/6.0_r_def
      cc2 = 0.5_r_def*(beta + 1.0_r_def)*(beta - 1.0_r_def)*(beta - 2.0_r_def)
      cc3 = -0.5_r_def*(beta + 1.0_r_def)*beta*(beta - 2.0_r_def)
      cc4 = (beta + 1.0_r_def)*beta*(beta - 1.0_r_def)/6.0_r_def
      psi_loc = cc1*psigg(jy-1) + cc2*psigg(jy) + cc3*psigg(jy+1) + cc4*psigg(jy+2)
    endif

    ! Change sign to return minus stream function
    psi_loc = -psi_loc*scaled_radius

    deallocate(hgg, psigg)

    ! Set geopotential perturbation

    ! Constants
    alpha = 1.0_r_def/3.0_r_def
    beta  = 1.0_r_def/15.0_r_def
    hpert = 120.0_r_def
    lat2  = 0.5_r_def*piby2

    ! Latitude and longitude
    l2   = lat
    clat = cos(l2)
    l1   = long
    if (l1 > pi) l1 = l1 - 2.0_r_def*pi

    ! Set perturbation as per eq(4) of Galewsky et al.
    e1 = exp(-(l1/alpha)**2)
    pow = max(-50.0_r_def,-((lat2 - l2)/beta)**2)
    e2 = exp(pow)
    h_loc = h_loc + gravity*hpert*clat*e1*e2

  ! ----------------------------------------------------------------

  end subroutine galewsky_profile

end module galewsky_test_case_mod
