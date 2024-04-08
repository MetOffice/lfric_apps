!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Helper functions for analytic orography.
!-------------------------------------------------------------------------------
module orography_helper_functions_mod

  use constants_mod, only : r_def, i_def
  use domain_mod,    only : domain_type
  use log_mod,       only : log_event,         &
                            log_scratch_space, &
                            LOG_LEVEL_INFO,    &
                            LOG_LEVEL_DEBUG

  implicit none

  private

  !> Domain length and domain width (horizontal).
  real(kind=r_def) :: domain_length, domain_width

  public :: eta2z_linear
  public :: eta2z_smooth
  public :: z2eta_linear
  public :: set_horizontal_domain_size
  public :: coord_transform_cart_biperiodic

contains
  !=============================================================================
  !> @brief Transforms nondimensional coordinate to physical height using
  !>        'smooth' method used in the UM
  !>
  !> @details Helper routine which calculates physical height from eta
  !>          (nondimensional) coordinate using 'smooth' method of the UM
  !>
  !> @param[in] eta            Nondimensional (terrain following) coordinate
  !> @param[in] surface_height Surface height (m)
  !> @param[in] domain_top     Height of the domain (m)
  !> @param[in] eta_c          Physical height at the interface to constant layers
  !> @return    z              Physical height (m)
  !=============================================================================
  real(kind=r_def) function eta2z_smooth(eta,                &
                                         surface_height,     &
                                         domain_top,         &
                                         stretching_height ) &
                                         result(z)

    implicit none

    ! Arguments
    real(kind=r_def), intent(in) :: eta
    real(kind=r_def), intent(in) :: surface_height
    real(kind=r_def), intent(in) :: domain_top
    real(kind=r_def), intent(in) :: stretching_height

    ! Value of eta at stretching_height
    real(kind=r_def) :: eta_c

    eta_c = stretching_height/domain_top

    ! If physical height is above the top of the domain, then set eta_c=1
    if ( eta_c > 1.0_r_def ) eta_c=1.0_r_def

    ! Calculate physical height from eta
    if ( eta < eta_c )then
      ! Quadratic
      z = eta*domain_top + (1.0_r_def - eta/eta_c)**2*surface_height
    else
      ! Linear
      z = eta*domain_top
    end if

    return
  end function eta2z_smooth

  !=============================================================================
  !> @brief Transforms nondimensional coordinate to physical height.
  !>
  !> @details Helper routine which calculates physical height from eta
  !>          (nondimensional) coordinate using linear transformation from
  !>          Wood et al. (2013), Section 7.
  !>
  !> @param[in] eta            Nondimensional (terrain following) coordinate
  !> @param[in] surface_height Surface height (m)
  !> @param[in] domain_top     Height of the domain (m)
  !> @return    z              Physical height (m)
  !=============================================================================
  real(kind=r_def) function eta2z_linear(eta,            &
                                         surface_height, &
                                         domain_top)     &
                                         result(z)

    implicit none

    ! Arguments
    real(kind=r_def), intent(in) :: eta
    real(kind=r_def), intent(in) :: surface_height
    real(kind=r_def), intent(in) :: domain_top

    ! Calculate physical height from eta
    z = eta*domain_top + (1.0_r_def - eta)*surface_height

    return
  end function eta2z_linear

  !=============================================================================
  !> @brief Transforms physical height to nondimensional coordinate.
  !>
  !> @details Helper routine which calculates nondimensional (terrain following)
  !>          coordinate eta from physical height using linear transformation
  !>          from Wood et al. (2013), Section 6.
  !>
  !> @param[in] z              Physical height (m)
  !> @param[in] surface_height Surface height (m)
  !> @param[in] domain_top     Height of the domain (m)
  !> @return    eta            Nondimensional (terrain following) coordinate
  !=============================================================================
  real(kind=r_def) function z2eta_linear(z,              &
                                         surface_height, &
                                         domain_top)     &
                                         result(eta)

    implicit none

    ! Arguments
    real(kind=r_def), intent(in) :: z
    real(kind=r_def), intent(in) :: surface_height
    real(kind=r_def), intent(in) :: domain_top

    ! Calculate eta from physical height
    eta = (z - surface_height)/(domain_top - surface_height)

    return
  end function z2eta_linear

  !=============================================================================
  !> @brief Set horizontal domain.
  !>
  !> @details Sets horizontal domain size for later use.
  !>
  !>          For now this is required only for the "biperiodic transforms" of
  !>          Cartesian coordinates for analytic orography profiles (Schar and
  !>          Witch-of-Agnesi mountains).
  !>
  !> @param[in] domain Global domain extents
  !=============================================================================
  subroutine set_horizontal_domain_size( domain )

    implicit none

    ! Arguments
    type(domain_type), intent(in) :: domain

    real(r_def) :: min_x, min_y
    real(r_def) :: max_x, max_y

    ! Calculate domain length and width
    if (domain%is_lonlat()) then
      min_x = domain%minimum_lonlat(axis=1)
      max_x = domain%maximum_lonlat(axis=1)
      min_y = domain%minimum_lonlat(axis=2)
      max_y = domain%maximum_lonlat(axis=2)
    else
      min_x = domain%minimum_xy(axis=1)
      max_x = domain%maximum_xy(axis=1)
      min_y = domain%minimum_xy(axis=2)
      max_y = domain%maximum_xy(axis=2)
    end if

      domain_length = max_x - min_x
      domain_width  = max_y - min_y

    write(log_scratch_space,'(A,A)') &
          "calc_domain_size_horizontal: Calculated horizontal domain size."
    call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    write(log_scratch_space,'(A,ES15.3E3)') "domain_length = ", domain_length
    call log_event(log_scratch_space, LOG_LEVEL_DEBUG)
    write(log_scratch_space,'(A,ES15.3E3)') "domain_width  = ", domain_width
    call log_event(log_scratch_space, LOG_LEVEL_DEBUG)

    return
  end subroutine set_horizontal_domain_size


  !=============================================================================
  !> @brief Transforms Cartesian coordinates to biperiodic domain.
  !>
  !> @details Helper routine which calculates "biperiodic transforms" of
  !>          Cartesian coordinates for analytic orography profiles (Schar and
  !>          Witch-of-Agnesi mountains). This enables continuation of mountains
  !>          around the boundary of Cartesian biperiodic domain.
  !>          Note: This works for analytic mountain profiles which use absolute
  !>          distance from the mountain centre. Minima of absolute distances
  !>          from the mountain centre are at x_centre and y_centre.
  !>          Periodicities are domain_length in x direction and domain_width in
  !>          y direction, so the maxima of absolute distances must be halfway
  !>          apart from the respective minima.
  !>
  !> @param[in]  x_coord     x coordinate (m)
  !> @param[in]  y_coord     y coordinate (m)
  !> @param[in]  x_centre    x coordinate centre of mountain function (m)
  !> @param[in]  y_centre    y coordinate centre of mountain function (m)
  !> @param[out] xper_coord  x coordinate biperiodic transform (m)
  !> @param[out] yper_coord  y coordinate biperiodic transform (m)
  !=============================================================================
  subroutine coord_transform_cart_biperiodic( x_coord,    &
                                              y_coord,    &
                                              x_centre,   &
                                              y_centre,   &
                                              xper_coord, &
                                              yper_coord )

    implicit none

    ! Arguments
    real(kind=r_def),    intent(in)  :: x_coord, y_coord, x_centre, y_centre
    real(kind=r_def),    intent(out) :: xper_coord, yper_coord
    ! Internal variables
    real(kind=r_def) :: xcen, ycen, lx, ly

    ! Note: Domain size is calculated in subroutine assign_orography_field
    ! before calling the selected analytic orography function.

    ! Calculate domain half-length (lx) for periodicity in x direction
    lx = domain_length/2.0_r_def

    ! Calculate domain half-width (ly) for periodicity in y direction
    ly = domain_width/2.0_r_def

    ! Calculate distance from the mountain centre in x and y directions
    xcen = x_coord - x_centre
    ycen = y_coord - y_centre

    ! Evaluate x and y coordinate biperiodic transforms
    xper_coord = lx-abs(lx-abs(xcen))
    yper_coord = ly-abs(ly-abs(ycen))

    return
  end subroutine coord_transform_cart_biperiodic

end module orography_helper_functions_mod
