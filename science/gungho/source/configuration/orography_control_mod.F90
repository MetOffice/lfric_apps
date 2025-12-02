!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!-----------------------------------------------------------------------
!> @brief Sets analytic orography parameters for the model.
!>
!> @details Contains subroutines for reading analytic orography parameters and
!>          constructing analytic orography types. The parameters are stored in
!>          relevant namelist files. Each namelist sets up the selected
!>          orography profile.
!-------------------------------------------------------------------------------
module orography_control_mod

  use constants_mod,          only : i_def, str_short
  use base_mesh_config_mod,   only : geometry,           &
                                     geometry_spherical, &
                                     geometry_planar
  use orography_config_mod,   only : profile,        &
                                     profile_schar,  &
                                     profile_agnesi, &
                                     profile_bell,   &
                                     profile_dcmip200
  use log_mod,                only : log_event,         &
                                     log_scratch_space, &
                                     LOG_LEVEL_INFO
  use analytic_orography_mod, only : orography_profile

  implicit none

  private

  public :: set_orography_option

contains

  !=============================================================================
  !> @brief Sets analytic orography options.
  !>
  !> @details Reads namelists with the parameters required for setting up
  !>          analytic orography profiles and initialises corresponding types.
  !>          There are currently four profiles to choose from:
  !>          1) Schar mountain,
  !>          2) Witch-of-Agnesi mountain,
  !>          3) DCMIP2.0.0.-type mountain,
  !>          4) Bell-shaped mountain.
  !>          The default option in no orography (flat planet surface).
  !=============================================================================
  subroutine set_orography_option()

    implicit none

    ! ----------- Deallocate abstract orography type if allocated -------------!
    if ( allocated (orography_profile) ) deallocate(orography_profile)

    ! ----------- Read orography namelist and set orography type --------------!
    select case( profile )
      ! Witch-of-Agnesi orography
      case( profile_agnesi )
        if ( geometry == geometry_spherical ) then
          ! Read parameters for Witch-of-Agnesi mountain in spherical
          ! coordinates and initialise the corresponding type
           call set_orography_agnesi_spherical()
        else
          ! Read parameters for Witch-of-Agnesi mountain in Cartesian
          ! coordinates and initialise the corresponding type
          call set_orography_agnesi_cartesian()
        end if
      ! Schar orography
      case( profile_schar )
        if ( geometry == geometry_spherical ) then
          ! Read parameters for Schar mountain in (lon,lat) coordinates and
          ! initialise the corresponding type
           call set_orography_schar_spherical()
        else
          ! Read parameters for Schar mountain in Cartesian coordinates and
          ! initialise the corresponding type
           call set_orography_schar_cartesian()
        end if
      ! DCMIP200 orography
      case( profile_dcmip200 )
        if ( geometry == geometry_spherical ) then
          ! Read parameters for dcmip200 mountain in spherical
          ! coordinates and initialise the corresponding type
          call set_orography_dcmip200_spherical()
        end if
      ! Bell-shaped orography
      case( profile_bell )
        if ( geometry == geometry_planar ) then
          ! Read parameters for bell-shaped mountain in Cartesian
          ! coordinates and initialise the corresponding type
          call set_orography_bell_cartesian()
        else
          call set_orography_bell_spherical()
        end if
      ! No orography (default)
      case default
        write(log_scratch_space,'(A,A)') &
              "set_orography_option: No analytic orography set."
        call log_event(log_scratch_space, LOG_LEVEL_INFO)
    end select

    return
  end subroutine set_orography_option

  !=============================================================================
  !> @brief Initialises analytic orography type for Witch-of-Agnesi mountain in
  !>        (lon,lat) coordinates using corresponding namelist parameters.
  !=============================================================================
  subroutine set_orography_agnesi_spherical()

    use agnesi_orography_spherical_mod,        only : agnesi_spherical_type
    use orography_agnesi_spherical_config_mod, only : mountain_height, &
                                                      half_width,      &
                                                      lambda_centre,   &
                                                      phi_centre,      &
                                                      lambda_focus,    &
                                                      phi_focus

    implicit none

    ! ----------- Initialise Witch-of-Agnesi spherical orography type ---------!
    allocate( orography_profile, &
              source = agnesi_spherical_type( mountain_height, &
                                              half_width,      &
                                              lambda_centre,   &
                                              phi_centre,      &
                                              lambda_focus,    &
                                              phi_focus ) )

    write(log_scratch_space,'(A,A)') "set_orography_agnesi_spherical: "// &
          "Set analytic orography type to spherical Witch-of-Agnesi mountain."
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    return
  end subroutine set_orography_agnesi_spherical

  !=============================================================================
  !> @brief Initialises analytic orography type for Witch-of-Agnesi mountain in
  !>        Cartesian coordinates using corresponding namelist parameters.
  !=============================================================================
  subroutine set_orography_agnesi_cartesian()

    use agnesi_orography_cartesian_mod,        only : agnesi_cartesian_type
    use orography_agnesi_cartesian_config_mod, only : mountain_height, &
                                                      half_width_x,    &
                                                      half_width_y,    &
                                                      x_centre,        &
                                                      y_centre,        &
                                                      direction

    implicit none

    ! ----------- Initialise Witch-of-Agnesi Cartesian orography type ---------!
    allocate( orography_profile,                               &
              source = agnesi_cartesian_type( mountain_height, &
                                              half_width_x,    &
                                              half_width_y,    &
                                              x_centre,        &
                                              y_centre,        &
                                              direction ) )

    write(log_scratch_space,'(A,A)') "set_orography_agnesi_cartesian: "// &
          "Set analytic orography type to Cartesian Witch-of-Agnesi mountain."
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    return
  end subroutine set_orography_agnesi_cartesian

  !=============================================================================
  !> @brief Initialises analytic orography type for Schar mountain in
  !>        (lon,lat) coordinates using corresponding namelist parameters.
  !=============================================================================
  subroutine set_orography_schar_spherical()

    use schar_orography_spherical_mod,        only : schar_spherical_type
    use orography_schar_spherical_config_mod, only : mountain_height, &
                                                     half_width,      &
                                                     wavelength,      &
                                                     lambda_centre,   &
                                                     phi_centre

    implicit none

    ! ----------- Initialise Schar spherical orography type -------------------!
    allocate( orography_profile,                              &
              source = schar_spherical_type( mountain_height, &
                                             half_width,      &
                                             wavelength,      &
                                             lambda_centre,   &
                                             phi_centre ) )

    write(log_scratch_space,'(A,A)') "set_orography_schar_spherical: "// &
          "Set analytic orography type to spherical Schar mountain."
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    return
  end subroutine set_orography_schar_spherical

  !=============================================================================
  !> @brief Initialises analytic orography type for Schar mountain in
  !>        Cartesian coordinates using corresponding namelist parameters.
  !=============================================================================
  subroutine set_orography_schar_cartesian()

    use schar_orography_cartesian_mod,        only : schar_cartesian_type
    use orography_schar_cartesian_config_mod, only : mountain_height, &
                                                     half_width_x,    &
                                                     half_width_y,    &
                                                     wavelength,      &
                                                     x_centre,        &
                                                     y_centre,        &
                                                     direction

    implicit none

    ! ----------- Initialise Schar Cartesian orography type -------------------!
    allocate( orography_profile,                              &
              source = schar_cartesian_type( mountain_height, &
                                             half_width_x,    &
                                             half_width_y,    &
                                             wavelength,      &
                                             x_centre,        &
                                             y_centre,        &
                                             direction ) )

    write(log_scratch_space,'(A,A)') "set_orography_schar_cartesian: "// &
          "Set analytic orography type to Cartesian Schar mountain."
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    return
  end subroutine set_orography_schar_cartesian

  !=============================================================================
  !> @brief Initialises analytic orography type for DCMIP200 mountain in
  !>        (lon,lat) coordinates using corresponding namelist parameters.
  !=============================================================================
  subroutine set_orography_dcmip200_spherical()

    use dcmip200_orography_spherical_mod,        only : dcmip200_spherical_type
    use orography_dcmip200_spherical_config_mod, only : mountain_height, &
                                                        radius,          &
                                                        osc_half_width,  &
                                                        lambda_centre,   &
                                                        phi_centre

    implicit none

    ! ----------- Initialise DCMIP200 spherical orography type ----------------!
    allocate( orography_profile,                                 &
              source = dcmip200_spherical_type( mountain_height, &
                                                radius,          &
                                                osc_half_width,  &
                                                lambda_centre,   &
                                                phi_centre ) )

    write(log_scratch_space,'(A,A)') "set_orography_dcmip200_spherical: "// &
          "Set analytic orography type to spherical dcmip200 mountain."
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    return
  end subroutine set_orography_dcmip200_spherical


  !=============================================================================
  !> @brief Initialises analytic orography type for bell-shaped mountain in
  !>        Cartesian coordinates using corresponding namelist parameters.
  !=============================================================================
  subroutine set_orography_bell_cartesian()

    use bell_orography_cartesian_mod,        only : bell_cartesian_type
    use orography_bell_cartesian_config_mod, only : mountain_height, &
                                                    half_width_x,    &
                                                    half_width_y,    &
                                                    x_centre,        &
                                                    y_centre,        &
                                                    direction

    implicit none

    ! ----------- Initialise bell-shaped Cartesian orography type ---------!
    allocate( orography_profile,                               &
              source = bell_cartesian_type( mountain_height, &
                                            half_width_x,    &
                                            half_width_y,    &
                                            x_centre,        &
                                            y_centre,        &
                                            direction ) )

    write(log_scratch_space,'(A,A)') "set_orography_bell_cartesian: "// &
          "Set analytic orography type to Cartesian bell-shaped mountain."
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    return
  end subroutine set_orography_bell_cartesian

  !=============================================================================
  !> @brief Initialises analytic orography type for two bell-shaped mountains in
  !>        Spherical coordinates using corresponding namelist parameters.
  !=============================================================================
  subroutine set_orography_bell_spherical()

    use bell_orography_spherical_mod,        only : bell_spherical_type
    use orography_bell_spherical_config_mod, only : mountain_height,&
                                                    radius_lon,     &
                                                    radius_lat,     &
                                                    lambda_centre1, &
                                                    lambda_centre2, &
                                                    phi_centre1,    &
                                                    phi_centre2

    implicit none

    ! ----------- Initialise bell-shaped Spherical orography type ---------!
    allocate( orography_profile,                             &
              source = bell_spherical_type( mountain_height, &
                                            radius_lon,      &
                                            radius_lat,      &
                                            lambda_centre1,  &
                                            lambda_centre2,  &
                                            phi_centre1,     &
                                            phi_centre2) )

    write(log_scratch_space,'(A,A)') "set_orography_bell_spherical: "// &
          "Set analytic orography type to Spherical bell-shaped mountain."
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    return
  end subroutine set_orography_bell_spherical


end module orography_control_mod

