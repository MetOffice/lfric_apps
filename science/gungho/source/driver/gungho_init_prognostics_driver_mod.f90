!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief This driver controls the initialisation of gungho's prognostics
module gungho_init_prognostics_driver_mod

  use log_mod,                          only: log_event,         &
                                              LOG_LEVEL_ERROR,   &
                                              LOG_LEVEL_INFO
  use constants_mod,                    only: r_def, i_def
  use init_gungho_prognostics_alg_mod,  only: init_u_field,      &
                                              init_theta_field,  &
                                              init_exner_field,  &
                                              init_rho_field,    &
                                              init_mr_fields
  use init_saturated_profile_alg_mod,   only: init_saturated_profile_alg
  use init_unsaturated_profile_alg_mod, only: init_unsaturated_profile_alg
  use init_thermo_profile_alg_mod,      only: init_thermo_profile_alg
  use idealised_config_mod,             only: test,                    &
                                              test_specified_profiles, &
                                              test_bryan_fritsch,      &
                                              test_grabowski_clark,    &
                                              perturb_init,            &
                                              perturb_magnitude
  use random_perturb_field_alg_mod,     only: random_perturb_field
  use sci_field_bundle_builtins_mod,    only: set_bundle_scalar
  use field_mod,                        only: field_type
  use field_collection_mod,             only: field_collection_type
  use mr_indices_mod,                   only: nummr, imr_v
  use moist_dyn_mod,                    only: num_moist_factors
  use moist_dyn_factors_alg_mod,        only: moist_dyn_factors_alg

  implicit none

  private
  public :: init_gungho_prognostics

contains

  !> @details A subroutine for initialising prognostic fields for gungho
  !> @param[in,out] prognostic_fields the collection of prognostics
  !> @param[in,out] mr Field bundle containing the moisture mixing ratios
  !> @param[in,out] moist_dyn Auxilliary fields for moist dynamics
  subroutine init_gungho_prognostics(prognostic_fields, mr, moist_dyn)

    implicit none

    ! Prognostic fields
    type( field_collection_type ), intent(inout) :: prognostic_fields
    type( field_type ),            intent(inout) :: mr(nummr), &
                                                    moist_dyn(num_moist_factors)

    ! Pointers to fields
    type( field_type ), pointer   :: u => null()
    type( field_type ), pointer   :: rho => null()
    type( field_type ), pointer   :: theta => null()
    type( field_type ), pointer   :: exner => null()

    real(kind=r_def)     :: initial_time

    call prognostic_fields%get_field('theta', theta)
    call prognostic_fields%get_field('u', u)
    call prognostic_fields%get_field('rho', rho)
    call prognostic_fields%get_field('exner', exner)

    !=== Initialise global prognostic fields ==================================!

    call log_event( "Gungho: Initialising prognostic fields", LOG_LEVEL_INFO )

    initial_time = 0.0_r_def

    call init_u_field( u , initial_time )

    ! Initialise potential temperature and water vapour
    call set_bundle_scalar(0.0_r_def, mr, nummr)
    if ( test == test_specified_profiles ) then
      call init_thermo_profile_alg( theta, mr(imr_v) )
    else if ( test == test_bryan_fritsch ) then
      call init_saturated_profile_alg( theta, mr, exner, rho, moist_dyn )
    else
      call init_theta_field( theta )
    end if
    call moist_dyn_factors_alg( moist_dyn, mr )

    ! Call subroutine to set up Exner and rho
    if ( test /= test_bryan_fritsch ) then
      call init_exner_field( exner, theta, moist_dyn, initial_time )
      call init_rho_field( rho, theta, exner, moist_dyn, initial_time )
    end if

    if ( test == test_grabowski_clark ) then
      call init_unsaturated_profile_alg( theta, mr, exner, rho, moist_dyn )
    else if (test /= test_bryan_fritsch .and. test /= test_specified_profiles) then
      call init_mr_fields( mr, theta, exner, rho, moist_dyn )
    end if

    if ( perturb_init ) then
      call log_event( "Gungho: Applying initial perturbation to theta",         &
                      LOG_LEVEL_INFO )
      call random_perturb_field(theta, perturb_magnitude)
    end if

    call log_event( "Gungho: Initialised prognostic fields", LOG_LEVEL_INFO )

    !==========================================================================!
    nullify( theta, rho, u, exner)

  end subroutine init_gungho_prognostics

end module gungho_init_prognostics_driver_mod
