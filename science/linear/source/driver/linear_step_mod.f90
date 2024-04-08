!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Steps the tangent linear app through one timestep

!> @details Handles the stepping (for a single timestep) of the
!>          linear app

module linear_step_mod

  use conservation_algorithm_mod,     only : conservation_algorithm
  use constants_mod,                  only : i_def, r_def
  use field_array_mod,                only : field_array_type
  use field_collection_mod,           only : field_collection_type
  use field_mod,                      only : field_type
  use formulation_config_mod,         only : moisture_formulation,     &
                                             moisture_formulation_dry, &
                                             use_physics
  use geometric_constants_mod,        only : get_da_at_w2
  use gungho_modeldb_mod,             only : modeldb_type
  use io_config_mod,                  only : write_conservation_diag, &
                                             write_minmax_tseries
  use log_mod,                        only : log_event, &
                                             log_scratch_space, &
                                             LOG_LEVEL_INFO,  &
                                             LOG_LEVEL_TRACE, &
                                             LOG_LEVEL_ERROR
  use mesh_mod,                       only : mesh_type
  use minmax_tseries_mod,             only : minmax_tseries
  use mr_indices_mod,                 only : nummr
  use model_clock_mod,                only : model_clock_type
  use moist_dyn_mod,                  only : num_moist_factors
  use moisture_conservation_alg_mod,  only : moisture_conservation_alg
  use moisture_fluxes_alg_mod,        only : moisture_fluxes_alg
  use tl_rk_alg_timestep_mod,         only : tl_rk_alg_step
  use tl_si_timestep_alg_mod,         only : tl_semi_implicit_alg_step
  use timestepping_config_mod,        only : method, &
                                             method_semi_implicit, &
                                             method_rk
  use field_minmax_alg_mod,           only : log_field_minmax

  implicit none

  private
  public linear_step

  contains

  !> @brief Steps the linear app through one timestep
  !> @param[in] mesh      The identifier of the primary mesh
  !> @param[in] twod_mesh The identifier of the two-dimensional mesh
  !> @param[inout] modeldbThe working data set for the model run
  !> @param[in] clock The model time
  subroutine linear_step( mesh,       &
                          twod_mesh,  &
                          modeldb,    &
                          model_clock )

    implicit none

    type( mesh_type ), pointer,      intent(in)    :: mesh
    type( mesh_type ), pointer,      intent(in)    :: twod_mesh
    type( modeldb_type ), target, intent(inout)    :: modeldb
    class(model_clock_type),         intent(in)    :: model_clock

    type( field_collection_type ), pointer :: prognostic_fields => null()
    type( field_collection_type ), pointer :: diagnostic_fields => null()
    type( field_type ),            pointer :: mr(:) => null()
    type( field_type ),            pointer :: moist_dyn(:) => null()
    type( field_collection_type ), pointer :: derived_fields => null()
    type( field_collection_type ), pointer :: ls_fields => null()
    type( field_type ),            pointer :: ls_mr(:) => null()
    type( field_type ),            pointer :: ls_moist_dyn(:) => null()

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()
    type( field_type), pointer :: ls_theta => null()
    type( field_type), pointer :: ls_u => null()
    type( field_type), pointer :: ls_rho => null()
    type( field_type), pointer :: ls_exner => null()
    type( field_type), pointer :: dA => null()  ! Areas of faces

    type(field_collection_type), pointer :: moisture_fields => null()
    type(field_array_type), pointer      :: mr_array => null()
    type(field_array_type), pointer      :: moist_dyn_array => null()
    type(field_array_type), pointer      :: ls_mr_array => null()
    type(field_array_type), pointer      :: ls_moist_dyn_array => null()

    write( log_scratch_space, '("/", A, "\ ")' ) repeat( "*", 76 )
    call log_event( log_scratch_space, LOG_LEVEL_TRACE )
    write( log_scratch_space, &
           '(A,I0)' ) 'Start of timestep ', model_clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    ! Get pointers to field collections for use downstream
    prognostic_fields => modeldb%fields%get_field_collection(&
                                          "prognostic_fields")
    diagnostic_fields => modeldb%fields%get_field_collection(&
                                          "diagnostic_fields")
    moisture_fields => modeldb%fields%get_field_collection("moisture_fields")
    call moisture_fields%get_field("mr", mr_array)
    call moisture_fields%get_field("moist_dyn", moist_dyn_array)
    mr => mr_array%bundle
    moist_dyn => moist_dyn_array%bundle
    derived_fields => modeldb%model_data%derived_fields

    ls_fields => modeldb%model_data%ls_fields
    call moisture_fields%get_field("ls_mr", ls_mr_array)
    call moisture_fields%get_field("ls_moist_dyn", ls_moist_dyn_array)
    ls_mr => ls_mr_array%bundle
    ls_moist_dyn => ls_moist_dyn_array%bundle

    ! Get pointers to fields in the prognostic/diagnostic field collections
    ! for use downstream
    call prognostic_fields%get_field('theta', theta)
    call prognostic_fields%get_field('u', u)
    call prognostic_fields%get_field('rho', rho)
    call prognostic_fields%get_field('exner', exner)
    call ls_fields%get_field('ls_theta', ls_theta)
    call ls_fields%get_field('ls_u', ls_u)
    call ls_fields%get_field('ls_rho', ls_rho)
    call ls_fields%get_field('ls_exner', ls_exner)
    dA => get_da_at_w2( mesh%get_id() )

    select case( method )
      case( method_semi_implicit )  ! Semi-Implicit
        call tl_semi_implicit_alg_step(u, rho, theta,                 &
                                       exner, mr, moist_dyn,          &
                                       ls_u, ls_rho, ls_theta,        &
                                       ls_exner, ls_mr, ls_moist_dyn, &
                                       derived_fields,                &
                                       model_clock, mesh, twod_mesh,  &
                                       modeldb%Configuration)
      case( method_rk )             ! RK
        call tl_rk_alg_step(u, rho, theta, moist_dyn, exner, mr,  &
                            ls_u, ls_rho, ls_theta,               &
                            ls_moist_dyn, ls_exner, ls_mr, model_clock)
    end select

    if ( write_conservation_diag ) then
      call conservation_algorithm( rho,              &
                                   u,                &
                                   theta,            &
                                   mr,               &
                                   exner )
      if ( moisture_formulation /= moisture_formulation_dry ) then
        call moisture_conservation_alg( rho,              &
                                        mr,               &
                                        'After timestep' )
        if ( use_physics ) then
          call log_event('use_physics not available for TL', LOG_LEVEL_ERROR)
        end if
      end if

      if (write_minmax_tseries) call minmax_tseries(u, 'u', mesh)

      call log_field_minmax( LOG_LEVEL_INFO, 'u', u )
      call log_field_minmax( LOG_LEVEL_INFO, 'theta', theta )

    end if

    write( log_scratch_space, &
           '(A,I0)' ) 'End of timestep ', model_clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '("\", A, "/ ")' ) repeat( "*", 76 )
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine linear_step

end module linear_step_mod
