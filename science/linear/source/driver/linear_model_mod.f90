!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Initialises and finalises the linear model numerical schemes.
module linear_model_mod

  use constants_mod,              only : i_def
  use field_array_mod,            only : field_array_type
  use field_mod,                  only : field_type
  use field_collection_mod,       only : field_collection_type
  use driver_modeldb_mod,         only : modeldb_type
  use tl_rk_alg_timestep_mod,     only : tl_rk_alg_init, &
                                         tl_rk_alg_final
  use tl_si_timestep_alg_mod,     only : tl_semi_implicit_alg_init, &
                                         tl_semi_implicit_alg_final
  use si_operators_alg_mod,       only : final_si_operators
  use semi_implicit_solver_alg_mod, &
                                  only : semi_implicit_solver_alg_final
  use timestepping_config_mod,    only : method,               &
                                         method_semi_implicit, &
                                         method_rk
  use log_mod,                    only : log_event, &
                                         log_scratch_space, &
                                         LOG_LEVEL_INFO,    &
                                         LOG_LEVEL_TRACE,   &
                                         LOG_LEVEL_ERROR
  use mesh_mod,                   only : mesh_type

  implicit none

  private
  public initialise_linear_model, &
         finalise_linear_model

contains

  !> @brief Completes the initialisation of the tangent linear model
  !> @param[in] model_clock Time wihtin the model.
  !> @param[in] mesh The primary mesh
  !> @param[in,out] modeldb The working data set for the model run
  !>
  subroutine initialise_linear_model( mesh,  &
                                      modeldb )
    implicit none

    type(mesh_type),     intent(in),   pointer :: mesh
    type(modeldb_type),  intent(inout), target :: modeldb

    type( field_collection_type ), pointer :: prognostic_fields => null()
    type( field_type ),            pointer :: mr(:) => null()
    type( field_collection_type ), pointer :: ls_fields
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

    type(field_collection_type), pointer :: moisture_fields => null()
    type(field_array_type), pointer      :: mr_array => null()
    type(field_array_type), pointer      :: ls_mr_array => null()
    type(field_array_type), pointer      :: ls_moist_dyn_array => null()

    ! Get pointers to field collections for use downstream
    prognostic_fields => modeldb%fields%get_field_collection(&
                                          "prognostic_fields")
    moisture_fields => modeldb%fields%get_field_collection("moisture_fields")
    call moisture_fields%get_field("mr", mr_array)
    mr => mr_array%bundle
    call moisture_fields%get_field("ls_mr", ls_mr_array)
    ls_mr => ls_mr_array%bundle
    call moisture_fields%get_field("ls_moist_dyn", ls_moist_dyn_array)
    ls_moist_dyn => ls_moist_dyn_array%bundle
    ls_fields => modeldb%fields%get_field_collection("ls_fields")

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

    select case( method )
      case( method_semi_implicit )  ! Semi-Implicit
        call semi_implicit_solver_alg_final()
        call final_si_operators()
        call tl_semi_implicit_alg_init(mesh, u, rho, theta, exner, &
                                       mr, ls_u, ls_rho, ls_theta, ls_exner, &
                                       ls_mr, ls_moist_dyn)

      case( method_rk )             ! RK
        ! Initialise and output initial conditions for first timestep

        call tl_rk_alg_init(mesh, u, rho, theta, exner, &
                            ls_u, ls_rho, ls_theta, ls_exner)

      case default
        call log_event("TL: Incorrect time stepping option chosen, "// &
                        "stopping program! ",LOG_LEVEL_ERROR)
    end select

  end subroutine initialise_linear_model

  !> @brief Finalises the remaining infrastructure and constants used by the linear model
  subroutine finalise_linear_model

    implicit none

    if ( method == method_rk )            call tl_rk_alg_final()
    if ( method == method_semi_implicit ) call tl_semi_implicit_alg_final()

  end subroutine finalise_linear_model

end module linear_model_mod
