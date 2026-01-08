!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Drives the execution of the adjoint miniapp.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module adjoint_test_driver_mod

  use adj_solver_lookup_cache_mod, only : adj_solver_lookup_cache_type
  use adj_trans_lookup_cache_mod,  only : adj_trans_lookup_cache_type
  use base_mesh_config_mod,        only : prime_mesh_name
  use extrusion_mod,               only : TWOD
  use fs_continuity_mod,           only : Wtheta, W3
  use field_mod,                   only : field_type
  use driver_modeldb_mod,          only : modeldb_type
  use log_mod,                     only : log_event, LOG_LEVEL_INFO
  use mesh_mod,                    only : mesh_type
  use mesh_collection_mod,         only : mesh_collection
  use sci_geometric_constants_mod, only : get_coordinates, &
                                          get_panel_id

  implicit none

  private

  public run

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief   Runs adjoint tests.
  !> @details Runs algorithm layer adjoint tests.
  !> @param[in,out]  modeldb  The model database
  subroutine run( modeldb )

    ! PSyAD generated tests
    use gen_adj_kernel_tests_mod,                   only : run_gen_adj_kernel_tests

    ! Handwritten kernel tests
    ! ./
    use adjt_sci_psykal_builtin_alg_mod,            only : adjt_convert_cart2sphere_vector_alg

    ! ./inter_function_space
    use adjt_sci_convert_hdiv_field_alg_mod,        only : adjt_sci_convert_hdiv_field_alg

    ! ./transport/mol
    use atlt_poly_adv_update_alg_mod,               only : atlt_poly_adv_update_alg
    use atlt_poly1d_vert_w3_recon_alg_mod,          only : atlt_poly1d_vert_w3_recon_alg
    use atlt_w3h_advective_update_alg_mod,          only : atlt_w3h_advective_update_alg
    !- Lookup table solutions.
    use adjt_poly1d_recon_lookup_alg_mod,           only : adjt_poly1d_recon_lookup_alg
    use adjt_poly2d_recon_lookup_alg_mod,           only : adjt_poly2d_recon_lookup_alg
    use adjt_poly_adv_upd_lookup_alg_mod,           only : adjt_poly_adv_upd_lookup_alg
    use adjt_w3h_adv_upd_lookup_alg_mod,            only : adjt_w3h_adv_upd_lookup_alg

    ! Handwritten algorithm tests
    ! ./interpolation
    use adjt_interpolation_alg_mod,                 only : adjt_interp_w3wth_to_w2_alg, &
                                                           adjt_interp_w2_to_w3wth_alg

    ! ./transport/common
    use adjt_flux_precomputations_alg_mod,          only : adjt_flux_precomputations_initialiser_alg, &
                                                           adjt_initialise_step_alg
    use adjt_wind_precomputations_alg_mod,          only : adjt_wind_precomputations_initialiser_alg
    use adjt_end_transport_step_alg_mod,            only : adjt_build_up_flux_alg
    use atlt_end_transport_step_alg_mod,            only : atlt_end_adv_step_alg, &
                                                           atlt_end_con_step_alg

    ! ./transport/mol
    use atlt_reconstruct_w3_field_alg_mod,          only : atlt_vert_w3_reconstruct_alg, &
                                                           atlt_reconstruct_w3_field_alg
    use adjt_reconstruct_w3_field_alg_mod,          only : adjt_hori_w3_reconstruct_alg
    use atlt_wt_advective_update_alg_mod,           only : atlt_hori_wt_update_alg, &
                                                           atlt_vert_wt_update_alg, &
                                                           atlt_wt_advective_update_alg
    use adjt_wt_advective_update_alg_mod,           only : adjt_hori_wt_update_alg
    use atlt_advective_and_flux_alg_mod,            only : atlt_advective_and_flux_alg
    use atlt_mol_conservative_alg_mod,              only : atlt_mol_conservative_alg
    use atlt_mol_advective_alg_mod,                 only : atlt_mol_advective_alg

    ! ./transport/control
    use atlt_transport_field_alg_mod,               only : atlt_transport_field_alg
    use atlt_wind_transport_alg_mod,                only : atlt_wind_transport_alg
    use atlt_moist_mr_transport_alg_mod,            only : atlt_moist_mr_transport_alg
    use atlt_theta_transport_alg_mod,               only : atlt_theta_transport_alg
    use adjt_transport_controller_alg_mod,          only : adjt_ls_wind_pert_rho_initialiser_alg, &
                                                           adjt_pert_wind_ls_rho_initialiser_alg
    use atlt_transport_controller_alg_mod,          only : atlt_transport_controller_initialiser_alg
    use atlt_transport_control_alg_mod,             only : atlt_transport_control_alg

    ! ./core_dynamics
    use atlt_pressure_gradient_bd_alg_mod,          only : atlt_pressure_gradient_bd_alg
    use atlt_rhs_alg_mod,                           only : atlt_rhs_alg
    use adjt_compute_vorticity_alg_mod,             only : adjt_compute_vorticity_alg
    use atlt_derive_exner_from_eos_alg_mod,         only : atlt_derive_exner_from_eos_alg
    use atlt_moist_dyn_factors_alg_mod,             only : atlt_moist_dyn_factors_alg

    ! ./solver
    use adjt_pressure_precon_alg_mod,               only : adjt_pressure_precon_alg
    use adjt_mixed_operator_alg_mod,                only : adjt_mixed_operator_alg
    use adjt_mixed_schur_preconditioner_alg_mod,    only : adjt_mixed_schur_preconditioner_alg
    use adjt_mixed_solver_alg_mod,                  only : adjt_mixed_solver_alg
    use adjt_semi_implicit_solver_step_alg_mod,     only : adjt_semi_implicit_solver_step_alg

    ! ./timestepping
    use atlt_si_timestep_alg_mod,                   only : atlt_si_timestep_alg

    implicit none

    ! Arguments
    type(modeldb_type), target, intent(inout) :: modeldb

    ! Internal variables
    type(field_type),                  pointer :: chi(:)
    type(field_type),                  pointer :: panel_id
    type(mesh_type),                   pointer :: mesh
    type(mesh_type),                   pointer :: twod_mesh
    type(adj_solver_lookup_cache_type)         :: adj_solver_lookup_cache
    type(adj_trans_lookup_cache_type)          :: adj_trans_lookup_cache

    mesh => mesh_collection%get_mesh( prime_mesh_name )
    chi => get_coordinates( mesh%get_id() )
    panel_id => get_panel_id( mesh%get_id() )
    twod_mesh => mesh_collection%get_mesh( mesh, TWOD )
    call adj_solver_lookup_cache%initialise(mesh)
    call adj_trans_lookup_cache%initialise(mesh)

    call log_event( "TESTING generated adjoint kernels", LOG_LEVEL_INFO )
    call run_gen_adj_kernel_tests( mesh, chi, panel_id )

    call log_event( "TESTING adjoint kernels", LOG_LEVEL_INFO )

    ! ./transport/mol
    call atlt_poly_adv_update_alg( mesh )
    call atlt_poly1d_vert_w3_recon_alg( mesh )
    call atlt_w3h_advective_update_alg( mesh )
    ! -- Lookup table solutions.
    call adjt_poly1d_recon_lookup_alg( mesh, adj_trans_lookup_cache )
    call adjt_poly2d_recon_lookup_alg( mesh, Wtheta, adj_trans_lookup_cache )
    call adjt_poly2d_recon_lookup_alg( mesh, W3, adj_trans_lookup_cache )
    call adjt_poly_adv_upd_lookup_alg( mesh, adj_trans_lookup_cache )
    call adjt_w3h_adv_upd_lookup_alg( mesh, adj_trans_lookup_cache )

    ! ./core_dynamics
    call atlt_pressure_gradient_bd_alg( mesh )

    ! ./inter_function_space
    call adjt_sci_convert_hdiv_field_alg( mesh, chi, panel_id )

    call log_event( "TESTING misc adjoints", LOG_LEVEL_INFO )
    ! ./
    call adjt_convert_cart2sphere_vector_alg( mesh )

    call log_event( "TESTING adjoint algorithms", LOG_LEVEL_INFO )
    ! ./interpolation
    call adjt_interp_w3wth_to_w2_alg( mesh )
    call adjt_interp_w2_to_w3wth_alg( mesh )

    ! ./transport/common
    call adjt_initialise_step_alg( mesh, modeldb%clock )
    call adjt_flux_precomputations_initialiser_alg( mesh, modeldb%clock )
    call adjt_wind_precomputations_initialiser_alg( mesh, modeldb%clock )
    call adjt_build_up_flux_alg( mesh, modeldb%clock )
    call atlt_end_adv_step_alg( mesh, modeldb%clock )
    call atlt_end_con_step_alg( mesh, modeldb%clock )

    ! ./transport/mol
    call adjt_hori_w3_reconstruct_alg( mesh, modeldb%clock )
    call atlt_vert_w3_reconstruct_alg( mesh, modeldb%clock )
    call atlt_reconstruct_w3_field_alg( mesh, modeldb%clock )
    call atlt_hori_wt_update_alg( mesh, modeldb%clock )
    call atlt_vert_wt_update_alg( mesh, modeldb%clock )
    call adjt_hori_wt_update_alg( mesh, modeldb%clock )
    call atlt_wt_advective_update_alg( mesh, modeldb%clock )
    call atlt_advective_and_flux_alg( mesh, modeldb%clock )
    call atlt_mol_conservative_alg( mesh, modeldb%clock )
    call atlt_mol_advective_alg( mesh, modeldb%clock )

    ! ./transport/control
    call atlt_transport_field_alg( mesh, modeldb%clock )
    call atlt_wind_transport_alg( mesh, modeldb%clock )
    call atlt_moist_mr_transport_alg( mesh, modeldb%clock )
    call atlt_theta_transport_alg( mesh, modeldb%clock )
    call adjt_ls_wind_pert_rho_initialiser_alg( mesh, modeldb%clock )
    call adjt_pert_wind_ls_rho_initialiser_alg( mesh, modeldb%clock )
    call atlt_transport_controller_initialiser_alg( mesh, modeldb%clock )
    call atlt_transport_control_alg( mesh, modeldb%clock )

    ! ./core_dynamics
    call atlt_rhs_alg( mesh, modeldb%clock )
    call adjt_compute_vorticity_alg( mesh )
    call atlt_derive_exner_from_eos_alg( mesh )
    call atlt_moist_dyn_factors_alg( mesh )

    ! ./solver
    call adjt_pressure_precon_alg( modeldb, mesh, modeldb%clock )
    call adjt_mixed_operator_alg( mesh, modeldb%clock )
    call adjt_mixed_schur_preconditioner_alg( modeldb,  mesh, modeldb%clock )
    call adjt_mixed_solver_alg( modeldb, mesh, modeldb%clock )
    call adjt_semi_implicit_solver_step_alg( modeldb, mesh, modeldb%clock )

    ! ./timestepping
    call atlt_si_timestep_alg( modeldb, mesh, twod_mesh, 1 )
    call atlt_si_timestep_alg( modeldb, mesh, twod_mesh, 2 )

    call log_event( "TESTING COMPLETE", LOG_LEVEL_INFO )

    call adj_solver_lookup_cache%finalise()
    call adj_trans_lookup_cache%finalise()

  end subroutine run

end module adjoint_test_driver_mod
