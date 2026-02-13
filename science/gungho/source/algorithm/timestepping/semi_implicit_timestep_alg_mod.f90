MODULE semi_implicit_timestep_alg_mod
  USE constants_mod, ONLY: i_def, r_def, l_def, str_def
  USE fs_continuity_mod, ONLY: Wtheta, W2
  USE log_mod, ONLY: log_event, log_scratch_space, LOG_LEVEL_INFO
  USE extrusion_mod, ONLY: TWOD
  USE namelist_mod, ONLY: namelist_type
  USE sci_fem_constants_mod, ONLY: get_mass_matrix_fe, get_mass_matrix_fv
  USE sci_field_bundle_builtins_mod, ONLY: clone_bundle, bundle_axpy, add_bundle, copy_bundle, set_bundle_scalar
  USE timestep_method_mod, ONLY: timestep_method_type
  USE section_choice_config_mod, ONLY: cloud, cloud_um, aerosol, aerosol_um
  USE physics_config_mod, ONLY: blayer_placement, blayer_placement_fast, convection_placement, convection_placement_fast, &
&stochastic_physics_placement, stochastic_physics_placement_fast, smagorinsky_placement, smagorinsky_placement_outer
  USE aerosol_config_mod, ONLY: glomap_mode, glomap_mode_dust_and_clim, glomap_mode_ukca
  USE formulation_config_mod, ONLY: use_physics, dlayer_on, use_wavedynamics, moisture_formulation, moisture_formulation_dry, &
&exner_from_eos
  USE io_config_mod, ONLY: write_conservation_diag, write_diag, use_xios_io, diagnostic_frequency, checkpoint_read
  USE initialization_config_mod, ONLY: init_option, init_option_checkpoint_dump, lbc_option_um2lfric_file
  USE mixed_solver_config_mod, ONLY: guess_np1, reference_reset_time
  USE timestepping_config_mod, ONLY: alpha, outer_iterations, inner_iterations, spinup_alpha
  USE transport_config_mod, ONLY: cheap_update, transport_ageofair
  USE derived_config_mod, ONLY: bundle_size
  USE boundaries_config_mod, ONLY: limited_area, blend_frequency, blend_frequency_inner, blend_frequency_outer, &
&blend_frequency_final
  USE finite_element_config_mod, ONLY: element_order_h, element_order_v
  USE field_array_mod, ONLY: field_array_type
  USE field_mod, ONLY: field_type
  USE field_parent_mod, ONLY: field_parent_type
  USE field_collection_mod, ONLY: field_collection_type
  USE r_tran_field_mod, ONLY: r_tran_field_type
  USE io_value_mod, ONLY: io_value_type, get_io_value
  USE driver_modeldb_mod, ONLY: modeldb_type
  USE mesh_mod, ONLY: mesh_type
  USE mesh_collection_mod, ONLY: mesh_collection
  USE model_clock_mod, ONLY: model_clock_type
  USE operator_mod, ONLY: operator_type
  USE rhs_alg_mod, ONLY: rhs_alg
  USE gungho_transport_control_alg_mod, ONLY: gungho_transport_control_alg_init, gungho_transport_control_alg
  USE si_operators_alg_mod, ONLY: create_si_operators, compute_si_operators, final_si_operators
  USE fast_physics_alg_mod, ONLY: fast_physics
  USE slow_physics_alg_mod, ONLY: slow_physics
  USE checks_and_balances_alg_mod, ONLY: check_fields
  USE semi_implicit_solver_alg_mod, ONLY: semi_implicit_solver_alg_init, semi_implicit_solver_alg_step, &
&semi_implicit_solver_alg_final
  USE derive_exner_from_eos_alg_mod, ONLY: derive_exner_from_eos
  USE sci_mass_matrix_solver_alg_mod, ONLY: mass_matrix_solver_alg
  USE moist_dyn_factors_alg_mod, ONLY: moist_dyn_factors_alg
  USE update_prognostic_scalars_alg_mod, ONLY: update_prognostic_scalars_alg
  USE mixing_alg_mod, ONLY: mixing_alg
  USE si_diagnostics_mod, ONLY: output_diags_for_si
  USE predictors_alg_mod, ONLY: predictors_alg
  USE limited_area_lbc_alg_mod, ONLY: lam_solver_lbc, lam_blend_lbc
  USE lam_rhs_alg_mod, ONLY: calc_rhs_lbc, apply_mask_rhs
  USE calc_phys_predictors_alg_mod, ONLY: calc_phys_predictors_alg
  USE map_physics_fields_alg_mod, ONLY: map_physics_fields_alg
  USE mr_indices_mod, ONLY: nummr
  USE moist_dyn_mod, ONLY: num_moist_factors, gas_law
  USE field_indices_mod, ONLY: igh_u, igh_t, igh_d, igh_p
  USE mixing_config_mod, ONLY: smagorinsky
  USE smagorinsky_alg_mod, ONLY: smagorinsky_alg
  USE cld_incs_mod, ONLY: cld_incs_init, cld_incs_output
  USE timing_mod, ONLY: start_timing, stop_timing, tik, LPROF
  USE ageofair_alg_mod, ONLY: ageofair_update
  IMPLICIT NONE
  PRIVATE
  TYPE, EXTENDS(timestep_method_type), PUBLIC :: semi_implicit_timestep_type
    PRIVATE
    LOGICAL(KIND = l_def) :: use_moisture
    LOGICAL(KIND = l_def) :: output_cld_incs
    TYPE(field_type), ALLOCATABLE :: state(:)
    TYPE(field_type), ALLOCATABLE :: state_n(:)
    TYPE(field_type), ALLOCATABLE :: state_after_slow(:)
    TYPE(field_type), ALLOCATABLE :: advected_state(:)
    TYPE(field_type), ALLOCATABLE :: mr_n(:), mr_after_adv(:), mr_after_slow(:)
    TYPE(field_type), ALLOCATABLE :: rhs_n(:), rhs_np1(:), rhs_adv(:)
    TYPE(field_type), ALLOCATABLE :: adv_inc_prev(:)
    TYPE(field_type), ALLOCATABLE :: rhs_phys(:), rhs_lbc(:)
    TYPE(field_type) :: dtheta, dtheta_cld
    TYPE(field_type) :: du
    TYPE(field_type) :: wind_prev
    TYPE(r_tran_field_type) :: total_dry_flux
    CONTAINS
    PRIVATE
    PROCEDURE, PUBLIC :: step => semi_implicit_alg_step
    PROCEDURE, PUBLIC :: finalise => semi_implicit_alg_final
    PROCEDURE, NOPASS :: run_init
    PROCEDURE, NOPASS :: run_step
    PROCEDURE, NOPASS :: conditional_collection_copy
  END TYPE semi_implicit_timestep_type
  INTERFACE semi_implicit_timestep_type
    MODULE PROCEDURE semi_implicit_alg_init
  END INTERFACE semi_implicit_timestep_type
  CONTAINS
  FUNCTION semi_implicit_alg_init(modeldb) RESULT(self)
    IMPLICIT NONE
    TYPE(semi_implicit_timestep_type) :: self
    TYPE(modeldb_type), INTENT(IN), TARGET :: modeldb
    TYPE(field_collection_type), POINTER :: prognostic_fields
    TYPE(field_collection_type), POINTER :: moisture_fields
    TYPE(field_type), POINTER :: u
    TYPE(field_type), POINTER :: rho
    TYPE(field_type), POINTER :: theta
    TYPE(field_type), POINTER :: exner
    TYPE(field_array_type), POINTER :: mr_array
    TYPE(field_type), POINTER :: mr(:)
    CLASS(model_clock_type), POINTER :: model_clock
    model_clock => modeldb % clock
    prognostic_fields => modeldb % fields % get_field_collection("prognostic_fields")
    CALL prognostic_fields % get_field('theta', theta)
    CALL prognostic_fields % get_field('u', u)
    CALL prognostic_fields % get_field('rho', rho)
    CALL prognostic_fields % get_field('exner', exner)
    moisture_fields => modeldb % fields % get_field_collection("moisture_fields")
    CALL moisture_fields % get_field("mr", mr_array)
    mr => mr_array % bundle
    CALL run_init(self, u, rho, theta, exner, mr, prognostic_fields, moisture_fields, model_clock)
    NULLIFY(prognostic_fields, moisture_fields, mr_array, model_clock, u, rho, theta, exner, mr)
  END FUNCTION semi_implicit_alg_init
  SUBROUTINE semi_implicit_alg_step(self, modeldb)
    IMPLICIT NONE
    CLASS(semi_implicit_timestep_type), INTENT(INOUT) :: self
    TYPE(modeldb_type), INTENT(IN), TARGET :: modeldb
    TYPE(mesh_type), POINTER :: mesh
    TYPE(mesh_type), POINTER :: twod_mesh
    CLASS(model_clock_type), POINTER :: model_clock
    TYPE(field_collection_type), POINTER :: prognostic_fields
    TYPE(field_collection_type), POINTER :: moisture_fields
    TYPE(field_type), POINTER :: u
    TYPE(field_type), POINTER :: rho
    TYPE(field_type), POINTER :: theta
    TYPE(field_type), POINTER :: exner
    TYPE(field_array_type), POINTER :: mr_array
    TYPE(field_type), POINTER :: mr(:)
    TYPE(field_array_type), POINTER :: moist_dyn_array
    TYPE(field_type), POINTER :: moist_dyn(:)
    TYPE(field_collection_type), POINTER :: adv_tracer_all_outer
    TYPE(field_collection_type), POINTER :: adv_tracer_last_outer
    TYPE(field_collection_type), POINTER :: con_tracer_all_outer
    TYPE(field_collection_type), POINTER :: con_tracer_last_outer
    TYPE(field_collection_type), POINTER :: derived_fields
    TYPE(field_collection_type), POINTER :: radiation_fields
    TYPE(field_collection_type), POINTER :: microphysics_fields
    TYPE(field_collection_type), POINTER :: electric_fields
    TYPE(field_collection_type), POINTER :: orography_fields
    TYPE(field_collection_type), POINTER :: turbulence_fields
    TYPE(field_collection_type), POINTER :: convection_fields
    TYPE(field_collection_type), POINTER :: cloud_fields
    TYPE(field_collection_type), POINTER :: surface_fields
    TYPE(field_collection_type), POINTER :: soil_fields
    TYPE(field_collection_type), POINTER :: snow_fields
    TYPE(field_collection_type), POINTER :: chemistry_fields
    TYPE(field_collection_type), POINTER :: aerosol_fields
    TYPE(field_collection_type), POINTER :: stph_fields
    TYPE(field_collection_type), POINTER :: lbc_fields
    TYPE(io_value_type), POINTER :: temp_corr_io_value
    REAL(KIND = r_def) :: dt
    REAL(KIND = r_def) :: dtemp_encorr
    prognostic_fields => modeldb % fields % get_field_collection("prognostic_fields")
    model_clock => modeldb % clock
    CALL prognostic_fields % get_field('theta', theta)
    CALL prognostic_fields % get_field('u', u)
    CALL prognostic_fields % get_field('rho', rho)
    CALL prognostic_fields % get_field('exner', exner)
    dt = REAL(model_clock % get_seconds_per_step(), r_def)
    moisture_fields => modeldb % fields % get_field_collection("moisture_fields")
    lbc_fields => modeldb % fields % get_field_collection("lbc_fields")
    radiation_fields => modeldb % fields % get_field_collection("radiation_fields")
    CALL moisture_fields % get_field("mr", mr_array)
    mr => mr_array % bundle
    CALL moisture_fields % get_field("moist_dyn", moist_dyn_array)
    moist_dyn => moist_dyn_array % bundle
    mesh => theta % get_mesh()
    twod_mesh => mesh_collection % get_mesh(mesh, TWOD)
    adv_tracer_all_outer => modeldb % fields % get_field_collection("adv_tracer_all_outer")
    adv_tracer_last_outer => modeldb % fields % get_field_collection("adv_tracer_last_outer")
    con_tracer_all_outer => modeldb % fields % get_field_collection("con_tracer_all_outer")
    con_tracer_last_outer => modeldb % fields % get_field_collection("con_tracer_last_outer")
    derived_fields => modeldb % fields % get_field_collection("derived_fields")
    microphysics_fields => modeldb % fields % get_field_collection("microphysics_fields")
    turbulence_fields => modeldb % fields % get_field_collection("turbulence_fields")
    convection_fields => modeldb % fields % get_field_collection("convection_fields")
    cloud_fields => modeldb % fields % get_field_collection("cloud_fields")
    surface_fields => modeldb % fields % get_field_collection("surface_fields")
    soil_fields => modeldb % fields % get_field_collection("soil_fields")
    snow_fields => modeldb % fields % get_field_collection("snow_fields")
    chemistry_fields => modeldb % fields % get_field_collection("chemistry_fields")
    aerosol_fields => modeldb % fields % get_field_collection("aerosol_fields")
    stph_fields => modeldb % fields % get_field_collection("stph_fields")
    electric_fields => modeldb % fields % get_field_collection("electric_fields")
    orography_fields => modeldb % fields % get_field_collection("orography_fields")
    temp_corr_io_value => get_io_value(modeldb % values, 'temperature_correction_io_value')
    dtemp_encorr = dt * temp_corr_io_value % data(1)
    CALL run_step(self, modeldb, u, rho, theta, exner, mr, moist_dyn, adv_tracer_all_outer, adv_tracer_last_outer, &
&con_tracer_all_outer, con_tracer_last_outer, prognostic_fields, moisture_fields, derived_fields, radiation_fields, &
&microphysics_fields, electric_fields, orography_fields, turbulence_fields, convection_fields, cloud_fields, surface_fields, &
&soil_fields, snow_fields, chemistry_fields, aerosol_fields, stph_fields, lbc_fields, model_clock, dtemp_encorr, mesh, twod_mesh)
    NULLIFY(model_clock, mesh, prognostic_fields, moisture_fields, u, rho, theta, exner, mr, moist_dyn, adv_tracer_all_outer, &
&adv_tracer_last_outer, con_tracer_all_outer, con_tracer_last_outer, derived_fields, radiation_fields, microphysics_fields, &
&electric_fields, orography_fields, turbulence_fields, convection_fields, cloud_fields, surface_fields, soil_fields, snow_fields, &
&chemistry_fields, aerosol_fields, stph_fields, lbc_fields, temp_corr_io_value, mr_array, moist_dyn_array)
  END SUBROUTINE semi_implicit_alg_step
  SUBROUTINE run_init(self, u, rho, theta, exner, mr, prognostic_fields, moisture_fields, model_clock)
    USE semi_implicit_timestep_alg_mod_psy, ONLY: invoke_0
    IMPLICIT NONE
    TYPE(semi_implicit_timestep_type), INTENT(INOUT) :: self
    TYPE(field_type), INTENT(IN) :: u, rho, theta, exner
    TYPE(field_type), INTENT(IN), OPTIONAL :: mr(:)
    TYPE(field_collection_type), INTENT(INOUT) :: prognostic_fields
    TYPE(field_collection_type), INTENT(INOUT) :: moisture_fields
    CLASS(model_clock_type), INTENT(IN) :: model_clock
    TYPE(mesh_type), POINTER :: mesh
    TYPE(field_type), POINTER :: rho_ref, theta_ref, exner_ref
    TYPE(field_array_type), POINTER :: moist_dyn_ref_array
    TYPE(field_type), POINTER :: moist_dyn_ref(:)
    REAL(KIND = r_def) :: cast_dt
    INTEGER(KIND = i_def) :: reference_reset_freq
    mesh => theta % get_mesh()
    ALLOCATE(self % state(bundle_size))
    ALLOCATE(self % state_n(bundle_size))
    ALLOCATE(self % state_after_slow(bundle_size))
    ALLOCATE(self % advected_state(bundle_size))
    ALLOCATE(self % rhs_n(bundle_size))
    ALLOCATE(self % rhs_np1(bundle_size))
    ALLOCATE(self % rhs_adv(bundle_size))
    ALLOCATE(self % adv_inc_prev(bundle_size))
    ALLOCATE(self % rhs_phys(bundle_size))
    ALLOCATE(self % rhs_lbc(bundle_size))
    ALLOCATE(self % mr_n(nummr))
    ALLOCATE(self % mr_after_slow(nummr))
    ALLOCATE(self % mr_after_adv(nummr))
    self % use_moisture = (moisture_formulation /= moisture_formulation_dry)
    self % output_cld_incs = (self % use_moisture .AND. write_diag .AND. use_xios_io)
    CALL u % copy_field_properties(self % state(igh_u))
    CALL theta % copy_field_properties(self % state(igh_t))
    CALL rho % copy_field_properties(self % state(igh_d))
    CALL exner % copy_field_properties(self % state(igh_p))
    CALL clone_bundle(self % state, self % state_n, bundle_size)
    CALL clone_bundle(self % state, self % state_after_slow, bundle_size)
    CALL clone_bundle(self % state, self % advected_state, bundle_size)
    CALL clone_bundle(self % state, self % rhs_n, bundle_size)
    CALL clone_bundle(self % state, self % rhs_np1, bundle_size)
    CALL clone_bundle(self % state, self % rhs_adv, bundle_size)
    CALL clone_bundle(self % state, self % rhs_phys, bundle_size)
    CALL clone_bundle(self % state, self % rhs_lbc, bundle_size)
    CALL clone_bundle(self % state, self % adv_inc_prev, bundle_size)
    CALL theta % copy_field_properties(self % dtheta)
    CALL theta % copy_field_properties(self % dtheta_cld)
    CALL invoke_0(self % dtheta_cld)
    CALL u % copy_field_properties(self % du)
    CALL u % copy_field_properties(self % wind_prev)
    CALL self % total_dry_flux % initialise(u % get_function_space())
    CALL clone_bundle(mr, self % mr_n, nummr)
    CALL clone_bundle(mr, self % mr_after_slow, nummr)
    IF (self % use_moisture) THEN
      CALL clone_bundle(mr, self % mr_after_adv, nummr)
    ELSE
      CALL set_bundle_scalar(0.0_r_def, self % mr_n, nummr)
      CALL set_bundle_scalar(0.0_r_def, self % mr_after_slow, nummr)
    END IF
    CALL set_bundle_scalar(0.0_r_def, self % rhs_phys, bundle_size)
    CALL create_si_operators(mesh)
    cast_dt = REAL(model_clock % get_seconds_per_step(), r_def)
    reference_reset_freq = NINT(reference_reset_time / cast_dt, i_def)
    IF (MOD(model_clock % get_first_step() - 1, reference_reset_freq) /= 0) THEN
      CALL moisture_fields % get_field("moist_dyn_ref", moist_dyn_ref_array)
      moist_dyn_ref => moist_dyn_ref_array % bundle
      CALL prognostic_fields % get_field('theta_ref', theta_ref)
      CALL prognostic_fields % get_field('rho_ref', rho_ref)
      CALL prognostic_fields % get_field('exner_ref', exner_ref)
      CALL compute_si_operators(theta_ref, rho_ref, exner_ref, model_clock, moist_dyn_ref)
      NULLIFY(theta_ref, rho_ref, exner_ref, moist_dyn_ref, moist_dyn_ref_array)
    END IF
    IF (use_wavedynamics) THEN
      CALL gungho_transport_control_alg_init(mesh)
      CALL semi_implicit_solver_alg_init(self % state)
    END IF
    NULLIFY(mesh)
    CALL log_event("semi_implicit_timestep: initialised timestepping algorithm", LOG_LEVEL_INFO)
  END SUBROUTINE run_init
  SUBROUTINE run_step(self, modeldb, u, rho, theta, exner, mr, moist_dyn, adv_tracer_all_outer, adv_tracer_last_outer, &
&con_tracer_all_outer, con_tracer_last_outer, prognostic_fields, moisture_fields, derived_fields, radiation_fields, &
&microphysics_fields, electric_fields, orography_fields, turbulence_fields, convection_fields, cloud_fields, surface_fields, &
&soil_fields, snow_fields, chemistry_fields, aerosol_fields, stph_fields, lbc_fields, model_clock, dtemp_encorr, mesh, twod_mesh)
    USE semi_implicit_timestep_alg_mod_psy, ONLY: invoke_11
    USE semi_implicit_timestep_alg_mod_psy, ONLY: invoke_10
    USE semi_implicit_timestep_alg_mod_psy, ONLY: invoke_9
    USE semi_implicit_timestep_alg_mod_psy, ONLY: invoke_update_rhs_phys_from_fast_physics
    USE semi_implicit_timestep_alg_mod_psy, ONLY: invoke_7
    USE semi_implicit_timestep_alg_mod_psy, ONLY: invoke_6
    USE semi_implicit_timestep_alg_mod_psy, ONLY: invoke_5
    USE semi_implicit_timestep_alg_mod_psy, ONLY: invoke_4
    USE semi_implicit_timestep_alg_mod_psy, ONLY: invoke_update_from_slow_physics
    USE semi_implicit_timestep_alg_mod_psy, ONLY: invoke_2
    USE semi_implicit_timestep_alg_mod_psy, ONLY: invoke_copy_init_fields_to_state
    IMPLICIT NONE
    TYPE(semi_implicit_timestep_type), INTENT(INOUT), TARGET :: self
    TYPE(modeldb_type), INTENT(IN), TARGET :: modeldb
    TYPE(field_type), INTENT(INOUT) :: u, rho, theta, exner
    TYPE(field_type), INTENT(INOUT) :: mr(nummr)
    TYPE(field_type), INTENT(INOUT) :: moist_dyn(num_moist_factors)
    TYPE(field_collection_type), INTENT(INOUT) :: adv_tracer_all_outer
    TYPE(field_collection_type), INTENT(INOUT) :: adv_tracer_last_outer
    TYPE(field_collection_type), INTENT(INOUT) :: con_tracer_all_outer
    TYPE(field_collection_type), INTENT(INOUT) :: con_tracer_last_outer
    TYPE(field_collection_type), INTENT(INOUT) :: prognostic_fields
    TYPE(field_collection_type), INTENT(INOUT) :: moisture_fields
    TYPE(field_collection_type), INTENT(INOUT) :: derived_fields
    TYPE(field_collection_type), INTENT(INOUT) :: radiation_fields
    TYPE(field_collection_type), INTENT(INOUT) :: microphysics_fields
    TYPE(field_collection_type), INTENT(INOUT) :: electric_fields
    TYPE(field_collection_type), INTENT(IN) :: orography_fields
    TYPE(field_collection_type), INTENT(INOUT) :: turbulence_fields
    TYPE(field_collection_type), INTENT(INOUT) :: convection_fields
    TYPE(field_collection_type), INTENT(INOUT) :: cloud_fields
    TYPE(field_collection_type), INTENT(INOUT) :: surface_fields
    TYPE(field_collection_type), INTENT(INOUT) :: soil_fields
    TYPE(field_collection_type), INTENT(INOUT) :: snow_fields
    TYPE(field_collection_type), INTENT(INOUT) :: chemistry_fields
    TYPE(field_collection_type), INTENT(INOUT) :: aerosol_fields
    TYPE(field_collection_type), INTENT(INOUT) :: stph_fields
    TYPE(field_collection_type), INTENT(INOUT) :: lbc_fields
    CLASS(model_clock_type), INTENT(IN) :: model_clock
    REAL(KIND = r_def), INTENT(IN) :: dtemp_encorr
    TYPE(mesh_type), INTENT(IN), POINTER :: mesh
    TYPE(mesh_type), INTENT(IN), POINTER :: twod_mesh
    TYPE(field_type), POINTER :: rho_ref
    TYPE(field_type), POINTER :: theta_ref
    TYPE(field_type), POINTER :: exner_ref
    TYPE(field_type), POINTER :: moist_dyn_ref(:)
    TYPE(field_array_type), POINTER :: moist_dyn_ref_array
    TYPE(field_type), POINTER :: mr_to_adv(:)
    TYPE(field_type) :: dcfl_tot, dcff_tot, dbcf_tot, dcfl_adv, dcff_adv, dbcf_adv
    CHARACTER(LEN = str_def), PARAMETER :: sec_tot = 'processed'
    CHARACTER(LEN = str_def), PARAMETER :: suffix_tot = 'tot'
    CHARACTER(LEN = str_def), PARAMETER :: sec_adv = 'advection'
    CHARACTER(LEN = str_def), PARAMETER :: suffix_adv = 'adv'
    REAL(KIND = r_def) :: cast_dt
    TYPE(field_type), POINTER :: ageofair
    TYPE(operator_type), POINTER :: mm_wt
    TYPE(operator_type), POINTER :: mm_vel
    INTEGER(KIND = i_def) :: outer, inner, reference_reset_freq
    REAL(KIND = r_def) :: varalpha, varbeta
    LOGICAL(KIND = l_def) :: write_moisture_diag
    TYPE(field_collection_type) :: adv_tracer_all_outer_after_slow
    TYPE(field_collection_type) :: adv_tracer_last_outer_after_slow
    TYPE(field_collection_type) :: con_tracer_all_outer_after_slow
    TYPE(field_collection_type) :: con_tracer_last_outer_after_slow
    LOGICAL(KIND = l_def) :: checkpoint_reference_fields
    TYPE(namelist_type), POINTER :: base_mesh_nml
    TYPE(namelist_type), POINTER :: initialization_nml
    TYPE(namelist_type), POINTER :: microphysics_nml
    TYPE(namelist_type), POINTER :: aerosol_nml
    TYPE(namelist_type), POINTER :: timestepping_nml
    CHARACTER(LEN = str_def) :: prime_mesh_name
    INTEGER(KIND = i_def) :: lbc_option
    LOGICAL(KIND = l_def) :: microphysics_casim
    LOGICAL(KIND = l_def) :: murk_lbc
    REAL(KIND = r_def) :: tau_r
    INTEGER(KIND = tik) :: id
    IF (LPROF) CALL start_timing(id, 'semi_implicit_timestep')
    cast_dt = REAL(model_clock % get_seconds_per_step(), r_def)
    IF (limited_area .AND. use_wavedynamics) THEN
      base_mesh_nml => modeldb % configuration % get_namelist('base_mesh')
      initialization_nml => modeldb % configuration % get_namelist('initialization')
      timestepping_nml => modeldb % configuration % get_namelist('timestepping')
      CALL base_mesh_nml % get_value('prime_mesh_name', prime_mesh_name)
      CALL initialization_nml % get_value('lbc_option', lbc_option)
      CALL timestepping_nml % get_value('tau_r', tau_r)
      IF (lbc_option == lbc_option_um2lfric_file) THEN
        aerosol_nml => modeldb % configuration % get_namelist('aerosol')
        CALL aerosol_nml % get_value('murk_lbc', murk_lbc)
      END IF
    END IF
    IF (lbc_option == lbc_option_um2lfric_file .OR. (use_physics .AND. cloud == cloud_um)) THEN
      microphysics_nml => modeldb % configuration % get_namelist('microphysics')
      CALL microphysics_nml % get_value('microphysics_casim', microphysics_casim)
    END IF
    IF (element_order_h == 0 .AND. element_order_v == 0) THEN
      mm_wt => get_mass_matrix_fv(Wtheta, mesh % get_id())
      mm_vel => get_mass_matrix_fv(W2, mesh % get_id())
    ELSE
      mm_wt => get_mass_matrix_fe(Wtheta, mesh % get_id())
      mm_vel => get_mass_matrix_fe(W2, mesh % get_id())
    END IF
    CALL invoke_copy_init_fields_to_state(self % state(igh_u), u, self % state(igh_t), theta, self % state(igh_d), rho, &
&self % state(igh_p), exner, self % total_dry_flux)
    reference_reset_freq = NINT(reference_reset_time / cast_dt, i_def)
    IF (MOD(model_clock % get_step() - 1_i_def, reference_reset_freq) == 0_i_def) THEN
      CALL compute_si_operators(self % state(igh_t), self % state(igh_d), self % state(igh_p), model_clock, moist_dyn)
      checkpoint_reference_fields = MOD(model_clock % get_first_step() - 1, reference_reset_freq) /= 0 .OR. MOD(model_clock % &
&get_last_step(), reference_reset_freq) /= 0
      IF (checkpoint_read .OR. init_option == init_option_checkpoint_dump) THEN
        IF (MOD(model_clock % get_first_step() - 1, reference_reset_freq) == 0 .AND. MOD(model_clock % get_last_step(), &
&reference_reset_freq) /= 0) THEN
          checkpoint_reference_fields = .FALSE.
        END IF
      END IF
      IF (checkpoint_reference_fields) THEN
        CALL prognostic_fields % get_field('theta_ref', theta_ref)
        CALL prognostic_fields % get_field('rho_ref', rho_ref)
        CALL prognostic_fields % get_field('exner_ref', exner_ref)
        CALL moisture_fields % get_field("moist_dyn_ref", moist_dyn_ref_array)
        moist_dyn_ref => moist_dyn_ref_array % bundle
        CALL copy_bundle(moist_dyn, moist_dyn_ref, num_moist_factors)
        CALL invoke_2(theta_ref, self % state(igh_t), rho_ref, self % state(igh_d), exner_ref, self % state(igh_p))
        NULLIFY(moist_dyn_ref_array, moist_dyn_ref, theta_ref, rho_ref, exner_ref)
      END IF
    END IF
    IF (spinup_alpha .AND. model_clock % is_spinning_up()) THEN
      varalpha = 1.0_r_def
    ELSE
      varalpha = alpha
    END IF
    varbeta = 1.0_r_def - varalpha
    CALL check_fields(self % state, cast_dt)
    IF (self % use_moisture) THEN
      CALL copy_bundle(mr, self % mr_n, nummr)
      CALL copy_bundle(mr, self % mr_after_slow, nummr)
    END IF
    CALL copy_bundle(self % state, self % state_n, bundle_size)
    CALL copy_bundle(self % state, self % state_after_slow, bundle_size)
    IF (use_physics) THEN
      IF (self % output_cld_incs) THEN
        CALL cld_incs_init(cloud_fields, dcfl_tot, dcff_tot, dbcf_tot, sec_tot, suffix_tot)
      END IF
      CALL slow_physics(modeldb, self % du, self % dtheta, self % mr_after_slow, self % state_n(igh_t), self % state_n(igh_u), &
&self % state_n(igh_d), self % state_n(igh_p), moist_dyn, self % mr_n, derived_fields, radiation_fields, microphysics_fields, &
&electric_fields, orography_fields, turbulence_fields, convection_fields, cloud_fields, surface_fields, soil_fields, snow_fields, &
&chemistry_fields, aerosol_fields, model_clock, cast_dt, dtemp_encorr, mesh, twod_mesh)
      CALL invoke_update_from_slow_physics(self % state_after_slow(igh_t), self % dtheta, self % state_after_slow(igh_u), self % du)
      IF (self % output_cld_incs) THEN
        CALL cld_incs_init(cloud_fields, dcfl_adv, dcff_adv, dbcf_adv, sec_adv, suffix_adv)
      END IF
    END IF
    CALL rhs_alg(self % rhs_n, varbeta * cast_dt, self % state_after_slow, self % state_n, moist_dyn, compute_eos = .FALSE., &
&compute_rhs_t_d = .TRUE., dlayer_rhs = .FALSE., model_clock = model_clock)
    CALL copy_bundle(self % state_after_slow, self % advected_state, bundle_size)
    mr_to_adv => self % mr_after_slow
    CALL predictors_alg(self % advected_state, self % state_n(igh_u), self % rhs_n(igh_u), varbeta, model_clock)
    CALL conditional_collection_copy(adv_tracer_all_outer_after_slow, generic_fields_to_copy = adv_tracer_all_outer, &
&field_list = adv_tracer_all_outer)
    CALL conditional_collection_copy(adv_tracer_last_outer_after_slow, generic_fields_to_copy = adv_tracer_last_outer, &
&field_list = adv_tracer_last_outer)
    CALL conditional_collection_copy(con_tracer_all_outer_after_slow, generic_fields_to_copy = con_tracer_all_outer, &
&field_list = con_tracer_all_outer)
    CALL conditional_collection_copy(con_tracer_last_outer_after_slow, generic_fields_to_copy = con_tracer_last_outer, &
&field_list = con_tracer_last_outer)
    CALL invoke_4(self % wind_prev, self % state(igh_u))
    outer_dynamics_loop:DO outer = 1, outer_iterations
      IF (use_wavedynamics) THEN
        CALL gungho_transport_control_alg(self % rhs_adv, self % advected_state, self % state(igh_u), self % state_n(igh_u), mr, &
&mr_to_adv, model_clock, outer, cheap_update, self % adv_inc_prev, self % wind_prev, self % state_after_slow(igh_d), &
&self % total_dry_flux, adv_tracer_all_outer, adv_tracer_all_outer_after_slow, adv_tracer_last_outer, &
&adv_tracer_last_outer_after_slow, con_tracer_all_outer, con_tracer_all_outer_after_slow, con_tracer_last_outer, &
&con_tracer_last_outer_after_slow)
        IF (cheap_update .AND. (outer < outer_iterations)) THEN
          CALL copy_bundle(self % rhs_adv, self % adv_inc_prev, bundle_size)
          CALL mass_matrix_solver_alg(self % du, self % rhs_adv(igh_u))
          CALL invoke_5(self % advected_state(igh_d), self % rhs_adv(igh_d), self % advected_state(igh_t), self % rhs_adv(igh_t), &
&self % advected_state(igh_u), self % du)
          CALL adv_tracer_all_outer_after_slow % clear
          CALL con_tracer_all_outer_after_slow % clear
          CALL conditional_collection_copy(adv_tracer_all_outer_after_slow, generic_fields_to_copy = adv_tracer_all_outer, &
&field_list = adv_tracer_all_outer)
          CALL conditional_collection_copy(con_tracer_all_outer_after_slow, generic_fields_to_copy = con_tracer_all_outer, &
&field_list = con_tracer_all_outer)
          IF (self % use_moisture) THEN
            mr_to_adv => self % mr_after_adv
          END IF
          CALL invoke_6(self % wind_prev, self % state(igh_u))
        END IF
        CALL invoke_7(self % dtheta, self % rhs_adv(igh_t), mm_wt)
        CALL rhs_alg(self % rhs_np1, - varalpha * cast_dt, self % state, self % state, moist_dyn, compute_eos = .TRUE., &
&compute_rhs_t_d = .TRUE., dlayer_rhs = dlayer_on, model_clock = model_clock)
      ELSE
        IF (self % use_moisture) CALL copy_bundle(self % mr_after_slow, mr, nummr)
      END IF
      IF (self % use_moisture) CALL copy_bundle(mr, self % mr_after_adv, nummr)
      IF (use_physics) THEN
        IF (blayer_placement == blayer_placement_fast .OR. convection_placement == convection_placement_fast .OR. &
&stochastic_physics_placement == stochastic_physics_placement_fast) THEN
          CALL calc_phys_predictors_alg(derived_fields, self % rhs_np1, self % rhs_adv, self % rhs_n, self % state, &
&self % state_after_slow, lbc_fields, model_clock)
        END IF
        IF (outer == outer_iterations .AND. self % output_cld_incs) THEN
          CALL cld_incs_output(cloud_fields, dcfl_adv, dcff_adv, dbcf_adv, sec_adv, suffix_adv)
        END IF
        CALL fast_physics(self % du, self % dtheta, mr, self % state_n(igh_t), self % state_n(igh_d), self % state_n(igh_u), &
&self % state_n(igh_p), self % mr_n, derived_fields, radiation_fields, microphysics_fields, orography_fields, turbulence_fields, &
&convection_fields, cloud_fields, surface_fields, soil_fields, snow_fields, chemistry_fields, aerosol_fields, stph_fields, outer, &
&model_clock, cast_dt)
        IF (smagorinsky .AND. smagorinsky_placement == smagorinsky_placement_outer) THEN
          CALL smagorinsky_alg(self % dtheta, self % du, mr, self % state(igh_t), self % state(igh_u), derived_fields, &
&self % state(igh_d), cast_dt)
        END IF
        IF (use_wavedynamics) THEN
          CALL set_bundle_scalar(0.0_r_def, self % rhs_phys, bundle_size)
          CALL invoke_update_rhs_phys_from_fast_physics(self % rhs_phys(igh_t), self % dtheta, mm_wt, self % rhs_phys(igh_u), &
&self % du, mm_vel)
        END IF
      END IF
      IF (use_wavedynamics) THEN
        IF (guess_np1) THEN
          IF (self % use_moisture) CALL moist_dyn_factors_alg(moist_dyn, mr)
          CALL update_prognostic_scalars_alg(self % state, self % rhs_n, self % rhs_adv, self % rhs_phys, moist_dyn(gas_law))
        END IF
        inner_dynamics_loop:DO inner = 1, inner_iterations
          WRITE(log_scratch_space, '(A,2I3)') 'loop indices (o, i): ', outer, inner
          CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
          IF (inner > 1) CALL rhs_alg(self % rhs_np1, - varalpha * cast_dt, self % state, self % state, moist_dyn, compute_eos = &
&.TRUE., compute_rhs_t_d = .FALSE., dlayer_rhs = dlayer_on, model_clock = model_clock)
          IF (limited_area .AND. inner == 1 .AND. outer == 1) THEN
            CALL lam_solver_lbc(self % state(igh_u), lbc_fields, prime_mesh_name)
            CALL calc_rhs_lbc(self % rhs_lbc, lbc_fields, model_clock, prime_mesh_name, tau_r)
          END IF
          CALL bundle_axpy(- 1.0_r_def, self % rhs_np1, self % rhs_n, self % rhs_np1, bundle_size)
          CALL add_bundle(self % rhs_np1, self % rhs_adv, self % rhs_np1, bundle_size)
          CALL add_bundle(self % rhs_np1, self % rhs_phys, self % rhs_np1, bundle_size)
          IF (limited_area) THEN
            IF (inner == 1 .AND. outer == 1) THEN
              CALL add_bundle(self % rhs_np1, self % rhs_lbc, self % rhs_np1, bundle_size)
            END IF
            CALL apply_mask_rhs(self % rhs_np1, prime_mesh_name)
          END IF
          IF (inner > 1) THEN
            CALL invoke_9(self % rhs_np1(igh_d), self % rhs_np1(igh_t))
          END IF
          write_moisture_diag = write_conservation_diag .AND. outer == outer_iterations .AND. inner == inner_iterations .AND. self &
&% use_moisture
          CALL semi_implicit_solver_alg_step(self % state, self % rhs_np1, moist_dyn(gas_law), mr, write_moisture_diag, &
&first_iteration = (inner == 1))
          IF (.NOT. guess_np1 .AND. self % use_moisture) CALL moist_dyn_factors_alg(moist_dyn, mr)
          IF (exner_from_eos) THEN
            CALL derive_exner_from_eos(self % state, moist_dyn(gas_law))
          END IF
          IF (limited_area) THEN
            IF ((blend_frequency == blend_frequency_inner) .OR. (blend_frequency == blend_frequency_outer .AND. inner == &
&inner_iterations) .OR. (blend_frequency == blend_frequency_final .AND. inner == inner_iterations .AND. outer == &
&outer_iterations)) THEN
              CALL lam_blend_lbc(self % state(igh_u), self % state(igh_p), self % state(igh_d), self % state(igh_t), mr, &
&lbc_fields, microphysics_fields, aerosol_fields, lbc_option, microphysics_casim, murk_lbc, moisture_formulation, prime_mesh_name)
            END IF
          END IF
        END DO inner_dynamics_loop
      ELSE
        CALL invoke_10(self % state(igh_t), self % state_after_slow(igh_t), self % dtheta, self % state(igh_u), &
&self % state_after_slow(igh_u), self % du)
      END IF
    END DO outer_dynamics_loop
    CALL adv_tracer_all_outer_after_slow % clear
    CALL adv_tracer_last_outer_after_slow % clear
    CALL con_tracer_all_outer_after_slow % clear
    CALL con_tracer_last_outer_after_slow % clear
    IF (transport_ageofair) THEN
      CALL con_tracer_last_outer % get_field('ageofair', ageofair)
      CALL ageofair_update(ageofair, model_clock)
    END IF
    CALL mixing_alg(mr, self % state(igh_t), self % state(igh_u), derived_fields, self % state(igh_d), cast_dt)
    IF (self % use_moisture) THEN
      CALL moist_dyn_factors_alg(moist_dyn, mr)
    END IF
    IF (use_physics) THEN
      CALL map_physics_fields_alg(self % state(igh_u), self % state(igh_p), self % state(igh_d), self % state(igh_t), moist_dyn, &
&derived_fields)
    END IF
    IF (use_physics .AND. self % output_cld_incs) THEN
      CALL cld_incs_output(cloud_fields, dcfl_tot, dcff_tot, dbcf_tot, sec_tot, suffix_tot)
    END IF
    IF (write_diag .AND. use_xios_io .AND. MOD(model_clock % get_step(), diagnostic_frequency) == 0) THEN
      CALL output_diags_for_si(self % state, self % state_n, self % state_after_slow, mr, self % mr_n, self % mr_after_slow, &
&self % mr_after_adv, derived_fields, self % du, self % dtheta, self % dtheta_cld)
    END IF
    CALL invoke_11(u, self % state(igh_u), theta, self % state(igh_t), rho, self % state(igh_d), exner, self % state(igh_p))
    NULLIFY(mm_wt, mm_vel)
    IF (LPROF) CALL stop_timing(id, 'semi_implicit_timestep')
  END SUBROUTINE run_step
  SUBROUTINE semi_implicit_alg_final(self)
    IMPLICIT NONE
    CLASS(semi_implicit_timestep_type), INTENT(INOUT) :: self
    CALL semi_implicit_solver_alg_final
    CALL final_si_operators
    IF (ALLOCATED(self % state)) DEALLOCATE(self % state)
    IF (ALLOCATED(self % state_n)) DEALLOCATE(self % state_n)
    IF (ALLOCATED(self % state_after_slow)) DEALLOCATE(self % state_after_slow)
    IF (ALLOCATED(self % advected_state)) DEALLOCATE(self % advected_state)
    IF (ALLOCATED(self % rhs_n)) DEALLOCATE(self % rhs_n)
    IF (ALLOCATED(self % rhs_np1)) DEALLOCATE(self % rhs_np1)
    IF (ALLOCATED(self % rhs_adv)) DEALLOCATE(self % rhs_adv)
    IF (ALLOCATED(self % rhs_phys)) DEALLOCATE(self % rhs_phys)
    IF (ALLOCATED(self % rhs_lbc)) DEALLOCATE(self % rhs_lbc)
    IF (ALLOCATED(self % adv_inc_prev)) DEALLOCATE(self % adv_inc_prev)
    IF (ALLOCATED(self % mr_n)) DEALLOCATE(self % mr_n)
    IF (ALLOCATED(self % mr_after_slow)) DEALLOCATE(self % mr_after_slow)
    IF (ALLOCATED(self % mr_after_adv)) DEALLOCATE(self % mr_after_adv)
    CALL self % dtheta % field_final
    CALL self % du % field_final
    CALL self % total_dry_flux % field_final
    RETURN
  END SUBROUTINE semi_implicit_alg_final
  SUBROUTINE conditional_collection_copy(generic_fields_copied, generic_fields_to_copy, field_list)
    USE semi_implicit_timestep_alg_mod_psy, ONLY: invoke_12
    USE field_collection_mod, ONLY: field_collection_type
    USE field_collection_iterator_mod, ONLY: field_collection_iterator_type
    IMPLICIT NONE
    TYPE(field_collection_type), INTENT(OUT) :: generic_fields_copied
    TYPE(field_collection_type), INTENT(IN) :: generic_fields_to_copy
    TYPE(field_collection_type), INTENT(IN) :: field_list
    TYPE(field_collection_iterator_type) :: iterator
    CLASS(field_parent_type), POINTER :: abstract_field_ptr
    TYPE(field_type), POINTER :: single_generic_field
    TYPE(field_type) :: copied_generic_field
    LOGICAL(KIND = l_def) :: l_copy
    NULLIFY(abstract_field_ptr)
    NULLIFY(single_generic_field)
    CALL generic_fields_copied % initialise(name = 'fields_copied')
    IF (generic_fields_to_copy % get_length() > 0) THEN
      CALL iterator % initialise(generic_fields_to_copy)
      DO
        IF (.NOT. iterator % has_next()) EXIT
        abstract_field_ptr => iterator % next()
        SELECT TYPE(abstract_field_ptr)
          TYPE IS (field_type)
          single_generic_field => abstract_field_ptr
        END SELECT
        l_copy = field_list % field_exists(single_generic_field % get_name())
        IF (l_copy) THEN
          CALL single_generic_field % copy_field_properties(copied_generic_field)
          CALL invoke_12(copied_generic_field, single_generic_field)
          CALL generic_fields_copied % add_field(copied_generic_field)
        END IF
      END DO
    END IF
  END SUBROUTINE conditional_collection_copy
END MODULE semi_implicit_timestep_alg_mod