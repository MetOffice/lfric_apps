!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Drives the execution of the GungHo model.
!>
!> This is a temporary solution until we have a proper driver layer.
!>
module gungho_driver_mod

  use base_mesh_config_mod,        only : prime_mesh_name
  use constants_mod,               only : r_def, l_def, str_def, i_def
  use derived_config_mod,          only : l_esm_couple
  use extrusion_mod,               only : TWOD
  use field_collection_mod,        only : field_collection_type
  use formulation_config_mod,      only : use_multires_coupling
  use gungho_diagnostics_driver_mod, &
                                   only : gungho_diagnostics_driver
  use iau_time_control_mod,        only : calc_iau_ts_end
  use gungho_init_fields_mod,      only : create_model_data, &
                                          initialise_model_data, &
                                          output_model_data, &
                                          finalise_model_data
  use driver_modeldb_mod,          only : modeldb_type
  use gungho_model_mod,            only : initialise_infrastructure, &
                                          initialise_model, &
                                          finalise_infrastructure, &
                                          finalise_model
  use gungho_step_mod,             only : gungho_step
  use gungho_time_axes_mod,        only : gungho_time_axes_type, &
                                          get_time_axes_from_collection
  use iau_multifile_io_mod,        only : init_multifile_io,       &
                                          setup_step_multifile_io, &
                                          finalise_multifile_io
  use initial_output_mod,          only : write_initial_output
  use io_config_mod,               only : checkpoint_read,      &
                                          write_diag,           &
                                          diagnostic_frequency, &
                                          multifile_io,         &
                                          nodal_output_on_w3
  use initialization_config_mod,   only : lbc_option,               &
                                          lbc_option_gungho_file,   &
                                          lbc_option_um2lfric_file, &
                                          ancil_option,             &
                                          ancil_option_updating,    &
                                          coarse_aerosol_ancil,     &
                                          coarse_ozone_ancil
  use init_gungho_lbcs_alg_mod,    only : update_lbcs_file_alg
  use log_mod,                     only : log_event,           &
                                          log_level_always,    &
                                          log_level_error,     &
                                          log_level_info,      &
                                          log_scratch_space
  use mesh_mod,                    only : mesh_type
  use mesh_collection_mod,         only : mesh_collection
  use multires_coupling_config_mod, &
                                   only : aerosol_mesh_name
  use remove_field_collection_mod, only : remove_field_collection
  use section_choice_config_mod,   only : iau,                   &
                                          iau_sst,               &
                                          iau_surf,              &
                                          stochastic_physics,    &
                                          stochastic_physics_um
  use io_value_mod,                only : io_value_type
  use time_config_mod,             only : timestep_start
  use timing_mod,                  only : start_timing, stop_timing, &
                                          tik, LPROF

#ifdef UM_PHYSICS
  use variable_fields_mod,         only : update_variable_fields
  use lfric_xios_time_axis_mod,    only : regridder
  use intermesh_mappings_alg_mod,  only : map_scalar_intermesh
  use update_ancils_alg_mod,       only : update_ancils_alg
  use gas_calc_all_mod,            only : gas_calc_all
  use multires_coupling_config_mod, &
                                   only : lowest_order_aero_flag
  use update_iau_inc_alg_mod,      only : update_iau_alg
  use update_iau_sst_alg_mod,      only : update_iau_sst_alg
  use update_iau_surf_alg_mod,     only : add_surf_inc_alg
  use iau_main_alg_mod,            only : iau_main_alg
  use iau_config_mod,              only : iau_mode,               &
                                          iau_mode_instantaneous, &
                                          iau_mode_time_mixed
  use stochastic_physics_config_mod, &
                                   only : use_random_parameters
  use stph_rp_main_alg_mod,        only : stph_rp_main_alg
  use flux_calc_all_mod,           only : flux_calc_init, &
                                          flux_calc_step
  use update_tile_temperature_alg_mod, &
                                   only : update_tile_temperature_alg
  use blayer_config_mod,           only : flux_bc_opt,                   &
                                          flux_bc_opt_specified_scalars, &
                                          flux_bc_opt_specified_scalars_tstar
  use specified_surface_config_mod, &
                                   only : function_name_sst,               &
                                          function_name_sst_time_interpolated
#endif
#ifdef COUPLED
  use esm_couple_config_mod,       only : l_esm_couple_test
  use coupler_mod,                 only : cpl_snd, cpl_rcv, cpl_fld_update
  use process_send_fields_2d_mod,  only : save_sea_ice_frac_previous
#endif

  implicit none

  private
  public initialise, step, finalise

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up required state in preparation for run.
  !> @param [in]     program_name An identifier given to the model being run
  !> @param [in,out] modeldb   The structure that holds model state
  subroutine initialise( program_name, modeldb )

    implicit none

    character(*),         intent(in)    :: program_name
    type(modeldb_type),   intent(inout) :: modeldb

    type(gungho_time_axes_type)     :: model_axes

    type(mesh_type),        pointer :: mesh              => null()
    type(mesh_type),        pointer :: twod_mesh         => null()
    type(mesh_type),        pointer :: aerosol_mesh      => null()
    type(mesh_type),        pointer :: aerosol_twod_mesh => null()

    type(io_value_type) :: temp_corr_io_value
    type(io_value_type) :: random_seed_io_value

    character(len=*), parameter :: io_context_name = "gungho_atm"
    integer(i_def) :: random_seed_size
    real(r_def), allocatable :: real_array(:)
    integer(tik)   :: id

#ifdef UM_PHYSICS
    type( field_collection_type ), pointer :: field_collection_ptr
    type( field_collection_type ), pointer :: soil_fields
    type( field_collection_type ), pointer :: snow_fields
    type( field_collection_type ), pointer :: surface_fields

    nullify( field_collection_ptr, soil_fields, snow_fields, surface_fields )
#endif

    call log_event('Initialising gungho', LOG_LEVEL_INFO)
    if ( LPROF ) call start_timing(id, 'gungho_driver.initialise')
    ! Initialise infrastructure and setup constants
    call initialise_infrastructure( io_context_name, modeldb )

    ! Add a place to store time axes in modeldb
    call modeldb%values%add_key_value('model_axes', model_axes)

    ! Get primary and 2D meshes for initialising model data
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)

    ! If aerosol data is on a different mesh, get this
    if (coarse_aerosol_ancil .or. coarse_ozone_ancil) then
      ! For now use the coarsest mesh
      aerosol_mesh => mesh_collection%get_mesh(aerosol_mesh_name)
      aerosol_twod_mesh => mesh_collection%get_mesh(aerosol_mesh, TWOD)
    else
      aerosol_mesh => mesh
      aerosol_twod_mesh => twod_mesh
    end if

    ! Rate of temperature adjustment for energy correction
    call temp_corr_io_value%init("temperature_correction_rate", [0.0_r_def])
    call modeldb%values%add_key_value( 'temperature_correction_io_value', &
                                       temp_corr_io_value)
    ! Total mass of dry atmosphere used for energy correction
    call modeldb%values%add_key_value( 'total_dry_mass', 0.0_r_def )
    ! Total energy of moist atmosphere for calculating energy correction
    call modeldb%values%add_key_value( 'total_energy', 0.0_r_def )
    ! Total energy of moist atmosphere at previous energy correction step
    call modeldb%values%add_key_value( 'total_energy_previous', 0.0_r_def )
    if ( stochastic_physics == stochastic_physics_um ) then
      ! Random seed for stochastic physics
      call random_seed(size = random_seed_size)
      allocate(real_array(random_seed_size))
      real_array(1:random_seed_size) = 0.0_r_def
      call random_seed_io_value%init("random_seed", real_array)
      call modeldb%values%add_key_value( 'random_seed_io_value', &
                                         random_seed_io_value )
      deallocate(real_array)
    end if

    ! Instantiate the fields stored in model_data
    call create_model_data( modeldb,         &
                            mesh, twod_mesh, &
                            aerosol_mesh, aerosol_twod_mesh )

    ! Set up io for multifile reading
    if ( multifile_io ) then
      call init_multifile_io( io_context_name, modeldb )
    end if

    ! Initialise the fields stored in the model_data
    call initialise_model_data( modeldb, mesh, twod_mesh )

    ! Initial output
    call write_initial_output( modeldb, mesh, twod_mesh, &
                               io_context_name, nodal_output_on_w3 )

    ! Model configuration initialisation
    call initialise_model( mesh, modeldb )

#ifdef COUPLED
    ! Placeholder for ESM coupling initialisation code.
    ! Check we have a value for related namelist control variable
    write(log_scratch_space,'(A,L1)') program_name//': Couple flag l_esm_couple_test: ', &
                                     l_esm_couple_test
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
#endif

#ifdef UM_PHYSICS
    ! If IAU is active and increments need to be added instantaneously, to the initial
    ! state, then do this now. The IAU should not be activated at this stage in
    ! the case of a checkpoint-restart.
    if ( ( iau ) .and.                               &
       ( ( iau_mode == iau_mode_instantaneous ) .OR. &
         ( iau_mode == iau_mode_time_mixed ) ) ) then
      if ( .not. checkpoint_read ) then
        call update_iau_alg( modeldb,                     &
                             twod_mesh,                   &
                             iau_ainc_active = .true.,    &
                             iau_addinf_active = .false., &
                             iau_bcorr_active = .false.,  &
                             iau_pertinc_active = .false. )
      end if

      ! IAU increment fields can now be cleared from the depository
      call remove_field_collection( modeldb, "iau_fields" )

    end if

    if ( iau_surf ) then
      field_collection_ptr => modeldb%fields%get_field_collection("iau_surf_fields")
      soil_fields => modeldb%fields%get_field_collection("soil_fields")
      snow_fields => modeldb%fields%get_field_collection("snow_fields")
      surface_fields => modeldb%fields%get_field_collection("surface_fields")

      if ( .not. checkpoint_read ) then
        call add_surf_inc_alg( field_collection_ptr, &
                               surface_fields,       &
                               soil_fields,          &
                               snow_fields )
      end if

      ! IAU surface increment fields can now be cleared from the depository
      call remove_field_collection( modeldb, "iau_surf_fields" )

    end if

    if ( iau_sst ) then
      field_collection_ptr => modeldb%fields%get_field_collection("iau_sst_fields")
      surface_fields => modeldb%fields%get_field_collection("surface_fields")

      if ( .not. checkpoint_read ) then
        call update_iau_sst_alg( field_collection_ptr, &
                                 surface_fields )
      end if

      ! IAU sst increment fields can now be cleared from the depository
      call remove_field_collection( modeldb, "iau_sst_fields" )

    end if

    ! Specified sensible and latent heat fluxes
    if ( flux_bc_opt == flux_bc_opt_specified_scalars .or. &
         flux_bc_opt == flux_bc_opt_specified_scalars_tstar ) then
      call flux_calc_init( modeldb )
    end if
#endif

    nullify(mesh, twod_mesh, aerosol_mesh, aerosol_twod_mesh)
    if ( LPROF ) call stop_timing(id, 'gungho_driver.initialise')

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Timestep the model, calling the desired timestepping algorithm
  !>        based upon the configuration.
  !> @param [in,out] modeldb   The structure that holds model state
  !>
  subroutine step( modeldb )

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    type(mesh_type), pointer :: mesh      => null()
    type(mesh_type), pointer :: twod_mesh => null()
    integer(kind=i_def)      :: ts_start, rc
    integer(tik)             :: tid_first, tid_rest

#if defined(COUPLED) || defined(UM_PHYSICS)
    type( field_collection_type ), pointer :: depository => null()
#endif

    type( field_collection_type ), pointer :: lbc_fields

    type(gungho_time_axes_type), pointer :: model_axes
    character(len=*), parameter :: io_context_name = "gungho_atm"

#ifdef UM_PHYSICS
    procedure(regridder), pointer :: regrid_operation => null()
    logical(l_def)                :: regrid_lowest_order

    type( field_collection_type ), pointer :: surface_fields
    type( field_collection_type ), pointer :: ancil_fields
#endif

    read(timestep_start,*,iostat=rc)  ts_start
    if ( LPROF ) then
      if ( modeldb%clock%get_step() == ts_start ) then
        call start_timing(tid_first, 'gungho_driver.first_timestep')
      else
        call start_timing(tid_rest, 'gungho_driver.timestep')
      end if
    end if
#ifdef UM_PHYSICS
    nullify( surface_fields, ancil_fields )

    regrid_operation => map_scalar_intermesh
    if (use_multires_coupling) then
      regrid_lowest_order = lowest_order_aero_flag
    else
      regrid_lowest_order = .false.
    end if

    ! Step multifile io if active. This sets up and shuts down XIOS contexts at the
    ! required time during the model run
    if( multifile_io ) then
      call setup_step_multifile_io( io_context_name, modeldb )
    end if

#endif
    ! Get model_axes out of modeldb
    model_axes => get_time_axes_from_collection(modeldb%values, "model_axes" )

    ! Get primary and 2D meshes
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)

#ifdef UM_PHYSICS
    ! Specified sensible and latent heat fluxes
    if ( flux_bc_opt == flux_bc_opt_specified_scalars .or. &
         flux_bc_opt == flux_bc_opt_specified_scalars_tstar ) then
      call flux_calc_step( modeldb )
    end if

    ! Specified surface temperatures
    if ( function_name_sst == function_name_sst_time_interpolated ) then
      surface_fields => modeldb%fields%get_field_collection("surface_fields")
      call update_tile_temperature_alg ( modeldb%clock, &
                                         surface_fields )
    end if
#endif

    lbc_fields => modeldb%fields%get_field_collection("lbc_fields")

#ifdef COUPLED
    if(l_esm_couple) then

       write(log_scratch_space, &
             '(A, I0)') 'Coupling timestep: ', modeldb%clock%get_step() - 1
       call log_event( log_scratch_space, LOG_LEVEL_INFO )

       depository => modeldb%fields%get_field_collection("depository")
       call save_sea_ice_frac_previous(depository)

       ! Receive all incoming (ocean/seaice fields) from the coupler
       call cpl_rcv( modeldb )

       ! Send all outgoing (ocean/seaice driving fields) to the coupler
       call cpl_snd( modeldb )

    endif
#endif

    if ( lbc_option == lbc_option_gungho_file .or. &
         lbc_option == lbc_option_um2lfric_file) then

      call update_lbcs_file_alg( modeldb%configuration,     &
                                 model_axes%lbc_times_list, &
                                 modeldb%clock, lbc_fields )
    endif

#ifdef UM_PHYSICS
    ! If IAU is active and increments need to be added over a time window, then do this
    ! at the start of every ts within the required time window
    if ( iau ) then
      call iau_main_alg( modeldb, twod_mesh )
    end if ! (iau)

    ! Apply RP scheme (stochastic perturbed parameters)
    if ( stochastic_physics == stochastic_physics_um .and. &
         use_random_parameters ) then
      call stph_rp_main_alg( modeldb )
    end if
#endif

    ! Perform a timestep
    call gungho_step( mesh, twod_mesh, modeldb, modeldb%clock )

    ! Use diagnostic output frequency to determine whether to write
    ! diagnostics on this timestep

    if ( ( mod(modeldb%clock%get_step(), diagnostic_frequency) == 0 ) &
         .and. ( write_diag ) ) then

      ! Calculation and output diagnostics
      call gungho_diagnostics_driver( modeldb, mesh, twod_mesh, &
                                      nodal_output_on_w3 )
    end if

#ifdef UM_PHYSICS
    ! Update time-varying ancillaries
    ! This is done last in the timestep, because the time data of the
    ! ancillaries needs to be that of the start of timestep, but the
    ! first thing that happens in a timestep is that the clock ticks to the
    ! end of timestep date.
    if ( ancil_option == ancil_option_updating ) then
      surface_fields => modeldb%fields%get_field_collection("surface_fields")
      ancil_fields => modeldb%fields%get_field_collection("ancil_fields")
      call update_variable_fields( model_axes%ancil_times_list, &
                                   modeldb%clock,                       &
                                   ancil_fields,                        &
                                   regrid_operation,                    &
                                   regrid_lowest_order )
      call update_ancils_alg( model_axes%ancil_times_list, &
                              modeldb%clock,                       &
                              ancil_fields,                        &
                              surface_fields)
    end if

    ! Update the time varying trace gases
    call gas_calc_all()
#endif

    ! Write out the model state
    call output_model_data( modeldb )

    nullify(mesh, twod_mesh)

    if ( LPROF ) then
      if ( modeldb%clock%get_step() == ts_start ) then
        call stop_timing(tid_first, 'gungho_driver.first_timestep')
      else
        call stop_timing(tid_rest, 'gungho_driver.timestep')
      end if
    end if

  end subroutine step

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Tidies up after a run.
  !> @param [in]     program_name An identifier given to the model begin run
  !> @param [in,out] modeldb      The structure that holds model state
  !>
  subroutine finalise( program_name, modeldb )

    implicit none

    character(*), intent(in)          :: program_name
    type(modeldb_type), intent(inout) :: modeldb
    integer(tik)                      :: id

#ifdef COUPLED
    type( field_collection_type ), pointer :: depository => null()
#endif

    if ( LPROF ) call start_timing(id, 'gungho_driver.finalise')

#ifdef COUPLED
    if (l_esm_couple) then
       depository => modeldb%fields%get_field_collection("depository")
       ! Ensure coupling fields are updated at the end of a cycle to ensure the values
       ! stored in and recovered from checkpoint dumps are correct and reproducible
       ! when (re)starting subsequent runs!
       call cpl_fld_update(modeldb)
    endif
#endif

    ! Multifile io finalisation
    if( multifile_io ) then
      call finalise_multifile_io( modeldb)
    end if

    ! Model configuration finalisation
    call finalise_model( modeldb,               &
                         program_name )

    ! Destroy the fields stored in model_data
    call finalise_model_data( modeldb )

    ! Finalise infrastructure and constants
    call finalise_infrastructure(modeldb)

    call log_event('gungho finalised', LOG_LEVEL_INFO)
    if ( LPROF ) call stop_timing(id, 'gungho_driver.finalise')

  end subroutine finalise

end module gungho_driver_mod
