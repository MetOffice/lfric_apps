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

  use base_mesh_config_mod,       only : prime_mesh_name
  use constants_mod,              only : r_def, l_def
  use derived_config_mod,         only : l_esm_couple
  use extrusion_mod,              only : TWOD
  use field_collection_mod,       only : field_collection_type
  use formulation_config_mod,     only : use_multires_coupling
  use gungho_diagnostics_driver_mod, &
                                  only : gungho_diagnostics_driver
  use iau_time_control_mod,       only : calc_iau_ts_end
  use gungho_init_fields_mod,     only : create_model_data, &
                                         initialise_model_data, &
                                         output_model_data, &
                                         finalise_model_data
  use gungho_modeldb_mod,         only : modeldb_type
  use gungho_model_mod,           only : initialise_infrastructure, &
                                         initialise_model, &
                                         finalise_infrastructure, &
                                         finalise_model
  use gungho_step_mod,            only : gungho_step
  use gungho_time_axes_mod,       only : gungho_time_axes_type, &
                                         get_time_axes_from_collection
  use initial_output_mod,         only : write_initial_output
  use io_config_mod,              only : write_diag, &
                                         diagnostic_frequency, &
                                         nodal_output_on_w3
  use initialization_config_mod,  only : lbc_option,               &
                                         lbc_option_gungho_file,   &
                                         lbc_option_um2lfric_file, &
                                         ancil_option,             &
                                         ancil_option_updating,    &
                                         coarse_aerosol_ancil
  use init_gungho_lbcs_alg_mod,   only : update_lbcs_file_alg
  use log_mod,                    only : log_event,           &
                                         log_level_always,    &
                                         log_level_error,     &
                                         log_level_info,      &
                                         log_scratch_space
  use mesh_mod,                   only : mesh_type
  use mesh_collection_mod,        only : mesh_collection
  use mpi_mod,                    only : mpi_type
  use multires_coupling_config_mod, &
                                  only : aerosol_mesh_name
  use section_choice_config_mod,  only : iau, &
                                         iau_surf

#ifdef UM_PHYSICS
  use variable_fields_mod,        only : update_variable_fields
  use lfric_xios_time_axis_mod,   only : regridder
  use intermesh_mappings_alg_mod, only : map_scalar_intermesh
  use update_ancils_alg_mod,      only : update_ancils_alg
  use gas_calc_all_mod,           only : gas_calc_all
  use multires_coupling_config_mod, &
                                  only : lowest_order_aero_flag
  use update_iau_inc_alg_mod,     only : update_iau_alg
  use update_iau_surf_alg_mod,    only : add_surf_inc_alg
  use iau_config_mod,             only : iau_mode, &
                                         iau_mode_instantaneous, &
                                         iau_mode_time_window, &
                                         iau_ts_start
  use section_choice_config_mod,  only: stochastic_physics, &
                                        stochastic_physics_um
  use stochastic_physics_config_mod, &
                                  only: use_random_parameters
  use stph_rp_main_alg_mod,       only: stph_rp_main_alg, &
                                        stph_rp_init_alg
#endif
#ifdef COUPLED
  use esm_couple_config_mod,      only : l_esm_couple_test
  use coupler_mod,                only : cpl_snd, cpl_rcv, cpl_fld_update
  use coupler_diagnostics_mod,    only : save_sea_ice_frac_previous
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

    character(len=*), parameter :: io_context_name = "gungho_atm"

#ifdef UM_PHYSICS
    ! For clearing IAU fields after use
    type( field_collection_type ), pointer :: field_collection_ptr => null()
#endif

    ! Initialise infrastructure and setup constants
    call initialise_infrastructure( io_context_name, modeldb )

    ! Add a place to store time axes in modeldb
    call modeldb%values%add_key_value('model_axes', model_axes)

    ! Get primary and 2D meshes for initialising model data
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)

    ! If aerosol data is on a different mesh, get this
    if (coarse_aerosol_ancil) then
      ! For now use the coarsest mesh
      aerosol_mesh => mesh_collection%get_mesh(aerosol_mesh_name)
      aerosol_twod_mesh => mesh_collection%get_mesh(aerosol_mesh, TWOD)
    else
      aerosol_mesh => mesh
      aerosol_twod_mesh => twod_mesh
    end if

    ! Rate of temperature adjustment for energy correction
    call modeldb%values%add_key_value( 'temperature_correction_rate', &
                                       0.0_r_def )
    ! Total mass of dry atmosphere used for energy correction
    call modeldb%values%add_key_value( 'total_dry_mass', 0.0_r_def )
    ! Total energy of moist atmosphere for calculating energy correction
    call modeldb%values%add_key_value( 'total_energy', 0.0_r_def )
    ! Total energy of moist atmosphere at previous energy correction step
    call modeldb%values%add_key_value( 'total_energy_previous', 0.0_r_def )

    ! Instantiate the fields stored in model_data
    call create_model_data( modeldb,         &
                            mesh, twod_mesh, &
                            aerosol_mesh, aerosol_twod_mesh )

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
    ! state, then do this now
    if ( ( iau ) .and. ( iau_mode == iau_mode_instantaneous ) ) then

      call update_iau_alg( modeldb, twod_mesh )
      field_collection_ptr => modeldb%fields%get_field_collection("iau_fields")
      call field_collection_ptr%clear()

    end if

    if ( iau_surf ) then
      field_collection_ptr => modeldb%fields%get_field_collection("iau_surf_fields")
      call add_surf_inc_alg( field_collection_ptr,              &
                             modeldb%model_data%surface_fields, &
                             modeldb%model_data%soil_fields,    &
                             modeldb%model_data%snow_fields )

      call field_collection_ptr%clear()

    end if

    ! Initialise RP scheme (stochastic perturbed parameters)
    if ( stochastic_physics == stochastic_physics_um .and. &
         use_random_parameters ) then
      call stph_rp_init_alg( modeldb )
    end if
#endif

    nullify(mesh, twod_mesh, aerosol_mesh, aerosol_twod_mesh)

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

#ifdef COUPLED
    type( field_collection_type ), pointer :: depository => null()
#endif

    type(gungho_time_axes_type), pointer :: model_axes

#ifdef UM_PHYSICS
    procedure(regridder), pointer :: regrid_operation => null()
    logical(l_def)                :: regrid_lowest_order

    ! For clearing IAU fields after use
    type( field_collection_type ), pointer :: field_collection_ptr => null()

    regrid_operation => map_scalar_intermesh
    if (use_multires_coupling) then
      regrid_lowest_order = lowest_order_aero_flag
    else
      regrid_lowest_order = .false.
    end if
#endif
    ! Get model_axes out of modeldb
    model_axes => get_time_axes_from_collection(modeldb%values, "model_axes" )

    ! Get primary and 2D meshes
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)

#ifdef COUPLED
    if(l_esm_couple) then

       write(log_scratch_space, &
             '(A, I0)') 'Coupling timestep: ', modeldb%clock%get_step() - 1
       call log_event( log_scratch_space, LOG_LEVEL_INFO )

       depository => modeldb%fields%get_field_collection("depository")

       call save_sea_ice_frac_previous(depository)

       ! Receive all incoming (ocean/seaice fields) from the coupler
       call cpl_rcv( modeldb%model_data%cpl_rcv,    &
                     depository, &
                     modeldb%clock )

       ! Send all outgoing (ocean/seaice driving fields) to the coupler
       call cpl_snd( modeldb%model_data%cpl_snd,    &
                     depository, &
                     modeldb%clock )

    endif
#endif

    if ( lbc_option == lbc_option_gungho_file .or. &
         lbc_option == lbc_option_um2lfric_file) then

      call update_lbcs_file_alg( modeldb%configuration,     &
                                 model_axes%lbc_times_list, &
                                 modeldb%clock, modeldb%model_data%lbc_fields )
    endif

#ifdef UM_PHYSICS
    ! If IAU is active and increments need to be added over a time window, then do this
    ! at the start of every ts within the required time window
    if ( ( iau ) .and. ( iau_mode == iau_mode_time_window ) ) then

      if ( ( modeldb%clock%get_step() >= iau_ts_start ) .and.                 &
           ( modeldb%clock%get_step() <= calc_iau_ts_end( modeldb%clock ) ) ) then

        call update_iau_alg( modeldb, twod_mesh )

      else if ( modeldb%clock%get_step() > calc_iau_ts_end( modeldb%clock ) ) then

        field_collection_ptr => modeldb%fields%get_field_collection("iau_fields")
        call field_collection_ptr%clear()

      end if

    end if

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
      call update_variable_fields( model_axes%ancil_times_list, &
                                   modeldb%clock,                       &
                                   modeldb%model_data%ancil_fields,     &
                                   regrid_operation,                    &
                                   regrid_lowest_order )
      call update_ancils_alg( model_axes%ancil_times_list, &
                              modeldb%clock,                       &
                              modeldb%model_data%ancil_fields,     &
                              modeldb%model_data%surface_fields)
    end if

    ! Update the time varying trace gases
    call gas_calc_all()
#endif

    nullify(mesh, twod_mesh)

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

#ifdef COUPLED
    type( field_collection_type ), pointer :: depository => null()

    if (l_esm_couple) then
       depository => modeldb%fields%get_field_collection("depository")
       ! Ensure coupling fields are updated at the end of a cycle to ensure the values
       ! stored in and recovered from checkpoint dumps are correct and reproducible
       ! when (re)starting subsequent runs!
       call cpl_fld_update(modeldb%model_data%cpl_snd,    &
                           depository,                    &
                           modeldb%clock)
    endif
#endif

    ! Write out the model state
    call output_model_data( modeldb )

    ! Model configuration finalisation
    call finalise_model( modeldb,               &
                         program_name )

    ! Destroy the fields stored in model_data
    call finalise_model_data( modeldb )

    ! Finalise infrastructure and constants
    call finalise_infrastructure(modeldb)

  end subroutine finalise

end module gungho_driver_mod
