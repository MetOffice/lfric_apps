!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>@brief Drives the execution of the (tangent) linear model.
module linear_driver_mod

  use constants_mod,              only : i_def, r_def, imdi, l_def, str_def, &
                                         i_medium
  use extrusion_mod,              only : TWOD
  use field_array_mod,            only : field_array_type
  use field_mod,                  only : field_type
  use field_collection_mod,       only : field_collection_type
  use io_value_mod,               only : io_value_type
  use section_choice_config_mod,  only : stochastic_physics, &
                                         stochastic_physics_um
  use gungho_diagnostics_driver_mod, &
                                  only : gungho_diagnostics_driver
  use gungho_model_mod,           only : initialise_infrastructure, &
                                         initialise_model, &
                                         finalise_infrastructure, &
                                         finalise_model
  use gungho_init_fields_mod,     only : create_model_data, &
                                         initialise_model_data, &
                                         output_model_data, &
                                         finalise_model_data
  use driver_modeldb_mod,         only : modeldb_type
  use gungho_step_mod,            only : gungho_step
  use gungho_time_axes_mod,       only : gungho_time_axes_type, &
                                         get_time_axes_from_collection
  use init_fd_prognostics_mod,    only : init_fd_prognostics_dump
  use initial_output_mod,         only : write_initial_output
  use initialization_config_mod,  only : init_option_fd_start_dump, &
                                         ls_option_file
  use io_context_mod,             only : io_context_type
  use log_mod,                    only : log_event,         &
                                         log_scratch_space, &
                                         LOG_LEVEL_ALWAYS,  &
                                         LOG_LEVEL_INFO
  use linear_model_data_mod,      only : linear_create_ls, &
                                         linear_init_ls,   &
                                         linear_init_pert
  use linear_model_mod,           only : initialise_linear_model, &
                                         finalise_linear_model
  use linear_diagnostics_driver_mod, &
                                  only : linear_diagnostics_driver
  use linear_step_mod,            only : linear_step
  use linear_data_algorithm_mod,  only : update_ls_file_alg
  use mesh_mod,                   only : mesh_type
  use mesh_collection_mod,        only : mesh_collection
  use namelist_mod,               only : namelist_type
  use create_tl_prognostics_mod,  only : create_tl_prognostics

  implicit none

  private
  public initialise, step, finalise

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up required state in preparation for run.
  !> @param [in]     program_name An identifier given to the model being run
  !> @param [in,out] modeldb   The structure that holds model state
  !>
  subroutine initialise( program_name, modeldb )

    implicit none

    character(*),         intent(in)    :: program_name
    type(modeldb_type),   intent(inout) :: modeldb

    type(gungho_time_axes_type)     :: model_axes


    type( mesh_type ),     pointer :: mesh
    type( mesh_type ),     pointer :: twod_mesh
    type( mesh_type ),     pointer :: aerosol_mesh
    type( mesh_type ),     pointer :: aerosol_twod_mesh

    type( namelist_type ), pointer :: base_mesh_nml
    type( namelist_type ), pointer :: multires_coupling_nml
    type( namelist_type ), pointer :: initialization_nml
    type( namelist_type ), pointer :: io_nml

    character( len=str_def )       :: prime_mesh_name
    character( len=str_def )       :: aerosol_mesh_name
    logical( kind=l_def )          :: coarse_aerosol_ancil
    logical( kind=l_def )          :: coarse_ozone_ancil
    integer( kind=i_def )          :: init_option
    logical( kind=l_def )          :: nodal_output_on_w3

    type( io_value_type ) :: temp_corr_io_value
    type( io_value_type ) :: random_seed_io_value
    type( field_collection_type ), pointer :: depository
    type( field_collection_type ), pointer :: fd_fields

    character(len=*), parameter :: io_context_name = "gungho_atm"
    integer(i_def) :: random_seed_size
    real(r_def), allocatable :: real_array(:)

    nullify( mesh, twod_mesh, aerosol_mesh, aerosol_twod_mesh, depository )
    nullify( base_mesh_nml, multires_coupling_nml, initialization_nml )

    depository => modeldb%fields%get_field_collection("depository")
    fd_fields => modeldb%fields%get_field_collection("fd_fields")

    call temp_corr_io_value%init('temperature_correction_rate', [0.0_r_def])
    call modeldb%values%add_key_value( 'temperature_correction_io_value', &
                                       temp_corr_io_value )
    call modeldb%values%add_key_value( 'total_dry_mass', 0.0_r_def )
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

    ! Initialise infrastructure and setup constants
    call initialise_infrastructure( io_context_name, modeldb )

    ! Add a place to store time axes in modeldb
    call modeldb%values%add_key_value('model_axes', model_axes)

    ! Get primary and 2D meshes for initialising model data
    base_mesh_nml => modeldb%configuration%get_namelist('base_mesh')
    call base_mesh_nml%get_value( 'prime_mesh_name', prime_mesh_name )
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)

    ! Get initialization configuration
    initialization_nml => modeldb%configuration%get_namelist('initialization')
    call initialization_nml%get_value( 'coarse_aerosol_ancil', &
                                        coarse_aerosol_ancil )
    call initialization_nml%get_value( 'coarse_ozone_ancil', &
                                       coarse_ozone_ancil )
    call initialization_nml%get_value( 'init_option', init_option )

    ! If aerosol data is on a different mesh, get this
    if (coarse_aerosol_ancil .or. coarse_ozone_ancil) then
      ! For now use the coarsest mesh
      multires_coupling_nml => &
        modeldb%configuration%get_namelist('multires_coupling')
      call multires_coupling_nml%get_value( 'aerosol_mesh_name', &
                                            aerosol_mesh_name )
      aerosol_mesh => mesh_collection%get_mesh(aerosol_mesh_name)
      aerosol_twod_mesh => mesh_collection%get_mesh(aerosol_mesh, TWOD)
      write( log_scratch_space,'(A,A)' ) "aerosol mesh name:", aerosol_mesh%get_mesh_name()
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
    else
      aerosol_mesh => mesh
      aerosol_twod_mesh => twod_mesh
    end if

    ! Instantiate the fields stored in model_data
    call create_model_data( modeldb,      &
                            mesh,         &
                            twod_mesh,    &
                            aerosol_mesh, &
                            aerosol_twod_mesh )

    ! Instantiate the fields required to read the initial
    ! conditions from a file.
    if ( init_option == init_option_fd_start_dump ) then
      call create_tl_prognostics( mesh, twod_mesh,   &
                                  fd_fields,         &
                                  depository)
    end if

    ! Instantiate the linearisation state
    call linear_create_ls( modeldb, mesh, twod_mesh)

    ! Initialise the fields stored in the model_data
    if ( init_option == init_option_fd_start_dump ) then
      call init_fd_prognostics_dump( fd_fields )
    else
      call initialise_model_data( modeldb, mesh, twod_mesh )
    end if

    ! Model configuration initialisation
    call initialise_model( mesh, &
                           modeldb )

    ! Initialise the linearisation state
    call linear_init_ls( mesh, twod_mesh, modeldb )

    ! Initialise the linear model perturbation state
    call linear_init_pert( mesh,      &
                           twod_mesh, &
                           modeldb )

    ! Get io configuration
    io_nml => modeldb%configuration%get_namelist('io')
    call io_nml%get_value( 'nodal_output_on_w3', nodal_output_on_w3 )

    ! Initial output
    call write_initial_output( modeldb, mesh, twod_mesh, &
                               io_context_name, nodal_output_on_w3 )

    ! Linear model configuration initialisation
    call initialise_linear_model( mesh,        &
                                  modeldb )

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Timestep the model, calling the desired timestepping algorithm
  !         based upon the configuration
  !> @param [in,out] model_data   The structure that holds model state
  subroutine step( modeldb )

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    type(field_collection_type), pointer :: moisture_fields
    type(field_collection_type), pointer :: ls_fields
    type(field_array_type), pointer      :: ls_mr_array
    type(field_array_type), pointer      :: ls_moist_dyn_array

    type( gungho_time_axes_type ), pointer :: model_axes
    type( mesh_type ),             pointer :: mesh
    type( mesh_type ),             pointer :: twod_mesh
    type( namelist_type ),         pointer :: base_mesh_nml
    type( namelist_type ),         pointer :: initialization_nml
    type( namelist_type ),         pointer :: io_nml
    character( len=str_def )               :: prime_mesh_name
    integer( kind=i_def )                  :: ls_option
    logical( kind=l_def )                  :: write_diag
    integer( kind=i_medium )               :: diagnostic_frequency
    logical( kind=l_def )                  :: nodal_output_on_w3

    nullify(mesh, twod_mesh, base_mesh_nml, initialization_nml, io_nml)
    nullify(moisture_fields, ls_mr_array, ls_moist_dyn_array)

    ! Get model_axes out of modeldb
    model_axes => get_time_axes_from_collection(modeldb%values, "model_axes" )

    moisture_fields => modeldb%fields%get_field_collection("moisture_fields")
    call moisture_fields%get_field("ls_mr", ls_mr_array)
    call moisture_fields%get_field("ls_moist_dyn", ls_moist_dyn_array)

    ls_fields => modeldb%fields%get_field_collection("ls_fields")

    ! Get initialization configuration
    initialization_nml => modeldb%configuration%get_namelist('initialization')
    call initialization_nml%get_value( 'ls_option', ls_option )

    if ( ls_option == ls_option_file ) then
      call update_ls_file_alg( model_axes%ls_times_list,     &
                               modeldb%clock,                &
                               ls_fields,                    &
                               ls_mr_array%bundle,           &
                               ls_moist_dyn_array%bundle )
    end if

    ! Get Mesh
    base_mesh_nml => modeldb%configuration%get_namelist('base_mesh')
    call base_mesh_nml%get_value( 'prime_mesh_name', prime_mesh_name )
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)

    call linear_step( mesh, twod_mesh, &
                      modeldb, modeldb%clock )

    ! Get io configuration
    io_nml => modeldb%configuration%get_namelist('io')
    call io_nml%get_value( 'diagnostic_frequency', diagnostic_frequency )
    call io_nml%get_value( 'write_diag', write_diag )
    call io_nml%get_value( 'nodal_output_on_w3', nodal_output_on_w3 )

    if ( ( mod(modeldb%clock%get_step(), diagnostic_frequency) == 0 ) &
         .and. ( write_diag ) ) then

      ! Calculation and output diagnostics
      call gungho_diagnostics_driver( modeldb, mesh, twod_mesh, &
                                      nodal_output_on_w3 )

      call linear_diagnostics_driver( mesh,    &
                                      modeldb, &
                                      nodal_output_on_w3 )
    end if

  end subroutine step

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Tidies up after a run.
  !> @param [in]     program_name An identifier given to the model being run
  !> @param [in,out] modeldb   The structure that holds model state
  subroutine finalise( program_name, modeldb )

    implicit none

    character(*), intent(in)          :: program_name
    type(modeldb_type), intent(inout) :: modeldb

    ! Write out the model state
    call output_model_data( modeldb )

    ! Model configuration finalisation
    call finalise_model( modeldb, &
                         program_name )

    call finalise_linear_model( )

    ! Destroy the fields stored in model_data
    call finalise_model_data( modeldb )

    ! Finalise infrastructure and constants
    call finalise_infrastructure(modeldb)

  end subroutine finalise

end module linear_driver_mod
