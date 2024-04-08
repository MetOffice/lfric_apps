!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>@brief Drives the execution of the (tangent) linear model.
module linear_driver_mod

  use base_mesh_config_mod,       only : prime_mesh_name
  use constants_mod,              only : i_def, r_def, imdi
  use extrusion_mod,              only : TWOD
  use field_array_mod,            only : field_array_type
  use field_mod,                  only : field_type
  use field_collection_mod,       only : field_collection_type
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
  use gungho_modeldb_mod,         only : modeldb_type
  use gungho_step_mod,            only : gungho_step
  use gungho_time_axes_mod,       only : gungho_time_axes_type, &
                                         get_time_axes_from_collection
  use init_fd_prognostics_mod,    only : init_fd_prognostics_dump
  use initial_output_mod,         only : write_initial_output
  use initialization_config_mod,  only : ls_option,                 &
                                         ls_option_file,            &
                                         init_option,               &
                                         init_option_fd_start_dump, &
                                         coarse_aerosol_ancil
  use io_config_mod,              only : write_diag, &
                                         diagnostic_frequency, &
                                         nodal_output_on_w3
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
  use model_clock_mod,            only : model_clock_type
  use multires_coupling_config_mod, &
                                  only : aerosol_mesh_name
  use create_tl_prognostics_mod,  only : create_tl_prognostics

  implicit none

  private
  public initialise, step, finalise

  type( mesh_type ), pointer :: mesh      => null()
  type( mesh_type ), pointer :: twod_mesh => null()

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

    type( mesh_type ),      pointer :: aerosol_mesh      => null()
    type( mesh_type ),      pointer :: aerosol_twod_mesh => null()

    type( field_collection_type ), pointer :: depository => null()

    character(len=*), parameter :: io_context_name = "gungho_atm"

    depository => modeldb%fields%get_field_collection("depository")

    call modeldb%values%add_key_value( 'temperature_correction_rate', &
                                       0.0_r_def )
    call modeldb%values%add_key_value( 'total_dry_mass', 0.0_r_def )
    call modeldb%values%add_key_value( 'total_energy_previous', 0.0_r_def )

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
      call create_tl_prognostics( mesh, twod_mesh,              &
                                  modeldb%model_data%fd_fields, &
                                  depository)
    end if

    ! Instantiate the linearisation state
    call linear_create_ls( modeldb, mesh )

    ! Initialise the fields stored in the model_data
    if ( init_option == init_option_fd_start_dump ) then
      call init_fd_prognostics_dump( modeldb%model_data%fd_fields )
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

    type(field_collection_type), pointer :: moisture_fields => null()
    type(field_array_type), pointer      :: ls_mr_array => null()
    type(field_array_type), pointer      :: ls_moist_dyn_array => null()

    type(gungho_time_axes_type), pointer :: model_axes

    ! Get model_axes out of modeldb
    model_axes => get_time_axes_from_collection(modeldb%values, "model_axes" )

    moisture_fields => modeldb%fields%get_field_collection("moisture_fields")
    call moisture_fields%get_field("ls_mr", ls_mr_array)
    call moisture_fields%get_field("ls_moist_dyn", ls_moist_dyn_array)

    if ( ls_option == ls_option_file ) then
      call update_ls_file_alg( model_axes%ls_times_list, &
                               modeldb%clock,            &
                               modeldb%model_data%ls_fields,     &
                               ls_mr_array%bundle,               &
                               ls_moist_dyn_array%bundle )
    end if

    call linear_step( mesh, twod_mesh, &
                      modeldb, modeldb%clock )

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
    call finalise_model( modeldb,               &
                         program_name )

    call finalise_linear_model( )

    ! Destroy the fields stored in model_data
    call finalise_model_data( modeldb )

    ! Finalise infrastructure and constants
    call finalise_infrastructure(modeldb)

  end subroutine finalise

end module linear_driver_mod
