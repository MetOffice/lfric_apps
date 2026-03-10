!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief    Module containing routines needed to initialise and finalise the
!>           modeldb object and also the linear model step method used by the
!>           tangent linear model from LFRic-JEDI
!>
!> @detail   The initialise and finalise methods set up/close down the modeldb
!>           as used by the JEDI interface to the linear model. This code uses
!>           the linear model as the basis to define the routines that need to
!>           be called to achieve this. The required code to setup the modeldb
!>           is contained in the linear program file and driver located in:
!>
!>             applications/linear_model/source/linear_model.f90
!>             science/linear/source/driver/linear_driver_mod.f90
!>
!>           A step method is also included that calls the linear model step
!>           method and ticks the model clock.
!>
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
! Notes for maintainer: This file is an attempt to follow the modelDB design
! outlined in:
!  https://code.metoffice.gov.uk/trac/lfric/wiki/ticket/3546/ModelDBnotes
! where each model has a method to initialise and finalise the modeldb. At
! present this design is not yet fully implemented in LFRic. For the linear
! model, the code that constitutes the full Initialisation and finalisation of
! modeldb is included in the program file and driver layer module. The code
! here is based on those files and includes the same sequence of calls that are
! related to initialising the modeldb. Due to differences in how JEDI handles
! the linearization state, the code here has modifications to allow for those
! differences. This module will therefore need to be updated in conjunction
! with the linear model to ensure the JEDI interface to the linear model does
! not break.
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------
module jedi_lfric_linear_modeldb_driver_mod

  use adjoint_model_mod,            only : initialise_adjoint_model, &
                                           adjoint_step,             &
                                           finalise_adjoint_model
  use atl_si_timestep_alg_mod,      only : atl_si_timestep_type
  use constants_mod,                only : r_def, l_def, str_def
  use driver_config_mod,            only : init_config
  use driver_time_mod,              only : init_time, final_time
  use extrusion_mod,                only : TWOD
  use gungho_init_fields_mod,       only : create_model_data, &
                                           finalise_model_data
  use gungho_mod,                   only : gungho_required_namelists
  use driver_modeldb_mod,           only : modeldb_type
  use gungho_model_mod,             only : initialise_infrastructure, &
                                           initialise_model, &
                                           finalise_infrastructure, &
                                           finalise_model
  use gungho_time_axes_mod,         only : gungho_time_axes_type
  use lfric_mpi_mod,                only : lfric_mpi_type
  use lfric_xios_context_mod,       only : lfric_xios_context_type
  use log_mod,                      only : log_event,         &
                                           log_scratch_space, &
                                           LOG_LEVEL_TRACE,   &
                                           LOG_LEVEL_ERROR
  use linear_model_data_mod,        only : linear_create_ls_analytic
  use linear_model_mod,             only : initialise_linear_model, &
                                           finalise_linear_model
  use linear_step_mod,              only : linear_step
  use mesh_mod,                     only : mesh_type
  use mesh_collection_mod,          only : mesh_collection
  use namelist_mod,                 only : namelist_type

  implicit none

  private
  public initialise_modeldb, finalise_modeldb
  public step_tl, step_ad, identity_step_tl

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up the modeldb and required infrastructure.
  !>
  !> @param [in]     modeldb_name    An identifier given to the modeldb
  !> @param [in]     filename        The configuration file path and name
  !> @param [in]     mpi_obj         The mpi object to be associated with the modeldb
  !> @param [in,out] modeldb         The structure that holds model state
  !> @param [in]     atl_si_timestep Object encapsulating adjoint model timestep routines
  subroutine initialise_modeldb( modeldb_name, filename, mpi_obj, modeldb, atl_si_timestep )

    implicit none

    character(*),                 intent(in)    :: modeldb_name
    character(len=*),             intent(in)    :: filename
    type(lfric_mpi_type), target, intent(in)    :: mpi_obj
    type(modeldb_type),           intent(inout) :: modeldb
    type(atl_si_timestep_type),   intent(inout) :: atl_si_timestep

    ! Local
    type( gungho_time_axes_type )            :: model_axes
    type( mesh_type ),               pointer :: mesh
    type( mesh_type ),               pointer :: twod_mesh
    type( mesh_type ),               pointer :: aerosol_mesh
    type( mesh_type ),               pointer :: aerosol_twod_mesh
    type( namelist_type ),           pointer :: base_mesh_nml
    type( namelist_type ),           pointer :: multires_coupling_nml
    type( namelist_type ),           pointer :: initialization_nml
    type( lfric_xios_context_type ), pointer :: io_context
    character( len=str_def )                 :: prime_mesh_name
    character( len=str_def )                 :: aerosol_mesh_name
    logical( kind=l_def )                    :: coarse_aerosol_ancil
    logical( kind=l_def )                    :: coarse_ozone_ancil

    character(len=*), parameter :: io_context_name = "gungho_atm"

    nullify( mesh, twod_mesh, aerosol_mesh, aerosol_twod_mesh )
    nullify( base_mesh_nml, multires_coupling_nml, initialization_nml )

    ! 1. Initialise modeldb field collections, configuration and mpi.
    modeldb%mpi => mpi_obj
    call modeldb%configuration%initialise( modeldb_name, table_len=10 )

    call modeldb%values%initialise('values', 5)

    ! Create the depository, prognostics and diagnostics field collections
    call modeldb%fields%add_empty_field_collection( "depository", &
                                                    table_len = 100 )
    call modeldb%fields%add_empty_field_collection( "prognostic_fields", &
                                                    table_len = 100)
    call modeldb%fields%add_empty_field_collection( "diagnostic_fields", &
                                                    table_len = 100 )
    call modeldb%fields%add_empty_field_collection( "lbc_fields", &
                                                    table_len = 100 )
    call modeldb%fields%add_empty_field_collection( "radiation_fields", &
                                                    table_len = 100 )
    call modeldb%fields%add_empty_field_collection( "fd_fields",        &
                                                    table_len = 100 )

    call modeldb%io_contexts%initialise(modeldb_name, table_len=100)

    call init_config( filename, gungho_required_namelists, &
                      modeldb%configuration )

    ! 2. Setup some model modeldb%values and initialise infrastructure
    call modeldb%values%add_key_value( 'temperature_correction_rate', 0.0_r_def )
    call modeldb%values%add_key_value( 'total_dry_mass', 0.0_r_def )
    call modeldb%values%add_key_value( 'total_energy_previous', 0.0_r_def )

    ! Initialise infrastructure
    call init_time( modeldb )
    call initialise_infrastructure( io_context_name, modeldb )

    ! Add a place to store time axes in modeldb
    call modeldb%values%add_key_value('model_axes', model_axes)

    ! 3. Setup the required mesh's and create the modeldb fields

    ! Get primary and 2D meshes for initialising model data
    base_mesh_nml => modeldb%configuration%get_namelist('base_mesh')
    call base_mesh_nml%get_value( 'prime_mesh_name', prime_mesh_name )
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)

    ! Get aerosol ancillary configuration logical
    initialization_nml => modeldb%configuration%get_namelist('initialization')
    call initialization_nml%get_value( 'coarse_aerosol_ancil', &
                                       coarse_aerosol_ancil )
    call initialization_nml%get_value( 'coarse_ozone_ancil', &
                                       coarse_ozone_ancil )
    if (coarse_aerosol_ancil .or. coarse_ozone_ancil) then
      ! For now use the coarsest mesh
      multires_coupling_nml => &
        modeldb%configuration%get_namelist('multires_coupling')
      call multires_coupling_nml%get_value( 'aerosol_mesh_name', &
                                            aerosol_mesh_name )
      aerosol_mesh => mesh_collection%get_mesh(aerosol_mesh_name)
      aerosol_twod_mesh => mesh_collection%get_mesh(aerosol_mesh, TWOD)
      write( log_scratch_space, '(A,A)' ) "aerosol mesh name:", &
                                         aerosol_mesh%get_mesh_name()
      call log_event( log_scratch_space, LOG_LEVEL_TRACE )
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

    ! Instantiate the linearisation state
    call linear_create_ls_analytic( modeldb, mesh, twod_mesh )

    ! 4. Initialise the model

    ! Model configuration initialisation
    call initialise_model( mesh, &
                           modeldb )

    ! Linear model configuration initialisation
    call initialise_linear_model( mesh, &
                                  modeldb )

    call initialise_adjoint_model( modeldb, &
                                   atl_si_timestep )

    ! Close IO context and clock so that it is not linked to XIOS
    call modeldb%io_contexts%get_io_context(io_context_name, io_context)
    call io_context%finalise_xios_context()   
    call final_time( modeldb )

    call log_event( "end of initialise_modeldb: initialise_linear_model", &
                    LOG_LEVEL_TRACE )

  end subroutine initialise_modeldb

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Increments the model state by a single timestep.
  !>
  !> @param [in,out] modeldb   The structure that holds model state
  subroutine step_tl( modeldb )

    implicit none

    type( modeldb_type ), intent(inout) :: modeldb

    ! Local
    logical( kind=l_def )          :: clock_running
    type( mesh_type ),     pointer :: mesh
    type( mesh_type ),     pointer :: twod_mesh
    type( namelist_type ), pointer :: base_mesh_nml
    character( len=str_def )       :: prime_mesh_name

    nullify(mesh, twod_mesh, base_mesh_nml)

    ! 1. Tick the clock and check its still running
    clock_running = modeldb%clock%tick()
    ! If the clock has finished then we need to abort and that may be due to
    ! incorrect model_clock configuration
    if ( .not. clock_running ) then
      write ( log_scratch_space, '(A)' ) &
              "jedi_lfric_linear_modeldb::The LFRic model clock has stopped."
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    ! Get mesh
    base_mesh_nml => modeldb%configuration%get_namelist('base_mesh')
    call base_mesh_nml%get_value( 'prime_mesh_name', prime_mesh_name )
    mesh => mesh_collection%get_mesh(prime_mesh_name)
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)

    ! 2. Step the linear model
    call linear_step( mesh,      &
                      twod_mesh, &
                      modeldb,   &
                      modeldb%clock )

  end subroutine step_tl

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Increments the model state by a single timestep (adjoint).
  !>
  !> @param [in,out] modeldb         The structure that holds model state
  !> @param [in]     atl_si_timestep Object encapsulating adjoint model timestep routines
  subroutine step_ad( modeldb, atl_si_timestep )

    implicit none

    type( modeldb_type ),         intent(inout) :: modeldb
    type( atl_si_timestep_type ), intent(inout) :: atl_si_timestep

    ! Local
    logical( kind=l_def ) :: clock_running

    ! 1. Tick the clock and check its still running
    clock_running = modeldb%clock%tick()

    ! If the clock has finished then we need to abort and that may be due to
    ! incorrect model_clock configuration
    if ( .not. clock_running ) then
      write ( log_scratch_space, '(A)' ) &
              "jedi_lfric_linear_modeldb::The LFRic model clock has stopped."
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    ! 2. Step the adjoint model
    call adjoint_step( modeldb, atl_si_timestep )

  end subroutine step_ad

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Increments the model state by a single timestep.
  !>
  !> @param [in,out] modeldb   The structure that holds model state
  subroutine identity_step_tl( modeldb )

    implicit none

    type( modeldb_type ), intent(inout) :: modeldb

    ! Local
    logical( kind=l_def ) :: clock_running

    ! Tick the clock and check its still running
    clock_running = modeldb%clock%tick()
    ! If the clock has finished then we need to abort and that may be due to
    ! incorrect model_clock configuration
    if ( .not. clock_running ) then
      write ( log_scratch_space, '(A)' ) &
              "jedi_lfric_linear_modeldb::The LFRic model clock has stopped."
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

  end subroutine identity_step_tl

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Finalise the modeldb
  !> @param [in,out] modeldb         The structure that holds model state
  !> @param [in]     atl_si_timestep Object encapsulating adjoint model timestep routines
  subroutine finalise_modeldb( modeldb, atl_si_timestep )

    implicit none

    type(modeldb_type),         intent(inout) :: modeldb
    type(atl_si_timestep_type), intent(inout) :: atl_si_timestep

    call log_event( 'Finalising linear model modeldb', LOG_LEVEL_TRACE )

    ! Model configuration finalisation
    call finalise_model( modeldb )

    ! Linear model configuration finalisation
    call finalise_linear_model()
    call finalise_adjoint_model(atl_si_timestep)

    ! Destroy the fields stored in the modeldb model_data
    call finalise_model_data( modeldb )

    ! Finalise infrastructure and constants
    call finalise_infrastructure( modeldb )

  end subroutine finalise_modeldb

end module jedi_lfric_linear_modeldb_driver_mod
