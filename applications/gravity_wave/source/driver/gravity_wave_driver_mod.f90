!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> Drives the execution of the gravity_wave miniapp.
!>
module gravity_wave_driver_mod

  use base_mesh_config_mod,     only: prime_mesh_name
  use driver_modeldb_mod,       only: modeldb_type
  use constants_mod,            only: i_def, r_def
  use gravity_wave_infrastructure_mod,                           &
                                only: initialise_infrastructure, &
                                       finalise_infrastructure
  use create_gravity_wave_prognostics_mod, &
                                only: create_gravity_wave_prognostics
  use field_mod,                only: field_type
  use function_space_chain_mod, only: function_space_chain_type
  use gravity_wave_constants_config_mod,             &
                                only: b_space,       &
                                      b_space_w0,    &
                                      b_space_w3,    &
                                      b_space_wtheta
  use gravity_wave_diagnostics_driver_mod, &
                                only: gravity_wave_diagnostics_driver
  use gw_init_fields_alg_mod,   only: gw_init_fields_alg
  use gravity_wave_alg_mod,     only: gravity_wave_alg_init, &
                                      gravity_wave_alg_step, &
                                      gravity_wave_alg_final
  use io_config_mod,            only: write_diag,            &
                                      checkpoint_read,       &
                                      checkpoint_write,      &
                                      diagnostic_frequency,  &
                                      use_xios_io,           &
                                      nodal_output_on_w3
  use io_context_mod,           only: io_context_type
  use sci_checksum_alg_mod,     only: checksum_alg

  use boundaries_config_mod,    only: limited_area
  use log_mod,                  only: log_event,        &
                                      log_level_always, &
                                      log_level_info,   &
                                      log_level_trace,  &
                                      log_level_error,  &
                                      log_scratch_space
  use mesh_mod,                 only: mesh_type
  use mesh_collection_mod,      only: mesh_collection
  use namelist_collection_mod,  only: namelist_collection_type
  use io_mod,                   only: ts_fname
  use files_config_mod,         only: checkpoint_stem_name

  implicit none

  private

  public initialise, step, finalise

  ! The prognostic fields
  type( field_type ), target :: wind
  type( field_type ), target :: pressure
  type( field_type ), target :: buoyancy

  type( mesh_type ), pointer :: mesh => null()

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Sets up required state in preparation for run.
  !>
  !> @param[in]    program_name  The name of the program being run
  !> @param[inout] modeldb       The modeldb object
  subroutine initialise( program_name, modeldb )

  implicit none

  character(*),       intent(in)    :: program_name
  type(modeldb_type), intent(inout) :: modeldb

  ! Initialise aspects of the infrastructure
  call initialise_infrastructure( program_name, modeldb)

  ! The limited area version is unable to work with buoyancy in W0 when using
  ! a biperiodic mesh. However, this could be solved by breaking the continuity
  ! at the edges - so that the W0 dofs on the boundary are not shared, e.g. by
  ! using a non-periodic mesh.
  if ( limited_area ) then
    select case( b_space )
      case( b_space_w0 )
        call log_event( 'Limited area version currently unable to run with buoyancy in W0', LOG_LEVEL_ERROR )
      case( b_space_w3 )
        call log_event( 'Limited area version has not been tested with buoyancy in W3', LOG_LEVEL_ERROR )
      case( b_space_wtheta )
        call log_event( 'Running limited area version with buoyancy in Wtheta', LOG_LEVEL_INFO )
      case default
        call log_event( 'Invalid buoyancy space', LOG_LEVEL_ERROR )
    end select
  endif

  ! Create the prognostic fields
  mesh => mesh_collection%get_mesh(prime_mesh_name)
  call create_gravity_wave_prognostics(mesh, wind, pressure, buoyancy)

  ! Initialise prognostic fields
  if (checkpoint_read) then                 ! R ecorded check point to start from
     call log_event( 'Reading checkpoint file to restart wind', LOG_LEVEL_INFO)
     call wind%read_checkpoint("restart_wind",                                &
          trim(ts_fname(checkpoint_stem_name,                                 &
          "","wind",modeldb%clock%get_step()-1,"")) )

     call log_event( 'Reading checkpoint file to restart pressure',           &
          LOG_LEVEL_INFO)
     call pressure%read_checkpoint("restart_pressure",                        &
          trim(ts_fname(checkpoint_stem_name,"",                              &
          "pressure",modeldb%clock%get_step()-1,"")) )

     call log_event( 'Reading checkpoint file to restart buoyancy',           &
          LOG_LEVEL_INFO)
     call buoyancy%read_checkpoint("restart_buoyancy",                        &
          trim(ts_fname(checkpoint_stem_name,"",                              &
          "buoyancy",modeldb%clock%get_step()-1,"")) )

  else                                      ! No check point to start from
     call gw_init_fields_alg(wind, pressure, buoyancy)
  end if

  ! Initialise the gravity-wave model
  call gravity_wave_alg_init(mesh, wind, pressure, buoyancy)

  ! Output initial conditions
  ! We only want these once at the beginning of a run
  if (modeldb%clock%is_initialisation() .and. write_diag) then
    call gravity_wave_diagnostics_driver( mesh,        &
                                          wind,        &
                                          pressure,    &
                                          buoyancy,    &
                                          modeldb%clock, &
                                          nodal_output_on_w3)
  end if

  end subroutine initialise

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Performs time step.
  !> @param[in]    program_name  The name of the program being run
  !> @param[inout] modeldb       The modeldb object
  subroutine step(program_name, modeldb)

    implicit none

    character(*),        intent(in)    :: program_name
    class(modeldb_type), intent(inout) :: modeldb


    ! Update XIOS calendar if we are using it for diagnostic output or
    ! checkpoint
    !
    if ( use_xios_io ) then
      call log_event( program_name//': Updating XIOS timestep', &
                      LOG_LEVEL_INFO )
    end if

    write( log_scratch_space, '("/",A,"\ ")' ) repeat('*', 76)
    call log_event( log_scratch_space, LOG_LEVEL_TRACE )
    write( log_scratch_space, &
           '(A,I0)' ) 'Start of timestep ', modeldb%clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    call gravity_wave_alg_step(wind, pressure, buoyancy, modeldb)

    write( log_scratch_space, &
           '(A,I0)' ) 'End of timestep ', modeldb%clock%get_step()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '("\",A,"/ ")' ) repeat('*', 76)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    if ( (mod(modeldb%clock%get_step(), diagnostic_frequency) == 0) &
         .and. (write_diag) ) then

      call log_event("Gravity Wave: writing diagnostic output", LOG_LEVEL_INFO)

      call gravity_wave_diagnostics_driver( mesh,        &
                                            wind,        &
                                            pressure,    &
                                            buoyancy,    &
                                            modeldb%clock, &
                                            nodal_output_on_w3)
    end if

  end subroutine step


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> Tidies up after a run.
  !> @param[in]    program_name  The name of the program being run
  !> @param[inout] modeldb       The modeldb object
  subroutine finalise(program_name, modeldb)

  implicit none

  character(*),        intent(in)    :: program_name
  class(modeldb_type), intent(inout) :: modeldb

  !--------------------------------------------------------------------------
  ! Model finalise
  !--------------------------------------------------------------------------

  ! Write checksums to file
  call checksum_alg( program_name, wind, 'wind', buoyancy, 'buoyancy', pressure, 'pressure')
  call gravity_wave_alg_final()

  ! Write checkpoint/restart files if required
  if( checkpoint_write ) then
     call wind%write_checkpoint("restart_wind",                               &
          trim(ts_fname(checkpoint_stem_name,                                 &
          "","wind",modeldb%clock%get_step(),"")) )

     call pressure%write_checkpoint("restart_pressure",                       &
          trim(ts_fname(checkpoint_stem_name,"",                              &
          "pressure",modeldb%clock%get_step(),"")) )

     call buoyancy%write_checkpoint("restart_buoyancy",                       &
          trim(ts_fname(checkpoint_stem_name,"",                              &
          "buoyancy",modeldb%clock%get_step(),"")) )

  end if

  !--------------------------------------------------------------------------
  ! Driver layer finalise
  !--------------------------------------------------------------------------

  call log_event( program_name//' completed.', LOG_LEVEL_ALWAYS )

  call finalise_infrastructure(modeldb)

  end subroutine finalise

end module gravity_wave_driver_mod
