!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!> @page gung_ho GungHo Program
!> This is a code that uses the LFRic infrastructure to build a model that
!> just includes the GungHo dynamical core.

!> @brief Main program used to illustrate gungho functionality.

!> @details This top-level code simply calls initialise, run and finalise
!>          routines that are required to run the model.

program gungho_model

  use cli_mod,                only: parse_command_line
  use derived_config_mod,     only: l_esm_couple
  use driver_collections_mod, only: init_collections, final_collections
  use driver_comm_mod,        only: init_comm, final_comm
  use driver_config_mod,      only: init_config, final_config
  use driver_counter_mod,     only: init_counters, final_counters
  use driver_log_mod,         only: init_logger, final_logger
  use driver_time_mod,        only: init_time, final_time
  use gungho_mod,             only: gungho_required_namelists
  use gungho_driver_mod,      only: initialise, step, finalise
  use driver_modeldb_mod,     only: modeldb_type
  use lfric_mpi_mod,          only: global_mpi
  use log_mod,                only: log_event,       &
                                    log_level_info,  &
                                    log_level_trace, &
                                    log_scratch_space
  use namelist_mod,           only: namelist_type
  use timing_mod,             only: init_timing, final_timing
  use io_config_mod,          only: timer_output_path


  implicit none

  ! Model run working data set
  type(modeldb_type) :: modeldb

  character(*), parameter      :: application_name = "gungho_model"
  character(:), allocatable    :: filename
  type(namelist_type), pointer :: io_nml
  logical                      :: lsubroutine_timers

  call parse_command_line( filename )

  modeldb%mpi => global_mpi

  call modeldb%configuration%initialise( application_name, &
                                         table_len=10 )
  call modeldb%values%initialise( 'values', 5 )

  ! Create the depository, prognostics and diagnostics field collections
  call modeldb%fields%add_empty_field_collection("depository"       , &
                                                 table_len = 100)
  call modeldb%fields%add_empty_field_collection("prognostic_fields", &
                                                 table_len = 100)
  call modeldb%fields%add_empty_field_collection("diagnostic_fields", &
                                                 table_len = 100)
  call modeldb%fields%add_empty_field_collection("lbc_fields",        &
                                                  table_len = 100)
  call modeldb%fields%add_empty_field_collection("radiation_fields",  &
                                                  table_len = 100)
  call modeldb%fields%add_empty_field_collection("ls_fields",         &
                                                  table_len = 100)
  call modeldb%fields%add_empty_field_collection("fd_fields",         &
                                                  table_len = 100)

  call modeldb%io_contexts%initialise(application_name, 100)

  call init_comm( application_name, modeldb )

  call init_config( filename, gungho_required_namelists, &
                    modeldb%configuration )
  call init_logger( modeldb%mpi%get_comm(), application_name )
  io_nml => modeldb%configuration%get_namelist('io')
  call io_nml%get_value('subroutine_timers', lsubroutine_timers)
  call init_timing( modeldb%mpi%get_comm(), lsubroutine_timers, application_name, timer_output_path )
  nullify( io_nml )
  call init_collections()
  call init_time( modeldb )
  call init_counters( application_name )
  deallocate( filename )

  write( log_scratch_space, '("Initialise ", A, " ...")' ) application_name
  call log_event( log_scratch_space, log_level_trace )
  call initialise( application_name, modeldb )

  if (l_esm_couple) then
    write(log_scratch_space,'(A)') "Configuration is coupled to one or more components"
  else
    write(log_scratch_space,'("No external coupler is used by this configuration")')
  end if
  call log_event( log_scratch_space, log_level_info )

  write(log_scratch_space,'("Running ", A, " ...")') application_name
  call log_event( log_scratch_space, log_level_trace )
  do while (modeldb%clock%tick())
    call step( modeldb )
  end do

  call log_event( 'Finalising '//application_name//' ...', log_level_trace )
  call finalise( application_name, modeldb )

  call final_counters( application_name )
  call final_time( modeldb )
  call final_collections()
  call final_timing( application_name )
  call final_logger( application_name )
  call final_config()
  call final_comm( modeldb )

end program gungho_model
