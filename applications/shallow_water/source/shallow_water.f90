!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @page shallow_water Shallow water equations miniapp
!> This is code that uses the LFRic infrastructure to build a shallow water
!> model that includes some of the GungHo routines.
!>
!> @brief Main program used to simulate shallow water equations.
!>
!> @details This top-level code simply calls initialise, run and finalise
!>          routines that are required to run the model.

program shallow_water

  use cli_mod,                   only: parse_command_line
  use driver_collections_mod,    only: init_collections, final_collections
  use driver_comm_mod,           only: init_comm, final_comm
  use driver_config_mod,         only: init_config, final_config
  use driver_counter_mod,        only: init_counters, final_counters
  use driver_log_mod,            only: init_logger, final_logger
  use driver_modeldb_mod,        only: modeldb_type
  use driver_time_mod,           only: init_time, final_time
  use driver_timer_mod,          only: init_timers, final_timers
  use lfric_mpi_mod,             only: global_mpi
  use log_mod,                   only: log_event,       &
                                       log_level_trace, &
                                       log_scratch_space
  use shallow_water_mod,         only: shallow_water_required_namelists
  use shallow_water_driver_mod,  only: initialise, &
                                       step,       &
                                       finalise

  implicit none

  character(*), parameter :: program_name = "shallow_water"

  ! Model run working data set
  type(modeldb_type)        :: modeldb
  character(:), allocatable :: filename

  call parse_command_line( filename )

  modeldb%mpi => global_mpi

  call modeldb%configuration%initialise( program_name, table_len=10 )

  ! Create the depository and prognostics field collections
  call modeldb%fields%add_empty_field_collection("depository", &
                                                 table_len = 100)
  call modeldb%fields%add_empty_field_collection("prognostics", &
                                                 table_len = 100)

  ! Initialise the key_value pair collection in the modeldb
  call modeldb%values%initialise( 'values', 5 )

  call modeldb%io_contexts%initialise(program_name, 100)

  call init_comm( program_name, modeldb )
  call init_config( filename, shallow_water_required_namelists, &
                    modeldb%configuration )
  call init_logger( global_mpi%get_comm(), program_name )
  call init_timers( program_name )
  call init_counters( program_name )
  call init_collections()
  call init_time( modeldb )
  deallocate( filename )

  call log_event( 'Initialising Infrastructure ...', log_level_trace )
  call initialise( modeldb, program_name, modeldb%calendar )
  write(log_scratch_space,'("Running ", A, "...")') program_name
  call log_event( log_scratch_space, log_level_trace )
  do while (modeldb%clock%tick())
    call step( modeldb )
  end do

  call log_event( 'Finalising ' // program_name // ' ...', log_level_trace )
  call finalise( modeldb, program_name )

  call final_time( modeldb )
  call final_collections()
  call final_counters( program_name )
  call final_timers( program_name )
  call final_logger( program_name )
  call final_config()
  call final_comm( modeldb )

end program shallow_water
