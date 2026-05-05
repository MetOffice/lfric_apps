!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @page gravity_wave Gravity Wave miniapp
!> Test program for the automatic generation of boundary condition enforcement
!> by PSyclone.
!>
!> @brief Main program used to simulate the linear gravity waves equations.

program gravity_wave

  use cli_mod,                 only: parse_command_line
  use driver_modeldb_mod,      only: modeldb_type
  use driver_collections_mod,  only: init_collections, final_collections
  use driver_comm_mod,         only: init_comm, final_comm
  use driver_config_mod,       only: init_config, final_config
  use driver_log_mod,          only: init_logger, final_logger
  use driver_time_mod,         only: init_time, final_time
  use driver_timer_mod,        only: init_timers, final_timers
  use gravity_wave_mod,        only: gravity_wave_required_namelists
  use gravity_wave_driver_mod, only: initialise, step, finalise
  use lfric_mpi_mod,           only: global_mpi
  use log_mod,                 only: log_event,       &
                                     log_level_trace, &
                                     log_scratch_space

  implicit none

  type(modeldb_type)        :: modeldb
  character(*), parameter   :: program_name = "gravity_wave"
  character(:), allocatable :: filename

  call parse_command_line( filename )

  call modeldb%configuration%initialise( program_name, table_len=10 )

  modeldb%mpi => global_mpi
  call init_comm( program_name, modeldb )
  call init_config( filename, gravity_wave_required_namelists, &
                    modeldb%configuration )
  deallocate( filename )

  call init_logger( modeldb%mpi%get_comm(), program_name )
  call init_collections()
  call init_timers( program_name )
  call init_time( modeldb )

  ! Create the depository field collection and place it in modeldb
  call modeldb%fields%add_empty_field_collection("depository")

  call modeldb%io_contexts%initialise(program_name, 100)
  call log_event( 'Initialising '//program_name//' ...', log_level_trace )
  call initialise( program_name, modeldb )

  write(log_scratch_space,'("Running ",A, " ...")') program_name
  call log_event( log_scratch_space, log_level_trace )
  do while (modeldb%clock%tick())
    call step( program_name, modeldb )
  end do

  call log_event( 'Finalising '//program_name//' ...', log_level_trace )
  call finalise( program_name, modeldb )

  call final_time( modeldb )
  call final_timers( program_name )
  call final_collections()
  call final_logger( program_name )
  call final_config()
  call final_comm( modeldb )

end program gravity_wave
