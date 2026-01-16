!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @page name_transport NAME Transport miniapp
!> Program file for running name transport miniapp. Subroutine calls include
!! initialise_name_transport(), run_name_transport() and finalise_name_transport().

program name_transport

  use cli_mod,                   only: parse_command_line
  use constants_mod,             only: i_def, r_def
  use driver_collections_mod,    only: init_collections, final_collections
  use driver_comm_mod,           only: init_comm, final_comm
  use driver_config_mod,         only: init_config, final_config
  use driver_log_mod,            only: init_logger, final_logger
  use driver_modeldb_mod,        only: modeldb_type
  use driver_time_mod,           only: init_time, final_time
  use driver_timer_mod,          only: init_timers, final_timers
  use lfric_mpi_mod,             only: global_mpi
  use log_mod,                   only: log_event,       &
                                       log_level_debug, &
                                       log_level_info,  &
                                       log_level_trace, &
                                       log_scratch_space
  use namelist_collection_mod,   only: namelist_collection_type

  use name_transport_mod,        only: name_transport_required_namelists
  use name_transport_driver_mod, only: initialise_name_transport, &
                                       step_name_transport,       &
                                       finalise_name_transport

  implicit none

  type(modeldb_type) :: modeldb
  character(*), parameter   :: program_name = "name_transport"
  character(:), allocatable :: filename

  call parse_command_line( filename )

  call modeldb%configuration%initialise( program_name, table_len=10 )
  modeldb%mpi => global_mpi
  call init_comm( program_name, modeldb )
  call init_config( filename, name_transport_required_namelists, &
                    modeldb%configuration )
  call init_logger( modeldb%mpi%get_comm(), program_name )

  call log_event( 'Miniapp will run with default precision set as:', &
    log_level_info )
  write(log_scratch_space, '("        r_def kind = ", I0)') kind(1.0_r_def)
  call log_event( log_scratch_space , log_level_info )
  write(log_scratch_space, '("        i_def kind = ", I0)') kind(1_i_def)
  call log_event( log_scratch_space , log_level_info )

  call init_timers( program_name )
  call init_collections()
  call init_time( modeldb )
  deallocate( filename )
  call modeldb%io_contexts%initialise(program_name, 100)

  ! Initialise name_transport app
  call log_event( 'Initialising ' // program_name // ' ...', log_level_debug )
  call initialise_name_transport( program_name, modeldb )

  ! Run time steps
  call log_event( 'Running ' // program_name // ' ...', log_level_trace )
  do while (modeldb%clock%tick())
    call step_name_transport( modeldb%clock )
  end do

  ! Finalise name_transport app
  call log_event( 'Finalising ' // program_name // ' ...', log_level_debug )
  call finalise_name_transport( program_name, modeldb )

  call final_time( modeldb )
  call final_collections()
  call final_timers( program_name )
  call final_logger( program_name )
  call final_config()
  call final_comm( modeldb )

end program name_transport
