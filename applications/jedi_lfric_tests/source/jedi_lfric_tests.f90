!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page Miniapp jedi_lfric program

!> @brief Tests running the linear model in the miniapp

!> @details Calls init, step and finalise routines from a
!>          DA specific linear driver module in the
!>          jedi-lfric component

program jedi_lfric_tests

  use cli_mod,                only: get_initial_filename
  use gungho_mod,             only: gungho_required_namelists
  use linear_driver_mod,      only: initialise, step, finalise
  use driver_collections_mod, only: init_collections, final_collections
  use driver_comm_mod,        only: init_comm, final_comm
  use driver_config_mod,      only: init_config, final_config
  use driver_log_mod,         only: init_logger, final_logger
  use gungho_modeldb_mod,     only: modeldb_type
  use driver_time_mod,        only: init_time, final_time
  use constants_mod,          only: i_def
  use log_mod,                only: log_event, log_level_trace
  use mpi_mod,                only: global_mpi

  implicit none

  type(modeldb_type) :: modeldb

  character(*), parameter :: program_name = "jedi_lfric_tests"

  character(:), allocatable :: filename

  modeldb%mpi => global_mpi

  call modeldb%configuration%initialise( program_name, table_len=10 )

  call init_comm( program_name, modeldb%mpi )
  call get_initial_filename( filename )
  call init_config( filename, gungho_required_namelists, &
                    modeldb%configuration )
  call init_logger( global_mpi%get_comm(), program_name )
  call init_collections()
  deallocate( filename )

  call modeldb%values%initialise( 'values', 5 )

  ! Create the depository, prognostics and diagnostics field collections
  call modeldb%fields%add_empty_field_collection("depository", table_len = 100)
  call modeldb%fields%add_empty_field_collection("prognostic_fields", &
                                                  table_len = 100)
  call modeldb%fields%add_empty_field_collection("diagnostic_fields", &
                                                  table_len = 100)

  call modeldb%io_contexts%initialise(program_name, 100)

  call log_event( 'Initialising ' // program_name // ' ...', log_level_trace )
  call init_time( modeldb%clock, modeldb%calendar )
  call initialise( program_name, modeldb )

  call log_event( 'Running ' // program_name // ' ...', log_level_trace )
  do while ( modeldb%clock%tick() )
    call step( modeldb )
  end do

  call log_event( 'Finalising ' // program_name // ' ...', log_level_trace )
  call finalise( program_name, modeldb )

  call final_time( modeldb%clock, modeldb%calendar )
  call final_collections()
  call final_logger( program_name )
  call final_config()
  call final_comm( modeldb%mpi )

end program jedi_lfric_tests
