!-----------------------------------------------------------------------------
! (C) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> This is a code that uses the LFRic infrastructure to build a model that
!> includes the GungHo dynamical core and physics parametrisation schemes
!> that are currently provided through the use of unified model code.

!> @brief Main program used to illustrate an atmospheric model built using
!>        LFRic infrastructure

!> @details This top-level code simply calls initialise, run and finalise
!>          routines that are required to run the atmospheric model.

program lfric_atm

  use cli_mod,                only: parse_command_line
  use driver_collections_mod, only: init_collections, final_collections
  use driver_comm_mod,        only: init_comm, final_comm
  use driver_config_mod,      only: init_config, final_config
  use driver_counter_mod,     only: init_counters, final_counters
  use driver_log_mod,         only: init_logger, final_logger
  use driver_time_mod,        only: init_time, final_time
  use driver_timer_mod,       only: init_timers, final_timers
  use gungho_mod,             only: gungho_required_namelists
  use driver_modeldb_mod,     only: modeldb_type
  use gungho_driver_mod,      only: initialise, step, finalise
  use lfric_mpi_mod,          only: global_mpi
  use namelist_mod,           only: namelist_type

  use timing_mod,             only: init_timing, start_timing, stop_timing, final_timing, tik, LPROF

  implicit none

  ! Model run working data set
  type(modeldb_type) :: modeldb

  character(*), parameter :: application_name = "lfric_atm"
  character(:), allocatable :: filename
  integer(tik)              :: timing_handle_global
  type(namelist_type), pointer :: io_nml
  logical :: lsubroutine_timers

  call parse_command_line( filename )

  modeldb%mpi => global_mpi

  call modeldb%configuration%initialise( application_name, &
                                         table_len=10 )
  call modeldb%values%initialise( 'values', 5 )

  ! Create the depository, prognostics and diagnostics field collections
  call modeldb%fields%add_empty_field_collection("depository", &
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
  call init_timers( application_name )

  io_nml => modeldb%configuration%get_namelist('io')
  call io_nml%get_value('subroutine_timers', lsubroutine_timers)
  call init_timing( modeldb%mpi%get_comm(), lsubroutine_timers )
  nullify( io_nml )
  if ( LPROF ) call start_timing( timing_handle_global, '__lfric_atm__ ')

  call init_collections()
  call init_time( modeldb )
  call init_counters( application_name )
  deallocate( filename )

  call initialise( application_name, modeldb )
  do while (modeldb%clock%tick())
    call step( modeldb )
  end do
  call finalise( application_name, modeldb )

  call final_counters( application_name )
  call final_time( modeldb )
  call final_collections()

  if ( LPROF ) call stop_timing( timing_handle_global )
  call final_timing()

  call final_timers( application_name )
  call final_logger( application_name )
  call final_config()
  call final_comm( modeldb )

end program lfric_atm
