!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @page linear Linear model Program
!> This is a code that uses the LFRic infrastructure to build a model that
!> is the tangent linear/ perturbation forecast for gungho.

!> @brief Main program used to illustrate linear model functionality.

!> @details This top-level code simply calls initialise, run and finalise
!>          routines that are required to run the model.

program linear_model

  use cli_mod,                only : get_initial_filename
  use driver_collections_mod, only : init_collections, final_collections
  use driver_comm_mod,        only : init_comm, final_comm
  use driver_config_mod,      only : init_config, final_config
  use driver_log_mod,         only : init_logger, final_logger
  use driver_time_mod,        only : init_time, final_time
  use driver_timer_mod,       only : init_timers, final_timers
  use gungho_mod,             only : gungho_required_namelists
  use driver_modeldb_mod,     only : modeldb_type
  use lfric_mpi_mod,          only : global_mpi
  use linear_driver_mod,      only : initialise, step, finalise
  use log_mod,                only : log_event,       &
                                     log_level_trace, &
                                     log_scratch_space

  implicit none

  ! Model run working data set
  type (modeldb_type) :: modeldb

  character(*), parameter   :: application_name = "linear_model"
  character(:), allocatable :: filename

  modeldb%mpi => global_mpi

  call modeldb%configuration%initialise( application_name, table_len=10 )

  call modeldb%values%initialise('values', 5)

  ! Create the depository, prognostics and diagnostics field collections
  call modeldb%fields%add_empty_field_collection("depository", table_len = 100)
  call modeldb%fields%add_empty_field_collection("prognostic_fields",         &
                                                    table_len = 100)
  call modeldb%fields%add_empty_field_collection("diagnostic_fields",         &
                                                    table_len = 100)
  call modeldb%fields%add_empty_field_collection("lbc_fields",                &
                                                    table_len = 100)
  call modeldb%fields%add_empty_field_collection("radiation_fields",          &
                                                    table_len = 100)
  call modeldb%fields%add_empty_field_collection("fd_fields",                 &
                                                     table_len = 100)

  call modeldb%io_contexts%initialise(application_name, 100)

  call init_comm( application_name, modeldb )
  call get_initial_filename( filename )
  call init_config( filename, gungho_required_namelists, &
                    modeldb%configuration )
  call init_logger( modeldb%mpi%get_comm(), application_name )
  call init_timers( application_name )
  call init_collections()
  call init_time( modeldb )
  deallocate( filename )

  call initialise( application_name, modeldb )

  write( log_scratch_space,'("Running ", A, " ...")' ) application_name
  call log_event( log_scratch_space, log_level_trace )
  do while (modeldb%clock%tick())
    call step( modeldb )
  end do

  call log_event( 'Finalising '//application_name//' ...', log_level_trace )
  call finalise( application_name, modeldb )

  call final_time( modeldb )
  call final_collections()
  call final_timers( application_name )
  call final_logger( application_name )
  call final_config()
  call final_comm( modeldb )

end program linear_model
