!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @page Miniapp lfric2lfric program

!> @brief Main program used to illustrate a horizontal regridding application
!!        built using LFRic infrastructure.

!> @details This top-level code simply calls initialise, run and finalise
!!          routines that are required to execute the regridding process.

program lfric2lfric

  use cli_mod,                only: parse_command_line
  use constants_mod,          only: precision_real
  use driver_collections_mod, only: init_collections, final_collections
  use driver_config_mod,      only: init_config, final_config
  use driver_comm_mod,        only: init_comm, final_comm
  use driver_log_mod,         only: init_logger, final_logger
  use driver_modeldb_mod,     only: modeldb_type
  use driver_time_mod,        only: init_time, final_time
#ifdef MCT
  use coupling_mod,           only: coupling_type
#endif
  use log_mod,                only: log_event,       &
                                    log_level_trace, &
                                    log_scratch_space
  use lfric_mpi_mod,          only: global_mpi

  use lfric2lfric_mod,        only: lfric2lfric_required_namelists
  use lfric2lfric_driver_mod, only: initialise, run, finalise

  implicit none

  ! The technical and scientific state
  type(modeldb_type) :: modeldb

  character(len=*), parameter   :: program_name = "lfric2lfric"
  character(:),     allocatable :: filename

#ifdef MCT
  ! Coupler objects
  type(coupling_type)          :: coupler
#endif

  call parse_command_line( filename )

  call modeldb%configuration%initialise( program_name, table_len=10 )

  write(log_scratch_space,'(A)')                          &
      'Application built with '// trim(precision_real) // &
      '-bit real numbers.'
  call log_event( log_scratch_space, log_level_trace )

  modeldb%mpi => global_mpi

  call modeldb%values%initialise( 'values', table_len=5 )

#ifdef MCT
  call modeldb%values%add_key_value('cpl_name', program_name)
  ! Create a second coupling object for single executable coupling
  ! We do not initialise or finalise it as initialising will set up
  ! MPI communication and this can only be done once
  call modeldb%values%add_key_value('coupling_dst', coupler)
#endif
  call init_comm( program_name, modeldb )
  call init_config( filename, lfric2lfric_required_namelists, &
                    modeldb%configuration                     )
  call init_logger( modeldb%mpi%get_comm(), program_name )
  call init_collections()
  call init_time( modeldb )
  deallocate( filename )

  ! Create the field collections and place in modeldb
  call modeldb%fields%add_empty_field_collection("depository")

  call modeldb%io_contexts%initialise(program_name, 100)

  call log_event( 'Initialising ' // program_name // ' ...', log_level_trace )
  call initialise( modeldb )

  call run( modeldb )

  call log_event( 'Finalising ' // program_name // ' ...', log_level_trace )
  call finalise( program_name, modeldb )

  call final_time( modeldb )
  call final_collections()
  call final_logger( program_name )
  call final_config()
  call final_comm( modeldb )

end program lfric2lfric
