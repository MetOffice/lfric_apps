!-----------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief   The top level program for the tangent linear tests.
!>@details The default is to run all available tests - which
!!         test whether the linear code is tangent linear to the
!!         corresponding nonlinear code.
program nwp_gal9

  use driver_collections_mod,  only: init_collections, final_collections
  use driver_time_mod,         only: init_time, final_time
  use driver_comm_mod,         only: init_comm, final_comm
  use driver_log_mod,          only: init_logger, final_logger
  use driver_config_mod,       only: init_config, final_config
  use driver_modeldb_mod,      only: modeldb_type
  use lfric_mpi_mod,           only: global_mpi
  use gungho_mod,              only: gungho_required_namelists
  use log_mod,                 only: log_event,       &
                                     LOG_LEVEL_ERROR, &
                                     LOG_LEVEL_INFO
  use linear_driver_mod,       only: initialise, finalise
  use tl_test_driver_mod,      only: run_timesteps,               &
                                     run_transport_control,       &
                                     run_semi_imp_alg,            &
                                     run_rhs_alg

  implicit none

  ! Model run working data set
  type(modeldb_type) :: modeldb
  character(*), parameter :: application_name = 'nwp_gal9'
  character(:), allocatable :: filename

  ! Variables used for parsing command line arguments
  integer :: length, status, nargs
  character(len=0) :: dummy
  character(len=:), allocatable :: program_name, test_flag

  ! Flags which determine the tests that will be carried out
  logical :: do_test_timesteps = .false.
  logical :: do_test_transport_control = .false.
  logical :: do_test_semi_imp_alg = .false.
  logical :: do_test_rhs_alg = .false.

  ! Usage message to print
  character(len=256) :: usage_message

  modeldb%mpi => global_mpi

  call modeldb%configuration%initialise( application_name, table_len=10 )
  call modeldb%values%initialise('values', 5)

  call log_event( 'TL testing running ...', LOG_LEVEL_INFO )

  ! Create the depository, prognostics and diagnostics field collections
  call modeldb%fields%add_empty_field_collection("depository", table_len = 100)
  call modeldb%fields%add_empty_field_collection("prognostic_fields", &
                                                    table_len = 100)
  call modeldb%fields%add_empty_field_collection("diagnostic_fields", &
                                                    table_len = 100)
  call modeldb%fields%add_empty_field_collection("lbc_fields",        &
                                                    table_len = 100)
  call modeldb%fields%add_empty_field_collection("radiation_fields",  &
                                                    table_len = 100)
  call modeldb%fields%add_empty_field_collection("fd_fields",         &
                                                    table_len = 100)


  call modeldb%io_contexts%initialise(application_name, 100)

  ! Parse command line parameters
  call get_command_argument( 0, dummy, length, status )
  allocate(character(length)::program_name)
  call get_command_argument( 0, program_name, length, status )
  nargs = command_argument_count()

  ! Print out usage message if wrong number of arguments is specified
  if (nargs /= 2) then
     write(usage_message,*) "Usage: ",trim(program_name), &
          " <namelist filename> "      // &
          " test_XXX with XXX in { "   // &
          " timesteps, "               // &
          " transport_control, "       // &
          " semi_imp_alg, "            // &
          " rhs_alg, "                 // &
          " } "
     call log_event( trim(usage_message), LOG_LEVEL_ERROR )
  end if

  call get_command_argument( 1, dummy, length, status )
  allocate( character(length) :: filename )
  call get_command_argument( 1, filename, length, status )

  call get_command_argument( 2, dummy, length, status )
  allocate(character(length)::test_flag)
  call get_command_argument( 2, test_flag, length, status )

  ! Choose test case depending on flag provided in the first command
  ! line argument
  select case (trim(test_flag))
  case ("test_timesteps")
     do_test_timesteps = .true.
  case ("test_transport_control")
     do_test_transport_control = .true.
  case ("test_semi_imp_alg")
     do_test_semi_imp_alg = .true.
  case ("test_rhs_alg")
     do_test_rhs_alg = .true.
  case default
     call log_event( "Unknown test", LOG_LEVEL_ERROR )
  end select

  call init_comm( application_name, modeldb )
  call init_config( filename, gungho_required_namelists, &
       modeldb%configuration )
  call init_logger( modeldb%mpi%get_comm(), application_name )
  call init_collections()
  call init_time( modeldb )
  call initialise( application_name, modeldb )

  if (do_test_timesteps) then
    call run_timesteps(modeldb)
  endif
  if (do_test_transport_control) then
    call run_transport_control(modeldb)
  endif
  if (do_test_rhs_alg) then
    call run_rhs_alg(modeldb)
  endif
  if (do_test_semi_imp_alg) then
    call run_semi_imp_alg(modeldb)
  endif

  call finalise( application_name, modeldb )
  call final_time( modeldb )
  call final_collections()
  call final_logger( application_name )
  call final_config()
  call final_comm( modeldb )

end program nwp_gal9
