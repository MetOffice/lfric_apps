!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief   The top level program for the da dev, jedi interface
!!         integrationf tests.
!>@details Sets up and runs the integration tests specified in
!!         jedi_lfric_time_test.py.
program jedi_lfric_time_test

  use configuration_mod,               only : final_configuration, &
                                              read_configuration
  use constants_mod,                   only : i_def, r_def, l_def
  use halo_comms_mod,                  only : initialise_halo_comms, &
                                              finalise_halo_comms
  use test_jedi_lfric_time_driver_mod, only : test_jedi_interface_init,             &
                                              test_jedi_interface_final,            &
                                              run_init_string_err,                  &
                                              run_copy_from_jedi_datetime_err,      &
                                              run_add_duration_to_datetime,         &
                                              run_duration_from_datetimes,          &
                                              run_YYYYMMDD_to_JDN,                  &
                                              run_JDN_to_YYYYMMDD_invalid,          &
                                              run_hhmmss_to_seconds,                &
                                              run_seconds_to_hhmmss_large,          &
                                              run_seconds_to_hhmmss_neg,            &
                                              run_duration_init_bad_string_err,     &
                                              run_duration_divide_zero_err,         &
                                              run_duration_divide_remainder_err,    &
                                              run_duration_divide_int_zero_err,     &
                                              run_duration_divide_int_remainder_err
  use lfric_mpi_mod,                   only : global_mpi, &
                                              create_comm, destroy_comm, &
                                              lfric_comm_type
  use log_mod,                         only : log_event,          &
                                              initialise_logging, &
                                              finalise_logging,   &
                                              LOG_LEVEL_ERROR,    &
                                              LOG_LEVEL_INFO
  use namelist_collection_mod,         only : namelist_collection_type

  implicit none

  ! MPI communicator
  type(lfric_comm_type) :: comm

  ! Number of processes and local rank
  integer(i_def) :: total_ranks, local_rank

  character(:), allocatable :: filename

  type(namelist_collection_type), save :: configuration

  ! Variables used for parsing command line arguments
  integer(i_def)   :: length, status, nargs
  character(len=0) :: dummy
  character(len=:), allocatable :: program_name, test_flag
  character(len=:), allocatable :: optional_arg

  ! Flags which determine the tests that will be carried out
  logical(l_def) :: do_test_init_string_err = .false.
  logical(l_def) :: do_test_copy_from_jedi_datetime_err = .false.
  logical(l_def) :: do_test_add_duration_to_datetime = .false.
  logical(l_def) :: do_test_duration_from_datetimes = .false.
  logical(l_def) :: do_test_YYYYMMDD_to_JDN = .false.
  logical(l_def) :: do_test_JDN_to_YYYYMMDD_invalid = .false.
  logical(l_def) :: do_test_hhmmss_to_seconds = .false.
  logical(l_def) :: do_test_seconds_to_hhmmss_large = .false.
  logical(l_def) :: do_test_seconds_to_hhmmss_neg = .false.
  logical(l_def) :: do_test_duration_init_bad_string_err = .false.
  logical(l_def) :: do_test_duration_divide_zero_err = .false.
  logical(l_def) :: do_test_duration_divide_remainder_err = .false.
  logical(l_def) :: do_test_duration_divide_int_zero_err = .false.
  logical(l_def) :: do_test_duration_divide_int_remainder_err = .false.

  ! Usage message to print
  character(len=512) :: usage_message

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Communicators and Logging Setup
  ! Initialise MPI communicatios and get a valid communicator
  call create_comm(comm)

  ! Save the communicator for later use
  call global_mpi%initialise(comm)

  ! Initialise halo functionality
  call initialise_halo_comms(comm)

  total_ranks = global_mpi%get_comm_size()
  local_rank  = global_mpi%get_comm_rank()

  call initialise_logging( comm%get_comm_mpi_val(), 'jedi-interface_test' )

  call log_event( 'jedi interface testing running ...', LOG_LEVEL_INFO )

  ! Parse command line parameters
  call get_command_argument( 0, dummy, length, status )
  allocate(character(length)::program_name)
  call get_command_argument( 0, program_name, length, status )
  nargs = command_argument_count()

  ! Print out usage message if wrong number of arguments is specified
  if (nargs /= 3) then
    write(usage_message,*) "Usage: ",trim(program_name), &
      " <namelist filename> "                // &
      " test_XXX "                           // &
      " optional_test__arg "                 // &
      " with XXX in { "                      // &
      " init_string_err, "                   // &
      " copy_from_jedi_datetime_err, "       // &
      " add_duration_to_datetime, "          // &
      " duration_from_datetimes, "           // &
      " YYYYMMDD_to_JDN, "                   // &
      " JDN_to_YYYYMMDD_invalid, "           // &
      " hhmmss_to_seconds, "                 // &
      " seconds_to_hhmmss_large, "           // &
      " seconds_to_hhmmss_neg, "             // &
      " duration_init_bad_string_err, "      // &
      " duration_divide_zero_err, "          // &
      " duration_divide_remainder_err, "     // &
      " duration_divide_int_zero_err, "      // &
      " duration_divide_int_remainder_err "  // &
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
  case ("test_init_string_err")
    do_test_init_string_err = .true.
  case ("test_copy_from_jedi_datetime_err")
    do_test_copy_from_jedi_datetime_err = .true.
  case ("test_add_duration_to_datetime")
    do_test_add_duration_to_datetime = .true.
  case ("test_duration_from_datetimes")
    do_test_duration_from_datetimes = .true.
  case ("test_YYYYMMDD_to_JDN")
    do_test_YYYYMMDD_to_JDN = .true.
  case ("test_JDN_to_YYYYMMDD_invalid")
    do_test_JDN_to_YYYYMMDD_invalid = .true.
  case ("test_hhmmss_to_seconds")
    do_test_hhmmss_to_seconds = .true.
  case ("test_seconds_to_hhmmss_large")
    do_test_seconds_to_hhmmss_large = .true.
  case ("test_seconds_to_hhmmss_neg")
    do_test_seconds_to_hhmmss_neg = .true.
  case ("test_duration_init_bad_string_err")
    do_test_duration_init_bad_string_err = .true.
  case ("test_duration_divide_zero_err")
    do_test_duration_divide_zero_err = .true.
  case ("test_duration_divide_remainder_err")
    do_test_duration_divide_remainder_err = .true.
  case ("test_duration_divide_int_zero_err")
    do_test_duration_divide_int_zero_err = .true.
  case ("test_duration_divide_int_remainder_err")
    do_test_duration_divide_int_remainder_err = .true.
  case default
    call log_event( "Unknown test", LOG_LEVEL_ERROR )
  end select

  ! Setup configuration, and initialise tests
  call configuration%initialise( program_name, table_len=10 )
  call read_configuration( filename, configuration )
  call test_jedi_interface_init()

  if ( do_test_init_string_err ) then
    call run_init_string_err()
  end if
  if ( do_test_copy_from_jedi_datetime_err ) then
    call run_copy_from_jedi_datetime_err()
  end if
  if ( do_test_add_duration_to_datetime ) then
    call run_add_duration_to_datetime()
  end if
  if ( do_test_duration_from_datetimes ) then
    call run_duration_from_datetimes()
  end if
  if ( do_test_YYYYMMDD_to_JDN ) then
    call run_YYYYMMDD_to_JDN()
  end if
  if ( do_test_JDN_to_YYYYMMDD_invalid ) then
    call run_JDN_to_YYYYMMDD_invalid()
  end if
  if ( do_test_hhmmss_to_seconds ) then
    call run_hhmmss_to_seconds()
  end if
  if ( do_test_seconds_to_hhmmss_large ) then
    call run_seconds_to_hhmmss_large()
  end if
  if ( do_test_seconds_to_hhmmss_neg ) then
    call run_seconds_to_hhmmss_neg()
  end if
  if ( do_test_duration_init_bad_string_err ) then

    call get_command_argument( 3, dummy, length, status )
    allocate(character(length)::optional_arg)
    call get_command_argument( 3, optional_arg, length, status )

    write(usage_message,*) "test_duration_init_bad_string_err " // &
                           "with bad string: ", trim(optional_arg)
    call log_event( trim(usage_message), LOG_LEVEL_INFO )

    call run_duration_init_bad_string_err( optional_arg )

    deallocate(optional_arg)

  end if
  if ( do_test_duration_divide_zero_err ) then
    call run_duration_divide_zero_err()
  end if
  if ( do_test_duration_divide_remainder_err ) then
    call run_duration_divide_remainder_err()
  end if
  if ( do_test_duration_divide_int_zero_err ) then
    call run_duration_divide_int_zero_err()
  end if
  if ( do_test_duration_divide_int_remainder_err ) then
    call run_duration_divide_int_remainder_err()
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Finalise and close down

  call test_jedi_interface_final()

  if (allocated(program_name)) deallocate(program_name)
  if (allocated(filename))     deallocate(filename)
  if (allocated(test_flag))    deallocate(test_flag)

  call log_event( 'jedi-interface functional testing completed ...', LOG_LEVEL_INFO )

  ! Finalise halo functionality
  call finalise_halo_comms()

  ! Finalise the logging system
  call global_mpi%finalise()

  call final_configuration

  call destroy_comm()

end program jedi_lfric_time_test
