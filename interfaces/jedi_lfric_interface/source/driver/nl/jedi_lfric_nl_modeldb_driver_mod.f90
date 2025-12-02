!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief    Module containing routines needed to initialise and finalise the
!>           modeldb object and also the model step method used by non-linear
!>           model from LFRic-JEDI
!>
!> @detail   The initialise and finalise methods set up/close down the modeldb
!>           as used by the JEDI interface to the non-linear model. This code uses
!>           the non-linear model as the basis to define the routines that need to
!>           be called to achieve this. The required code to setup the modeldb
!>           is contained in the non-linear program file located in:
!>
!>             applications/gungho_model/source/gungho.f90
!>
!>
!>           A step method is also included that calls the non-linear model step
!>           method and ticks the model clock.
!>
module jedi_lfric_nl_modeldb_driver_mod

  use constants_mod,                only : l_def
  use driver_config_mod,            only : init_config
  use driver_time_mod,              only : init_time, final_time
  use gungho_init_fields_mod,       only : finalise_model_data
  use gungho_mod,                   only : gungho_required_namelists
  use driver_modeldb_mod,           only : modeldb_type
  use gungho_model_mod,             only : finalise_infrastructure, &
                                           finalise_model
  use lfric_mpi_mod,                only : lfric_mpi_type
  use log_mod,                      only : log_event,         &
                                           log_scratch_space, &
                                           LOG_LEVEL_TRACE,   &
                                           LOG_LEVEL_ERROR
  use namelist_mod,                 only : namelist_type

  implicit none

  private
  public initialise_modeldb, step_nl, finalise_modeldb

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Sets up the modeldb and required infrastructure.
  !>
  !> @param [in]     modeldb_name An identifier given to the modeldb
  !> @param [in]     filename     The configuration file path and name
  !> @param [in]     mpi_obj      The mpi object to be associated with the modeldb
  !> @param [in,out] modeldb      The structure that holds model state
  subroutine initialise_modeldb( modeldb_name, filename, mpi_obj, modeldb )

    use gungho_driver_mod, only : gh_initialise => initialise

    implicit none

    character(*),                 intent(in) :: modeldb_name
    character(len=*),             intent(in) :: filename
    type(lfric_mpi_type), target, intent(in) :: mpi_obj
    type(modeldb_type),        intent(inout) :: modeldb
    ! Local
    character(len=*), parameter :: io_context_name = "gungho_atm"

    ! 1. Initialise modeldb field collections, configuration and mpi.
    modeldb%mpi => mpi_obj
    call modeldb%configuration%initialise( modeldb_name, table_len=10 )

    call modeldb%values%initialise('values', table_len = 5)

    ! 2. Create the depository, prognostics and diagnostics field collections
    call modeldb%fields%add_empty_field_collection( "depository",        &
                                                    table_len = 100 )
    call modeldb%fields%add_empty_field_collection( "prognostic_fields", &
                                                    table_len = 100)
    call modeldb%fields%add_empty_field_collection( "diagnostic_fields", &
                                                    table_len = 100 )
    call modeldb%fields%add_empty_field_collection("lbc_fields",         &
                                                   table_len = 100)
    call modeldb%fields%add_empty_field_collection("radiation_fields",   &
                                                    table_len = 100)
    call modeldb%fields%add_empty_field_collection("ls_fields",          &
                                                    table_len = 100)
    call modeldb%fields%add_empty_field_collection("fd_fields",          &
                                                    table_len = 100)

    call modeldb%io_contexts%initialise(modeldb_name, table_len=100)

    call init_config( filename, gungho_required_namelists, &
                      modeldb%configuration )

    ! 3. Initialise the clock and calendar
    call init_time( modeldb )

    ! 4. Call the gungho driver initialise
    call gh_initialise(modeldb_name, modeldb)

  end subroutine initialise_modeldb

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Increments the model state by a single timestep.
  !>
  !> @param [in,out] modeldb   The structure that holds model state
  subroutine step_nl( modeldb )

    use gungho_driver_mod, only : step

    implicit none

    type( modeldb_type ), intent(inout) :: modeldb

    ! Local
    logical( kind=l_def )          :: clock_running

    !. 1. Tick the clock and check its still running
    clock_running = modeldb%clock%tick()
    ! If the clock has finished then we need to abort and that may be due to
    ! incorrect model_clock configuration
    if ( .not. clock_running ) then
      write ( log_scratch_space, '(A)' ) &
              "jedi_lfric_linear_modeldb::The LFRic model clock has stopped."
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    ! 2. Step the non-linear model
    call step( modeldb )

  end subroutine step_nl

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Finalise the modeldb
  !> @param [in,out] modeldb      The structure that holds model state
  subroutine finalise_modeldb( modeldb )

    implicit none

    type(modeldb_type), intent(inout) :: modeldb

    call log_event( 'Finalising linear model modeldb', LOG_LEVEL_TRACE )

    ! Model configuration finalisation
    call finalise_model( modeldb )

    ! Destroy the fields stored in the modeldb model_data
    call finalise_model_data( modeldb )

    ! Finalise infrastructure and constants
    call finalise_infrastructure( modeldb )

    ! Finalise time
    call final_time( modeldb )

  end subroutine finalise_modeldb

end module jedi_lfric_nl_modeldb_driver_mod
