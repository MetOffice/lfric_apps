!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing an class that handles LFRic initialisation.
!>
!> @details This class handles the initialisation and finalisation of LFRic
!
module jedi_run_mod

  use constants_mod,           only : i_def, str_def
  use namelist_collection_mod, only : namelist_collection_type

  implicit none

  private

type, public :: jedi_run_type
  private
  character(str_def)             :: jedi_run_name
  type(namelist_collection_type) :: configuration

contains

  !> Run initialiser.
  procedure, public :: initialise

  !> LFRic initialiser.
  procedure, public :: initialise_infrastructure

  !> Finalizer
  final             :: jedi_run_destructor

end type jedi_run_type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @brief Initialiser for jedi_run_type
!>
!> @param [in]  program_name     The model name
!> @param [out] out_communicator The communicator to be used by the application
subroutine initialise( self, program_name, out_communicator )

  use mpi_mod,                only: create_comm
  use jedi_lfric_comm_mod,    only: init_external_comm

  implicit none

  class( jedi_run_type ), intent(inout) :: self
  character(len=*),       intent(in)    :: program_name
  integer(i_def),         intent(out)   :: out_communicator

  ! Local
  integer(i_def) :: world_communicator

  self%jedi_run_name = program_name

  ! JEDI will initialise MPI so calling it here to enforce that behaviour.
  ! It will be called outside the scope of the model interface.
  call create_comm( world_communicator )

  ! Call to initialise external dependencies like XIOS that require the world
  ! comm
  call init_external_comm( program_name, world_communicator, out_communicator)

end subroutine initialise

!> @brief    Initialiser for LFRic infrastructure
!>
!> @param [in]    filename           A character that contains the
!>                                   location of the namelist file.
!> @param [in]    model_communicator The communicator used by the model.
subroutine initialise_infrastructure( self, filename, model_communicator )

  use jedi_lfric_comm_mod,           only: init_internal_comm
  use driver_collections_mod,        only: init_collections
  use driver_config_mod,             only: init_config
  use driver_log_mod,                only: init_logger
  use jedi_lfric_tests_mod,          only: jedi_lfric_tests_required_namelists
  use namelist_collection_mod,       only: namelist_collection_type

  implicit none

  class( jedi_run_type ),         intent(inout) :: self
  character(len=*),               intent(in)    :: filename
  integer(i_def),                 intent(in)    :: model_communicator

  ! Initialise the configuration
  call self%configuration%initialise( self%jedi_run_name, table_len=10 )

  ! Initialise the model communicator to setup global_mpi
  call init_internal_comm( model_communicator )

  ! Setup the config which is curently global
  call init_config( filename, jedi_lfric_tests_required_namelists, &
                    self%configuration )

  ! Initialise the logger
  call init_logger( model_communicator, self%jedi_run_name )

  ! Initialise collections
  call init_collections()

end subroutine initialise_infrastructure

!> @brief    Finalizer for jedi_run_type
!>
subroutine jedi_run_destructor(self)

  use driver_collections_mod,        only: final_collections
  use driver_config_mod,             only: final_config
  use driver_log_mod,                only: final_logger
  use jedi_lfric_comm_mod,           only: final_external_comm, &
                                           final_internal_comm
  use mpi_mod,                       only: destroy_comm

  implicit none

  type(jedi_run_type), intent(inout) :: self

  ! Finalise collections
  call final_collections()

  ! Finalise logger
  call final_logger(self%jedi_run_name)

  ! Finalise the config
  call final_config()

  ! Finalise internal communicator groups
  call final_internal_comm()

  ! Finalise external communicator groups
  call final_external_comm()

  ! Finalise the communicator
  call destroy_comm()

end subroutine jedi_run_destructor

end module jedi_run_mod
