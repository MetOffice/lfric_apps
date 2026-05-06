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

  use constants_mod,           only : i_def, l_def, str_def
  use namelist_collection_mod, only : namelist_collection_type

  implicit none

  private

type, public :: jedi_run_type
  private
  character(str_def)             :: jedi_run_name
  type(namelist_collection_type) :: configuration
  logical(kind=l_def)            :: timers_finalised

contains

  !> Run initialiser.
  procedure, public :: initialise

  !> LFRic initialiser.
  procedure, public :: initialise_infrastructure

  !> Get a pointer to the stored configuration.
  procedure, public ::  get_configuration

  !> Just finalise subroutine timing; to get useful timing statistics from failed adjoint tests
  procedure, public ::  finalise_timers

  !> Finalizer
  procedure, public :: finalise

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

  use lfric_mpi_mod,          only: create_comm, &
                                    lfric_comm_type
  use jedi_lfric_comm_mod,    only: init_external_comm

  implicit none

  class( jedi_run_type ), intent(inout) :: self
  character(len=*),       intent(in)    :: program_name
  integer(i_def),         intent(out)   :: out_communicator
  ! Local
  type(lfric_comm_type) :: lfric_comm
  integer(i_def) :: world_communicator

  self%jedi_run_name = program_name

  ! JEDI will initialise MPI so calling it here to enforce that behaviour.
  ! It will be called outside the scope of the model interface.
  call create_comm( lfric_comm )
  world_communicator = lfric_comm%get_comm_mpi_val()

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
  use namelist_mod,                  only: namelist_type
  use timing_mod,                    only: init_timing
  use io_config_mod,                 only: timer_output_path
  use jedi_lfric_tests_mod,          only: jedi_lfric_tests_required_namelists
  use lfric_mpi_mod,                 only: lfric_comm_type

  implicit none

  class( jedi_run_type ),         intent(inout) :: self
  character(len=*),               intent(in)    :: filename
  integer(i_def),                 intent(in)    :: model_communicator

  type(lfric_comm_type)                         :: lfric_comm
  type(namelist_type), pointer                  :: io_nml
  logical                                       :: lsubroutine_timers


  ! Initialise the configuration
  call self%configuration%initialise( self%jedi_run_name, table_len=10 )

  ! Initialise the model communicator to setup global_mpi
  call init_internal_comm( model_communicator )

  ! Setup the config which is curently global
  call init_config( filename, jedi_lfric_tests_required_namelists, &
                    self%configuration )

  ! Initialise the logger
  call lfric_comm%set_comm_mpi_val(model_communicator)
  call init_logger( lfric_comm, self%jedi_run_name )

  ! Initialise timing wrapper
  io_nml => self%configuration%get_namelist('io')
  call io_nml%get_value('subroutine_timers', lsubroutine_timers)
  call init_timing( lfric_comm, lsubroutine_timers, trim(self%jedi_run_name), timer_output_path )
  self%timers_finalised = .false.

  ! Initialise collections
  call init_collections()

end subroutine initialise_infrastructure

!> @brief    Get pointer to the stored configuration
!>
!> @return  configuration A pointer to the configuration
function get_configuration(self) result(configuration)

  class( jedi_run_type ), target, intent(inout) :: self
  type( namelist_collection_type ),     pointer :: configuration

  configuration => self%configuration

end function get_configuration

!> @brief    Just finalise subroutine timing; to get useful timing statistics from failed adjoint tests
!>
subroutine finalise_timers(self)

  use timing_mod,       only: final_timing

  implicit none

  class(jedi_run_type), intent(inout) :: self

  call final_timing( self%jedi_run_name )
  self%timers_finalised = .true.

end subroutine finalise_timers

!> @brief    Finalizer for jedi_run_type
!>
subroutine finalise(self)

  use driver_collections_mod,        only: final_collections
  use driver_config_mod,             only: final_config
  use driver_log_mod,                only: final_logger
  use timing_mod,                    only: final_timing
  use jedi_lfric_comm_mod,           only: final_external_comm, &
                                           final_internal_comm
  use lfric_mpi_mod,                 only: destroy_comm

  implicit none

  class(jedi_run_type), intent(inout) :: self

  ! Finalise collections
  call final_collections()

  ! Finalise subroutine timers
  if (.not. self%timers_finalised) call final_timing( self%jedi_run_name )

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

end subroutine finalise

end module jedi_run_mod
