!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief A dummy method for use when running without evolving time-step. As the
!>        object holds no data, we can rely on the default constructor.

module no_timestep_alg_mod

  use log_mod,                          only: log_event,          &
                                              log_scratch_space,  &
                                              LOG_LEVEL_WARNING
  use timestep_method_mod,              only: timestep_method_type
  use driver_modeldb_mod,               only: modeldb_type

  implicit none

  type, extends(timestep_method_type), public :: no_timestep_type
    private

  contains
    private

    procedure, public  :: step
    procedure, public  :: finalise

  end type no_timestep_type

contains

  !> @brief Fulfils API without doing anything
  !> @param[in] modeldb Holds the model state
  subroutine step( self, modeldb )

    implicit none

    class(no_timestep_type), intent(inout) :: self

    type(modeldb_type),              intent(in), target  :: modeldb
    write( log_scratch_space, &
                    '(A, A)' ) 'CAUTION: Running with no timestepping. ' // &
                    ' Prognostic fields not evolved'
    call log_event( log_scratch_space, LOG_LEVEL_WARNING )

  end subroutine step

  !> @brief Fulfils API without doing anything
  subroutine finalise( self )
    implicit none

    class(no_timestep_type), intent(inout) :: self

  end subroutine finalise


end module no_timestep_alg_mod
