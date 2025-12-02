!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief A module providing the trajectory implementation of the post processor
!>
!> @details This module includes a class that handles the post processing of
!>          state objects. This post processor stores a pointer to the linear
!>          model instance and calls the set_trajectory method to store the
!>          linear state. The linear state is created by running the non-linear
!>          forecast model.
!>
module jedi_post_processor_traj_mod

  use constants_mod,              only : i_def
  use jedi_base_linear_model_mod, only : jedi_base_linear_model_type
  use jedi_post_processor_mod,    only : jedi_post_processor_type
  use jedi_state_mod,             only : jedi_state_type
  use jedi_linear_model_mod,      only : jedi_linear_model_type
  use jedi_id_linear_model_mod,   only : jedi_id_linear_model_type
  use log_mod,                    only : log_event,          &
                                         log_scratch_space,  &
                                         LOG_LEVEL_ERROR

  implicit none

  private

  ! Enumerator types to define the linear model to use
  integer(kind=i_def), parameter :: type_full = 5
  integer(kind=i_def), parameter :: type_identity = 8
  integer(kind=i_def), parameter :: type_none = 43

  type, public, extends(jedi_post_processor_type) :: jedi_post_processor_traj_type
    private

    !> Pointer to a linear model instance
    type(jedi_linear_model_type),    pointer :: jedi_linear_model => null()
    type(jedi_id_linear_model_type), pointer :: jedi_id_linear_model => null()

    !> Enumerator that holds the enrolled model type
    integer(kind=i_def) :: type_of_linear_model = type_none

  contains
    private

    !> Initialise method that creates the post processor
    procedure, public :: initialise

    !> Methods to process the data
    procedure, public :: pp_init
    procedure, public :: process
    procedure, public :: pp_final

    !> Finalizer
    final             :: post_processor_traj_destructor

  end type jedi_post_processor_traj_type

contains

  !> @brief    Initialiser for the jedi_post_processor_traj_type.
  !>
  !> @param [inout] jedi_linear_model The linear model instance that the
  !>                                  post-processor will use
  subroutine initialise( self, jedi_linear_model )

  implicit none

  class( jedi_post_processor_traj_type ),       intent(inout) :: self
  class( jedi_base_linear_model_type ), target, intent(inout) :: jedi_linear_model


  ! Select the correct extened class
  select type(jedi_linear_model)
    type is ( jedi_linear_model_type )
      self%jedi_linear_model => jedi_linear_model
      self%type_of_linear_model = type_full
    type is ( jedi_id_linear_model_type )
      self%jedi_id_linear_model => jedi_linear_model
      self%type_of_linear_model = type_identity
    class default
      call log_event( "The supplied linear model is not supported", &
                      LOG_LEVEL_ERROR )
    end select

  end subroutine initialise

  !> @brief    Calls the post processor initialise method.
  !>
  !> @param [inout] jedi_state The state to post process.
  subroutine pp_init( self, jedi_state )

    implicit none

    class(jedi_post_processor_traj_type), intent(inout) :: self
    type(jedi_state_type),                   intent(in) :: jedi_state

  end subroutine pp_init

  !> @brief    Calls the post processor process method to store current state
  !>           in trajectory.
  !>
  !> @param [inout] jedi_state The state to post process.
  subroutine process( self, jedi_state )

    implicit none

    class(jedi_post_processor_traj_type), intent(inout) :: self
    type(jedi_state_type),                intent(inout) :: jedi_state

    select case(self%type_of_linear_model)
    case(type_full)
      call self%jedi_linear_model%set_trajectory( jedi_state )
    case(type_identity)
      call self%jedi_id_linear_model%set_trajectory( jedi_state )
    case default
      call log_event( 'The linear model has not been set-up', LOG_LEVEL_ERROR )
    end select

  end subroutine process

  !> @brief    Calls the post processor finalise method.
  !>
  !> @param [inout] jedi_state The state to post process.
  subroutine pp_final( self, jedi_state )

    implicit none

    class(jedi_post_processor_traj_type), intent(inout) :: self
    type(jedi_state_type),                   intent(in) :: jedi_state

  end subroutine pp_final

  !> @brief    Finalize the jedi_post_processor_traj_type
  !>
  subroutine post_processor_traj_destructor(self)

    implicit none

    type(jedi_post_processor_traj_type), intent(inout) :: self

    ! Nullify the linear model pointers
    nullify(self%jedi_linear_model)
    nullify(self%jedi_id_linear_model)
    ! Reset the enumerator
    self%type_of_linear_model = type_none

  end subroutine post_processor_traj_destructor

end module jedi_post_processor_traj_mod
