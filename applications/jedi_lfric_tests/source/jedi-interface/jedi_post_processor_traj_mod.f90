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

  use jedi_post_processor_mod, only : jedi_post_processor_type
  use jedi_state_mod,          only : jedi_state_type
  use jedi_linear_model_mod,   only : jedi_linear_model_type

  implicit none

  private

  type, public, extends(jedi_post_processor_type) :: jedi_post_processor_traj_type
    private

    !> Pointer to a linear model instance so the trajectory can be updated
    type(jedi_linear_model_type), pointer :: jedi_linear_model => null()

  contains
    private

    !> Initialise method that creates the post processor
    procedure, public :: initialise

    !> Methods to process the data
    procedure, public :: pp_init
    procedure, public :: process
    procedure, public :: pp_final

    !> Finalizer
    final              :: post_processor_traj_destructor

  end type jedi_post_processor_traj_type

contains

  !> @brief    Initialiser for the jedi_post_processor_traj_type.
  !>
  !> @param [inout] jedi_linear_model The linear model instance that the
  !>                                  post-processor will use
  subroutine initialise( self, jedi_linear_model )

  implicit none

  class( jedi_post_processor_traj_type ), intent(inout) :: self
  type( jedi_linear_model_type ), target, intent(inout) :: jedi_linear_model

    self%jedi_linear_model => jedi_linear_model

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

    call self%jedi_linear_model%set_trajectory( jedi_state )

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

    ! Nullify the linear model pointer
    nullify(self%jedi_linear_model)

  end subroutine post_processor_traj_destructor

end module jedi_post_processor_traj_mod
