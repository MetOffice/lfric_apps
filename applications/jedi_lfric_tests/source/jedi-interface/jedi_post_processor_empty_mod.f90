!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief A module providing an empty implementation of the post processor
!>
!> @details This module includes a class that handles the post processing of
!>          state objects. This post post processor is an empty implementation
!>          that does nothing.
!>
module jedi_post_processor_empty_mod

  use jedi_post_processor_mod,  only : jedi_post_processor_type
  use jedi_state_mod,           only : jedi_state_type

  implicit none

  private

  type, public, extends(jedi_post_processor_type) :: &
                                          jedi_post_processor_empty_type
    private

  contains
    private

    !> Methods to process the data
    procedure, public :: pp_init
    procedure, public :: process
    procedure, public :: pp_final

  end type jedi_post_processor_empty_type

contains

  !> @brief    Calls the post processor initialise method.
  !>
  !> @param [inout] jedi_state The state to post process.
  subroutine pp_init( self, jedi_state )

    implicit none

    class(jedi_post_processor_empty_type), intent(inout) :: self
    type(jedi_state_type),                    intent(in) :: jedi_state

  end subroutine pp_init

  !> @brief    Calls the post processor process method with an empty method.
  !>
  !> @param [inout] jedi_state The state to post process.
  subroutine process( self, jedi_state )

    implicit none

    class(jedi_post_processor_empty_type), intent(inout) :: self
    type(jedi_state_type),                 intent(inout) :: jedi_state

  end subroutine process

  !> @brief    Calls the post processor finalise method.
  !>
  !> @param [inout] jedi_state The state to post process.
  subroutine pp_final( self, jedi_state )

    implicit none

    class(jedi_post_processor_empty_type), intent(inout) :: self
    type(jedi_state_type),                    intent(in) :: jedi_state

  end subroutine pp_final

end module jedi_post_processor_empty_mod
