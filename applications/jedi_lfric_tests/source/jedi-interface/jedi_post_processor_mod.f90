!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief A module providing an abstract base class for post processors
!>
!> @details This module includes an abstract base class that handles the post
!>          processing of state objects. The post processor will accesses but
!>          not update the state.
!>
module jedi_post_processor_mod

  use jedi_state_mod,        only : jedi_state_type

  implicit none

  private

  type, public, abstract :: jedi_post_processor_type
    private
  contains
    private
    procedure(pp_init_interface),  public, deferred :: pp_init
    procedure(process_interface),  public, deferred :: process
    procedure(pp_final_interface), public, deferred :: pp_final
  end type jedi_post_processor_type

  abstract interface
    subroutine pp_init_interface( self, jedi_state )
      import jedi_post_processor_type, jedi_state_type
      implicit none
      class(jedi_post_processor_type), intent(inout) :: self
      type(jedi_state_type),              intent(in) :: jedi_state
    end subroutine pp_init_interface

    subroutine process_interface( self, jedi_state )
      import jedi_post_processor_type, jedi_state_type
      implicit none
      class(jedi_post_processor_type), intent(inout) :: self
      type(jedi_state_type),           intent(inout) :: jedi_state
    end subroutine process_interface

    subroutine pp_final_interface( self, jedi_state )
        import jedi_post_processor_type, jedi_state_type
        implicit none
        class(jedi_post_processor_type), intent(inout) :: self
        type(jedi_state_type),              intent(in) :: jedi_state
    end subroutine pp_final_interface

  end interface

end module jedi_post_processor_mod
