!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Abstract class for time-discretisation methods
!> @description An application can call a choice of time-stepping methods,
!>              and these methods store internal state. This module provides
!>              an abstract parent that can be included in the application's
!>              model database and point to the chosen time-stepping method.

module timestep_method_mod

  use field_mod,                    only: field_type
  use key_value_collection_mod,     only: key_value_collection_type
  use key_value_mod,                only: abstract_value_type
  use driver_modeldb_mod,           only: modeldb_type

  implicit none

  private
  public :: get_timestep_method_from_collection

  type, extends(abstract_value_type), public, abstract :: timestep_method_type
    private
  contains
    private
    procedure(step_interface), public,  deferred :: step
    procedure(final_interface), public, deferred :: finalise

  end type timestep_method_type

  abstract interface

    !> The API for the time-stepping procedure
    !> @param[in] modeldb Holds the model state

    subroutine step_interface(self, modeldb)
      import :: timestep_method_type
      import :: modeldb_type

      class(timestep_method_type),     intent(inout)       :: self
      type(modeldb_type),              intent(in), target  :: modeldb

    end subroutine step_interface

    !> Deallocates internal data structures
    subroutine final_interface(self)
      import :: timestep_method_type

      class(timestep_method_type), intent(inout) :: self

    end subroutine final_interface
  end interface

contains

  !-----------------------------------------------------------------------------
  ! Non-type-bound helper function
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> @brief Helper function to extract a timestep alg object from a
  !>        key-value collection
  !> @param[in] collection The key-value collection to extract from
  !> @param[in] name       The name of the timestep method object to extract
  !> @return    timestep_method   The requested timestep object
  function get_timestep_method_from_collection(collection, name) &
                                               result(timestep_method)

  implicit none

    type(key_value_collection_type), intent(in) :: collection
    character(*),                    intent(in) :: name

    class(timestep_method_type), pointer :: timestep_method
    class(abstract_value_type), pointer :: abstract_value

    call collection%get_value(trim(name), abstract_value)
    select type(abstract_value)
      class is (timestep_method_type)
        timestep_method => abstract_value
    end select

  end function get_timestep_method_from_collection

end module timestep_method_mod
