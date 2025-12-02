!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Container for time axes.
!>
module gungho_time_axes_mod

  use constants_mod,                only : l_def
  use log_mod,                      only : log_event, LOG_LEVEL_ERROR
  use key_value_collection_mod,     only : key_value_collection_type
  use key_value_mod,                only : abstract_value_type
  use linked_list_mod,              only : linked_list_type
  use lfric_xios_time_axis_mod,     only : time_axis_type

  implicit none

  private
  public :: get_time_axes_from_collection

  !> Collection of time axes.
  !>
  type, extends(abstract_value_type), public :: gungho_time_axes_type

    private

    !> Time varying ancillaries time axis.
    type(linked_list_type), public :: ancil_times_list

    !> Time varying LBC time axis.
    type(linked_list_type), public :: lbc_times_list

    !> Pointer to the LBC time axis object
    type(time_axis_type), pointer, public :: lbc_time_axis

    !> Time varying linearisation state time axis.
    !>
    !> @todo Is this part of the linear model?
    !>
    type(linked_list_type), public :: ls_times_list

    !> Pointer to the linearisation state time axis (currently unused)
    type(time_axis_type), pointer, public :: ls_time_axis

  contains
    private

    procedure, public :: initialise
    procedure, public :: make_lbc_time_axis
    procedure, public :: save_lbc_time_axis

  end type gungho_time_axes_type

contains

  !-----------------------------------------------------------------------------
  ! Type-bound helper function
  !-----------------------------------------------------------------------------

  !> @brief Initialise time axes object
  !> @param[in] self      Time axes object
  subroutine initialise(self)
    implicit none
    class(gungho_time_axes_type), intent(inout) :: self

    self%lbc_time_axis => null()
    self%ls_time_axis => null()
  end subroutine initialise

  !> @brief Create an lbc time axis
  !> @param[in] self      Time axes object
  subroutine make_lbc_time_axis(self)
    implicit none
    class(gungho_time_axes_type), intent(inout) :: self

    logical(l_def),   parameter   :: cyclic=.false.
    logical(l_def),   parameter   :: interp_flag=.true.

    if (associated(self%lbc_time_axis)) &
      call log_event('attempt to recreate LBC time axis', LOG_LEVEL_ERROR)
    allocate(self%lbc_time_axis)
    call self%lbc_time_axis%initialise( &
            "lbc_time", file_id="lbc", yearly=cyclic, &
            interp_flag = interp_flag )
  end subroutine make_lbc_time_axis

  !> @brief Add the lbc time axis if any to the lbc times list
  !> @param[in,out] self      Time axes object
  subroutine save_lbc_time_axis(self)
    implicit none
    class(gungho_time_axes_type), intent(inout) :: self

    if (associated(self%lbc_time_axis)) &
      call self%lbc_times_list%insert_item(self%lbc_time_axis)
  end subroutine save_lbc_time_axis

  !-----------------------------------------------------------------------------
  ! Non-type-bound helper function
  !-----------------------------------------------------------------------------

  !-----------------------------------------------------------------------------
  !> @brief Helper function to extract a concrete time axes object from a
  !>        key-value collection
  !> @param[in] collection The key-value collection to extract from
  !> @param[in] name       The name of the time axes object to extract
  !> @return    fs_chain   The requested time_axes object
  function get_time_axes_from_collection(collection, name) result(time_axes)

  implicit none

    type(key_value_collection_type), intent(in) :: collection
    character(*),                    intent(in) :: name

    type(gungho_time_axes_type), pointer    :: time_axes

    class(abstract_value_type), pointer :: abstract_value

    call collection%get_value(trim(name), abstract_value)
    select type(abstract_value)
      type is (gungho_time_axes_type)
      time_axes => abstract_value
    end select

  end function get_time_axes_from_collection

end module gungho_time_axes_mod
