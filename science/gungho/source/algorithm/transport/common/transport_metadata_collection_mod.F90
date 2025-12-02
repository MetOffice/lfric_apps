!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!>
!> @brief Holds and manages transport metadata objects created during a model run.
!>
!> @details A container which holds type definition of a collection of
!!          transport metadata objects. Provides functionality to return
!!          a metadata for a given field name.
module transport_metadata_collection_mod

  use constants_mod,          only : i_def
  use transport_metadata_mod, only : transport_metadata_type
  use log_mod,                only : log_event, log_scratch_space,    &
                                     LOG_LEVEL_ERROR
  use linked_list_mod,        only : linked_list_type,                &
                                     linked_list_item_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Type which is is a collection of metadata objects held in a linked list
  !-----------------------------------------------------------------------------
  type, public :: transport_metadata_collection_type
    private
    type(linked_list_type) :: metadata_list

    !> An unused allocatable integer that prevents an intenal compiler error
    !> with the Gnu Fortran compiler. Adding an allocatable forces the compiler
    !> to accept that the object has a finaliser. It gets confused without it.
    !> This is a workaround for GCC bug id 61767 - when this bug is fixed, the
    !> integer can be removed.
    integer(i_def), allocatable :: dummy_for_gnu

  contains
    procedure, public :: get_transport_metadata
    procedure, public :: set_transport_metadata
    procedure, public :: clear
    final             :: transport_metadata_collection_destructor
  end type transport_metadata_collection_type
  !-----------------------------------------------------------------------------

  interface transport_metadata_collection_type
    module procedure transport_metadata_collection_constructor
  end interface

  ! Module level variable to make the transport metadata collection
  ! globally available
  type(transport_metadata_collection_type), public, allocatable :: &
                                                transport_metadata_collection

contains
  !-----------------------------------------------------------------------------
  ! Construct the transport metadata collection
  !-----------------------------------------------------------------------------
  !> Function to construct a transport metadata collection

  function transport_metadata_collection_constructor() result(self)

    implicit none

    type(transport_metadata_collection_type) :: self

    self%metadata_list = linked_list_type()

  end function transport_metadata_collection_constructor


  !-----------------------------------------------------------------------------
  ! Get a metadata object
  !-----------------------------------------------------------------------------
  !> Function to get an instance of a transport_metadata from the linked list
  !> @param[in] fname Name of field (or group of fields) to get
  function get_transport_metadata( self,          &
                                   fname ) result(transport_metadata)

    use transport_config_mod, only: field_names

    implicit none

    class(transport_metadata_collection_type), intent(inout) :: self
    character(len=*),                          intent(in)    :: fname

    type(transport_metadata_type), pointer :: transport_metadata

    integer(i_def) :: i_name

    nullify(transport_metadata)

    if (.not.(any(trim(fname) == field_names))) then
      write(log_scratch_space, "(5A,*(A,:,', '))")                        &
      "Transport metadata (", trim(fname), ") is not defined in LFRic.",  &
      new_line("A"), "Available metadata types are: ",                    &
      (trim(field_names(i_name)), i_name = 1, size(field_names))
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

    transport_metadata => get_existing_metadata( self, fname )

  end function get_transport_metadata

  !-----------------------------------------------------------------------------
  ! Create a metadata object
  !-----------------------------------------------------------------------------
  !> Function to set an instance of a transport_metadata type in the linked list
  !> or create it if it doesn't exist
  !> @param[in] transport_metadata   Metadata object to store
  subroutine set_transport_metadata( self,  transport_metadata )

    implicit none

    class(transport_metadata_collection_type), intent(inout) :: self
    type(transport_metadata_type),             intent(in)    :: transport_metadata

    call self%metadata_list%insert_item(transport_metadata)

  end subroutine set_transport_metadata

  !-----------------------------------------------------------------------------
  ! Clear the transport metadata collection
  !-----------------------------------------------------------------------------
  !> Function to clear all items from the transport metadata collection
  !> linked list
  subroutine clear(self)

    implicit none

    class(transport_metadata_collection_type), intent(inout) :: self

    call self%metadata_list%clear()
    if (allocated(self%dummy_for_gnu)) deallocate(self%dummy_for_gnu)

  end subroutine clear

  !-----------------------------------------------------------------------------
  ! transport_metadata_collection destructor
  !-----------------------------------------------------------------------------

  subroutine transport_metadata_collection_destructor(self)

    implicit none

    type (transport_metadata_collection_type), intent(inout) :: self

    call self%clear()

  end subroutine transport_metadata_collection_destructor


  !------------------------------------------------------------------------------
  ! Private function (not accessible through the API - and only called from within
  ! this module) to scan the transport metadata collection for transport metadata with the
  ! given name and return a pointer to it. A null pointer is returned if the
  ! requested transport metadata does not exist.
  !
  ! param[in] fname Name of the field (or field group) to find
  ! return <pointer> Pointer to transport metadata object or null()
  function get_existing_metadata( self, &
                                  fname ) result(instance)

    implicit none

    class(transport_metadata_collection_type) :: self
    character(len=*),              intent(in) :: fname

    type(transport_metadata_type), pointer :: instance

    type(linked_list_item_type), pointer :: loop

    ! Point to head of the metadata linked list
    loop => self%metadata_list%get_head()

    ! Loop through the linked list
    do
      if ( .not. associated(loop) ) then
        ! Have reach the end of the list so either
        ! the list is empty or at the end of list.
        instance => null()

        loop => self%metadata_list%get_tail()
        call log_event('Unable to match metadata name', LOG_LEVEL_ERROR)
      end if

      ! Need to 'cast' the payload as the specific
      ! linked list data type, i.e. transport_metadata_type,
      ! before we can use it.
      select type(listmetadata => loop%payload)
      type is (transport_metadata_type)
        ! Check the properties of the payload in the current item
        ! to see if its the one being requested.
        if ( trim(fname) == trim(listmetadata%get_name()) ) then
          instance => listmetadata
          exit
        end if
      end select

      loop => loop%next
    end do

    nullify(loop)

  end function get_existing_metadata

end module transport_metadata_collection_mod
