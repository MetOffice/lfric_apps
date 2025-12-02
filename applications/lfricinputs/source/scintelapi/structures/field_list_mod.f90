! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module field_list_mod
!
! This module defines and provides access to the internal global field list of
! the API. It also provides routines for initialising the field list and
! creating a pointer to the field in the list.
!

use field_mod,         only: field_type
use log_mod,           only: log_event, LOG_LEVEL_ERROR
use constants_def_mod, only: field_name_len, max_no_fields

implicit none

! Field array containing all fields
type(field_type), target :: field_list(max_no_fields)

! Field write id array
character(len=field_name_len) :: field_io_name_list(max_no_fields)

! Actual number of fields in field list
integer :: no_fields

contains

subroutine init_field_list()
!
! Routine that initialises the global field list
!

implicit none

!
! Local variables
!
! Iterable
integer :: l

! Set number of field to zero
no_fields = 0

! Set associated field array that contains the field dump write names to blank
! strings.
do l = 1, max_no_fields
  field_io_name_list(l) = repeat(' ', field_name_len)
end do

end subroutine init_field_list


function get_field_pointer(field_id)
!
! Function that produces a pointer to a field with a given field id
!

implicit none

!
! Arguments
!
! Field identifier
character(len=*), intent(in) :: field_id

!
! Local variables
!
! Field index
integer :: field_index

! pointer to field in the in the field list
type(field_type), pointer :: get_field_pointer

! Iterable
integer :: l

! Loop over field list items to find field id and store its position index in
! the global field list.
field_index = 0
do l = 1, no_fields
  if (trim(field_id) == trim(field_list(l)%get_name())) then
    field_index = l
    exit
  end if
end do

! Create field pointer to the field in the global field list with the given id.
! If no field in global field list exist with the given id, report issue to
! user and abort the API.
if (field_index == 0) then ! No field in global field list has the given id

  get_field_pointer => null()
  call log_event('Field ' // trim(field_id) // ' is not in global field list', &
                 LOG_LEVEL_ERROR)

 else ! Field with given id has been found in the global field list

  get_field_pointer => field_list(field_index)

end if

end function get_field_pointer


function get_field_index(field_id)
!
! Function that produces the index in the global field list a given field_id
! corresponds to.
!

implicit none

!
! Arguments
!
! Field identifier
character(len=*), intent(in) :: field_id

!
! Local variables
!
! Field index
integer :: get_field_index

! Iterable
integer :: l

! Loop over field list items to find field id and store its position index in
! the global field list.
get_field_index = 0
do l = 1, no_fields
  if (trim(field_id) == trim(field_list(l)%get_name())) then
    get_field_index = l
    exit
  end if
end do

! If no field in global field list exist with the given id, report issue to
! user and abort the API.
if (get_field_index == 0) then ! No field in global field list has the given id

  call log_event('Field ' // trim(field_id) // ' is not in global field list', &
                 LOG_LEVEL_ERROR)

end if

end function get_field_index

end module field_list_mod
