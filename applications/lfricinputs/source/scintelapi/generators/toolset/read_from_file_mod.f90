! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module read_from_file_mod
!
! This module contains a generator to read a field from file.
!

use dependency_graph_mod, only: dependency_graph

implicit none

contains

subroutine read_from_file(dep_graph)
!
! This generator reads the output field in the dependency graph from file
!

use gen_io_check_mod,       only: gen_io_check
use lfric_xios_read_mod,    only: read_field_generic
use field_mod,              only: field_proxy_type
use field_parent_mod,       only: read_interface
use log_mod,                only: log_event, log_scratch_space, LOG_LEVEL_ERROR
use constants_def_mod,      only: genpar_len, field_name_len

implicit none

!
! Argument definitions:
!
! Dependency graph to be processed
class(dependency_graph), intent(in out) :: dep_graph

!
! Local variables
!
! IO procedure pointers
procedure(read_interface), pointer :: tmp_read_ptr
!
! Parameter list
character(len=genpar_len) :: parlist
!
! Id of field to read from file
character(len=field_name_len) field_io_name

! Error code for reading parameter list
integer :: ioerr

!
! Perform some initial input checks
!
call gen_io_check(                                                             &
                  dep_graph=dep_graph,                                         &
                  input_field_no=0,                                            &
                  output_field_no=1,                                           &
                  parameter_no=1                                               &
                 )
!
! Done with initial input checks
!

! Get id of field to read from file
parlist = dep_graph % genpar
read(parlist,'(A)',iostat=ioerr) field_io_name
if (ioerr /= 0 ) then
  write(log_scratch_space,'(A)') 'Error occured when parsing parameter ' //    &
                                 'list in routine READ_FROM_FILE. ' //         &
                                 'Input should be a single character string.'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

tmp_read_ptr => read_field_generic

! Read field from file
call dep_graph % output_field(1) % field_ptr %                                 &
                                   set_read_behaviour(tmp_read_ptr)
call dep_graph % output_field(1) % field_ptr %                                 &
                                   read_field('read_' // trim(field_io_name))

! Nullify read pointer
nullify(tmp_read_ptr)

end subroutine read_from_file

end module read_from_file_mod
