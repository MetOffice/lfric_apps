! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module scintelapi_interface_mod
!
! This module provides the necessary routines for initialising, configuring and
! finalising the API.
!

use log_mod,                   only: log_event, log_scratch_space,             &
                                     LOG_LEVEL_INFO, LOG_LEVEL_ERROR
use scintelapi_namelist_mod,   only: scintelapi_nl, required_lfric_namelists
use constants_def_mod,         only: field_kind_name_len, field_name_len,      &
                                     gen_id_len, genpar_len,                   &
                                     field_id_list_max_size, empty_string, rmdi
use lfricinp_setup_io_mod,     only: io_config
use lfricinp_datetime_mod,     only: datetime

implicit none

contains


subroutine scintelapi_initialise()
!
! This routine initialises all the necessary LFRic, XIOS, YAXT and API
! infrastructure.
!

use field_list_mod,            only: init_field_list
use generator_library_mod,     only: init_generator_lib
use dependency_graph_list_mod, only: init_dependency_graph_list
use lfricinp_lfric_driver_mod, only: lfricinp_initialise_lfric, model_clock,   &
                                     lfric_nl_fname
use lfricinp_setup_io_mod,     only: io_fname
use lfricinp_read_command_line_args_mod, only: lfricinp_read_command_line_args

implicit none

logical                             :: l_advance

! Read namelist file names from command line
call lfricinp_read_command_line_args(scintelapi_nl, lfric_nl_fname, io_fname)

! Set up IO file configuration
call io_config%load_namelist()

! Load date and time information
call datetime % initialise()

! Initialise LFRic infrastructure
call lfricinp_initialise_lfric(program_name_arg="scintelapi",                  &
     required_lfric_namelists = required_lfric_namelists,                      &
     start_date = datetime % first_validity_time,                              &
     time_origin = datetime % first_validity_time,                             &
     first_step = datetime % first_step,                                       &
     last_step = datetime % last_step,                                         &
     spinup_period = datetime % spinup_period,                                 &
     seconds_per_step = datetime % seconds_per_step)

! Advance clock to first time step, so output can be written to file
l_advance = model_clock % tick()
if (.not. l_advance) then
  call log_event('Failed to advance clock on initialisation', LOG_LEVEL_ERROR)
end if

! Initialise the field list
call init_field_list()

! Initialise the generator library
call init_generator_lib()

! Initialise the dependency graph list
call init_dependency_graph_list()

end subroutine scintelapi_initialise


subroutine scintelapi_finalise()
!
! This routine finalises all the used infrastructure
!

use lfricinp_lfric_driver_mod, only: lfricinp_finalise_lfric

implicit none

! Finalise LFRic infrastructure.
call lfricinp_finalise_lfric()

end subroutine scintelapi_finalise


subroutine scintelapi_add_field(field_id, field_kind, n_data, write_name)
!
! This routine is used to add a field to the internally stored global field
! list. It also performs several checks on the input for validity and
! consistency.
!

use finite_element_config_mod,      only: element_order_h, element_order_v
use function_space_collection_mod , only: function_space_collection
use fs_continuity_mod,              only: W3, Wtheta
use lfricinp_lfric_driver_mod,      only: mesh, twod_mesh
use field_list_mod,                 only: no_fields, field_list,               &
                                          field_io_name_list
use field_mod,                      only: field_proxy_type
use mesh_mod,                       only: mesh_type

implicit none

!
! Arguments
!
! Field identifier of new field to be added to list
character(len=*), optional, intent(in) :: field_id

! Field type identifier
character(len=*), optional, intent(in) :: field_kind

! Field non-spatial dimension size
integer, optional, intent(in)          :: n_data

! Field proxy
type(field_proxy_type)                 :: field_proxy

! XIOS write identifier to use as defined in iodef.xml file
character(len=*), optional, intent(in) :: write_name

!
! Local variables
!
! Mesh to use
type(mesh_type), pointer :: tmp_mesh => null()

! Function space to use
integer :: Fspace

! ndata
integer :: ndata

! Iterable
integer :: l

! Logicals used
logical :: l_field_id_exists, l_write_name_exists, l_field_id_present,         &
           l_write_name_present, l_field_kind_present
logical :: ndata_first

! User feedback
write(log_scratch_space,'(A)') 'Attempt to add ' // trim(field_id) //          &
                               ' to global field list.'
call log_event(log_scratch_space, LOG_LEVEL_INFO)

!
! START input checks ...
!

! Check field id is present and not already in global field list
l_field_id_present = .false.
if (present(field_id)) then
  if (trim(field_id) /= empty_string) l_field_id_present = .true.
end if

if (l_field_id_present) then

  l_field_id_exists = .false.
  do l = 1, no_fields
    if (trim(field_list(l)%get_name()) == trim(field_id)) then
      l_field_id_exists = .true.
      exit
    end if
  end do

  if (l_field_id_exists) then
    write(log_scratch_space,'(A)') trim(field_id) // ' already in field list'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

else

call log_event('No field id provided', LOG_LEVEL_ERROR)

end if

! Check if any fields in global field list already has the same write name
l_write_name_present = .false.
if (present(write_name)) then
  if (trim(write_name) /= empty_string) l_write_name_present = .true.
end if

if (l_write_name_present) then

  l_write_name_exists = .false.
  do l = 1, no_fields
    if (trim(field_io_name_list(l)) == trim(write_name)) then
      l_write_name_exists = .true.
      exit
    end if
  end do

  if (l_write_name_exists) then
    write(log_scratch_space,'(A)') 'The supplied write_name ' //               &
                                   trim(write_name) // ' for field ' //        &
                                   trim(field_id) // ' already exist for ' //  &
                                   'another field in the field list'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

end if

! Set size of non-spatial dimension based on input provided
if (present(n_data)) then
  ndata = n_data
else
  ndata = 1
end if

! Set mesh type, function space, etc based on field type provided
l_field_kind_present = .false.
if (present(field_kind)) then
  if (trim(field_kind) /= empty_string) l_field_kind_present = .true.
end if

if (l_field_kind_present) then

  select case (trim(field_kind))

    case('W3_field')
      tmp_mesh => mesh
      Fspace = W3
      ndata_first = .false.

    case('Wtheta_field')
      tmp_mesh => mesh
      Fspace = Wtheta
      ndata_first = .false.

    case('W3_field_2d')
      tmp_mesh => twod_mesh
      Fspace = W3
      ndata_first = .false.

    case('W3_soil_field')
      tmp_mesh => twod_mesh
      Fspace = W3
      ndata_first = .true.

    case DEFAULT
      write(log_scratch_space, '(A,A,A)')                                      &
         "Field type ", trim(field_kind), " not recognised"
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)

  end select

else

call log_event('No field type provided', LOG_LEVEL_ERROR)

end if

!
! Input checks DONE ...
!

! Set index of new field in the list
l = no_fields + 1

! Set new field
call field_list(l) % initialise(vector_space =                                  &
                                  function_space_collection%get_fs(             &
                                                               tmp_mesh,        &
                                                               element_order_h, &
                                                               element_order_v, &
                                                               Fspace,          &
                                                               ndata = ndata,   &
                                                               ndata_first =    &
                                                                 ndata_first    &
                                                                   ),           &
                                name = field_id)

! Initialise new field data to rmdi
field_proxy = field_list(l) % get_proxy()
field_proxy % data(:) = rmdi

! Set write id of the new field, if defined
if (l_write_name_present) then
  field_io_name_list(l) = write_name
end if

! Update the number of defined fields
no_fields = no_fields + 1

call log_event('Field successfully added', LOG_LEVEL_INFO)

nullify (tmp_mesh)

end subroutine scintelapi_add_field


subroutine scintelapi_add_dependency_graph(input_fields, output_fields,        &
                                          generator, genpar)
!
! This routine is used to add a dependency graph to the internally stored global
! dependency graph list. It also performs several checks on the input for
! validity and consistency.
!

use dependency_graph_list_mod, only: no_dependency_graphs, dependency_graph_list
use field_list_mod,            only: get_field_pointer, field_list, no_fields
use generator_library_mod,     only: generator_index, generator_list,          &
                                     no_generators

implicit none

!
! Arguments
!
! Identifier of generator to run
character(len=*), optional, intent(in) :: generator

! Input field id array
character(len=*), optional, intent(in) :: input_fields(:)

! Output field id array
character(len=*), optional, intent(in) :: output_fields(:)

! Optional input parameter list to the generator
character(len=*), optional, intent(in) :: genpar

!
! Local variables
!
! Sizes of input and output field arrays for the dep graph
integer :: no_input_fields, no_output_fields

! Iterable
integer :: i, j, l

! Logical to check if input and output fields are defined in global field list
logical :: l_field_does_not_exist

! Logical to check if generator is defined in generator library
logical :: l_generator_does_not_exist

! Other logicals used in checking input
logical :: l_input_fields_present, l_output_fields_present, l_generator_present

call log_event('Attempting to ADD DEPENDENCY GRAPH', LOG_LEVEL_INFO)

!
! START input checks ...
!

! Check the input and output field lists are minimal, i.e. no repeated field
! ids, the fields are already in the global field list, and the output field
! list has at least one field id.
l_input_fields_present = .false.
if (present(input_fields)) then
  if (size(input_fields) /= 0) l_input_fields_present = .true.
end if

if (l_input_fields_present) then
  write(log_scratch_space,'(100(A,1X))')                                       &
        'INPUT FIELDS:', (trim(input_fields(i)), i=1,size(input_fields))
  call log_event(log_scratch_space, LOG_LEVEL_INFO)
  do i = 1, size(input_fields)

    do j = 1, size(input_fields)
      if ( (trim(input_fields(j)) == trim(input_fields(i))) .and. j /= i ) then
        write(log_scratch_space, '(A)') 'Repeated field ids in input list'
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end if
    end do

    l_field_does_not_exist = .true.
    do j = 1, no_fields
      if (trim(field_list(j)%get_name()) == trim(input_fields(i))) then
        l_field_does_not_exist = .false.
        exit
      end if
    end do

    if (l_field_does_not_exist) then
      write(log_scratch_space, '(A)') 'Input field not defined'
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

  end do
end if
!
l_output_fields_present = .false.
if (present(output_fields)) then
  if (size(output_fields) /= 0) l_output_fields_present = .true.
end if

if (l_output_fields_present) then
  write(log_scratch_space,'(100(A,1X))')                                       &
        'OUTPUT FIELDS:', (trim(output_fields(i)), i=1,size(output_fields))
  call log_event(log_scratch_space, LOG_LEVEL_INFO)
  do i = 1, size(output_fields)

    do j = 1, size(output_fields)
      if ( (trim(output_fields(j)) == trim(output_fields(i)))                  &
          .and. j /= i ) then
        write(log_scratch_space, '(A)') 'Repeated field ids in output list'
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)
      end if
    end do

    l_field_does_not_exist = .true.
    do j = 1, no_fields
      if (trim(field_list(j)%get_name()) == trim(output_fields(i))) then
        l_field_does_not_exist = .false.
        exit
      end if
    end do

    if (l_field_does_not_exist) then
      write(log_scratch_space, '(A)') 'Output field not defined'
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if

  end do
else
  write(log_scratch_space, '(A)') 'No output fields in dependency graph!'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

! Check generator id has been provided and that it exists in the generator
! library
l_generator_present = .false.
if (present(generator)) then
  if (trim(generator) /= empty_string) l_generator_present = .true.
end if

if (l_generator_present) then
  write(log_scratch_space,'(2(A,1X))') 'GENERATOR:', trim(generator)
  call log_event(log_scratch_space, LOG_LEVEL_INFO)

  l_generator_does_not_exist = .true.
  do j = 1, no_generators
    if (trim(generator_list(j)%identifier) == trim(generator)) then
      l_generator_does_not_exist = .false.
      exit
    end if
  end do

  if (l_generator_does_not_exist) then
    write(log_scratch_space, '(A)') 'Generator not defined!'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

else
  write(log_scratch_space, '(A)') 'No generator id provided!'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

!
! Input checks DONE ...
!

! Set index of new dep graph in the list
l = no_dependency_graphs + 1

! Set up the pointers to the input fields
if (l_input_fields_present) then
  no_input_fields = size(input_fields)
  allocate(dependency_graph_list(l)%input_field(no_input_fields))
  do i = 1, no_input_fields
    dependency_graph_list(l)%input_field(i)%field_ptr =>                       &
                                            get_field_pointer(input_fields(i))
  end do
end if

! Set up the pointers to the output fields
no_output_fields = size(output_fields)
allocate(dependency_graph_list(l)%output_field(no_output_fields))
do i = 1, no_output_fields
  dependency_graph_list(l)%output_field(i)%field_ptr =>                        &
                                           get_field_pointer(output_fields(i))
end do

! Set up generator for this dep graph
dependency_graph_list(l)%gen = generator_list(generator_index(generator))

! Set up generator parameter list
if (present(genpar)) then
  dependency_graph_list(l)%genpar = genpar
else
  dependency_graph_list(l)%genpar = empty_string
end if

! Update the number of defined dep graphs
no_dependency_graphs = no_dependency_graphs + 1

call log_event('DEPENDENCY GRAPH successfully added', LOG_LEVEL_INFO)

end subroutine scintelapi_add_dependency_graph


subroutine scintelapi_add_fields_from_nl()
!
! This routine is used to add a list of fields to the internally stored global
! field list. Input is read from the scintelapi_nl namelist file as set in the
! module scintelapi_namelist_mod.
!

use lfricinp_unit_handler_mod,          only: get_free_unit

implicit none

! Unit number to use in namelist file reading.
integer :: nml_unit

! Input string declarations used in namelists
character(len=field_name_len)      :: field_id
character(len=field_kind_name_len) :: field_kind
integer                            :: n_data
character(len=field_name_len)      :: write_name
character(len=1)                   :: write_to_dump

! Field definition namelist
namelist /field_definitions/ field_id, field_kind, n_data, write_name,      &
                             write_to_dump

! Read namelist file for field definitions, and add said fields to internal
! field list
call get_free_unit(nml_unit)
open(unit=nml_unit, file=trim(scintelapi_nl))

do ! Loop over all field_definitions namelists

  ! Initialise field definition items
  field_id   = empty_string
  field_kind = empty_string
  n_data     = 1
  write_name = empty_string

  ! Read namelist items. Exit loop if EOF is reached
  read(unit=nml_unit, nml=field_definitions, end=101)

  ! Report warning to user if write_name item and write_to_dump settings are in
  ! in conflict with each other. In case of such a conflict, default position is
  ! not to write field to dump, i.e. write_name is an empty string.
  if ((trim(write_name) /= empty_string).and.(write_to_dump == 'n')) then
    write(log_scratch_space,'(A)') 'WARNING: Dump write name set, but '    //  &
                                   'writing to dump option turned off. '   //  &
                                   'Field will not be written to dump'
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
    write_name = empty_string
  end if

  if ((trim(write_name) == empty_string).and.(write_to_dump == 'y')) then
    write(log_scratch_space,'(A)') 'WARNING: Dump write name not set, but ' // &
                                   'writing to dump option turned on. '     // &
                                   'Field will not be written to dump'
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
  end if

  ! Add field to global field list
  call scintelapi_add_field(field_id      = field_id,                          &
                            field_kind    = field_kind,                        &
                            n_data        = n_data,                            &
                            write_name    = write_name)

end do ! End of loop over namelists

101 close(unit=nml_unit)

end subroutine scintelapi_add_fields_from_nl


subroutine scintelapi_add_dependency_graphs_from_nl()
!
! This routine is used to add a list of dependency graphs to the internally
! stored global dependency graph list. Input is read from the scintelapi_nl
! namelist file as set in the module scintelapi_namelist_mod.
!

use lfricinp_unit_handler_mod,          only: get_free_unit

implicit none

! Unit number to use in namelist file reading.
integer :: nml_unit

! Iterable(s)
integer :: i, k

! Allocatable arrays needed
character(len=field_name_len), allocatable :: input_fields_trimmed(:)
character(len=field_name_len), allocatable :: output_fields_trimmed(:)

! Input string declarations used in namelists
character(len=field_name_len)      :: output_fields(field_id_list_max_size)
character(len=field_name_len)      :: input_fields(field_id_list_max_size)
character(len=gen_id_len)          :: generator
character(len=genpar_len)          :: genpar

! Dependency graph definition namelist
namelist /dependency_graphs/ input_fields, output_fields, generator, genpar

! Read namelist file for dependency graph definitions, and add said dependency
! graphs to internal list
call get_free_unit(nml_unit)
open(unit=nml_unit, file=trim(scintelapi_nl))

do ! Loop over all dependency_graphs namelists

  ! Initialise dependency graph items. Assume input and output field name arrays
  ! are of maximum size, which will be trimmed down later.
  input_fields  = [(empty_string, i = 1,field_id_list_max_size)]
  output_fields = [(empty_string, i = 1,field_id_list_max_size)]
  generator     = empty_string
  genpar        = empty_string

  ! Read namelist items. Exit loop if EOF is reached
  read(unit=nml_unit, nml=dependency_graphs, end=102)

  ! Trim down input field name array to minimum size.
  !
  ! First find minimum size and allocate necessary array with that size.
  k = 0
  do i = 1, field_id_list_max_size
    if (trim(input_fields(i)) /= empty_string) k = k + 1
  end do
  allocate(input_fields_trimmed(k))
  !
  ! Second fill allocated minimun sized array with the necessary data.
  k = 0
  do i = 1, field_id_list_max_size
    if (trim(input_fields(i)) /= empty_string) then
      k = k + 1
      input_fields_trimmed(k) = input_fields(i)
    end if
  end do

  ! Trim down output field name array to minimum size.
  !
  ! First find minimum size and allocate necessary array with that size.
  k = 0
  do i = 1, field_id_list_max_size
    if (trim(output_fields(i)) /= empty_string) k = k + 1
  end do
  allocate(output_fields_trimmed(k))
  !
  ! Second fill allocated minimun sized array with the necessary data.
  k = 0
  do i = 1, field_id_list_max_size
    if (trim(output_fields(i)) /= empty_string) then
      k = k + 1
      output_fields_trimmed(k) = output_fields(i)
    end if
  end do

  ! Add dependendency graph to global dependency graph list
  call scintelapi_add_dependency_graph(input_fields  = input_fields_trimmed,   &
                                       output_fields = output_fields_trimmed,  &
                                       generator     = generator,              &
                                       genpar        = genpar)

  deallocate(input_fields_trimmed)
  deallocate(output_fields_trimmed)

end do ! End of loop over namelists

102 close(unit=nml_unit)

end subroutine scintelapi_add_dependency_graphs_from_nl

end module scintelapi_interface_mod
