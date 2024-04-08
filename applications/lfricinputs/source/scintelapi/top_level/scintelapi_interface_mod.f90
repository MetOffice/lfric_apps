! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE scintelapi_interface_mod
!
! This module provides the necessary routines for initialising, configuring and
! finalising the API.
!

USE log_mod,                   ONLY: log_event, log_scratch_space,             &
                                     LOG_LEVEL_INFO, LOG_LEVEL_ERROR
USE scintelapi_namelist_mod,   ONLY: scintelapi_nl, required_lfric_namelists
USE constants_def_mod,         ONLY: field_kind_name_len, field_name_len,      &
                                     gen_id_len, genpar_len,                   &
                                     field_id_list_max_size, empty_string, rmdi
USE lfricinp_setup_io_mod,     ONLY: io_config
USE lfricinp_datetime_mod,     ONLY: datetime

IMPLICIT NONE

CONTAINS


SUBROUTINE scintelapi_initialise()
!
! This routine initialises all the necessary LFRic, XIOS, YAXT and API
! infrastructure.
!

USE field_list_mod,            ONLY: init_field_list
USE generator_library_mod,     ONLY: init_generator_lib
USE dependency_graph_list_mod, ONLY: init_dependency_graph_list
USE lfricinp_lfric_driver_mod, ONLY: lfricinp_initialise_lfric, model_clock,   &
                                     lfric_nl_fname
USE lfricinp_setup_io_mod,     ONLY: io_fname
USE lfricinp_read_command_line_args_mod, ONLY: lfricinp_read_command_line_args

IMPLICIT NONE

LOGICAL                             :: l_advance

! Read namelist file names from command line
CALL lfricinp_read_command_line_args(scintelapi_nl, lfric_nl_fname, io_fname)

! Set up IO file configuration
CALL io_config%load_namelist()

! Load date and time information
CALL datetime % initialise()

! Initialise LFRic infrastructure
CALL lfricinp_initialise_lfric(program_name_arg="scintelapi",                  &
     required_lfric_namelists = required_lfric_namelists,                      &
     start_date = datetime % first_validity_time,                              &
     time_origin = datetime % first_validity_time,                             &
     first_step = datetime % first_step,                                       &
     last_step = datetime % last_step,                                         &
     spinup_period = datetime % spinup_period,                                 &
     seconds_per_step = datetime % seconds_per_step)

! Advance clock to first time step, so output can be written to file
l_advance = model_clock % tick()
IF (.NOT. l_advance) THEN
  CALL log_event('Failed to advance clock on initialisation', LOG_LEVEL_ERROR)
END IF

! Initialise the field list
CALL init_field_list()

! Initialise the generator library
CALL init_generator_lib()

! Initialise the dependency graph list
CALL init_dependency_graph_list()

END SUBROUTINE scintelapi_initialise


SUBROUTINE scintelapi_finalise()
!
! This routine finalises all the used infrastructure
!

USE lfricinp_lfric_driver_mod, ONLY: lfricinp_finalise_lfric

IMPLICIT NONE

! Finalise LFRic infrastructure.
CALL lfricinp_finalise_lfric()

END SUBROUTINE scintelapi_finalise


SUBROUTINE scintelapi_add_field(field_id, field_kind, n_data, write_name)
!
! This routine is used to add a field to the internally stored global field
! list. It also performs several checks on the input for validity and
! consistency.
!

USE finite_element_config_mod,      ONLY: element_order
USE function_space_collection_mod , ONLY: function_space_collection
USE fs_continuity_mod,              ONLY: W3, Wtheta
USE lfricinp_lfric_driver_mod,      ONLY: mesh, twod_mesh
USE field_list_mod,                 ONLY: no_fields, field_list,               &
                                          field_io_name_list
USE field_mod,                      ONLY: field_proxy_type
USE mesh_mod,                       ONLY: mesh_type

IMPLICIT NONE

!
! Arguments
!
! Field identifier of new field to be added to list
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: field_id

! Field type identifier
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: field_kind

! Field non-spatial dimension size
INTEGER, OPTIONAL, INTENT(IN)          :: n_data

! Field proxy
TYPE(field_proxy_type)                 :: field_proxy

! XIOS write identifier to use as defined in iodef.xml file
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: write_name

!
! Local variables
!
! Mesh to use
TYPE(mesh_type), pointer :: tmp_mesh => null()

! Function space to use
INTEGER :: Fspace

! ndata
INTEGER :: ndata

! Iterable
INTEGER :: l

! Logicals used
LOGICAL :: l_field_id_exists, l_write_name_exists, l_field_id_present,         &
           l_write_name_present, l_field_kind_present
LOGICAL :: ndata_first

! User feedback
WRITE(log_scratch_space,'(A)') 'Attempt to add ' // TRIM(field_id) //          &
                               ' to global field list.'
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

!
! START input checks ...
!

! Check field id is present and not already in global field list
l_field_id_present = .false.
IF (PRESENT(field_id)) THEN
  IF (TRIM(field_id) /= empty_string) l_field_id_present = .true.
END IF

IF (l_field_id_present) THEN

  l_field_id_exists = .false.
  DO l = 1, no_fields
    IF (TRIM(field_list(l)%get_name()) == TRIM(field_id)) THEN
      l_field_id_exists = .true.
      EXIT
    END IF
  END DO

  IF (l_field_id_exists) THEN
    WRITE(log_scratch_space,'(A)') TRIM(field_id) // ' already in field list'
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END IF

ELSE

CALL log_event('No field id provided', LOG_LEVEL_ERROR)

END IF

! Check if any fields in global field list already has the same write name
l_write_name_present = .false.
IF (PRESENT(write_name)) THEN
  IF (TRIM(write_name) /= empty_string) l_write_name_present = .true.
END IF

IF (l_write_name_present) THEN

  l_write_name_exists = .false.
  DO l = 1, no_fields
    IF (TRIM(field_io_name_list(l)) == TRIM(write_name)) THEN
      l_write_name_exists = .true.
      EXIT
    END IF
  END DO

  IF (l_write_name_exists) THEN
    WRITE(log_scratch_space,'(A)') 'The supplied write_name ' //               &
                                   TRIM(write_name) // ' for field ' //        &
                                   TRIM(field_id) // ' already exist for ' //  &
                                   'another field in the field list'
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END IF

END IF

! Set size of non-spatial dimension based on input provided
IF (PRESENT(n_data)) THEN
  ndata = n_data
ELSE
  ndata = 1
END if

! Set mesh type, function space, etc based on field type provided
l_field_kind_present = .false.
IF (PRESENT(field_kind)) THEN
  IF (TRIM(field_kind) /= empty_string) l_field_kind_present = .true.
END IF

IF (l_field_kind_present) THEN

  SELECT CASE (TRIM(field_kind))

    CASE('W3_field')
      tmp_mesh => mesh
      Fspace = W3
      ndata_first = .FALSE.

    CASE('Wtheta_field')
      tmp_mesh => mesh
      Fspace = Wtheta
      ndata_first = .FALSE.

    CASE('W3_field_2d')
      tmp_mesh => twod_mesh
      Fspace = W3
      ndata_first = .FALSE.

    CASE('W3_soil_field')
      tmp_mesh => twod_mesh
      Fspace = W3
      ndata_first = .TRUE.

    CASE DEFAULT
      WRITE(log_scratch_space, '(A,A,A)')                                      &
         "Field type ", TRIM(field_kind), " not recognised"
      CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)

  END SELECT

ELSE

CALL log_event('No field type provided', LOG_LEVEL_ERROR)

END IF

!
! Input checks DONE ...
!

! Set index of new field in the list
l = no_fields + 1

! Set new field
CALL field_list(l) % initialise(vector_space =                                 &
                                function_space_collection%get_fs(              &
                                                                tmp_mesh,      &
                                                                element_order, &
                                                                Fspace,        &
                                                                ndata = ndata, &
                                                                ndata_first =  &
                                                                  ndata_first  &
                                                                ),             &
                                name         = field_id)

! Initialise new field data to rmdi
field_proxy = field_list(l) % get_proxy()
field_proxy % data(:) = rmdi

! Set write id of the new field, if defined
IF (l_write_name_present) THEN
  field_io_name_list(l) = write_name
END IF

! Update the number of defined fields
no_fields = no_fields + 1

CALL log_event('Field successfully added', LOG_LEVEL_INFO)

NULLIFY (tmp_mesh)

END SUBROUTINE scintelapi_add_field


SUBROUTINE scintelapi_add_dependency_graph(input_fields, output_fields,        &
                                          generator, genpar)
!
! This routine is used to add a dependency graph to the internally stored global
! dependency graph list. It also performs several checks on the input for
! validity and consistency.
!

USE dependency_graph_list_mod, ONLY: no_dependency_graphs, dependency_graph_list
USE field_list_mod,            ONLY: get_field_pointer, field_list, no_fields
USE generator_library_mod,     ONLY: generator_index, generator_list,          &
                                     no_generators

IMPLICIT NONE

!
! Arguments
!
! Identifier of generator to run
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: generator

! Input field id array
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: input_fields(:)

! Output field id array
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: output_fields(:)

! Optional input parameter list to the generator
CHARACTER(LEN=*), OPTIONAL, INTENT(IN) :: genpar

!
! Local variables
!
! Sizes of input and output field arrays for the dep graph
INTEGER :: no_input_fields, no_output_fields

! Iterable
INTEGER :: i, j, l

! Logical to check if input and output fields are defined in global field list
LOGICAL :: l_field_does_not_exist

! Logical to check if generator is defined in generator library
LOGICAL :: l_generator_does_not_exist

! Other logicals used in checking input
LOGICAL :: l_input_fields_present, l_output_fields_present, l_generator_present

CALL log_event('Attempting to ADD DEPENDENCY GRAPH', LOG_LEVEL_INFO)

!
! START input checks ...
!

! Check the input and output field lists are minimal, i.e. no repeated field
! ids, the fields are already in the global field list, and the output field
! list has at least one field id.
l_input_fields_present = .false.
IF (PRESENT(input_fields)) THEN
  IF (SIZE(input_fields) /= 0) l_input_fields_present = .true.
END IF

IF (l_input_fields_present) THEN
  WRITE(log_scratch_space,'(100(A,1X))')                                       &
        'INPUT FIELDS:', (TRIM(input_fields(i)), i=1,SIZE(input_fields))
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
  DO i = 1, SIZE(input_fields)

    DO j = 1, SIZE(input_fields)
      IF ( (TRIM(input_fields(j)) == TRIM(input_fields(i))) .AND. j /= i ) THEN
        WRITE(log_scratch_space, '(A)') 'Repeated field ids in input list'
        CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
      END IF
    END DO

    l_field_does_not_exist = .true.
    DO j = 1, no_fields
      IF (TRIM(field_list(j)%get_name()) == TRIM(input_fields(i))) THEN
        l_field_does_not_exist = .false.
        EXIT
      END IF
    END DO

    IF (l_field_does_not_exist) THEN
      WRITE(log_scratch_space, '(A)') 'Input field not defined'
      CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
    END IF

  END DO
END IF
!
l_output_fields_present = .false.
IF (PRESENT(output_fields)) THEN
  IF (SIZE(output_fields) /= 0) l_output_fields_present = .true.
END IF

IF (l_output_fields_present) THEN
  WRITE(log_scratch_space,'(100(A,1X))')                                       &
        'OUTPUT FIELDS:', (TRIM(output_fields(i)), i=1,SIZE(output_fields))
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
  DO i = 1, SIZE(output_fields)

    DO j = 1, SIZE(output_fields)
      IF ( (TRIM(output_fields(j)) == TRIM(output_fields(i)))                  &
          .AND. j /= i ) THEN
        WRITE(log_scratch_space, '(A)') 'Repeated field ids in output list'
        CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
      END IF
    END DO

    l_field_does_not_exist = .true.
    DO j = 1, no_fields
      IF (TRIM(field_list(j)%get_name()) == TRIM(output_fields(i))) THEN
        l_field_does_not_exist = .false.
        EXIT
      END IF
    END DO

    IF (l_field_does_not_exist) THEN
      WRITE(log_scratch_space, '(A)') 'Output field not defined'
      CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
    END IF

  END DO
ELSE
  WRITE(log_scratch_space, '(A)') 'No output fields in dependency graph!'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

! Check generator id has been provided and that it exists in the generator
! library
l_generator_present = .false.
IF (PRESENT(generator)) THEN
  IF (TRIM(generator) /= empty_string) l_generator_present = .true.
END IF

IF (l_generator_present) THEN
  WRITE(log_scratch_space,'(2(A,1X))') 'GENERATOR:', TRIM(generator)
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

  l_generator_does_not_exist = .true.
  DO j = 1, no_generators
    IF (TRIM(generator_list(j)%identifier) == TRIM(generator)) THEN
      l_generator_does_not_exist = .false.
      EXIT
    END IF
  END DO

  IF (l_generator_does_not_exist) THEN
    WRITE(log_scratch_space, '(A)') 'Generator not defined!'
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END IF

ELSE
  WRITE(log_scratch_space, '(A)') 'No generator id provided!'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

!
! Input checks DONE ...
!

! Set index of new dep graph in the list
l = no_dependency_graphs + 1

! Set up the pointers to the input fields
IF (l_input_fields_present) THEN
  no_input_fields = SIZE(input_fields)
  ALLOCATE(dependency_graph_list(l)%input_field(no_input_fields))
  DO i = 1, no_input_fields
    dependency_graph_list(l)%input_field(i)%field_ptr =>                       &
                                            get_field_pointer(input_fields(i))
  END DO
END IF

! Set up the pointers to the output fields
no_output_fields = SIZE(output_fields)
ALLOCATE(dependency_graph_list(l)%output_field(no_output_fields))
DO i = 1, no_output_fields
  dependency_graph_list(l)%output_field(i)%field_ptr =>                        &
                                           get_field_pointer(output_fields(i))
END DO

! Set up generator for this dep graph
dependency_graph_list(l)%gen = generator_list(generator_index(generator))

! Set up generator parameter list
IF (PRESENT(genpar)) THEN
  dependency_graph_list(l)%genpar = genpar
ELSE
  dependency_graph_list(l)%genpar = empty_string
END IF

! Update the number of defined dep graphs
no_dependency_graphs = no_dependency_graphs + 1

CALL log_event('DEPENDENCY GRAPH successfully added', LOG_LEVEL_INFO)

END SUBROUTINE scintelapi_add_dependency_graph


SUBROUTINE scintelapi_add_fields_from_nl()
!
! This routine is used to add a list of fields to the internally stored global
! field list. Input is read from the scintelapi_nl namelist file as set in the
! module scintelapi_namelist_mod.
!

USE lfricinp_unit_handler_mod,          ONLY: get_free_unit

IMPLICIT NONE

! Unit number to use in namelist file reading.
INTEGER :: nml_unit

! Input string declarations used in namelists
CHARACTER(LEN=field_name_len)      :: field_id
CHARACTER(LEN=field_kind_name_len) :: field_kind
INTEGER                            :: n_data
CHARACTER(LEN=field_name_len)      :: write_name
CHARACTER(LEN=1)                   :: write_to_dump

! Field definition namelist
NAMELIST /field_definitions/ field_id, field_kind, n_data, write_name,      &
                             write_to_dump

! Read namelist file for field definitions, and add said fields to internal
! field list
CALL get_free_unit(nml_unit)
OPEN(unit=nml_unit, file=TRIM(scintelapi_nl))

DO ! Loop over all field_definitions namelists

  ! Initialise field definition items
  field_id   = empty_string
  field_kind = empty_string
  n_data     = 1
  write_name = empty_string

  ! Read namelist items. Exit loop if EOF is reached
  READ(unit=nml_unit, nml=field_definitions, END=101)

  ! Report warning to user if write_name item and write_to_dump settings are in
  ! in conflict with each other. In case of such a conflict, default position is
  ! not to write field to dump, i.e. write_name is an empty string.
  IF ((TRIM(write_name) /= empty_string).AND.(write_to_dump == 'n')) THEN
    WRITE(log_scratch_space,'(A)') 'WARNING: Dump write name set, but '    //  &
                                   'writing to dump option turned off. '   //  &
                                   'Field will not be written to dump'
    CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
    write_name = empty_string
  END IF

  IF ((TRIM(write_name) == empty_string).AND.(write_to_dump == 'y')) THEN
    WRITE(log_scratch_space,'(A)') 'WARNING: Dump write name not set, but ' // &
                                   'writing to dump option turned on. '     // &
                                   'Field will not be written to dump'
    CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
  END IF

  ! Add field to global field list
  CALL scintelapi_add_field(field_id      = field_id,                          &
                            field_kind    = field_kind,                        &
                            n_data        = n_data,                            &
                            write_name    = write_name)

END DO ! End of loop over namelists

101 CLOSE(unit=nml_unit)

END SUBROUTINE scintelapi_add_fields_from_nl


SUBROUTINE scintelapi_add_dependency_graphs_from_nl()
!
! This routine is used to add a list of dependency graphs to the internally
! stored global dependency graph list. Input is read from the scintelapi_nl
! namelist file as set in the module scintelapi_namelist_mod.
!

USE lfricinp_unit_handler_mod,          ONLY: get_free_unit

IMPLICIT NONE

! Unit number to use in namelist file reading.
INTEGER :: nml_unit

! Iterable(s)
INTEGER :: i, k

! Allocatable arrays needed
CHARACTER(LEN=field_name_len), ALLOCATABLE :: input_fields_trimmed(:)
CHARACTER(LEN=field_name_len), ALLOCATABLE :: output_fields_trimmed(:)

! Input string declarations used in namelists
CHARACTER(LEN=field_name_len)      :: output_fields(field_id_list_max_size)
CHARACTER(LEN=field_name_len)      :: input_fields(field_id_list_max_size)
CHARACTER(LEN=gen_id_len)          :: generator
CHARACTER(LEN=genpar_len)          :: genpar

! Dependency graph definition namelist
NAMELIST /dependency_graphs/ input_fields, output_fields, generator, genpar

! Read namelist file for dependency graph definitions, and add said dependency
! graphs to internal list
CALL get_free_unit(nml_unit)
OPEN(unit=nml_unit, file=TRIM(scintelapi_nl))

DO ! Loop over all dependency_graphs namelists

  ! Initialise dependency graph items. Assume input and output field name arrays
  ! are of maximum size, which will be trimmed down later.
  input_fields  = [(empty_string, i = 1,field_id_list_max_size)]
  output_fields = [(empty_string, i = 1,field_id_list_max_size)]
  generator     = empty_string
  genpar        = empty_string

  ! Read namelist items. Exit loop if EOF is reached
  READ(unit=nml_unit, nml=dependency_graphs, END=102)

  ! Trim down input field name array to minimum size.
  !
  ! First find minimum size and allocate necessary array with that size.
  k = 0
  DO i = 1, field_id_list_max_size
    IF (TRIM(input_fields(i)) /= empty_string) k = k + 1
  END DO
  ALLOCATE(input_fields_trimmed(k))
  !
  ! Second fill allocated minimun sized array with the necessary data.
  k = 0
  DO i = 1, field_id_list_max_size
    IF (TRIM(input_fields(i)) /= empty_string) THEN
      k = k + 1
      input_fields_trimmed(k) = input_fields(i)
    END IF
  END DO

  ! Trim down output field name array to minimum size.
  !
  ! First find minimum size and allocate necessary array with that size.
  k = 0
  DO i = 1, field_id_list_max_size
    IF (TRIM(output_fields(i)) /= empty_string) k = k + 1
  END DO
  ALLOCATE(output_fields_trimmed(k))
  !
  ! Second fill allocated minimun sized array with the necessary data.
  k = 0
  DO i = 1, field_id_list_max_size
    IF (TRIM(output_fields(i)) /= empty_string) THEN
      k = k + 1
      output_fields_trimmed(k) = output_fields(i)
    END IF
  END DO

  ! Add dependendency graph to global dependency graph list
  CALL scintelapi_add_dependency_graph(input_fields  = input_fields_trimmed,   &
                                       output_fields = output_fields_trimmed,  &
                                       generator     = generator,              &
                                       genpar        = genpar)

  DEALLOCATE(input_fields_trimmed)
  DEALLOCATE(output_fields_trimmed)

END DO ! End of loop over namelists

102 CLOSE(unit=nml_unit)

END SUBROUTINE scintelapi_add_dependency_graphs_from_nl

END MODULE scintelapi_interface_mod
