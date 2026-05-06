!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Subroutines for determining fields to initialise from iodef and the
!>        routine for field creation.
module lfric2lfric_field_init_mod

  use constants_mod,                      only: i_def, str_def
  use field_collection_mod,               only: field_collection_type
  use field_mod,                          only: field_type
  use field_parent_mod,                   only: read_interface,            &
                                                write_interface,           &
                                                checkpoint_read_interface, &
                                                checkpoint_write_interface
  use file_mod,                           only: FILE_MODE_READ, &
                                                FILE_OP_OPEN
  use function_space_mod,                 only: function_space_type
  use function_space_collection_mod,      only: function_space_collection
  use lfric_ncdf_file_mod,                only: lfric_ncdf_file_type
  use lfric_xios_diag_mod,                only: field_is_valid
  use lfric_xios_read_mod,                only: read_field_generic, &
                                                checkpoint_read_xios
  use lfric_xios_write_mod,               only: write_field_generic, &
                                                checkpoint_write_xios
  use lfric2lfric_config_mod,             only: mode_ics, mode_lbc
  use log_mod,                            only: log_event,         &
                                                log_scratch_space, &
                                                log_level_info,    &
                                                log_level_trace,   &
                                                log_level_error
  use mesh_mod,                           only: mesh_type
  use netcdf,                             only: nf90_inquire_variable, &
                                                nf90_inquire,          &
                                                nf90_max_name
  use space_from_metadata_mod,            only: space_from_metadata

  implicit none

  private

  public :: get_field_list, field_maker

contains

  !> @brief   Finds all the valid fields in a dump/checkpoint/ancillary file.
  !> @details It is probably not the best place for this subroutine...
  !!          It would be better to use a linked list rather than a vector
  !!          to store field names, but as it requires quite a bit of coding
  !!          it was left for later. The only disadvantage of using a vector
  !!          is that we will be wasting some memory, or having two loops
  !!          instead of one and therefore wasting some computational time.
  !!          TODO : #325 Define fields to initialise for lfric2lfric miniapp
  !> @param [out] num_fields    Number of valid fields in the input file
  !> @param [out] config_list   Vector of valid field names in the input file
  !> @param [in]  file_name     Input file to be checked
  !> @param [in]  prefix        Variable name prefix
  subroutine get_field_list(num_fields, config_list, file_name, prefix)

    implicit none

    integer(kind=i_def),                 intent(out) :: num_fields
    character(len=str_def), allocatable, intent(out) :: config_list(:)
    character(len=str_def),              intent(in)  :: file_name
    character(len=nf90_max_name),        intent(in)  :: prefix
    ! Local variables
    type(lfric_ncdf_file_type)       :: input_file
    integer(kind=i_def), allocatable :: varid_list(:)
    integer(kind=i_def)              :: i, nvar, ierr
    character(len=nf90_max_name)     :: var_name

    ! Create a list of all fields contained in the input file
    input_file = lfric_ncdf_file_type(trim(file_name) // '.nc',             &
                             open_mode=FILE_OP_OPEN, io_mode=FILE_MODE_READ )

    ! Get the number of fields held in the file for looping over
    ierr = nf90_inquire(input_file%get_id(), nvariables=nvar)

    ! For referencing the actual fields in the file
    allocate(varid_list(nvar))
    varid_list = input_file%get_all_varids()

    ! Allocate the config_list array length to the number of variables
    if (allocated(config_list)) deallocate(config_list)
    allocate(config_list(nvar))

    ! Loop over the number of fields and extract field name from the input file
    num_fields = 0
    do i = 1, nvar

        ! Read from file set to 'ierr' because it is an unused variable
        ierr = nf90_inquire_variable(input_file%get_id(), &
                                     varid_list(i),       &
                                     name=var_name        )

        if (field_is_valid(trim(prefix)//trim(var_name))) then
            num_fields = num_fields + 1                         ! Increment counter
            config_list(num_fields) = trim(var_name)            ! Add field to the config_list

            ! Log the field name to the console
            write(log_scratch_space, '(A)') &
              "Found variable: " // trim(var_name)
            call log_event(log_scratch_space, log_level_info)
        else
          write(log_scratch_space, '(A)') &
            "Field not found in iodef file, skipping: "//&
            trim(prefix)//trim(var_name)
          call log_event(log_scratch_space, log_level_trace)
        endif
    end do

    deallocate(varid_list)
    call input_file%close_file()

    ! Check whether we got 0 fields in the file and throw exception
    if (num_fields == 0) then
      log_scratch_space = 'The input file ' // trim(file_name) // &
                          ' does not contain any valid fields'
      call log_event(log_scratch_space, log_level_error)
    end if

  end subroutine get_field_list

  !> @brief   Creates a field in a given field collection.
  !> @details Calls lfric_xios_space_from_metadata_mod to create a vector space
  !!          and uses this to create a field in the given field collection
  !!          with the given name. The desired 3D and 2D meshes to be
  !!          initialised on must be given but space_from_metadata will
  !!          determine which to use.
  !> @param [in,out]  field_collection  The field collection to which the field
  !!                                    will be added
  !> @param [in]      field_name        The name the field will be given
  !> @param [in]      mesh              The 3D mesh the field could be
  !!                                    initalised on
  !> @param [in]      twod_mesh         The 2D mesh the field could be
  !!                                    intialised on
  subroutine field_maker(field_collection, field_name, mesh, twod_mesh, prefix)

    type(field_collection_type), pointer, intent(inout) :: field_collection
    character(len=*),                     intent(in)    :: field_name
    type(mesh_type),             pointer, intent(in)    :: mesh
    type(mesh_type),             pointer, intent(in)    :: twod_mesh
    character(len=nf90_max_name),         intent(in)    :: prefix

    ! Field object to initialise
    type(field_type) :: field

    ! For field creation
    type( function_space_type ),           pointer :: vector_space
    procedure(write_interface),            pointer :: write_behaviour
    procedure(read_interface),             pointer :: read_behaviour
    procedure(checkpoint_write_interface), pointer :: cp_write_behaviour
    procedure(checkpoint_read_interface),  pointer :: cp_read_behaviour

    ! Set write behaviour to generic
    write_behaviour    => write_field_generic
    read_behaviour     => read_field_generic
    cp_read_behaviour  => checkpoint_read_xios
    cp_write_behaviour => checkpoint_write_xios

    if ( .NOT. field_collection%field_exists(field_name) ) then
      ! Get function space from metadata
      vector_space => space_from_metadata(trim(prefix)//trim(field_name), &
                                          'Regridding',                   &
                                          mesh_3d=mesh,                   &
                                          mesh_2d=twod_mesh               )

      ! Initialise the field
      call field%initialise(vector_space=vector_space, &
                            name=field_name)

      ! Set the generic write behaviours
      call field%set_read_behaviour(read_behaviour)
      call field%set_write_behaviour(write_behaviour)
      call field%set_checkpoint_read_behaviour(cp_read_behaviour)
      call field%set_checkpoint_write_behaviour(cp_write_behaviour)

      call log_event("Adding " // trim(field_name) // &
                     " to field collection " //       &
                    field_collection%get_name(),      &
                    LOG_LEVEL_INFO)

      ! Add to field_colleciton
      call field_collection%add_field(field)
    endif

  end subroutine field_maker

end module lfric2lfric_field_init_mod
