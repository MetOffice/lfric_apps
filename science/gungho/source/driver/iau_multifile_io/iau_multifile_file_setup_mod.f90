!-------------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @brief Sets up files for use by multifile IO for the IAU
!>
module iau_multifile_file_setup_mod

  use constants_mod,         only: i_def, &
                                   str_def
  use driver_modeldb_mod,    only : modeldb_type
  use field_collection_mod,  only: field_collection_type
  use file_mod,              only: file_type, FILE_MODE_READ
  use lfric_xios_file_mod,   only: lfric_xios_file_type, OPERATION_ONCE
  use linked_list_mod,       only: linked_list_type

  ! Configuration modules
  use io_config_mod,       only: use_xios_io
  use time_config_mod,     only: timestep_start, &
                                 timestep_end

  implicit none

  private
  public :: init_iau_inc_files

contains

!> @brief Sets up increment files for multifile IO
!> @param[out]   files_list Structure holding file information
!> @param[inout] modeldb    Structure holding model data
!> @param[in]    iau_incs   Type of IAU increment
!> @param[in]    filename   Name of file to be added to files_list
!>
  subroutine init_iau_inc_files(files_list, modeldb, &
                                iau_incs, filename)

    implicit none

    type(linked_list_type), intent(out)    :: files_list
    type(modeldb_type),     intent(inout)  :: modeldb
    character(*),           intent(in)     :: iau_incs
    character(str_def),     intent(in)     :: filename

    integer(i_def) :: ts_start, ts_end
    integer(i_def) :: rc

    type(field_collection_type), pointer :: multifile_fields

    multifile_fields  => modeldb%fields%get_field_collection(iau_incs)

    ! Get time configuration in integer form
    read(timestep_start,*,iostat=rc)  ts_start
    read(timestep_end,*,iostat=rc)  ts_end

    if ( use_xios_io) then

      call files_list%insert_item( &
        lfric_xios_file_type( filename, &
        xios_id=iau_incs, &
        io_mode=FILE_MODE_READ, &
        freq=1, &
        operation=OPERATION_ONCE, &
        fields_in_file=multifile_fields))

    end if ! use_xios_io

  end subroutine init_iau_inc_files

end module iau_multifile_file_setup_mod
