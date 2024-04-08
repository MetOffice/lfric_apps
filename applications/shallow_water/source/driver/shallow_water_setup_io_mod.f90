!-------------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief   Sets up I/O configuration from within shallow_water.
!> @details Collects configuration information relevant for the I/O subsystem
!!          and formats it so that it can be passed to the infrastructure.
module shallow_water_setup_io_mod

  use constants_mod,                 only: i_def, &
                                           str_def, str_max_filename
  use driver_modeldb_mod,            only: modeldb_type
  use file_mod,                      only: FILE_MODE_READ, &
                                           FILE_MODE_WRITE
  use files_config_mod,              only: checkpoint_stem_name
  use io_config_mod,                 only: diagnostic_frequency, &
                                           checkpoint_write,     &
                                           checkpoint_read,      &
                                           write_dump
  use lfric_xios_file_mod,           only: lfric_xios_file_type
  use linked_list_mod,               only: linked_list_type, &
                                           linked_list_item_type
  use time_config_mod,               only: timestep_start, &
                                           timestep_end

  implicit none

  private
  public :: init_shallow_water_files

  contains

  !> @brief  Sets up I/O configuration.
  !> @param[out] files_list The list of I/O files
  subroutine init_shallow_water_files(files_list, modeldb)

    implicit none

    type(linked_list_type),        intent(out) :: files_list
    type(modeldb_type), optional, intent(inout)  :: modeldb

    character(len=str_max_filename) :: checkpoint_write_fname, &
                                       checkpoint_read_fname
    integer(i_def)                  :: ts_start, ts_end
    integer(i_def)                  :: rc

    ! Get time configuration in integer form
    read(timestep_start,*,iostat=rc)  ts_start
    read(timestep_end,*,iostat=rc)    ts_end

    ! Setup diagnostic output file
    call files_list%insert_item( lfric_xios_file_type( "lfric_diag",            &
                                                       xios_id="lfric_diag",    &
                                                       io_mode=FILE_MODE_WRITE, &
                                                       freq=diagnostic_frequency) )

    ! Setup checkpoint writing context information
    if ( checkpoint_write ) then
      ! Create checkpoint filename from stem and end timestep
      write(checkpoint_write_fname,'(A,A,I6.6)') &
                           trim(checkpoint_stem_name),"_", ts_end
      call files_list%insert_item( lfric_xios_file_type( checkpoint_write_fname,           &
                                                         xios_id="lfric_checkpoint_write", &
                                                         freq=ts_end,                      &
                                                         io_mode=FILE_MODE_WRITE,          &
                                                         field_group_id="checkpoint_fields") )
    end if

    ! Setup checkpoint reading context information
    if ( checkpoint_read ) then
      ! Create checkpoint filename from stem and (start - 1) timestep
      write(checkpoint_read_fname,'(A,A,I6.6)') &
                   trim(checkpoint_stem_name),"_", (ts_start - 1)
      call files_list%insert_item( lfric_xios_file_type( checkpoint_read_fname,           &
                                                         xios_id="lfric_checkpoint_read", &
                                                         freq=ts_start - 1,               &
                                                         io_mode=FILE_MODE_READ,          &
                                                         field_group_id="checkpoint_fields") )
    end if

  end subroutine init_shallow_water_files

end module shallow_water_setup_io_mod
