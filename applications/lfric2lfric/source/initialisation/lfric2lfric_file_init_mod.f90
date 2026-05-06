!-------------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief   Sets up I/O configuration from within lfric2lfric.
!> @details Collects configuration information relevant for the I/O subsystem
!!          and formats it so that it can be passed to the infrastructure.
module lfric2lfric_file_init_mod

  use constants_mod,          only: i_def,                &
                                    str_max_filename
  use driver_modeldb_mod,     only: modeldb_type
  use file_mod,               only: FILE_MODE_READ,       &
                                    FILE_MODE_WRITE
  use io_config_mod,          only: diagnostic_frequency, &
                                    checkpoint_write,     &
                                    checkpoint_read,      &
                                    write_diag,           &
                                    use_xios_io
  use lfric_xios_file_mod,    only: lfric_xios_file_type, &
                                    OPERATION_TIMESERIES
  use lfric2lfric_config_mod, only: dst_ancil_directory,           &
                                    dst_orography_mean_ancil_path, &
                                    src_ancil_directory,           &
                                    src_orography_mean_ancil_path, &
                                    mode_ics, mode_lbc
  use linked_list_mod,        only: linked_list_type
  use namelist_mod,           only: namelist_type
  use orography_config_mod,   only: orog_init_option,    &
                                    orog_init_option_ancil

  implicit none

  private
  public :: init_lfric2lfric_src_files, init_lfric2lfric_dst_files

  contains

  !> @brief   Sets up source I/O configuration.
  !> @details Initialises the file list for the source I/O context, using the
  !!          start_dump_filename extracted from the `files` namelist, or the
  !!          source_file_lbc extracted from the `lfric2lfric` namelist.
  !> @param [out]        files_list    The list of I/O files.
  !> @param [in,out]     modeldb       Required by init_io.
  subroutine init_lfric2lfric_src_files( files_list, modeldb )

    implicit none

    type(linked_list_type),        intent(out)   :: files_list
    type(modeldb_type), optional,  intent(inout) :: modeldb

    type(namelist_type),    pointer :: lfric2lfric_nml
    type(namelist_type),    pointer :: files_nml

    integer(kind=i_def)             :: mode
    character(len=str_max_filename) :: start_dump_filename
    character(len=str_max_filename) :: source_file_lbc
    character(len=str_max_filename) :: src_ancil_fname

    if( use_xios_io ) then

      lfric2lfric_nml    => modeldb%configuration%get_namelist('lfric2lfric')
      call lfric2lfric_nml%get_value( 'mode', mode )

      files_nml    => modeldb%configuration%get_namelist('files')
      if (mode == mode_ics) then
        call files_nml%get_value( 'start_dump_filename', start_dump_filename )

        ! Setup checkpoint reading context information
        call files_list%insert_item(                               &
            lfric_xios_file_type( trim(start_dump_filename),       &
                                  xios_id="lfric_checkpoint_read", &
                                  io_mode=FILE_MODE_READ )         )

      else if (mode == mode_lbc) then
        call lfric2lfric_nml%get_value( 'source_file_lbc', source_file_lbc )

        ! Setup lbc source reading context information
        call files_list%insert_item(                               &
            lfric_xios_file_type( trim(source_file_lbc),           &
                                  xios_id="lfric_lbc_read",        &
                                  io_mode=FILE_MODE_READ,          &
                                  operation=OPERATION_TIMESERIES,  &
                                  freq=diagnostic_frequency) )
      endif

      ! Setup orography ancillary file
      if ( orog_init_option == orog_init_option_ancil ) then
        ! Set orography ancil filename from namelist
        write(src_ancil_fname,'(A)') trim(src_ancil_directory)//'/'// &
                                     trim(src_orography_mean_ancil_path)
        call files_list%insert_item( lfric_xios_file_type(                     &
                                           trim(src_ancil_fname),              &
                                           xios_id="src_orography_mean_ancil", &
                                           io_mode=FILE_MODE_READ ) )
      end if

    endif

  end subroutine init_lfric2lfric_src_files

  !> @brief   Sets up destination I/O configuration.
  !> @details Initialises the file list for the destination I/O context, using
  !!          the checkpoint_stem_name extracted from the `files` namelist.
  !> @param [out]        files_list    The list of I/O files.
  !> @param [in,out]     modeldb       Required by init_io.
  subroutine init_lfric2lfric_dst_files( files_list, modeldb )

    implicit none

    type(linked_list_type),        intent(out)   :: files_list
    type(modeldb_type), optional,  intent(inout) :: modeldb

    ! Local variables
    type(namelist_type),    pointer :: lfric2lfric_nml
    type(namelist_type),    pointer :: files_nml
    integer(kind=i_def),  parameter :: checkpoint_frequency = 1_i_def

    integer(kind=i_def)             :: mode
    character(len=str_max_filename) :: dst_ancil_fname
    character(len=str_max_filename) :: checkpoint_stem_name
    character(len=str_max_filename) :: diag_stem_name

    if( use_xios_io ) then

      lfric2lfric_nml    => modeldb%configuration%get_namelist('lfric2lfric')
      call lfric2lfric_nml%get_value( 'mode', mode )

      files_nml    => modeldb%configuration%get_namelist('files')
      ! Set up diagnostic writing info
      if( write_diag ) then
        call files_nml%get_value( 'diag_stem_name', diag_stem_name )

        ! Setup diagnostic output file
        call files_list%insert_item(                         &
            lfric_xios_file_type( trim( diag_stem_name ),    &
                                  xios_id="lfric_diag",      &
                                  io_mode=FILE_MODE_WRITE,   &
                                  freq=diagnostic_frequency) )
      endif

      if (mode == mode_ics) then
        ! Setup checkpoint writing context information
        if ( checkpoint_write ) then
          call files_nml%get_value( 'checkpoint_stem_name', checkpoint_stem_name )

          call files_list%insert_item(                                &
              lfric_xios_file_type( trim( checkpoint_stem_name ),     &
                                    xios_id="lfric_checkpoint_write", &
                                    io_mode=FILE_MODE_WRITE,          &
                                    freq=checkpoint_frequency )       )
        end if
      else if (mode == mode_lbc) then
        ! Setup lbc writing context information
        call files_list%insert_item(                                &
            lfric_xios_file_type( "lfric2lfric_lbc",                &
                                  xios_id="lfric_lbc_write",        &
                                  io_mode=FILE_MODE_WRITE,          &
                                  operation=OPERATION_TIMESERIES,   &
                                  freq=diagnostic_frequency ) )
      endif

      ! Setup orography ancillary file
      if ( orog_init_option == orog_init_option_ancil ) then
        ! Set orography ancil filename from namelist
        write(dst_ancil_fname,'(A)') trim(dst_ancil_directory)//'/'// &
                                     trim(dst_orography_mean_ancil_path)
        call files_list%insert_item( lfric_xios_file_type(                     &
                                           trim(dst_ancil_fname),              &
                                           xios_id="dst_orography_mean_ancil", &
                                           io_mode=FILE_MODE_READ ) )
      end if

    endif

  end subroutine init_lfric2lfric_dst_files

end module lfric2lfric_file_init_mod
