! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!> @brief     Module containing io_config TYPE
!> @details   Holds details from the io namelist including input and output
!!            file details and any ancillary files required

MODULE lfricinp_setup_io_mod

USE constants_mod,                 ONLY: str_max_filename
USE file_mod,                      ONLY: FILE_MODE_READ,                       &
                                         FILE_MODE_WRITE
USE log_mod,                       ONLY: log_event, LOG_LEVEL_ERROR
USE lfric_xios_file_mod,           ONLY: lfric_xios_file_type,                 &
                                         OPERATION_TIMESERIES
USE linked_list_mod,               ONLY: linked_list_type
USE lfricinp_um_parameters_mod,    ONLY: fnamelen
IMPLICIT NONE
PRIVATE
PUBLIC :: io_config, io_fname

INTEGER, PARAMETER              :: max_number_ancfiles = 20

! Type to contain the details from the io namelist, along with namelist reading
! and initialisation procedures
TYPE :: config
  CHARACTER(LEN=str_max_filename) :: checkpoint_read_file  = 'unset'
  CHARACTER(LEN=str_max_filename) :: checkpoint_write_file = 'unset'
  CHARACTER(LEN=str_max_filename) :: ancil_file_map(max_number_ancfiles) = 'unset'
  LOGICAL :: checkpoint_write, checkpoint_read, ancil_read

  CONTAINS

  PROCEDURE :: load_namelist
  PROCEDURE :: init_lfricinp_files

END TYPE config

TYPE(config) :: io_config
CHARACTER(LEN=fnamelen) :: io_fname

CONTAINS

!> @brief   Reads details from the io namelist and store the output in io_config
SUBROUTINE load_namelist(self)

USE lfricinp_unit_handler_mod, ONLY: get_free_unit

IMPLICIT NONE

CLASS(config), INTENT(IN OUT) :: self

CHARACTER(LEN=512) :: message = 'No namelist read'
INTEGER            :: unit_number
INTEGER            :: status = -1

! Initialise local variables to read the namelist into
CHARACTER(LEN=str_max_filename) :: checkpoint_read_file  = 'unset'
CHARACTER(LEN=str_max_filename) :: checkpoint_write_file = 'unset'
CHARACTER(LEN=str_max_filename) :: ancil_file_map(max_number_ancfiles) = 'unset'
LOGICAL :: checkpoint_write, checkpoint_read, ancil_read

NAMELIST /iofiles/ checkpoint_read,                                            &
                   checkpoint_write,                                           &
                   ancil_read,                                                 &
                   checkpoint_read_file,                                       &
                   checkpoint_write_file,                                      &
                   ancil_file_map

! Read the namelist
CALL get_free_unit(unit_number)

OPEN(UNIT=unit_number, FILE=io_fname, IOSTAT=status, IOMSG=message)
IF (status /= 0) CALL log_event(message, LOG_LEVEL_ERROR)

READ(unit_number, NML=iofiles, IOSTAT=status, IOMSG=message)
IF (status /= 0) CALL log_event(message, LOG_LEVEL_ERROR)

! Load namelist variables into object
self%checkpoint_read = checkpoint_read
self%checkpoint_write = checkpoint_write
self%ancil_read = ancil_read
self%checkpoint_read_file = checkpoint_read_file
self%checkpoint_write_file = checkpoint_write_file
self%ancil_file_map = ancil_file_map

CLOSE(unit_number)

END SUBROUTINE load_namelist

!> @brief Takes the namelist variables and sets up a list of files to be used.
SUBROUTINE init_lfricinp_files(self, files_list)

IMPLICIT NONE

CLASS(config),                            INTENT(IN OUT) :: self
TYPE(linked_list_type),                   INTENT(OUT)    :: files_list

INTEGER, PARAMETER                     :: checkpoint_frequency = 1
CHARACTER(LEN=str_max_filename)        :: ancil_xios_file_id,                  &
                                          ancil_file_path,                     &
                                          afm
INTEGER                                :: split_idx, i

IF (self%ancil_read) THEN
  ! Set ancil file reading context information for all required ancil files
  DO i = 1, max_number_ancfiles
    ! Exit loop if entry in acil file map is unset
    afm = self%ancil_file_map(i)
    IF( TRIM(afm) == 'unset') EXIT
    ! From ancil file map string extract ancil xios file id and ancil file path
    split_idx = INDEX(afm, ':')
    ancil_xios_file_id = afm(1:split_idx-1)
    ancil_file_path = afm(split_idx+1:)
    ! Initial ancil file and insert file in file list
    CALL files_list%insert_item( lfric_xios_file_type(                          &
                                              TRIM(ancil_file_path),            &
                                              xios_id=TRIM(ancil_xios_file_id), &
                                              io_mode=FILE_MODE_READ) )
  END DO
END IF

! Setup checkpoint writing context information
IF (self%checkpoint_write) THEN
  ! Create checkpoint filename from stem and first timestep
  CALL files_list%insert_item( lfric_xios_file_type(                            &
                                              TRIM(self%checkpoint_write_file), &
                                              xios_id="lfric_checkpoint_write", &
                                              io_mode=FILE_MODE_WRITE,          &
                         ! For some reason LI outputs checkpoints as a timeseries
                                              operation=OPERATION_TIMESERIES,   &
                                              freq=checkpoint_frequency ) )
END IF

! Setup checkpoint reading context information
IF (self%checkpoint_read) THEN
  ! Create checkpoint filename from stem
  CALL files_list%insert_item( lfric_xios_file_type(                            &
                                               TRIM(self%checkpoint_read_file), &
                                               xios_id="lfric_checkpoint_read", &
                                               io_mode=FILE_MODE_READ ) )
END IF

END SUBROUTINE init_lfricinp_files

END MODULE lfricinp_setup_io_mod
