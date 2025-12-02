! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!> @brief     Module containing io_config type
!> @details   Holds details from the io namelist including input and output
!!            file details and any ancillary files required

module lfricinp_setup_io_mod

use constants_mod,                 only: str_max_filename
use file_mod,                      only: FILE_MODE_READ,                       &
                                         FILE_MODE_WRITE
use log_mod,                       only: log_event, LOG_LEVEL_ERROR
use lfric_xios_file_mod,           only: lfric_xios_file_type,                 &
                                         OPERATION_TIMESERIES
use linked_list_mod,               only: linked_list_type
use lfricinp_um_parameters_mod,    only: fnamelen
implicit none
private
public :: io_config, io_fname

integer, parameter              :: max_number_ancfiles = 20

! Type to contain the details from the io namelist, along with namelist reading
! and initialisation procedures
type :: config
  character(len=str_max_filename) :: checkpoint_read_file  = 'unset'
  character(len=str_max_filename) :: checkpoint_write_file = 'unset'
  character(len=str_max_filename) :: ancil_file_map(max_number_ancfiles) = 'unset'
  logical :: checkpoint_write, checkpoint_read, ancil_read

  contains

  procedure :: load_namelist
  procedure :: init_lfricinp_files

end type config

type(config) :: io_config
character(len=fnamelen) :: io_fname

contains

!> @brief   Reads details from the io namelist and store the output in io_config
subroutine load_namelist(self)

use lfricinp_unit_handler_mod, only: get_free_unit

implicit none

class(config), intent(in out) :: self

character(len=512) :: message = 'No namelist read'
integer            :: unit_number
integer            :: status = -1

! Initialise local variables to read the namelist into
character(len=str_max_filename) :: checkpoint_read_file  = 'unset'
character(len=str_max_filename) :: checkpoint_write_file = 'unset'
character(len=str_max_filename) :: ancil_file_map(max_number_ancfiles) = 'unset'
logical :: checkpoint_write, checkpoint_read, ancil_read

namelist /iofiles/ checkpoint_read,                                            &
                   checkpoint_write,                                           &
                   ancil_read,                                                 &
                   checkpoint_read_file,                                       &
                   checkpoint_write_file,                                      &
                   ancil_file_map

! Read the namelist
call get_free_unit(unit_number)

open(unit=unit_number, file=io_fname, iostat=status, iomsg=message)
if (status /= 0) call log_event(message, LOG_LEVEL_ERROR)

read(unit_number, nml=iofiles, iostat=status, iomsg=message)
if (status /= 0) call log_event(message, LOG_LEVEL_ERROR)

! Load namelist variables into object
self%checkpoint_read = checkpoint_read
self%checkpoint_write = checkpoint_write
self%ancil_read = ancil_read
self%checkpoint_read_file = checkpoint_read_file
self%checkpoint_write_file = checkpoint_write_file
self%ancil_file_map = ancil_file_map

close(unit_number)

end subroutine load_namelist

!> @brief Takes the namelist variables and sets up a list of files to be used.
subroutine init_lfricinp_files(self, files_list)

implicit none

class(config),                            intent(in out) :: self
type(linked_list_type),                   intent(out)    :: files_list

integer, parameter                     :: checkpoint_frequency = 1
character(len=str_max_filename)        :: ancil_xios_file_id,                  &
                                          ancil_file_path,                     &
                                          afm
integer                                :: split_idx, i

if (self%ancil_read) then
  ! Set ancil file reading context information for all required ancil files
  do i = 1, max_number_ancfiles
    ! Exit loop if entry in acil file map is unset
    afm = self%ancil_file_map(i)
    if( trim(afm) == 'unset') exit
    ! From ancil file map string extract ancil xios file id and ancil file path
    split_idx = index(afm, ':')
    ancil_xios_file_id = afm(1:split_idx-1)
    ancil_file_path = afm(split_idx+1:)
    ! Initial ancil file and insert file in file list
    call files_list%insert_item( lfric_xios_file_type(                          &
                                              trim(ancil_file_path),            &
                                              xios_id=trim(ancil_xios_file_id), &
                                              io_mode=FILE_MODE_READ) )
  end do
end if

! Setup checkpoint writing context information
if (self%checkpoint_write) then
  ! Create checkpoint filename from stem and first timestep
  call files_list%insert_item( lfric_xios_file_type(                            &
                                              trim(self%checkpoint_write_file), &
                                              xios_id="lfric_checkpoint_write", &
                                              io_mode=FILE_MODE_WRITE,          &
                         ! For some reason LI outputs checkpoints as a timeseries
                                              operation=OPERATION_TIMESERIES,   &
                                              freq=checkpoint_frequency ) )
end if

! Setup checkpoint reading context information
if (self%checkpoint_read) then
  ! Create checkpoint filename from stem
  call files_list%insert_item( lfric_xios_file_type(                            &
                                               trim(self%checkpoint_read_file), &
                                               xios_id="lfric_checkpoint_read", &
                                               io_mode=FILE_MODE_READ ) )
end if

end subroutine init_lfricinp_files

end module lfricinp_setup_io_mod
