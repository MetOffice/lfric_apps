! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!> @brief   Module containing lfric2um configuration type
!> @details Hold information on input and output files and fields to be
!!          regridded, including proceedures to read this information from
!!          the relevant namelist.
module lfric2um_namelists_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int64, real64

! lfricinputs modules
use lfricinp_um_parameters_mod,     only: fnamelen, um_imdi

implicit none
private

public :: lfric2um_config, required_lfric_namelists

type :: config
  character(len=fnamelen) :: output_filename = 'unset'
  character(len=fnamelen) :: target_grid_namelist = 'unset'
  character(len=fnamelen) :: stashmaster_file = 'unset'
  character(len=fnamelen) :: weights_file_face_centre_to_p_bilinear = 'unset'
  character(len=fnamelen) :: weights_file_face_centre_to_p_neareststod = 'unset'
  character(len=fnamelen) :: weights_file_face_centre_to_u_bilinear = 'unset'
  character(len=fnamelen) :: weights_file_face_centre_to_v_bilinear = 'unset'
  integer(kind=int64), allocatable :: stash_list(:)
  integer(kind=int64), allocatable :: lbtim_list(:)
  integer(kind=int64), allocatable :: lbproc_list(:)
  integer(kind=int64) :: um_version_int = um_imdi
  integer(kind=int64) :: dump_validity_time(6) = um_imdi
  integer :: num_fields

  integer :: status = -1
  character(len=512) :: message = 'No namelist read'
  integer :: unit_number
contains

  procedure :: load_namelists

end type config

! Input namelist configuration
type(config) :: lfric2um_config

! Namelist filenames read from command line
character(len=fnamelen), public :: lfric2um_nl_fname

integer(kind=int64), parameter :: max_stash_list = 999

character(*), parameter  :: required_lfric_namelists(6) = ['logging         ', &
                                                           'finite_element  ', &
                                                           'base_mesh       ', &
                                                           'planet          ', &
                                                           'extrusion       ', &
                                                           'io              ']

contains

subroutine load_namelists(self)

! Descriptions:
!  Reads in lfric2um namelists. Performs checking on namelist values.
!  Populates UM grid object and prints out grid diagnostics.

! LFRic modules
use log_mod,                    only: log_event,       &
                                      LOG_LEVEL_ERROR, &
                                      LOG_LEVEL_INFO
use constants_mod,              only: imdi, rmdi

! lfricinputs modules
use lfricinp_grid_namelist_mod, only: grid,                                &
                                      lambda_origin_targ, phi_origin_targ, &
                                      phi_pole, lambda_pole,               &
                                      delta_lambda_targ, delta_phi_targ,   &
                                      points_lambda_targ, points_phi_targ, &
                                      igrid_targ, rotated
use lfricinp_unit_handler_mod,  only: get_free_unit

! lfric2um modules
use lfricinp_um_grid_mod,       only: um_grid

implicit none
class(config) :: self


! Local variables
integer :: i_stash

! Namelist variables
character(len=fnamelen) :: output_filename = 'unset'
character(len=fnamelen) :: stashmaster_file = 'unset'
character(len=fnamelen) :: target_grid_namelist = 'unset'
character(len=fnamelen) :: weights_file_face_centre_to_p_bilinear = 'unset'
character(len=fnamelen) :: weights_file_face_centre_to_p_neareststod = 'unset'
character(len=fnamelen) :: weights_file_face_centre_to_u_bilinear = 'unset'
character(len=fnamelen) :: weights_file_face_centre_to_v_bilinear = 'unset'
integer(kind=int64) :: stash_list(max_stash_list)
integer(kind=int64) :: lbtim_list(max_stash_list)
integer(kind=int64) :: lbproc_list(max_stash_list)
integer(kind=int64) :: um_version_int = um_imdi
integer(kind=int64) :: dump_validity_time(6) = um_imdi

namelist /configure_lfric2um/ output_filename,                                 &
                              target_grid_namelist,                            &
                              stashmaster_file,                                &
                              weights_file_face_centre_to_p_bilinear,          &
                              weights_file_face_centre_to_p_neareststod,       &
                              weights_file_face_centre_to_u_bilinear,          &
                              weights_file_face_centre_to_v_bilinear,          &
                              stash_list,                                      &
                              lbtim_list,                                      &
                              lbproc_list,                                     &
                              um_version_int,                                  &
                              dump_validity_time

stash_list(:) = um_imdi
lbtim_list(:) = um_imdi
lbproc_list(:) = um_imdi

self%status = 0
self%message = 'Reading namelist from ' // trim(lfric2um_nl_fname)

call get_free_unit(self%unit_number)

open(unit=self%unit_number, file=lfric2um_nl_fname, iostat=self%status,                    &
                            iomsg=self%message)
if (self%status /= 0) call log_event(self%message, LOG_LEVEL_ERROR)

read(self%unit_number, nml=configure_lfric2um, iostat=self%status,             &
                       iomsg=self%message)
if (self%status /= 0) call log_event(self%message, LOG_LEVEL_ERROR)

if (trim(output_filename) == 'unset') then
  self%status = 1
  self%message='Target filename is unset'
  call log_event(self%message, LOG_LEVEL_ERROR)
end if

if (trim(stashmaster_file) == 'unset') then
  self%status = 1
  self%message='Stashmaster filename is unset'
  call log_event(self%message, LOG_LEVEL_ERROR)
end if

if (trim(weights_file_face_centre_to_p_bilinear) == 'unset') then
  self%status = 1
  self%message='weights_file_face_centre_to_p_bilinear filename is unset'
  call log_event(self%message, LOG_LEVEL_ERROR)
end if

if (trim(weights_file_face_centre_to_p_neareststod) == 'unset') then
  self%status = 1
  self%message='weights_file_face_centre_to_p_neareststod filename is unset'
  call log_event(self%message, LOG_LEVEL_ERROR)
end if

if (trim(weights_file_face_centre_to_u_bilinear) == 'unset') then
  self%status = 1
  self%message='weights_file_face_centre_to_u_bilinear filename is unset'
  call log_event(self%message, LOG_LEVEL_ERROR)
end if

if (trim(weights_file_face_centre_to_v_bilinear) == 'unset') then
  self%status = 1
  self%message='weights_file_face_centre_to_v_bilinear filename is unset'
  call log_event(self%message, LOG_LEVEL_ERROR)
end if

if (trim(target_grid_namelist) == 'unset') then
  self%status = 1
  self%message='Target grid namelist is unset'
  call log_event(self%message, LOG_LEVEL_ERROR)
end if

close(self%unit_number)

! Now read grid namelist to define UM grid
open(unit=self%unit_number, file=target_grid_namelist, iostat=self%status,     &
                            iomsg=self%message)
if (self%status /= 0) call log_event(self%message, LOG_LEVEL_ERROR)
read(self%unit_number, nml=grid, iostat=self%status,                           &
                       iomsg=self%message)
if (self%status /= 0) call log_event(self%message, LOG_LEVEL_ERROR)

! Load namelist variables into objects
self%output_filename = output_filename
self%target_grid_namelist = target_grid_namelist
self%stashmaster_file = stashmaster_file
self%weights_file_face_centre_to_p_bilinear =                                  &
                                 weights_file_face_centre_to_p_bilinear
self%weights_file_face_centre_to_p_neareststod =                               &
                                 weights_file_face_centre_to_p_neareststod
self%weights_file_face_centre_to_u_bilinear =                                  &
                                 weights_file_face_centre_to_u_bilinear
self%weights_file_face_centre_to_v_bilinear =                                  &
                                 weights_file_face_centre_to_v_bilinear
self%um_version_int = um_version_int
self%dump_validity_time(:) = dump_validity_time(:)

self%num_fields=0
! Count how many fields have been requested
do i_stash = 1, max_stash_list
  if (stash_list(i_stash) == um_imdi) then
    exit
  else
    self%num_fields = self%num_fields + 1
  end if
end do

if (self%num_fields <= 0) then
  call log_event('No fields selected in stash_list namelist variable',         &
       LOG_LEVEL_ERROR)
end if

! Can now allocate type variable
allocate(self%stash_list(self%num_fields))
allocate(self%lbtim_list(self%num_fields))
allocate(self%lbproc_list(self%num_fields))
self%stash_list(:) = stash_list(1:self%num_fields)
self%lbtim_list(:) = lbtim_list(1:self%num_fields)
self%lbproc_list(:) = lbproc_list(1:self%num_fields)

call um_grid%set_grid_coords(                                                  &
     grid_staggering = igrid_targ ,                                            &
     num_p_points_x = points_lambda_targ,                                      &
     num_p_points_y = points_phi_targ,                                         &
     grid_spacing_x = delta_lambda_targ,                                       &
     grid_spacing_y = delta_phi_targ,                                          &
     ! namelist provides p grid values and routine
     ! wants base grid origin, so apply offset
     grid_origin_x = lambda_origin_targ - (0.5_real64 * delta_lambda_targ),           &
     grid_origin_y = phi_origin_targ - (0.5_real64 * delta_phi_targ))

call um_grid%print_grid_coords()

self%status = 0
self%message = 'Successfully read namelists from ' // trim(lfric2um_nl_fname)
call log_event(self%message, LOG_LEVEL_INFO)
close(self%unit_number)

end subroutine load_namelists

end module lfric2um_namelists_mod
