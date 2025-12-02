! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_regrid_options_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int64

! UM2LFRic modules
use lfricinp_um_parameters_mod,     only: fnamelen, um_imdi

implicit none
private
public :: lfricinp_init_regrid_options, interp_method, winds_on_w3, &
          regrid_type, specify_nearest_neighbour, nn_fields

character(len=fnamelen) :: interp_method = 'bilinear'
character(len=fnamelen) :: regrid_type = 'global_to_global'
logical :: winds_on_w3 = .true.

integer(kind=int64), parameter :: max_stash_list = 999
integer(kind=int64)     :: specify_nearest_neighbour(max_stash_list)
integer :: nn_fields ! Number of nearest neighbour interpolation fields

contains

subroutine lfricinp_init_regrid_options(fname)

use lfricinp_unit_handler_mod, only: get_free_unit
use log_mod,                   only: log_event, log_scratch_space,             &
                                     LOG_LEVEL_ERROR

implicit none

character(len=fnamelen) :: fname
integer                 :: status = -1
character(len=512)      :: message = 'No namelist read'
integer                 :: unit_number
integer                 :: i_stash

namelist /regrid_options/ interp_method, regrid_type, winds_on_w3, &
     specify_nearest_neighbour

specify_nearest_neighbour(:) = um_imdi

call get_free_unit(unit_number)

open(unit=unit_number, file=fname, iostat=status, iomsg=message)
if (status /= 0) call log_event(message, LOG_LEVEL_ERROR)

read(unit_number, nml=regrid_options, iostat=status, iomsg=message)
if (status /= 0) call log_event(message, LOG_LEVEL_ERROR)

if ( (trim(interp_method) /= 'bilinear') .and.                                 &
     (trim(interp_method) /= 'copy') ) then
  write(log_scratch_space, '(A)') 'Interpolation method: '                     &
                                  // trim(interp_method) // ' not recognised.'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

if ( (trim(regrid_type) /= 'global_to_global') .and.                           &
     (trim(regrid_type) /= 'lam_to_lam') .and.                                 &
     (trim(regrid_type) /= 'lbc_to_lbc') ) then
  write(log_scratch_space, '(A)') 'Regridding type: ' // trim(regrid_type) //  &
                                  ' not recognised.'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

if ( (trim(interp_method) == 'copy') .and.                                     &
     (trim(regrid_type) /= 'lam_to_lam' .and.                                  &
      trim(regrid_type) /= 'lbc_to_lbc') ) then
  write(log_scratch_space, '(A)') 'The copy interpolation method is only ' //  &
                                  'valid for LAM to LAM or LBC to LBC ' //     &
                                  'regridding'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

if ( (trim(interp_method) == 'bilinear') .and.                                 &
     (.not. winds_on_w3) ) then
  write(log_scratch_space, '(A)') 'For bilinear interpolation method winds ' //&
                                  'must be placed on W3'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

if ( (trim(interp_method) == 'copy') .and. winds_on_w3 ) then
  write(log_scratch_space, '(A)')                                              &
                            'For copying method winds must not be placed on W3'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

nn_fields=0
! Count how many fields in specified lists
do i_stash = 1, max_stash_list
  if (specify_nearest_neighbour(i_stash) == um_imdi) then
    exit
  else
    nn_fields = nn_fields + 1
  end if
end do

close(unit_number)

end subroutine lfricinp_init_regrid_options

end module lfricinp_regrid_options_mod
