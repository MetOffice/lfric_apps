! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module um2lfric_post_process_fields_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : real64, int64

implicit none

private

public :: um2lfric_post_process_fields

contains

subroutine um2lfric_post_process_fields(field, stashcode)
! Description: Applies transformations to fields. Should be applied
! post regridding but before the fields are transferred into the
! LFRic field types.

! lfric modules
use log_mod,                         only: log_event,      &
                                           LOG_LEVEL_INFO, &
                                           log_scratch_space
! um2lfric modules
use lfricinp_um_grid_mod,            only: um_grid
! lfricinputs modules
use lfricinp_stashmaster_mod,        only: stashcode_area_cf,    &
                                           stashcode_v,          &
                                           get_stashmaster_item, &
                                           pseudL,               &
                                           snow_layers_and_tiles
use lfricinp_reorder_snow_field_mod, only: lfricinp_reorder_snow_field
use lfricinp_add_bottom_level_mod,   only: lfricinp_add_bottom_level
use lfricinp_regrid_options_mod,     only: winds_on_w3
use lfricinp_um_parameters_mod,      only: um_rmdi

implicit none

! Arguments
real(kind=real64), allocatable, intent(in out) :: field(:,:)
integer(kind=int64), intent(in) :: stashcode

! Any post procesing steps determined by STASHcode

select case (stashcode)
case (stashcode_area_cf)
  write(log_scratch_space, '(A,I0)') "Adding bottom level for stashcode: ",    &
       stashcode
  call log_event(log_scratch_space, LOG_LEVEL_INFO)
  call lfricinp_add_bottom_level(field)

case (stashcode_v)
  if (.not. winds_on_w3) then
    write(log_scratch_space,'(A)') "Change sign of 'V' component of " //       &
                                   "wind on W2H function space"
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
    ! Switch sign of data, but leave RMDI values unchanged
    field(:,:) = field(:,:) * merge(-1.0_real64, 1.0_real64, field /= um_rmdi)
  end if
end select

! Post processing steps determined by pseudo level code

select case (get_stashmaster_item(stashcode, pseudL))
case(snow_layers_and_tiles)
  write(log_scratch_space, '(A,I0)') "Reordering snow layers for stashcode: ", &
       stashcode
  call log_event(log_scratch_space, LOG_LEVEL_INFO)
  call lfricinp_reorder_snow_field(field, um_grid)
end select

end subroutine um2lfric_post_process_fields

end module um2lfric_post_process_fields_mod
