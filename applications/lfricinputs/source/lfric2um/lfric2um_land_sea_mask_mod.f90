!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some content in this file was generated or refactored with assistance from
! - Met Office Github Copilot Enterprise (Claude Opus 5.5)
!> @brief   Optional UM land-sea mask functionality
module lfric2um_land_sea_mask_mod

use, intrinsic :: iso_fortran_env, only: int64, real64

implicit none

private

public :: lfric2um_read_land_sea_mask, &
          lfric2um_apply_land_sea_mask

! UM P-grid land-sea mask (land = .true.), only allocated on rank 0
! (Passing the mask as a module-level variable is unpleasant,
! but it is the lfricinputs way)
logical, allocatable :: um_land_sea_mask(:,:)

contains

!> @brief Read the UM land-sea mask (stash 30) from a UM file
!> @param[in] fname UM fieldsfile/ancillary containing the land-sea mask
subroutine lfric2um_read_land_sea_mask(fname)

! lfricinputs modules
use lfricinp_check_shumlib_status_mod, only: shumlib
use lfricinp_stashmaster_mod,          only: stashcode_lsm
use lfricinp_um_grid_mod,              only: um_grid
use lfricinp_um_parameters_mod,        only: fnamelen

! LFRic modules
use log_mod, only: log_event,         &
                   log_scratch_space, &
                   LOG_LEVEL_INFO,    &
                   LOG_LEVEL_ERROR

! Shumlib modules
use f_shum_file_mod,  only: shum_file_type
use f_shum_field_mod, only: shum_field_type

implicit none

character(len=fnamelen), intent(in) :: fname

character(len=*), parameter :: routinename = 'lfric2um_read_land_sea_mask'
type(shum_file_type) :: mask_file
type(shum_field_type), allocatable :: mask_fields(:)

call log_event('Loading UM land-sea mask from '//trim(fname), LOG_LEVEL_INFO)
call shumlib(routinename//'::open_file', mask_file%open_file(fname))
call shumlib(routinename//'::read_header', mask_file%read_header())
call shumlib(routinename//'::find_fields_in_file',      &
             mask_file%find_fields_in_file(mask_fields, &
                                           stashcode=stashcode_lsm))

! Mask may be held as real or integer data depending on the file
! Convert to a logical mask
if (allocated(mask_fields(1)%rdata)) then
  um_land_sea_mask = (mask_fields(1)%rdata /= 0.0_real64)
else if (allocated(mask_fields(1)%idata)) then
  um_land_sea_mask = (mask_fields(1)%idata /= 0_int64)
else
  call log_event('No data found for UM land-sea mask in '//trim(fname), &
                 LOG_LEVEL_ERROR)
end if

call shumlib(routinename//'::close_file', mask_file%close_file())

if (size(um_land_sea_mask, 1) /= um_grid%num_p_points_x .or. &
    size(um_land_sea_mask, 2) /= um_grid%num_p_points_y) then
  write(log_scratch_space, '(A,4(I0,A))')                          &
       'UM land-sea mask shape (', size(um_land_sea_mask, 1), ',', &
       size(um_land_sea_mask, 2), ') does not match UM P grid (',  &
       um_grid%num_p_points_x, ',', um_grid%num_p_points_y, ')'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

write(log_scratch_space, '(A,I0,A)') 'UM land-sea mask has ', &
     count(um_land_sea_mask), ' land points'
call log_event(log_scratch_space, LOG_LEVEL_INFO)

end subroutine lfric2um_read_land_sea_mask

!> @brief Set sea points of land-compressed fields to missing data
!> @param[in]     stashcode Stashcode of the field
!> @param[in,out] field     Regridded field on the UM P grid
subroutine lfric2um_apply_land_sea_mask(stashcode, field)

! lfricinputs modules
use lfricinp_stashmaster_mod,   only: get_stashmaster_item, grid, &
                                      land_compressed
use lfricinp_um_parameters_mod, only: um_rmdi

implicit none

integer(kind=int64), intent(in)   :: stashcode
real(kind=real64),   intent(in out) :: field(:,:)

! If there is no mask, or this isn't a land-only field, just return
if (.not. allocated(um_land_sea_mask)) return
if (get_stashmaster_item(stashcode, grid) /= land_compressed) return

! Having an explicit loop would be more OpenMP friendly, but as we don't
! use OpenMP in lfricinputs, this is probably the most efficient method
where (.not. um_land_sea_mask) field = um_rmdi

end subroutine lfric2um_apply_land_sea_mask

end module lfric2um_land_sea_mask_mod
