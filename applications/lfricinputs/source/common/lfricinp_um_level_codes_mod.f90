! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_um_level_codes_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only: int64

! Shumlib modules
use f_shum_file_mod, only: shum_file_type

! lfricinputs modules
use lfricinp_check_shumlib_status_mod, only: shumlib
use lfricinp_grid_type_mod, only: lfricinp_grid_type

! lfric modules
use log_mod,  only : log_event, log_scratch_space,         &
                     LOG_LEVEL_ERROR

implicit none

private

public :: lfricinp_get_num_levels, lfricinp_get_first_level_num, &
     lfricinp_get_last_level_num, lfricinp_get_num_pseudo_levels


contains

!------------------------------------------------------------------

function lfricinp_get_num_levels(um_file, stashcode) result(num_levels)
! Description:
!  Returns the number of levels expected in a field as defined by
!  first and last level codes in the stashmaster

implicit none

type(shum_file_type), intent(in) :: um_file
integer(kind=int64), intent(in) :: stashcode
! Result
integer(kind=int64) :: num_levels

num_levels = lfricinp_get_last_level_num(um_file, stashcode) - &
             lfricinp_get_first_level_num(stashcode) + 1

end function lfricinp_get_num_levels

!------------------------------------------------------------------

function lfricinp_get_first_level_num(stashcode) result(first_level_num)
! Description:
!  Takes stashcode as input and interogates the stashmaster first
!  level code to determine the first/bottom level number for the field

! lfricinp modules
use lfricinp_stashmaster_mod, only: get_stashmaster_item, levelf

implicit none
integer(kind=int64), intent(in) :: stashcode
! Result
integer(kind=int64) :: first_level_num
! Local variables
integer(kind=int64) :: first_level_code

first_level_code =  get_stashmaster_item(stashcode, levelf)

select case(first_level_code)
  case (-1) ! Unset /single level
    first_level_num = 1
  case (1) ! First atmos level
    first_level_num = 1
  case (8) ! First soil level
    first_level_num = 1
  case (38,40) ! Zeroth atmos level
    first_level_num = 0
case DEFAULT
  write(log_scratch_space, '(A,I0,A)') &
     "First level code ", first_level_code, " not supported"
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end select

end function lfricinp_get_first_level_num

!------------------------------------------------------------------

function lfricinp_get_last_level_num(um_file, stashcode) result(last_level_num)
! Description:
!  Takes um_file and stashcode as input and interogates the stashmaster
!  last level code to determine the last/top level number for the field

! lfricinp modules
use lfricinp_stashmaster_mod, only: get_stashmaster_item, levell
use lfricinp_um_parameters_mod, only: ih_model_levels, ih_soilT_levels
implicit none
type(shum_file_type), intent(in) :: um_file
integer(kind=int64), intent(in) :: stashcode
! Result
integer(kind=int64) :: last_level_num
! Local variables
integer(kind=int64) :: last_level_code
integer(kind=int64) :: model_levels
character(len=*), parameter :: routinename='lfricinp_get_last_level_num'


last_level_code =  get_stashmaster_item(stashcode, levell)
! Get model_levels
call shumlib(routinename //'::get_integer_constants_by_index', &
     um_file%get_integer_constants_by_index(ih_model_levels,   &
     model_levels))

select case(last_level_code)
  case (-1) ! Unset /single level
    last_level_num = 1
  case (2,3) ! Top atmos level
    last_level_num = model_levels
  case (9) ! Last soil level
    call shumlib(routinename //'::get_integer_constants_by_index', &
     um_file%get_integer_constants_by_index(ih_soilT_levels,   &
     last_level_num))
  case(19) ! Number of exner levels includes level above model top
    last_level_num = model_levels
  case(23) ! Top ozone level
    last_level_num = model_levels
case DEFAULT
  write(log_scratch_space, '(A,I0,A)') &
     "Last level code ", last_level_code, " not supported"
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end select

end function lfricinp_get_last_level_num

!------------------------------------------------------------------

function lfricinp_get_num_pseudo_levels(um_grid, stashcode) &
     result(num_pseudo_levels)
! Description:
!  Returns the number of pseudo levels expected in a field as defined by
!  first and last pseudo level codes in the stashmaster

implicit none

type(lfricinp_grid_type), intent(in):: um_grid
integer(kind=int64), intent(in) :: stashcode
! Result
integer(kind=int64) :: num_pseudo_levels

num_pseudo_levels = lfricinp_get_last_pseudo_level_num(um_grid, stashcode) - &
             lfricinp_get_first_pseudo_level_num(stashcode) + 1

end function lfricinp_get_num_pseudo_levels

!------------------------------------------------------------------

function lfricinp_get_first_pseudo_level_num(stashcode) &
     result(first_pseudo_level_num)
! Description:
!  Takes stashcode as input and interogates the stashmaster first
!  pseudo_level code to determine the first/bottom pseudo_level
!  number for the field

! lfricinp modules
use lfricinp_stashmaster_mod, only: get_stashmaster_item, pseudf

implicit none
integer(kind=int64), intent(in) :: stashcode
! Result
integer(kind=int64) :: first_pseudo_level_num
! Local variables
integer(kind=int64) :: first_pseudo_level_code

first_pseudo_level_code =  get_stashmaster_item(stashcode, pseudf)

select case(first_pseudo_level_code)
  case (1) ! Dimension starts at 1
    first_pseudo_level_num = 1
case DEFAULT
  write(log_scratch_space, '(A,I0,A)') &
     "First pseudo_level code ", first_pseudo_level_code, " not supported"
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end select

end function lfricinp_get_first_pseudo_level_num

!------------------------------------------------------------------

function lfricinp_get_last_pseudo_level_num(um_grid, stashcode) &
     result(last_pseudo_level_num)
! Description:
!  Takes um_grid and stashcode as input and interogates the stashmaster
!  last pseudo_level code to determine the last/top pseudo_level number
!  for the field

! lfricinp modules
use lfricinp_stashmaster_mod, only: get_stashmaster_item, pseudl

implicit none
type(lfricinp_grid_type), intent(in) :: um_grid
integer(kind=int64), intent(in) :: stashcode
! Result
integer(kind=int64) :: last_pseudo_level_num
! Local variables
integer(kind=int64) :: last_pseudo_level_code
character(len=*), parameter :: routinename = &
     'lfricinp_get_last_pseudo_level_num'

last_pseudo_level_code =  get_stashmaster_item(stashcode, pseudl)

select case(last_pseudo_level_code)
case (7,9) ! ntypes == ntiles (lfricinputs doesn't support aggregate tile)
  last_pseudo_level_num = um_grid % num_surface_types
case (10)
  last_pseudo_level_num = um_grid%num_ice_cats
case (11)
  last_pseudo_level_num = um_grid%num_snow_layers * um_grid%num_surface_types
case DEFAULT
  write(log_scratch_space, '(A,I0,A,I0)') &
       "Last pseudo_level code ", last_pseudo_level_code, &
       " not supported for stashcode ",stashcode
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end select

end function lfricinp_get_last_pseudo_level_num


end module lfricinp_um_level_codes_mod
