! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_read_um_time_data_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int64, real64
use, intrinsic :: iso_c_binding, only: c_bool

! LFRic modules
use constants_mod, only: i_def

! LFRic Inputs modules
use lfricinp_datetime_mod, only: datetime

implicit none

contains

subroutine lfricinp_read_um_time_data(um_file, stash_list)
!
! This routine reads, for a list of stash codes, the forecast and validity times
! associated with those stash codes from a UM dump/fieldsfile and forms a list
! of unique forecast and validity times and stores the results in a LFRic Inputs
! datetime type

! LFRic modules
use log_mod,       only: log_event, LOG_LEVEL_ERROR, LOG_LEVEL_INFO,           &
                         log_scratch_space
! Shumlib modules
use f_shum_field_mod, only: shum_field_type
use f_shum_file_mod,  only: shum_file_type
use f_shum_fixed_length_header_indices_mod, only: calendar
! lfricinp modules
use lfricinp_check_shumlib_status_mod, only: shumlib

implicit none

type(shum_file_type),    intent(in out) :: um_file
integer(kind=int64),     intent(in)     :: stash_list(:)

! Local variables
integer(kind=int64), parameter :: calendar_gregorian = 1, calendar_360 = 2,    &
                                  calendar_365 = 3

! Array of shumlib field objects that will be returned from UM file
integer(kind=int64), allocatable  :: um_input_field_indices(:)

integer(kind=int64) :: stashcode, calendar_type

! Error reporting
character(len=*), parameter :: routinename='datetime%read_um_time_data'

! Variables needed for forecast time checking
real(kind=real64), parameter :: tol_fct = 1.0e-6_real64 ! fctime tolerance in Shumlib
real(kind=real64) :: fctime, period, pdiff
character(len=16) :: timestring
integer :: i_field, level, time_idx, time_idx_insert, t_idx, time_idx_max
integer :: errorstatus
logical :: l_new_fct, l_periodic
logical(kind=c_bool) :: true_cbool

true_cbool = logical(.true., kind=c_bool)

! Get forecast and validity times present in the um file for requested
! stashlist, sorting them from smallest forecast to greatest forecast time
! as we go along.
!
time_idx = 0
l_new_fct = .false.
do i_field = 1, size(stash_list)

  stashcode = stash_list(i_field)

  call shumlib("um2lfric::find_field_indices_in_file",                         &
                um_file%find_field_indices_in_file(um_input_field_indices,     &
                stashcode = stashcode, lbproc = 0_int64),                      &
                ignore_warning = true_cbool, errorstatus = errorstatus)

  ! If stashcode is not present in dump, move onto next one
  if (errorstatus /= 0 ) then
    if (allocated(um_input_field_indices)) deallocate(um_input_field_indices)
    cycle
  end if

  do level = 1, size(um_input_field_indices)

    ! Get forecast and validity time for this specific field
    call shumlib("um2lfric::get_real_fctime",                                  &
                  um_file%fields(um_input_field_indices(level)) % get_real_fctime(fctime))
    call shumlib("um2lfric::get_timestring",                                   &
                  um_file%fields(um_input_field_indices(level)) % get_timestring(timestring))

    ! Check if forecast time is already in current stored list of forecast
    ! times. If not insert new forecast time into array, preserving ordering.
    l_new_fct = (.not. any(abs(datetime%fctimes - fctime) <= tol_fct))
    if (l_new_fct .or. time_idx == 0) then
      write(log_scratch_space,*) ' STASH code: ', stashcode,                   &
                                 ' has a forecast time of: ', fctime, ' and',  &
                                 ' validity time of: ', timestring
      call log_event(log_scratch_space, LOG_LEVEL_INFO)

      ! Note: For very first entry the following two loops will be skipped by
      ! default for the initial values of time_idx = 0 and time_idx_insert = 1

      ! Find location where to insert new forecast time in the array.
      time_idx_insert = 1
      do t_idx = 1, time_idx
        if (fctime > datetime % fctimes(t_idx)) time_idx_insert = t_idx + 1
      end do

      ! Shift existing longer forecast time entries one index up in the array
      do t_idx = time_idx+1, time_idx_insert+1,-1
        datetime % fctimes(t_idx) = datetime % fctimes(t_idx-1)
        datetime % validity_times(t_idx) = datetime % validity_times(t_idx-1)
      end do

      ! Now insert existing entry into the array
      datetime % fctimes(time_idx_insert) = fctime
      datetime % validity_times(time_idx_insert) =                             &
            timestring(1:4)//'-'//timestring(5:6)//'-'//timestring(7:8)//' '// &
            timestring(10:11)//':'//timestring(12:13)//':'//timestring(14:15)

      ! Update number of newly found forecast times
      time_idx = time_idx + 1
    ENDIF

  end do

  if (allocated(um_input_field_indices)) deallocate(um_input_field_indices)

end do

! Set number of unique forecast/validity times
time_idx_max = time_idx
datetime % num_times = int(time_idx_max,kind=int64)

! Set first fctime and the corresponding first validity time
datetime % first_fctime = datetime % fctimes(1)
datetime % first_validity_time = datetime % validity_times(1)

! Get calendar type from fixed length header
call shumlib(routinename//'::get_fixed_length_header_by_index',                &
     um_file % get_fixed_length_header_by_index(calendar, calendar_type))
select case(calendar_type)
  case(calendar_gregorian)
    datetime % calendar = 'Gregorian'
  case(calendar_360)
    datetime % calendar = '360'
  case(calendar_365)
    datetime % calendar = '365'
  case DEFAULT
    call log_event('Unrecognised calendar type', LOG_LEVEL_ERROR)
end select

if (datetime % num_times > 1) then
  ! Check forecast times are periodic. Currently LFRic Inputs only assume
  ! timeseries with a fix output period/frequency
  period = datetime % fctimes(2) - datetime % fctimes(1)
  do t_idx = 3, time_idx_max
    pdiff = abs(datetime%fctimes(t_idx) - datetime%fctimes(t_idx-1) - period)
    l_periodic = (pdiff < tol_fct)
    if (.not. l_periodic) then
      call log_event('Only periodic timeseries are valid', LOG_LEVEL_ERROR)
    end if
  end do
  write(log_scratch_space, *) 'Forecast time period is ', period, ' hours'
  call log_event(log_scratch_space, LOG_LEVEL_INFO)
else
  ! ... else set period to one hour by default
  period = 1.0_real64
  write(log_scratch_space, *) 'Only single validity time detected. Setting '// &
                              'forecast period to one hour by default'
  call log_event(log_scratch_space, LOG_LEVEL_INFO)
end if

! Set time step information
datetime % seconds_per_step = period * 3600.0_real64
datetime % first_step = 1
datetime % last_step = int(datetime % num_times, kind=i_def)

end subroutine lfricinp_read_um_time_data

end module lfricinp_read_um_time_data_mod
