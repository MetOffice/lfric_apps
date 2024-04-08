! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_read_um_time_data_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int64, real64
USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_BOOL

! LFRic modules
USE constants_mod, ONLY: i_def

! LFRic Inputs modules
USE lfricinp_datetime_mod, ONLY: datetime

IMPLICIT NONE

CONTAINS

SUBROUTINE lfricinp_read_um_time_data(um_file, stash_list)
!
! This routine reads, for a list of stash codes, the forecast and validity times
! associated with those stash codes from a UM dump/fieldsfile and forms a list
! of unique forecast and validity times and stores the results in a LFRic Inputs
! datetime type

! LFRic modules
USE log_mod,       ONLY: log_event, LOG_LEVEL_ERROR, LOG_LEVEL_INFO,           &
                         log_scratch_space
! Shumlib modules
USE f_shum_field_mod, ONLY: shum_field_type
USE f_shum_file_mod,  ONLY: shum_file_type
USE f_shum_fixed_length_header_indices_mod, ONLY: calendar
! lfricinp modules
USE lfricinp_check_shumlib_status_mod, ONLY: shumlib

IMPLICIT NONE

TYPE(shum_file_type),    INTENT(IN OUT) :: um_file
INTEGER(KIND=int64),     INTENT(IN)     :: stash_list(:)

! Local variables
INTEGER(KIND=int64), PARAMETER :: calendar_gregorian = 1, calendar_360 = 2,    &
                                  calendar_365 = 3

! Array of shumlib field objects that will be returned from UM file
INTEGER(KIND=int64), ALLOCATABLE  :: um_input_field_indices(:)

INTEGER(KIND=int64) :: stashcode, calendar_type

! Error reporting
CHARACTER(LEN=*), PARAMETER :: routinename='datetime%read_um_time_data'

! Variables needed for forecast time checking
REAL(KIND=real64), PARAMETER :: tol_fct = 1.0e-6 ! fctime tolerance in Shumlib
REAL(KIND=real64) :: fctime, period, pdiff
CHARACTER(LEN=16) :: timestring
INTEGER :: i_field, level, time_idx, time_idx_insert, t_idx, time_idx_max
INTEGER :: errorstatus
LOGICAL :: l_new_fct, l_periodic
LOGICAL(KIND=C_BOOL) :: true_cbool

true_cbool = LOGICAL(.TRUE., KIND=C_BOOL)

! Get forecast and validity times present in the um file for requested
! stashlist, sorting them from smallest forecast to greatest forecast time
! as we go along.
!
time_idx = 0
l_new_fct = .FALSE.
DO i_field = 1, SIZE(stash_list)

  stashcode = stash_list(i_field)

  CALL shumlib("um2lfric::find_field_indices_in_file",                         &
                um_file%find_field_indices_in_file(um_input_field_indices,     &
                stashcode = stashcode, lbproc = 0_int64),                      &
                ignore_warning = true_cbool, errorstatus = errorstatus)

  ! If stashcode is not present in dump, move onto next one
  IF (errorstatus /= 0 ) THEN
    IF (ALLOCATED(um_input_field_indices)) DEALLOCATE(um_input_field_indices)
    CYCLE
  END IF

  DO level = 1, SIZE(um_input_field_indices)

    ! Get forecast and validity time for this specific field
    CALL shumlib("um2lfric::get_real_fctime",                                  &
                  um_file%fields(um_input_field_indices(level)) % get_real_fctime(fctime))
    CALL shumlib("um2lfric::get_timestring",                                   &
                  um_file%fields(um_input_field_indices(level)) % get_timestring(timestring))

    ! Check if forecast time is already in current stored list of forecast
    ! times. If not insert new forecast time into array, preserving ordering.
    l_new_fct = (.NOT. ANY(ABS(datetime%fctimes - fctime) <= tol_fct))
    IF (l_new_fct .OR. time_idx == 0) THEN
      WRITE(log_scratch_space,*) ' STASH code: ', stashcode,                   &
                                 ' has a forecast time of: ', fctime, ' and',  &
                                 ' validity time of: ', timestring
      CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

      ! Note: For very first entry the following two loops will be skipped by
      ! default for the initial values of time_idx = 0 and time_idx_insert = 1

      ! Find location where to insert new forecast time in the array.
      time_idx_insert = 1
      DO t_idx = 1, time_idx
        IF (fctime > datetime % fctimes(t_idx)) time_idx_insert = t_idx + 1
      END DO

      ! Shift existing longer forecast time entries one index up in the array
      DO t_idx = time_idx+1, time_idx_insert+1,-1
        datetime % fctimes(t_idx) = datetime % fctimes(t_idx-1)
        datetime % validity_times(t_idx) = datetime % validity_times(t_idx-1)
      END DO

      ! Now insert existing entry into the array
      datetime % fctimes(time_idx_insert) = fctime
      datetime % validity_times(time_idx_insert) =                             &
            timestring(1:4)//'-'//timestring(5:6)//'-'//timestring(7:8)//' '// &
            timestring(10:11)//':'//timestring(12:13)//':'//timestring(14:15)

      ! Update number of newly found forecast times
      time_idx = time_idx + 1
    ENDIF

  END DO

  IF (ALLOCATED(um_input_field_indices)) DEALLOCATE(um_input_field_indices)

END DO

! Set number of unique forecast/validity times
time_idx_max = time_idx
datetime % num_times = INT(time_idx_max,KIND=int64)

! Set first fctime and the corresponding first validity time
datetime % first_fctime = datetime % fctimes(1)
datetime % first_validity_time = datetime % validity_times(1)

! Get calendar type from fixed length header
CALL shumlib(routinename//'::get_fixed_length_header_by_index',                &
     um_file % get_fixed_length_header_by_index(calendar, calendar_type))
SELECT CASE(calendar_type)
  CASE(calendar_gregorian)
    datetime % calendar = 'Gregorian'
  CASE(calendar_360)
    datetime % calendar = '360'
  CASE(calendar_365)
    datetime % calendar = '365'
  CASE DEFAULT
    CALL log_event('Unrecognised calendar type', LOG_LEVEL_ERROR)
END SELECT

IF (datetime % num_times > 1) THEN
  ! Check forecast times are periodic. Currently LFRic Inputs only assume
  ! timeseries with a fix output period/frequency
  period = datetime % fctimes(2) - datetime % fctimes(1)
  DO t_idx = 3, time_idx_max
    pdiff = ABS(datetime%fctimes(t_idx) - datetime%fctimes(t_idx-1) - period)
    l_periodic = (pdiff < tol_fct)
    IF (.NOT. l_periodic) THEN
      CALL log_event('Only periodic timeseries are valid', LOG_LEVEL_ERROR)
    END IF
  END DO
  WRITE(log_scratch_space, *) 'Forecast time period is ', period, ' hours'
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
ELSE
  ! ... else set period to one hour by default
  period = 1.0
  WRITE(log_scratch_space, *) 'Only single validity time detected. Setting '// &
                              'forecast period to one hour by default'
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
END IF

! Set time step information
datetime % seconds_per_step = period * 3600.0
datetime % first_step = 1
datetime % last_step = INT(datetime % num_times, KIND=i_def)

END SUBROUTINE lfricinp_read_um_time_data

END MODULE lfricinp_read_um_time_data_mod
