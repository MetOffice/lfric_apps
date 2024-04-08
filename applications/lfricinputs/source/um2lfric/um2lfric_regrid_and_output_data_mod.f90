! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE um2lfric_regrid_and_output_data_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY: real64, int64

! lfricinputs modules
USE lfricinp_lfric_driver_mod,               ONLY: model_clock, lfric_fields
USE lfricinp_datetime_mod,                   ONLY: datetime_type

! um2lfric modules
USE um2lfric_regrid_fields_mod,    ONLY: um2lfric_regrid_fields

! LFRic modules
USE log_mod,                       ONLY: log_event, log_scratch_space,         &
                                         LOG_LEVEL_INFO, LOG_LEVEL_ERROR
USE lfric_xios_write_mod,          ONLY: write_state

! External libraries
USE xios,                          ONLY: xios_date_convert_to_string,          &
                                         xios_get_current_date, xios_date,     &
                                         xios_context_finalize
USE mod_wait,                      ONLY: init_wait

IMPLICIT NONE

PRIVATE

PUBLIC :: um2lfric_regrid_and_output_data

CONTAINS

SUBROUTINE um2lfric_regrid_and_output_data(datetime)
!
! This routine, for each forecast time, regrids a set of UM fields to a lfric
! field collection followed by a write to the data output file/stream. In a
! final step it calls the time axis adjustment routine, which will correct the
! time axis values if needed.
!

TYPE(datetime_type),    INTENT(IN)    :: datetime

TYPE(xios_date)              :: xios_current_date
CHARACTER(LEN=32)            :: xios_current_date_str
REAL(KIND=real64)            :: fctime
INTEGER(KIND=int64)          :: time_step
LOGICAL                      :: l_advance

! NOTE: For the logic below the calendar/clock time step is one unit  ahead of
! the actual forecast time. This will be corrected in a post-processing step.
! This is a current work around the fact that XIOS does not appear allow fields
! to be written before time step 1.
DO time_step = datetime % first_step, datetime % last_step

  l_advance = model_clock % tick()
  IF (.NOT. l_advance) THEN
    WRITE(log_scratch_space,'(A,I0)') 'Failed to advance clock on time step ', &
                                      time_step
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END IF

  CALL xios_get_current_date(xios_current_date)
  CALL xios_date_convert_to_string(xios_current_date, xios_current_date_str)
  fctime = datetime % fctimes(time_step)

  WRITE(log_scratch_space,*) 'Regrid fields for forecast time: ',              &
                              fctime,                                          &
                              '... where current xios date is: ',              &
                              TRIM(xios_current_date_str)
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
  CALL um2lfric_regrid_fields(fctime = fctime)
  CALL write_state(lfric_fields)

END DO

! Finalizes XIOS file context thus forcing data out of IO buffers
CALL xios_context_finalize()
! We have closed the context on our end, but we need to make sure that XIOS
! has closed the files for all servers before we process them
CALL init_wait()
CALL log_event('Finalise XIOS context', LOG_LEVEL_INFO)

! Post process correct output file by offsetting time axis by one time step
CALL adjust_time_axis(time_axis_offset = datetime % seconds_per_step)

END SUBROUTINE um2lfric_regrid_and_output_data


SUBROUTINE adjust_time_axis(time_axis_offset)
!
! This routine adjusts the time-axis values in the output file by the
! given input time duration. It is assumed the time duration has the same units
! as the time axis values, and the time axis is named 'time'. If no time axis
! variable with that name is found a warning message will be reported to that
! effect, and the output file will remain unaltered.
!

! External libraries
USE netcdf,             ONLY: nf90_noerr, nf90_open, nf90_write, nf90_close,   &
                              nf90_inq_dimid, nf90_inquire_dimension,          &
                              nf90_inq_varid, nf90_get_var, nf90_put_var
USE mpi_mod,            ONLY: global_mpi

! LFRic modules
USE constants_mod,      ONLY: i_def, r_second

! NetCDF wrapper routine
USE lfricinp_check_stat_ncdf_mod, ONLY: check_stat_ncdf

! Path to output file
USE lfricinp_setup_io_mod, ONLY: io_config

  IMPLICIT NONE

  REAL(KIND=r_second),     INTENT(IN) :: time_axis_offset

  INTEGER(KIND=i_def)                 :: local_rank
  INTEGER                             :: ncid, varid, dimid, dim_size,         &
                                         check_time_axis_status
  CHARACTER(LEN=:),       ALLOCATABLE :: file_path
  REAL(KIND=r_second),    ALLOCATABLE :: temp_time(:), temp_time_bounds(:,:)

  local_rank = global_mpi%get_comm_rank()

  IF (local_rank == 0) THEN

    file_path = TRIM(ADJUSTL(io_config%checkpoint_write_file)) // '.nc'

    WRITE(log_scratch_space,*) 'NetCDF file to be processed: ', file_path
    CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

    ! Open NetCDF file
    CALL check_stat_ncdf(nf90_open(path=file_path, mode=nf90_write, ncid=ncid))

    ! Check is time variable exist in file
    check_time_axis_status = nf90_inq_varid(ncid, 'time', varid)

    ! Based on whether time axis variable exist, update axis or do nothing
    IF (check_time_axis_status == nf90_noerr) THEN ! Time axis found

      ! Get size of dimension named "time"
      CALL check_stat_ncdf(nf90_inq_dimid(ncid, 'time', dimid))
      CALL check_stat_ncdf(nf90_inquire_dimension(ncid, dimid, len=dim_size))

      ! Allocate temporary time data array
      ALLOCATE(temp_time(dim_size))

      ! Get size of dimension named "axis_nbounds"
      CALL check_stat_ncdf(nf90_inq_dimid(ncid, 'axis_nbounds', dimid))
      CALL check_stat_ncdf(nf90_inquire_dimension(ncid, dimid, len=dim_size))

      ! Allocate temporary time_bounds data array
      ALLOCATE(temp_time_bounds(dim_size, SIZE(temp_time)))

      ! Get "time" variable id
      CALL check_stat_ncdf(nf90_inq_varid(ncid, 'time', varid))

      ! Get data for "time" variable using variable id
      CALL check_stat_ncdf(nf90_get_var(ncid, varid, temp_time))

      ! Shift time data values
      temp_time = temp_time - time_axis_offset

      ! Write shifted values back to file
      CALL check_stat_ncdf(nf90_put_var(ncid, varid, temp_time))

      ! Get "time_bounds" variable id
      CALL check_stat_ncdf(nf90_inq_varid(ncid, 'time_bounds', varid))

      ! Get data for "time_bounds" variable using variable id
      CALL check_stat_ncdf(nf90_get_var(ncid, varid, temp_time_bounds))

      ! Shift time_bounds data values
      temp_time_bounds = temp_time_bounds - time_axis_offset

      ! Write shifted values back to file
      CALL check_stat_ncdf(nf90_put_var(ncid, varid, temp_time_bounds))

      ! Deallocate temporary arrays
      DEALLOCATE(temp_time, temp_time_bounds)

    ELSE ! Time axis variable not found

      WRITE(log_scratch_space,*) 'No time variable found in file'
      CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

    END IF

    ! Close NetCDF file
    CALL check_stat_ncdf(nf90_close(ncid))

  END IF

  END SUBROUTINE adjust_time_axis

END MODULE um2lfric_regrid_and_output_data_mod
