! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module um2lfric_regrid_and_output_data_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only: real64, int64

! lfricinputs modules
use lfricinp_lfric_driver_mod,               only: model_clock, lfric_fields
use lfricinp_datetime_mod,                   only: datetime_type

! um2lfric modules
use um2lfric_regrid_fields_mod,    only: um2lfric_regrid_fields

! LFRic modules
use log_mod,                       only: log_event, log_scratch_space,         &
                                         LOG_LEVEL_INFO, LOG_LEVEL_ERROR
use lfric_xios_write_mod,          only: write_state

! External libraries
use xios,                          only: xios_date_convert_to_string,          &
                                         xios_get_current_date, xios_date,     &
                                         xios_context_finalize
use mod_wait,                      only: init_wait

implicit none

private

public :: um2lfric_regrid_and_output_data

contains

subroutine um2lfric_regrid_and_output_data(datetime)
!
! This routine, for each forecast time, regrids a set of UM fields to a lfric
! field collection followed by a write to the data output file/stream. In a
! final step it calls the time axis adjustment routine, which will correct the
! time axis values if needed.
!

type(datetime_type),    intent(in)    :: datetime

type(xios_date)              :: xios_current_date
character(len=32)            :: xios_current_date_str
real(kind=real64)            :: fctime
integer(kind=int64)          :: time_step
logical                      :: l_advance

! NOTE: For the logic below the calendar/clock time step is one unit  ahead of
! the actual forecast time. This will be corrected in a post-processing step.
! This is a current work around the fact that XIOS does not appear allow fields
! to be written before time step 1.
do time_step = datetime % first_step, datetime % last_step

  l_advance = model_clock % tick()
  if (.not. l_advance) then
    write(log_scratch_space,'(A,I0)') 'Failed to advance clock on time step ', &
                                      time_step
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

  call xios_get_current_date(xios_current_date)
  call xios_date_convert_to_string(xios_current_date, xios_current_date_str)
  fctime = datetime % fctimes(time_step)

  write(log_scratch_space,*) 'Regrid fields for forecast time: ',              &
                              fctime,                                          &
                              '... where current xios date is: ',              &
                              trim(xios_current_date_str)
  call log_event(log_scratch_space, LOG_LEVEL_INFO)
  call um2lfric_regrid_fields(fctime = fctime)
  call write_state(lfric_fields)

end do

! Finalizes XIOS file context thus forcing data out of IO buffers
call xios_context_finalize()
! We have closed the context on our end, but we need to make sure that XIOS
! has closed the files for all servers before we process them
call init_wait()
call log_event('Finalise XIOS context', LOG_LEVEL_INFO)

! Post process correct output file by offsetting time axis by one time step
call adjust_time_axis(time_axis_offset = datetime % seconds_per_step)

end subroutine um2lfric_regrid_and_output_data


subroutine adjust_time_axis(time_axis_offset)
!
! This routine adjusts the time-axis values in the output file by the
! given input time duration. It is assumed the time duration has the same units
! as the time axis values, and the time axis is named 'time'. If no time axis
! variable with that name is found a warning message will be reported to that
! effect, and the output file will remain unaltered.
!

! External libraries
use netcdf,             only: nf90_noerr, nf90_open, nf90_write, nf90_close,   &
                              nf90_inq_dimid, nf90_inquire_dimension,          &
                              nf90_inq_varid, nf90_get_var, nf90_put_var
use lfric_mpi_mod,      only: global_mpi

! LFRic modules
use constants_mod,      only: i_def, r_second

! NetCDF wrapper routine
use lfricinp_check_stat_ncdf_mod, only: check_stat_ncdf

! Path to output file
use lfricinp_setup_io_mod, only: io_config

  implicit none

  real(kind=r_second),     intent(in) :: time_axis_offset

  integer(kind=i_def)                 :: local_rank
  integer                             :: ncid, varid, dimid, dim_size,         &
                                         check_time_axis_status
  character(len=:),       allocatable :: file_path
  real(kind=r_second),    allocatable :: temp_time(:), temp_time_bounds(:,:)

  local_rank = global_mpi%get_comm_rank()

  if (local_rank == 0) then

    file_path = trim(adjustl(io_config%checkpoint_write_file)) // '.nc'

    write(log_scratch_space,*) 'NetCDF file to be processed: ', file_path
    call log_event(log_scratch_space, LOG_LEVEL_INFO)

    ! Open NetCDF file
    call check_stat_ncdf(nf90_open(path=file_path, mode=nf90_write, ncid=ncid))

    ! Check is time variable exist in file
    check_time_axis_status = nf90_inq_varid(ncid, 'time', varid)

    ! Based on whether time axis variable exist, update axis or do nothing
    if (check_time_axis_status == nf90_noerr) then ! Time axis found

      ! Get size of dimension named "time"
      call check_stat_ncdf(nf90_inq_dimid(ncid, 'time', dimid))
      call check_stat_ncdf(nf90_inquire_dimension(ncid, dimid, len=dim_size))

      ! Allocate temporary time data array
      allocate(temp_time(dim_size))

      ! Get size of dimension named "axis_nbounds"
      call check_stat_ncdf(nf90_inq_dimid(ncid, 'axis_nbounds', dimid))
      call check_stat_ncdf(nf90_inquire_dimension(ncid, dimid, len=dim_size))

      ! Allocate temporary time_bounds data array
      allocate(temp_time_bounds(dim_size, size(temp_time)))

      ! Get "time" variable id
      call check_stat_ncdf(nf90_inq_varid(ncid, 'time', varid))

      ! Get data for "time" variable using variable id
      call check_stat_ncdf(nf90_get_var(ncid, varid, temp_time))

      ! Shift time data values
      temp_time = temp_time - time_axis_offset

      ! Write shifted values back to file
      call check_stat_ncdf(nf90_put_var(ncid, varid, temp_time))

      ! Get "time_bounds" variable id
      call check_stat_ncdf(nf90_inq_varid(ncid, 'time_bounds', varid))

      ! Get data for "time_bounds" variable using variable id
      call check_stat_ncdf(nf90_get_var(ncid, varid, temp_time_bounds))

      ! Shift time_bounds data values
      temp_time_bounds = temp_time_bounds - time_axis_offset

      ! Write shifted values back to file
      call check_stat_ncdf(nf90_put_var(ncid, varid, temp_time_bounds))

      ! Deallocate temporary arrays
      deallocate(temp_time, temp_time_bounds)

    else ! Time axis variable not found

      write(log_scratch_space,*) 'No time variable found in file'
      call log_event(log_scratch_space, LOG_LEVEL_INFO)

    end if

    ! Close NetCDF file
    call check_stat_ncdf(nf90_close(ncid))

  end if

  end subroutine adjust_time_axis

end module um2lfric_regrid_and_output_data_mod
