! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_check_stat_ncdf_mod

implicit none

private

public :: check_stat_ncdf

contains

! Subroutine name doesn't follow naming convention 'lfricinp_' as we want
! it to be short due to its use as a wrapper to netcdf function calls
subroutine check_stat_ncdf(stat)
! Description:
!  Wrapper routine to check return status of any calls to
!  netcdf library and call an abort if necessary

! External libraries
use netcdf, only: NF90_NOERR, NF90_STRERROR
! LFRic modules
use log_mod, only: LOG_LEVEL_INFO, LOG_LEVEL_ERROR, log_event

implicit none

integer, intent(in) :: stat

if (stat /= NF90_NOERR) then
  call log_event("Issue with call to netcdf library", LOG_LEVEL_INFO)
  call log_event(trim(NF90_STRERROR(stat)), LOG_LEVEL_ERROR)
end if

end subroutine check_stat_ncdf

end module lfricinp_check_stat_ncdf_mod
