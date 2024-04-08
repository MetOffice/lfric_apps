! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_check_stat_ncdf_mod

IMPLICIT NONE

PRIVATE

PUBLIC :: check_stat_ncdf

CONTAINS

! Subroutine name doesn't follow naming convention 'lfricinp_' as we want
! it to be short due to its use as a wrapper to netcdf function calls
SUBROUTINE check_stat_ncdf(stat)
! Description:
!  Wrapper routine to check return status of any calls to
!  netcdf library and call an abort if necessary

! External libraries
USE netcdf, ONLY: NF90_NOERR, NF90_STRERROR
! LFRic modules
USE log_mod, ONLY: LOG_LEVEL_INFO, LOG_LEVEL_ERROR, log_event

IMPLICIT NONE

INTEGER, INTENT(IN) :: stat

IF (stat /= NF90_NOERR) THEN
  CALL log_event("Issue with call to netcdf library", LOG_LEVEL_INFO)
  CALL log_event(TRIM(NF90_STRERROR(stat)), LOG_LEVEL_ERROR)
END IF

END SUBROUTINE check_stat_ncdf

END MODULE lfricinp_check_stat_ncdf_mod
