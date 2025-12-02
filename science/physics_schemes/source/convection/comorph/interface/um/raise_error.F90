! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module raise_error_mod

implicit none

contains

! Subroutines to handle errors in the CoMorph convection scheme.
! This is the version for CoMorph convection runs within the UM,
! so just contains wrappers for the UM's error-handling routines.


!----------------------------------------------------------------
! Routine to raise a fatal error, killing the run.
!----------------------------------------------------------------
subroutine raise_fatal( routinename_in, error_message )

use ereport_mod, only: ereport

implicit none

! Routine which generated the error
character(len=*), intent(in) :: routinename_in

! Error message to be printed
character(len=*), intent(in) :: error_message

! Error code to pass into the UM's ereport
integer :: error_code


! Call UM's ereport to output the error message and kill the run.
error_code = 1
call ereport( routinename_in, error_code, error_message )


return
end subroutine raise_fatal



!----------------------------------------------------------------
! Routine to raise a non-fatal error, just prints a warning.
!----------------------------------------------------------------
subroutine raise_warning( routinename_in, error_message )

use ereport_mod, only: ereport

implicit none

! Routine which generated the error
character(len=*), intent(in) :: routinename_in

! Error message to be printed
character(len=*), intent(in) :: error_message

! Error code to pass into the UM's ereport
integer :: error_code


! Call UM's ereport to output the error message.
error_code = -1
call ereport( routinename_in, error_code, error_message )


return
end subroutine raise_warning


end module raise_error_mod
