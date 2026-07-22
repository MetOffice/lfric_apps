! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

#if !defined(LFRIC)
module raise_error_mod

implicit none

contains

! Subroutines to handle errors in the CoMorph convection scheme.
! This version is for stand-alone CoMorph runs, where there are
! no error-handling routines to borrow from the host-model.


!----------------------------------------------------------------
! Routine to raise a fatal error, killing the run.
!----------------------------------------------------------------
subroutine raise_fatal( routinename_in, error_message )

use comorph_constants_mod, only: newline
use, intrinsic :: iso_fortran_env, only: error_unit

implicit none

! Routine which generated the error
character(len=*), intent(in) :: routinename_in

! Error message to be printed
character(len=*), intent(in) :: error_message


! Print the error message to standard error
write(error_unit,"(A)")                                                        &
  "Error in subroutine " // routinename_in // ":" // newline //                &
  error_message

! Kill the run
stop


return
end subroutine raise_fatal


!----------------------------------------------------------------
! Routine to raise a non-fatal error, just prints a warning.
!----------------------------------------------------------------
subroutine raise_warning( routinename_in, error_message )

use comorph_constants_mod, only: newline
use, intrinsic :: iso_fortran_env, only: error_unit

implicit none

! Routine which generated the error
character(len=*), intent(in) :: routinename_in

! Error message to be printed
character(len=*), intent(in) :: error_message


! Print the warning message to standard error
write(error_unit,"(A)")                                                        &
  "Warning in subroutine " // routinename_in // ":" // newline //              &
  error_message


return
end subroutine raise_warning


end module raise_error_mod
#endif
