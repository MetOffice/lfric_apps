! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation
module mphys_die

use errormessagelength_mod, only: errormessagelength

implicit none

character(len=*), parameter, private :: ModuleName='MPHYS_DIE'

! Integer constants for various CASIM errors
! These are used within the CASIM repository

integer, parameter :: incorrect_opt = 1
integer, parameter :: bad_values    = 2
integer, parameter :: warn          = -1

! Length of a standard CASIM message
! This is deliberately shorter than errormessagelength
! to allow for additional CASIM information to be added.
integer, parameter :: std_msg_len = errormessagelength - 150

! Standard CASIM message, used exclusively within CASIM
! The same variable is available in mphys_die on the
! CASIM repository for CASIM runs not made with the UM.
! It is preferable to use std_msg within the CASIM code
! but revert to using cmessage or ummessage within
! UM routines.
character(len=std_msg_len) :: std_msg = ''

contains

subroutine throw_mphys_error(itype, casim_routine, info)

! If modifying the subroutine or argument list, ensure that the
! CASIM version of mphys_die (on the CASIM repository)
! is also modified to give the same result

use ereport_mod, only: ereport
use umprintmgr,  only: newline

use errormessagelength_mod, only: errormessagelength

implicit none

integer, intent(in) :: itype
! Type of error:
! 1:  Incorrect specification of options
! 2:  Bad values found
! 3+:  Unknown error
! <0: Warning, code prints out a warning to ereport and continues

character(len=*), intent(in) :: casim_routine
! CASIM Routine causing the error or warning

character(len=std_msg_len), intent(in) :: info
! Additional error or warning information

! Local variables
integer, parameter :: um_error_flag = 100
integer, parameter :: um_warning   = -100
! Positive value for call to ereport, negative for warning

character(len=errormessagelength) :: message_str = ''
! Error or Warning message string

! character(len=*), parameter   :: RoutineName='THROW_MPHYS_ERROR'

integer :: errcode ! Error code

!----------------------------------------------------------------------------

if ( itype > 0 ) then

  !  --------------------------------------------
  !                    error
  !  Code must output an error message and abort
  !  --------------------------------------------

  if ( itype == incorrect_opt ) then
    message_str = 'error in CASIM microphysics:' //newline//                   &
      'Incorrect specification of options.'      //newline//                   &
      'Additional information:'                  //newline//                   &
      trim(info)

  else if ( itype == bad_values ) then

    message_str = 'error in CASIM microphysics:' //newline//                   &
      'Bad values found.'                        //newline//                   &
      'Additional information:'                  //newline//                   &
      trim(info)

  else

    !   Unknown error. Just report a general error message

    message_str = 'error in CASIM microphysics:' //newline//                   &
      'Additional information:'                  //newline//                   &
      trim(info)

  end if ! itype

  errcode = um_error_flag

  call ereport(casim_routine, errcode, message_str)

else if ( itype < 0 ) then

  !  --------------------------------------------
  !                  WARNING
  !  Output a warning via ereport and continue
  !  --------------------------------------------

  message_str = 'WARNING from CASIM microphysics:' //newline//                 &
                 trim(info)

  errcode = um_warning

  call ereport(casim_routine, errcode, message_str)

end if

end subroutine throw_mphys_error


subroutine mphys_message(casim_routine, casim_message)

use umprintmgr,  only: umPrint, umMessage

implicit none

character( * ), intent(in) :: casim_routine
! CASIM Routine causing the error

character( * ), intent(in) :: casim_message
! Additional error information

!----------------------------------------------------------------------------

write(umMessage, '(A)') 'CASIM Message | '//trim(casim_routine)//              &
                         trim(casim_message)

call umPrint( umMessage, src='mphys_message' )

end subroutine mphys_message

end module mphys_die
