! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_check_shumlib_status_mod

implicit none

private
public :: shumlib

contains

subroutine shumlib(routinename, status,                                        &
                   print_on_success, ignore_warning, errorstatus)

! Each routine in the shumlib fieldsfile API returns a status object.
! This routine checks the status of such an object. It can act
! as a wrapper routine if the shumlib function is called directly
! in the argument list, for example:
!
! call shumlib("my_shumlib_func", my_shumlib_func(var1, var2))
!
! Here my_shumlib_func is called first, the arguments var1 and
! var2 are passed to the routine. They remain in scope in the calling
! routine and are available to be set if their intent allows.
! my_shumlib_func then returns a status object that is passed into this
! routine for checking.

! Intrinsic modules
use, intrinsic :: iso_c_binding, only: c_bool

! LFRic modules
use log_mod, only: log_event, LOG_LEVEL_ERROR, LOG_LEVEL_INFO

! Shumlib modules
use f_shum_ff_status_mod, only: shum_ff_status_type, operator(==),            &
                                operator(<), SHUMLIB_SUCCESS

implicit none

! Arguments
type(shum_ff_status_type), intent(in) :: status
character(len=*),          intent(in) :: routinename
logical(kind=c_bool), optional,  intent(in) :: print_on_success, ignore_warning
integer,              optional, intent(out) :: errorstatus
! Internal variables

! Message - set to be the maximum SHUMlib message length (1024) plus a
! reasonable size for a routine name (128)
character(len=1152) :: message
! Set message
write(message, '(A,A,A,A)') '[', trim(routinename), '] ',trim(status%message)

if (status < SHUMLIB_SUCCESS) then
  ! This is potentially a warning
  if (present(errorstatus)) errorstatus = -1
  if (present(ignore_warning)) then
    if (ignore_warning) then
      ! Ignore warning is true, print and carry on
      call log_event(message, LOG_LEVEL_INFO)
    else
      ! Ignore warning is false, print and abort
      call log_event(message, LOG_LEVEL_ERROR)
    end if
  else
    ! Ignore warning is not set, print and abort
    call log_event(message, LOG_LEVEL_ERROR)
  end if
else if (status == SHUMLIB_SUCCESS) then
  ! This is a success; if print_on_success is provided, and true, print message
  ! otherwise carry on as normal
  if (present(errorstatus)) errorstatus = 0
  if (present(print_on_success)) then
    if (print_on_success) then
      call log_event(message, LOG_LEVEL_INFO)
    end if
  end if
else
  if (present(errorstatus)) errorstatus = 1
  ! This is a definite failure, abort whatever else is provided
  call log_event(message, LOG_LEVEL_ERROR)
end if

end subroutine shumlib

end module lfricinp_check_shumlib_status_mod
