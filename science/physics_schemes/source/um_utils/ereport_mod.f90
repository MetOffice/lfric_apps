!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!> @brief UM interface to LFRic error reporting
!>
!> @details call log_event to signal errors or warnings
module ereport_mod

  implicit none

  public ereport

  contains

    subroutine ereport( routine_name, error_status, message )

      use log_mod, only: log_event,        &
                         LOG_LEVEL_INFO,   &
                         LOG_LEVEL_ERROR,  &
                         LOG_LEVEL_WARNING

      implicit none

      character (len=*), intent(in) :: routine_name
      integer, intent(inout)        :: error_status
      character (len=*), intent(in) :: message

      ! call log_event with different flag depending on UM errorstatus
      ! convention

      if (error_status > 0) then
        call log_event(routine_name//message,LOG_LEVEL_ERROR)
      else if (error_status < 0) then
        call log_event(routine_name//message,LOG_LEVEL_WARNING)
      else if (error_status == 0) then
        call log_event(routine_name//message,LOG_LEVEL_INFO)
      end if

      ! reset error_status
      error_status = 0

    end subroutine ereport

end module ereport_mod
