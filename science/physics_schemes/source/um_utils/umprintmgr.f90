!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!> @brief UM interface to LFRic printing
!>
!> @details call log_event for info only
module umprintmgr

  implicit none

  public

! The following public variables/parameters are all used by
! the UM physics code, so there names/values should not change

! ========================
! PARAMETERS
! ========================
  INTEGER, PARAMETER            :: PrStatus_Min    = 1  ! Minimum output
  INTEGER, PARAMETER            :: PrStatus_Normal = 2  ! Short info output
  INTEGER, PARAMETER            :: PrStatus_Oper   = 3  ! Full info output
  INTEGER, PARAMETER            :: PrStatus_Diag   = 4  ! Diagnostic output
!Shorter names for future use to avoid code clutter
  INTEGER, PARAMETER            :: PrMin=PrStatus_Min
  INTEGER, PARAMETER            :: PrNorm=PrStatus_Normal
  INTEGER, PARAMETER            :: PrOper=PrStatus_Oper
  INTEGER, PARAMETER            :: PrDiag=PrStatus_Diag
! parameterisation of storage/buffering ammounts
  INTEGER, PARAMETER            :: maxLineLen=1024
! Declare newline character
  CHARACTER(LEN=1), PARAMETER   :: newline=NEW_LINE('a')
! ========================
! A buffer that clients
! can use for list
! directed writes
! ========================
  CHARACTER (LEN=maxLineLen)    :: umMessage
!$OMP THREADPRIVATE (umMessage)
! ========================
! Runtime variables and
! initial defaults.
! ========================
  INTEGER                       :: PrintStatus     = PrStatus_Normal

  contains

    subroutine umprint(line, level, pe, src, model,  &
                       UsrPrefix, HangIndent, stdErrorToo)

      use log_mod, only: log_event, LOG_LEVEL_INFO

      implicit none
      CHARACTER(LEN=*)            :: line
      INTEGER, OPTIONAL           :: level
      INTEGER, OPTIONAL           :: pe
      CHARACTER(LEN=*), OPTIONAL  :: src
      CHARACTER(LEN=*), OPTIONAL  :: model
      CHARACTER(LEN=*), OPTIONAL  :: UsrPrefix       ! User prefix for each line
      CHARACTER(LEN=*), OPTIONAL  :: HangIndent      ! Hanging indent for new lines
      LOGICAL, OPTIONAL           :: stdErrorToo

      ! call log_event with info only
      call log_event(line,LOG_LEVEL_INFO)

    end subroutine umprint

    ! Dummy routine
    subroutine umprintflush
    implicit none
    end subroutine umprintflush

end module umprintmgr
