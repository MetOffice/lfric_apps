! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module wtrac_move_phase_mod

use um_types,           only: real_64, real_32
use water_tracers_mod,  only: l_wtrac_limit_chg

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim


implicit none

! Description:
!    Move water tracer between phases due to a single phase change.
!
!    If the phase change amount (qchange) is positive, the phase change
!    is from  type 1 to 2.  If it is negative, the phase change is from
!    type 2 to 1.
!
!    **There are two versions for different precisions of the
!        inputs/outputs.**
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

public :: wtrac_move_phase
interface wtrac_move_phase
  module procedure wtrac_move_phase_real64,                                    &
                   wtrac_move_phase_real32
end interface

character(len=*), parameter, private :: ModuleName='WTRAC_MOVE_PHASE_MOD'

contains

!===========================================================================

! Version for real_64 precision

subroutine wtrac_move_phase_real64(npnts, ni, idx, l_limit_chg, qchange,       &
                                   qratio1_wtrac, qratio2_wtrac,               &
                                   q1_wtrac, q2_wtrac, qchange_wtrac)

implicit none

integer, parameter :: field_kind = real_64

character(len=*), parameter :: RoutineName='WTRAC_MOVE_PHASE_REAL64'

#include "wtrac_move_phase.h"

return
end subroutine wtrac_move_phase_real64

!============================================================================

! Version for real_32 precision

subroutine wtrac_move_phase_real32(npnts, ni, idx, l_limit_chg, qchange,       &
                                   qratio1_wtrac, qratio2_wtrac,               &
                                   q1_wtrac, q2_wtrac, qchange_wtrac)
implicit none

integer, parameter :: field_kind = real_32

character(len=*), parameter :: RoutineName='WTRAC_MOVE_PHASE_REAL32'

#include "wtrac_move_phase.h"

return
end subroutine wtrac_move_phase_real32

end module wtrac_move_phase_mod
