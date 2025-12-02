! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module wtrac_all_phase_chg_mod

use um_types,                only: real_64, real_32
use wtrac_calc_ratio_mod,    only: wtrac_calc_ratio_fn
use wtrac_move_phase_mod,    only: wtrac_move_phase

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Description:
!   General routine to update water tracers due to a single phase change.
!   Note, the difference between q1_old and q1 must be due to a single phase
!   change between q1 and q2.
!
!   **There are two versions of the function for different input/output
!     precisions.**
!
! Method:
!   The amount of water changing phase is input to the routine and the
!   water tracers are updated using this input plus the ratio of water
!   tracer to normal water.
!   Fractionation will occur here for water isotopes under certain
!   phase changes.
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

public :: wtrac_all_phase_chg
interface wtrac_all_phase_chg
  module procedure wtrac_all_phase_chg_real64,                                 &
                   wtrac_all_phase_chg_real32
end interface

character(len=*), parameter, private :: ModuleName='WTRAC_ALL_PHASE_CHG_MOD'

contains

!==========================================================================

! Version for real_64 precision

subroutine wtrac_all_phase_chg_real64(npnts, n_wtrac, q1_old, q2_old,          &
                                      q_change, q1, q2,                        &
                                      q1_type, q2_type, direction,             &
                                      q1_wtrac, q2_wtrac,                      &
                                      l_no_fract, l_limit_chg, qchange_wtrac)

implicit none

integer, parameter :: field_kind = real_64

character(len=*), parameter :: RoutineName='WTRAC_ALL_PHASE_CHG_REAL64'

#include "wtrac_all_phase_chg.h"

return
end subroutine wtrac_all_phase_chg_real64

!==========================================================================

! Version for real_32 precision

subroutine wtrac_all_phase_chg_real32(npnts, n_wtrac, q1_old, q2_old,          &
                                      q_change, q1, q2,                        &
                                      q1_type, q2_type, direction,             &
                                      q1_wtrac, q2_wtrac,                      &
                                      l_no_fract, l_limit_chg, qchange_wtrac)

implicit none

integer, parameter :: field_kind = real_32

character(len=*), parameter :: RoutineName='WTRAC_ALL_PHASE_CHG_REAL32'

#include "wtrac_all_phase_chg.h"

return
end subroutine wtrac_all_phase_chg_real32

end module wtrac_all_phase_chg_mod
