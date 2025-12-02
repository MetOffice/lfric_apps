! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module wtrac_calc_ratio_mod

use um_types,          only: real_64, real_32
use water_tracers_mod, only: l_wtrac_no_neg_ratio, min_q_ratio, wtrac_info

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim


implicit none

! Description:
!   Subroutine and function to calculate the ratio of water tracer to
!   normal water.
!   If there is no water, the ratio is set to a standard value
!
!   **There are two versions of the subroutine and function for different
!   input/output precisions.**
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

public :: wtrac_calc_ratio_fn
interface wtrac_calc_ratio_fn

  module procedure wtrac_calc_ratio_fn_real64,                                 &
                   wtrac_calc_ratio_fn_real32
end interface

integer(kind=jpim), private, parameter :: zhook_in  = 0
integer(kind=jpim), private, parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle


character(len=*), parameter, private :: ModuleName='WTRAC_CALC_RATIO_MOD'

contains

!===========================================================================

! Function for real_64 precision

real function wtrac_calc_ratio_fn_real64(i_wt, q_wtrac, q, l_rate_in)

implicit none

integer, parameter :: field_kind = real_64


#include "wtrac_calc_ratio_fn.h"


wtrac_calc_ratio_fn_real64 = q_ratio

end function wtrac_calc_ratio_fn_real64

!===========================================================================

! Function for real_32 precision

real function wtrac_calc_ratio_fn_real32(i_wt, q_wtrac, q, l_rate_in)

implicit none

integer, parameter :: field_kind = real_32


#include "wtrac_calc_ratio_fn.h"


wtrac_calc_ratio_fn_real32 = q_ratio

end function wtrac_calc_ratio_fn_real32

end module wtrac_calc_ratio_mod
