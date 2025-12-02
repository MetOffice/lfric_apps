! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer

! PURPOSE: To calculate buoyancy parameters on p,T,q-levels

! CODE DESCRIPTION:
!   LANGUAGE: FORTRAN 95
!   THIS CODE is WRITTEN to UMDP 3 PROGRAMMING STANDARDS.

module buoy_tq_mod

use um_types, only: real_64, real_32

implicit none

private
public :: buoy_tq

interface buoy_tq
  module procedure buoy_tq_64b, buoy_tq_32b
end interface

character(len=*), parameter, private :: ModuleName = 'BUOY_TQ_MOD'
contains

subroutine buoy_tq_64b (                                                       &
! in dimensions/logicals
 bl_levels,                                                                    &
! in fields
 p,t,q,qcf,qcl,cf_bulk,                                                        &
! out fields
 bt,bq,bt_cld,bq_cld,bt_gb,bq_gb,a_qs,a_dqsdt,dqsdt                            &
 )

use atm_fields_bounds_mod, only: tdims, tdims_l
use bl_option_mod, only: l_noice_in_turb
use gen_phys_inputs_mod, only: l_mr_physics
use planet_constants_mod, only: r, repsilon, c_virtual, etar, lcrcp, ls, lsrcp
use water_constants_mod, only: lc, tm
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use qsat_mod, only: qsat, qsat_mix, qsat_wat, qsat_wat_mix

implicit none

character(len=*), parameter :: RoutineName='BUOY_TQ_64B'
integer, parameter :: prec = real_64
#include "buoy_tq.h"

return
end subroutine buoy_tq_64b

subroutine buoy_tq_32b (                                                       &
! in dimensions/logicals
 bl_levels,                                                                    &
! in fields
 p,t,q,qcf,qcl,cf_bulk,                                                        &
! out fields
 bt,bq,bt_cld,bq_cld,bt_gb,bq_gb,a_qs,a_dqsdt,dqsdt                            &
 )

use atm_fields_bounds_mod, only: tdims, tdims_l
use bl_option_mod, only: l_noice_in_turb
use gen_phys_inputs_mod, only: l_mr_physics
use planet_constants_mod, only: r => r_32b, repsilon => repsilon_32b,          &
     c_virtual => c_virtual_32b, etar => etar_32b, lcrcp => lcrcp_32b,         &
     ls => ls_32b, lsrcp => lsrcp_32b
use water_constants_mod, only: lc => lc_32b, tm => tm_32b
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use qsat_mod, only: qsat, qsat_mix, qsat_wat, qsat_wat_mix

implicit none

character(len=*), parameter :: RoutineName='BUOY_TQ_32B'
integer, parameter :: prec = real_32
#include "buoy_tq.h"

return
end subroutine buoy_tq_32b

end module buoy_tq_mod
