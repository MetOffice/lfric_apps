! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module subregion_mod

use comorph_constants_mod, only: name_length

implicit none

! Module storing indices to be used in a super-array of
! properties of different sub-regions of the grid-box, as
! defined below...

! Integer indicators for each of the 4 sub-regions
! described below...
integer, parameter :: n_regions = 4
integer, parameter :: i_dry = 1  ! Clear-sky
integer, parameter :: i_liq = 2  ! Liquid-cloud only
integer, parameter :: i_mph = 3  ! Mixed-phase cloud
integer, parameter :: i_icr = 4  ! Ice / rain (no liq cloud)

! Names of the regions stored as character strings, used
! for naming diagnostics or error reporting
character(len=name_length), parameter :: region_names(n_regions)               &
                                = [ "dry", "liq", "mph", "icr" ]

! The 4 sub-regions calculated are:
! - clear-sky         (dry)
! - liquid-only cloud (liq)
! - mixed-phase cloud (mph)
! - ice and rain only (icr)
! These 4 regions are all assumed to have the same virtual
! temperature, but differing water vapour.  The liquid and
! mixed-phase cloud regions are assumed to be saturated w.r.t.
! liquid water, whereas the clear-sky and ice/rain regions
! are subsaturated.  Each condensed water species is assumed
! to be confined within certain regions:
! - q_cl    (liq,mph    )
! - q_rain  (liq,mph,icr)
! - q_cf    (    mph,icr)
! - q_snow  (    mph,icr)
! - q_graup (    mph,icr)
! The local values of T, q_vap and condensed water mixing ratios
! are calculated in each region so-as to satisfy the above
! assumptions, and the conservation laws applying to the
! fractional area f of each region:
! f_dry       + f_liq       + f_mph       + f_icr       = 1
! f_dry T_dry + f_liq T_liq + f_mph T_mph + f_icr T_icr = T_env
! f_dry q_dry + f_liq q_liq + f_mph q_mph + f_icr q_icr = q_env

end module subregion_mod
