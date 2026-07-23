! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module sublevs_mod

implicit none

save

! Module contains addresses for a super-array which stores the
! parcel buoyancy and other properties at several intermediate heights
! within a given level-step, in height order.
! This can be used to calculate detrainment
! and/or vertical momentum source terms, with accurate
! interpolation between the different levels
! (e.g. to the height where the parcel first hits saturation)

! Flag for whether this module has been initialised
logical :: l_init_sublevs_mod = .false.

! Total number of fields stored
integer :: n_sublev_vars = 0

! Fields which maybe stored in the sub-levels super-array:
! (note these are named starting with j rather than i, to avoid
!  variable-name clashes with the addresses of some of the same
!  variables in other super-arrays, e.g. i_height in grid super-arrays,
!  and i_massflux_d in parcel super-arrays)

! Height of each sub-level
integer :: j_height = 0

! Mass-flux at sub-level height
integer :: j_massflux_d = 0

! Environment virtual temperature at sub-level heights
integer :: j_env_tv = 0
! Increment to env Tv due to compensating subsidence
integer :: j_delta_tv = 0

! Environment vertical velocity excess at sub-level heights
integer :: j_env_w = 0

! Buoyancy of the parcel mean and core
integer :: j_mean_buoy = 0
integer :: j_core_buoy = 0

! Vertical velocity excess of the parcel mean and core
integer :: j_mean_wex = 0
integer :: j_core_wex = 0

! Max number of sub-level heights at which they maybe stored
integer :: max_sublevs = 0

! Possible different sub-level heights where these fields
! maybe stored:
! a) Previous model-level interface
! b) Next model-level interface
! c) Height where parcel mean properties first hit saturation
! d) Height where parcel core first hit saturation

! Note, these get arranged in height order in the array,
! so they don't have fixed addresses (address depends on
! whether the parcel mean and core hits saturation, and
! if so where the saturation heights c,d are relative to the
! model-grid heights a,b,c.

! Address of prev is fixed, so declare here
integer, parameter :: i_prev = 1


contains

! Subroutine to set the addreses of the fields in the sub-levels super-array
subroutine sublevs_set_addresses()

use comorph_constants_mod, only: l_par_core

implicit none

! Set flag to indicate that we don't need to call this routine again!
l_init_sublevs_mod = .true.

! Height, mass-flux, saturated fraction, env Tv and w and
! subsidence Tv increment are always stored
j_height     = 1
j_massflux_d = 2
j_env_tv     = 3
j_delta_tv   = 4
j_env_w      = 5

n_sublev_vars = 5

! Buoyancies and w-excesses always stored at the end
j_mean_buoy = n_sublev_vars + 1
j_mean_wex  = n_sublev_vars + 2
n_sublev_vars = n_sublev_vars + 2

! Core buoyancy and w-excess are stored if using parcel core properties
if ( l_par_core ) then
  j_core_buoy = n_sublev_vars + 1
  j_core_wex  = n_sublev_vars + 2
  n_sublev_vars = n_sublev_vars + 2
end if

! Allowed heights are previous level, next level, saturation height,
! and edge saturation height if the parcel core is used
if ( l_par_core ) then
  max_sublevs = 4
else
  max_sublevs = 3
end if

return
end subroutine sublevs_set_addresses


end module sublevs_mod
