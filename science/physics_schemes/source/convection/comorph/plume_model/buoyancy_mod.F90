! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module buoyancy_mod

implicit none

! Module contains addresses for a super-array
! which stores the parcel buoyancy at each of several heights
! within a given level-step, and the relevant heights,
! in height order.  This can be used to calculate detrainment
! and/or vertical momentum source terms, with accurate
! interpolation between the different levels
! (and the height where the parcel first hits saturation)

! The first field this array stores is height;
! don't need to declare this, as can just use i_height from
! grid_type_mod which is similarly set to 1.

! Second field this array stores is the parcel mean buoyancy
integer, parameter :: i_mean_buoy = 2

! Third field is the buoyancy of the parcel core
integer, parameter :: i_core_buoy = 3

! Possible different sub-level heights where these fields
! maybe stored:
! a) Previous model-level interface
! b) Centre of current level k
! c) Next model-level interface
! d) Height where parcel mean properties first hit saturation
! e) Height where parcel core first hit saturation

! Note, these get arranged in height order in the array,
! so they don't have fixed addresses (address depends on
! whether the parcel mean and core hits saturation, and
! if so where the saturation heights d,e are relative to the
! model-grid heights a,b,c.

! Address of prev is fixed, so declare here
integer, parameter :: i_prev = 1


end module buoyancy_mod
