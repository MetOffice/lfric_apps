! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module surface_flux_mod

! Purpose:
!   Contains information on surface flux forcing  and SST forcing for the
!   idealised model.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

use missing_data_mod, only: rmdi, imdi

implicit none

save

! Parameters
! ======================
! Maximum size of arrays that can be read through namelist
integer, parameter :: num_surf_max  = 100

integer, parameter :: zero_flux       = 1
integer, parameter :: diurnal_flux    = 2
integer, parameter :: constant_flux   = 3
integer, parameter :: hot_spot        = 4
integer, parameter :: time_varying    = 5

integer, parameter :: sst_fixed       = 1
integer, parameter :: sst_varying     = 2

! Variables from idealised namelist
! ------------------------------------

integer :: IdlSurfFluxSeaOption = imdi

real    :: IdlSurfFluxseaParams(4) = rmdi

! Time varying surface fluxes
integer :: num_surface_flux_times   = imdi

real    :: surface_flux_time(num_surf_max) = rmdi
real    :: sh_flux(num_surf_max)           = rmdi
real    :: lh_flux(num_surf_max)           = rmdi

! SST
integer :: IdlSSTOption = imdi

! Time varying SST
integer :: num_sst_times   = imdi

real    :: sst_time(num_surf_max) = rmdi
real    :: sst_data(num_surf_max) = rmdi

end module surface_flux_mod
