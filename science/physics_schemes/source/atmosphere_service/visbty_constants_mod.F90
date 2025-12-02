! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!   A module containing constants used in visibility diagnosis

! Description:
!   This module contains declarations for constants used to diagnose
!   horizontal visibility.

!   Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_visibility

! Code description:
!   This code is written to UMDP3 standards.

module visbty_constants_mod

use conversions_mod, only: pi
use  missing_data_mod, only: rmdi

use um_types, only: real_umphys

implicit none

! General constants (these should be promoted to a more general module)

real(kind=real_umphys), parameter :: onethird   = 1.0/3.0
real(kind=real_umphys), parameter :: fourthirds = 4*onethird

! Information about the diagnostics themselves

integer, parameter :: n_vis_thresh = 2
real(kind=real_umphys), parameter    :: visfog = 1000.0
                                           ! Visibility defining fog
real(kind=real_umphys), parameter    :: vismist = 5000.0
                                           ! Visibility defining mist
real(kind=real_umphys), parameter    ::                                        &
                               vis_thresh(n_vis_thresh)=[ visfog, vismist ]
real(kind=real_umphys)               :: calc_prob_of_vis = rmdi
                  ! tunable parameter: the cumulative prob value at
                  ! which vis is estimated, set in RUN_BL namelist

! Parameters used in visibility calculations

real(kind=real_umphys), parameter :: liminalcontrast = 0.05
real(kind=real_umphys), parameter :: lnliminalcontrast = -2.99573
                                      ! Natural log of Liminal contrast
real(kind=real_umphys), parameter :: recipvisair = 1.0e-5
                                        ! Reciprocal of the clean air
                                      ! visibility (100km)

! Microphysical parameters used in visibility calculations

real(kind=real_umphys), parameter :: n0 = 200.0e7
                                        ! Standard number density
                                        ! of murk (/m3)
real(kind=real_umphys), parameter :: b0= 0.14             ! Activation parameter
real(kind=real_umphys), parameter :: radius0 = 0.11e-6
                                        ! Radius of standard murk
                                        ! particle (m)
real(kind=real_umphys), parameter :: rho_aerosol = 1700.0
                                        ! Density of murk (kg/m3)
real(kind=real_umphys), parameter :: rho_air = 1.0
                                        ! Density of air (kg/m3)

real(kind=real_umphys), parameter :: beta0  = 1.5 * pi
                                        ! Scattering coefficient
                                        ! normalisation
real(kind=real_umphys), parameter :: a0 = 1.2e-9
                                        ! Constant involving surface
                                        ! energy of water
real(kind=real_umphys), parameter ::                                           &
                   m0 = fourthirds * pi * radius0 * radius0 * radius0          &
                        * (rho_aerosol/rho_air) * n0
                                      ! Standard aerosol mass mixing
                                      ! ratio (kg/kg)
real(kind=real_umphys), parameter :: power = 1.0/6.0
                                        ! Murk particle radius/mass
                                        ! loading power
real(kind=real_umphys), parameter :: visfactor = -lnliminalcontrast / beta0
                                      ! transformation to visibility
                                      ! ( = ln(liminal contrast) / Beta0 )

real(kind=real_umphys), parameter :: aero0 = 0.1
                                        ! Minimum allowed murk aerosol
real(kind=real_umphys), parameter :: aeromax = 200.0
                                        ! Maximum allowed murk aerosol

end module visbty_constants_mod
