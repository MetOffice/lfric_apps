! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: gravity_wave_drag
module c_gwave_mod

use um_types, only: real_umphys

implicit none
!
!  Description: This module defines the constants for the 4A and 5A
!               versions of the orographic gravity wave drag scheme. These are
!               tuneable parameters but are unlikely to be changed.

! Orographic drag scheme parameters:

! 4A scheme only
      ! Switch to determine which wave saturation test is used
integer,parameter :: Amplitude_saturation = 1
integer,parameter :: Stress_saturation    = 0

! 5A scheme only
      ! Fixed value for group velocity angle (used when l_nonhydro true
      ! and l_dynbeta false)
real(kind=real_umphys),parameter     ::  beta_fix     = 100.0
! Fraction of vert wavelength to deposit acceleration over
! (when l_smooth true)
real(kind=real_umphys),parameter     ::  frac_wl      = 1.0
! Minimum value for local vert wavelength (U/N) (m)
real(kind=real_umphys),parameter     ::  lambdaz_min  = 250.0
! Maximum value for local vert wavelength (U/N) (m)
real(kind=real_umphys),parameter     ::  lambdaz_max  = 3000.0
! Threshold for brunt-vaisala frequency squared (ensq)
! to define a neutral layer (used in gw_block and gw_wave)
real(kind=real_umphys), parameter    ::  nsq_neutral  = 1.0e-6
! Threshold for % agreement for Zav convergence in gw_block
real(kind=real_umphys), parameter    ::  zav_converge = 0.05
! Num of iterations to allow Zav to converge in gw_block
integer, parameter ::  zav_iterate  = 5
!

end module c_gwave_mod
