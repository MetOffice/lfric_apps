! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Purpose: Holds constants required for the electric scheme

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Electric

module electric_constants_mod

use conversions_mod, only: zerodegc

use um_types, only: real_umphys

implicit none

! Conversion terms
! minutes in a second
real(kind=real_umphys), parameter :: min2sec = 1.0_real_umphys/60.0_real_umphys

! five minutes into a second
real(kind=real_umphys), parameter :: fivemin2sec = min2sec / 5.0_real_umphys

! To go from kg to g.
real(kind=real_umphys), parameter :: kg_to_g = 1000.0_real_umphys

! Graupel Water Path method terms
! Minus 5 level
real(kind=real_umphys), parameter :: minus5 = zerodegc - 5.0_real_umphys

! McCaul et al (2009) parameters
real(kind=real_umphys), parameter :: minus15 = zerodegc - 15.0_real_umphys
                                             ! Minus 15 C level for McCaul
real(kind=real_umphys), parameter :: mccaul_r1 = 0.95_real_umphys
                                    ! Factor for inclusion of flash1 in McCaul
                                    ! scheme
real(kind=real_umphys), parameter :: mccaul_r2 = 0.05_real_umphys
                                    ! Factor for inclusion of flash2 in McCaul
                                    ! scheme
                                    ! Note: mccaul_r1 + mccaul_r2 should always
                                    !       equal 1.0 exactly.

! Constants for mixed-precision use.
real(kind=real_umphys), parameter :: zero = 0.0_real_umphys
real(kind=real_umphys), parameter :: one  = 1.0_real_umphys

end module electric_constants_mod

