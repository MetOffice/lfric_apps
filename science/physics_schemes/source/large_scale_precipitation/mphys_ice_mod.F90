! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module mphys_ice_mod

! Description:
! Holds Ice constants required by the large-scale
! precipitation scheme
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

use mphys_constants_mod, only: mprog_min

use um_types, only: real_umphys

implicit none

!-----------------------------------------------------------------
! Nucleation of ice
!-----------------------------------------------------------------

! Note that the assimilation scheme uses temperature thresholds
! in its calculation of qsat.

! Nucleation mass
real(kind=real_umphys), parameter :: m0    =  1.0e-12

! Maximum temperature for homogeneous nucleation
real(kind=real_umphys), parameter :: thomo = -40.0

!  1.0/Scaling quantity for ice in crystals
real(kind=real_umphys), parameter :: qcf0  = 1.0e4
                                    ! This is an inverse quantity

! Minimum allowed QCF after microphysics
real(kind=real_umphys), parameter:: qcfmin  = mprog_min

! 1/scaling temperature in aggregate fraction calculation
real(kind=real_umphys), parameter :: t_scaling = 0.0384

!  Minimum temperature limit in calculation  of N0 for ice (deg C)
real(kind=real_umphys), parameter :: t_agg_min = -45.0

!-----------------------------------------------------------------
! Turbulent generation of mixed phase cloud
!-----------------------------------------------------------------
real(kind=real_umphys), parameter ::  rhoi = 100.0 ! Ice density [kg m-3]

end module mphys_ice_mod
