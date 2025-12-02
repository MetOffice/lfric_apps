! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lsp_dif_mod

! Description:
! The apb and tw parameters represent diffusional growth constants
! and wet bulb temperature parameters for use in LSP.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

use um_types, only: real_umphys

implicit none


! Values of reference variables
real(kind=real_umphys),parameter:: air_density0=1.0             ! kg m-3
real(kind=real_umphys),parameter:: air_viscosity0=1.717e-5      ! kg m-1 s-1
real(kind=real_umphys),parameter:: air_conductivity0=2.40e-2    ! J m-1 s-1 K-1
real(kind=real_umphys),parameter:: air_diffusivity0=2.21e-5     ! m2 s-1
real(kind=real_umphys),parameter:: air_pressure0=1.0e5          ! Pa

! Values of constants used in correction terms
real(kind=real_umphys), parameter :: tcor1 = 393.0
real(kind=real_umphys), parameter :: tcor2 = 120.0
real(kind=real_umphys), parameter :: cpwr  = 1.5

! Values of diffusional growth parameters (set in lspcon)
! Terms in deposition and sublimation
real(kind=real_umphys) :: apb1 ! =(lc+lf)**2 * repsilon /(r*air_conductivity0)
real(kind=real_umphys) :: apb2 ! =(lc+lf) / air_conductivity0
real(kind=real_umphys) :: apb3 ! =r/(repsilon*air_pressure0*air_diffusivity0)
! Terms in evap of melting snow and rain
real(kind=real_umphys) :: apb4 ! =lc**2*repsilon/(r*air_conductivity0)
real(kind=real_umphys) :: apb5 ! =lc /air_conductivity0
real(kind=real_umphys) :: apb6 ! =r/(repsilon*air_pressure0*air_diffusivity0)

! Values of numerical approximation to wet bulb temperature
! Numerical fit to wet bulb temperature
real(kind=real_umphys),parameter:: tw1=1329.31
real(kind=real_umphys),parameter:: tw2=0.0074615
real(kind=real_umphys),parameter:: tw3=0.85e5
! Numerical fit to wet bulb temperature
real(kind=real_umphys),parameter:: tw4=40.637
real(kind=real_umphys),parameter:: tw5=275.0

! Ventilation parameters
real(kind=real_umphys),parameter:: sc=0.6
! f(v)  =  vent_ice1 + vent_ice2  Sc**(1/3) * Re**(1/2)
real(kind=real_umphys),parameter:: vent_ice1=0.65
real(kind=real_umphys),parameter:: vent_ice2=0.44
! f(v)  =  vent_rain1 + vent_rain2  Sc**(1/3) * Re**(1/2)
real(kind=real_umphys),parameter:: vent_rain1=0.78
real(kind=real_umphys),parameter:: vent_rain2=0.31

end module lsp_dif_mod
