! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Global data module for variables concerned with solar incidence.

module solinc_data

use um_types, only: real_umphys

implicit none
save

! Description:
!   Global data necessary for calculating the angle of solar incidence
!   on sloping terrain.
!
! Method:
!   Provides global data.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: radiation_control
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v7.4 programming standards.
!
! Declarations:
!
! Global variables (#include statements etc):

real(kind=real_umphys), allocatable :: slope_aspect(:,:), slope_angle(:,:)
real(kind=real_umphys), allocatable :: f_orog(:,:), orog_corr(:,:)
real(kind=real_umphys), allocatable :: horiz_ang(:,:,:,:), horiz_aspect(:,:,:)
real(kind=real_umphys), allocatable :: sky(:,:)
logical :: l_orog = .false.
logical :: l_skyview = .false.
integer :: horiz_limit = 30
integer :: n_horiz_layer = 1
integer, parameter :: n_horiz_ang = 16

! slope_aspect: The direction faced by the mean slope - i.e. the
!               bearing of the slope normal projected on the surface
!               measured in radians clockwise from true north.
!
! slope_angle:  Angle of the mean slope measured in radians from
!               the horizontal.
!
! orog_corr:    correction factor for the direct solar flux
!               reaching the surface for sloping terrain.
!
! f_orog:       The extra direct solar flux at the surface due to
!               the orography correction. This is used in the
!               correction to the sw_incs calculation and to
!               correct net_atm_flux.
!
! l_orog:       model switch for orography scheme
!
! horiz_ang:    Angle in radians measured from the local zenith to
!               the obscuring terrain.
!
! horiz_aspect: The local bearing for each horizon angle measured
!               in radians clockwise from true north.
!
! sky:          Sky-view correction factor for net surface LW.
!
! l_skyview:    Model switch for skyview scheme.
!
! horiz_limit:  Number of grid-lengths to furthest point considered
!               in horizon calculation.
!
! n_horiz_layer: Number of layers where horizon angles are calculated.
!
! n_horiz_ang:  Number of horizon angles calculated.
!
!- End of header

end module solinc_data
