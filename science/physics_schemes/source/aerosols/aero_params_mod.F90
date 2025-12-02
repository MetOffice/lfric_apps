! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Description:
!   Set median radii and standard deviation (sigma) of different aerosol types.
!
! Method:
!   Declare and set the two parameters of the logarithmic size distribution
!   (i.e. median radius, also called geometric mean radius, and geometric
!   standard deviation) for different aerosol types. Note that this is
!   valid for dry aerosols. All values here are exactly the same ones
!   as those in the standalone program calc_pm_params.f90.
!   If needed this module can be extended in the future to include other
!   parameters such as density.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: aerosols
!
! Code description:
!   Language: FORTRAN 2003.
!   This code is written to UMDP3 programming standards.
! --------------------------------------------------------------------------

module aero_params_mod

use um_types, only: real_umphys

implicit none

! Variables for ammonium sulphate (Aitken and accumulation)
real(kind=real_umphys), parameter :: r_bar_su_ait  = 0.0065e-6
                                              ! all median radii are in metres
real(kind=real_umphys), parameter :: sigma_su_ait  = 1.3
real(kind=real_umphys), parameter :: r_bar_su_acc  = 0.095e-6
real(kind=real_umphys), parameter :: sigma_su_acc  = 1.4

! Variables for black carbon (fresh and aged)
real(kind=real_umphys), parameter :: r_bar_bc_fr   = 0.04e-6
real(kind=real_umphys), parameter :: sigma_bc_fr   = 2.0
real(kind=real_umphys), parameter :: r_bar_bc_ag   = 0.04e-6
real(kind=real_umphys), parameter :: sigma_bc_ag   = 2.0

! Variables for biomass burning aerosol (fresh and aged)
real(kind=real_umphys), parameter :: r_bar_bb_fr   = 0.1e-6
real(kind=real_umphys), parameter :: sigma_bb_fr   = 1.3
real(kind=real_umphys), parameter :: r_bar_bb_ag   = 0.12e-6
real(kind=real_umphys), parameter :: sigma_bb_ag   = 1.3

! Variables for OCFF (fresh and aged)
real(kind=real_umphys), parameter :: r_bar_ocff_fr = 0.1e-6
real(kind=real_umphys), parameter :: sigma_ocff_fr = 1.3
real(kind=real_umphys), parameter :: r_bar_ocff_ag = 0.12e-6
real(kind=real_umphys), parameter :: sigma_ocff_ag = 1.3

! Variables for biogenic SOA
real(kind=real_umphys), parameter :: r_bar_soa     = 0.095e-6
real(kind=real_umphys), parameter :: sigma_soa     = 1.5

! Variables for sea salt aerosol (film and jet)
real(kind=real_umphys), parameter :: r_bar_ss_fi   = 0.1e-6
real(kind=real_umphys), parameter :: sigma_ss_fi   = 1.9
real(kind=real_umphys), parameter :: r_bar_ss_je   = 1.0e-6
real(kind=real_umphys), parameter :: sigma_ss_je   = 2.0

! Variables for ammonium nitrate (accumulation mode)
real(kind=real_umphys), parameter :: r_bar_ni_acc  = 0.095e-6
real(kind=real_umphys), parameter :: sigma_ni_acc  = 1.4

end module aero_params_mod
