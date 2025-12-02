! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Global module for storing parameters used in parametrization of
! fractional standard deviation (FSD) of subgrid water content.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: radiation_control
!
!----------------------------------------------------------------------

module fsd_parameters_mod

use um_types, only: real_umphys

implicit none
save

! Parameters for FSD parametrization

real(kind=real_umphys), allocatable :: f_arr(:,:,:,:)
!   Array of regime and layer thickness dependent
!   parameters required for parametrization of FSD.

real(kind=real_umphys), allocatable :: f_arr_c(:,:)
!   Array of regime and layer thickness dependent
!   parameters required for parametrization of FSD.
!   Gathered version for use in microphysics.

real(kind=real_umphys) :: f_cons(3)
!   Constant parameters required for parametrization of FSD

real(kind=real_umphys) :: fsd_eff_lam
                    ! effective delta_lambda for fsd parametrization
real(kind=real_umphys) :: fsd_eff_phi
                    ! effective delta_phi for fsd parametrization

integer, parameter :: ip_fsd_constant     = 0
!   flag to use constant value of FSD/scaling factor
integer, parameter :: ip_fsd_param        = 1
!   flag to use FSD parametrisation described in Hill et al (2012)
!   doi: 10.1002/qj.1893
integer, parameter :: ip_fsd_regime        = 2
!   flag to use regime dependent FSD parametrisation
integer, parameter :: ip_fsd_regime_no_sh  = 3
!   flag to use regime dependent FSD parametrisation
!   doesn't include shallow convection
integer, parameter :: ip_fsd_boutle        = 4
!   flag to use FSD parametrisation described in Boutle et al (2013)
!   doi: 10.1002/qj.2140
integer, parameter :: ip_fsd_regime_smooth = 5
!   flag to use regime dependent FSD parametrisation on smoothed
!   cca field
integer, parameter :: ip_fsd_regime_smooth_no_sh  = 6
!   flag to use regime dependent FSD parametrisation on smoothed
!   cca field. Doesn't include shallow convection
integer, parameter :: ip_fsd_regime_cca    = 7
!   flag to use regime dependent FSD parametrisation based on cca only
!   not shallow, mid, deep
integer, parameter :: ip_fsd_regime_cca_v8 = 8
!   flag to use more extreme regime dependent FSD
!   parametrisation based on local cca > 0.02

end module fsd_parameters_mod
