! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Constants used by convective diagnosis routines.
!
module cv_diag_param_mod

!-----------------------------------------------------------------------
! Description:
!   Module containing parameters used by the convective diagnosis routines
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN90
!  This code is written to UMDP3  v8 programming standards.
!
!-----------------------------------------------------------------------

use um_types, only: r_bl

implicit none

! Parameters used in determining the thetav perturbation given to the parcel
! for the undilute ascent.

real(kind=r_bl), parameter ::                                                  &
  a_plume = 0.2_r_bl         & ! Minimum initial parcel dthetav
 ,b_plume = 3.26_r_bl        & ! Used in initial parcel perturbation cal.
 ,max_t_grad = 1.0e-3_r_bl     ! Used in initial parcel perturbation cal.


! Parameter used in cumulus testing to see if a cloud layer present

real(kind=r_bl), parameter ::                                                  &
  sc_cftol   = 0.1_r_bl        ! Cloud fraction required for a cloud layer
                          ! to be diagnosed.

! Parameters for the calculation of T at the lifting condensation level.
! See Bolton 1980 Mon Wea Rev 108, P1046-1053 for details of empirical
! relationship

real(kind=r_bl), parameter ::                                                  &
  a_bolton = 55.0_r_bl                                                         &
, b_bolton = 2840.0_r_bl                                                       &
, c_bolton = 3.5_r_bl                                                          &
, d_bolton = 4.805_r_bl

end module cv_diag_param_mod
