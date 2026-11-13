! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Module mym_const_mod----------------------------------------------

!  Purpose: To define constants including the closure constants used
!           in the MY model.

!  Programming standard : UMDP 3

!  Documentation: UMDP 025

!  This code is based on the code provided by the authors who wrote
!  the following papers.
!   * Nakanishi, M. and H. Niino, 2009: Development of an improved
!      turbulence closure model for the atmospheric boundary layer.
!      J. Meteor. Soc. Japan, 87, 895-912.
!   * Nakanishi, M. and H. Niino, 2006: An improved Mellor-Yamada
!      Level-3 model: Its numerical stability and application to
!      a regional prediction of advection fog.
!      Boundary-Layer Meteor., 119, 397-407.
!   * Nakanishi, M. and H. Niino, 2004: An improved Mellor-Yamada
!      Level-3 model with condensation physics: Its design and
!      verification.
!      Boundary-Layer Meteor., 112, 1-31.
!   * Nakanishi, M., 2001: Improvement of the Mellor-Yamada
!      turbulence closure model based on large-eddy simulation data.
!      Boundary-Layer Meteor., 99, 349-378.
!   The web site publicising their code:
!    http://www.nda.ac.jp/~naka/MYNN/index.html

!---------------------------------------------------------------------
module mym_const_mod
use um_types, only: r_bl

implicit none
save

 ! For the meaing of the variables, see the papers above.
 ! N2001: Nakanishi, M., 2001
 ! NN2004: Nakanishi, M. and H. Niino, 2004
 ! NN2006: Nakanishi, M. and H. Niino, 2006

 ! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
real(kind=r_bl) ::                                                      &
   g1,                                                                         &
                  ! gamma_1 = 1/3 - 2A_1 / B_1
   g2,                                                                         &
                  ! gamma_2 defined at Eq. (B4) in N2001
   a1,                                                                         &
                  ! closure constant A_1 defined at Eq.(7) in N2001
   a2,                                                                         &
                  ! closure constant A_2 defined at Eq.(8) in N2001
   b1,                                                                         &
                  ! closure constant B_1 defined at Eq.(5) in N2001
   b2,                                                                         &
                  ! closure constant B_2 defined at Eq.(6) in N2001
   c1,                                                                         &
                  ! closure constant C_1 defined at Eq.(7) in N2001
   c2,                                                                         &
                  ! closure constant C_2 defined at Eq.(7) in N2001
   c3,                                                                         &
                  ! closure constant C_3 defined at Eq.(7) in N2001
   c4,                                                                         &
                  ! closure constant C_1 defined at Eq.(7) in N2001
   c5,                                                                         &
                  ! closure constant C_1 defined at Eq.(7) in N2001
   a1_2,                                                                       &
                  ! A_1 / A_2
   pr,                                                                         &
                  ! Prandtl Number
   rfc,                                                                        &
                  ! Critical flux Richardson number
   f1,                                                                         &
                  ! F_1 defined at Eq. (B4) in N2001
   f2,                                                                         &
                  ! F_2 defined at Eq. (B4) in N2001
   rf1,                                                                        &
                  ! R_f1 defined at Eq. (B4) in N2001
   rf2,                                                                        &
                  ! R_f2 defined at Eq. (B4) in N2001
   smc,                                                                        &
                  ! A_1 F_1 / (A_2 F_2) appeared at Eq. (B2) in N2001
   shc,                                                                        &
                  ! 3 A_2 (gamma_1 + gamma_2) appeared at Eq. (B2)
                  ! in N2001
   ri1,                                                                        &
                  ! Ri_1 defined at Eq. (B7) in N2001
   ri2,                                                                        &
                  ! Ri_2 defined at Eq. (B7) in N2001
   ri3,                                                                        &
                  ! Ri_3 defined at Eq. (B7) in N2001
   ri4,                                                                        &
                  ! Ri_4 defined at Eq. (B7) in N2001
   cc2,                                                                        &
                  ! 1-C_2
   cc3,                                                                        &
                  ! 1-C_3
   e1c,                                                                        &
                  ! constant appeared in phi1 at Eq.(3a) in NN2006
   e2c,                                                                        &
                  ! constant appeared in phi2 at Eq.(3b) in NN2006
   e3c,                                                                        &
                  ! constant appeared in phi3 and phi3'
                  ! at Eq. (3c,d) in NN2006
   e4c,                                                                        &
                  ! constant appeared in phi4 and phi4'
                  ! at Eq. (3e,f) in NN2006
   e5c,                                                                        &
                  ! constant appeared in D and D'
                  ! at Eq.(2a,b) in NN2006
   my_alpha1,                                                                  &
                  ! alpha_1 defined at Eq.(40) in N2001
   my_alpha2,                                                                  &
                  ! alpha_2 defined at Eq.(41) in N2001
   my_alpha3,                                                                  &
                  ! alpha_3 defined at Eq.(41) in N2001
   my_alpha4,                                                                  &
                  ! alpha_4 defined at Eq.(39) in N2001
   elt_min,                                                                    &
                  ! lower limit of elt
   one_third,                                                                  &
                  ! 1.0 / 3.0
   two_thirds,                                                                 &
                  ! 2.0 / 3.0
   coef_trbvar_diff_tke,                                                       &
                  ! factor of the diffusion for TKE to the one
                  ! for momentum
   coef_trbvar_diff,                                                           &
                  ! factor of the diffusion for the other prognostic
                  ! variables appeared in MY  to the one
                  ! for momentum
   qke_max,                                                                    &
                  ! upper limit for qke (twice of TKE) for safety.
   e_trb_max
                  ! upper limit for e_trb for safety.
end module mym_const_mod
