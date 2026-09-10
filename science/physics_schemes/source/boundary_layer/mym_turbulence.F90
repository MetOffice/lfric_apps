! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Purpose: To calculate the diffusion coefficients and counter
!           gradient term for momentum, heat and moisture, and
!           integrate the prognostic variables appearing
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

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!---------------------------------------------------------------------
module mym_turbulence_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'MYM_TURBULENCE_MOD'

contains

subroutine mym_turbulence(                                                     &
! IN levels/switches
      bl_levels, levflag, BL_diag,                                             &
! IN fields
      z_uv, z_tq,                                                              &
      vq, vt, gtr, fqw, ftl, wb_ng,                                            &
      dbdz, dtldz, dqwdz, dvdzm, dudz, dvdz, delta_smag,                       &
      r_mosurf, u_s, fb_surf, pmz, phh,                                        &
! INOUT fields
      qke, tsq, qsq, cov, dfm, dfh,                                            &
! OUT fields
      dfu_cg, dfv_cg, dft_cg, dfq_cg, leonard_kl_tke)

use atm_fields_bounds_mod, only: tdims, pdims, tdims_l, tdims_s
use bl_diags_mod, only: strnewbldiag
use conversions_mod, only: pi
use mym_const_mod, only: e1c,e2c,e3c,e4c,e5c,a1,a2,c1,c5,b2,qke_max,           &
      coef_trbvar_diff,coef_trbvar_diff_tke,two_thirds,a1_2,                   &
      b1,one_third,cc3
use mym_option_mod, only:                                                      &
      my_lowest_pd_surf, l_my_extra_level, my_z_extra_fact,                    &
      l_my_prod_adj, my_prod_adj_fact, tke_levels,                             &
      l_my_lowest_pd_surf_tqc, l_leonard_tke,                                  &
      deardorff, mymodel2, mymodel25, mymodel3

use model_domain_mod, only: model_type, mt_single_column

use parkind1, only: jprb, jpim
use planet_constants_mod, only: vkman
use yomhook, only: lhook, dr_hook

use mym_length_mod, only: mym_length
use mym_level2_mod, only: mym_level2
use mym_update_covariance_mod, only: mym_update_covariance
use mym_update_fields_mod, only: mym_update_fields

implicit none

! Intent IN Variables
integer, intent(in) ::                                                         &
   bl_levels,                                                                  &
                 ! Max. no. of "boundary" levels
   levflag
                 ! to indicate the level of the MY model
                 ! 0: Deardorff
                 ! 1: level 2
                 ! 2: level 2.5
                 ! 3: level 3

real(kind=r_bl), intent(in) ::                                          &
   z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                   &
        bl_levels+1),                                                          &
                 ! Z_UV(*,k) is height of u level k
   z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        bl_levels),                                                            &
                 ! Z_TQ(*,k) is height of theta level k.
   vq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                 ! A buoyancy param on theta level k-1
   vt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                 ! A buoyancy param on theta level k-1
   gtr(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                 ! g/thetav on theta level k-1
   fqw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                 ! Moisture flux between layers
                 !   (kg per square metre per sec).
                 ! FQW(,1) is total water flux
                 ! from surface, 'E'.
                 ! Defined on rho levels.
   ftl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                 ! FTL(,k) contains net turbulent
                 ! sensible heat flux into layer k
                 ! from below; so FTL(,1) is the
                 ! surface sensible heat, H. (W/m2)
                 ! Defined on rho levels.
   wb_ng(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         tke_levels),                                                          &
                 ! buoyancy flux related to the skewness
                 ! on theta k-1 levels
   dbdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        2:tke_levels),                                                         &
                 ! Buoyancy gradient across layer
                 ! interface interpolated to theta levels.
                 ! (:,:,k) represents the value on theta level k-1
   dtldz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         2:bl_levels),                                                         &
                 ! gradient of TL across layer
                 ! interface interpolated to theta levels.
                 ! (:,:,k) represents the value on theta level k-1
   dqwdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         2:bl_levels),                                                         &
                 ! gradient of QW across layer
                 ! interface interpolated to theta levels.
                 ! (:,:,k) represents the value on theta level k-1
   dvdzm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                  &
         2:bl_levels),                                                         &
                 ! Modulus of wind shear at theta levels.
                 ! (:,:,k) represents the value on theta level k-1
   dudz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        2:bl_levels),                                                          &
                 ! Gradient of u at theta levels.
                 !(:,:,k) represents the value on theta level k-1
   dvdz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        2:bl_levels),                                                          &
                 ! Gradient of v at theta levels.
                 !(:,:,k) rerserents the value on theta level k-1
   delta_smag(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),            &
                 ! IN delta_x used by Smagorinsky
   r_mosurf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),              &
                 ! reciprocal of Monin-Obukhov length
   u_s(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
                 ! Surface friction velocity
   fb_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),               &
                 ! Surface buoyancy flux
   pmz(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                   &
                 ! gradient function for momentum at surface
   phh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                 ! gradient function for scalars at surface

! Intent INOUT Variables
real(kind=r_bl), intent(in out) ::                                      &
   qke(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,            &
       bl_levels),                                                             &
                 ! twice of TKE (denoted to q**2) on theta level k-1
   tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       bl_levels),                                                             &
                 ! Self covariance of liquid potential temperature
                 ! (thetal'**2) defined on theta levels k-1
   qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       bl_levels),                                                             &
                 ! Self covariance of total water
                 ! (qw'**2) defined on theta levels k-1
   cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       bl_levels),                                                             &
                 ! Correlation between thetal and qw
                 ! (thetal'qw') defined on theta levels k-1
   dfm(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end,            &
       bl_levels),                                                             &
                 ! diffusion coefficient for momentum
                 ! on theta level k-1
   dfh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels)
                 ! diffusion coefficient for scalars
                 ! on theta level k-1

!  Declaration of BL diagnostics.
type (strnewbldiag), intent(in out) :: BL_diag

! Intent OUT Variables
real(kind=r_bl), intent(out) ::                                         &
   dfu_cg(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end,         &
          2:bl_levels),                                                        &
                 ! counter gradient term for u
                 ! on theta level k-1
   dfv_cg(tdims_s%i_start:tdims_s%i_end,tdims_s%j_start:tdims_s%j_end,         &
          2:bl_levels),                                                        &
                 ! counter gradient term for v
                 ! on theta level k-1
   dft_cg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          2:bl_levels),                                                        &
                 ! counter gradient term for TL
                 ! on theta level k-1
   dfq_cg(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          2:bl_levels)
                 ! counter gradient term for QW
                 ! on theta level k-1
real(kind=r_bl), intent(out) ::                                                &
         leonard_kl_tke( tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end,                            &
                         1:bl_levels, 2 )

! Local variables
! Scalar
integer ::                                                                     &
   i, j, k, k_start, k_start_cor
                 ! Loop indexes

real(kind=r_bl) ::                                                      &
   e1,                                                                         &
                 ! a variable denoted to E1 in the papers
   e3,                                                                         &
                 ! a variable denoted to E3 in the papers
   e4,                                                                         &
                 ! a variable denoted to E4 in the paper
   q2sq,                                                                       &
                 ! qke derived by level 2 scheme
   t2sq,                                                                       &
                 ! tsq derived by level 2 scheme
   r2sq,                                                                       &
                 ! qsq derived by level 2 scheme
   t3sq,                                                                       &
                 ! tsq derived by level 2.5 or 3 scheme
   r3sq,                                                                       &
                 ! qsq derived by level 2.5 or 3 scheme
   c3sq,                                                                       &
                 ! cov deribed by level 2.5 or 3 scheme
   wden,                                                                       &
                 ! work variable
   eden,                                                                       &
                 ! work variable
   reden,                                                                      &
                 ! reciprocal of eden
   e6c,                                                                        &
                 ! work variable
   coef,                                                                       &
                 ! work variable
   elq,                                                                        &
                 ! mixing length  times qkw appeared
                 ! in the production term of qke
   elh,                                                                        &
                 ! mixing length times qkw appeared
                 ! in the production terms of tsq, qsq and cov.
   phm,                                                                        &
                 ! work variable
   disp_coef,                                                                  &
   b1l,                                                                        &
                 ! B1 (closure constant) times mixing length
   b2l,                                                                        &
                 ! work variable
   clow,                                                                       &
                 ! lower limit for difference between cov in level 3
                 ! and level 2
   cupp
                 ! upper limit for difference between cov in level 3
                 ! and level 2

real(kind=r_bl) ::                                                      &
   gm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                 ! square of wind shear on theta level k-1
                 ! (a denominator of gradient Richardson number)
   gh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                 ! - buoyancy gradient on theta level k-1
                 ! (a numerator of gradient Richardson number)
   sm(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                 ! Non-dimensional diffusion coefficients for
                 ! momentum derived by level 2 scheme
                 ! defined on theta level k-1
   sh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                 ! Non-dimensional diffusion coefficients for
                 ! scalars derived by level 2 scheme
                 ! define on theta level k-1
   qkw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
        tke_levels),                                                           &
                 ! q=sqrt(qke) on theta level k-1
   elsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        tke_levels),                                                           &
                 ! square of mixing length
   gmel(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        tke_levels),                                                           &
                 ! GM times the mixing length squared
   ghel(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        tke_levels),                                                           &
                 ! GH times the mixing length squared
   qdiv(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        tke_levels),                                                           &
                 ! factor for flux correction: sqrt(q3sq/q2sq)
   gamv(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        tke_levels),                                                           &
                 ! counter gradient correction for production of qke
   gamv_coef(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,              &
             tke_levels),                                                      &
                 ! coefficient of (c3sq-c2sq) in gamv
   gamt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        tke_levels),                                                           &
                 ! counter gradient term for flux of TL
                 ! gamt = gamt_tsq * tsq + gamt_cov * cov + gamt_res
   gamt_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            tke_levels),                                                       &
                 ! a linear part to tsq in gamt
   gamt_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            tke_levels),                                                       &
                 ! a linear part to cov in gamt
   gamt_res(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            tke_levels),                                                       &
                 ! a residual part in gamt
   gamt_factor(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
               tke_levels),                                                    &
                 ! stability factor for gamt
   gamq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        tke_levels),                                                           &
                 ! counter gradient term for flux of QW
                 ! gamq = gamq_qsq * qsq + gamq_cov * cov + gamq_res
   gamq_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            tke_levels),                                                       &
                 ! a linear part to qsq in gamq
   gamq_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            tke_levels),                                                       &
                 ! a linear part to cov in gamq
   gamq_res(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            tke_levels),                                                       &
                 ! a residual part in gamq
   gamq_factor(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,            &
               tke_levels),                                                    &
                 ! stability factor for gamq
   pdc_factor(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,             &
              tke_levels),                                                     &
                 ! stability factor for pdc
   smd(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                 ! counter gradient correction for SM
   smd_coef(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            tke_levels),                                                       &
                 ! coefficient of (c3sq-c2sq) in smd
   e2(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                     &
      tke_levels),                                                             &
                 ! a variable denoted to E2 in the papers
   cu(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                 ! <u'**2> / qke in level 3
   cv(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                 ! <v'**2> / qke  in level 3
   cu25(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        tke_levels),                                                           &
                 ! <u'**2> / qke  in level 2.5
   cv25(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        tke_levels),                                                           &
                 ! <v'**2> / qke in level 2.5
   cw25(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        tke_levels),                                                           &
                 ! <w'**2> / qke = 1 - <u'**2> - <v'**2> in level 2.5
   pdk(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
        tke_levels),                                                           &
                 ! production term of qke
   pdt_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                 ! a linear part to tsq in a production term of tsq
   pdt_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                 ! a linear part to cov in a production term of tsq
   pdt(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                 ! production term of tsq excluding a linear part
                 ! of tsq and cov
                 ! (Note that it is the production term itself
                 ! in level 2.5)
   pdq_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                 ! a linear part to qsq in a production term of qsq
   pdq_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                 ! a linear part to cov in a production term of qsq
   pdq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                 ! production term of qsq excluding a linear part
                 ! of qsq and cov
                 ! (Note that it is the production term itself
                 ! in level 2.5)
   pdc_tsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                 ! a linear part to tsq in a production term of cov
   pdc_qsq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                 ! a linear part to qsq in a production term of cov
   pdc_cov(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                &
           tke_levels),                                                        &
                 ! a linear part to cov in a production term of cov
   pdc(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                    &
       tke_levels),                                                            &
                 ! production term of cov excluding a linear part
                 ! of cov, tsq and qsq
                 ! (Note that it is the production term itself
                 ! in level 2.5)
   bp(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                 ! coefficients of qke in the dissipation term
   rp(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                 ! production term of qke
   el(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,tke_levels),         &
                 ! mixing length
   q3sq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
        tke_levels),                                                           &
                 ! qke derived by level 2.5 or level 3
   c2sq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                   &
         tke_levels)
                 ! cov derived by level 2

real(kind=r_bl), allocatable ::                                         &
   ! These variables are required only when imp_mode /= FULL_IMPL
   ! So usually they are not used.
   ! (That is why they have an "allocatable" attribute.)
   bp_tsq(:, :, :),                                                            &
                 ! coefficients of tsq in the dissipation term
   rp_tsq(:, :, :),                                                            &
                 ! production term of tsq
   bp_qsq(:, :, :),                                                            &
                 ! coefficients of qsq in the dissipation term
   rp_qsq(:, :, :),                                                            &
                 ! production term of qsq
   bp_cov(:, :, :),                                                            &
                 ! coefficients of cov in the dissipation term
   rp_cov(:, :, :)
                 ! production term of cov

integer, parameter ::                                                          &
   ! Symbols for a switch
   full_impl = 0,                                                              &
   half_impl = 1,                                                              &
   expl      = 2

integer, parameter ::                                                          &
   imp_mode = full_impl
      ! mode to integrate covariances

character(len=*), parameter ::  RoutineName = 'MYM_TURBULENCE'

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
if (l_my_extra_level) then
  k_start = 1
else
  k_start = 2
end if

if (l_leonard_tke) then
  do k = 1, bl_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        leonard_kl_tke(i,j,k,1) = 1.0
        leonard_kl_tke(i,j,k,2) = 1.0
      end do
    end do
  end do
end if

call mym_level2(                                                               &
      bl_levels, dbdz, dvdzm, gm, gh, sm, sh)

call mym_length(                                                               &
      tdims%i_end,tdims%j_end,tdims_l%halo_i,tdims_l%halo_j,bl_levels,         &
      qke, z_uv, z_tq, dbdz, delta_smag, r_mosurf, fb_surf, qkw, el)


! This section calculates 'growing turbulence' correction
! qdiv = 1-alpha in UMDP 25,
! limited sm and sh
! and diagnostic Cu, Cv Cw, relative u,v,w contributions to TKE at level 2.5.
! Note, these appear not to be used apart from Cw which is used to bound C3sq.

do k = 2, tke_levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      elsq(i, j, k) = el(i, j, k) ** 2
      q2sq = b1 * elsq(i, j, k)                                                &
           * (sm(i, j, k) * gm(i, j, k) + sh(i, j, k) * gh(i, j, k))
      q3sq(i, j, k) = qkw(i, j, k) ** 2
      gmel(i, j, k) = gm(i, j, k) * elsq(i, j, k)
      ghel(i, j, k) = gh(i, j, k) * elsq(i, j, k)

      ! adjust SM and SH by SQRT(q3sq / q2sq)
      ! qdiv is 1-alpha in Sec 2.7.2
      ! q2sq is L2 (equilibrium) q2 (or u_t in PAC notation) so
      ! check is q2 is less than equilibrium, i.e. growing.
      if ( q3sq(i, j, k) < q2sq ) then
        ! Growing turbulence corrections as per Sec. 2.7.2
        qdiv(i, j, k) = sqrt(q3sq(i, j, k) / q2sq)
        sm(i, j, k) = sm(i, j, k) * qdiv(i, j, k)
        sh(i, j, k) = sh(i, j, k) * qdiv(i, j, k)

        ! Constants needed for final SM, SH
        ! Note gm is just shear and gh bouyancy gradient at this point, so
        ! GM = gm L^2 / u_t^2 =  gmel / u_t^2
        ! Hence these are L3 computation if qkw prognostic.
        ! e1 etc. have units of u_t^2 (i.e. q3sq)
        ! reden (the denominator) units of (u_t^2)^2
        ! sm, sh are properly dimensionless.
                ! e1= u_t^2[ 1 - (1-alpha^2) 9.0 a1 a2 (1-c2) G_H
                !     + (1-alpha^2) 6.0 a1^2 G_M ] = u_t^2 E1
        e1   = q3sq(i, j, k)                                                   &
             - e1c * ghel(i, j, k) * qdiv(i, j, k) ** 2
                ! e2 = u_t^2[ 1 - (1-alpha^2) 9.0 a1 a2 (1-c2) G_H]
        e2(i, j, k)   = q3sq(i, j, k)                                          &
             - e2c * ghel(i, j, k) * qdiv(i, j, k) ** 2
        e3 = e1 + e3c * ghel(i, j, k) * qdiv(i, j, k) ** 2
        e4 = e1 - e4c * ghel(i, j, k) * qdiv(i, j, k) ** 2
        eden = e2(i, j, k) * e4                                                &
             + e3 * e5c * gmel(i, j, k) * qdiv(i, j, k) ** 2
        eden = max(eden, 1.0e-20)
        reden = 1.0 / eden
      else
        ! Constants needed for final SM, SH
        e1 = q3sq(i, j, k) - e1c * ghel(i, j, k)
        e2(i, j, k) = q3sq(i, j, k) - e2c * ghel(i, j, k)
        e3 = e1 + e3c * ghel(i, j, k)
        e4 = e1 - e4c * ghel(i, j, k)
        eden = e2(i, j, k) * e4 + e3 * e5c * gmel(i, j, k)
        eden = max(eden, 1.0e-20)
        reden = 1.0 / eden
        qdiv(i, j, k) = 1.0
        ! L3 sm, sh if L3 qkw
        sm(i, j, k) = q3sq(i, j, k) * a1 * (e3 - 3.0 * c1 *e4)                 &
                           * reden
        sh(i, j, k) = q3sq(i, j, k)                                            &
             * a2 * (e2(i, j, k) + 3.0 * c1 * e5c * gmel(i, j, k))             &
             * reden
      end if ! test if q3sq < q2sq
      cu25(i, j, k) =(e2(i, j, k)                                              &
           + 3.0 * c1 * e5c * gmel(i, j, k)                                    &
           * qdiv(i, j, k) ** 2) * one_third * reden
      cv25(i, j, k) = cu25(i, j, k)                                            &
           * (e4 - 0.5 * e4c * ghel(i, j, k) * qdiv(i, j, k) ** 2)
      cw25(i, j, k) = cu25(i, j, k) * e1
      cu25(i, j, k) = 1.0 - cv25(i, j, k) - cw25(i, j, k)

    end do
  end do
end do
if (l_leonard_tke) then
  do k = 2, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        e1 = (el(i, j, k) / delta_smag(i,j))**2
        leonard_kl_tke(i,j,k,1) = 36.0 * a2 *  e1                              &
                                  * (sm(i, j, k) + (1.0 - c5) * sh(i, j, k) )
        leonard_kl_tke(i,j,k,2) = 72.0 * a1 * e1  * (1.0 - c5) * sm(i, j, k)
      end do
    end do
  end do
endif
! This section primarily computes L3 non-local terms, gamma_theta and gamma_q
! and corresponding perturbations to diffusion coefficients, dfm, dfh.

if ( levflag == mymodel3 ) then
  do k = 2, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end

        ! At L2, L2.5 C_theta  = b2 * sh

        t2sq = qdiv(i, j, k) * b2 * elsq(i, j, k)                              &
                     * sh(i, j, k) * dtldz(i, j, k) ** 2
        r2sq = qdiv(i, j, k) * b2 * elsq(i, j, k)                              &
                     * sh(i, j, k) * dqwdz(i, j, k) ** 2
        c2sq(i, j, k) = qdiv(i, j, k) * b2 * elsq(i, j, k)                     &
              * sh(i, j, k) * dtldz(i, j, k) * dqwdz(i, j, k)
        t3sq = max(tsq(i, j, k), 0.0)
        r3sq = max(qsq(i, j, k), 0.0)
        c3sq = cov(i, j, k)

        c3sq = sign( min( abs(c3sq), sqrt(t3sq*r3sq) ), c3sq )

        ! t2sq = vt <theta_L'^2> + vq <theta_L' q_t'> = <theta_L' theta_v'>_2.5
        t2sq = vt(i, j, k) * t2sq + vq(i, j, k) * c2sq(i, j, k)
        ! r2sq = vt <theta_L' q_t'>  + vq <q_t'^2> = <q_t' theta_v'>_2.5
        r2sq = vt(i, j, k) * c2sq(i, j, k) + vq(i, j, k) * r2sq
        ! c2sq = vt <theta_L' theta_v'>  + vq <q_t'theta_v'> = <theta_v'^2>_2.5
        c2sq(i, j, k) = max(vt(i, j, k) * t2sq + vq(i, j, k) * r2sq,           &
                            0.0)
        ! t3sq = vt <theta_L'^2> + vq <theta_L' q_t'> = <theta_L' theta_v'>
        t3sq = vt(i, j, k) * t3sq + vq(i, j, k) * c3sq
        ! r3sq = vt <theta_L' q_t'>  + vq <q_t'^2> = <q_t' theta_v'>
        r3sq = vt(i, j, k) * c3sq + vq(i, j, k) * r3sq
        ! c3sq = vt <theta_L' theta_v'>  + vq <q_t'theta_v'> = <theta_v'^2>
        c3sq = max(vt(i, j, k) * t3sq + vq(i, j, k) * r3sq, 0.0)

        !  Limitation on q, instead of L/q
        ! Note, gh = -db/dz so this applies when db/dz > 0.
        ! We could use -ghel here.
        if ( q3sq(i, j, k) < -gh(i, j, k) * elsq(i, j, k)) then
          q3sq(i, j, k) = -elsq(i, j, k) * gh(i, j, k)
        end if

        ! Some duplication from above.

        ! Limitation on c3sq (0.12 =< cw =< 0.76)
        ! e2 = q^2 * phi2'
        e2(i, j, k)   = q3sq(i, j, k)                                          &
                            - e2c*ghel(i, j, k) * qdiv(i, j, k)**2
        ! e3 = q^2 * phi3'
        e3   = q3sq(i, j, k) + e3c*ghel(i, j, k) * qdiv(i, j, k)**2
        ! e4 = q^2 * phi4'
        e4   = q3sq(i, j, k) - e4c*ghel(i, j, k) * qdiv(i, j, k)**2
        ! eden = q^4 D'
        eden = e2(i, j, k) * e4                                                &
                        + e3 *e5c*gmel(i, j, k) * qdiv(i, j, k)**2

        ! wden = numerator in the square bracket in (10a) in NN2006
        !        times (1-c3) * (g/thetav)**2 * GH
        wden = cc3*gtr(i, j, k) **2                                            &
                 * elsq(i, j, k)**2 / elsq(i, j, k)                            &
                 * qdiv(i, j, k)**2                                            &
                 *( e2(i, j, k)*e4c                                            &
                       - e3c*e5c*gmel(i, j, k) * qdiv(i, j, k)**2 )

        if ( wden /= 0.0 ) then
          clow = q3sq(i, j, k) * ( 0.12-cw25(i, j, k) )*eden/wden
          cupp = q3sq(i, j, k) *( 0.76-cw25(i, j, k) )*eden/wden

          if ( wden > 0.0 ) then
            c3sq  = min( max( c3sq, c2sq(i, j, k) + clow),                     &
                                     c2sq(i, j, k) + cupp)
          else
            c3sq  = max( min( c3sq, c2sq(i, j, k) + clow),                     &
                                     c2sq(i, j, k) + cupp)
          end if
        end if

        ! e1 has dimensions q^2
        e1   = e2(i, j, k) + e5c*gmel(i, j, k) * qdiv(i, j, k) ** 2
        ! eden has dimensions q^4
        eden = max( eden, 1.0e-20 )

        ! reden has dimensions q^-4
        reden = 1.0 / eden

        ! Surely elsq cancels?
!        e6c  = 3.0 * a2 *cc3 * gtr(i, j, k)                                    &
!                                   * elsq(i, j, k) / elsq(i, j, k)
        e6c  = 3.0 * a2 *cc3 * gtr(i, j, k)

        ! Calculate each term in  Gamma_theta
        ! coef is E_H g/theta_0 /ut^2 in UMDP 25
        coef = - e1 * qdiv(i, j, k) * e6c * reden
        ! Hence all Gamma coefficients are actually ~ut^(-2)
        gamt_tsq(i, j, k) = coef * vt(i, j, k)
        gamt_cov(i, j, k) = coef * vq(i, j, k)
        ! This appears to be missing a vt but its built into t2sq
        ! -Gamma_t_2.5  or T_theta
        gamt_res(i, j, k) = - coef * t2sq

        ! Calculate each term in  Gamma_q
        gamq_qsq(i, j, k) = coef * vq(i, j, k)
        gamq_cov(i, j, k) = coef * vt(i, j, k)
        ! -Gamma_q_2.5
        gamq_res(i, j, k) = - coef * r2sq

        ! for Sm' and Sh'd(Theta_V)/dz
        smd_coef(i, j, k)  = elsq(i, j, k) * qdiv(i, j, k) * e6c               &
             * gtr(i, j, k) * reden * qdiv(i, j, k) ** 2                       &
             * (e3c + e4c) * a1_2
        gamv_coef(i, j, k) = e1 * qdiv(i, j, k) * e6c * gtr(i, j, k)           &
                              * reden
        smd(i, j, k) = smd_coef(i, j, k) * (c3sq - c2sq(i, j, k))
        gamv(i, j, k) = gamv_coef(i, j, k) * (c3sq - c2sq(i, j,k))

        ! For elh (see below), qdiv in Level 3 is reset to 1.0.
        qdiv(i, j, k) = 1.0

        ! Calculate diffusion coefficients
        elq = el(i, j, k) * qkw(i, j, k)
        dfm(i, j, k) = elq * sm(i, j, k)
        dfh(i, j, k) = elq * sh(i, j, k)

      end do
    end do
  end do

  ! Adjustment for Gamma_theta and Gamma_q
  ! After the adjustment, Gamma_theta and Gamma_q are calculated
  if (l_my_prod_adj .and.                                                      &
        (imp_mode == half_impl .or. imp_mode == full_impl)) then
    do k = 2, tke_levels
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end

          elq = el(i, j, k) * qkw(i, j, k)
          elh = elq * qdiv(i, j, k)
          disp_coef = qkw(i, j, k) / (b2 * el(i, j, k))                        &
                   + 0.5 * coef_trbvar_diff * dfm(i, j, k)                     &
                      * (2.0 * pi * my_prod_adj_fact(k)                        &
                        / (z_uv(i, j, k) - z_uv(i, j, k - 1))) ** 2

          pdt_tsq(i, j, k) = elh * gamt_tsq(i, j, k) * dtldz(i, j, k)
          if (disp_coef < pdt_tsq(i, j, k)) then
            gamt_factor(i, j, k) = disp_coef / pdt_tsq(i, j, k)
          else
            gamt_factor(i, j, k) = 1.0
          end if

          pdq_qsq(i, j, k) = elh * gamq_qsq(i, j, k) * dqwdz(i, j, k)
          if (disp_coef < pdq_qsq(i, j, k)) then
            gamq_factor(i, j, k) = disp_coef / pdq_qsq(i, j, k)
          else
            gamq_factor(i, j, k) = 1.0
          end if

          gamt_tsq(i, j, k) = gamt_factor(i, j, k) * gamt_tsq(i, j, k)
          gamt_cov(i, j, k) = gamt_factor(i, j, k) * gamt_cov(i, j, k)
          gamt_res(i, j, k) = gamt_factor(i, j, k) * gamt_res(i, j, k)

          gamq_qsq(i, j, k) = gamq_factor(i, j, k) * gamq_qsq(i, j, k)
          gamq_cov(i, j, k) = gamq_factor(i, j, k) * gamq_cov(i, j, k)
          gamq_res(i, j, k) = gamq_factor(i, j, k) * gamq_res(i, j, k)

          pdc_cov(i, j, k) = elh                                               &
                        * (gamt_cov(i, j, k) * dqwdz(i, j, k)                  &
                         + gamq_cov(i, j, k) * dtldz(i, j, k)) * 0.5
          if (disp_coef < pdc_cov(i, j, k)) then
            pdc_factor(i, j, k) = disp_coef / pdc_cov(i, j, k)
          else
            pdc_factor(i, j, k) = 1.0
          end if
          gamt_tsq(i, j, k) = pdc_factor(i, j, k) * gamt_tsq(i, j, k)
          gamt_cov(i, j, k) = pdc_factor(i, j, k) * gamt_cov(i, j, k)
          gamt_res(i, j, k) = pdc_factor(i, j, k) * gamt_res(i, j, k)

          gamq_qsq(i, j, k) = pdc_factor(i, j, k) * gamq_qsq(i, j, k)
          gamq_cov(i, j, k) = pdc_factor(i, j, k) * gamq_cov(i, j, k)
          gamq_res(i, j, k) = pdc_factor(i, j, k) * gamq_res(i, j, k)

          gamt(i, j, k) = gamt_tsq(i, j, k) * tsq(i, j, k)                     &
                        + gamt_cov(i, j, k) * cov(i, j, k)                     &
                        + gamt_res(i, j, k)

          gamq(i, j, k) = gamq_qsq(i, j, k) * qsq(i, j, k)                     &
                        + gamq_cov(i, j, k) * cov(i, j, k)                     &
                        + gamq_res(i, j, k)

        end do
      end do
    end do
  else
    do k = 2, tke_levels
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          ! Non-local difference term:
          ! E_H g/theta_0 /ut^2 (<theta_L' theta_v'> - <theta_L' theta_v'>_2.5)
          gamt(i, j, k) = gamt_tsq(i, j, k) * tsq(i, j, k)                     &
                        + gamt_cov(i, j, k) * cov(i, j, k)                     &
                        + gamt_res(i, j, k)
          ! E_H g/theta_0 /ut^2 (<q_t' theta_v'> - <q_t' theta_v'>_2.5)
          gamq(i, j, k) = gamq_qsq(i, j, k) * qsq(i, j, k)                     &
                        + gamq_cov(i, j, k) * cov(i, j, k)                     &
                        + gamq_res(i, j, k)

          gamt_factor(i, j, k) = 1.0
          gamq_factor(i, j, k) = 1.0
          pdc_factor(i, j, k) = 1.0
        end do
      end do
    end do
  end if ! IF L_MY_PROD_ADJ

  ! gamt appears to be E_H/q^2 g/theta_v (beta_t (thetaL'_3^2 - thetaL'_2^2)
  !                                      +beta_q (theta_L'qt'_3 - thetaL'qt'_2))
  ! as per UMDP 25 eq.2.257
  ! Calculate production terms
  do k = 2, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end

        elq = el(i, j, k) * qkw(i, j, k)
        elh = elq * qdiv(i, j, k)

        ! TKE production
        ! Using Level 2.5 coefficients
        ! wb_ng is buoyancy flux related to the skewness (input)
        ! Level 2 production
        pdk(i, j, k) = elq * (sm(i, j, k) * gm(i, j, k)                        &
                                + sh(i, j, k) * gh(i, j, k))                   &
                                + wb_ng(i,j,k)

                ! pdt = elh * (sh * dtldz - Gamma_t_2.5) * dtl/dz
                !     = (P_theta^0 + R_theta)/2
        pdt(i, j, k) = elh                                                     &
                * (sh(i, j, k) * dtldz(i, j, k) + gamt_res(i, j, k))           &
                * dtldz(i, j, k)
                ! pdt_tsq = (1-alpha) L u_t E_H g/theta_0 beta_t / ut^2 dtl/dz
                !         = T_theta/2
        pdt_tsq(i, j, k) = elh * gamt_tsq(i, j, k) * dtldz(i, j, k)
                ! pdt_cov = (1-alpha) L u_t E_H g/theta_0 beta_q / ut^2 dtl/dz = p^tc/2
        pdt_cov(i, j, k) = elh * gamt_cov(i, j, k) * dtldz(i, j, k)


                ! pdq = elh * (sh * dqt/dz - Gamma_q_2.5) * dqt/dz
                !     = (P_q^0 + R_q)/2
        pdq(i, j, k) = elh                                                     &
             * (sh(i, j, k) * dqwdz(i, j, k) + gamq_res(i, j, k))              &
             * dqwdz(i, j, k)
                ! pdq_qsq = (1-alpha) L u_t E_H g/theta_0  beta_q /ut^2 dqt/dz = Q_q/2
        pdq_qsq(i, j, k) = elh * gamq_qsq(i, j, k) * dqwdz(i, j, k)
                ! pdq_cov = (1-alpha) L u_t E_H g/theta_0  beta_t /ut^2 dqt/dz = p^qc/2
        pdq_cov(i, j, k) = elh * gamq_cov(i, j, k) * dqwdz(i, j, k)

                ! pdc = 0.5 * elh * ((sh * dtldz - Gamma_t_2.5) * dqt/dz
                !                  + (sh * dqt/dz - Gamma_q_2.5) * dtl/dz)
                !     = (P_c^0 + R_c)/2
        pdc(i, j, k) = 0.5 * elh                                               &
             * ((sh(i, j, k) * dtldz(i, j, k)                                  &
                     + gamt_res(i, j, k)) * dqwdz(i, j, k)                     &
              + (sh(i, j, k) * dqwdz(i, j, k)                                  &
                     + gamq_res(i, j, k)) * dtldz(i, j, k))

                ! pdc_tsq = (1-alpha) L u_t E_H g/theta_0 beta_t / ut^2 * dqt/dz
                !         = p^ct/2
        pdc_tsq(i, j, k) = elh                                                 &
                         * gamt_tsq(i, j, k) * dqwdz(i, j, k) * 0.5
                ! pdc_qsq = (1-alpha) L u_t E_H g/theta_0  beta_q /ut^2 dtl/dz = p^cq/2
        pdc_qsq(i, j, k) = elh                                                 &
                         * gamq_qsq(i, j, k) * dtldz(i, j, k) * 0.5
                ! pdc_cov = 0.5 * (1-alpha) L u_t (
                !              E_H g/theta_0 beta_q / ut^2 dqt/dz
                !            + E_H g/theta_0  beta_t /ut^2 dtl/dz) = C_c/2
        pdc_cov(i, j, k) = 0.5 * elh                                           &
                         * (gamt_cov(i, j, k) * dqwdz(i, j, k)                 &
                          + gamq_cov(i, j, k) * dtldz(i, j, k))

        dfu_cg(i, j, k) = elq * smd(i, j, k) * dudz(i, j, k)
        dfv_cg(i, j, k) = elq * smd(i, j, k) * dvdz(i, j, k)
        dft_cg(i, j, k) = elq * gamt(i, j, k)
        dfq_cg(i, j, k) = elq * gamq(i, j, k)
      end do
    end do
  end do
else  ! level 2.5
  do k = 2, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        !     In Level 2.5, qdiv is not reset.
        gamt(i, j, k) = 0.0
        gamq(i, j, k) = 0.0
        gamv(i, j, k) = 0.0
        smd(i, j, k) = 0.0
        cu(i, j, k) = cu25(i, j, k)
        cv(i, j, k) = cv25(i, j, k)

        elq = el(i, j, k) * qkw(i, j, k)
        elh = elq * qdiv(i, j, k)

        pdk(i, j, k) = elq                                                     &
                    * (sm(i, j, k) * gm(i, j, k)                               &
                     + sh(i, j, k) * gh(i, j, k)) + wb_ng(i,j,k)
        pdt(i, j, k) = elh                                                     &
                    * (sh(i, j, k) * dtldz(i, j, k)) * dtldz(i, j, k)
        pdq(i, j, k) = elh                                                     &
                    * (sh(i, j, k) * dqwdz(i, j, k)) * dqwdz(i, j, k)
        pdc(i, j, k) = elh                                                     &
             * (sh(i, j, k) * dtldz(i, j, k)) * dqwdz(i, j, k) * 0.5           &
             + elh                                                             &
             * (sh(i, j, k) * dqwdz(i, j, k)) * dtldz(i, j, k) * 0.5

        dfm(i, j, k) = elq * sm(i, j, k)
        dfh(i, j, k) = elq * sh(i, j, k)
        dfu_cg(i, j, k) = 0.0
        dfv_cg(i, j, k) = 0.0
        dft_cg(i, j, k) = 0.0
        dfq_cg(i, j, k) = 0.0
      end do
    end do
  end do
end if  ! test if levflag == mymodel3

! Overwrite production terms by ones calculated with surface fluxes
if (my_lowest_pd_surf > 0) then
  if (l_my_extra_level) then
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        pdk(i, j, 1) = 1.0 * u_s(i, j) ** 3 * pmz(i, j)                        &
                    / (vkman * z_tq(i, j, 1) * my_z_extra_fact)
      end do
    end do
    if (l_my_lowest_pd_surf_tqc) then
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          phm  = 1.0 / u_s(i, j) * phh(i, j)                                   &
             / (vkman * z_tq(i, j, 1) * my_z_extra_fact)
          pdt(i, j, 1) = phm * ftl(i, j, 1) ** 2
          pdq(i, j, 1) = phm * fqw(i, j, 1) ** 2
          pdc(i, j, 1) = phm * ftl(i, j, 1) * fqw(i, j, 1)
        end do
      end do
    end if
  else    ! NOT L_MY_Extra_level
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        pdk(i, j, 2) = 1.0 * u_s(i, j) ** 3 * pmz(i, j)                        &
                             / (vkman * z_tq(i, j, 1))
        pdk(i, j, 1) = 0.0
        pdt(i, j, 1) = 0.0
        pdq(i, j, 1) = 0.0
        pdc(i, j, 1) = 0.0
      end do
    end do
    if (l_my_lowest_pd_surf_tqc) then
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          phm  = 1.0 / u_s(i, j)* phh(i, j)                                    &
                            / (vkman * z_tq(i, j, 1))
          pdt(i, j, 2) = phm * ftl(i, j, 1) ** 2
          pdq(i, j, 2) = phm * fqw(i, j, 1) ** 2
          pdc(i, j, 2) = phm * ftl(i, j, 1) * fqw(i, j, 1)
          pdt_tsq(i, j, 2) = 0.0
          pdt_cov(i, j, 2) = 0.0
          pdq_qsq(i, j, 2) = 0.0
          pdq_cov(i, j, 2) = 0.0
          pdc_tsq(i, j, 2) = 0.0
          pdc_qsq(i, j, 2) = 0.0
          pdc_cov(i, j, 2) = 0.0
        end do
      end do
    end if  ! IF L_MY_lowest_pd_surf_tqc
  end if ! IF L_MY_EXTRA_LEVEL
else  ! MY_lowest_pd_surf = 0
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      pdk(i, j, 1) = 0.0
      pdt(i, j, 1) = 0.0
      pdq(i, j, 1) = 0.0
      pdc(i, j, 1) = 0.0
    end do
  end do
end if  ! IF MY_lowest_pd_surf

! for diagnostics
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    gamt(i, j, 1) = 0.0
    gamq(i, j, 1) = 0.0
    gamv(i, j, 1) = 0.0
    smd(i, j, 1) = 0.0
    dfh(i, j, 1) = 0.0
    pdt_tsq(i, j, 1) = 0.0
    pdt_cov(i, j, 1) = 0.0
    pdq_qsq(i, j, 1) = 0.0
    pdq_cov(i, j, 1) = 0.0
    pdc_tsq(i, j, 1) = 0.0
    pdc_qsq(i, j, 1) = 0.0
    pdc_cov(i, j, 1) = 0.0
  end do
end do

if (BL_diag%l_tke_shr_prod) then
  do k = 2, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        elq = el(i, j, k) * qkw(i, j, k)
        BL_diag%tke_shr_prod(i, j, k) = el(i, j, k) * qkw(i, j, k)             &
                       * (sm(i, j, k) + smd(i, j, k)) * gm(i, j, k)
      end do
    end do
  end do
end if

if (BL_diag%l_tke_boy_prod) then
  do k = 2, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        BL_diag%tke_boy_prod(i, j, k) = el(i, j, k) * qkw(i, j, k)             &
                       * (sh(i, j, k) * gh(i, j, k)                            &
                                    + gamv(i, j, k)) + wb_ng(i,j,k)
      end do
    end do
  end do
end if

if (BL_diag%l_tke_dissp) then
  do k = k_start, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        BL_diag%tke_dissp(i, j, k) = qkw(i, j, k) ** 3                         &
                                              / (b1 * el(i, j, k))
      end do
    end do
  end do
end if

if (levflag == mymodel3) then
  ! Integrate the covariances

  if (imp_mode == full_impl) then
    call mym_update_covariance(                                                &
    ! IN levels
                bl_levels,                                                     &
    ! IN fields
                z_uv, z_tq,                                                    &
                qkw, el, dfm, pdt_tsq, pdt_cov, pdt,                           &
                pdq_qsq, pdq_cov, pdq, pdc_cov, pdc_tsq, pdc_qsq, pdc,         &
    ! INOUT fields
                tsq, qsq, cov)
  else   ! half implict or explicit
    allocate(bp_tsq(tdims%i_start:tdims%i_end,                                 &
                    tdims%j_start:tdims%j_end, tke_levels))
    allocate(rp_tsq(tdims%i_start:tdims%i_end,                                 &
                    tdims%j_start:tdims%j_end, tke_levels))
    allocate(bp_qsq(tdims%i_start:tdims%i_end,                                 &
                    tdims%j_start:tdims%j_end, tke_levels))
    allocate(rp_qsq(tdims%i_start:tdims%i_end,                                 &
                    tdims%j_start:tdims%j_end, tke_levels))
    allocate(bp_cov(tdims%i_start:tdims%i_end,                                 &
                    tdims%j_start:tdims%j_end, tke_levels))
    allocate(rp_cov(tdims%i_start:tdims%i_end,                                 &
                    tdims%j_start:tdims%j_end, tke_levels))

    if (imp_mode == half_impl) then
      do k = k_start, tke_levels
        do j = tdims%j_start, tdims%j_end
          do i = tdims%i_start, tdims%i_end
            pdt(i, j, k) = pdt(i, j, k)                                        &
                                + pdt_cov(i, j, k) * cov(i, j, k)

            pdq(i, j, k) = pdq(i, j, k)                                        &
                                + pdq_cov(i, j, k) * cov(i, j, k)

            pdc(i, j, k) = pdc(i, j, k)                                        &
                                + pdc_tsq(i, j, k) * tsq(i, j, k)              &
                                + pdc_qsq(i, j, k) * qsq(i, j, k)

            b2l = 2.0 * qkw(i, j, k) / (b2 * el(i, j, k))

            bp_tsq(i, j, k) = b2l - 2.0 * pdt_tsq(i, j, k)
            rp_tsq(i, j, k) = 2.0 * pdt(i, j, k)

            bp_qsq(i, j, k) = b2l - 2.0 * pdq_qsq(i, j, k)
            rp_qsq(i, j, k) = 2.0 * pdq(i, j, k)

            bp_cov(i, j, k) = b2l - 2.0 * pdc_cov(i, j, k)
            rp_cov(i, j, k) = 2.0 * pdc(i, j, k)
          end do
        end do
      end do
    else if (imp_mode == expl) then
      do k = k_start, tke_levels
        do j = tdims%j_start, tdims%j_end
          do i = tdims%i_start, tdims%i_end
            pdt(i, j, k) = pdt(i, j, k)                                        &
                                + pdt_tsq(i, j, k) * tsq(i, j, k)              &
                                + pdt_cov(i, j, k) * cov(i, j, k)
            pdq(i, j, k) = pdq(i, j, k)                                        &
                                + pdq_qsq(i, j, k) * qsq(i, j, k)              &
                                + pdq_cov(i, j, k) * cov(i, j, k)
            pdc(i, j, k) = pdc(i, j, k)                                        &
                                + pdc_cov(i, j, k) * cov(i, j, k)              &
                                + pdc_tsq(i, j, k) * tsq(i, j, k)              &
                                + pdc_qsq(i, j, k) * qsq(i, j, k)

            b2l = 2.0 * qkw(i, j, k) / (b2 * el(i, j, k))

            bp_tsq(i, j, k) = b2l
            rp_tsq(i, j, k) = 2.0 * pdt(i, j, k)

            bp_qsq(i, j, k) = b2l
            rp_qsq(i, j, k) = 2.0 * pdq(i, j, k)

            bp_cov(i, j, k) = b2l
            rp_cov(i, j, k) = 2.0 * pdc(i, j, k)
          end do
        end do
      end do
    end if
    call mym_update_fields(                                                    &
          bl_levels, coef_trbvar_diff, z_uv, z_tq, dfm, rp_tsq, bp_tsq,tsq)

    call mym_update_fields(                                                    &
          bl_levels, coef_trbvar_diff, z_uv, z_tq, dfm, rp_qsq, bp_qsq,qsq)

    call mym_update_fields(                                                    &
          bl_levels, coef_trbvar_diff, z_uv, z_tq, dfm, rp_cov, bp_cov,cov)

    deallocate(rp_cov)
    deallocate(bp_cov)
    deallocate(rp_qsq)
    deallocate(bp_qsq)
    deallocate(rp_tsq)
    deallocate(bp_tsq)

  end if  ! if imp_mode == FULL_IMPL
else  ! level 2.5
  ! In level 2.5, tsq, qsq, cov are diagnosed assuming balance between
  ! production and dissipation.
  do k = k_start, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        if (qkw(i, j, k) <= 1.0e-4) then
          b2l = 0.0
        else
          b2l = b2 * el(i, j, k) / qkw(i, j, k)
        end if
        tsq(i, j, k) = b2l * 2.0 * pdt(i, j, k)
        qsq(i, j, k) = b2l * 2.0 * pdq(i, j, k)
        cov(i, j, k) = b2l * 2.0 * pdc(i, j, k)
      end do
    end do
  end do
end if

if (levflag >= mymodel25) then
  ! predict qke
  if (my_lowest_pd_surf > 0) then
    k_start_cor = k_start + 1
  else
    k_start_cor = k_start
  end if

  if (levflag == mymodel3 .and.                                                &
         (imp_mode == half_impl .or. imp_mode == full_impl)) then
    ! add correction terms evaluated with integrated tsq, qsq and cov
    do k = k_start_cor, tke_levels
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end

          t3sq = max(tsq(i, j, k), 0.0)
          r3sq = max(qsq(i, j, k), 0.0)
          c3sq = cov(i, j, k)

          c3sq = sign( min( abs(c3sq), sqrt(t3sq*r3sq) ), c3sq )

          t3sq = vt(i, j, k) * t3sq + vq(i, j, k) * c3sq
          r3sq = vt(i, j, k) * c3sq + vq(i, j, k) * r3sq
          c3sq = max(vt(i, j, k) * t3sq + vq(i, j, k) * r3sq, 0.0)

          elq = el(i, j, k) * qkw(i, j, k)
          ! Non-local TKE production from shear.
          ! S_M' G_M = E_M (L/u_t^2 g/theta)^2
          smd(i, j, k) = smd_coef(i, j, k) * (c3sq - c2sq(i, j, k))

          pdk(i, j, k) = pdk(i, j, k) + elq                                    &
               * (smd(i, j, k) * gm(i, j, k)                                   &
                  + gamv_coef(i, j, k) * (c3sq- c2sq(i, j, k)))
        end do
      end do
    end do
  else
    do k = k_start_cor, tke_levels
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          pdk(i, j, k) = pdk(i, j, k)                                          &
                + el(i, j, k) * qkw(i, j, k)                                   &
                      * (smd(i, j, k) * gm(i, j, k) + gamv(i, j, k))
        end do
      end do
    end do
  end if ! if test levflag == mymodel3

  do k = k_start, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        b1l = b1 * el(i, j, k)
        bp(i, j, k) = 2.0 * qkw(i, j, k) / b1l
        rp(i, j, k) = 2.0 * pdk(i, j, k)
      end do
    end do
  end do

  call mym_update_fields(                                                      &
        bl_levels, coef_trbvar_diff_tke, z_uv, z_tq, dfm, rp, bp, qke)

else
   ! level 2
   ! diagnose qke
  do k = 2, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        b2l = b2 * el(i, j, k)
        qke(i, j, k) = (max(b2l * 2.0 * pdk(i, j, k), 0.0))                    &
                                                     ** two_thirds
      end do
    end do
  end do
end if  ! test if levflag >= mymodel25

do k = 1, tke_levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      qke(i, j, k) = min(max(qke(i, j, k), 1.0e-20), qke_max)
      tsq(i, j, k) = max(tsq(i, j, k), 0.0)
      qsq(i, j, k) = max(qsq(i, j, k), 0.0)
    end do
  end do
end do

do k = tke_levels + 1, bl_levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      qke(i, j, k) = 0.0
      tsq(i, j, k) = 0.0
      qsq(i, j, k) = 0.0
      cov(i, j, k) = 0.0
      dfm(i, j, k) = 0.0
      dfh(i, j, k) = 0.0
      dfu_cg(i, j, k) = 0.0
      dfv_cg(i, j, k) = 0.0
      dft_cg(i, j, k) = 0.0
      dfq_cg(i, j, k) = 0.0
    end do
  end do
end do

if (BL_diag%l_elm) then
  do k = 2, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        BL_diag%elm(i, j, k) = el(i, j, k)
      end do
    end do
  end do
end if

if (BL_diag%l_sm) then
  do k = 2, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        BL_diag%sm(i, j, k) = sm(i, j, k)
      end do
    end do
  end do
end if

if (BL_diag%l_sh) then
  do k = 2, tke_levels
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        BL_diag%sh(i, j, k) = sh(i, j, k)
      end do
    end do
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine mym_turbulence
end module mym_turbulence_mod
