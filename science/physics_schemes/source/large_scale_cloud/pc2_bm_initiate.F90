! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

!  PC2 Cloud Scheme: Initiation

module pc2_bm_initiate_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'PC2_BM_INITIATE_MOD'
contains

subroutine pc2_bm_initiate(                                                    &
!      Pressure related fields
 p_theta_levels, cumulus,                                                      &
!      Mulit-modal related fields
 tgrad_bm, bl_w_var,tau_dec_bm,tau_hom_bm,tau_mph_bm,                          &
 ri_bm, mix_len_bm, zh,zhsc,dzh,bl_type_7,                                     &
!      Array dimensions
 nlevels,zlcl_mixed,r_theta_levels,z_theta,                                    &
!      Prognostic Fields
   t, cf, cfl, cff, q, qcl, qcf,sskew, svar_turb, svar_bm, entzone,            &
   sl_modes, qw_modes, rh_modes, sd_modes,                                     &
   l_calc_diag, l_mixing_ratio)

use conversions_mod,       only: zerodegc
use water_constants_mod,   only: lc
use planet_constants_mod,  only: g,r,lcrcp,kappa,repsilon,cp
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim
use atm_fields_bounds_mod, only: pdims, tdims, pdims_l
use pc2_constants_mod,     only: bm_tiny,                                      &
                                 pc2init_logic_original,                       &
                                 pc2init_logic_simplified,                     &
                                 pc2init_logic_smooth
use cloud_inputs_mod,      only: i_pc2_init_logic, cloud_pc2_tol,              &
                                 l_bm_sigma_s_grad,                            &
                                 i_bm_ez_opt, i_bm_ez_orig, i_bm_ez_subcrit,   &
                                 i_bm_ez_entpar, turb_var_fac_bm, max_sigmas,  &
                                 min_sigx_ft, min_sigx_fac
use qsat_mod,              only: qsat_wat, qsat_wat_mix, qsat, qsat_mix

use pc2_total_cf_mod,    only: pc2_total_cf
use bm_cld_mod,          only: bm_cld
use bm_ql_mean_mod,      only: bm_ql_mean
use bm_ez_diagnosis_mod, only: bm_ez_diagnosis
use bm_entrain_parcel_mod, only: bm_entrain_parcel

implicit none

! Purpose:
!   Initiate liquid and total cloud fraction and liquid water content

! Method:
!   Use bimodal scheme with gaussian pdfs, linked to turbulence-based
!   variances. Define an entrainment zone based on the vertical profile
!   of liquid potential temperature and apply a mixture of two pdfs within
!   this entrainment zone: one mode from the free troposphere above the
!   inversion, and one mode from the bottom of the entrainment zone. Apply
!   weights to the modes so that the scheme is conservative for saturation
!   departure.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.
!
!   Documentation: PC2 Cloud Scheme Documentation

!  Subroutine Arguments:------------------------------------------------
integer ::                                                                     &
                      !, intent(in)
 nlevels
!       No. of levels being processed.

logical ::                                                                     &
                      !, intent(in)
 cumulus(       tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end),                                    &
!       Is this a boundary layer cumulus point
   l_mixing_ratio,                                                             &
! Use mixing ratio formulation
   l_calc_diag
! Flag to turn on calculation of extra diagnostics

! height of lcl in a well-mixed BL (types 3 or 4), 0 otherwise
real(kind=real_umphys) :: zlcl_mixed( pdims%i_start:pdims%i_end,               &
                                      pdims%j_start:pdims%j_end)


real(kind=real_umphys) ::                                                      &
                      !, intent(in)
 p_theta_levels(pdims%i_start:pdims%i_end,                                     &
                pdims%j_start:pdims%j_end,                                     &
                nlevels),                                                      &
!       Pressure at all points (Pa)
   r_theta_levels(pdims_l%i_start:pdims_l%i_end,                               &
                pdims_l%j_start:pdims_l%j_end,                                 &
                0:nlevels),                                                    &
!       Pressure at all points (Pa)
   bl_w_var(    tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,nlevels),                            &
!       Vertical velocity variance
   tgrad_bm(    tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                nlevels),                                                      &
!       Gradient of liquid potential temperature
   ri_bm(       tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                nlevels),                                                      &
!       Richardson Number for bimodal cloud scheme
   mix_len_bm(  tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                            1:tdims%k_end),                                    &
!       Turbulent mixing length for bimodal cloud scheme
   zh(          tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end),                                    &
!       Boundary-layer height for bimodal cloud scheme
   zhsc(        tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end),                                    &
!       Decoupled layer height for bimodal cloud scheme
   dzh(         tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end),                                    &
!       Inversion Thickness for bimodal cloud scheme
   bl_type_7(   tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end),                                    &
!       Shear-driven boundary layer indicator
   cff(         tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                nlevels),                                                      &
   z_theta(     tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                nlevels),                                                      &
!       Height of levels
   tau_dec_bm(  tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                nlevels),                                                      &
!       Decorrelation time scale
   tau_mph_bm(  tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                nlevels),                                                      &
!       Phase-relaxation time scale
   tau_hom_bm(  tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                nlevels),                                                      &
!       Turbulence homogenisation time scale
   qcf(         tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                nlevels)
!       Frozen content (kg water per kg air)

integer ::                                                                     &
 kez_inv(       tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                nlevels),                                                      &
!       if greater than 1, this is a flag indicating that this level belongs to
!       an entrainment zone (EZ) and its value is set to the k-level of the
!       inversion above this EZ.
 kez_top(       tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                nlevels),                                                      &
!       if greater than 1, this is the k-level identified as the mode
!       representative for the air above the inversion for the EZ that this
!       level belongs to
 kez_bottom(    tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                nlevels)
!       if greater than 1, this is the k-level identified as the mode from the
!       bottom of the EZ for the EZ that this level belongs to


real(kind=real_umphys) ::                                                      &
                      !, intent(INOUT)
 t(             tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                nlevels),                                                      &
!       Temperature (K)
   cf(          tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end  ,                                   &
                nlevels),                                                      &
!       Total cloud fraction (no units)
   cfl(         tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                nlevels),                                                      &
!       Liquid cloud fraction (no units)
   q(           tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                nlevels),                                                      &
!       Vapour content (kg water per kg air)
   qcl(         tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                nlevels)
!       Liquid content (kg water per kg air)


real(kind=real_umphys) ::                                                      &
                      !, intent(out)
   sskew(       tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,nlevels),                            &
!       Skewness of mixture of two s-distributions on levels
   svar_turb(   tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,nlevels),                            &
!       Variance of turbulence-based uni-modal s-distribution on levels
   svar_bm(     tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,nlevels),                            &
!       Variance of mixture of two s-distributions in the bi-modal scheme
   entzone(     tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,nlevels),                            &
!       Diagnostic indicating where entrainment zones are diagnosed
   sl_modes(    tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,nlevels,3),                          &
   qw_modes(    tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,nlevels,3),                          &
   rh_modes(    tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,nlevels,3),                          &
   sd_modes(    tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,nlevels,3)
!       Diagnostics of SL = T + g/cp z - Lc/cp qcl,  qw = q + qcl,
!       RHt = qw / qsat(Tl), and local turbulent standard-deviation of RHt
!       for the current level plus modes from above and below


!  Local scalars--------------------------------------------------------

integer :: k,i,j,                                                              &
!       Loop counters: K   - vertical level index
!                      I,J - horizontal position index
          qc_points
!                      Number of points to iterate over
!                      Indices for RHcrit array

!  Local scalars--------------------------------------------------------
real(kind=real_umphys) ::                                                      &
   alphl,                                                                      &
                      ! Local gradient of clausius-clapeyron
   alphal,                                                                     &
                      ! repsilon*lc/r
   mux,                                                                        &
                      ! Local first moment of the s-distribution
   sigx ( tdims%i_start:tdims%i_end, 3 ),                                      &
                      ! Local turbulence-based variance
   alx,                                                                        &
                      ! Local latent-heat correction term
   tlx,                                                                        &
                      ! Local liquid potential temperature
   qsl,                                                                        &
                      ! Local liquid saturation specific humidity
   qsi,                                                                        &
                      ! Local ice saturation specific humidity
   qnx_max ( tdims%i_start:tdims%i_end ),                                      &
   qnx_min ( tdims%i_start:tdims%i_end ),                                      &
                      ! Maximum and minimum value ofq qcl/bs
   qc,                                                                         &
                      ! Qc = al ( q + qcl - qsl(Tl) )
   frac_init,                                                                  &
                      ! Fraction of final liquid water initiated this timestep
   tmp
                      ! Work variable for calculating limit on turb variance

!  Local dynamic arrays-------------------------------------------------
real(kind=real_umphys) ::                                                      &
   cf_c(     tdims%i_len*tdims%j_len),                                         &
!       Total cloud fraction on compressed points (no units)
   cfl_c(    tdims%i_len*tdims%j_len),                                         &
!       Liquid cloud fraction on compressed points (no units)
   cff_c(    tdims%i_len*tdims%j_len),                                         &
!       Ice cloud fraction on compressed points (no units)
   deltacl_c(tdims%i_len*tdims%j_len),                                         &
!       Change in liquid cloud fraction (no units)
   deltaql_c(tdims%i_len*tdims%j_len),                                         &
!       Change in liquid cloud water (kg kg-1)
   deltacf_c(tdims%i_len*tdims%j_len)
!       Change in ice cloud fraction (no units)


real(kind=real_umphys) ::                                                      &
   qt_in(         tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,nlevels),                          &
!       Copy of initial q+qcl
   tl_in(         tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,nlevels),                          &
!       Copy of initial tliquid
   tl3d(          tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,nlevels),                          &
!       Copy of initial tliquid
   qt3d(          tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,nlevels),                          &
!       Copy of initial q+qcl
   cfl3d(         tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,nlevels),                          &
!       working copy for liquid cloud fraction
   qcl3d(         tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,nlevels),                          &
!       working copy for liquid cloud water
   cff3d(         tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,nlevels),                          &
!       working copy for ice cloud fraction
   qcf3d(         tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,nlevels),                          &
!       working copy for ice cloud water
   cf3d(         tdims%i_start:tdims%i_end,                                    &
                  tdims%j_start:tdims%j_end,nlevels),                          &
!       working copy for total cloud fraction
   cfl_max(        tdims%i_start:tdims%i_end,                                  &
                  tdims%j_start:tdims%j_end),                                  &
!       max cfl in column working downwards
   qsl_lay(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       Liquid saturation specific humidity for bimodal cloud scheme on layers
!       (level 1 is bottom-of-EZ mode, level 2 is k-level to calculate cloud
!       for and level 3 is mode from above inversion)
   qsi_lay(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       Ice saturation specific humidity for bimodal cloud scheme on layers
!       (level 1 is bottom-of-EZ mode, level 2 is k-level to calculate cloud
!       for and level 3 is mode from above inversion)
   ql_lay(        tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       Total water (vapour + liquid) for bimodal cloud scheme on layers
!       (level 1 is bottom-of-EZ mode, level 2 is k-level to calculate cloud
!       for and level 3 is mode from above inversion)
   tl_lay(        tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       Liquid water temperature for bimodal cloud scheme on layers
!       (level 1 is bottom-of-EZ mode, level 2 is k-level to calculate cloud
!       for and level 3 is mode from above inversion)
   wvar_lay(      tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       Local vertical velocity variance for bimodal cloud scheme on layers
!       (level 1 is bottom-of-EZ mode, level 2 is k-level to calculate cloud
!       for and level 3 is mode from above inversion)
   tdc_lay(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       Turbulence decorrelation time scale  for bimodal cloud scheme on layers
!       (level 1 is bottom-of-EZ mode, level 2 is k-level to calculate cloud
!       for and level 3 is mode from above inversion)
   inv_thm_lay(   tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       Inverse of turbulence homogenisation time scale for bimodal cloud
!       scheme on layers (level 1 is bottom-of-EZ mode, level 2 is k-level to
!       calculate cloud for and level 3 is mode from above inversion)
   inv_tmp_lay(   tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       Inverse of phase-relaxation time scale for bimodal cloud scheme
!       on layers (level 1 is bottom-of-EZ mode, level 2 is k-level to
!       calculate cloud for and level 3 is mode from above inversion)
   svar_ini(      tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
      ! Part of variance of turbulence-based uni-modal s-distribution on level
   dtldz_lay(     tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
   dqtdz_lay(     tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,3),                                &
!       Local vertical gradients dTl/dz and dqt/dz
   ql_mean(        tdims%i_start:tdims%i_end,                                  &
                   tdims%j_start:tdims%j_end, nlevels)
!       Mean ql in column below, used for tapering minimum-allowed variance
!       within the surface mixed-layer.

! Entraining parcel properties from above
real(kind=real_umphys) :: tl_above ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys) :: qt_above ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys) :: tau_dec_above ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys) :: tau_hom_above ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys) :: wvar_above ( tdims%i_start:tdims%i_end,              &
                                       tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys) :: dtldz_above ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys) :: dqtdz_above ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end, nlevels )
! Entraining parcel properties from below
real(kind=real_umphys) :: tl_below ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys) :: qt_below ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys) :: tau_dec_below ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys) :: tau_hom_below ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys) :: wvar_below ( tdims%i_start:tdims%i_end,              &
                                       tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys) :: dtldz_below ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end, nlevels )
real(kind=real_umphys) :: dqtdz_below ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end, nlevels )

logical ::                                                                     &
 lqc_max(       tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end)
!       True for points non-zero clouds for initiation

logical ::                                                                     &
 lqc_min(       tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end)
!       True for points with non-unity clouds for removal

! Flag for points where we calculate modes from above and below
logical :: l_set_modes ( tdims%i_start:tdims%i_end )

integer ::                                                                     &
 idx(tdims%j_len*tdims%i_len,2)

integer  :: kk, km1, kp1, kkm1, kkp1
!       indices for bimodal cloud scheme

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real   (kind=jprb)            :: zhook_handle

character(len=*), parameter :: RoutineName='PC2_BM_INITIATE'

!- End of Header

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set up a flag to state whether RHcrit is a single parameter or defined
! on all points.

! ==Main Block==--------------------------------------------------------

! ----------------------------------------------------------------------------
! --  Section 1 - initialisations                                           --
! --  Copy the initial q and t arrays to qt_in and tl_in arrays and         --
! --  initialise working arrays                                             --
! ----------------------------------------------------------------------------

alphl=repsilon*lc/r

!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED(nlevels,tdims,qt_in,tl_in,qt3d,tl3d,q,qcl,t,lcrcp,qcl3d,cfl3d,    &
!$OMP qcf,qcf3d,cff3d,cf3d,cff,cfl_max)                                        &
!$OMP private(k,j,i)
!$OMP do SCHEDULE(STATIC)
do k = 1, nlevels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      qt_in(i,j,k)     = q(i,j,k)+qcl(i,j,k)
      tl_in(i,j,k)     = t(i,j,k)-lcrcp*qcl(i,j,k)
      qt3d(i,j,k)      = q(i,j,k)+qcl(i,j,k)
      tl3d(i,j,k)      = t(i,j,k)-lcrcp*qcl(i,j,k)
      cff3d(i,j,k)     = cff(i,j,k)
      qcf3d(i,j,k)     = qcf(i,j,k)
      cfl3d(i,j,k)     = 0.0
      qcl3d(i,j,k)     = 0.0
      cf3d(i,j,k)      = 0.0
    end do
  end do
end do
!$OMP end do
!$OMP do SCHEDULE(STATIC)
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    cfl_max(i,j) = 0.0
  end do
end do
!$OMP end do
!$OMP end PARALLEL

! Find mean qt in each column, used for tapering the minimum-allowed
! variance in the surface mixed-layer.
call bm_ql_mean( nlevels, qt_in, p_theta_levels, ql_mean )

! ----------------------------------------------------------------------------
! --  Section 2 - Detection of entrainment zones                            --
! ----------------------------------------------------------------------------

select case ( i_bm_ez_opt )
case ( i_bm_ez_orig, i_bm_ez_subcrit )
  ! Options using entrainment zone diagnosis on model-levels

  call  bm_ez_diagnosis( p_theta_levels,tgrad_bm,z_theta,ri_bm,zh,zhsc,dzh,    &
                         bl_type_7,nlevels,tl_in,qt_in,                        &
                         l_mixing_ratio,kez_inv,kez_bottom,kez_top)

case ( i_bm_ez_entpar )
  ! Construct modes based on entraining parcels from above and below

  call bm_entrain_parcel(     nlevels,                                         &
                              zh, zhsc, dzh, bl_type_7,                        &
                              z_theta, bl_w_var, tau_dec_bm, tau_hom_bm,       &
                              mix_len_bm, tl_in, qt_in,                        &
                              tl_below, qt_below, wvar_below, tau_dec_below,   &
                              tau_hom_below, dtldz_below, dqtdz_below,         &
                              tl_above, qt_above, wvar_above, tau_dec_above,   &
                              tau_hom_above, dtldz_above, dqtdz_above )

end select  ! ( i_bm_ez_opt )

! ----------------------------------------------------------------------------
! -- Section 3 - Calculate cloud scheme input variables                      -
! --  Calculate the cloud scheme input values for the bottom of the EZ (1),  -
! --  the k-level of interest (2) and the mode from above the inversion (3). -
! --  Bring the modes from bottom of EZ and above the EZ dry adiabatically   -
! --  to the k-level to calculate cloud for.                                 -
! ----------------------------------------------------------------------------

    ! Loop over levels to prepare fields for the cloud scheme

!$OMP  PARALLEL                                                                &
!$OMP  DEFAULT(none)                                                           &
!$OMP  SHARED(nlevels,tdims,repsilon,r,g,cp,q,t,qt3d,tl3d,cfl3d,qcl3d,alphl,   &
!$OMP  lcrcp,kappa,sskew,svar_turb,svar_bm,tl_in,qt_in,p_theta_levels,         &
!$OMP  svar_ini,l_mixing_ratio,lqc_max,lqc_min,bl_w_var,                       &
!$OMP  kez_top,kez_bottom,kez_inv,tau_dec_bm,tau_hom_bm,tau_mph_bm,            &
!$OMP  cfl,cf,cff,qcl,qcf,i_pc2_init_logic,r_theta_levels,                     &
!$OMP  zlcl_mixed,cumulus,qcf3d,cff3d,cf3d,cfl_max,cloud_pc2_tol,              &
!$OMP  z_theta,l_bm_sigma_s_grad,                                              &
!$OMP  i_bm_ez_opt, turb_var_fac_bm, mix_len_bm, max_sigmas,                   &
!$OMP  min_sigx_ft, min_sigx_fac,                                              &
!$OMP  tl_above, qt_above, wvar_above, tau_dec_above, tau_hom_above,           &
!$OMP  dtldz_above, dqtdz_above,                                               &
!$OMP  tl_below, qt_below, wvar_below, tau_dec_below, tau_hom_below,           &
!$OMP  dtldz_below, dqtdz_below,                                               &
!$OMP  ql_mean, l_calc_diag,                                                   &
!$OMP  sl_modes, qw_modes, rh_modes, sd_modes, entzone)                        &
!$OMP  private(j,i,kk,qsl,qsi,alphal,alx,tlx,mux,sigx,deltacl_c,qc_points,idx, &
!$OMP  deltacf_c,cf_c,cfl_c,cff_c,deltaql_c,qnx_min,qnx_max,                   &
!$OMP  tl_lay,ql_lay,qsl_lay,qsi_lay,tdc_lay,inv_thm_lay,inv_tmp_lay,wvar_lay, &
!$OMP  qc, frac_init, km1,kp1,kkm1,kkp1,dtldz_lay,dqtdz_lay,tmp,l_set_modes )

do k = nlevels, 1, -1   ! need to work down for cfl_max
  ! Indices of levels above and below, but not allowed to go out-of-bounds!
  km1 = max( k - 1, 1 )
  kp1 = min( k + 1, nlevels )
!$OMP do SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end

      cfl3d(i,j,k)     = 0.0
      qcl3d(i,j,k)     = 0.0
      sskew(i,j,k)     = 0.0
      svar_turb(i,j,k) = 0.0
      svar_bm(i,j,k)   = 0.0

      tl_lay(i,j,1)      = 0.0
      ql_lay(i,j,1)      = 0.0
      qsl_lay(i,j,1)     = 0.0
      qsi_lay(i,j,1)     = 0.0
      tdc_lay(i,j,1)     = bm_tiny
      inv_thm_lay(i,j,1) = bm_tiny
      inv_tmp_lay(i,j,1) = bm_tiny
      wvar_lay(i,j,1)    = bm_tiny
      dtldz_lay(i,j,1)   = 0.0
      dqtdz_lay(i,j,1)   = 0.0

      tl_lay(i,j,3)      = 0.0
      ql_lay(i,j,3)      = 0.0
      qsl_lay(i,j,3)     = 0.0
      qsi_lay(i,j,3)     = 0.0
      tdc_lay(i,j,3)     = bm_tiny
      inv_thm_lay(i,j,3) = bm_tiny
      inv_tmp_lay(i,j,3) = bm_tiny
      wvar_lay(i,j,3)    = bm_tiny
      dtldz_lay(i,j,3)   = 0.0
      dqtdz_lay(i,j,3)   = 0.0

      !-----------------------------------------------------------------------
      ! Prepare layer fields to be passed on to the cloud scheme:
      ! - tl_lay      : liquid temperature, adiabatically brought to level k
      ! - q_lay       : total water
      ! - qsl_lay     : liquid saturation specific humidity
      ! - qsi_lay     : ice saturation specific humidity
      ! - tdc_lay     : turbulence decorrelation time scale
      ! - inv_thm_lay : inverse homogenisation time scale
      ! - inv_tmp_lay : inverse phase-relaxation time scale
      ! - wvar_lay    : vertical velocity variance
      !-----------------------------------------------------------------------

      tl_lay(i,j,2)      = tl_in(i,j,k)
      ql_lay(i,j,2)      = qt_in(i,j,k)
      tdc_lay(i,j,2)     = max(tau_dec_bm(i,j,k),bm_tiny)
      inv_thm_lay(i,j,2) = 1.0/max(tau_hom_bm(i,j,k),bm_tiny)
      inv_tmp_lay(i,j,2) = 1.0/max(tau_mph_bm(i,j,k),bm_tiny)
      wvar_lay(i,j,2)    = max(bl_w_var(i,j,k),bm_tiny)
      dtldz_lay(i,j,2)   = ( tl_in(i,j,kp1) - tl_in(i,j,km1) )                 &
                         / ( z_theta(i,j,kp1) - z_theta(i,j,km1) )
      dqtdz_lay(i,j,2)   = ( qt_in(i,j,kp1) - qt_in(i,j,km1) )                 &
                         / ( z_theta(i,j,kp1) - z_theta(i,j,km1) )
      svar_ini(i,j,2)    = max((0.5*bl_w_var(i,j,k)*tau_dec_bm(i,j,k))/        &
                               (1.0/tau_hom_bm(i,j,k)),bm_tiny)

      if ( l_mixing_ratio ) then
        call qsat_wat_mix(qsl,tl_in(i,j,k),p_theta_levels(i,j,k))
        call qsat_mix(qsi,tl_in(i,j,k),p_theta_levels(i,j,k))
      else
        call qsat_wat(qsl,tl_in(i,j,k),p_theta_levels(i,j,k))
        call qsat(qsi,tl_in(i,j,k),p_theta_levels(i,j,k))
      end if

      qsl_lay(i,j,2) = qsl
      qsi_lay(i,j,2) = qsi
      alphal = alphl * qsl / (tl_in(i,j,k) * tl_in(i,j,k))
      alx  = 1.0 / (1.0 + (lcrcp * alphal))

      if ( l_bm_sigma_s_grad ) then
        ! Account for local gradients in calculation of sigma_S
        sigx(i,2) = alx * sqrt(svar_ini(i,j,2))                                &
                        * abs( alphal*( g/cp + dtldz_lay(i,j,2) )              &
                             - dqtdz_lay(i,j,2) ) * turb_var_fac_bm
      else
        ! Don't account for local gradients
        sigx(i,2) = alx * sqrt(svar_ini(i,j,2))                                &
                        * alphal*g/cp * turb_var_fac_bm
      end if
      ! Impose minimum variance
      tmp = min_sigx_fac * max( ql_mean(i,j,k) - ql_lay(i,j,2), 0.0 )
      sigx(i,2) = max( sigx(i,2), max( 0.01*alx*qsl/max_sigmas,                &
                                       alx*qsl*min_sigx_ft*tmp/(tmp+qsl) ) )

      mux        = alx*(qt_in(i,j,k) - qsl)
      qnx_max(i) = mux/(max_sigmas*sigx(i,2))
      qnx_min(i) = mux/(max_sigmas*sigx(i,2))
      ! store first guess of the uni-modal turbulence-based variance. This
      ! will be updated in bm_cld in case adjustment for the latent heat
      ! relase term need to be made within the iteration loop.
      ! (note this actually stores standard deviation, not variance)
      svar_turb(i,j,k) = sigx(i,2)

    end do  ! i = tdims%i_start, tdims%i_end

      !-----------------------------------------------------------------------
      !-- If this level is identified as being part of an entrainment zone
      !-- with its inversion level stored in the kez_inv array, prepare
      !-- the variables for the bimodal scheme for the mode from above the
      !-- inversion and the mode from the bottom of the EZ.
      !-- Otherwise, these levels remain zero and the bimodal scheme will
      !-- revert to a uni-modal scheme with a gaussian pdf and
      !-- turbulence-based variances
      !-----------------------------------------------------------------------

    ! Set properties of modes from above and below, depending on i_bm_ez_opt...
    ! NOTE: there is no real reason why the qsat calls need to be duplicated
    ! identically in the 2 loops below, but moving them into another
    ! loop separate from the setting of the mode 1 and 3 Tl causes CCE
    ! to change answers, presumably due to optimisations behaving differently.
    ! NOTE: Cray fortran likes to split this loop up and vectorise it in a way
    ! that changes answers only if not using OpenMP.  NOFISSION directive is
    ! used to suppress this behaviour in order to preserve bit-comparison
    ! between OMP and NOOMP builds.
    if ( i_bm_ez_opt == i_bm_ez_entpar ) then
      ! Using entraining parcel option; set modes from above and below
      ! wherever the depth-scale is not negligible
#if defined (CRAYFTN_VERSION)
!DIR$ NOFISSION
#endif
      do i = tdims%i_start, tdims%i_end
        l_set_modes(i) = mix_len_bm(i,j,k) >                                   &
                         0.0001 * ( z_theta(i,j,kp1) - z_theta(i,j,km1) )
        if ( l_set_modes(i) ) then
          ! Copy properties of mode from above from subsided parcel arrays
          tlx = tl_above(i,j,k)
          if ( l_mixing_ratio ) then
            call qsat_wat_mix(qsl_lay(i,j,3),tlx,p_theta_levels(i,j,k))
            call qsat_mix(qsi_lay(i,j,3),tlx,p_theta_levels(i,j,k))
          else
            call qsat_wat(qsl_lay(i,j,3),tlx,p_theta_levels(i,j,k))
            call qsat(qsi_lay(i,j,3),tlx,p_theta_levels(i,j,k))
          end if
          tl_lay(i,j,3)      = tlx
          ql_lay(i,j,3)      = qt_above(i,j,k)
          wvar_lay(i,j,3)    = max(wvar_above(i,j,k),bm_tiny)
          tdc_lay(i,j,3)     = max(tau_dec_above(i,j,k),bm_tiny)
          inv_thm_lay(i,j,3) = 1.0/max(tau_hom_above(i,j,k),bm_tiny)
          dtldz_lay(i,j,3)   = dtldz_above(i,j,k)
          dqtdz_lay(i,j,3)   = dqtdz_above(i,j,k)
          svar_ini(i,j,3)    = max((0.5*wvar_above(i,j,k)*tau_dec_above(i,j,k))&
                                  /(1.0/tau_hom_above(i,j,k)),bm_tiny)
          ! Copy properties of mode from below from lifted parcel arrays
          tlx = tl_below(i,j,k)
          if ( l_mixing_ratio ) then
            call qsat_wat_mix(qsl_lay(i,j,1),tlx,p_theta_levels(i,j,k))
            call qsat_mix(qsi_lay(i,j,1),tlx,p_theta_levels(i,j,k))
          else
            call qsat_wat(qsl_lay(i,j,1),tlx,p_theta_levels(i,j,k))
            call qsat(qsi_lay(i,j,1),tlx,p_theta_levels(i,j,k))
          end if
          tl_lay(i,j,1)      = tlx
          ql_lay(i,j,1)      = qt_below(i,j,k)
          wvar_lay(i,j,1)    = max(wvar_below(i,j,k),bm_tiny)
          tdc_lay(i,j,1)     = max(tau_dec_below(i,j,k),bm_tiny)
          inv_thm_lay(i,j,1) = 1.0/max(tau_hom_below(i,j,k),bm_tiny)
          dtldz_lay(i,j,1)   = dtldz_below(i,j,k)
          dqtdz_lay(i,j,1)   = dqtdz_below(i,j,k)
          svar_ini(i,j,1)    = max((0.5*wvar_below(i,j,k)*tau_dec_below(i,j,k))&
                                  /(1.0/tau_hom_below(i,j,k)),bm_tiny)
        end if
      end do  ! i = tdims%i_start, tdims%i_end
    else  ! ( i_bm_ez_opt )
      ! Using discrete entrainment zone diagnosis; only set modes from
      ! above and below when inside an entrainment zone
#if defined (CRAYFTN_VERSION)
!DIR$ NOFISSION
#endif
      do i = tdims%i_start, tdims%i_end
        l_set_modes(i) = (k > 2) .and. (k < nlevels-2)                         &
                                 .and. (kez_inv(i,j,k) > 1)                    &
                                 .and. (kez_bottom(i,j,kez_inv(i,j,k)) < k)
        if ( l_set_modes(i) ) then
          ! Copy Properties of mode from above from model-level kez_top
          kk = kez_top(i,j,kez_inv(i,j,k))
          tlx = tl_in(i,j,kk)*(p_theta_levels(i,j,k)/                          &
                               p_theta_levels(i,j,kk))**kappa
          if ( l_mixing_ratio ) then
            call qsat_wat_mix(qsl_lay(i,j,3),tlx,p_theta_levels(i,j,k))
            call qsat_mix(qsi_lay(i,j,3),tlx,p_theta_levels(i,j,k))
          else
            call qsat_wat(qsl_lay(i,j,3),tlx,p_theta_levels(i,j,k))
            call qsat(qsi_lay(i,j,3),tlx,p_theta_levels(i,j,k))
          end if
          tl_lay(i,j,3)      = tlx
          ql_lay(i,j,3)      = qt_in(i,j,kk)
          wvar_lay(i,j,3)    = max(bl_w_var(i,j,kk),bm_tiny)
          tdc_lay(i,j,3)     = max(tau_dec_bm(i,j,kk),bm_tiny)
          inv_thm_lay(i,j,3) = 1.0/max(tau_hom_bm(i,j,kk),bm_tiny)
          kkm1 = max( kk - 1, 1 )
          kkp1 = min( kk + 1, nlevels )
          dtldz_lay(i,j,3)   = ( tl_in(i,j,kkp1) - tl_in(i,j,kkm1) )           &
                             / ( z_theta(i,j,kkp1) - z_theta(i,j,kkm1) )
          dqtdz_lay(i,j,3)   = ( qt_in(i,j,kkp1) - qt_in(i,j,kkm1) )           &
                             / ( z_theta(i,j,kkp1) - z_theta(i,j,kkm1) )
          svar_ini(i,j,3)    = max((0.5*bl_w_var(i,j,kk)*tau_dec_bm(i,j,kk))/  &
                                   (1.0/tau_hom_bm(i,j,kk)),bm_tiny)
          ! Copy Properties of mode from below from model-level kez_bottom
          kk = kez_bottom(i,j,kez_inv(i,j,k))
          tlx = tl_in(i,j,kk)*(p_theta_levels(i,j,k)/                          &
                               p_theta_levels(i,j,kk))**kappa
          if ( l_mixing_ratio ) then
            call qsat_wat_mix(qsl_lay(i,j,1),tlx,p_theta_levels(i,j,k))
            call qsat_mix(qsi_lay(i,j,1),tlx,p_theta_levels(i,j,k))
          else
            call qsat_wat(qsl_lay(i,j,1),tlx,p_theta_levels(i,j,k))
            call qsat(qsi_lay(i,j,1),tlx,p_theta_levels(i,j,k))
          end if
          tl_lay(i,j,1)      = tlx
          ql_lay(i,j,1)      = qt_in(i,j,kk)
          wvar_lay(i,j,1)    = max(bl_w_var(i,j,kk),bm_tiny)
          tdc_lay(i,j,1)     = max(tau_dec_bm(i,j,kk),bm_tiny)
          inv_thm_lay(i,j,1) = 1.0/max(tau_hom_bm(i,j,kk),bm_tiny)
          kkm1 = max( kk - 1, 1 )
          kkp1 = min( kk + 1, nlevels )
          dtldz_lay(i,j,1)   = ( tl_in(i,j,kkp1) - tl_in(i,j,kkm1) )           &
                             / ( z_theta(i,j,kkp1) - z_theta(i,j,kkm1) )
          dqtdz_lay(i,j,1)   = ( qt_in(i,j,kkp1) - qt_in(i,j,kkm1) )           &
                             / ( z_theta(i,j,kkp1) - z_theta(i,j,kkm1) )
          svar_ini(i,j,1)    = max((0.5*bl_w_var(i,j,kk)*tau_dec_bm(i,j,kk))/  &
                                   (1.0/tau_hom_bm(i,j,kk)),bm_tiny)
        end if
      end do  ! i = tdims%i_start, tdims%i_end
    end if  ! ( i_bm_ez_opt )

    do i = tdims%i_start, tdims%i_end

      if ( l_set_modes(i) ) then

        !---------------------------------------------------------------------
        ! Mode from above the inversion
        !---------------------------------------------------------------------

        ! Calculate dqs/dT
        qsl = qsl_lay(i,j,3)
        alphal = alphl*qsl / (tl_lay(i,j,3)*tl_lay(i,j,3))
        alx = 1.0 / (1.0 + (lcrcp * alphal))

        if ( l_bm_sigma_s_grad ) then
          ! Account for local gradients in calculation of sigma_S
          sigx(i,3) = alx * sqrt(svar_ini(i,j,3))                              &
                          * abs( alphal*( g/cp + dtldz_lay(i,j,3) )            &
                               - dqtdz_lay(i,j,3) ) * turb_var_fac_bm
        else
          ! Don't account for local gradients
          sigx(i,3) = alx * sqrt(svar_ini(i,j,3))                              &
                          * alphal*g/cp * turb_var_fac_bm
        end if
        ! Impose minimum variance
        tmp = min_sigx_fac * max( ql_mean(i,j,k) - ql_lay(i,j,3), 0.0 )
        sigx(i,3) = max( sigx(i,3), max( 0.01*alx*qsl/max_sigmas,              &
                                         alx*qsl*min_sigx_ft*tmp/(tmp+qsl) ) )

        ! Use local phase-relaxation time scale
        inv_tmp_lay(i,j,3) = 1.0/max(tau_mph_bm(i,j,k),bm_tiny)
        mux = alx*(ql_lay(i,j,3) - qsl)
        qnx_max(i) = max(qnx_max(i),mux/(max_sigmas*sigx(i,3)))
        qnx_min(i) = min(qnx_min(i),mux/(max_sigmas*sigx(i,3)))

        !---------------------------------------------------------------------
        ! Mode from entrainment zone bottom
        !---------------------------------------------------------------------

        ! Calculate dqs/dT
        qsl = qsl_lay(i,j,1)
        alphal = alphl*qsl / (tl_lay(i,j,1)*tl_lay(i,j,1))
        alx = 1.0 / (1.0 + (lcrcp * alphal))

        if ( l_bm_sigma_s_grad ) then
          ! Account for local gradients in calculation of sigma_S
          sigx(i,1) = alx * sqrt(svar_ini(i,j,1))                              &
                          * abs( alphal*( g/cp + dtldz_lay(i,j,1) )            &
                               - dqtdz_lay(i,j,1) ) * turb_var_fac_bm
        else
          ! Don't account for local gradients
          sigx(i,1) = alx * sqrt(svar_ini(i,j,1))                              &
                          * alphal*g/cp * turb_var_fac_bm
        end if
        ! Impose minimum variance
        tmp = min_sigx_fac * max( ql_mean(i,j,k) - ql_lay(i,j,1), 0.0 )
        sigx(i,1) = max( sigx(i,1), max( 0.01*alx*qsl/max_sigmas,              &
                                         alx*qsl*min_sigx_ft*tmp/(tmp+qsl) ) )

        ! Use local phase-relaxation time scale
        inv_tmp_lay(i,j,1) = 1.0/max(tau_mph_bm(i,j,k),bm_tiny)
        mux = alx*(ql_lay(i,j,1) - qsl)
        qnx_max(i) = max(qnx_max(i),mux/(max_sigmas*sigx(i,1)))
        qnx_min(i) = min(qnx_min(i),mux/(max_sigmas*sigx(i,1)))

      end if  ! ( l_set_modes(i) )

      lqc_max(i,j) = (qnx_max(i)  >   -1.0)
      lqc_min(i,j) = (qnx_min(i)  <    1.0)

    end do  ! i = tdims%i_start, tdims%i_end

    if ( l_calc_diag ) then
      ! Store properties of the modes for diagnostics...
      do i = tdims%i_start, tdims%i_end

        ! Copy properties of middle mode, from current level
        alphal = alphl * qsl_lay(i,j,2) / (tl_lay(i,j,2) * tl_lay(i,j,2))
        alx  = 1.0 / (1.0 + (lcrcp * alphal))
        sl_modes(i,j,k,2) = tl_lay(i,j,2) + (g/cp) * z_theta(i,j,k)
        qw_modes(i,j,k,2) = ql_lay(i,j,2)
        rh_modes(i,j,k,2) = ql_lay(i,j,2) / qsl_lay(i,j,2)
        sd_modes(i,j,k,2) = sigx(i,2) / ( alx * qsl_lay(i,j,2) )

        if ( l_set_modes(i) ) then
          ! Modes from above and below defined;
          ! Copy properties of mode from below
          alphal = alphl * qsl_lay(i,j,1) / (tl_lay(i,j,1) * tl_lay(i,j,1))
          alx  = 1.0 / (1.0 + (lcrcp * alphal))
          sl_modes(i,j,k,1) = tl_lay(i,j,1) + (g/cp) * z_theta(i,j,k)
          qw_modes(i,j,k,1) = ql_lay(i,j,1)
          rh_modes(i,j,k,1) = ql_lay(i,j,1) / qsl_lay(i,j,1)
          sd_modes(i,j,k,1) = sigx(i,1) / ( alx * qsl_lay(i,j,1) )
          ! Copy properties of mode from above
          alphal = alphl * qsl_lay(i,j,3) / (tl_lay(i,j,3) * tl_lay(i,j,3))
          alx  = 1.0 / (1.0 + (lcrcp * alphal))
          sl_modes(i,j,k,3) = tl_lay(i,j,3) + (g/cp) * z_theta(i,j,k)
          qw_modes(i,j,k,3) = ql_lay(i,j,3)
          rh_modes(i,j,k,3) = ql_lay(i,j,3) / qsl_lay(i,j,3)
          sd_modes(i,j,k,3) = sigx(i,3) / ( alx * qsl_lay(i,j,3) )
          ! Set entrainment zone indicator to 1
          entzone(i,j,k) = 1.0
        else
          ! Modes from above and below not defined; set diags to current level
          sl_modes(i,j,k,1) = sl_modes(i,j,k,2)
          qw_modes(i,j,k,1) = qw_modes(i,j,k,2)
          rh_modes(i,j,k,1) = rh_modes(i,j,k,2)
          sd_modes(i,j,k,1) = sd_modes(i,j,k,2)
          sl_modes(i,j,k,3) = sl_modes(i,j,k,2)
          qw_modes(i,j,k,3) = qw_modes(i,j,k,2)
          rh_modes(i,j,k,3) = rh_modes(i,j,k,2)
          sd_modes(i,j,k,3) = sd_modes(i,j,k,2)
          ! Set entrainment zone indicator to 0
          entzone(i,j,k) = 0.0
        end if

      end do  ! i = tdims%i_start, tdims%i_end
    end if  ! ( l_calc_diag )

  end do
!$OMP end do NOWAIT

  ! --------------------------------------------------------------------------
  ! -- Section 4 - Store the i,j indices of points that require cloud       --
  ! --             initiation or remgoval and fill compression arrays.      --
  ! --------------------------------------------------------------------------

  if ( i_pc2_init_logic == pc2init_logic_simplified ) then

    ! Copy points into compressed arrays
    qc_points = 0
!$OMP do SCHEDULE(STATIC)
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        if ( ( cfl(i,j,k) < cloud_pc2_tol                                      &
            .and. (r_theta_levels(i,j,k)-r_theta_levels(i,j,0))                &
             > zlcl_mixed(i,j)  .and.  lqc_max(i,j)  ) .or.                    &
            (cfl(i,j,k) > 1.0 - cloud_pc2_tol .and. lqc_min(i,j) )             &
             ) then

          qc_points = qc_points+1
          idx(qc_points,1) = i
          idx(qc_points,2) = j

        end if
      end do
    end do
!$OMP end do NOWAIT

  else if ( i_pc2_init_logic == pc2init_logic_original ) then

      ! Copy points into compressed arrays
    qc_points = 0
!$OMP do SCHEDULE(STATIC)
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        if ( ( ( ( .not. cfl(i,j,k) > 0.0_real_umphys )                        &
           .or. (cfl(i,j,k)  <   0.05_real_umphys .and. t(i,j,k) < zerodegc))  &
           .and. (.not. cumulus(i,j))                                          &
           .and. (r_theta_levels(i,j,k)-r_theta_levels(i,j,0))                 &
           > zlcl_mixed(i,j) .and.  lqc_max(i,j)) .or.                         &
           ( ( .not. cfl(i,j,k) < 1.0_real_umphys ) .and. lqc_min(i,j) )) then

          qc_points = qc_points+1
          idx(qc_points,1) = i
          idx(qc_points,2) = j

        end if
      end do
    end do
!$OMP end do NOWAIT

  else if ( i_pc2_init_logic == pc2init_logic_smooth ) then
    ! Smooth initiation logic; its possible we might want to initiate
    ! at least a little bit of extra liquid water at any point where
    ! we expect the bimodal scheme to diagnose non-zero qcl,
    ! i.e. wherever   Al (qt - qs(Tl)) > -3 sigx
    ! The result of this test has been stored in lqc_max:

    qc_points = 0
!$OMP do SCHEDULE(STATIC)
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        if ( lqc_max(i,j) ) then
          qc_points = qc_points+1
          idx(qc_points,1) = i
          idx(qc_points,2) = j
        end if
      end do
    end do
!$OMP end do NOWAIT

  end if  ! ( i_pc2_init_logic )

  ! --------------------------------------------------------------------------
  ! -- Section 5 - Call BM_CLD to calculate cloud water content,         --
  ! --             specific humidity, water cloud fraction and determine    --
  ! --             temperature.                                             --
  ! --------------------------------------------------------------------------
  ! Qc_points_if:


  if (qc_points  >   0) then

    call bm_cld(  p_theta_levels(1,1,k),                                       &
                  qsl_lay(:,:,:),                                              &
                  qsi_lay(:,:,:),                                              &
                  ql_lay(:,:,:),                                               &
                  tl_lay(:,:,:),                                               &
                  inv_thm_lay(:,:,:),                                          &
                  tdc_lay(:,:,:),                                              &
                  inv_tmp_lay(:,:,:),                                          &
                  wvar_lay(:,:,:),                                             &
                  dtldz_lay, dqtdz_lay, ql_mean(:,:,k),                        &
                  qcl3d(1,1,k),qcf3d(1,1,k),cfl3d(1,1,k),cff3d(1,1,k),         &
                  cf3d(1,1,k),cfl_max(1,1),qt3d(1,1,k),tl3d(1,1,k),            &
                  sskew(1,1,k),svar_turb(1,1,k),svar_bm(1,1,k),                &
                  idx,qc_points,l_mixing_ratio)


    ! ------------------------------------------------------------------------
    ! -- Section 6 - Calculate increments of cloud fraction and water       --
    ! --             content, calculate total cloud fraction and update     --
    ! --             fields                                                 --
    ! ------------------------------------------------------------------------

!DIR$ IVDEP
    do i = 1, qc_points
      ! cloud fraction and cloud liquid increments
      deltacl_c(i) = cfl3d(idx(i,1),idx(i,2),k)-cfl(idx(i,1),idx(i,2),k)
      deltacf_c(i) = 0.0
      deltaql_c(i) = qcl3d(idx(i,1),idx(i,2),k)-qcl(idx(i,1),idx(i,2),k)

      ! Gather variables
      cf_c (i) = cf(idx(i,1),idx(i,2),k)
      cfl_c(i) = cfl(idx(i,1),idx(i,2),k)
      cff_c(i) = cff(idx(i,1),idx(i,2),k)
    end do !i

    ! Calculate change in total cloud fraction.
    call pc2_total_cf(                                                         &
          qc_points,cfl_c,cff_c,deltacl_c,deltacf_c,cf_c)

    if ( i_pc2_init_logic == pc2init_logic_smooth ) then
      ! Smooth initiation logic; only use updated values
      ! where the diagnosed qcl exceeds the prognosed qcl
!DIR$ IVDEP
      do i = 1, qc_points
        if ( deltaql_c(i) > 0.0 ) then

          ! Calculate Qc (corresponds to the qcl we would have with
          ! no sub-grid moisture variability)

          ! qc = al ( q + qcl - qsl(Tl) )
          ! sd = al ( qsl(T) - q )
          ! qc + sd = al ( qcl + qsl(T) - qsl(Tl) )
          !         = al ( qcl + qsl(T) - qsl(T) + alpha lcrcp qcl )
          !         = al qcl ( 1 + alpha lcrcp )
          !         = qcl
          ! => sd = qcl - qc

          tlx = tl_in(idx(i,1),idx(i,2),k)
          if ( l_mixing_ratio ) then
            call qsat_wat_mix(qsl,tlx,p_theta_levels(idx(i,1),idx(i,2),k))
          else
            call qsat_wat(qsl,tlx,p_theta_levels(idx(i,1),idx(i,2),k))
          end if
          alphal = alphl * qsl / (tlx * tlx)
          alx  = 1.0 / (1.0 + (lcrcp * alphal))
          qc = alx * ( qt_in(idx(i,1),idx(i,2),k) - qsl )

          if ( qc < 0.0 ) then
            ! If qc<0 (total-water subsaturation)
            ! Find fraction of final qcl that was just created by initiation
            frac_init = deltaql_c(i) / qcl3d(idx(i,1),idx(i,2),k)
          else
            ! If qc>0 (total-water supersaturation)
            ! Find fraction of final sd = qcl-qc that was created by initiation
            frac_init = deltaql_c(i) / ( qcl3d(idx(i,1),idx(i,2),k)            &
                                       - min( qc, qcl(idx(i,1),idx(i,2),k) ) )
            ! (safety check avoids getting frac > 1 if qcl < qc, i.e. sd < 0
            !  which shouldn't really be possible!)
          end if

          ! Use the updated q, qcl and T at these points
          q  (idx(i,1),idx(i,2),k) = qt3d(idx(i,1),idx(i,2),k)
          qcl(idx(i,1),idx(i,2),k) = qcl3d(idx(i,1),idx(i,2),k)
          t  (idx(i,1),idx(i,2),k) = tl3d(idx(i,1),idx(i,2),k)

          ! New cloud fraction is weighted towards the value from the diagnostic
          ! scheme, in proportion to the fraction of qcl or sd created
          cf (idx(i,1),idx(i,2),k) =      frac_init  * cf_c(i)                 &
                                   + (1.0-frac_init) * cf (idx(i,1),idx(i,2),k)
          cfl(idx(i,1),idx(i,2),k) =   frac_init * cfl3d(idx(i,1),idx(i,2),k)  &
                                   + (1.0-frac_init) * cfl(idx(i,1),idx(i,2),k)

        end if
      end do
    else  ! ( i_pc2_init_logic )
      ! Other initiation logic options
!DIR$ IVDEP
      do i = 1, qc_points
        ! Update fields
        if ((cfl(idx(i,1),idx(i,2),k) <  0.5 .and. deltaql_c(i) > 0.0) .or.    &
            (cfl(idx(i,1),idx(i,2),k) >= 0.5 .and. deltaql_c(i) < 0.0)) then

          q  (idx(i,1),idx(i,2),k) = qt3d(idx(i,1),idx(i,2),k)
          qcl(idx(i,1),idx(i,2),k) = qcl(idx(i,1),idx(i,2),k)+deltaql_c(i)
          t  (idx(i,1),idx(i,2),k) = tl3d(idx(i,1),idx(i,2),k)
          cf (idx(i,1),idx(i,2),k) = cf_c(i)
          cfl(idx(i,1),idx(i,2),k) = cfl_c(i) + deltacl_c(i)
        end if
      end do !i
    end if  ! ( i_pc2_init_logic )

  end if ! Qc_points_if

end do !k loop
!$OMP end PARALLEL

! End of the subroutine

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine pc2_bm_initiate
end module pc2_bm_initiate_mod
