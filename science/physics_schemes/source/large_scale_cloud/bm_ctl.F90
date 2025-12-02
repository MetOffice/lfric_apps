! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module bm_ctl_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='BM_CTL_MOD'

contains
! Bimodal cloud Scheme.
! Subroutine Interface:
subroutine bm_ctl(                                                             &
!      Pressure related fields
 p_theta_levels,                                                               &
!      Time scales for variances
 tgrad_bm, bl_w_var, tau_dec_bm, tau_hom_bm, tau_mph_bm, z_theta,              &
 ri_bm, mix_len_bm, zh, zhsc, dzh, bl_type_7,                                  &
!      Array dimensions
 levels,                                                                       &
!      From convection diagnosis (only used if A05_4A)
 l_mixing_ratio,                                                               &
!      Prognostic Fields
 t, cf, q, qcf, qcl,                                                           &
!      Liquid and frozen ice cloud fractions
 cfl, cff,                                                                     &
 sskew, svar_turb, svar_bm, entzone,                                           &
 sl_modes, qw_modes, rh_modes, sd_modes,                                       &
 l_calc_diag, errorstatus)

use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim
use atm_fields_bounds_mod, only: pdims, tdims
use cloud_inputs_mod,      only: ice_fraction_method, cloud_top_temp,          &
                                 l_bm_sigma_s_grad,                            &
                                 i_bm_ez_opt, i_bm_ez_orig, i_bm_ez_subcrit,   &
                                 i_bm_ez_entpar, turb_var_fac_bm, max_sigmas,  &
                                 min_sigx_ft, min_sigx_fac
use planet_constants_mod,  only: g,r,lcrcp, kappa, repsilon,cp
use water_constants_mod,   only: lc
use pc2_constants_mod,     only: bm_tiny
use ereport_mod,           only: ereport
use errormessagelength_mod,only: errormessagelength
use umPrintMgr,            only: newline
use qsat_mod,              only: qsat_wat, qsat_wat_mix, qsat_mix, qsat
use bm_cld_mod,            only: bm_cld
use bm_ql_mean_mod,        only: bm_ql_mean
use bm_ez_diagnosis_mod,   only: bm_ez_diagnosis
use bm_entrain_parcel_mod, only: bm_entrain_parcel


implicit none

! Purpose:
!   This subroutine calculates liquid and ice cloud fractional cover
!   for use with the enhanced precipitation microphysics scheme, using the
!   diagnostic bimodal cloud scheme

! Method:
!   Diagnostic scheme using gaussian pdfs, linked to turbulence-based
!   variances. Define an entrainment zone based on the vertical profile
!   of liquid potential temperature and apply a mixture of two pdfs
!   within this entrainment zone: one mode from the free troposphere above
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

!   Documentation: UMDP No. 39


!  Global Variables:----------------------------------------------------

!  Subroutine Arguments:------------------------------------------------
integer ::                                                                     &
                      !, intent(in)
 levels
!       No. of levels being processed.

real(kind=real_umphys) ::                                                      &
                      !, intent(in)
 qcf(           tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       Cloud ice content at processed levels (kg water per kg air).
   p_theta_levels(pdims%i_start:pdims%i_end,                                   &
                  pdims%j_start:pdims%j_end,levels),                           &
!       pressure at all points (Pa).
   bl_w_var(    tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       Vertical velocity variance.
   z_theta(     tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       Height of levels
   tgrad_bm(    tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       Gradient of liquid potential temperature
   ri_bm(       tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
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
   tau_dec_bm(  tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       Decorrelation time scale
   tau_mph_bm(  tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       Phase-relaxation time scale
   tau_hom_bm(  tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels)
!       Turbulence homogenisation time scale

logical ::                                                                     &
 l_mixing_ratio,                                                               &
! true if using mixing ratios
   l_calc_diag
! Flag to turn on calculation of extra diagnostics

real(kind=real_umphys) ::                                                      &
                      !, intent(INOUT)
 q(             tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       On input : Total water content (QW) (kg per kg air).
!       On output: Specific humidity at processed levels
!                  (kg water per kg air).
   t(           tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels)
!       On input : Liquid/frozen water temperature (TL) (K).
!       On output: Temperature at processed levels (K).

real(kind=real_umphys) ::                                                      &
                      !, intent(out)
 cf(            tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       Cloud fraction at processed levels (decimal fraction).
   qcl(         tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       Cloud liquid water content at processed levels (kg per kg air).
   cfl(         tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       Liquid cloud fraction at processed levels (decimal fraction).
   cff(         tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       Frozen cloud fraction at processed levels (decimal fraction).
   sskew(       tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       Skewness of mixture of two s-distributions on levels
   svar_turb(   tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       Variance of turbulence-based uni-modal s-distribution on levels
   svar_bm(     tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       Variance of mixture of two s-distributions in the bi-modal scheme
   entzone(     tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       Diagnostic indicating where entrainment zones are diagnosed
   sl_modes(    tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels,3),                           &
   qw_modes(    tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels,3),                           &
   rh_modes(    tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels,3),                           &
   sd_modes(    tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels,3)
!       Diagnostics of SL = T + g/cp z - Lc/cp qcl,  qw = q + qcl,
!       RHt = qw / qsat(Tl), and local turbulent standard-deviation of RHt
!       for the current level plus modes from above and below

!     Error Status:
integer :: errorstatus     !, intent(out)  0 if OK; 1 if bad arguments.

!  Local scalars--------------------------------------------------------
real(kind=real_umphys) ::                                                      &
 alphal,                                                                       &
                      ! Local gradient of clausius-clapeyron
 alphl,                                                                        &
                      ! repsilon*lc/r
 mux,                                                                          &
                      ! Local first moment of the s-distribution
 sigx ( tdims%i_start:tdims%i_end, 3 ),                                        &
                      ! Local unimodal turbulence-based variance
 alx,                                                                          &
                      ! Local latent-heat correction term
 tlx,                                                                          &
                      ! Local liquid potential temperature
 qsl,                                                                          &
                      ! Local liquid saturation specific humidity
 tmp
                      ! Work variable for calculating limit on turb variance

!  (b) Others.
integer :: k,i,j      ! Loop counters: K - vertical level index.
!                                      I,J - horizontal field indices.

integer :: qc_points
                      ! No. points with non-zero cloud

!  Local dynamic arrays-------------------------------------------------
integer ::                                                                     &
 kez_inv(       tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                levels),                                                       &
!       If greater than 1, this is a flag indicating that this level belongs to
!       an entrainment zone (EZ) and its value is set to the k-level of the
!       inversion above this EZ.
 kez_top(       tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                levels),                                                       &
!       If greater than 1, this is the k-level identified as the mode
!       representative for the air above the inversion for the EZ that this
!       level belongs to.
 kez_bottom(    tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,                                     &
                levels)
!       If greater than 1, this is the k-level identified as the mode from the
!       bottom of the EZ for the EZ that this level belongs to.

real ::                                                                        &
   q_in(        tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels),                             &
!       Copy of initial q
   t_in(        tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,levels)
!       Copy of initial t

real(kind=real_umphys) ::                                                      &
   qsl_lay(     tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,3),                                  &
!       Liquid saturation specific humidity for bimodal cloud scheme on layers
!       (level 1 is bottom-of-EZ mode, level 2 is k-level to calculate cloud
!       for and level 3 is mode from above inversion)
   qsi_lay(     tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,3),                                  &
!       Ice saturation specific humidity for bimodal cloud scheme on layers
!       (level 1 is bottom-of-EZ mode, level 2 is k-level to calculate cloud
!       for and level 3 is mode from above inversion)
   ql_lay(      tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,3),                                  &
!       Total water (vapour + liquid) for bimodal cloud scheme on layers
!       (level 1 is bottom-of-EZ mode, level 2 is k-level to calculate cloud
!       for and level 3 is mode from above inversion)
   tl_lay(      tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,3),                                  &
!       Liquid water temperature for bimodal cloud scheme on layers
!       (level 1 is bottom-of-EZ mode, level 2 is k-level to calculate cloud
!       for and level 3 is mode from above inversion)
   wvar_lay(    tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,3),                                  &
!       Local vertical velocity variance for bimodal cloud scheme on layers
!       (level 1 is bottom-of-EZ mode, level 2 is k-level to calculate cloud
!       for and level 3 is mode from above inversion)
   tdc_lay(     tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,3),                                  &
!       Turbulence decorrelation time scale  for bimodal cloud scheme on layers
!       (level 1 is bottom-of-EZ mode, level 2 is k-level to calculate cloud
!       for and level 3 is mode from above inversion)
   inv_thm_lay( tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,3),                                  &
!       Inverse of turbulence homogenisation time scale for bimodal cloud
!       scheme on layers (level 1 is bottom-of-EZ mode, level 2 is k-level to
!       calculate cloud for and level 3 is mode from above inversion)
   inv_tmp_lay( tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end,3),                                  &
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
                   tdims%j_start:tdims%j_end, levels),                         &
!       Mean ql in column below, used for tapering minimum-allowed variance
!       within the surface mixed-layer.
   qnx(tdims%i_start:tdims%i_end),                                             &
!       Local cloud water normalised with BS.
   qsi_v(tdims%i_start:tdims%i_end),                                           &
!       Vector version for qsi along index i
   qsl_v(tdims%i_start:tdims%i_end)
!       Vector version for qsl along index i

real(kind=real_umphys) ::                                                      &
   cfl_max( tdims%i_start:tdims%i_end,                                         &
            tdims%j_start:tdims%j_end)
!       Max liquid cloud fraction in column, working down (decimal fraction).

! Entraining parcel properties from above
real(kind=real_umphys) :: tl_above ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end, levels )
real(kind=real_umphys) :: qt_above ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end, levels )
real(kind=real_umphys) :: tau_dec_above ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, levels )
real(kind=real_umphys) :: tau_hom_above ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, levels )
real(kind=real_umphys) :: wvar_above ( tdims%i_start:tdims%i_end,              &
                                       tdims%j_start:tdims%j_end, levels )
real(kind=real_umphys) :: dtldz_above ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end, levels )
real(kind=real_umphys) :: dqtdz_above ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end, levels )
! Entraining parcel properties from below
real(kind=real_umphys) :: tl_below ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end, levels )
real(kind=real_umphys) :: qt_below ( tdims%i_start:tdims%i_end,                &
                                     tdims%j_start:tdims%j_end, levels )
real(kind=real_umphys) :: tau_dec_below ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, levels )
real(kind=real_umphys) :: tau_hom_below ( tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end, levels )
real(kind=real_umphys) :: wvar_below ( tdims%i_start:tdims%i_end,              &
                                       tdims%j_start:tdims%j_end, levels )
real(kind=real_umphys) :: dtldz_below ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end, levels )
real(kind=real_umphys) :: dqtdz_below ( tdims%i_start:tdims%i_end,             &
                                        tdims%j_start:tdims%j_end, levels )

logical ::                                                                     &
 lcld(          tdims%i_start:tdims%i_end,                                     &
                tdims%j_start:tdims%j_end)
!       True for points with non-zero cloud

! Flag for points where we calculate modes from above and below
logical :: l_set_modes ( tdims%i_start:tdims%i_end )

integer ::                                                                     &
 idx(tdims%j_len*tdims%i_len,2)
!       Index for points with non-zero cloud


integer  :: kk, km1, kp1, kkm1, kkp1
!       indices for bimodal cloud scheme


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='BM_CTL'

character(len=errormessagelength) :: comments

!- End of Header

! ----------------------------------------------------------------------
!  Check input arguments for potential over-writing problems.
! ----------------------------------------------------------------------
errorstatus=0

alphl=repsilon*lc/r

! ==Main Block==--------------------------------------------------------
! Subroutine structure :
! Loop round levels to be processed.
! ----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------------
! --  Section 1 - initialisations                                           --
! --  Copy the initial q and t arrays to q_in and t_in arrays               --
! ----------------------------------------------------------------------------

!$OMP PARALLEL do                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(none)                                                            &
!$OMP SHARED(levels,tdims,q_in,t_in,q,t)                                       &
!$OMP private(k,j,i)
do k = 1, levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      q_in(i,j,k)       = q(i,j,k)
      t_in(i,j,k)       = t(i,j,k)
    end do
  end do
end do
!$OMP end PARALLEL do

! Find mean qt in each column, used for tapering the minimum-allowed
! variance in the surface mixed-layer.
call bm_ql_mean( levels, q_in, p_theta_levels, ql_mean )

! ----------------------------------------------------------------------------
! --  Section 2 - Detection of entrainment zones                            --
! ----------------------------------------------------------------------------

select case ( i_bm_ez_opt )
case ( i_bm_ez_orig, i_bm_ez_subcrit )
  ! Options using entrainment zone diagnosis on model-levels

  call  bm_ez_diagnosis( p_theta_levels,tgrad_bm,z_theta,ri_bm,zh,zhsc,dzh,    &
                         bl_type_7,levels,t_in,q_in,                           &
                         l_mixing_ratio,kez_inv,kez_bottom,kez_top)

case ( i_bm_ez_entpar )
  ! Construct modes based on entraining parcels from above and below

  call bm_entrain_parcel(     levels,                                          &
                              zh, zhsc, dzh, bl_type_7,                        &
                              z_theta, bl_w_var, tau_dec_bm, tau_hom_bm,       &
                              mix_len_bm, t_in, q_in,                          &
                              tl_below, qt_below, wvar_below, tau_dec_below,   &
                              tau_hom_below, dtldz_below, dqtdz_below,         &
                              tl_above, qt_above, wvar_above, tau_dec_above,   &
                              tau_hom_above, dtldz_above, dqtdz_above )

end select  ! ( i_bm_ez_opt )

! -----------------------------------------------------------------------------
! --  Section 3 - Check that we are not using ice_fraction_method = 2        --
! -----------------------------------------------------------------------------

if (ice_fraction_method  ==  cloud_top_temp) then
  errorstatus = 1

  comments = ' When using the bimodal cloud scheme,'//newline//                &
             ' ice_fraction_method = 2 (cloud_top_temp) is not '//newline//    &
             ' implemented. Please select option 1 (smith_orig), '//newline//  &
             ' 3 (min_liq_overlap) or 4 (min_liq_overlap_col) or '//newline//  &
             ' 5 (all-or-nothing) instead.'

  call ereport(RoutineName,errorstatus,comments)
end if

! ----------------------------------------------------------------------------
! -- Section 4 - Calculate cloud scheme input variables                     --
! --  Calculate the cloud scheme input values for the bottom of the EZ (1), --
! --  the k-level of interest (2) and the mode from above the inversion (3).--
! --  Bring the modes from bottom of EZ and above the EZ dry adiabatically  --
! --  to the k-level to calculate cloud for.                                --
! ----------------------------------------------------------------------------

!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED(levels,tdims,qcl,cfl,sskew,svar_turb,svar_bm,t_in,q_in,z_theta,   &
!$OMP        p_theta_levels,svar_ini,l_mixing_ratio,g,cp,r,qcf,cff,cf,cfl_max, &
!$OMP        lcrcp,kappa,repsilon,lcld,alphl,kez_top,kez_bottom,kez_inv,t,q,   &
!$OMP        tau_dec_bm,tau_hom_bm,tau_mph_bm,bl_w_var,l_bm_sigma_s_grad,      &
!$OMP        i_bm_ez_opt, turb_var_fac_bm, mix_len_bm, max_sigmas,             &
!$OMP        min_sigx_ft, min_sigx_fac,                                        &
!$OMP        tl_above, qt_above, wvar_above, tau_dec_above, tau_hom_above,     &
!$OMP        dtldz_above, dqtdz_above,                                         &
!$OMP        tl_below, qt_below, wvar_below, tau_dec_below, tau_hom_below,     &
!$OMP        dtldz_below, dqtdz_below,                                         &
!$OMP        ql_mean, l_calc_diag,                                             &
!$OMP        sl_modes, qw_modes, rh_modes, sd_modes, entzone)                  &
!$OMP private(j,i,k,kk,qsl,alphal,alx,qnx,tlx,mux,sigx,qc_points,idx,          &
!$OMP         tl_lay,ql_lay,qsl_lay,qsi_lay,tdc_lay,inv_thm_lay,inv_tmp_lay,   &
!$OMP         wvar_lay, qsl_v, qsi_v, km1, kp1, kkm1, kkp1,                    &
!$OMP         dtldz_lay, dqtdz_lay, tmp, l_set_modes)

!$OMP do SCHEDULE(STATIC)
do k = 1, levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end

      ! Initialise values

      qcl(i,j,k)         = 0.0
      cfl(i,j,k)         = 0.0
      cff(i,j,k)         = 0.0
      cf(i,j,k)          = 0.0
      sskew(i,j,k)       = 0.0
      svar_bm(i,j,k)     = 0.0
    end do
  end do
end do
!$OMP end do

!$OMP do SCHEDULE(STATIC)
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    ! Initialise values
    cfl_max(i,j) = 0.0
  end do
end do
!$OMP end do

do k = levels, 1, -1  ! this loop needs to work downwards for the calculation
                      ! of cfl_max in bm_cld to work correctly
  ! Indices of levels above and below, but not allowed to go out-of-bounds!
  km1 = max( k - 1, 1 )
  kp1 = min( k + 1, levels )

!$OMP do SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end

      ! Initialise values

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

      tl_lay(i,j,2)      = t_in(i,j,k)
      ql_lay(i,j,2)      = q_in(i,j,k)
      tdc_lay(i,j,2)     = max(tau_dec_bm(i,j,k),bm_tiny)
      inv_thm_lay(i,j,2) = 1.0/max(tau_hom_bm(i,j,k),bm_tiny)
      inv_tmp_lay(i,j,2) = 1.0/max(tau_mph_bm(i,j,k),bm_tiny)
      wvar_lay(i,j,2)    = max(bl_w_var(i,j,k),bm_tiny)
      dtldz_lay(i,j,2)   = ( t_in(i,j,kp1) - t_in(i,j,km1) )                   &
                         / ( z_theta(i,j,kp1) - z_theta(i,j,km1) )
      dqtdz_lay(i,j,2)   = ( q_in(i,j,kp1) - q_in(i,j,km1) )                   &
                         / ( z_theta(i,j,kp1) - z_theta(i,j,km1) )
      svar_ini(i,j,2) = max((0.5*bl_w_var(i,j,k)*tau_dec_bm(i,j,k))/           &
                            (1.0/tau_hom_bm(i,j,k)),bm_tiny)

    end do

    do i = tdims%i_start, tdims%i_end
      if ( l_mixing_ratio ) then
        call qsat_wat_mix(qsl_v(i),t_in(i,j,k),p_theta_levels(i,j,k))
        call qsat_mix(qsi_v(i),t_in(i,j,k),p_theta_levels(i,j,k))
      else
        call qsat_wat(qsl_v(i),t_in(i,j,k),p_theta_levels(i,j,k))
        call qsat(qsi_v(i),t_in(i,j,k),p_theta_levels(i,j,k))
      end if
    end do

    do i = tdims%i_start, tdims%i_end

      qsl_lay(i,j,2) = qsl_v(i)
      qsi_lay(i,j,2) = qsi_v(i)

      alphal = alphl * qsl_v(i) / (t_in(i,j,k) * t_in(i,j,k))
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
      sigx(i,2) = max( sigx(i,2), max( 0.01*alx*qsl_v(i)/max_sigmas,           &
                                alx*qsl_v(i)*min_sigx_ft*tmp/(tmp+qsl_v(i)) ) )

      mux  = alx*(q_in(i,j,k) - qsl_v(i))
      qnx(i)  = mux/(max_sigmas*sigx(i,2))

      ! store first guess of the uni-modal turbulence-based variance. This
      ! will be updated in bm_cld for cloudy grid points, taking into account
      ! iterative latent heat release.
      ! (note this actually stores standard deviation, not variance)
      svar_turb(i,j,k) = sigx(i,2)
    end do

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
        l_set_modes(i) = (k > 2) .and. (k < levels-2)                          &
                                 .and. (kez_inv(i,j,k) > 1)                    &
                                 .and. (kez_bottom(i,j,kez_inv(i,j,k)) < k)
        if ( l_set_modes(i) ) then
          ! Copy Properties of mode from above from model-level kez_top
          kk = kez_top(i,j,kez_inv(i,j,k))
          tlx = t_in(i,j,kk)*(p_theta_levels(i,j,k)/                           &
                              p_theta_levels(i,j,kk))**kappa
          if ( l_mixing_ratio ) then
            call qsat_wat_mix(qsl_lay(i,j,3),tlx,p_theta_levels(i,j,k))
            call qsat_mix(qsi_lay(i,j,3),tlx,p_theta_levels(i,j,k))
          else
            call qsat_wat(qsl_lay(i,j,3),tlx,p_theta_levels(i,j,k))
            call qsat(qsi_lay(i,j,3),tlx,p_theta_levels(i,j,k))
          end if
          tl_lay(i,j,3)      = tlx
          ql_lay(i,j,3)      = q_in(i,j,kk)
          wvar_lay(i,j,3)    = max(bl_w_var(i,j,kk),bm_tiny)
          tdc_lay(i,j,3)     = max(tau_dec_bm(i,j,kk),bm_tiny)
          inv_thm_lay(i,j,3) = 1.0/max(tau_hom_bm(i,j,kk),bm_tiny)
          kkm1 = max( kk - 1, 1 )
          kkp1 = min( kk + 1, levels )
          dtldz_lay(i,j,3)   = ( t_in(i,j,kkp1) - t_in(i,j,kkm1) )             &
                             / ( z_theta(i,j,kkp1) - z_theta(i,j,kkm1) )
          dqtdz_lay(i,j,3)   = ( q_in(i,j,kkp1) - q_in(i,j,kkm1) )             &
                             / ( z_theta(i,j,kkp1) - z_theta(i,j,kkm1) )
          svar_ini(i,j,3)    = max((0.5*bl_w_var(i,j,kk)*tau_dec_bm(i,j,kk))/  &
                                   (1.0/tau_hom_bm(i,j,kk)),bm_tiny)
          ! Copy Properties of mode from below from model-level kez_bottom
          kk = kez_bottom(i,j,kez_inv(i,j,k))
          tlx = t_in(i,j,kk)*(p_theta_levels(i,j,k)/                           &
                              p_theta_levels(i,j,kk))**kappa
          if ( l_mixing_ratio ) then
            call qsat_wat_mix(qsl_lay(i,j,1),tlx,p_theta_levels(i,j,k))
            call qsat_mix(qsi_lay(i,j,1),tlx,p_theta_levels(i,j,k))
          else
            call qsat_wat(qsl_lay(i,j,1),tlx,p_theta_levels(i,j,k))
            call qsat(qsi_lay(i,j,1),tlx,p_theta_levels(i,j,k))
          end if
          tl_lay(i,j,1)      = tlx
          ql_lay(i,j,1)      = q_in(i,j,kk)
          wvar_lay(i,j,1)    = max(bl_w_var(i,j,kk),bm_tiny)
          tdc_lay(i,j,1)     = max(tau_dec_bm(i,j,kk),bm_tiny)
          inv_thm_lay(i,j,1) = 1.0/max(tau_hom_bm(i,j,kk),bm_tiny)
          kkm1 = max( kk - 1, 1 )
          kkp1 = min( kk + 1, levels )
          dtldz_lay(i,j,1)   = ( t_in(i,j,kkp1) - t_in(i,j,kkm1) )             &
                             / ( z_theta(i,j,kkp1) - z_theta(i,j,kkm1) )
          dqtdz_lay(i,j,1)   = ( q_in(i,j,kkp1) - q_in(i,j,kkm1) )             &
                             / ( z_theta(i,j,kkp1) - z_theta(i,j,kkm1) )
          svar_ini(i,j,1)    = max((0.5*bl_w_var(i,j,kk)*tau_dec_bm(i,j,kk))/  &
                                   (1.0/tau_hom_bm(i,j,kk)),bm_tiny)
        end if
      end do  ! i = tdims%i_start, tdims%i_end
    end if  ! ( i_bm_ez_opt )

    do i = tdims%i_start, tdims%i_end

      if ( l_set_modes(i) ) then

        !--------------------------------------------------------------------
        ! Mode from above the inversion
        !--------------------------------------------------------------------

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

        ! Use local phase-relaxation time scale for level k
        inv_tmp_lay(i,j,3) = 1.0/max(tau_mph_bm(i,j,k),bm_tiny)
        mux = alx*(ql_lay(i,j,3) - qsl)
        qnx(i)             = max(qnx(i),mux/(max_sigmas*sigx(i,3)))

        !------------------------------ 2 -----------------------------------
        ! Mode from entrainment zone bottom
        !--------------------------------------------------------------------

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

        ! Use local phase-relaxation time scale for level k
        inv_tmp_lay(i,j,1) = 1.0/max(tau_mph_bm(i,j,k),bm_tiny)
        mux = alx*(ql_lay(i,j,1) - qsl)
        qnx(i)             = max(qnx(i),mux/(max_sigmas*sigx(i,1)))

      end if  ! ( l_set_modes(i) )

      lcld(i,j) = (qnx(i) > -1.0 .or. qcf(i,j,k) > 0.0)

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
  ! -- Section 5 - Store the i,j indices of k-levels that potentially       --
  ! --             contain cloud.                                           --
  ! --------------------------------------------------------------------------

  qc_points = 0

!$OMP do SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      if (lcld(i,j)) then
        qc_points = qc_points + 1
        idx(qc_points,1) = i
        idx(qc_points,2) = j
      end if
    end do ! i
  end do ! j
!$OMP end do NOWAIT


  ! --------------------------------------------------------------------------
  ! -- Section 6 - Call BM_CLD to calculate cloud water content,            --
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
                  qcl(1,1,k),qcf(1,1,k),cfl(1,1,k),cff(1,1,k),cf(1,1,k),       &
                  cfl_max(1,1),q(1,1,k),t(1,1,k),                              &
                  sskew(1,1,k),svar_turb(1,1,k),svar_bm(1,1,k),                &
                  idx,qc_points,l_mixing_ratio )

  end if ! Qc_points_if

end do !k loop
!$OMP end PARALLEL


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine bm_ctl
end module bm_ctl_mod
