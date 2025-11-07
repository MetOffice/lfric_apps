! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module bm_calc_tau_mod

implicit none

character(len=*), parameter, private ::                                        &
  ModuleName='BM_CALC_TAU_MOD'

contains

subroutine bm_calc_tau( q, theta, exner_theta, qcf, bl_levels,                 &
                        cff, p_theta_levels, bl_w_var,                         &
                        elm, mix_len_bm, rho_dry, rho_moist,                   &
                        icenumber, snownumber,                                 &
                        tau_dec, tau_hom, tau_mph, qcf2)

!Microphysics modules
use mphys_constants_mod,   only: cx, constp, mp_smallnum, mp_one_third
use lsp_moments_mod,       only: lsp_moments
use mphys_ice_mod,         only: rhoi
use lsp_dif_mod,           only: air_conductivity0, air_diffusivity0, tcor1,   &
                                 tcor2, cpwr
use cloud_inputs_mod,      only: i_bm_ez_opt, i_bm_ez_entpar
! Stochastic physics
use stochastic_physics_run_mod, only:  l_rp2, rp_idx, mp_czero_rp

!General and constants modules
use gen_phys_inputs_mod,   only: l_mr_physics
use conversions_mod,       only: pi
use water_constants_mod,   only: lc, lf, tm
use planet_constants_mod,  only: cp, repsilon, rv, pref

! Grid bounds module
use atm_fields_bounds_mod, only: tdims,tdims_l,tdims_s

use qsat_mod, only: qsat_mix

! Dr Hook modules
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim


use mphys_inputs_mod, only: l_casim, mp_czero, mp_tau_lim

use mphys_parameters,     only: snow_params, ice_params
use thresholds,           only: qi_small, ni_small, qs_small, ns_small

implicit none

! Purpose:
! Calculate time scales and variances for the bi-modal cloud scheme:
!     - decorrelation time scale              (tau_dec_bm)
!     - homogenisation time scale             (tau_hom_bm)
!     - phase-relaxation time scale           (tau_mph_bm)
! This routine is similar to mphys_turb_gen_mixed_phase, but uses the actual
! blended mixing length scale as defined in the boundary layer scheme. It also
! just calculates the time scales and is called in atmos_physics2 to make sure
! the time scales and vertical velocity variance are updated before calling
! the cloud scheme.

! Method and paper reference:
! Uses the stochastic model from Field et al (2013), Q. J. R. Met. Soc
! "Mixed-phase clouds in a turbulent environment. Part 2: Analytic treatment"

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud

! Description of Code:
!  Fortran90
!  This code is written to UMDP3 programming standards.

! Documentation: UMDP 39

! Subroutine arguments

integer, intent(in) :: bl_levels
!                    Number of boundary layer levels

real, intent(in) ::  qcf(tdims_l%i_start:tdims_l%i_end,                        &
                         tdims_l%j_start:tdims_l%j_end,                        &
                         tdims_l%k_start:tdims_l%k_end)
!                    snow content (kg per kg air)
real, intent(in) ::  qcf2(tdims_l%i_start:tdims_l%i_end,                       &
                         tdims_l%j_start:tdims_l%j_end,                        &
                         tdims_l%k_start:tdims_l%k_end)
!                    Cloud ice content (kg per kg air)

real, intent(in) ::  cff(tdims_l%i_start:tdims_l%i_end,                        &
                         tdims_l%j_start:tdims_l%j_end,                        &
                         tdims_l%k_start:tdims_l%k_end)
!                    Cloud ice fraction

real, intent(in) ::  icenumber(tdims_l%i_start:tdims_l%i_end,                  &
                               tdims_l%j_start:tdims_l%j_end,                  &
                               tdims_l%k_start:tdims_l%k_end)
!                    Ice number concentration from Casim
real, intent(in) ::  snownumber(tdims_l%i_start:tdims_l%i_end,                 &
                                tdims_l%j_start:tdims_l%j_end,                 &
                                tdims_l%k_start:tdims_l%k_end)
!                    Snow number concentration from Casim

real, intent(in) ::  q(tdims_l%i_start:tdims_l%i_end,                          &
                       tdims_l%j_start:tdims_l%j_end,                          &
                       tdims_l%k_start:tdims_l%k_end)
!                    Specific humidity (kg/kg)

real, intent(in) ::  theta(tdims_s%i_start:tdims_s%i_end,                      &
                           tdims_s%j_start:tdims_s%j_end,                      &
                           tdims_s%k_start:tdims_s%k_end)
!                    Potential temperature (K)

real, intent(in) ::  exner_theta(tdims_s%i_start:tdims_s%i_end,                &
                                 tdims_s%j_start:tdims_s%j_end,                &
                                 tdims_s%k_start:tdims_s%k_end)
!                    Exner on theta levels

real, intent(in) ::  bl_w_var(tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end,                     &
                                          1 : tdims%k_end )
!                    BL w-variance diagnostic (m/s)^2

real, intent(in) ::  elm(tdims%i_start : tdims%i_end,                          &
                         tdims%j_start : tdims%j_end,                          &
                                     1 : bl_levels )
!                    Mixing length for momentum (on levels k+1)

real, intent(in) ::  mix_len_bm(tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end,                   &
                                            1 : tdims%k_end )
!                    Alternative mixing-length optionally used (on levels k)

real, intent(in) ::  p_theta_levels(tdims_s%i_start:tdims_s%i_end,             &
                                    tdims_s%j_start:tdims_s%j_end,             &
                                    tdims_s%k_start:tdims_s%k_end)
!                    Pressure on theta levels (Pa).

real, intent(in) ::  rho_dry(tdims%i_start : tdims%i_end,                      &
                             tdims%j_start : tdims%j_end,                      &
                                             tdims%k_end-1)
!                    Dry density on theta levels (not top) (kg/m3)


real, intent(in) ::  rho_moist(tdims%i_start : tdims%i_end,                    &
                               tdims%j_start : tdims%j_end,                    &
                                               tdims%k_end-1)
!                    Wet density on theta levels (not top) (kg/m3)

real, intent(out) :: tau_dec(tdims%i_start : tdims%i_end,                      &
                             tdims%j_start : tdims%j_end,                      &
                                         1 : tdims%k_end )
!                    Decorrelation Timescale (s)

real, intent(out) :: tau_hom(tdims%i_start : tdims%i_end,                      &
                             tdims%j_start : tdims%j_end,                      &
                                         1 : tdims%k_end )
!                    Homogenisation Timescale (s)

real, intent(out) :: tau_mph(tdims%i_start : tdims%i_end,                      &
                             tdims%j_start : tdims%j_end,                      &
                                         1 : tdims%k_end )
!                    Phase-Relaxation Timescale (s)

real :: mix_len          ! Turbulent mixing-length (m)
real :: ei               ! Saturated vapour pressure wrt ice (Pa)
real :: b0               ! Thermodynamic term (for definition see Field (2013))
real :: ai               ! Thermodynamic term (for definition see Field (2013))
real :: bi               ! Thermodynamic term (for definition see Field (2013))
real :: ka               ! temperature-corrected conductivity
real :: dv               ! temperature and pressure-corrected diffusivity
real :: t_corr           ! temperature correction for constants
real :: p_corr           ! pressure correction for diffusivity
real :: qsi
real :: qkw              ! 3 x vertical velocity variance
real :: dissp            ! Eddy dissipation rate

integer :: i             ! Loop counter in x direction
integer :: j             ! Loop counter in y direction
integer :: k             ! Loop counter in z direction

real, parameter :: b1 = 24.0

real :: q_local2d(tdims%i_start : tdims%i_end,                                 &
                  tdims%j_start : tdims%j_end)
! Local value of humidity (kg/kg)
real :: t_local2d(tdims%i_start : tdims%i_end,                                 &
                  tdims%j_start : tdims%j_end)
! Local value of temperature (K)
real :: qcf_local2d(tdims%i_start : tdims%i_end,                               &
                    tdims%j_start : tdims%j_end)

real :: mom1(tdims%i_start : tdims%i_end)
! First moment of the ice particle size distribution

real :: rho_air(tdims%i_start : tdims%i_end,                                   &
                tdims%j_start : tdims%j_end )
! Air density (kg m-3)

real :: cff_inv(tdims%i_start : tdims%i_end,                                   &
                tdims%j_start : tdims%j_end,                                   &
                                bl_levels )
! 1/frozen cloud fraction

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

real :: ipcx, ipdx, spcx, spdx, imu, smu, lami, lams  ! casim parameters
real :: Gam_1_imu_id, Gam_1_imu, Gam_1_smu_sd, Gam_1_smu

real :: icenumber_cas(tdims%i_start : tdims%i_end,                             &
                      tdims%j_start : tdims%j_end,                             &
                                  1 :  tdims%k_end )
real :: snownumber_cas(tdims%i_start : tdims%i_end,                            &
                       tdims%j_start : tdims%j_end,                            &
                                  1 :  tdims%k_end )


character(len=*), parameter :: RoutineName='BM_CALC_TAU'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!=============================================================================
! START OF PHYSICS
!=============================================================================

! If RP scheme is in use, set parameters to their perturbed values
if ( l_rp2 ) then
  mp_czero = mp_czero_rp(rp_idx)
end if

!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED(tdims,tau_mph,tau_dec,tau_hom)                                    &
!$OMP private(k,j,i)
do k = 1, tdims%k_end
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      ! Initialise time scales
      tau_mph(i,j,k) = 1.0e-6
      tau_dec(i,j,k) = 1.0e-6
      tau_hom(i,j,k) = 1.0e-6
    end do
  end do
end do
!$OMP end PARALLEL do

if (l_casim) then
  ipcx=ice_params%c_x
  ipdx=ice_params%d_x
  imu=ice_params%fix_mu
  spcx=snow_params%c_x
  spdx=snow_params%d_x
  smu=snow_params%fix_mu
  Gam_1_imu_id=gamma(1+imu+ipdx)
  Gam_1_imu=gamma(1+imu)
  Gam_1_smu_sd=gamma(1+smu+spdx)
  Gam_1_smu=gamma(1+smu)

!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none) private(i, j, k)              &
!$OMP SHARED(tdims, icenumber_cas, snownumber_cas, icenumber, snownumber)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        icenumber_cas(i,j,k)=icenumber(i,j,k)
        snownumber_cas(i,j,k)=snownumber(i,j,k)
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  ipcx=0.0
  ipdx=0.0
  imu=0.0
  spcx=0.0
  spdx=0.0
  smu=0.0
  lami=0.0
  lams=0.0
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none) private(i, j, k)              &
!$OMP SHARED(tdims, icenumber_cas, snownumber_cas)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        icenumber_cas(i,j,k)=0.0
        snownumber_cas(i,j,k)=0.0
      end do
    end do
  end do
!$OMP end PARALLEL do
end if


!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP private(k,j,i,t_corr,ka,bi,ai,ei,dv,p_corr,b0,qsi,qkw,dissp,             &
!$OMP         lami,lams,mix_len)                                               &
!$OMP SHARED(tdims,q,theta,rho_dry,rho_moist,cff_inv,l_mr_physics,tau_mph,     &
!$OMP        cff,qcf,q_local2d,qcf_local2d,t_local2d,rho_air,repsilon,         &
!$OMP        cp,constp,pref,p_theta_levels,cx,bl_w_var,                        &
!$OMP        elm,mix_len_bm,tau_dec,tau_hom,exner_theta,bl_levels,             &
!$OMP        l_casim, qcf2, icenumber_cas, snownumber_cas, spcx,spdx,ipcx,ipdx,&
!$OMP        mp_czero,mp_tau_lim, ni_small, gam_1_imu_id, gam_1_imu, imu,      &
!$OMP        ns_small, gam_1_smu_sd, gam_1_smu, smu, qi_small, qs_small,       &
!$OMP        i_bm_ez_opt, mom1)
do k = 1, bl_levels
  do j = tdims%j_start, tdims%j_end
!$OMP do SCHEDULE(STATIC)
    do i = tdims%i_start, tdims%i_end

      ! Calculate decorrelation and homogenisation time scales

      if (bl_w_var(i,j,k) > 1.0e-12 .and. k < bl_levels) then

        qkw = sqrt(3.0*bl_w_var(i,j,k))

        if ( i_bm_ez_opt == i_bm_ez_entpar ) then
          ! Use bespoke mixing-length for bimodal scheme when available
          mix_len = mix_len_bm(i,j,k)
        else
          ! NOTE: elm is defined on levels k+1
          mix_len = elm(i,j,k+1)
        end if

        ! Calculate eddy dissipation rate as in Mellor-Yamada(1982) scheme.

        dissp = (qkw**3.0)/(b1*mix_len)

        ! Calculate the decorrelation and homogenisation time scales
        tau_dec(i,j,k) = 2.0*bl_w_var(i,j,k)/(dissp*mp_czero)
        tau_hom(i,j,k) = ((mix_len**2.0)/dissp)**(mp_one_third)

        if ( .not. i_bm_ez_opt == i_bm_ez_entpar ) then
          ! Impose max limit on timescales only if not using the entraining
          ! parcel-based mixing-length.  These limits are disabled for the
          ! entraining parcel model because they break the scaling of the
          ! turbulent fluctuations with mix_len in an ad-hoc way.
          tau_dec(i,j,k) = min( tau_dec(i,j,k), mp_tau_lim )
          tau_hom(i,j,k) = min( tau_hom(i,j,k), mp_tau_lim )
        end if

      end if

      ! Store and modify values of q and T for mixed phase calculation

      q_local2d(i,j) = q(i,j,k)
      t_local2d(i,j) = theta(i,j,k)*exner_theta(i,j,k)
      qcf_local2d(i,j) = qcf(i,j,k)
      if (l_casim) qcf_local2d(i,j) = qcf_local2d(i,j) + qcf2(i,j,k)

      ! dry air density
      if (l_mr_physics) then
        ! rho is the dry density
        rho_air(i,j) = rho_dry(i,j,k)
      else
        ! rho is the moist density
        rho_air(i,j) = rho_moist(i,j,k)

        ! Subsequent liquid cloud calculation uses mixing
        ! ratios. Convert moisture variable to a mixing ratio
        q_local2d(i,j) = rho_air(i,j)*q_local2d(i,j)/rho_dry(i,j,k)

      end if

      cff_inv(i,j,k) = 1.0 / max( cff(i,j,k), 0.001 )

    end do
!$OMP end do
  end do

  do j = tdims%j_start, tdims%j_end

    !Call lsp_moments row by row to help it play nicely with the interface
    if (l_casim) then
      ! For CASIM derive ice moment 1 by combining the snow and ice moments
!$OMP do SCHEDULE(STATIC)
      do i=tdims%i_start, tdims%i_end
        mom1(i)=0.0
        if ((icenumber_cas(i,j,k) > ni_small) .and.                            &
                                                 (qcf2(i,j,k) > qi_small)) then
          lami=(ipcx*icenumber_cas(i,j,k)/qcf2(i,j,k)*                         &
                                      Gam_1_imu_id/Gam_1_imu)**(1/ipdx)
          mom1(i)=rho_air(i,j)*(icenumber_cas(i,j,k)/lami) * (1.0+imu)
        end if
        if ((snownumber_cas(i,j,k) > ns_small) .and.                           &
                                                  (qcf(i,j,k) > qs_small)) then
          lams=(spcx*snownumber_cas(i,j,k)/qcf(i,j,k)*                         &
                                      Gam_1_smu_sd/Gam_1_smu)**(1/spdx)
          mom1(i)=mom1(i) +                                                    &
                rho_air(i,j)*(snownumber_cas(i,j,k)/lams) * (1.0+smu)
        end if

      end do
!$OMP end do
    else
!$OMP SINGLE
      call lsp_moments( tdims%i_len, rho_air(:,j), t_local2d(:,j),             &
           qcf_local2d(:,j), cff_inv(:,j,k), cx(84), mom1 )
!$OMP END SINGLE
    end if !l_casim
!$OMP do SCHEDULE(STATIC)
    do i = tdims%i_start, tdims%i_end

      ! First statement in if test catches any -ve qv from pc2
      if (q_local2d(i,j) > mp_smallnum .and. qcf_local2d(i,j) > 0.0 ) then

        call qsat_mix(qsi, t_local2d(i,j), p_theta_levels(i,j,k))

        ei    = qsi * p_theta_levels(i,j,k) /                                  &
              ( qsi + repsilon )

        ! Pressure and temperature corrections:
        t_corr = ( (t_local2d(i,j)/tm)**cpwr) *                                &
                   ( tcor1/( t_local2d(i,j)+tcor2 ) )
        p_corr = ( pref / p_theta_levels(i,j,k) )

        dv = air_diffusivity0 * t_corr * p_corr
        ka = air_conductivity0 * t_corr

        bi = 1.0 / q_local2d(i,j) + (Lc+Lf)**2 /                               &
             (cp * rv * t_local2d(i,j) ** 2)

        ai = 1.0 / (rhoi *(Lc+Lf)**2 / (ka*rv*t_local2d(i,j)**2) +             &
                  rhoi * rv * t_local2d(i,j) / (ei*dv))

        if (.not. l_casim) then
          b0 = 4.0 * pi * constp(35) * rhoi * Ai / rho_air(i,j)
        else
          ! 1.0 is a representation of capacitance and used
          ! to replace constp(35), which is not available with
          ! CASIM
          b0 = 4.0 * pi * 1.0 * rhoi * Ai / rho_air(i,j)
        end if

        tau_mph(i,j,k)  = min(1.0/ max(bi * b0 * mom1(i), 1.0e-6), mp_tau_lim)

      end if       ! qlocal2d > 0

    end do         ! i
!$OMP end do
  end do           ! j
end do             ! k
!$OMP end PARALLEL
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine bm_calc_tau
end module bm_calc_tau_mod
