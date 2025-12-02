! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Large-scale precipitation scheme. Generation of mixed-phase cloud
! by turbulent processes
module mphys_turb_gen_mixed_phase_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private ::                                        &
  ModuleName='MPHYS_TURB_GEN_MIXED_PHASE_MOD'

contains

subroutine mphys_turb_gen_mixed_phase( q_work, t_work, qcl_work, qcf_work,     &
                                       q_inc, qcl_inc, cfl_inc,  cf_inc,       &
                                       t_inc,  dqcl_mp, bl_levels,             &
                                       bl_w_var, cff_work, cfl_work, cf_work,  &
                                       q_n, cfl_n, cf_n, p_layer_centres,      &
                                       rhodz_dry, rhodz_moist, deltaz,         &
                                       qcl_mpt, tau_d, inv_prt, disprate,      &
                                       inv_mt, si_avg, dcfl_mp, sigma2_s ,     &
                                       qcf2_work, icenumber, snownumber  )

! Microphysics modules
use mphys_inputs_mod,      only: mp_dz_scal
use mphys_constants_mod,   only: cx, constp, mp_smallnum, mp_one_third,        &
     mprog_min
use thresholds,            only: ni_small, ns_small
use lsp_moments_mod,       only: lsp_moments
use mphys_ice_mod,         only: rhoi
use lsp_dif_mod,           only: air_conductivity0, air_diffusivity0, tcor1,   &
                                 tcor2, cpwr

! Stochastic physics
use stochastic_physics_run_mod, only:  l_rp2, rp_idx, mp_czero_rp

! General and constants modules
use gen_phys_inputs_mod,   only: l_mr_physics
use conversions_mod,       only: pi, zerodegc
use water_constants_mod,   only: lc, lf
use planet_constants_mod,  only: cp, r, repsilon, pref, rv, g

! Grid bounds module
use atm_fields_bounds_mod, only: tdims, tdims_l

use qsat_mod, only: qsat_mix, qsat_wat_mix

use mphys_inputs_mod,      only: l_casim, mp_czero, mp_tau_lim
use mphys_parameters,      only: snow_params, ice_params

use free_tracers_inputs_mod, only: l_wtrac
use wtrac_pc2_mod,           only: wtrac_pc2

! Dr Hook modules
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim

implicit none

! Purpose:
! Produces mixed phase cloud by subgrid turbulence, in order to
! improve the radiation budget and account for missing microphysical
! processes

! Method and paper reference:
! Uses the stochastic model from Field et al (2013), Q. J. R. Met. Soc
! "Mixed-phase clouds in a turbulent environment. Part 2: Analytic treatment"

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Description of Code:
!  Fortran90
!  This code is written to UMDP3 programming standards.

! Documentation: UMDP 29A.

! Subroutine arguments

integer, intent(in) :: bl_levels

real(kind=real_umphys), intent(in out) ::                                      &
                        q_work(   tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                  &
                                             1 : tdims%k_end )

real(kind=real_umphys), intent(in out) ::                                      &
                         qcl_work(tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                  &
                                             1 : tdims%k_end )

real(kind=real_umphys), intent(in out) ::                                      &
                         qcf_work(tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                  &
                                             1 : tdims%k_end )

real(kind=real_umphys), intent(in out) ::                                      &
                        t_work(  tdims%i_start : tdims%i_end,                  &
                                tdims%j_start : tdims%j_end,                   &
                                          1 : tdims%k_end )

real(kind=real_umphys), intent(in out) ::                                      &
                         qcf2_work(tdims%i_start : tdims%i_end,                &
                                 tdims%j_start : tdims%j_end,                  &
                                             1 : tdims%k_end )

real, intent(in) ::  icenumber(tdims_l%i_start:tdims_l%i_end,                  &
                               tdims_l%j_start:tdims_l%j_end,                  &
                               tdims_l%k_start:tdims_l%k_end)
!                    Ice number concentration from Casim
real, intent(in) ::  snownumber(tdims_l%i_start:tdims_l%i_end,                 &
                                tdims_l%j_start:tdims_l%j_end,                 &
                                tdims_l%k_start:tdims_l%k_end)
!                    Snow number concentration from Casim

real(kind=real_umphys), intent(in out) ::                                      &
                        cff_work( tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                  &
                                             1 : tdims%k_end )

real(kind=real_umphys), intent(in out) ::                                      &
                        cfl_work( tdims%i_start : tdims%i_end,                 &
                                 tdims%j_start : tdims%j_end,                  &
                                             1 : tdims%k_end )

real(kind=real_umphys), intent(in out) ::                                      &
                        cf_work( tdims%i_start : tdims%i_end,                  &
                                tdims%j_start : tdims%j_end,                   &
                                            1 : tdims%k_end )

real(kind=real_umphys), intent(in) :: q_n   ( tdims_l%i_start : tdims_l%i_end, &
                                              tdims_l%j_start : tdims_l%j_end, &
                                              tdims_l%k_start : tdims_l%k_end )

real(kind=real_umphys), intent(in) :: cfl_n ( tdims_l%i_start : tdims_l%i_end, &
                                              tdims_l%j_start : tdims_l%j_end, &
                                              tdims_l%k_start : tdims_l%k_end )

real(kind=real_umphys), intent(in) :: cf_n  ( tdims_l%i_start : tdims_l%i_end, &
                                              tdims_l%j_start : tdims_l%j_end, &
                                              tdims_l%k_start : tdims_l%k_end )

real(kind=real_umphys), intent(in) ::                                          &
                    bl_w_var(    tdims%i_start : tdims%i_end,                  &
                                 tdims%j_start : tdims%j_end,                  &
                                             1 : tdims%k_end )

real(kind=real_umphys), intent(in) ::                                          &
                    p_layer_centres( tdims%i_start : tdims%i_end,              &
                                     tdims%j_start : tdims%j_end,              &
                                                 0 : tdims%k_end )

real(kind=real_umphys), intent(in) ::                                          &
                    rhodz_dry(  tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end,                   &
                                            1 : tdims%k_end )

real(kind=real_umphys), intent(in) ::                                          &
                      rhodz_moist(tdims%i_start : tdims%i_end,                 &
                                  tdims%j_start : tdims%j_end,                 &
                                              1 : tdims%k_end )

real(kind=real_umphys), intent(in) ::                                          &
                      deltaz( tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end,                     &
                                          1 : tdims%k_end )

real(kind=real_umphys), intent(in out) ::                                      &
                        q_inc( tdims%i_start : tdims%i_end,                    &
                               tdims%j_start : tdims%j_end,                    &
                               tdims%k_start : tdims%k_end )

real(kind=real_umphys), intent(in out) ::                                      &
                        qcl_inc( tdims%i_start : tdims%i_end,                  &
                                 tdims%j_start : tdims%j_end,                  &
                                 tdims%k_start : tdims%k_end )

real(kind=real_umphys), intent(in out) ::                                      &
                        cfl_inc( tdims%i_start : tdims%i_end,                  &
                                 tdims%j_start : tdims%j_end,                  &
                                 tdims%k_start : tdims%k_end )

real(kind=real_umphys), intent(in out) ::                                      &
                        cf_inc( tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end,                   &
                                tdims%k_start : tdims%k_end )

real(kind=real_umphys), intent(in out) ::                                      &
                        t_inc( tdims%i_start : tdims%i_end,                    &
                               tdims%j_start : tdims%j_end,                    &
                                           1 : tdims%k_end )

real(kind=real_umphys), intent(in out) ::                                      &
                        dqcl_mp( tdims%i_start : tdims%i_end,                  &
                                tdims%j_start : tdims%j_end,                   &
                                            1 : tdims%k_end )

! Diagnostics
!                    qcl generate by turbulent mixed_phase
real(kind=real_umphys), intent(out) ::                                         &
                     qcl_mpt( tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end,                     &
                                          1 : tdims%k_end )

!                    Turbulent decorrelation timescale [s]
real(kind=real_umphys), intent(out) ::                                         &
                     tau_d(   tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end,                     &
                                          1 : tdims%k_end )

!                    Inverse Phase-Relaxation Timescale [s-1]
real(kind=real_umphys), intent(out) ::                                         &
                     inv_prt( tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end,                     &
                                          1 : tdims%k_end )

!                    Diagnosed turbulent dissipation rate [m2 s-3]
real(kind=real_umphys), intent(out) ::                                         &
                     disprate(tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end,                     &
                                          1 : tdims%k_end )

!                    Inverse Mixing Timescale [s-1]
real(kind=real_umphys), intent(out) ::                                         &
                     inv_mt(  tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end,                     &
                                          1 : tdims%k_end )

!                    Mean of subgrid supersaturation PDF
real(kind=real_umphys), intent(out) ::                                         &
                     si_avg(  tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end,                     &
                                          1 : tdims%k_end )

!                    Liquid cloud fraction diagnosed by mixed phase scheme
real(kind=real_umphys), intent(out) ::                                         &
                     dcfl_mp( tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end,                     &
                                          1 : tdims%k_end )

!                    Variance of subgrid supersaturation PDF
real(kind=real_umphys), intent(out) ::                                         &
                     sigma2_s(tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end,                     &
                                          1 : tdims%k_end )

real(kind=real_umphys) :: Ei
                         ! Saturated vapour pressure wrt ice [Pa]
real(kind=real_umphys) :: rhice            ! Relative humidity of ice [0-1]
real(kind=real_umphys) :: siw              ! The value of ice supersaturation at
                         ! water saturation []
real(kind=real_umphys) :: siw_lim          ! Limit of ice supersaturation
real(kind=real_umphys) :: aa
                         ! Thermodynamic term (for definition see Field (2013))
real(kind=real_umphys) :: b0
                         ! Thermodynamic term (for definition see Field (2013))
real(kind=real_umphys) :: Ai
                         ! Thermodynamic term (for definition see Field (2013))
real(kind=real_umphys) :: bi
                         ! Thermodynamic term (for definition see Field (2013))
real(kind=real_umphys) :: ka               ! temperature-corrected conductivity
real(kind=real_umphys) :: dv
                         ! temperature and pressure-corrected diffusivity
real(kind=real_umphys) :: fdist            ! frequency distribution of
real(kind=real_umphys) :: qv_excess        ! Excess moisture mixing ratio
real(kind=real_umphys) :: Sice             ! Saturation
real(kind=real_umphys) :: deltas           ! Local change of sigma_s
real(kind=real_umphys) :: fac              ! Local factors ...
real(kind=real_umphys) :: fac2             ! ... used in calculations
real(kind=real_umphys) :: four_root_sigmas ! local variable 4.0 * sqrt(sigma_s)
real(kind=real_umphys) :: dz_scal
                         ! Scaling factor applied to layer thickness
real(kind=real_umphys) :: t_limit
                         ! local temperature limit applied to mixed phase
                         ! calculation
real(kind=real_umphys) :: t_corr
                         ! temperature correction for constants
real(kind=real_umphys) :: p_corr           ! pressure correction for diffusivity
real(kind=real_umphys) :: tau_d_work
                         ! local working value of the turbulent decorrelation
                         ! timescale

! Number of bins in the mixed phase calculation (balance cost vs accuracy)
integer, parameter :: nbins_mp = 100 ! hard-wired to 100

integer :: grd_pts       ! Number of points in decomposed domain

integer :: ibin          ! Bin number

integer :: i             ! Loop counter in x direction
integer :: j             ! Loop counter in y direction
integer :: k             ! Loop counter in z direction

real(kind=real_umphys) ::                                                      &
        t_local2d(tdims%i_start : tdims%i_end,                                 &
                  tdims%j_start : tdims%j_end)
! Local value of temperature (K)

real(kind=real_umphys) ::                                                      &
        q_local2d(tdims%i_start : tdims%i_end,                                 &
                  tdims%j_start : tdims%j_end)
! Local value of humidity (kg/kg)

real(kind=real_umphys) ::                                                      &
        qsi_2d( tdims%i_start : tdims%i_end,                                   &
                tdims%j_start : tdims%j_end )

real(kind=real_umphys) ::                                                      &
        qsw_2d( tdims%i_start : tdims%i_end,                                   &
                tdims%j_start : tdims%j_end )

real(kind=real_umphys) :: mom1(   tdims%i_start : tdims%i_end)
! First moment of the ice particle size distribution

real(kind=real_umphys) ::                                                      &
        rho_air( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end )
! Air density (kg m-3)

real(kind=real_umphys) ::                                                      &
        rho_dry( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end )
! Dry air density (kg m-3)

real(kind=real_umphys) ::                                                      &
        cff_inv(  tdims%i_start : tdims%i_end,                                 &
                  tdims%j_start : tdims%j_end,                                 &
                              1 : tdims%k_end )
! 1/frozen cloud fraction

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

logical, parameter :: l_subgrid_cfl_mp_by_erosion = .false.
real :: ipcx, ipdx, spcx, spdx, imu, smu, lami, lams  ! casim parameters
real :: Gam_1_imu_id, Gam_1_imu, Gam_1_smu_sd, Gam_1_smu
real :: icenumber_cas(tdims%i_start : tdims%i_end,                             &
                      tdims%j_start : tdims%j_end,                             &
                                  1 :  tdims%k_end )
real :: snownumber_cas(tdims%i_start : tdims%i_end,                            &
                       tdims%j_start : tdims%j_end,                            &
                                  1 :  tdims%k_end )

character(len=*), parameter :: RoutineName='MPHYS_TURB_GEN_MIXED_PHASE'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!=============================================================================
! START OF PHYSICS
!=============================================================================

! If RP scheme is in use, set parameters to their perturbed values
if ( l_rp2 ) then
  mp_czero = mp_czero_rp(rp_idx)
end if

! Set the temperature and saturated ice water limits
t_limit = zerodegc - 0.01
siw_lim = mp_smallnum

grd_pts =  tdims%j_len  *                                                      &
           tdims%i_len

if (l_casim) then
  ipcx = ice_params%c_x
  ipdx = ice_params%d_x
  imu=ice_params%fix_mu
  spcx = snow_params%c_x
  spdx = snow_params%d_x
  smu=snow_params%fix_mu
  Gam_1_imu_id=gamma(1+imu+ipdx)
  Gam_1_imu=gamma(1+imu)
  Gam_1_smu_sd=gamma(1+smu+spdx)
  Gam_1_smu=gamma(1+smu)
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(k,j,i ) SHARED(tdims, icenumber_cas, snownumber_cas, icenumber,  &
!$OMP                        snownumber)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        icenumber_cas(i,j,k) = icenumber(i,j,k)
        snownumber_cas(i,j,k) = snownumber(i,j,k)
      end do
    end do
  end do
!$OMP end PARALLEL do
else
  ipcx = 0.0
  ipdx = 0.0
  spcx = 0.0
  spdx = 0.0
  imu=0.0
  smu=0.0
  lami=0.0
  lams=0.0
!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP private(k,j,i) SHARED(tdims, icenumber_cas, snownumber_cas)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        icenumber_cas(i,j,k) = 0.0
        snownumber_cas(i,j,k) = 0.0
      end do
    end do
  end do
!$OMP end PARALLEL do
end if

!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP private(k,j,i,q_local2d,t_local2d,rho_dry,rho_air,qsi_2d,qsw_2d,mom1,    &
!$OMP         rhice,siw,ei,t_corr,p_corr,dv,ka,bi,Ai,b0,aa,dz_scal,            &
!$OMP         tau_d_work,fac,four_root_sigmas,fac2,deltas,ibin,sice,           &
!$OMP         qv_excess,fdist,lami,lams)                                       &
!$OMP SHARED(tdims,dqcl_mp,qcl_mpt,tau_d,inv_prt,disprate,inv_mt,si_avg,       &
!$OMP        dcfl_mp,sigma2_s,bl_levels,q_work,t_work,rhodz_dry,deltaz,        &
!$OMP        l_mr_physics,rhodz_moist,cff_inv,cff_work,                        &
!$OMP        p_layer_centres,grd_pts,qcf_work,cx,repsilon,t_limit,bl_w_var,    &
!$OMP        pref,cp,constp,g,r,mp_dz_scal,siw_lim,cfl_work,q_n,cfl_n,cf_n,    &
!$OMP        qcl_inc,q_inc,t_inc,cfl_inc,cf_inc,qcl_work,cf_work, l_casim,     &
!$OMP        icenumber_cas, snownumber_cas, qcf2_work, ipcx, ipdx, spcx,spdx,  &
!$OMP        l_wtrac, wtrac_pc2, mp_czero, mp_tau_lim, Gam_1_imu_id, Gam_1_imu,&
!$OMP        Gam_1_smu_sd, Gam_1_smu, ni_small, ns_small)
!$OMP do SCHEDULE(STATIC)
do k=1, tdims%k_end
  do j=tdims%j_start, tdims%j_end
    do i=tdims%i_start, tdims%i_end
      dqcl_mp(i,j,k)  =  0.0
      qcl_mpt(i,j,k)  =  0.0
      tau_d(i,j,k)    = -2.0 ! To isolate points where the scheme is not used.
      inv_prt(i,j,k)  =  0.0
      disprate(i,j,k) =  0.0
      inv_mt(i,j,k)   =  0.0
      si_avg(i,j,k)   =  0.0
      dcfl_mp(i,j,k)  =  0.0
      sigma2_s(i,j,k) =  0.0
    end do
  end do
end do
!$OMP end do

! Loop over boundary layer levels
!$OMP do SCHEDULE(DYNAMIC)
do k = 1, bl_levels-1

  ! Store and modify values of q and T for mixed phase calculation
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end

      q_local2d(i,j) = q_work(i,j,k)
      t_local2d(i,j) = t_work(i,j,k)

      ! dry air density
      rho_dry(i,j) = rhodz_dry(i,j,k) / deltaz(i,j,k)

      if (l_mr_physics) then
         ! rho is the dry density
        rho_air(i,j) = rho_dry(i,j)
      else
         ! rho is the moist density
        rho_air(i,j) = rhodz_moist(i,j,k) / deltaz(i,j,k)

        ! Subsequent liquid cloud calculation uses mixing
        ! ratios. Convert moisture variable to a mixing ratio
        q_local2d(i,j) = rho_air(i,j)*q_local2d(i,j)/rho_dry(i,j)

      end if

      cff_inv(i,j,k) = 1.0 / max( cff_work(i,j,k), 0.001 )

    end do
  end do

  !  Calls to qsat_mix, qsat_wat_mix and lsp_moments.
  !  Done level-by-level for efficiency and generates the values
  !  of mom1, qsi_2d and qsw_2d

  ! Always request saturated values as mixing ratios
  ! for subsequent calculation of liquid cloud
  call qsat_mix(qsi_2d, t_local2d, p_layer_centres(:,:,k),                     &
                    tdims%j_len,tdims%i_len)

  call qsat_wat_mix(qsw_2d, t_local2d, p_layer_centres(:,:,k),                 &
                        tdims%j_len,tdims%i_len)

  ! Assuming the Field 'Generic' PSD and no ventilation effects here,
  ! to keep things simple.
  ! Note: this may not be consistent with assumptions in the microphysical
  ! deposition rate calculation in lsp_deposition

  do j = tdims%j_start, tdims%j_end


    if (l_casim) then
       ! For CASIM derive ice moment 1 by combining the snow and ice moments
      do i=tdims%i_start, tdims%i_end
        mom1(i) = 0.0
        if (icenumber_cas(i,j,k) > ni_small .and.                              &
             qcf2_work(i,j,k) > mprog_min) then
          lami=(ipcx*icenumber_cas(i,j,k)/qcf2_work(i,j,k)*                    &
                                        Gam_1_imu_id/Gam_1_imu)**(1/ipdx)
          mom1(i)=rho_air(i,j)*(icenumber_cas(i,j,k)/lami)
        end if
        if (snownumber_cas(i,j,k) > ns_small .and.                             &
             qcf_work(i,j,k) > mprog_min) then
          lams=(spcx*snownumber_cas(i,j,k)/qcf_work(i,j,k)*                    &
                                         Gam_1_smu_sd/Gam_1_smu)**(1/spdx)
          mom1(i)=mom1(i) + rho_air(i,j)*(snownumber_cas(i,j,k)/lams)
        end if
      end do
    else
       !Call lsp_moments row by row to help it play nicely with the interface
      call lsp_moments( tdims%i_len, rho_air(:,j), t_local2d(:,j),             &
           qcf_work(:,j,k), cff_inv(:,j,k), cx(84), mom1 )
    end if

    do i = tdims%i_start, tdims%i_end

       ! First statement in if test catches any -ve qv from pc2
      if (q_local2d(i,j)    > mp_smallnum .and.                                &
          t_local2d(i,j)    < t_limit  .and.                                   &
          bl_w_var(i,j,k)   > 1e-12 ) then

          ! N.B. Using approximate mixing ratio of
          ! vapour to get an RH
        rhice = q_local2d(i,j) / qsi_2d(i,j)

        ! N.B. This way of obtaining Siw and Ei
        ! is for qsw, qsi mixing ratios
        siw   = (qsw_2d(i,j) / qsi_2d(i,j)) - 1.0

        ei    = qsi_2d(i,j) * p_layer_centres(i,j,k) /                         &
              ( qsi_2d(i,j) + repsilon )

        ! Pressure and temperature corrections:
        t_corr = ( (t_local2d(i,j) / zerodegc)**cpwr) *                        &
                 ( tcor1 / ( t_local2d(i,j) + tcor2 ) )

        p_corr = ( pref / p_layer_centres(i,j,k) )

        dv = air_diffusivity0 * t_corr * p_corr

        ka = air_conductivity0 * t_corr

        bi = 1.0 / q_local2d(i,j) + (Lc+Lf)**2 /                               &
             (cp * rv * t_local2d(i,j) ** 2)

        Ai = 1.0 / (rhoi *(Lc+Lf)**2 / (ka*rv*t_local2d(i,j)**2) +             &
                  rhoi * rv * t_local2d(i,j) / (Ei*Dv))

        if (.not. l_casim) then
          b0 = 4.0 * pi * constp(35) * rhoi * Ai / rho_air(i,j)
        else
           ! 1.0 is a representation of capacitance and used
           ! to replace constp(35), which is not available with
           ! CASIM
          b0 = 4.0 * pi * 1.0 * rhoi * Ai / rho_air(i,j)
        end if

        aa = ( g / (r*t_local2d(i,j) ) *                                       &
             ( (Lc+Lf)*r / (cp*rv*t_local2d(i,j))-1.0))

        dz_scal = mp_dz_scal * deltaz(i,j,k)

        disprate(i,j,k) = 2.0 * bl_w_var(i,j,k)**1.5 /                         &
                          (mp_czero * dz_scal)

        ! Limit tau_d to specified limit
        tau_d_work = 2.0 * bl_w_var(i,j,k) / (disprate(i,j,k)*mp_czero)

        inv_prt(i,j,k)  = bi * b0 * mom1(i)

        inv_mt(i,j,k) = (disprate(i,j,k) / dz_scal**2)**mp_one_third

        ! Add factor to speed up calculations
        fac = 1.0 / (inv_prt(i,j,k) + inv_mt(i,j,k))

        sigma2_s(i,j,k) = 0.5 *aa**2 * sqrt(bl_w_var(i,j,k))                   &
                        * dz_scal * fac

        si_avg(i,j,k) = (rhice-1.0) * inv_mt(i,j,k) * fac

        if ( siw > siw_lim .and. aa > mp_smallnum .and.                        &
             tau_d_work <=  mp_tau_lim   ) then

          tau_d(i,j,k) = tau_d_work

          four_root_sigmas = 4 * sqrt(sigma2_s(i,j,k))
          fac2   = sqrt(2 * pi * sigma2_s(i,j,k) )
          deltas = four_root_sigmas / nbins_mp

          do ibin = 0, nbins_mp-1

            sice      = siw + ibin * deltas
            qv_excess = qsi_2d(i,j) * ibin * deltas

            if ( abs(sice-si_avg(i,j,k)) <= four_root_sigmas ) then

              fdist = exp(-(sice-si_avg(i,j,k))**2/(2*sigma2_s(i,j,k)))        &
                      / fac2
            else
              fdist     = 0.0
            end if

            qcl_mpt(i,j,k) = qcl_mpt(i,j,k) + qv_excess * fdist * deltas
            dcfl_mp(i,j,k) = dcfl_mp(i,j,k) + fdist * deltas

          end do ! ibin

          ! The calculated qcl_mpt is a mixing ratio.
          ! Convert to specific quantity if necessary
          if (.not. l_mr_physics) then
            qcl_mpt(i,j,k) = rho_dry(i,j)*qcl_mpt(i,j,k)/rho_air(i,j)
          end if

          if ( cfl_work(i,j,k) < mp_smallnum .or.                              &
              .not. l_subgrid_cfl_mp_by_erosion ) then

            ! Ensure that qcl_inc + qcl_mpt < q_n
            qcl_mpt(i,j,k) = min( qcl_mpt(i,j,k),                              &
                                  q_n(i,j,k) - qcl_inc(i,j,k) )

            qcl_inc(i,j,k) = qcl_inc(i,j,k) + qcl_mpt(i,j,k)
            q_inc(i,j,k)   = q_inc(i,j,k)   - qcl_mpt(i,j,k)
            t_inc(i,j,k)   = t_inc(i,j,k)   + Lc * qcl_mpt(i,j,k) / cp

            cfl_inc(i,j,k) = min( cfl_inc(i,j,k) + dcfl_mp(i,j,k),             &
                                       1.0 - cfl_n(i,j,k)       )

            cf_inc(i,j,k)  = min( cf_inc(i,j,k) + dcfl_mp(i,j,k),              &
                                       1.0 - cf_n(i,j,k))

            qcl_work(i,j,k) = qcl_work(i,j,k) + qcl_mpt(i,j,k)
            q_work(i,j,k)   = q_work(i,j,k)   - qcl_mpt(i,j,k)
            t_work(i,j,k)   = t_work(i,j,k)   + Lc * qcl_mpt(i,j,k) / cp

            cfl_work(i,j,k) = min( cfl_work(i,j,k) + dcfl_mp(i,j,k), 1.0)

            cf_work(i,j,k)  = min( cf_work(i,j,k) + dcfl_mp(i,j,k), 1.0)

            if (l_wtrac) wtrac_pc2%q_cond(i,j,k) = qcl_mpt(i,j,k)

          else ! cfl_work == 0 etc

            dqcl_mp(i,j,k) = qcl_mpt(i,j,k) - qcl_work(i,j,k)

          end if  ! cfl_work == 0 etc

        else if (tau_d_work >  mp_tau_lim) then

          tau_d(i,j,k) = -1.0 ! Flag for tau_d limit exceeded

        end if    ! Siw > 0 etc

      end if       ! qlocal2d > 0

    end do         ! i
  end do           ! j
end do             ! k
!$OMP end do NOWAIT
!$OMP end PARALLEL

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine mphys_turb_gen_mixed_phase
end module mphys_turb_gen_mixed_phase_mod
