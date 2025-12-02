! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
module parcel_ascent_mod
use um_types, only: r_bl

implicit none
private
! Calculates parcel ascent
public parcel_ascent
character(len=*), parameter, private :: ModuleName='PARCEL_ASCENT_MOD'

contains
!
! Subroutine Interface:
subroutine parcel_ascent (npnts, nSCMDpkgs,                                    &
                          nlcl, k_plume,                                       &
                          l_dilute, L_SCMDiags,                                &
                          sl_plume, qw_plume,                                  &
                          t, Tl, qcl, qcf, qw, t_dens_env,                     &
                          p_theta_lev, exner_theta_levels,                     &
                          z_theta, entrain_fraction,                           &
                          T_parc, t_parc_dil, ql_parc,                         &
                          buoyancy, buoyancy_dil,                              &
                          env_svl, par_svl, par_svl_dil,                       &
                          denv_bydz, dpar_bydz, dqsatdz, qv_parc_out)

use cv_run_mod, only:                                                          &
  plume_water_load, dil_plume_water_load, qlmin, fac_qsat

use cv_param_mod, only:                                                        &
  qlcrit

use planet_constants_mod, only:                                                &
  r => rd_bl, repsilon => repsilon_bl, c_virtual => c_virtual_bl, g => g_bl,   &
  ls => ls_bl, lsrcp => lsrcp_bl, lcrcp => lcrcp_bl, gamma_dry => grcp_bl

use water_constants_mod, only: lc => lc_bl, tm => tm_bl
use bl_option_mod, only: zero, one, one_half
use model_domain_mod,      only: model_type, mt_single_column


use s_scmop_mod,           only: default_streams,                              &
                                 t_inst,d_wet,scmdiag_conv
use nlsizes_namelist_mod,  only: model_levels

use qsat_mod, only: qsat, qsat_mix

use gen_phys_inputs_mod, only: l_mr_physics
!$ use omp_lib,          only: omp_get_max_threads

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! ------------------------------------------------------------------------------
! Description:
!   This routine calculates a parcel ascent from k_plume.
!   An undilute ascent is always calculated. If the option l_dilute is set
!   to .true. a dilute ascent is also calculated. The dilute ascent uses
!   the entrainment rate held in entrain_fraction to mix in environmental air.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.
!
! ------------------------------------------------------------------------------

! Subroutine arguments

integer, intent(in) ::                                                         &
  npnts                & ! Number of points
 ,nSCMDpkgs              ! SCM  - No of diagnostics packages

integer, intent(in) ::                                                         &
  nlcl(npnts)          & ! Lifting condensation level
 ,k_plume(npnts)         ! Starting model level for plume ascent

logical, intent(in) ::                                                         &
  l_dilute             & ! .true. if a dilute parcel ascent also required.
 ,L_SCMDiags(nSCMDpkgs)  ! SCM - Logicals for diagnostics packages

real(kind=r_bl), intent(in) ::                                                 &
  t(npnts,model_levels)   & ! Temperature on model levels (K)
 ,Tl(npnts,model_levels)  & ! Liquid water temperature on model lev (K)
 ,qcl(npnts,model_levels) & ! cloud liquid water on model levels (kg/kg)
 ,qcf(npnts,model_levels) & ! cloud ice water on model levels (kg/kg)
 ,qw(npnts,model_levels)  & ! total water (kg/kg)
 ,t_dens_env(npnts,model_levels)
                             ! Density potential temperature of environment (K)

real(kind=r_bl), intent(in) ::                                                 &
  p_theta_lev(npnts,model_levels)        & ! Pressure on theta levels (Pa)
 ,exner_theta_levels(npnts,model_levels) & ! Exner Pressure on theta levels
 ,z_theta(npnts,model_levels)            & ! Height of theta levels  (m)
 ,entrain_fraction(npnts,model_levels)     ! fraction of environmental air
                                               ! to mix with parcel

real(kind=r_bl), intent(in out) ::                                             &
  sl_plume(npnts)      & ! SL at start of plume
 ,qw_plume(npnts)      & ! total water at plume start (kg/kg)
 ,t_parc_dil(npnts,model_levels)     ! Dilute Parcel temperature (K)

real(kind=r_bl), intent(out) ::                                                &
  t_parc(npnts,model_levels)       & ! Parcel temperature  (K)
 ,ql_parc(npnts,model_levels)      & ! Parcel water content
 ,buoyancy(npnts,model_levels)     & ! Parcel buoyancy  (K)
 ,buoyancy_dil(npnts,model_levels) & ! Dilute parcel buoyancy (K)
 ,env_svl(npnts,model_levels)      & ! Density (virtual) static energy
                                         ! over CP for layer.
 ,par_svl(npnts,model_levels)      & ! Density (virtual) static energy
                                         ! over CP of parcel for level.
 ,par_svl_dil(npnts,model_levels)  & ! Density (virtual) static energy
                                         ! over CP of parcel for level.
 ,denv_bydz(npnts, model_levels)   & ! Gradient of density potential
                                         ! temperature in the environment.
 ,dpar_bydz(npnts, model_levels)   & ! Gradient of density potential
 ,dqsatdz(npnts, model_levels)       ! dqsat/dz along an adiabat undilute
                                         ! parcel

real(kind=r_bl), intent(out),   optional          ::                           &
 qv_parc_out(npnts, model_levels)        ! parcel water undilute plume

! Local variables

integer ::                                                                     &
  ii,k            &  ! loop counters
 ,omp_block       &  ! for openmp blocking
 ,ij              &  ! for indexing over openmp block
 ,nij             &  ! final loop value in openmp block
 ,npnts_omp          ! number of points in openmp block

real(kind=r_bl) ::                                                             &
  q_liq_env       &  ! Condensed water content of environment.
 ,dq_sat_env      &  ! DQSAT/DT for environment
 ,lrcp_const      &  ! lc or lc+lf over cp
 ,lrcp_const_env  &  ! lc or lc+lf over cp
 ,lrcp_const_parc &  ! lc or lc+lf over cp
 ,l_const         &  ! lc or lc+lf
 ,l_const_env     &  ! lc or lc+lf
 ,dz              &  ! layer depth
 ,dtdz            &  ! temperature gradient along undilute parcel ascent
 ,z_pr            &  ! used in estimating th_ref at next level
 ,th_par          &  ! theta value for parcel
 ,dq_sat_par      &  ! dqsat/dT for parcel
 ,dq_sat_par_dil  &  ! dqsat/dT for parcel
 ,temp_parc       &  ! average temperature of parcel after entrainment
 ,q_vap_parc      &  ! Vapour content of undilute parcel
 ,q_liq_parc      &  ! Liquid water content of undilute parcel
 ,qcl_parc        &  ! parcel qcl - dilute parcel cal
 ,qcf_parc        &  ! parcel qcf - dilute parcel cal
 ,ql_remove          ! ql removed from plume


! parcel calculation

real(kind=r_bl) ::                                                             &
  t_ref(npnts)                         & ! reference temperature
 ,th_ref(npnts)                        & ! reference potential temperature
 ,th_par_km1(npnts)                    & ! theta refernce for level below
 ,qsat_lev(npnts, model_levels)        & ! qsat for reference temperature
 ,qsat_env(npnts, model_levels)        & ! qsat for environment temperature
 ,t_dens_parc(npnts, model_levels)       ! Density potential temperature
                                         ! of parcel.
real(kind=r_bl) ::                                                             &
  max_qw_plume                           ! maximum plume water content

! Arrays added for dilute parcel calculation

real(kind=r_bl) ::                                                             &
  t_ref_dil(npnts)                    & ! dilute parcel reference temperature
 ,th_ref_dil(npnts)                   & ! reference potential temperature
 ,th_par_km_dil(npnts)                & ! reference potential temperature 2nd
 ,qsat_lev_dil(npnts)                 & ! qsat for dilute parcel
 ,qw_parc(npnts,model_levels)         & ! parcel total water undilute plume
 ,qv_parc(npnts,model_levels)         & ! parcel water undilute plume (array not
                                        ! essential but may require in future)
 ,sl_parc(npnts,model_levels)         & ! parcel SL undilute
 ,t_dens_parc_dil(npnts,model_levels) & ! dilute parcel t_dens
 ,ql_parc_dil(npnts,model_levels)       ! dilute parcel liquid water

logical ::                                                                     &
  l_keep_water    ! if true keeps water loading in plume
                  ! false removed if water exceeds 1g/kg


character(len=*), parameter ::  RoutineName = 'PARCEL_ASCENT'



! Model constants

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! 1.0 Initialisation
!-----------------------------------------------------------------------

l_keep_water  = .false. ! water loading restricted to 1g/kg

! Initialise parcel reference theta and value for level below

omp_block = npnts
!$ omp_block = ceiling(real(npnts)/omp_get_max_threads())

!$OMP PARALLEL DEFAULT(SHARED)                                                 &
!$OMP private(ii, ij, nij, k, npnts_omp, lrcp_const, lrcp_const_env,           &
!$OMP         l_const, dq_sat_env, dq_sat_par, q_liq_parc, q_liq_env,          &
!$OMP         z_pr, lrcp_const_parc, ql_remove, max_qw_plume,                  &
!$OMP         q_vap_parc, th_par, temp_parc, qcl_parc, qcf_parc,               &
!$OMP         l_const_env, dq_sat_par_dil, dz, dtdz)
!$OMP do SCHEDULE(STATIC)
do ii=1, npnts
  th_ref(ii) = tl(ii,k_plume(ii))/exner_theta_levels(ii,k_plume(ii))
  th_par_km1(ii) = th_ref(ii)
end do
!$OMP end do NOWAIT

if (l_dilute) then     ! dilute parcel ascent
!$OMP do SCHEDULE(STATIC)
  do ii=1, npnts
    th_ref_dil(ii)    = th_ref(ii)
    th_par_km_dil(ii) = th_ref_dil(ii)
  end do
!$OMP end do NOWAIT
end if

!-----------------------------------------------------------------------
! 2.0 Parcel ascent
!-----------------------------------------------------------------------
! Parcel starts from level k_plume and is lifted upwards to top of model.
! Dilute parcel ascent - mix in environmental air above lifting
! condensation level
!
!-----------------------------------------------------------------------
! Calculate parcel water by linearising qsat about the Parcel's
! temperature extrapolated up to the next grid_level.
!-----------------------------------------------------------------------

! Level loop calculating parcel ascent

if ( l_mr_physics ) then
!$OMP do SCHEDULE(STATIC)
  do  k = 1,model_levels
    call qsat_mix(qsat_env(:,k),t(:,k),p_theta_lev(:,k),npnts)
  end do
!$OMP end do
else
!$OMP do SCHEDULE(STATIC)
  do  k = 1,model_levels
    call qsat(qsat_env(:,k),t(:,k),p_theta_lev(:,k),npnts)
  end do
!$OMP end do
end if

!$OMP do SCHEDULE(STATIC)
do ij = 1, npnts, omp_block
  nij = min(ij+omp_block-1,npnts)
  npnts_omp = nij - ij + 1
  do  k = 1,model_levels

    ! Require t_ref on all point for qsat call

    do ii=ij, nij
      t_ref(ii)     = th_ref(ii)*exner_theta_levels(ii,k)
    end do

    if ( l_mr_physics ) then
      call qsat_mix(qsat_lev(ij:nij,k),t_ref(ij:nij),                          &
                        p_theta_lev(ij:nij,k),npnts_omp)
    else
      call qsat(qsat_lev(ij:nij,k),t_ref(ij:nij),p_theta_lev(ij:nij,k),        &
                    npnts_omp)
    end if

    if (l_dilute) then     ! dilute parcel ascent

      do ii=ij, nij
        t_ref_dil(ii) = th_ref_dil(ii)*exner_theta_levels(ii,k)
      end do

      if ( l_mr_physics ) then
        call qsat_mix(qsat_lev_dil(ij:nij),t_ref_dil(ij:nij),                  &
                          p_theta_lev(ij:nij,k),npnts_omp)
      else
        call qsat(qsat_lev_dil(ij:nij),t_ref_dil(ij:nij),                      &
                      p_theta_lev(ij:nij,k),npnts_omp)
      end if

    end if                 ! dilute parcel

    ! Undilute parcel calculation always required

    do ii=ij, nij

      if (T_ref(ii) >  tm) then
        lrcp_const = lcrcp
        l_const    = lc
      else
        lrcp_const = lsrcp
        l_const    = ls
      end if
      if (t(ii,k) >  tm) then
        lrcp_const_env = lcrcp
        l_const_env    = lc
      else
        lrcp_const_env = lsrcp
        l_const_env    = ls
      end if

      dq_sat_env = repsilon*l_const_env*qsat_env(ii,k)/(r*t(ii,k)**2)
      dq_sat_par = repsilon*l_const    *qsat_lev(ii,k)/(r*T_ref(ii)**2)

      q_liq_parc = max( zero, ( qw_plume(ii) - qsat_lev(ii,k)                  &
                -dq_sat_par*( sl_plume(ii)-gamma_dry*z_theta(ii,k)-T_ref(ii) ) &
                                   ) / (one+lrcp_const*dq_sat_par) )

      q_liq_env  = max( zero, ( qw(ii,k) - qsat_env(ii,k)                      &
                  -dq_sat_env*( tl(ii,k)               - t(ii,k) )             &
                                   ) / (one+Lrcp_const_env*dq_sat_env) )
      !
      ! add on the difference in the environment's ql as calculated by the
      ! UM cloud scheme (using some RH_CRIT value) and what it
      ! would be If RH_CRIT=1. This then imitates partial condensation
      ! in the parcel.
      !
      ql_parc(ii,k) = q_liq_parc + qcl(ii,k) + qcf(ii,k) - q_liq_env

      t_parc(ii,k) = sl_plume(ii)-gamma_dry*z_theta(ii,k) + lrcp_const         &
                    *ql_parc(ii,k)


      ! May need to recalculate if T_parc is > Tm and T_ref < Tm

      if (T_ref(ii) <= tm .and. T_parc(ii,k) >  tm) then

        ! recalculate using corrected latent heats
        lrcp_const_parc = lcrcp

        q_liq_parc = max( zero, ( qw_plume(ii) - qsat_lev(ii,k)                &
             -dq_sat_par*( sl_plume(ii)-gamma_dry*z_theta(ii,k)-T_ref(ii) )    &
                            ) / (one+lrcp_const_parc*dq_sat_par) )

        ql_parc(ii,k) = q_liq_parc + qcl(ii,k) + qcf(ii,k)- q_liq_env

        ! revised at parcel calculation

        t_parc(ii,k)=sl_plume(ii)-gamma_dry*z_theta(ii,k)                      &
                                         +lrcp_const_parc*ql_parc(ii,k)

      end if       ! test on T_ref and t_parc

      ! Limit this from going to zero or negative
      t_parc(ii,k) = max(t_parc(ii,k), 1.0e-12_r_bl)

      ! Add water removal from undilute plume
      if (plume_water_load == 1) then

        ! water removed from parcel after condensation if > 0.001 kg/kg
        if (ql_parc(ii,k) > 0.001_r_bl) then
          ql_remove = ql_parc(ii,k) -  0.001_r_bl
          ql_parc(ii,k) = 0.001_r_bl
          qw_plume(ii) = qw_plume(ii) - ql_remove
          ! Also adjust sl_plume as well as altered energy
          sl_plume(ii) = sl_plume(ii) + lrcp_const*ql_remove
        end if

      else if (plume_water_load == 2) then

        ! water remaining depends on qsat environment
        max_qw_plume = fac_qsat*qsat_env(ii,k)
        max_qw_plume = min (qlcrit,max_qw_plume)  ! max
        max_qw_plume = max (qlmin,max_qw_plume)   ! min value
        if (ql_parc(ii,k) > max_qw_plume) then
          ql_remove = ql_parc(ii,k) - max_qw_plume
          ql_parc(ii,k) = max_qw_plume
          qw_plume(ii) = qw_plume(ii) - ql_remove
          ! Also adjust sl_plume as well as altered energy
          sl_plume(ii) = sl_plume(ii) + lrcp_const*ql_remove

        end if

      end if    ! end test on plume_water_load

      q_vap_parc=qw_plume(ii)-ql_parc(ii,k)
      qv_parc(ii,k) = q_vap_parc

      t_dens_parc(ii,k)=t_parc(ii,k)*(one+                                     &
           c_virtual*q_vap_parc-ql_parc(ii,k))

      ! calculate t_ref for next level
      if (k >  1 .and. k <   model_levels-1) then
        z_pr = (z_theta(ii,k+1)-z_theta(ii,k))/(z_theta(ii,k)-z_theta(ii,k-1))
        th_par = t_parc(ii,k)/exner_theta_levels(ii,k)
        th_ref(ii) = th_par*(one+z_pr) - th_par_km1(ii)*z_pr

        ! Check sensible value otherwise set to previous reference value
        ! Problems can occur near top of model where calculation are nolonger
        ! important.
        if (th_ref(ii) <= zero) then
          th_ref(ii) = th_par_km1(ii)
        end if
        if (th_par > zero) then
          th_par_km1(ii) = th_par
        end if
      end if

      !-----------------------------------------------------------------------
      ! Dilute parcel ascent
      !-----------------------------------------------------------------------
      if (l_dilute) then

        if (k <= nlcl(ii)) then  ! Dilute parcel same as undilute
                                 ! parcel ascent (no entrainment)

          sl_parc(ii,k) = sl_plume(ii)
          qw_parc(ii,k) = qw_plume(ii)
          t_parc_dil(ii,k)     = t_parc(ii,k)
          ql_parc_dil(ii,k)    = ql_parc(ii,k)
          t_dens_parc_dil(ii,k)= t_dens_parc(ii,k)
          th_ref_dil(ii)       = th_ref(ii)
          th_par_km_dil(ii)    = th_par_km1(ii)

        else                      ! Dilute parcel ascent now required

          if (t_ref_dil(ii) >  tm) then
            lrcp_const = lcrcp
            l_const    = lc
          else
            lrcp_const = lsrcp
            l_const    = ls
          end if

          !---------------------------------------------------------------------
          ! Dilute parcel
          !---------------------------------------------------------------------
          ! Mix in entrain_fraction from environmental air from level below and
          ! raise this to current level.
          ! Assume mix in fraction of mass from environment.
          ! Estimate parcel properties after mixing air from environment with
          ! parcel. Temperature given approximately by average

          temp_parc = (t_parc_dil(ii,k-1)                                      &
                                 + entrain_fraction(ii,k)*t(ii,k-1))           &
                          /(one+entrain_fraction(ii,k))


          qw_parc(ii,k) = (qw_parc(ii,k-1) +                                   &
                                entrain_fraction(ii,k)*qw(ii,k-1))             &
                           /(one+entrain_fraction(ii,k))

          qcl_parc = (ql_parc_dil(ii,k-1)   +                                  &
                                entrain_fraction(ii,k)*qcl(ii,k-1))            &
                           /(one+entrain_fraction(ii,k))

          qcf_parc = (zero     +                                               &
                                entrain_fraction(ii,k)*qcf(ii,k-1))            &
                           /(one+entrain_fraction(ii,k))

          ! All condensed water either ice or liquid based on t_ref
          sl_parc(ii,k) = temp_parc - lrcp_const*(qcl_parc+qcf_parc)           &
                                  +gamma_dry*z_theta(ii,k-1)

          dq_sat_par_dil = repsilon*l_const*qsat_lev_dil(ii)                   &
                                      /(r*t_ref_dil(ii)**2)

          q_liq_parc = max( zero, ( qw_parc(ii,k) - qsat_lev_dil(ii)           &
        -dq_sat_par_dil*( sl_parc(ii,k)-gamma_dry*z_theta(ii,k)-t_ref_dil(ii)) &
                                     ) / (one+lrcp_const*dq_sat_par_dil) )

          ! add on the dIfference in the environment's ql as calculated by the
          ! UM cloud scheme (using some RH_CRIT value) and what it
          ! would be If RH_CRIT=1. This then imitates partial condensation
          ! in the parcel.
          !
          ql_parc_dil(ii,k) = q_liq_parc + qcl(ii,k) + qcf(ii,k) - q_liq_env

          t_parc_dil(ii,k) = sl_parc(ii,k)-gamma_dry*z_theta(ii,k)             &
                                        +lrcp_const*ql_parc_dil(ii,k)

          ! May need to recalculate if T_parc is > Tm and T_ref < Tm

          if (t_ref_dil(ii) <= tm .and. t_parc_dil(ii,k) >  tm) then

            ! recalculate using corrected latent heats
            lrcp_const_parc = lcrcp

            q_liq_parc = max( zero, ( qw_parc(ii,k) - qsat_lev_dil(ii)         &
         -dq_sat_par_dil*(sl_parc(ii,k)-gamma_dry*z_theta(ii,k)-t_ref_dil(ii)) &
                           ) / (one+lrcp_const_parc*dq_sat_par_dil) )

            ql_parc_dil(ii,k) = q_liq_parc + qcl(ii,k)                         &
                                                + qcf(ii,k)- q_liq_env

            ! revised at parcel calculation

            t_parc_dil(ii,k)=sl_parc(ii,k)-gamma_dry*z_theta(ii,k)             &
                                     +lrcp_const_parc*ql_parc_dil(ii,k)

          end if   ! test on t_ref

          q_vap_parc=qw_parc(ii,k)-ql_parc_dil(ii,k)

          ! Water loading in plume  ( option 0 water remains in plume)

          if (dil_plume_water_load == 1) then

            ! water removed from parcel after condesation if > 0.001 kg/kg
            if (ql_parc_dil(ii,k) > 0.001_r_bl) then
              ql_parc_dil(ii,k) = 0.001_r_bl
              qw_parc(ii,k) = q_vap_parc + ql_parc_dil(ii,k)
            end if

          else if (dil_plume_water_load == 2) then

            ! water remaining depends on qsat environment
            max_qw_plume = one_half*qsat_env(ii,k)
            max_qw_plume = min (qlcrit,max_qw_plume)  ! max
            max_qw_plume = max (qlmin,max_qw_plume)   ! min value
            if (ql_parc_dil(ii,k) > max_qw_plume) then
              ql_parc_dil(ii,k) = max_qw_plume
              qw_parc(ii,k) = q_vap_parc + ql_parc_dil(ii,k)
            end if

          end if    ! end test on dil_plume_water_load

          t_dens_parc_dil(ii,k)=t_parc_dil(ii,k)*                              &
                          (one+c_virtual*q_vap_parc-ql_parc_dil(ii,k))

          ! calculate dilute t_ref for next level
          if (k >  1 .and. k <   model_levels-1) then
            z_pr = (z_theta(ii,k+1)-z_theta(ii,k))                             &
                                        /(z_theta(ii,k)-z_theta(ii,k-1))

            th_par = t_parc_dil(ii,k)/exner_theta_levels(ii,k)
            th_ref_dil(ii) = th_par*(one+z_pr) - th_par_km_dil(ii)*z_pr
            ! Check new reference sensible
            if (th_ref_dil(ii) <= zero) then
              th_ref_dil(ii) = th_par_km_dil(ii)
            end if
            if (th_par > zero) then
              th_par_km_dil(ii) = th_par
            end if
          end if      ! k level test

        end if   ! test on LCL
      end if    ! test on L_dilute

    end do     ! Loop over ii

  end do      ! level loop
end do
!$OMP end do

!$OMP do SCHEDULE(STATIC)
do  k = 1,model_levels
  do ii=1, npnts
    buoyancy(ii,k) = t_dens_parc(ii,k) - t_dens_env(ii,k)

    env_svl(ii,k)  = t_dens_env(ii,k)  + gamma_dry*z_theta(ii,k)

    par_svl(ii,k)  = t_dens_parc(ii,k) + gamma_dry*z_theta(ii,k)

    if (L_dilute) then
      buoyancy_dil(ii,k) = t_dens_parc_dil(ii,k) - t_dens_env(ii,k)
      par_svl_dil(ii,k)  = t_dens_parc_dil(ii,k) + gamma_dry*z_theta(ii,k)
    end if
  end do     ! Loop over ii
end do      ! level loop
!$OMP end do

!$OMP do SCHEDULE(STATIC)
do  k = 1,model_levels
  ! Gradient calculations

  if (k >= 2) then

    if ( l_mr_physics ) then
      call qsat_mix(qsat_lev(:,k),T_parc(:,k),p_theta_lev(:,k),npnts)
    else
      call qsat(qsat_lev(:,k),T_parc(:,k),p_theta_lev(:,k),npnts)
    end if

    do ii=1,npnts
      !-------------------------------------------------------------
      ! Find vertical gradients in parcel and environment SVL
      ! (using values from level below (i.e. K-1)).
      !-------------------------------------------------------------

      dz = z_theta(ii,k) - z_theta(ii,k-1)

      dpar_bydz(ii,k) = (par_svl(ii,k) - par_svl(ii,k-1))/dz

      denv_bydz(ii,k) = (env_svl(ii,k) - env_svl(ii,k-1))/dz

      !-----------------------------------------------------------------
      ! Temperature gradient and qsat/dz along an undilute parcel ascent
      !-----------------------------------------------------------------

      dtdz = (T_parc(ii,k) - T_parc(ii,k-1))/dz

      if (T_parc(ii,k) >  tm) then
        l_const=lc
      else
        l_const=ls
      end if

      dqsatdz(ii,k) = (qsat_lev(ii,k)/(r*T_parc(ii,k)))                        &
                            * (dtdz*repsilon*l_const/T_parc(ii,k) + g)

    end do    ! ii loop

  else

    do ii=1,npnts
      dpar_bydz(ii,k) = zero
      denv_bydz(ii,k) = zero
      dqsatdz(ii,k) = zero
    end do    ! ii loop

  end if   ! test on k

end do      ! level loop
!$OMP end do NOWAIT

if (present (qv_parc_out)) then
!$OMP do SCHEDULE(STATIC)
  do  k = 1,model_levels
    qv_parc_out(:,k) = qv_parc(:,k)
  end do
!$OMP end do
end if
!$OMP end PARALLEL

!-----------------------------------------------------------------------
! SCM diagnostics from convective diagnosis - only if required
!-----------------------------------------------------------------------
! Dilute parcel information not passed back to higher level routine
!-----------------------------------------------------------------------
!-------------------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine parcel_ascent
end module parcel_ascent_mod
