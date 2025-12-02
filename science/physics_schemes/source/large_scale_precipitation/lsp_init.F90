! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Initialisation of variables
! Subroutine Interface:
module lsp_init_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_INIT_MOD'

contains

subroutine lsp_init(                                                           &
  points,                                                                      &
                                     ! Number of points
  timestep, timestep_mp,                                                       &
                                     ! Timesteps
  t, p, cttemp,                                                                &
                                     ! Temperature and pressure
  deltaz, rhodz_dry, rhodz_moist,                                              &
                                           ! Air density information
  rho, rhor,                                                                   &
  dhi, dhir,                                                                   &
                                     ! Tstep and layer thick. ratios
  qcf, qcf2, qrain, qgraup, rainrate,                                          &
                                                   ! Water contents
  qcf_agg, qcf_cry, qcf_tot, frac_agg,                                         &
  qs, qsl, esi, esw,                                                           &
                                     ! Saturation water quantities
  cf_transfer_diag, cfl_transfer_diag,                                         &
  cff_transfer_diag, rf_transfer_diag,                                         &
                                     ! Cloud and rain fraction
                                     ! transfer diagnostics
  snow_cry, snow_agg, snowt_cry,                                               &
                                     ! Precipitation rates
  snowt_agg, rainratet, graupratet,                                            &
  lheat_correc_liq, lheat_correc_ice,                                          &
                                           ! Latent heat corrections
  corr, corr2, rocor,                                                          &
                                     ! Fall speed and diffusivity
                                     ! corrections
  agg_nofall, cry_nofall, rain_nofall, graup_nofall,                           &
                                     ! Fraction of precip not falling out
  tcg,tcgi,tcgc,tcgci,tcgg,tcggi,                                              &
                                           ! Ice PSD intercepts
  cficekeep, vm_cry, vm_agg, l_use_agg_vt, vm_used,                            &
       ! Variables for different fallspeeds with generic psd
  graut_psdep, graut_psacw, psacw_prev,                                        &
       ! Variables required for graupel autoconversion
  rainfrac, precfrac_k )

!Use in reals in lsprec precision, both microphysics related and general atmos
use lsprec_mod,         only: t_scaling, m0, qcf0, t_agg_min, cx, constp,      &
                              tcor1, tcor2, zerodegc, lc, lf, cp, r,           &
                              repsilon, recip_epsilon, cpwr
use lsprec_mod,         only: zero, half, one

  ! Microphysics modules- logicals and integers
use mphys_inputs_mod,    only: l_psd, l_mcr_qrain, l_mcr_qcf2, l_mcr_qgraup,   &
                               l_diff_icevt, l_mcr_precfrac
use mphys_bypass_mod,    only: l_crystals

! General and Atmospheric Modules- logicals and integers
use gen_phys_inputs_mod,  only: l_mr_physics

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,             only: real_lsprec

use qsat_mod, only: qsat, qsat_mix, qsat_wat, qsat_wat_mix

! Dr Hook Modules
use yomhook,             only: lhook, dr_hook
use parkind1,            only: jprb, jpim

! Mathematical modules
use vectlib_mod,   only: powr_v => powr_v_interface,                           &
                         exp_v  => exp_v_interface

! Large scale precip modules
use lsp_moments_mod, only: lsp_moments

implicit none

! Purpose:
!   Perform initialization of variables required for the
!   microphysics scheme.

! Method:
!   Set variables to their initial values.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation


! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.



! Subroutine Arguments


integer, intent(in) ::                                                         &
  points             ! Number of points

real (kind=real_lsprec), intent(in) ::                                         &
  timestep,                                                                    &
                     ! Model physics timestep / s
  timestep_mp,                                                                 &
                     ! Timestep of each microphysics iteration / s
  t(points),                                                                   &
                     ! Temperature / K
  p(points),                                                                   &
                     ! Pressure / N m-2
  qcf(points),                                                                 &
                     ! Ice aggregate content / kg kg-1
  qcf2(points),                                                                &
                     ! Ice crystal content / kg kg-1
  cttemp(points),                                                              &
                     ! Cloud top temperature / K
  rainrate(points),                                                            &
                       ! Rainrate into layer / kg m-2 s-1
  deltaz(points),                                                              &
                     ! Layer thickness / m
  rhodz_dry(points),                                                           &
                         ! Air density*deltaz for dry air based on
                         ! non-hydrostatic assumption / kg m-2
  rhodz_moist(points),                                                         &
                         ! Air density*deltaz for moist air based on
  cficekeep(points)
                     ! input ice cloud fraction

real (kind=real_lsprec), intent(in out) ::                                     &
  qgraup(points),                                                              &
                     ! Graupel content / kg kg-1
  snow_cry(points),                                                            &
                        ! Ice crystal precip into layer / kg m-2 s-1
  snow_agg(points)  ! Ice agg. precip into layer / kg m-2 s-1

real (kind=real_lsprec), intent(in out) ::                                     &
  frac_agg(points)
                     ! Fraction of ice that is aggregates (dimensionless)

real (kind=real_lsprec), intent(out) ::                                        &
  snowt_cry(points),                                                           &
                        ! Ice crystal precip out of layer / kg m-2 s-1
  snowt_agg(points),                                                           &
                        ! Ice agg. precip out of layer / kg m-2 s-1
  rainratet(points),                                                           &
                        ! Rain precip rate out of layer / kg m-2 s-1
  graupratet(points),                                                          &
                        ! Graupel precip rate out of layer/ kg m-2 s-1
  qcf_agg(points),                                                             &
                        ! Ice aggregate mixing ratio / kg kg-1
  qcf_cry(points),                                                             &
                        ! Ice crystal mixing ratio / kg kg-1
  qcf_tot(points)
                        ! Total ice content / kg kg-1

real (kind=real_lsprec), intent(in out) ::                                     &
  qrain(points)
                        ! Rain water content / kg kg-1
real (kind=real_lsprec), intent(out) ::                                        &
  qs(points),                                                                  &
                        ! Saturated humidity wrt ice / kg kg-1
  qsl(points),                                                                 &
                        ! Saturated humidity wrt liquid / kg kg-1
  rho(points),                                                                 &
                     ! Air density / kg m-3
  rhor(points),                                                                &
                     ! 1 / air density / m-3 kg
  esi(points),                                                                 &
                     ! Vapour pressure wrt ice / N m-2
  esw(points),                                                                 &
                     ! Vapour pressure wrt liquid / N m-2
  lheat_correc_liq(points),                                                    &
                               ! Latent heat correction for
                               ! liquid (no units)
  lheat_correc_ice(points),                                                    &
                               ! Latent heat correction for
                               ! ice (no units)
  dhi(points),                                                                 &
                        ! microphysics timestep / deltaz / s m-1
  dhir(points),                                                                &
                        ! deltaz / microphysic timestep / m s-1
  corr(points),                                                                &
                        ! Fall speed correction factor (no units)
  corr2(points),                                                               &
                        ! Diffusivity correction factor (no units)
  rocor(points),                                                               &
                        ! sqrt(rho*corr*corr2) (no units)
  agg_nofall(points),                                                          &
                        ! Fraction of the ice aggregate mass not falling out
  cry_nofall(points),                                                          &
                        ! Fraction of the ice crystal mass not falling out
  rain_nofall(points),                                                         &
                        ! Fraction of the rain-mass that is not falling out
  graup_nofall(points),                                                        &
                        ! Fraction of the graupel-mass that is not falling out
  tcg(points),                                                                 &
                        ! Temperature dependent aggregate PSD
                        ! intercept factor (no units)
  tcgi(points),                                                                &
                        ! 1 / tcg (no units)
  tcgc(points),                                                                &
                        ! Temperature dependent crystal PSD
                        ! intercept factor (no units)
  tcgci(points),                                                               &
                        ! 1 / tcgc (no units)
  tcgg(points),                                                                &
                        ! Temperature dependent graupel PSD
                        ! intercept factor (no units)
  tcggi(points),                                                               &
                        ! 1 / tcgg (no units)

  cf_transfer_diag(points),                                                    &
                        ! Diagnostic for change in cloud frac
  cfl_transfer_diag(points),                                                   &
                        ! Diagnostic for change in liq cloud frac
  cff_transfer_diag(points),                                                   &
                        ! Diagnostic for change in ice cloud frac
  rf_transfer_diag(points),                                                    &
                          ! Diagnostic for change in rain fraction
  vm_cry(points),                                                              &
                   ! Mass-weighted mean fallspeed crystal parameters
  vm_agg(points),                                                              &
                   ! Mass-weighted mean fallspeed aggregate parameters
  vm_used(points),                                                             &
                   ! Selected mass-weighted mean fallspeed
  graut_psdep(points),                                                         &
                   ! Snow deposition rate for graupel autoconversion
  graut_psacw(points),                                                         &
                   ! Aggregate-Snow riming rate for graupel autoconversion
  psacw_prev(points)
                   ! Saved aggregate-snow riming rate for use in rate
                   ! calculations


logical, intent(out) ::                                                        &
  l_use_agg_vt(points)
                 ! At each point defines which branch of the
! fallspeed relation is used. If .true. then use aggregate fallspeed
! parameters; else use crystal parameters.

! Sub-grid fraction of rain (and optionally graupel)
real (kind=real_lsprec), intent(in out) :: rainfrac(points)

! Prognostic precipitation fraction at level k;
! initial value used to initialise rain-fraction if used.
real (kind=real_lsprec), intent(in) :: precfrac_k(points)

! Local Variables

real (kind=real_lsprec) ::                                                     &
  corr_in(points),                                                             &
  corr2_in(points),                                                            &
  rocor_in(points),                                                            &
  frac_agg_out(points),                                                        &
  tcg_in(points),                                                              &
  tcgc_in(points),                                                             &
  tcgi_in(points),                                                             &
  tcgci_in(points)

integer ::                                                                     &
  i              ! Loop counter

real (kind=real_lsprec), parameter ::                                          &
  rho1 = one
                                   ! Ref. air density / kg m-3

! Variables for use when different fallspeed relations are used
! for crystals and aggregates with the generic psd
real (kind=real_lsprec) ::                                                     &
  m_bic_dic(points),                                                           &
         ! Psd moment that determines vertical ice mass flux when
! crystal fallspeed parameters are used
    m_bi_di(points)
           ! Psd moment that determines vertical ice mass flux when
! aggregate fallspeed parameters are used

real (kind=real_lsprec) ::                                                     &
  cfice(points),                                                               &
         ! ice cloud fraction
  cficei(points)
         ! reciprocal ice cloud fraction

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_INIT'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i = 1, points

    !-----------------------------------------------
    ! Set fluxes to zero
    !-----------------------------------------------
  snowt_cry(i)  = zero
  snowt_agg(i)  = zero
  rainratet(i)  = zero
  graupratet(i) = zero

   !--------------------------------------------------------------
   ! Set transfer variables (only used in lsp_ice and not actually
   ! coupled to STASH diagnostics) to zero
   !--------------------------------------------------------------
  cf_transfer_diag(i)  = zero
  cfl_transfer_diag(i) = zero
  cff_transfer_diag(i) = zero
  rf_transfer_diag(i)  = zero

   !----------------------------------------------
   ! Initialise fallspeed variables
   !----------------------------------------------
  vm_cry(i) = zero
  vm_agg(i) = zero
  vm_used(i) = zero
  l_use_agg_vt(i) = .false.

   !----------------------------------------------
   ! Initialise process rates required for
   ! graupel autoconversion
   !----------------------------------------------
  graut_psdep(i) = zero
  graut_psacw(i) = zero
  psacw_prev(i)  = zero

end do

    !-----------------------------------------------
    ! Set mixing ratios of graupel to zero if not used.
    ! If not done then (with high optimisation) this can be any
    ! old rubbish from memory and you get d_qgraup_dt being
    ! calculated as nan-nan in later routines, which causes a crash
    !-----------------------------------------------

if (.not. l_mcr_qgraup) then

  do i = 1, points

    qgraup(i) = zero

  end do

end if  ! .not. l_mcr_qgraup

if (l_mcr_qcf2) then
      ! If l_mcr_qcf2 is true then there are two ice/snow
      ! prognostics (qcf and qcf2), so copy them to
      ! qcf_cry and qcf_agg

  do i = 1, points

    qcf_cry(i) = qcf2(i)
    qcf_agg(i) = qcf(i)

  end do

else if (.not. l_crystals) then
      ! All ice is placed in a single ice category (aggregates).

  do i = 1, points

    frac_agg(i) = one
    qcf_cry(i)  = zero
    qcf_agg(i)  = qcf(i)
    snow_cry(i) = zero

  end do

else
      ! Split the one ice/snow
      ! prognostic (qcf) diagnostically into ice crystals
      ! (qcf_cry) and snow aggregates (qcf_agg)

  do i = 1, points

    frac_agg(i) = -t_scaling*max((t(i)-cttemp(i)),zero)                        &
                  *max(qcf(i)*qcf0,zero)
  end do
  call exp_v(points,frac_agg,frac_agg_out)
  do i = 1, points
    frac_agg(i) = max(one-frac_agg_out(i) , zero)
        ! Allocate ice content to crystals and aggregates
    qcf_cry(i) = qcf(i) * (one-frac_agg(i))
    qcf_agg(i) = qcf(i) * frac_agg(i)

        ! Assume falling snow is partitioned into crystals and
        ! aggregates. Snow_agg contains total snow on input
    snow_cry(i) = snow_agg(i) * (one-frac_agg(i))
    snow_agg(i) = snow_agg(i) * frac_agg(i)

  end do

end if ! l_mcr_qcf2

    !-----------------------------------------------
    ! Calculate total ice content
    !-----------------------------------------------
do i = 1, points

  qcf_tot(i) = qcf_cry(i) + qcf_agg(i)

end do

    !-----------------------------------------------
    ! If rain is a diagnostic, convert flux (kg m-2 s-1)
    ! to mass (kg kg-1)
    !-----------------------------------------------
if (.not. l_mcr_qrain) then

  ! Rain is a diagnostic quantity

  do i = 1, points

    if (rainrate(i)  >   zero) then

      if (l_mr_physics) then

        ! Mixing ratio formulation
        qrain(i) = rainrate(i) * timestep / rhodz_dry(i)

      else ! l_mr_physics

        ! Specific humidity formulation
        qrain(i) = rainrate(i) * timestep / rhodz_moist(i)

      end if  ! l_mr_physics

    else    ! rainrate > 0

      qrain(i) = zero

    end if  ! rainrate > 0

  end do ! points

end if  ! .not. l_mcr_qrain

    !-----------------------------------------------
    ! Calculate saturation specific humidities
    !-----------------------------------------------

if (l_mr_physics) then
  ! Qsat with respect to ice
  call qsat_mix(qs,t,p,points)
  ! Qsat with respect to liquid water
  call qsat_wat_mix(qsl,t,p,points)
else
  ! Qsat with respect to ice
  call qsat(qs,t,p,points)
  ! Qsat with respect to liquid water
  call qsat_wat(qsl,t,p,points)
end if

!When working in 32-bit, a Cray compiler bug breaks PROC comparability.
!Conditionally using NOVECTOR makes this go away. Review upon upgrade.
#if defined(LSPREC_32B) && defined (CRAYFTN_VERSION) && (CRAYFTN_VERSION <8004000)
!DIR$ NOVECTOR
#endif
do i = 1, points

    !-----------------------------------------------
    ! Calculate saturation vapour pressures
    !-----------------------------------------------
  esi(i) = qs (i) * p(i) * recip_epsilon
  esw(i) = qsl(i) * p(i) * recip_epsilon

    !-----------------------------------------------
    ! Calculate density of air
    !-----------------------------------------------
  if (l_mr_physics) then

    ! rho is the dry density
    rho(i) = rhodz_dry(i) / deltaz(i)

  else

    ! rho is the moist density
    rho(i) = rhodz_moist(i) / deltaz(i)

  end if  ! l_mr_physics

      ! Calculate the inverse of the air density
  rhor(i) = one / rho(i)

      !-----------------------------------------------
      ! Estimate latent heat correction to rate of evaporation
      !-----------------------------------------------
  lheat_correc_liq(i) = one/(one+repsilon*lc**2*qsl(i)                         &
                           /(cp*r*t(i)**2))
  lheat_correc_ice(i) = one/(one+repsilon*(lc+lf)**2*qs(i)                     &
                           /(cp*r*t(i)**2))

      !-----------------------------------------------
      ! Calculate CFL timestep divided by level separation
      !-----------------------------------------------

      ! Use the formulation based on the heights of the levels
  dhi(i) = timestep_mp/deltaz(i)

      ! Calculate the inverse
  dhir(i)       = one/dhi(i)

end do ! points

    !-----------------------------------------------
    ! Correction factors due to air density and temperature
    !-----------------------------------------------

do i = 1, points

      ! Correction of fall speeds
  if (l_mr_physics) then

    corr_in(i) = rho1*deltaz(i) / rhodz_moist(i)

  else

    corr_in(i) = rho1*rhor(i)

  end if

      ! Correction factor in viscosity etc. due to temperature

  corr2_in(i) = (t(i) * (one / zerodegc))

end do

call powr_v( points, corr_in,  0.4_real_lsprec,   corr  )
call powr_v( points, corr2_in, cpwr , corr2 )

! Determine which set of fallspeed parameters
! to use in process rate calculations, based on
! mass-weighted mean fallspeed

if ( l_diff_icevt .and. l_psd ) then
  do i = 1, points

    ! -------------------------------------------------------------
    ! Check that ice cloud fraction is sensible. This should mimic
    ! exactly what is done in lsp_ice prior to the process
    ! rate calculations.
    ! -------------------------------------------------------------

    cfice(i)  = max( cficekeep(i), 0.001_real_lsprec )
    cficei(i) = one / cfice(i)

  end do

  !-------------------------------------------------------------
  ! Calculate the psd moments required for crystal and
  ! aggregate fallspeeds
  !-------------------------------------------------------------


  call lsp_moments(points,rho,t,qcf,cficei,cx(182),m_bic_dic)
    ! ice mass flux moment with crystal parameters

  call lsp_moments(points,rho,t,qcf,cficei,cx(82),m_bi_di)
    ! ice mass flux moment with aggregate parameters

  do i = 1, points

     !-------------------------------------------------------------
     ! Calculate mass-weighted mean ice fallspeeds for the
     ! two vt-D relations
     !-------------------------------------------------------------
    if (qcf(i) >  m0) then
        ! Only calculate a fallspeed if qcf exceeds minimum value
        ! set by the nucleation mass. This is how the
        ! fallspeeds are calculated in lsp_fall
      vm_cry(i) = constp(182)* corr(i) * m_bic_dic(i)                          &
                / (rho(i) * qcf(i) * cficei(i))
      vm_agg(i) = constp(82) * corr(i) * m_bi_di(i)                            &
                / (rho(i) * qcf(i) * cficei(i))
      !-------------------------------------------------------------
      ! Determine which set of fallspeed parameters to use
      ! in process rate calculations
      !-------------------------------------------------------------
      if ( vm_cry(i)  <=  vm_agg(i) ) then
        l_use_agg_vt(i) = .false.
                   ! Use crystal vt-D parameters
      else
        l_use_agg_vt(i) = .true.
                   ! Use aggregate vt-D parameters
      end if
    end if
    if (qcf(i)  <=   m0) then
        ! Use crystal fallspeeds for very small qcf values
      l_use_agg_vt(i) = .false.
    end if

  end do

end if !  l_diff_icevt .and. l_psd

do i = 1, points

  corr2(i) = corr2(i)*( tcor1 /(t(i)+ tcor2))

  ! Combined correction factor

  if (l_mr_physics) then

    rocor_in(i) = rhodz_moist(i) / deltaz(i) * corr(i) * corr2(i)

  else

    rocor_in(i) = rho(i)*corr(i)*corr2(i)

  end if

end do

call powr_v( points, rocor_in, half, rocor )

    !-----------------------------------------------
    ! Calculate ice particle size distributions
    !-----------------------------------------------

do i = 1, points

  ! Calculate a temperature factor for N0aggregates
  tcg_in(i)  = -cx(32)*max(t(i)-zerodegc,t_agg_min)

  ! Define inverse of TCG values
  tcgi_in(i)  = -tcg_in(i)

end do

call exp_v( points, tcg_in,   tcg   )
call exp_v( points, tcgi_in,  tcgi  )

if ( l_crystals ) then

  do i = 1, points
    ! Calculate a temperature factor for N0crystals
    tcgc_in(i) = -cx(12)*max(t(i)-zerodegc,t_agg_min)

    ! Define inverse of TCGC values
    tcgci_in(i) = -tcgc_in(i)

  end do

  call exp_v( points, tcgc_in,  tcgc  )
  call exp_v( points, tcgci_in, tcgci )

end if ! l_crystals

    !-----------------------------------------------
    ! Calculate graupel size distributions
    !-----------------------------------------------
do i = 1, points

  tcgg(i)  = one
  tcggi(i) = one ! (Equivalent to 1.0/tcgg, so will always be 1.0)

end do

! If using prognostic precipitation fraction...
if ( l_mcr_precfrac ) then

  ! Initialise the rain fraction equal to the prognostic precip fraction
  ! at level k (the first thing rainfrac is used for is computing the
  ! fall-speed of the pre-existing rain at level k, so need to set to
  ! fraction consistent with existing precip at level k)
  do i = 1, points
    rainfrac(i) = precfrac_k(i)
  end do

end if  ! ( l_mcr_precfrac )

! Initialise the not falling-out fractions of ice, rain and graupel to 1
do i = 1, points
  agg_nofall(i) = one
  cry_nofall(i) = one
end do
if ( l_mcr_qrain ) then
  do i = 1, points
    rain_nofall(i) = one
  end do
end if
if ( l_mcr_qgraup ) then
  do i = 1, points
    graup_nofall(i) = one
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_init
end module lsp_init_mod
