! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Nucleation of ice particles
module lsp_nucleation_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_NUCLEATION_MOD'

contains

subroutine lsp_nucleation(                                                     &
  points, timestep,                                                            &
                                          ! Number of points and tstep
  q, qcl, tnuc_new, qrain, qcf, qgraup, t,                                     &
                                          ! Water contents, temperature
  qs, qsl,                                                                     &
                                          ! Saturated quantities
  cfliq,                                                                       &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
  area_liq, area_mix, cf, cfl, cff, rainfrac,                                  &
                                          ! Current cloud and rain
                                          ! fractions for updating
  rain_liq, rain_mix, rain_ice, rain_clear, rainfraci,                         &
                                          ! Overlaps of cloud with rain
  rho, rhor, lheat_correc_ice,                                                 &
                                          ! Parametrization information
  corr,  dhir, rain_nofall,                                                    &
                                          ! Parametrization information
  lfrcp, lsrcp,                                                                &
                                          ! Microphysical information
  hettransfer, hettransfer2, homtransfer, homtransfer2,                        &
                                          ! Mass transfer diagnostics
  one_over_tsi,                                                                &
                                          ! 1/(timestep*iterations)
  cftransfer, cfltransfer, cfftransfer, rf_transfer_diag,                      &
                                          ! Cloud transfer diagnostics
  dqprec_liq, dqprec_mix, dqprec_ice, dqprec_clear,                            &
                                          ! Total precip incs within partitions
  wtrac_mp_cpr_old                                                             &
                                          ! Store used for water tracers
  )

use lsprec_mod, only: thomo, m0, tnuc, zerodegc,                               &
                      zero, one
use mphys_inputs_mod,     only: l_het_freezing_rain, l_mcr_precfrac,           &
                                l_subgrid_graupel_frac, i_update_precfrac,     &
                                i_homog_areas, l_progn_tnuc

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,             only: real_lsprec

! Dr Hook modules
use yomhook,         only: lhook, dr_hook
use parkind1,        only: jprb, jpim

use lsp_het_freezing_rain_mod,   only: lsp_het_freezing_rain
use free_tracers_inputs_mod,     only: l_wtrac
use wtrac_mphys_mod,             only: mp_cpr_old_wtrac_type

implicit none

! Purpose:
!   Update cloud prognostics as a result of homogeneous and
!   heterogeneous ice nucleation

! Method:
!   Homogeneous nucleation converts all liquid to ice with a temperature
!   threshold. Heterogeneous nucleation converts a seed amount of liquid
!   or vapour to ice. Immersion freezing of raindrops is based on Bigg (1953).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

! Subroutine Arguments

integer, intent(in) ::                                                         &
  points
                        ! Number of points to calculate
real (kind=real_lsprec), intent(in) ::                                         &
  timestep,                                                                    &
                        ! Timestep / s
  qs(points),                                                                  &
                        ! Saturated humidity wrt ice / kg kg-1
  qsl(points),                                                                 &
                        ! Saturated humidity wrt liquid / kg kg-1
  cfliq(points),                                                               &
                        ! Fraction of gridbox with liquid cloud
  area_liq(points),                                                            &
                        ! Fraction of gridbox with only liquid cloud
  area_mix(points),                                                            &
                        ! Fraction of gridbox with mixed phase cloud
  rho(points),                                                                 &
                          ! Air density / kg m-3
  rhor(points),                                                                &
                        ! 1/Air density / kg m-3
  lheat_correc_ice(points),                                                    &
                               ! Ice latent heat correction factor
  corr(points),                                                                &
                        ! Air density fall speed correction (no units)
  dhir(points),                                                                &
                        ! Thickness of model layer / timestep / m s-1
  rain_nofall(points),                                                         &
                        ! Fraction of the rain-mass that is not falling out
  lfrcp,                                                                       &
                          ! Latent heat of fusion
                          ! /heat capacity of air (cP) / K
  lsrcp,                                                                       &
                          ! Latent heat of sublimation/cP / K
  one_over_tsi          ! 1/(timestep*iterations)

real (kind=real_lsprec), intent(in out) ::                                     &
  q(points),                                                                   &
                        ! Vapour content / kg kg-1
  qcl(points),                                                                 &
                        ! Liquid water content / kg kg-1
  qrain(points),                                                               &
                        ! Rain mixing ratio / kg kg-1
  qcf(points),                                                                 &
                        ! Ice water content in ice category to be
!                           updated    / kg kg-1
  qgraup(points),                                                              &
                        ! Graupel mixing ration / kg kg-1
  t(points),                                                                   &
                          ! Temperature / K
  cf(points),                                                                  &
                          ! Current cloud fraction
  cfl(points),                                                                 &
                          ! Current liquid cloud fraction
  cff(points),                                                                 &
                          ! Current ice cloud fraction
  rainfrac(points),                                                            &
                          ! Current rain fraction
  rain_liq(points),                                                            &
                          ! Frac. of gridbox with rain and liquid cld.
  rain_mix(points),                                                            &
                          ! Frac. of grbx. with rain and mixed ph. cld
  rain_ice(points),                                                            &
                          ! Frac. of gridbox with rain and ice
  rain_clear(points),                                                          &
                          ! Frac. of gridbox with rain and no cloud
  rainfraci(points),                                                           &
                          ! 1.0 / rainfrac
  homtransfer(points),                                                         &
                          ! Mass homog. nucleated this ts / kg kg-1
  homtransfer2(points),                                                        &
                          ! Mass homog. (rain)
                          !  nucleated this ts / kg kg-1
  hettransfer(points),                                                         &
                          ! Mass het. nucleated this ts / kg kg-1
  hettransfer2(points)
                          ! Mass het. (rain)
                          !  nucleated this ts / kg kg-1
real (kind=real_lsprec), intent(in) ::                                         &
    tnuc_new(points)
                          ! tnuc as function of dust / deg cel

real (kind=real_lsprec), intent(in out) ::                                     &
  cftransfer(points),                                                          &
                           ! Cloud fraction increment this tstep
  cfltransfer(points),                                                         &
                           ! Liquid cloud fraction inc this tstep
  cfftransfer(points),                                                         &
                           ! Ice cloud fraction inc this tstep
  rf_transfer_diag(points)
                           ! Rain fraction inc this tstep

! Precipitation mass increment within each sub-region of the precipitation
! fraction (bits that overlap liquid-cloud, ice-cloud, clear-sky).
! Used to update the prognostic precipitation fraction.
real (kind=real_lsprec), intent(in out) :: dqprec_liq(points)
real (kind=real_lsprec), intent(in out) :: dqprec_mix(points)
real (kind=real_lsprec), intent(in out) :: dqprec_ice(points)
real (kind=real_lsprec), intent(in out) :: dqprec_clear(points)

type(mp_cpr_old_wtrac_type), intent(in out) :: wtrac_mp_cpr_old

! Local Variables

integer ::                                                                     &
  i                 ! Loop counter for points

real (kind=real_lsprec) ::                                                     &
  dqi,                                                                         &
                        ! Temporary in calculating ice transfers
  dqil,                                                                        &
                        ! Temporary in calculating liquid transfers
  rhnuc
                        ! Nucleation relative humidity threshold
real (kind=real_lsprec) :: tnuc_kelvin
                        ! Nucleation temperature in Kelvin

! Total increment to qrain (+qgraup) from all ice nucleation processes
real (kind=real_lsprec) :: dqprec(points)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_NUCLEATION'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if ( l_mcr_precfrac .and. i_update_precfrac == i_homog_areas ) then
  ! If using prognostic precip fraction, save value of qrain before
  ! nucleation so that we can calculate the increment afterwards
  if ( l_subgrid_graupel_frac ) then
    ! If using sub-grid graupel fraction, want total inc to rain + graupel
    do i = 1, points
      dqprec(i) = qrain(i) + qgraup(i)
    end do
  else
    ! Otherwise only want rain increment
    do i = 1, points
      dqprec(i) = qrain(i)
    end do
  end if
end if

do i = 1, points

      !-----------------------------------------------
      ! Homogeneous nucleation - freeze all liquid if T < threshold
      !-----------------------------------------------
  if (t(i) <  (zerodegc+thomo) .and. qcl(i) >  zero) then

        !-----------------------------------------------
        ! Update cloud fractions
        !-----------------------------------------------
    cfftransfer(i) = cfftransfer(i) + (cf(i) - cff(i) )* one_over_tsi
    cfltransfer(i) = cfltransfer(i) -  cfl(i)* one_over_tsi

    cff(i) = cf(i)
    cfl(i) = zero

        !-----------------------------------------------
        ! Update water contents
        !-----------------------------------------------
    homtransfer(i) = homtransfer(i) + qcl(i)*one_over_tsi
    if (l_wtrac) wtrac_mp_cpr_old%qchange(i) = qcl(i)

    qcf(i) = qcf(i) + qcl(i)
    t(i)   = t(i)   + lfrcp * qcl(i)
    qcl(i) = zero

  end if ! T lt 0+thomo etc.

      !-----------------------------------------------
      ! Homogeneous nucleation - freeze all rain if T < threshold
      !-----------------------------------------------
  if (t(i) <  (zerodegc+thomo) .and. qrain(i) >  zero) then

        !-----------------------------------------------
        ! Update rain fractions
        !-----------------------------------------------
    cfftransfer(i)      = cfftransfer(i) +                                     &
                        ( cf(i) - cff(i) ) * one_over_tsi
    cff(i)      = cf(i)

    if ( .not. l_mcr_precfrac ) then
      ! Update rain-fraction
      ! (not doing this here if using prognostic precip fraction,
      !  as in that case rainfrac is updated later)
      rf_transfer_diag(i) = rf_transfer_diag(i) -                              &
                            rainfrac(i) * one_over_tsi
      rainfrac(i) = zero
    end if

        !-----------------------------------------------
        ! Update water contents
        !-----------------------------------------------
    homtransfer2(i) = homtransfer2(i) + qrain(i)*one_over_tsi

    if (l_wtrac) wtrac_mp_cpr_old%hom_qr(i) = qrain(i)      ! qrain -> qcf

    qcf(i)   = qcf(i) + qrain(i)
    t(i)     = t(i)   + lfrcp * qrain(i)
    qrain(i) = zero

  end if ! T lt 0+thomo etc.

      !-----------------------------------------------
      ! Heterogeneous nucleation for T < tnuc
      !-----------------------------------------------
  if (l_progn_tnuc) then
      ! Use prognostic-dust-based nucleation temperature
    tnuc_kelvin = zerodegc+tnuc_new(i)
  else
      ! Use default nucleation temperature
    tnuc_kelvin = zerodegc+tnuc
  end if

  if (t(i) <  tnuc_kelvin .and. area_liq(i) >  zero                            &
      .and. cfliq(i)  >   zero ) then

        !-----------------------------------------------
        ! Calculate number and mixing ratio of nucleated crystals
        !-----------------------------------------------
    dqi=min(0.01_real_lsprec*exp(-0.6_real_lsprec*(t(i)-zerodegc)),            &
               1.0e5_real_lsprec)
    dqi=m0 * dqi * rhor(i)

        !-----------------------------------------------
        ! How much moisture is available for ice formation
        !-----------------------------------------------
        ! Firstly, calculate the threshold relative humidity
    rhnuc=(188.92_real_lsprec+2.81_real_lsprec*(t(i)-zerodegc)                 &
                +0.013336_real_lsprec*(t(i)-zerodegc)**2)*0.01_real_lsprec
    rhnuc=min(rhnuc,one)-0.1_real_lsprec
    rhnuc=max(qsl(i)*rhnuc,qs(i))

        ! Next calculate the available water
    dqil = (qcl(i)/cfliq(i)+qsl(i)-rhnuc)
    dqi = max(area_liq(i) * min(dqi,dqil*lheat_correc_ice(i)),zero)
    qcf(i) = qcf(i)+dqi

        !-----------------------------------------------
        ! Store nucleation rate
        !-----------------------------------------------
    hettransfer(i) = hettransfer(i)+dqi*one_over_tsi

        !-----------------------------------------------
        ! Calculate mass transfers
        !-----------------------------------------------

          ! Firstly take mass from liquid water

    if (area_mix(i) > zero) then

      dqil = min(dqi,qcl(i)*area_liq(i)/cfliq(i))

    else  ! area_liq/cfliq must be 1

      dqil = min(dqi,qcl(i))

    end if

    qcl(i)  = qcl(i)-dqil
    t(i)    = t(i)+lfrcp*dqil
          ! If more ice is needed then the mass comes from vapour
    dqi  = dqi-dqil
    t(i) = t(i)+lsrcp*dqi
    q(i) = q(i)-dqi

    ! Store phase changes for water tracer use
    if (l_wtrac) then
      ! qcl -> qcf
      wtrac_mp_cpr_old%qchange(i) = wtrac_mp_cpr_old%qchange(i) + dqil
      ! q -> qcf
      wtrac_mp_cpr_old%het_q(i)   = dqi
    end if

        !-----------------------------------------------
        ! Udate cloud fractions.
        !-----------------------------------------------

    cfftransfer(i) = cfftransfer(i) + (cf(i) - cff(i))*one_over_tsi
    cff(i)      = cf(i)


  end if  ! On temperature threshold

end do  ! Points

if (l_het_freezing_rain) then
      ! Call heterogeneous freezing rain
  call lsp_het_freezing_rain(points, timestep,                                 &
                qrain, qcf, qgraup, t,                                         &
                cf, cff, rainfrac,                                             &
                rain_liq, rain_mix, rain_ice,                                  &
                rho, rhor, corr, dhir, rain_nofall,                            &
                lfrcp, hettransfer2, one_over_tsi,                             &
                cftransfer, cfftransfer, rf_transfer_diag                      &
               )
end if  ! l_het_freezing_rain

if ( l_mcr_precfrac .and. i_update_precfrac == i_homog_areas ) then
  ! If using prognostic precip fraction...

  ! Compute increment to precip from all nucleation processes,
  ! by differencing latest value with saved value from start of this routine
  if ( l_subgrid_graupel_frac ) then
    ! If using sub-grid graupel fraction, want total inc to rain + graupel
    do i = 1, points
      dqprec(i) = qrain(i) + qgraup(i) - dqprec(i)
    end do
  else
    ! Otherwise only want rain increment
    do i = 1, points
      dqprec(i) = qrain(i) - dqprec(i)
    end do
  end if

  ! Increment applies in all of the rain area partitions;
  ! divvy it up between them in proportion to area
  do i = 1, points
    ! rainfraci stores 1 / (rain_liq + rain_mix + rain_ice + rain_clear)
    dqprec_liq(i) = dqprec(i) * rain_liq(i) * rainfraci(i)
    dqprec_mix(i) = dqprec(i) * rain_mix(i) * rainfraci(i)
    dqprec_ice(i) = dqprec(i) * rain_ice(i) * rainfraci(i)
    dqprec_clear(i) = dqprec(i) * rain_clear(i) * rainfraci(i)
  end do

end if  ! ( l_mcr_precfrac .and. i_update_precfrac == i_homog_areas )

! Note: for the other option for updating precfrac, there is no change to
! precfrac from this process, since it affects the whole rain / graupel
! region homogeneously.

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_nucleation
end module lsp_nucleation_mod
