! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Heterogeneous nucleation of rain
module lsp_het_freezing_rain_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_HET_FREEZING_RAIN_MOD'

contains

subroutine lsp_het_freezing_rain(                                              &
  points, timestep,                                                            &
                                          ! Number of points and tstep
  qrain, qcf, qgraup, t,                                                       &
                                          ! Water contents, temperature
  cf, cff, rainfrac,                                                           &
                                          ! Current cloud and rain
                                          ! fractions for updating
  rain_liq, rain_mix, rain_ice,                                                &
                                          ! Overlaps of cloud with rain
  rho, rhor,                                                                   &
                                          ! Parametrization information
  corr,  dhir, rain_nofall,                                                    &
                                          ! Parametrization information
  lfrcp,                                                                       &
                                          ! Microphysical information
  hettransfer2,                                                                &
                                          ! Mass transfer diagnostic
  one_over_tsi,                                                                &
                                          ! 1/(timestep*iterations)
  cftransfer, cfftransfer, rf_transfer_diag                                    &
                                          ! Cloud transfer diagnostics
  )

use lsprec_mod,           only: zerodegc, zero,                                &
                                one, qcfmin, cx, constp
use mphys_inputs_mod,     only: l_mcr_qrain,                                   &
                                l_mcr_qgraup, l_mcr_precfrac

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,             only: real_lsprec

! Dr Hook modules
use yomhook,         only: lhook, dr_hook
use parkind1,        only: jprb, jpim

implicit none

! Purpose:
!   Update cloud prognostics as a result of heterogeneous rain freezing

! Method:
!   Immersion freezing of raindrops is based on Bigg (1953).
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
  rho(points),                                                                 &
                        ! Air density / kg m-3
  rhor(points),                                                                &
                        ! 1/Air density / kg m-3
  corr(points),                                                                &
                        ! Air density fall speed correction (no units)
  dhir(points),                                                                &
                        ! Thickness of model layer / timestep / m s-1
  rain_nofall(points),                                                         &
                        ! Fraction of the rain-mass that is not falling out
  lfrcp,                                                                       &
                        ! Latent heat of fusion
                        ! /heat capacity of air (cP) / K
  one_over_tsi          ! 1/(timestep*iterations)

real (kind=real_lsprec), intent(in out) ::                                     &
  qrain(points),                                                               &
                        ! Rain mixing ratio / kg kg-1
  qcf(points),                                                                 &
                        ! Ice water content in ice category to be
                        ! updated    / kg kg-1
  qgraup(points),                                                              &
                        ! Graupel mixing ration / kg kg-1
  t(points),                                                                   &
                        ! Temperature / K
  cf(points),                                                                  &
                        ! Current cloud fraction
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
  hettransfer2(points)
                        ! Mass het. (rain)
                        !  nucleated this ts / kg kg-1

real (kind=real_lsprec), intent(in out) ::                                     &
  cftransfer(points),                                                          &
                           ! Cloud fraction increment this tstep
  cfftransfer(points),                                                         &
                        ! Ice cloud fraction inc this tstep
  rf_transfer_diag(points)
                        ! Rain fraction inc this tstep

! Local Variables

integer ::                                                                     &
  i                     ! Loop counter for points

real (kind=real_lsprec) ::                                                     &
  dqir(points),                                                                &
                        ! Temporary in calculating rain transfers
  cff_final(points),                                                           &
                        ! Final ice cloud fraction
  delta_cff(points),                                                           &
                        ! Ice cloud fraction transferred
  lamr1(points),                                                               &
                        ! Reciprocal of slope parameter in raindrop
                        ! size distribution  / m
  lamr3
                        ! lamr1**3  /m3

real (kind=real_lsprec), parameter :: a_bigg = 0.66_real_lsprec
                        ! Bigg parameter A  / K-1
real (kind=real_lsprec), parameter :: b_bigg = 100.0_real_lsprec
                        ! Bigg parameter B  / m3 s-1

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_HET_FREEZING_RAIN'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
do i = 1, points

  !-----------------------------------------------
  ! Immersion freezing of raindrops - allowed if T < -4 C
  !-----------------------------------------------

  if (t(i) <  (zerodegc - 4.0_real_lsprec) .and.                               &
      qrain(i)*rain_nofall(i) >  qcfmin) then

    !-----------------------------------------------
    ! Calculate reciprocal of lambda for rain
    !-----------------------------------------------

    if (l_mcr_qrain) then
      ! rain is a mixing ratio (kg kg-1)
      lamr1(i) = (rho(i) * constp(50) * qrain(i)*rain_nofall(i) / rainfrac(i))
    else
      ! rain is a flux (kg m-2 s-1)
      lamr1(i) = (qrain(i) * rho(i) * dhir(i)                                  &
               /( rainfrac(i)*constp(42)*corr(i)) )**(cx(42))
    end if  ! l_mcr_qrain

    !------------------------------------------------
    ! Calculate transfer
    !------------------------------------------------

    ! Calculate lambda^-3
    lamr3 = lamr1(i)*lamr1(i)*lamr1(i)

    ! constp(58)=20.0*pi*pi*rho_water*x1r
    dqir(i) = constp(58) * b_bigg * lamr1(i)**(-cx(46)) * rhor(i) *            &
               (exp(a_bigg*(zerodegc-t(i)))-one) *                             &
               lamr3 * lamr3 * lamr1(i) * rainfrac(i) * timestep

    ! Limit transfer to the mass of rain that is available
    dqir(i) = min(dqir(i), qrain(i))

    ! Store raindrop immersion freezing rate
    hettransfer2(i) = hettransfer2(i) + dqir(i) * one_over_tsi

    !------------------------------------------------
    ! Adjust graupel (or snow/ice) and rain contents
    !------------------------------------------------

    if (l_mcr_qgraup) then
      ! Move contents to graupel
      qgraup(i) = qgraup(i) + dqir(i)
    else
      ! Put transfer quantities as snow/ice
      qcf(i) = qcf(i) + dqir(i)
    end if  ! l_mcr_qgraup

    ! Now remove the appropriate transfer amount
    ! from the rain category and adjust the
    ! temperature

    qrain(i) = qrain(i) - dqir(i)
    t(i)     = t(i)     + dqir(i) * lfrcp

    !------------------------------------------------
    ! Update cloud fractions
    !------------------------------------------------

    !Update ice cloud fractions
    if (.not. l_mcr_qgraup) then

      cff_final(i) = max(rainfrac(i),cff(i))
      delta_cff(i) = cff_final(i) - cff(i)

      !Calculate transfer rate diagnostics
      cfftransfer(i) = cfftransfer(i) + delta_cff(i) * one_over_tsi
      cftransfer(i)  = cftransfer(i)  + delta_cff(i) * one_over_tsi

      cf(i)  = cf(i)  + delta_cff(i)
      cff(i) = cff(i) + delta_cff(i)

    end if  ! l_mcr_graupel

    !-----------------------------------------------
    ! Update rain fractions
    !-----------------------------------------------
    !   Assume no rain fraction change except if all
    !   rain mass is frozen
    !
    if (qrain(i) <= zero .and. (.not. l_mcr_precfrac) ) then
      ! (not doing this here if using prognostic precip fraction,
      !  as in that case rainfrac is updated later)

      rf_transfer_diag(i) = rf_transfer_diag(i) -                              &
                          rainfrac(i) * one_over_tsi

      rainfrac(i)  = zero
      rain_liq(i)  = zero
      rain_mix(i)  = zero
      rain_ice(i)  = zero

    end if  !qrain = 0

  end if  ! tc < -4, qrain(i) >  qcfmin

end do  ! Points

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_het_freezing_rain
end module lsp_het_freezing_rain_mod
