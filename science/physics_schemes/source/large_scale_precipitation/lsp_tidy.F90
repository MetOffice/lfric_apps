! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Riming of ice particles
! Subroutine Interface:
module lsp_tidy_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_TIDY_MOD'

contains

subroutine lsp_tidy(                                                           &
  points,                                                                      &
                                          ! Number of points
  one_over_tsi,                                                                &
                                          ! 1/(timestep*iterations)
  q, qcl, qcf, qcf2, qrain, qgraup, t,                                         &
                                          ! Water contents and temp
  area_liq, area_mix, area_ice,                                                &
                                          ! Cloud fraction information
  cfice, cficei,                                                               &
                                          ! at start of microphysics ts
  cf, cfl, cff,                                                                &
                                          ! Current cloud fractions for
                                          ! updating
  rainfrac, rain_liq, rain_mix,                                                &
                                          ! Rain fractions
  rain_ice, rain_clear,                                                        &
  q_ice, qs, qsl, snow_agg, snow_cry, rainrate, grauprate,                     &
                                          ! Other water contents
  rho, rhor, p,                                                                &
                                          ! Other model prognostics
  cttemp,dhi,dhir,frac_ice_fall,                                               &
                                          ! Other information
  lcrcp, lfrcp, lsrcp,                                                         &
                                          ! Microphysical information
  psdep, pidep, psmlt, pimlt, prevp,                                           &
                                          ! Mass transfer diagnostics
  cftransfer,cfltransfer,cfftransfer,                                          &
                                          ! Cloud transfer diagnostics
  rftransfer,                                                                  &
                                          ! Rain transfer diagnostics
  precfrac_k,                                                                  &
                                          ! Prog fraction of rain/graup on level
  precfrac_fall,                                                               &
                                          ! Prog fraction of falling rain/graup
  wtrac_mp_cpr, wtrac_mp_cpr_old                                               &
                                          ! Water tracers
  )

!Use in reals in lsprec precision, both microphysics related and general atmos
use lsprec_mod,       only: qcfmin, zerodegc, tw1, tw2, tw3, tw4, tw5,         &
                            zero, one

! Microphysics modules- logicals and integers
use mphys_inputs_mod, only: l_mcr_qcf2, l_mcr_precfrac, l_subgrid_graupel_frac,&
                            l_proc_fluxes

use lsp_combine_precfrac_mod, only: lsp_combine_precfrac

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,         only: real_lsprec

use science_fixes_mod, only: l_fix_tidy_rainfracs

! Water tracers
use free_tracers_inputs_mod, only: l_wtrac, n_wtrac
use wtrac_mphys_mod,         only: mp_cpr_wtrac_type, mp_cpr_old_wtrac_type
use wtrac_calc_ratio_mod,    only: wtrac_calc_ratio_fn

! Dr Hook Modules
use yomhook,          only: lhook, dr_hook
use parkind1,         only: jprb, jpim

implicit none

! Purpose:
!   Tidy up small numerical values after large-scale precipitation.
!   Ideally, this code would be unnecessary, but in reality there
!   are always going to be small values left over numerically that
!   should be reset.

! Method:
!   1. Evaporate rain amounts
!   2. Evaporate small ice amounts
!   3. Reset cloud fractions
!   4. Melt small snow amounts
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
                        ! Number of points

real (kind=real_lsprec), intent(in) ::                                         &
  one_over_tsi,                                                                &
                        ! 1/(timestep*iterations)
  qcl(points),                                                                 &
                        ! Liquid water content / kg kg-1
  area_liq(points),                                                            &
                        ! Fraction of gridbox with liquid-only cloud
  area_mix(points),                                                            &
                        ! Fraction of gridbox with mixed phase cloud
  area_ice(points),                                                            &
                        ! Fraction of gridbox with ice-only cloud
  cfice(points),                                                               &
                        ! Fraction of gridbox with ice cloud
  cficei(points),                                                              &
                        ! 1/Fraction of gridbox with ice cloud
  q_ice(points),                                                               &
                        ! Vapour content in ice cloud / kg kg-1
  qs(points),                                                                  &
                        ! Saturated water content wrt ice / kg kg-1
  qsl(points),                                                                 &
                        ! Saturated water content wrt liquid / kg kg-1
  rho(points),                                                                 &
                        ! Air density / kg m-3
  rhor(points),                                                                &
                        ! 1 / air density / m3 kg-1
  p(points),                                                                   &
                        ! Air pressure / N m-2
  dhi(points),                                                                 &
                        ! Timestep/layer thickness / s m-1
  dhir(points)
                        ! 1/dhi / m s-1
real (kind=real_lsprec), intent(in out) ::                                     &
  frac_ice_fall(points)
                        ! Ice cloud fraction associated with the ice fall-flux
real (kind=real_lsprec), intent(in) ::                                         &
  lcrcp,                                                                       &
                        ! Latent heat of condensation/cP / K
  lfrcp,                                                                       &
                        ! Latent heat of fusion/cP / K
  lsrcp             ! Latent heat of sublimation/cP / K

real (kind=real_lsprec), intent(in out) ::                                     &
  q(points),                                                                   &
                        ! Vapour mixing ratio / kg kg-1
  qcf(points),                                                                 &
                        ! Aggregates mixing ratio / kg kg-1
  qcf2(points),                                                                &
                        ! Crystals mixing ratio / kg kg-1
  qrain(points),                                                               &
                        ! Rain mixing ratio / kg kg-1
  qgraup(points),                                                              &
                        ! Graupel mixing ratio / kg kg-1
  t(points),                                                                   &
                        ! Temperature / K
  cttemp(points),                                                              &
                        ! Ice-cloud top temperature / K
  cf(points),                                                                  &
                        ! Current cloud fraction
  cfl(points),                                                                 &
                        ! Current liquid cloud fraction
  cff(points),                                                                 &
                        ! Current ice cloud fraction
  rainfrac(points),                                                            &
                        ! Rain fraction (no units)
  rain_liq(points),                                                            &
                        ! Overlap of rain with liquid (no units)
  rain_mix(points),                                                            &
                        ! Overlap of rain with mixed phase region
  rain_ice(points),                                                            &
                        ! Overlap of rain with ice
  rain_clear(points),                                                          &
                        ! Overlap of rain with clear sky
  snow_agg(points),                                                            &
                        ! Aggregate snowfall rate / kg m-2 s-1
  snow_cry(points),                                                            &
                        ! Crystal snowfall rate / kg m-2 s-1
  rainrate(points),                                                            &
                        ! Rainfall rate / kg m-2 s-1
  grauprate(points),                                                           &
                        ! Fall-out flux of graupel / kg m-2 s-1
  psdep(points),                                                               &
                        ! Deposition of aggregates diag. / kg kg s-1
  pidep(points),                                                               &
                        ! Deposition of crystals diag. / kg kg s-1
  psmlt(points),                                                               &
                        ! Melting of aggregates diagnostic / kg kg s-1
  pimlt(points),                                                               &
                        ! Melting of crystals diagnostic / kg kg s-1
  prevp(points),                                                               &
                        ! Evaporation of rain diagnostic / kg kg s-1
  cftransfer(points),                                                          &
                         ! Rate of change of bulk cloud frac / s-1
  cfltransfer(points),                                                         &
                         ! Rate of change of liquid cloud frac / s-1
  cfftransfer(points),                                                         &
                         ! Rate of change of ice cloud frac / s-1
  rftransfer(points),                                                          &
                         ! Rate of change of rain fraction / s-1
  precfrac_k(points),                                                          &
                         ! Prognostic precip fraction on current level
  precfrac_fall(points)
                         ! Value of prognostic precip fraction for fall-out

! Water tracer fields
type(mp_cpr_wtrac_type), intent(in out)     :: wtrac_mp_cpr(n_wtrac)
type(mp_cpr_old_wtrac_type), intent(in out) :: wtrac_mp_cpr_old

! Local Variables


real (kind=real_lsprec), parameter :: smallnum = 2.2e-14_real_lsprec
                        ! Small value used in if tests

integer ::                                                                     &
  i, i_wt               ! Loop counter for points

real (kind=real_lsprec) ::                                                     &
  dpr,                                                                         &
                        ! Mass transfer / kg kg-1
  dpr2,                                                                        &
                        ! Equivalent to dpr, but with units kg m-2 s-1
  temp7,                                                                       &
                        ! Wet bulb temperature / deg C
  tempw,                                                                       &
                        ! Saturation excess / kg kg-1
  cfnew,                                                                       &
                        ! New cloud fraction
  cflnew,                                                                      &
                        ! New liquid cloud fraction
  cffnew                ! New ice cloud fraction

! Rain-rate increment due to emergency melting; used for updating
! the prognostic precip fraction, if used.
real (kind=real_lsprec) :: drainrate_melt(points)

! Precipitation rate of rain + graupel
real (kind=real_lsprec) :: precrate(points)

! Water tracer to normal water ratio for snow
real (kind=real_lsprec) :: snow_ratio_wtrac

real (kind=real_lsprec) ::                                                     &
  dpr_wt,                                                                      &
                        ! Water mass transfer / kg kg-1
  dpr2_wt
                        ! Equivalent to dpr_wt, but with units kg m-2 s-1


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_TIDY'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i = 1, points

  !RJJG- should we use epsilon?
  if ((qrain(i)  <=  qcfmin .and. qcl(i) < 1.0e-16_real_lsprec)                &
       .or. qrain(i)  <   zero) then
        !-----------------------------------------------
        ! 1. If there is a tiny rain amount and no liquid then
        ! evaporate the rain
        !-----------------------------------------------
    dpr = qrain(i)

    if ( l_proc_fluxes ) then
      ! If handling fluxes consistently with values at current level,
      ! evaporate any remaining flux as well
      dpr = dpr + rainrate(i) * dhi(i)*rhor(i)
      rainrate(i) = zero
    end if

        ! Update prognostics
    q(i)   = q(i) + dpr
    t(i)   = t(i) - dpr * lcrcp
    qrain(i) = zero

        ! Update water tracers consistently
    if (l_wtrac) then
      do i_wt = 1, n_wtrac
        wtrac_mp_cpr(i_wt)%q(i) = wtrac_mp_cpr(i_wt)%q(i)                      &
                                  + wtrac_mp_cpr(i_wt)%qrain(i)
        wtrac_mp_cpr(i_wt)%qrain(i) = zero
      end do
    end if

        ! Update rain fractions
    if ( .not. l_mcr_precfrac ) then
      ! Don't bother with this update if using prognostic precip fraction
      rainfrac(i) = zero
      rain_liq(i) = zero
      rain_mix(i) = zero
      rain_ice(i) = zero
      rain_clear(i) = zero
      rftransfer(i) = rftransfer(i)                                            &
                    - rainfrac(i) * one_over_tsi
    end if

        ! Update evaporation rate
    prevp(i) = prevp(i) + dpr * one_over_tsi

  end if  ! rain lt qcfmin etc.

      !-----------------------------------------------
      ! 2. Evaporate small ice amounts
      !-----------------------------------------------
  if (l_mcr_qcf2) then
        ! Evaporate small aggregate and crystal amounts separately

        !-----------------------------------------------
        ! 2a. Aggregates
        !-----------------------------------------------
    if ((qcf(i)  <   qcfmin) .and.                                             &
      (t(i)  >   zerodegc .or.                                                 &
      (q_ice(i)  <=  qs(i) .and. area_mix(i)  <  smallnum)                     &
      .or. qcf(i) <  zero)) then
          ! Ice is very small and T>0 and ice is not growing by
          ! deposition, so evaporate it.
      dpr = qcf(i)

      if ( l_proc_fluxes ) then
        ! If handling fluxes consistently with values at current level,
        ! evaporate any remaining flux as well
        dpr = dpr + snow_agg(i) * dhi(i)*rhor(i)
        snow_agg(i) = zero
      end if

          ! Update prognostics

      q(i)   = q(i) + dpr
      t(i)   = t(i) - lsrcp * dpr
      qcf(i) = zero


          ! Update deposition rate
      psdep(i) = psdep(i) - dpr * one_over_tsi

    end if  ! qcf < qcfmin etc

        !-----------------------------------------------
        ! 2b. Crystals
        !-----------------------------------------------
    if ((qcf2(i)  <   qcfmin) .and.                                            &
      (t(i)  >   zerodegc .or.                                                 &
      (q_ice(i)  <=  qs(i) .and. area_mix(i)  <  smallnum)                     &
      .or. qcf2(i)  <   zero)) then
          ! Ice is very small and T>0 and ice is not growing by
          ! deposition, so evaporate it.
      dpr = qcf2(i)

      if ( l_proc_fluxes ) then
        ! If handling fluxes consistently with values at current level,
        ! evaporate any remaining flux as well
        dpr = dpr + snow_cry(i) * dhi(i)*rhor(i)
        snow_cry(i) = zero
      end if

          ! Update prognostics
      q(i)   = q(i) + dpr
      t(i)   = t(i) - lsrcp * dpr
      qcf2(i)= zero

          ! Update deposition rate
      pidep(i) = pidep(i) - dpr * one_over_tsi

    end if  ! qcf < 0 etc.

    if (qcf(i)  <=  zero .and. qcf2(i)  <=  zero) then

          !----------------------------------------------
          ! 2c. Update ice cloud amounts
          !----------------------------------------------

      cfftransfer(i) = -cff(i) * one_over_tsi
      cftransfer(i)  = ( cfl(i) - cff(i) ) * one_over_tsi
      cff(i)         = zero
      cf(i)          = cfl(i)


          !----------------------------------------------
          ! 2d. Update cloud top temperature
          !----------------------------------------------
      cttemp(i) = t(i)
    end if  ! qcf eq 0 etc

  else  ! l_mcr_qcf2

        !----------------------------------------------
        ! 2e. One prognostic ice content is active
        !----------------------------------------------
    if ((qcf(i)  <   qcfmin) .and.                                             &
      (t(i)  >   zerodegc .or.                                                 &
      (q_ice(i)  <=  qs(i) .and. area_mix(i)  <  smallnum)                     &
      .or. qcf(i) <= zero) ) then
          ! Ice is very small and T>0 and ice is not growing by
          ! deposition, so evaporate it.
      dpr = qcf(i)

      if ( l_proc_fluxes ) then
        ! If handling fluxes consistently with values at current level,
        ! evaporate any remaining flux as well
        dpr = dpr + snow_agg(i) * dhi(i)*rhor(i)
        snow_agg(i) = zero
      end if

          ! Update prognostics

      q(i)   = q(i) + dpr
      t(i)   = t(i) - lsrcp * dpr
      qcf(i) = zero

         ! Update water tracers consistently
      if (l_wtrac) then
        do i_wt = 1, n_wtrac
          wtrac_mp_cpr(i_wt)%q(i) = wtrac_mp_cpr(i_wt)%q(i)                    &
                                    + wtrac_mp_cpr(i_wt)%qcf(i)
          wtrac_mp_cpr(i_wt)%qcf(i) = zero
        end do
      end if

          ! Update deposition rate
      psdep(i) = psdep(i) - dpr * one_over_tsi

          ! Update ice cloud amounts
      cfftransfer(i) = -cff(i) * one_over_tsi
      cftransfer(i)  = (cfl(i)-cff(i)) * one_over_tsi

      cff(i) = zero
      cf(i)  = cfl(i)

          ! Update cloud top temperature
      cttemp(i) = t(i)

    end if  ! qcf lt qcfmin etc

  end if  ! l_mcr_qcf2

      !------------------------------------------------
      ! 3. Limit cloud fractions to physically reasonable values
      !------------------------------------------------
      ! It is not clear that this ought to be a parallel calculation
      ! but it is allowed to be for the present time.

        ! Calculate new cloud fraction values
  cffnew = max(min(cff(i),one),zero)
  cflnew = max(min(cfl(i),one),zero)
  cfnew  = max( max(cflnew,cffnew) , cf(i) )
  cfnew  = min( min(cflnew+cffnew,one), cfnew )

      ! Calculate transfer rates
  cfftransfer(i) = cfftransfer(i)                                              &
                 + (cffnew - cff(i)) * one_over_tsi
  cfltransfer(i) = cfltransfer(i)                                              &
                 + (cflnew - cfl(i)) * one_over_tsi
  cftransfer(i)  = cftransfer(i)                                               &
                 + (cfnew  - cf (i)) * one_over_tsi

      ! Update cloud fractions

  cff(i) = cffnew
  cfl(i) = cflnew
  cf(i)  = cfnew

end do  ! i = 1, points

    !------------------------------------------------
    ! 4. Emergency melting of snow to avoid excess snow at
    !    warm temperatures, i.e. > zerodegc
    !------------------------------------------------
if ( .not. l_proc_fluxes ) then
  ! Only need to do emergency melting of the falling-out snow fluxes
  ! if not applying process rates to fluxes in the main process calculations.
  ! (melting the fluxes again here would be double-counting)

  if ( l_mcr_precfrac ) then
    ! Initialise rain-rate increment due to emergency melting
    do i = 1, points
      drainrate_melt(i) = zero
    end do
  end if

  do i = 1, points

        !------------------------------------------------
        ! Melt snow_agg first
        !------------------------------------------------
    if (snow_agg(i)  >   zero .and. t(i)  >   zerodegc) then

          ! Numerical approximation of wet bulb temperature excess
      tempw = area_ice(i) * max(qsl(i)-q_ice(i),zero)*cficei(i)
      temp7 = t(i) - zerodegc -tempw*(tw1+tw2*(p(i)-tw3)-tw4*(t(i)-tw5))
      temp7 = max(temp7,zero)

          ! Calculate transfer rate
      dpr  = temp7 / lfrcp ! Rate based on Tw excess
      dpr2 = dpr*rho(i)*dhir(i)

          ! Limit to the amount of snow available
      dpr  = min(dpr , snow_agg(i) * dhi(i)*rhor(i) )
      dpr2 = min(dpr2, snow_agg(i))

          ! Add to melting rate
      psmlt(i) = psmlt(i) + dpr * one_over_tsi

      ! Update water tracers (before snow_agg is updated as need value
      !   before update in ratio calculation)
      if (l_wtrac) then
        do i_wt = 1, n_wtrac
          snow_ratio_wtrac = wtrac_calc_ratio_fn(i_wt,                         &
                             wtrac_mp_cpr(i_wt)%lssnow(i), snow_agg(i))
          ! Limit to the amount of water tracer in snow available
          dpr_wt  = min(dpr*snow_ratio_wtrac ,                                 &
                             wtrac_mp_cpr(i_wt)%lssnow(i) * dhi(i)*rhor(i) )
          dpr2_wt = min(dpr2*snow_ratio_wtrac, wtrac_mp_cpr(i_wt)%lssnow(i))

          wtrac_mp_cpr(i_wt)%lssnow(i) = wtrac_mp_cpr(i_wt)%lssnow(i)          &
                                         - dpr2_wt
          wtrac_mp_cpr(i_wt)%qrain(i)  = wtrac_mp_cpr(i_wt)%qrain(i)           &
                                         + dpr_wt
        end do
      end if

          ! Update values of snow and rain

      snow_agg(i) = snow_agg(i) - dpr2
      t(i)        = t(i)        - dpr * lfrcp

      if ( l_mcr_precfrac ) then
        ! If using prognostic precip fraction, store increment for
        ! updating the rain-rate and fall-out fraction later
        drainrate_melt(i) = drainrate_melt(i) + dpr2
      else
        ! Otherwise, just add the melted ice-flux onto the current level qrain
        qrain(i)    = qrain(i)    + dpr

        ! Update rain-fraction at level k
        if ( l_fix_tidy_rainfracs ) then

            ! Rain fraction will take on the value of the ice fraction
          rftransfer(i) = rftransfer(i) +                                      &
            max( cfice(i)-rainfrac(i) , zero)                                  &
            * one_over_tsi

            ! Update rain fractions
          rainfrac(i) = max(rainfrac(i),cfice(i))
          rain_liq(i) = min(area_liq(i),rainfrac(i))
          rain_mix(i) = min(area_mix(i),rainfrac(i)-rain_liq(i))
          rain_ice(i) =                                                        &
                 min(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i))
          rain_clear(i) =                                                      &
                 rainfrac(i)-rain_liq(i)-rain_mix(i)-rain_ice(i)

        end if  !  l_fix_tidy_rainfracs

      end if  ! not l_mcr_precfrac

    end if  !  snow_agg gt 0 etc

        !------------------------------------------------
        ! Melt snow_cry next
        !------------------------------------------------
    if (l_mcr_qcf2 .and.                                                       &
        snow_cry(i)  >   zero .and. t(i)  >   zerodegc) then

          ! Numerical approximation of wet bulb temperature excess
      tempw = area_ice(i) * max(qsl(i)-q_ice(i),zero)*cficei(i)
      temp7 = t(i) - zerodegc - tempw*(tw1+tw2*(p(i)-tw3)-tw4*(t(i)-tw5))
      temp7 = max(temp7,zero)

          ! Calculate transfer rate
      dpr  = temp7 / lfrcp ! Rate based on Tw excess
      dpr2 = dpr*rho(i)*dhir(i)

          ! Limit to the amount of snow available
      dpr  = min(dpr , snow_cry(i) * dhi(i)*rhor(i) )
      dpr2 = min(dpr2, snow_cry(i))

          ! Add to melting rate
      pimlt(i) = pimlt(i) + dpr * one_over_tsi



          ! Update values of snow and rain

      snow_cry(i) = snow_cry(i) - dpr2
      t(i)        = t(i)        - dpr * lfrcp

      if ( l_mcr_precfrac ) then
        ! If using prognostic precip fraction, store increment for
        ! updating the rain-rate and fall-out fraction later
        drainrate_melt(i) = drainrate_melt(i) + dpr2
      else
        ! Otherwise, just add the melted ice-flux onto the current level qrain
        qrain(i)    = qrain(i)    + dpr

            ! Update rain fractions

          ! Rain fraction will take on the value of the ice fraction
        ! Note: not sure why this is done for melting of crystals but
        ! not melting of aggregates??
        rftransfer(i) = rftransfer(i) +                                        &
          max( cfice(i)-rainfrac(i), zero)                                     &
          * one_over_tsi

        rainfrac(i) = max(rainfrac(i),cfice(i))
        rain_liq(i) = min(area_liq(i),rainfrac(i))
        rain_mix(i) = min(area_mix(i),rainfrac(i)-rain_liq(i))
        rain_ice(i) =                                                          &
               min(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i))
        rain_clear(i) =                                                        &
               rainfrac(i)-rain_liq(i)-rain_mix(i)-rain_ice(i)
      end if

    end if  ! snow_cry lt 0 etc

  end do  ! Points

  if ( l_mcr_precfrac ) then
    ! If using prognostic precip fraction...

    ! Update the fraction of precip falling out of the current level
    ! with the increment due to emergency melting of falling out ice...

    if ( l_subgrid_graupel_frac ) then
      ! If including graupel in the prognostic precip fraction,
      ! the precip-rate includes the graupel flux:
      do i = 1, points
        precrate(i) = rainrate(i) + grauprate(i)
      end do
    else
      ! Otherwise it is just the rainrate
      do i = 1, points
        precrate(i) = rainrate(i)
      end do
    end if

    ! Calculate combined fraction of already falling out precip
    ! and rain-flux just created by melting the ice flux
    ! (precfrac_fall is updated)
    call lsp_combine_precfrac( points,                                         &
                               precrate, drainrate_melt,                       &
                               precfrac_fall, frac_ice_fall )

    do i = 1, points
      ! Update the rain-rate with the source from emergency melting
      rainrate(i) = rainrate(i) + drainrate_melt(i)
    end do

    ! If emergency melting has removed all of the ice fall-flux,
    ! reset the fraction of falling-out ice to zero as well.
    if ( l_mcr_qcf2 ) then
      do i = 1, points
        if ( snow_cry(i) + snow_agg(i) <= zero )  frac_ice_fall(i) = zero
      end do
    else
      do i = 1, points
        if ( snow_agg(i) <= zero )  frac_ice_fall(i) = zero
      end do
    end if

  end if  ! ( l_mcr_precfrac )

end if  ! ( .not. l_proc_fluxes )


if ( l_mcr_precfrac ) then
  ! Tidy-ups for the prognostic precip fraction...

  ! If the precip mass or flux has gone to zero, reset the associated
  ! precip fraction to zero as well
  if ( l_subgrid_graupel_frac ) then
    ! Precip fraction includes graupel
    do i = 1, points
      if ( qrain(i) + qgraup(i) <= zero )        precfrac_k(i) = zero
      if ( rainrate(i) + grauprate(i) <= zero )  precfrac_fall(i) = zero
    end do
  else
    ! Precip fraction doesn't include graupel
    do i = 1, points
      if ( qrain(i) <= zero )     precfrac_k(i) = zero
      if ( rainrate(i) <= zero )  precfrac_fall(i) = zero
    end do
  end if

  ! Reset the output rainfrac equal to the prognostic precip fraction
  ! for diagnostics.  Note: it is possible for all the rain left on
  ! the current model-level to evaporate, but still have rain falling
  ! through to the level below.  In this case, we do want a nonzero
  ! rainfrac for the diagnostic, even though qrain has gone to zero.
  ! So set rainfrac to max of precfrac_k and precfrac_fall.
  do i = 1, points
    rainfrac(i) = max( precfrac_k(i), precfrac_fall(i) )
  end do
  if ( l_subgrid_graupel_frac ) then
    ! If prognostic precip fraction includes graupel, reset rainfrac
    ! to zero where there is graupel but no rain
    do i = 1, points
      if ( qrain(i) <= zero .and. rainrate(i) <= zero )  rainfrac(i) = zero
    end do
  end if

end if  ! ( l_mcr_precfrac )

! Deallocate water tracer arrays
if (l_wtrac) then
  deallocate(wtrac_mp_cpr_old%depos_l)
  deallocate(wtrac_mp_cpr_old%depos_v)
  deallocate(wtrac_mp_cpr_old%hom_qr)
  deallocate(wtrac_mp_cpr_old%het_q)
  deallocate(wtrac_mp_cpr_old%qchange)
  deallocate(wtrac_mp_cpr_old%droplet_flux)
  deallocate(wtrac_mp_cpr_old%t)
  deallocate(wtrac_mp_cpr_old%qrain)
  deallocate(wtrac_mp_cpr_old%qcf)
  deallocate(wtrac_mp_cpr_old%qcl)
  deallocate(wtrac_mp_cpr_old%q)
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_tidy
end module lsp_tidy_mod
