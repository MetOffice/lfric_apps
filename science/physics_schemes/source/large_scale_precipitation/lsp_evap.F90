! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Evaporation of rain
! Subroutine Interface:
module lsp_evap_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_EVAP_MOD'

contains

subroutine lsp_evap(                                                           &
  points, timestep,                                                            &
                                          ! Number of points and tstep
  p, q, qrain, t, qgraup, q_ice, q_clear,                                      &
                                          ! Water contents and temp
  rainfrac, rain_liq, rain_mix,                                                &
                                          ! Rain fractions for updating
  rain_ice, rain_clear,                                                        &
  rho, corr, corr2,rocor,dhir,                                                 &
                                          ! Parametrization information
  lcrcp, lheat_correc_liq, qsl, esw, rain_nofall,                              &
  ptransfer, rftransfer,                                                       &
                                          ! Mass transfer diagnostic
  one_over_tsi,                                                                &
                                          ! 1/(timestep*iterations)
  dqprec_ice, dqprec_clear,                                                    &
                                          ! Total precip incs within partitions
  precfrac_k,                                                                  &
                                          ! Prognostic precip fraction
  wtrac_mp_cpr_old                                                             &
                                          ! Store used for water tracers
  )

!Use in reals in lsprec precision, both microphysics related and general atmos
use lsprec_mod,         only: apb4, apb5, apb6, qcfmin, m0, cx, constp,        &
                              rho_q_veloc, lam_evap_enh, max_as_enh,           &
                              zero, one

  ! Microphysics modules
use mphys_inputs_mod,    only: l_warm_new, l_mcr_qrain, l_mcr_precfrac,        &
                               l_subgrid_graupel_frac, heavy_rain_evap_fac

use lsp_evap_precfrac_mod, only: lsp_evap_precfrac

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,             only: real_lsprec

! Water tracers
use free_tracers_inputs_mod, only: l_wtrac
use wtrac_mphys_mod,         only: mp_cpr_old_wtrac_type

! Dr Hook Modules
use yomhook,             only: lhook, dr_hook
use parkind1,            only: jprb, jpim

implicit none

! Purpose:
!   Update rain and water vapour due to evaporation

! Method:
!   Integrate evaporation rates of a single raindrop over the
!   raindrop size distribution
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

integer, intent(in) ::                                                         &
  points
                        ! Number of points to calculate

real (kind=real_lsprec), intent(in) ::                                         &
  timestep,                                                                    &
                        ! Timestep / s
  p(points),                                                                   &
                        ! Air pressure / N m-2
  qgraup(points),                                                              &
                        ! Graupel content / kg kg-1
  q_ice(points),                                                               &
                        ! Local vapour in ice partition / kg kg-1
  q_clear(points),                                                             &
                        ! Local vapour in clear partition / kg kg-1
  lcrcp,                                                                       &
                        ! Latent heat of condensation/cP / K
  dhir(points),                                                                &
                        ! layer thickness/timestep / m

  corr(points),                                                                &
                        ! Fall speed correction factor due to
                        !   air density changes (no units)
  corr2(points),                                                               &
                        ! Air diffusivity correction (no units)
  rho(points),                                                                 &
                        ! Air density / kg m-3
  lheat_correc_liq(points),                                                    &
                                ! Correction of subsaturation due
                                ! to latent heat
  qsl(points),                                                                 &
                        ! Saturated spec humidity wrt liquid / kg kg-1
  esw(points),                                                                 &
                        ! Saturated vapour pressure wrt liquid / N m-2
  rocor(points),                                                               &
                        ! Air density and viscosity correction factor
                        !   (no units)
  rain_nofall(points),                                                         &
                        ! Fraction of the rain-mass that is not falling out
  one_over_tsi
                        ! 1/(timestep*iterations)

real (kind=real_lsprec), intent(in out) ::                                     &
  qrain(points),                                                               &
                        ! Rain water content / kg kg-1
  q(points),                                                                   &
                        ! Vapour content / kg kg-1
  t(points),                                                                   &
                        ! Temperature / K
  ptransfer(points),                                                           &
                        ! Evaporation rate / kg kg-1 s-1
  rftransfer(points),                                                          &
                        ! Rate of change of rain fraction / s-1
  rainfrac(points),                                                            &
                        ! Rain fraction (no units)
  rain_liq(points),                                                            &
                        ! Overlap of rain with liquid (no units)
  rain_mix(points),                                                            &
                        ! Overlap of rain with mixed phase region
  rain_ice(points),                                                            &
                        ! Overlap of rain with ice
  rain_clear(points)! Overlap of rain with clear sky

! Stuff for updating the prognostic precip fraction;
! increments to rain+graupel in the ice-only cloud and cloud-free regions
real(kind=real_lsprec), intent(in out) :: dqprec_ice(points)
real(kind=real_lsprec), intent(in out) :: dqprec_clear(points)

! Prognostic precipitation fraction
real(kind=real_lsprec), intent(in out) :: precfrac_k(points)

type(mp_cpr_old_wtrac_type), intent(in out) :: wtrac_mp_cpr_old
                        ! Water tracers store

!     Real, Intent(Out) ::

! Local Variables

integer ::                                                                     &
  i                 ! Loop counter

real (kind=real_lsprec) ::                                                     &
  dpr(points),                                                                 &
                        ! Amount of mass evaporated / kg kg-1
  pr04(points),                                                                &
                        ! Temporary in evaporation rate calc.
  lamr1,                                                                       &
                        ! Powers of drop size distribution slope
  lamr2,                                                                       &
                        !   lambda.
  lamr3,                                                                       &
  temp7,                                                                       &
                        ! Subsaturation in gridbox / kg kg-1

! Local variables for Abel & Shipway-style diagnostic rain evaporation
! Enhancement

    prelamnew,                                                                 &
                          ! Precursor to new lambda calculated for A/S
    lamnew,                                                                    &
                          ! Lambda for Abel and Shipway calculations
    standard_as_veloc,                                                         &
                          ! Standard velocity for A/S
                          ! determined from rho_q_veloc input
    true_as_veloc,                                                             &
                          ! Actual Abel and Shipway velocity
                          ! determined from lamnew
    evap_enh_factor
                          ! Factor with which to increase (enhance)
                          ! the evaporation rate

! Local compression variable
integer ::                                                                     &
  npts,                                                                        &
                        ! Number of compressed points to compute
  c,                                                                           &
                        ! Compressed point index
  ix(points)
                        ! Original index


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_EVAP'

! End of header and all declarations

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! When working in 32-bit, a Cray compiler bug on power breaks PROC
! comparability.  Conditionally using NOVECTOR makes this go
! away. This directive will be valid until the end of the program unit
! or when a VECTOR directive (without any extra string) is used.
! Review upon upgrade.
#if defined(LSPREC_32B) && defined (CRAYFTN_VERSION) && (CRAYFTN_VERSION <8004000)
!DIR$ NOVECTOR
#endif

do i = 1, points
      !-----------------------------------------------
      ! If there is only a small amount of rain present, evaporate
      ! it completely to avoid numerical problems.
      !-----------------------------------------------
  if (qrain(i)*rain_nofall(i)  <   qcfmin .and. .not. l_warm_new .or.          &
      qrain(i)*rain_nofall(i) < m0 .and. l_warm_new) then
    ! original code allowed significant evaporation of rain in-cloud,
    ! which is unphysical. l_warm_new prevents this by lowering
    ! the threshold to define a "small amount"

        ! Evaporate all this rain
    dpr(i) = qrain(i)

    t(i)   = t(i) - lcrcp * dpr(i)
    q(i)   = q(i) + dpr(i)
    qrain(i) = zero

      ! Store evaporation rate
    ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

    if (l_wtrac) wtrac_mp_cpr_old%qchange(i) = dpr(i) !Store for water tracers

    if ( .not. l_mcr_precfrac ) then
      ! Update rain fractions
      ! (not doing this here if using prognostic precip fraction,
      !  as in that case rainfrac is updated later)
      rftransfer(i) = rftransfer(i)                                            &
                    - rainfrac(i) * one_over_tsi
      rainfrac(i)  = zero
      rain_liq(i)  = zero
      rain_mix(i)  = zero
      rain_ice(i)  = zero
      rain_clear(i)= zero
    end if

  end if  ! qrain lt qcfmin

end do

! Identify the points that need to be calculated.
npts = 0
do i = 1, points

  if ( qrain(i)  >   zero .and. rainfrac(i)  >   zero .and.                    &
       ( rain_clear(i) > zero .or. rain_ice(i) > zero ) ) then
    ! Can only evaporate if rain fraction overlaps with subsaturated air.

    npts = npts + 1
    ix(npts) = i

  end if

end do

! In all the following loops, the 'i' variable is used for the
! original index while the 'c' variable is used for the compressed
! index.

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)

      !-----------------------------------------------
      ! Calculate evaporation parameters
      !-----------------------------------------------
  pr04(i) = ((apb4-apb5*t(i))*esw(i)+apb6*p(i)*t(i)**3)

      !-----------------------------------------------
      ! Calculate powers of drop size distribution slopes
      !-----------------------------------------------
  if (l_mcr_qrain) then
        ! Rain is a mixing ratio - use fall speeds as parametrized
    lamr1 = qrain(i)*rain_nofall(i) * constp(50) * rho(i) / (rainfrac(i))
    lamr2 = lamr1 ** (cx(47)*cx(52))
    lamr3 = lamr1 ** (cx(49)*cx(52))
  else
        ! Rain is a diagnostic quantity
        ! - use fall speeds to allow sedimentation into next layer
    lamr1 = qrain(i) * rho(i) * dhir(i) /                                      &
               (constp(42) * corr(i) * rainfrac(i))
    lamr2 = lamr1 ** (cx(47)*cx(48))
    lamr3 = lamr1 ** (cx(49)*cx(48))
  end if  ! l_mcr_qrain

      !-----------------------------------------------
      ! Calculate evaporation rate
      !-----------------------------------------------
  dpr(i) = constp(46) * t(i)**2 * esw(i) * timestep
  dpr(i) = dpr(i) * ( (real(heavy_rain_evap_fac,real_lsprec)                   &
                                  * rocor(i) * lamr1)                          &
                    + (constp(47) * corr2(i) * lamr2)                          &
                    + (constp(48) * rocor(i) * lamr3) )


end do

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)

      !-----------------------------------------------
      ! Simulate  Abel & Shipway Section for diagnostic
      ! rain only
      !-----------------------------------------------

      ! Idea is that droplets are all small for big
      ! lambda and so for
      ! diagnostic rain they will not fall as far before
      ! evaporating.

  if (.not. l_mcr_qrain ) then

        ! Calculate new lambda for rain based on A/S (07)

    prelamnew = qrain(i)*rho(i)*dhir(i)                                        &
               /(rainfrac(i)*constp(42)*corr(i))

    lamnew =  one / (prelamnew**cx(48))

    if (lamnew  >  lam_evap_enh) then

          !Drops are smaller than mean of those at the
          !point at where the Abel and Shipway curve
          !diverges from the standard UM curve.

          !-----------------------------------------------------
          ! Calculate standard velocity
          !-----------------------------------------------------
      standard_as_veloc = rho_q_veloc / (rho(i)*qrain(i))

          !-----------------------------------------------------
          ! Calculate actual Abel and Shipway Bulk Velocity
          !-----------------------------------------------------
      true_as_veloc = (constp(57)/(rho(i)*qrain(i)))*                          &
                         (lamnew**cx(46)) *(                                   &
     ( constp(54) / ((lamnew+cx(56))**cx(59) ) ) +                             &
     ( constp(55) / ((lamnew+cx(57))**cx(60) ) )   )

          !-----------------------------------------------------
          ! Calculate evaporation enhancement factor
          ! and
          ! Limit enhancement factor to maximum value
          ! (Set in module)
          !-----------------------------------------------------
      evap_enh_factor = min((standard_as_veloc/true_as_veloc), max_as_enh)

          !-----------------------------------------------------
          ! Finally, enhance evaporation rate
          !-----------------------------------------------------
      dpr(i) = dpr(i)*evap_enh_factor

    end if   !lamnew >  Lam_evap_enh

  end if ! .not. l_mcr_qrain

end do

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)
      !-----------------------------------------------
      ! Calculate transfers
      !-----------------------------------------------
      ! Calculate gridbox mean supersaturation
  temp7 = (q_ice(i) - qsl(i)) * rain_ice(i)                                    &
         +(q_clear(i) - qsl(i)) * rain_clear(i)

      ! Limit on the gridbox mean supersaturation
  dpr(i) = dpr(i) * max(-temp7*lheat_correc_liq(i),zero)                       &
         /(qsl(i) * rho(i) * pr04(i) + dpr(i))

      ! Limit on the amount of rain available
  if ( rain_mix(i) > zero .or. rain_liq(i) > zero) then
    dpr(i) = min( dpr(i), qrain(i) *                                           &
            (rain_ice(i) + rain_clear(i)) / rainfrac(i))
  else
        ! No overlap between rain and liquid cloud, so
        ! (rain_ice + rain_clear) / rainfrac must be 1.
    dpr(i) = min( dpr(i), qrain(i))
  end if

end do

if ( l_mcr_precfrac ) then
  ! If using prognostic precipitation fraction, calculate precfrac increments
  call lsp_evap_precfrac( points, npts, ix, l_subgrid_graupel_frac,            &
                          rainfrac, rain_clear, rain_ice,                      &
                          dqprec_clear, dqprec_ice,                            &
                          qrain, qgraup, dpr, precfrac_k )
end if  ! ( l_mcr_precfrac )

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)

      !-----------------------------------------------
      ! Store process rate (kg/kg/s)
      !-----------------------------------------------
  ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

      !-----------------------------------------------
      ! Update values of rain, vapour and temperature
      !-----------------------------------------------
  qrain(i) = qrain(i) - dpr(i)
  q(i)     = q(i)     + dpr(i)
  t(i)     = t(i)     - dpr(i) * lcrcp

end do

if (l_wtrac) then     ! Store phase change amount for water tracer use
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do c = 1, npts
    i = ix(c)
    wtrac_mp_cpr_old%qchange(i) = dpr(i)
  end do
end if


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_evap
end module lsp_evap_mod
