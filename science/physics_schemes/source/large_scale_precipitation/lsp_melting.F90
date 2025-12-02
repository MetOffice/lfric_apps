! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Melting of ice particles
! Subroutine Interface:
module lsp_melting_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_MELTING_MOD'

contains

subroutine lsp_melting(                                                        &
  points, timestep,                                                            &
                                          ! Number of points and tstep
  q, q_ice, qgraup, qcf, qcft, qrain, qsl,                                     &
                                          ! Water contents
  t, p,                                                                        &
                                          ! Temperature and pressure
  area_liq, area_mix, area_ice, area_clear,                                    &
                                          ! Partition information
  cficei, frac_ice_above,                                                      &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
  rainfrac, rain_liq, rain_mix,                                                &
                                          ! Rain fraction information
  rain_ice, rain_clear, rain_new,                                              &
  cf, cff,                                                                     &
                                          ! Current cloud fractions for
                                          ! updating
  rho, rhor, m0, tcg, tcgi,                                                    &
                                          ! Parametrization information
  corr2, rocor, ice_nofall,                                                    &
  lfrcp, ice_type,                                                             &
                                          ! Microphysical information
  l_psd,                                                                       &
                                          ! Code options
  ptransfer,                                                                   &
                                          ! Mass transfer diagnostic
  one_over_tsi,                                                                &
                                          ! 1/(timestep*iterations)
  cftransfer,cfftransfer,                                                      &
                                          ! Cloud transfer diagnostics
  rftransfer,                                                                  &
                                          ! Rain fraction transfer
  l_use_agg_vt,                                                                &
  dqprec_ice, dqprec_mix, dqprec_new,                                          &
                                          ! Total precip incs within ice-cloud
  precfrac_k,                                                                  &
                                          ! Prognostic precipitation fraction
  wtrac_mp_cpr_old                                                             &
                                          ! Store used for water tracers
  )

!Use in reals in lsprec precision, both microphysics related and general atmos
use lsprec_mod, only: cx, constp, tw1, tw2, tw3, tw4, tw5, zerodegc,           &
                      zero, one, small_number

  ! Microphysics modules- logicals and integers
use mphys_constants_mod, only: ice_type_offset
use mphys_inputs_mod,    only: l_diff_icevt,                                   &
                               l_mcr_precfrac, l_subgrid_graupel_frac,         &
                               i_update_precfrac, i_homog_areas, i_sg_correl

! Water tracers
use free_tracers_inputs_mod, only: l_wtrac
use wtrac_mphys_mod,         only: mp_cpr_old_wtrac_type

! Dr Hook Modules
use yomhook,             only: lhook, dr_hook
use parkind1,            only: jprb, jpim

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,             only: real_lsprec

! Large scale precipitation modules
use lsp_moments_mod, only: lsp_moments
use lsp_combine_precfrac_mod, only: lsp_combine_precfrac

implicit none

! Purpose:
!   Update cloud prognostics as a result of melting of ice particles

! Method:
!   Solve implicitly the microphysical transfer equation for a
!   specified distribution of ice particles melting at an
!   approximated wet-bulb temperature
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.


! Subroutine Arguments

integer, intent(in) ::                                                         &
  points,                                                                      &
                        ! Number of points to calculate
  ice_type          ! Type of ice (0 - crystals, 1 - aggregates
                        !              3 - graupel)

real (kind=real_lsprec), intent(in) ::                                         &
  timestep,                                                                    &
                        ! Timestep / s
  q(points),                                                                   &
                        ! Gridbox mean vapour content / kg kg-1
  q_ice(points),                                                               &
                        ! Vapour content in ice partition / kg kg-1
  qgraup(points),                                                              &
                        ! Graupel content / kg kg-1
  qcft(points),                                                                &
                        ! Total ice water content / kg kg-1
  qsl(points),                                                                 &
                        ! Saturated mixing ratio wrt liquid / kg kg-1
  p(points),                                                                   &
                        ! Pressure / N m-2
  area_liq(points),                                                            &
                        ! Fraction of gridbox with liquid-only cloud
  area_mix(points),                                                            &
                        ! Fraction of gridbox with mixed phase cloud
  area_ice(points),                                                            &
                        ! Fraction of gridbox with ice-only cloud
  area_clear(points),                                                          &
                        ! Frac of gridbox with no cloud
  cficei(points),                                                              &
                        ! 1/Fraction of gridbox with ice cloud
  frac_ice_above(points),                                                      &
                             ! Ice fraction in layer above this layer
  rho(points),                                                                 &
                        ! Air density / kg m-3
  rhor(points),                                                                &
                        ! 1/Air density / m3 kg-1
  corr2(points),                                                               &
                        ! Temperature dependency of diffusion params
  rocor(points),                                                               &
                        ! Combination of corr and corr2
  m0,                                                                          &
                        ! Seed ice water content / kg kg-1
  tcg(points),                                                                 &
                        ! T dependent function in ice size distribution
  tcgi(points),                                                                &
                        ! 1/tcg (no units)
  ice_nofall(points),                                                          &
                        ! Fraction of qcf that is not falling out
  lfrcp,                                                                       &
                        ! Latent heat of fusion
                        ! / heat capacity of air / K
  one_over_tsi          ! 1/(timestep*iterations)

real (kind=real_lsprec), intent(in out) ::                                     &
  qcf(points),                                                                 &
                        ! Ice water content of this category / kg kg-1
  qrain(points),                                                               &
                        ! Rain mixing ratio / kg kg-1
  t(points),                                                                   &
                        ! Temperature / K
  cf(points),                                                                  &
                        ! Current cloud fraction
  cff(points),                                                                 &
                        ! Current ice cloud fraction
  rainfrac(points),                                                            &
                        ! Rain fraction
  rain_liq(points),                                                            &
                        ! Overlap fraction of rain and liquid
  rain_mix(points),                                                            &
                        ! Overlap fraction of rain and mixed phase
  rain_ice(points),                                                            &
                        ! Overlap fraction of rain and ice
  rain_clear(points),                                                          &
                        ! Overlap fraction of rain and clear sky
  rain_new(points),                                                            &
                        ! Area where rain produced in air that had no rain
  ptransfer(points),                                                           &
                        ! Mass rimed in this timestep / kg kg-1
  rftransfer(points),                                                          &
                        ! Transfer of rain fraction this timestep
  cftransfer(points),                                                          &
                         ! Cumulative cloud fraction increment
  cfftransfer(points)! Cumulative liquid cloud fraction increment

logical, intent(in) ::                                                         &
  l_psd,                                                                       &
                        ! Use generic ice particle size distribution
  l_use_agg_vt(points)
                        ! Determines which vt-D parameters to use

! Stuff for updating the prognostic precip fraction;
! increments to rain+graupel in the mixed-phase and ice-only cloud regions
real(kind=real_lsprec), intent(in out) :: dqprec_mix(points)
real(kind=real_lsprec), intent(in out) :: dqprec_ice(points)
! Increment in air which didn't have any rain already
real(kind=real_lsprec), intent(in out) :: dqprec_new(points)
! Prognostic precip fraction
real(kind=real_lsprec), intent(in out) :: precfrac_k(points)

type(mp_cpr_old_wtrac_type), intent(in out) :: wtrac_mp_cpr_old
                        ! Water tracers store

! Local Variables

integer ::                                                                     &
  i,                                                                           &
                        ! Loop counter for points
  cry_offset        ! Index offset for ice crystals

real (kind=real_lsprec) ::                                                     &
  tempw,                                                                       &
                        ! Average supersaturation in gridbox
  temp7(points),                                                               &
                        ! Wet bulb temperature minus 0 C / deg C
  pr02,                                                                        &
                        ! Temporary in melting mass calculation
  dpr(points),                                                                 &
                        ! Mass melted this timestep / kg kg-1
  delta_cf,                                                                    &
                        ! Cloud fraction transfered (no units)
  delta_cff,                                                                   &
                        ! Ice cloud fraction transfered (no units)
  rf_final,                                                                    &
                        ! Final rain fraction (no units)
  m_1(points),                                                                 &
                        ! 1st moment of generic particle size dist.
  m_0p5_dip3(points),                                                          &
                        ! 1+(di+1)/2 moment of generic PSD
  m_0p5_dicp3(points)
                        ! 1+(dic+1)/2 moment of generic PSD

! Ice mass that is not falling out
real (kind=real_lsprec) :: qcf_nofall(points)

! Area fraction of rain produced by melting where no pre-existing rain
real(kind=real_lsprec) :: area_new
! Work variable for computing fractions of the increment in different areas
real(kind=real_lsprec) :: tmp
! Area fraction where rain created
real(kind=real_lsprec) :: area_inc(points)
! Precip content (rain+graupel) before the increment is applied
real(kind=real_lsprec) :: qprec(points)

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

character(len=*), parameter :: RoutineName='LSP_MELTING'


    !-----------------------------------------------
    ! Select appropriate ice parametrization (see c_lspsiz)
    !-----------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! When working in 32-bit, a Cray compiler bug on power breaks PROC
! comparability.  Conditionally using NOVECTOR makes this go
! away. This directive will be valid until the end of the program unit
! or when a VECTOR directive (without any extra string) is used.
! Review upon upgrade.
#if defined(LSPREC_32B) && defined (CRAYFTN_VERSION) && (CRAYFTN_VERSION <8004000)
!DIR$ NOVECTOR
#endif

cry_offset=ice_type*ice_type_offset

! Set amount of qcf that is not falling ouot this timestep
do i = 1, points
  qcf_nofall(i) = qcf(i) * ice_nofall(i)
end do

    ! Use the generic ice particle size distribution
    ! Calculate the 1st (cx(84)) moment and the 1+0.5(di+1) (cx(85))
    ! moment of the ice particle size distribution
if (l_psd) then

  call lsp_moments(points,rho,t,qcf_nofall,cficei,cx(84),m_1)
  if (.not. l_diff_icevt) then
          ! Only one set of vt-D parameters

    call lsp_moments(points,rho,t,qcf_nofall,cficei,cx(85),m_0p5_dip3)
  else
    ! The vt-D relation which gives the
    ! least mass weighted mean fallspeed will be used so
    ! calculate both required moments
    call lsp_moments(points,rho,t,qcf_nofall,cficei,cx(185),m_0p5_dicp3)
                        ! psd moment with crystal parameters
    call lsp_moments(points,rho,t,qcf_nofall,cficei,cx(85),m_0p5_dip3)
                        ! psd moment with aggregate parameters
  end if
end if


! Identify the points that need to be calculated.
npts = 0
do i = 1, points

  if (qcf_nofall(i) >  m0 .and. t(i) >  zerodegc) then

    npts = npts + 1
    ix(npts) = i

  else

    dpr(i) = zero

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
      ! Calculate the average wet-bulb temperature (tempw)
      !-----------------------------------------------
      ! First calculate the average supersaturation
  if (ice_type  ==  3) then
        ! Graupel does not use ice cloud fraction
    tempw = max(qsl(i)-q(i),zero)
  else
        ! Use the full ice cloud fraction information
    tempw=area_ice(i)*max(qsl(i)-q_ice(i),zero)*cficei(i)
  end if

      ! Now calculate the average amount (temp7) by which the
      ! wet-bulb temperature exceeds 0 degrees C.
  temp7(i) = t(i) - zerodegc - tempw                                           &
           * (tw1 + tw2 * (p(i)-tw3) - tw4*(t(i)-tw5) )
  temp7(i) = max(temp7(i),zero)


end do

if (l_psd) then
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do c = 1, npts

    i = ix(c)

    !-----------------------------------------------
    ! Calculate melting rate
    !-----------------------------------------------
    ! Use generic particle size distribution
    if (.not. l_diff_icevt) then
      ! constp(84) = ventilation coefficient 1
      ! constp(85) = ventilation coefficient 2
      ! constp(90) = 2 pi axial_ratio_correction
      !              * air_conductivity / Lf
      dpr(i) = rhor(i) * constp(90) * timestep *                               &
                 ( constp(84)*m_1(i)*corr2(i) +                                &
                 constp(85)*rocor(i)*m_0p5_dip3(i) )
    else
      if (l_use_agg_vt(i)) then
        ! Use aggregate parameters for ventilation
        dpr(i) = rhor(i) * constp(90) * timestep *                             &
                 ( constp(84)*m_1(i)*corr2(i) +                                &
                 constp(85)*rocor(i)*m_0p5_dip3(i) )
      else
        ! Use crystal parameters for ventilation
        dpr(i) = rhor(i) * constp(90) * timestep *                             &
                 ( constp(84)*m_1(i)*corr2(i) +                                &
                 constp(185)*rocor(i)*m_0p5_dicp3(i) )
      end if
    end if
  end do
else
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do c = 1, npts

    i = ix(c)

    ! Use particle size distribution based on intercepts
    if (ice_type  ==  3) then
      ! Graupel does not use ice cloud fraction
      pr02 = rho(i)*qcf_nofall(i)*constp(5+cry_offset)*tcgi(i)
    else
      ! Use the full ice cloud information
      pr02 = rho(i)*qcf_nofall(i)*cficei(i) * constp(5+cry_offset)*tcgi(i)
    end if  ! ice_type==3

    dpr(i)  = tcg(i) * constp(14+cry_offset) * timestep *                      &
              (constp(7+cry_offset) * corr2(i) *                               &
              pr02**cx(4+cry_offset)                                           &
            + constp(8+cry_offset) * rocor(i) *                                &
              pr02**cx(5+cry_offset)) * rhor(i)

  end do
end if  ! l_psd

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)

      !-----------------------------------------------
      ! Calculate transfer
      !-----------------------------------------------
      ! Solve implicitly in terms of temperature
  dpr(i) = temp7(i) * (one-one/(one+dpr(i)*lfrcp))/lfrcp

      ! Limit mass transfer to the mass available
  dpr(i) = min(dpr(i),qcf(i))

  ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

      !------------------------------------------------
      ! Adjust ice and rain contents
      !------------------------------------------------
  qcf(i)   = qcf(i)   - dpr(i)
  qrain(i) = qrain(i) + dpr(i)
  t(i)     = t(i)     - dpr(i) * lfrcp

end do

if (l_wtrac) then     ! Store phase change amount for water tracer use
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do c = 1, npts
    i = ix(c)
    wtrac_mp_cpr_old%qchange(i) = dpr(i)
  end do
end if

if ( .not. l_mcr_precfrac ) then
  !------------------------------------------------
  ! Update rain fractions
  !------------------------------------------------
  ! (not doing this here if using prognostic precip fraction,
  !  as in that case rainfrac is updated later)

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do c = 1, npts

    i = ix(c)


    if (dpr(i)  >   zero) then
      rf_final = max(max(rainfrac(i),cff(i)),frac_ice_above(i))

      rftransfer(i) = rftransfer(i) + (rf_final - rainfrac(i))

      rainfrac(i)= rf_final
      rain_liq(i)= min(area_liq(i),rainfrac(i))
      rain_mix(i)= min(area_mix(i),rainfrac(i)-rain_liq(i))
      rain_ice(i)=                                                             &
           min(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i))
      rain_clear(i)=                                                           &
             rainfrac(i) - rain_liq(i)-rain_mix(i)-rain_ice(i)
      rain_clear(i) = min(area_clear(i),rain_clear(i))

    end if  ! dpr  >   0

  end do

end if  ! ( .not. l_mcr_precfrac )

!------------------------------------------------
! Update cloud fractions
!------------------------------------------------
if (ice_type /= 3) then
  ! ( ice_type /=3 as graupel does not change
  !   cloud fractions )
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do c = 1, npts

    i = ix(c)

    ! Assume that ice cloud fraction is reduced in
    ! proportion to the mass removed
    delta_cff = -cff(i) * dpr(i) / qcft(i)

        ! Assume a random overlap between the removed ice
        ! cloud and the liquid cloud
    delta_cf  = delta_cff * area_ice(i) * cficei(i)

        ! Calculate transfer rate diagnostics
    cftransfer(i)  = cftransfer(i)                                             &
                     + delta_cf  * one_over_tsi
    cfftransfer(i) = cfftransfer(i)                                            &
                     + delta_cff * one_over_tsi

    cf(i)  = cf (i) + delta_cf
    cff(i) = cff(i) + delta_cff

  end do

end if  ! ice_type /= 3

if ( l_mcr_precfrac .and. ( .not. (ice_type==3 .and. l_subgrid_graupel_frac) ) &
   ) then
  ! If using prognostic precipitation fraction...
  ! Need to store contribution of melting of ice to the
  ! total precip increment in the ice-cloud partition.
  ! (not needed if this is melting of graupel, and including graupel
  !  as part of the prognostic precip fraction, since in that case
  !  converting graupel to rain makes no difference to the precip mass
  !  within the ice-cloud partition).

  if ( i_update_precfrac == i_homog_areas ) then

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do c = 1, npts
      i = ix(c)
      ! Assuming that all the melting source occurs within the ice-cloud.
      ! If part of the ice-cloud area doesn't already contain rain/graupel,
      ! need to assign part of the precip mass increment to "new" area...

      ! Compute area of newly-created rain
      area_new = max( zero, area_mix(i) - rain_mix(i) )                        &
               + max( zero, area_ice(i) - rain_ice(i) )

      ! Compute 1.0 / total area producing rain by melting
      tmp = one / max( rain_mix(i) + rain_ice(i) + area_new, small_number )

      ! Increment applies in the rain_ice, rain_mix  and new area partitions;
      ! divvy it up between them in proportion to area
      dqprec_mix(i) = dqprec_mix(i) + dpr(i) * rain_mix(i) * tmp
      dqprec_ice(i) = dqprec_ice(i) + dpr(i) * rain_ice(i) * tmp
      dqprec_new(i) = dqprec_new(i) + dpr(i) * area_new * tmp

      ! Update area of new precip passed back up
      rain_new(i) = max( rain_new(i), area_new )

    end do

  else if ( i_update_precfrac == i_sg_correl ) then

    ! Set precip mass before the increment was applied
    ! (ignoring negative values), and area fraction of the increment
    do i = 1, points
      qprec(i) = max( qrain(i) - dpr(i), zero )
      dpr(i) = max( qrain(i), zero ) - qprec(i)
      area_inc(i) = area_mix(i) + area_ice(i)
    end do
    if ( l_subgrid_graupel_frac ) then
      do i = 1, points
        qprec(i) = qprec(i) + max( qgraup(i), zero )
      end do
    end if

    ! Calculate combined fraction of the existing precip and increment
    ! (precfrac_k is updated)
    call lsp_combine_precfrac( points,                                         &
                               qprec, dpr, precfrac_k, area_inc )

  end if

end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_melting
end module lsp_melting_mod
