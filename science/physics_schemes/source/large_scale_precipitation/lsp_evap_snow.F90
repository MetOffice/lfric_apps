! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Evaporation of melting snow
! Subroutine Interface:
module lsp_evap_snow_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_EVAP_SNOW_MOD'

contains

subroutine lsp_evap_snow(                                                      &
  points, timestep,                                                            &
                                          ! Number of points and tstep
  q, q_ice, qcf, qcft, t, p,                                                   &
                                          ! Water contents, temp, pres
  esw, qsl,                                                                    &
                                          ! Saturated quantities
  area_ice, cficei,                                                            &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
  cf, cff,                                                                     &
                                          ! Current cloud fractions for
                                          ! updating
  rho, tcg, tcgi,                                                              &
                                          ! Parametrization information
  corr2, rocor, ice_nofall, lheat_correc_liq,                                  &
  lsrcp, ice_type,                                                             &
                                          ! Microphysical information
  l_psd,                                                                       &
                                          ! Code options
  ptransfer,                                                                   &
                                          ! Mass transfer diagnostic
  one_over_tsi,                                                                &
                                          ! 1/(timestep*iterations)
  cftransfer,cfftransfer,                                                      &
                                          ! Cloud transfer diagnostics
  l_use_agg_vt,                                                                &
                                          ! Fallspeed branch choice
  wtrac_mp_cpr_old                                                             &
                                          ! Store used for water tracers
  )

!Use in reals in lsprec precision, both microphysics related and general atmos
use lsprec_mod, only: apb4, apb5, apb6, m0, cx, constp, zerodegc,              &
                      zero

  ! Microphysics Modules- logicals and integers
use mphys_constants_mod, only: ice_type_offset
use mphys_inputs_mod,    only: l_diff_icevt

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,             only: real_lsprec

! Water tracers
use free_tracers_inputs_mod, only: l_wtrac
use wtrac_mphys_mod,         only: mp_cpr_old_wtrac_type

! Dr Hook Modules
use yomhook,             only: lhook, dr_hook
use parkind1,            only: jprb, jpim

! Large scale precipitation modules
use lsp_moments_mod,     only: lsp_moments

implicit none

! Purpose:
!   Update cloud prognostics as a result of sublimation of melting
!   snow

! Method:
!   Solve the microphysical transfer equation for a specified
!   distribution of ice particles in a distribution of vapour
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

! Source for vapour. Sink for ice.

! Subroutine Arguments

integer, intent(in) ::                                                         &
  points,                                                                      &
                        ! Number of points to calculate
  ice_type          ! Type of ice (0 - crystals, 1 - aggregates)

real (kind=real_lsprec), intent(in) ::                                         &
  timestep,                                                                    &
                        ! Timestep / s
  p(points),                                                                   &
                        ! Air pressure / N m-2
  q_ice(points),                                                               &
                        ! Vapour content in ice partition / kg kg-1
  esw(points),                                                                 &
                        ! Saturated vapour pres. over liquid / N m-2
  qsl(points),                                                                 &
                        ! Saturated humidity wrt liquid / kg kg-1
  area_ice(points),                                                            &
                        ! Fraction of gridbox with ice but no liquid
  cficei(points),                                                              &
                        ! 1/Fraction of gridbox with ice cloud
  rho(points),                                                                 &
                        ! Air density / kg m-3
  tcg(points),                                                                 &
                        ! T dependent function in ice size distribution
  tcgi(points),                                                                &
                        ! 1/tcg (no units)
  corr2(points),                                                               &
                        ! Temperature correction factor (no units)
  rocor(points),                                                               &
                        ! Combined fall and corr2 correction factor
  ice_nofall(points),                                                          &
                        ! Fraction of qcf that is not falling out
  lheat_correc_liq(points),                                                    &
                               ! Liquid latent heat correction factor
  lsrcp,                                                                       &
                        ! Latent heat of sublimation
                        ! /heat capacity of air (cP) / K
  one_over_tsi          ! 1/(timestep*iterations)

real (kind=real_lsprec), intent(in out) ::                                     &
  q(points),                                                                   &
                        ! Vapour content / kg kg-1
  qcf(points),                                                                 &
                        ! Ice water content in ice category to be
!                           updated    / kg kg-1
    qcft(points),                                                              &
                          ! Ice water in all ice categories
!                           (for cloud fraction calculations)
    t(points),                                                                 &
                          ! Temperature / K
    cf(points),                                                                &
                          ! Current cloud fraction
    cff(points),                                                               &
                          ! Current ice cloud fraction
    ptransfer(points)  ! Mass deposited in this timestep / kg kg-1

real (kind=real_lsprec), intent(in out) ::                                     &
  cftransfer(points),                                                          &
                         ! Cloud fraction increment this tstep
  cfftransfer(points)! Ice cloud fraction inc this tstep

logical, intent(in) ::                                                         &
  l_psd,                                                                       &
                        ! Use generic ice particle size distribution
  l_use_agg_vt(points)
                        ! Determines which vt-D parameters to use

type(mp_cpr_old_wtrac_type), intent(in out) :: wtrac_mp_cpr_old
                        ! Water tracers store

! Local Variables

integer ::                                                                     &
  i,                                                                           &
                        ! Loop counter for points
  cry_offset        ! Index offset for ice crystals

real (kind=real_lsprec) ::                                                     &
  pr02,                                                                        &
                        ! Temporary in calculation of PSD slopes
  pr04,                                                                        &
                        ! Temporary in calculation of PSD slopes
  dpr,                                                                         &
                        ! Temporary in calculating ice transfers
  deltacf,                                                                     &
                        ! Change in cf across timestep
  tempw,                                                                       &
                        ! Available subsaturation for evaporation
  m_1(points),                                                                 &
                        ! 1st moment of the generic ice size distribtn
  m_0p5_dip3(points),                                                          &
                        ! 1+(di+1)/2 moment of the generic ice PSD
  m_0p5_dicp3(points)
                        ! 1+(dic+1)/2 moment of the generic ice PSD

! Amount of qcf that is not falling out
real (kind=real_lsprec) :: qcf_nofall(points)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_EVAP_SNOW'


    !-----------------------------------------------
    ! Select appropriate ice parametrization (see c_lspsiz)
    !-----------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

cry_offset=ice_type*ice_type_offset

! Set qcf amount that is not falling out
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

do i = 1, points

  if (qcf_nofall(i) >  m0 .and. t(i) >  zerodegc) then

        !-----------------------------------------------
        ! Diffusional growth parameters
        !-----------------------------------------------
    pr04 = ((apb4-apb5*t(i))*esw(i)+apb6*p(i)*t(i)**3)

        !-----------------------------------------------
        ! Calculate transfer rates
        !-----------------------------------------------
    if (l_psd) then
          ! Use generic particle size distribution
      if (.not. l_diff_icevt) then
         ! constp(83) = 2 pi axial_ratio_correction
         ! constp(84) = ventilation coefficient 1
         ! constp(85) = ventilation coefficient 2
         !        * Sc^(1/3)*ci^0.5/viscosity0^0.5
        dpr = constp(83) * timestep * t(i)**2 * esw(i)                         &
                        * (constp(84)*m_1(i)*corr2(i)                          &
                        +  constp(85)*rocor(i)*m_0p5_dip3(i))                  &
                        / (qsl(i) * rho(i) * pr04)
      else
        if (l_use_agg_vt(i)) then
           ! Use aggregate parameters for ventilation
          dpr = constp(83) * timestep * t(i)**2 * esw(i)                       &
                       * (constp(84)*m_1(i)*corr2(i)                           &
                       +  constp(85)*rocor(i)*m_0p5_dip3(i))                   &
                       / (qsl(i) * rho(i) * pr04)
        else
           ! Use crystal parameters for ventilation
          dpr = constp(83) * timestep * t(i)**2 * esw(i)                       &
                       * (constp(84)*m_1(i)*corr2(i)                           &
                       +  constp(185)*rocor(i)*m_0p5_dicp3(i))                 &
                       / (qsl(i) * rho(i) * pr04)
        end if
      end if
    else
          ! Use particle size distribution based on intercepts
      pr02=rho(i)*qcf_nofall(i)*cficei(i)*constp(5+cry_offset)*tcgi(i)
      dpr=tcg(i)*constp(6+cry_offset)*t(i)**2*esw(i)*timestep*                 &
      (constp(7+cry_offset)*corr2(i)*pr02**cx(4+cry_offset)                    &
       +constp(8+cry_offset)*rocor(i)*pr02**cx(5+cry_offset))                  &
       /(qsl(i)*rho(i)*pr04)

    end if  ! l_psd

        !-----------------------------------------------
        ! Limit transfers
        !-----------------------------------------------
        ! tempw is the subsaturation that is available
    tempw = area_ice(i) * (qsl(i) - q_ice(i))
    dpr = dpr * tempw
    dpr = max(min(dpr,tempw*lheat_correc_liq(i)),zero)

        ! Limit on the amount of ice available
    dpr = min(dpr,qcf(i))

        !-----------------------------------------------
        ! Store process rate
        !-----------------------------------------------
    ptransfer(i) = ptransfer(i) + dpr * one_over_tsi

        !-----------------------------------------------
        ! Store and update ice cloud fraction and total fractions
        !-----------------------------------------------

    deltacf     = ( -cff(i) * dpr/qcft(i) )
    cftransfer(i)  = cftransfer(i) + deltacf * one_over_tsi
    cfftransfer(i) = cfftransfer(i)+ deltacf * one_over_tsi

    cf(i)  = cf(i)  + deltacf
    cff(i) = cff(i) + deltacf


        !-----------------------------------------------
        ! Update values of ice and vapour
        !-----------------------------------------------

    qcf(i) = qcf(i) - dpr
    q(i)   = q(i)   + dpr
    t(i)   = t(i)   - dpr*lsrcp
    if (l_wtrac)  wtrac_mp_cpr_old%qchange(i) = dpr

  end if  ! qcf_nofall gt m0 etc.

end do  ! Points

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_evap_snow
end module lsp_evap_snow_mod
