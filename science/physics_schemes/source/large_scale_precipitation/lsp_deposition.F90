! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Deposition of ice particles
! Subroutine Interface:
module lsp_deposition_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_DEPOSITION_MOD'

contains

subroutine lsp_deposition(                                                     &
  points, timestep,                                                            &
                                          ! Number of points and tstep
  q, qcl, qcf, qcft, t, p,                                                     &
                                          ! Water contents, temp, pres
  q_ice_1, q_ice_2,                                                            &
                                          ! Subgrid-scale water contents
  area_ice_1, area_ice_2,                                                      &
                                          ! Subgrid-scale areas
  esi, qs, qsl,                                                                &
                                          ! Saturated quantities
  area_mix, cfliq, cfice, cficei, areamix_over_cfliq,                          &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
  cf, cfl, cff,                                                                &
                                          ! Current cloud fractions for
                                          ! updating
  rho, tcg, tcgi,                                                              &
                                          ! Parametrization information
  corr2, rocor, lheat_correc_ice, ice_nofall,                                  &
  lfrcp, lsrcp, ice_type,                                                      &
                                          ! Microphysical information
  l_psd,                                                                       &
                                          ! Code options
  ptransfer,                                                                   &
                                          ! Mass transfer diagnostic
  graut_psdep,                                                                 &
                                          ! Deposition rate on a single
                                          ! iteration for graupel
                                          ! autoconversion
  one_over_tsi,                                                                &
                                          ! 1/(timestep*iterations)
  cftransfer,cfltransfer,cfftransfer,                                          &
                                          ! Cloud transfer diagnostics
  l_use_agg_vt,                                                                &
                                          ! Vt-D branch to use
  wtrac_mp_cpr_old                                                             &
                                          ! Water tracers structure
  )

!Use in reals in lsprec precision, both microphysics related and general atmos
use lsprec_mod, only: apb1, apb2, apb3, m0, cx, constp, zerodegc, zero, half,  &
                      one

  ! Microphysics modules- logicals and integers
use mphys_constants_mod, only: ice_type_offset
use mphys_inputs_mod,    only: l_diff_icevt, l_proc_fluxes

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,             only: real_lsprec

! Dr Hook Modules
use yomhook,             only: lhook, dr_hook
use parkind1,            only: jprb, jpim

! Large scale precipitation modules
use lsp_moments_mod,     only: lsp_moments

use science_fixes_mod,   only: l_fix_gr_autoc

use free_tracers_inputs_mod, only: l_wtrac
use wtrac_mphys_mod,         only: mp_cpr_old_wtrac_type

implicit none

! Purpose:
!   Update cloud prognostics as a result of ice deposition and
!   sublimation

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

! Source for ice. Sink for liquid water and vapour.

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
  esi(points),                                                                 &
                        ! Saturated vapour pressure over ice / N m-2
  qs(points),                                                                  &
                        ! Saturated humidity wrt ice / kg kg-1
  qsl(points),                                                                 &
                        ! Saturated humidity wrt liquid / kg kg-1
  q_ice_1(points),                                                             &
                        ! Mean vapour in ice only region / kg kg-1
  q_ice_2(points),                                                             &
                        ! Mean vapour in clear region / kg kg-1
  area_ice_1(points),                                                          &
                        ! Ice only area that is growing by deposition
  area_ice_2(points),                                                          &
                        ! Ice only area that is subliming
  area_mix(points),                                                            &
                        ! Fraction of gridbox with mixed phase cloud
  cfliq(points),                                                               &
                        ! Fraction of gridbox with liquid cloud
  cfice(points),                                                               &
                        ! Fraction of gridbox with ice cloud
  cficei(points),                                                              &
                        ! 1/Fraction of gridbox with ice cloud
  areamix_over_cfliq(points),                                                  &
                        ! area_mix(points)/cfliq(points)
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
  lheat_correc_ice(points),                                                    &
                        ! Ice latent heat correction factor
  ice_nofall(points),                                                          &
                        ! Fraction of qcf that is not falling out
  lfrcp,                                                                       &
                        ! Latent heat of fusion
                        ! /heat capacity of air (cP) / K
  lsrcp,                                                                       &
                        ! Latent heat of sublimation/cP / K
  one_over_tsi          ! 1/(timestep*iterations)
logical, intent(in) ::                                                         &
   l_use_agg_vt(points)
                         ! Determines which vt-D parameters to use

real (kind=real_lsprec), intent(in out) ::                                     &
  q(points),                                                                   &
                        ! Vapour content / kg kg-1
  qcl(points),                                                                 &
                        ! Liquid water content / kg kg-1
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
    cfl(points),                                                               &
                          ! Current liquid cloud fraction
    cff(points),                                                               &
                          ! Current ice cloud fraction
    ptransfer(points),                                                         &
                          ! Mass deposited in this timestep / kg kg-1
    graut_psdep(points)
                          ! Deposition rate on a single
                          ! iteration for graupel
                          ! autoconversion

real (kind=real_lsprec), intent(in out) ::                                     &
  cftransfer(points),                                                          &
                         ! Cloud fraction increment this tstep
  cfltransfer(points),                                                         &
                         ! Liquid cloud fraction inc this tstep
  cfftransfer(points)! Ice cloud fraction inc this tstep

logical, intent(in) ::                                                         &
  l_psd
                        ! Use generic ice particle size distribution

! Water tracer structure
type(mp_cpr_old_wtrac_type), intent(in out) :: wtrac_mp_cpr_old

! Local Variables

integer ::                                                                     &
  i,                                                                           &
                        ! Loop counter for points
  cry_offset        ! Index offset for ice crystals

real (kind=real_lsprec) ::                                                     &
  aplusb,                                                                      &
                        ! A+B terms in diffusional growth equation
  tempw_dep(points),                                                           &
                        ! Temporary for available moisture
  tempw_sub(points),                                                           &
                        ! Temporary for saturation deficit
  pr02,                                                                        &
                        ! Temporary in calculation of PSD slopes
  lamr1,                                                                       &
                        ! Power of PSD slope
  lamr2,                                                                       &
                        ! Power of PSD slope
  dqi,                                                                         &
                        ! Temporary in calculating ice transfers
  dqil,                                                                        &
                        ! Temporary in calculating liquid transfers
  dqi_dep(points),                                                             &
                        ! Deposition amount / kg kg-1
  dqi_sub(points),                                                             &
                        ! Sublimation amount / kg kg-1
  cfltemp,                                                                     &
                        ! Temporary in calc of cfl and qcl change
  deltacf,                                                                     &
                        ! Change in cf across timestep
  deltacfl          ! Change in cfl across timestep

real (kind=real_lsprec) ::                                                     &
  m_1(points),                                                                 &
                        ! 1st moment of generic particle size dist.
  m_0p5_dicp3(points),                                                         &
                        ! 1+(dic+1)/2 moment of generic PSD
  m_0p5_dip3(points)
                        ! 1+(di+1)/2 moment of generic PSD

! Ice-mass that is not falling out this microphysics timestep
real (kind=real_lsprec) :: qcf_nofall(points)

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

character(len=*), parameter :: RoutineName='LSP_DEPOSITION'


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

! Set qcf amount that is not falling out
do i = 1, points
  qcf_nofall(i) = qcf(i) * ice_nofall(i)
end do

    ! Use the generic ice particle size distribution
if (l_psd) then

   ! Calculate the 1st (cx(84)) moment of the
   ! ice particle size distribution

  call lsp_moments(points,rho,t,qcf_nofall,cficei,cx(84),m_1)

  if (.not. l_diff_icevt ) then
    ! Use only one vt-D relation
    ! Calculate the 1+0.5(di+1) (cx(85)) moment of the
    ! ice particle size distribution

    call lsp_moments(points,rho,t,qcf_nofall,cficei,cx(85),m_0p5_dip3)
  else
     ! The vt-D relation which gives the
     ! least mass weighted mean fallspeed will be used so
     ! calculate both required moments

    call lsp_moments(points,rho,t,qcf_nofall,cficei,cx(185),m_0p5_dicp3)
             ! ice ventilation moment with crystal parameters

    call lsp_moments(points,rho,t,qcf_nofall,cficei,cx(85),m_0p5_dip3)
             ! ice ventilation moment with aggregate parameters
  end if ! l_diff_icevt


end if


! Identify the points that need to be calculated.
npts = 0
do i = 1, points

  if (qcf_nofall(i) >  m0 .and. t(i) <  zerodegc) then

    npts = npts + 1
    ix(npts) = i

  end if ! qcf_nofall(i) >  m0 etc.

end do


! In all the following loops, the 'i' variable is used for the
! original index while the 'c' variable is used for the compressed
! index.

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)


      !-----------------------------------------------
      ! Diffusional growth parameters
      !-----------------------------------------------
  aplusb = (apb1-apb2*t(i)) * esi(i)
  aplusb = aplusb + (t(i)**3) * p(i) * apb3

      !-----------------------------------------------
      ! Moisture available from subgrid scale calculation
      !-----------------------------------------------
  tempw_dep(i) = qsl(i) * area_mix(i)                                          &
               + min(q_ice_1(i),qsl(i)) * area_ice_1(i)                        &
               - qs(i) * (area_mix(i) + area_ice_1(i))
  tempw_sub(i) = (q_ice_2(i) - qs(i)) * area_ice_2(i)

      !-----------------------------------------------
      ! Calculate transfer rates
      !-----------------------------------------------
  if (l_psd) then
        ! Use generic particle size distribution
    if (.not. l_diff_icevt) then
       ! Only one choice of vt-D parameters for ventilation
       ! constp(83) = 2 pi axial_ratio_correction
       ! constp(84) = ventilation coefficient 1
       ! constp(85) = ventilation coefficient 2
       !        * Sc^(1/3)*ci^0.5/viscosity0^0.5
      dqi = constp(83) * t(i)**2 * esi(i)                                      &
                       * (constp(84)*m_1(i)*corr2(i)                           &
                       +  constp(85)*rocor(i)*m_0p5_dip3(i))                   &
                       / (qs(i) * aplusb * rho(i))
    else
      if (l_use_agg_vt(i)) then
         ! Use aggregate parameters for ventilation
        dqi = constp(83) * t(i)**2 * esi(i)                                    &
                      * (constp(84)*m_1(i)*corr2(i)                            &
                      +  constp(85)*rocor(i)*m_0p5_dip3(i))                    &
                      / (qs(i) * aplusb * rho(i))
      else
         ! Use crystal parameters for ventilation
        dqi = constp(83) * t(i)**2 * esi(i)                                    &
                      * (constp(84)*m_1(i)*corr2(i)                            &
                      +  constp(185)*rocor(i)*m_0p5_dicp3(i))                  &
                      / (qs(i) * aplusb * rho(i))
      end if

    end if ! l_diff_icevt

  else
        ! Use particle size distribution based on intercepts
    pr02=rho(i)*qcf_nofall(i)*cficei(i)*constp(5+cry_offset)*tcgi(i)
    lamr1 = pr02**cx(4+cry_offset)
    lamr2 = pr02**cx(5+cry_offset)
    dqi = tcg(i) * constp(6+cry_offset) * t(i)**2 *                            &
          esi(i) * (constp(7+cry_offset) * corr2(i) *                          &
          lamr1 + constp(8+cry_offset) * rocor(i) *                            &
          lamr2) / (qs(i) * aplusb * rho(i))

  end if  ! l_psd

  dqi_dep(i) = dqi * tempw_dep(i)
  dqi_sub(i) = dqi * tempw_sub(i)

end do

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)


  if (dqi_dep(i) >  zero) then  ! Limits depend on whether
                                ! deposition or sublimation

        !-----------------------------------------------
        ! Deposition is occuring
        !-----------------------------------------------

        !-----------------------------------------------
        ! Molecular diffusion to/from a surface is more efficient
        ! when a particle is at a molecular step. This is more
        ! likely for sublimation. For growth, reduce rate by 10%
        ! ----------------------------------------------
    dqi_dep(i)=0.9_real_lsprec*dqi_dep(i)

        !-----------------------------------------------
        ! Latent heat correction (equivalent to aL in LS cloud)
        !-----------------------------------------------
    tempw_dep(i)=tempw_dep(i)*lheat_correc_ice(i)

        !-----------------------------------------------
        ! Calculate available moisture and transfer amount.
        !-----------------------------------------------
    if (cfliq(i) >  zero) then
        ! Include liquid water contribution. Ignore latent
        ! heat correction from freezing liquid.
      tempw_dep(i)=tempw_dep(i)+qcl(i)*areamix_over_cfliq(i)                   &
                  +max((q_ice_1(i)-qsl(i)),zero)*area_ice_1(i)
    end if
    dqi_dep(i)=min(dqi_dep(i)*timestep,tempw_dep(i))

  end if  ! dqi_dep gt 0.0

end do


!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)

  if (dqi_sub(i)  <   zero) then
        !-----------------------------------------------
        ! Sublimation is occuring
        !-----------------------------------------------
        ! Limits are spare moisture capacity and QCF
        ! outside liquid cloud
    dqi_sub(i) = max( max( dqi_sub(i) * timestep,                              &
                           tempw_sub(i) * lheat_correc_ice(i) ),               &
                 -( qcf(i) * area_ice_2(i) * cficei(i) ))
    ! Due to rounding errors, the limit ( qcf(i) * area_ice_2(i) * cficei(i) )
    ! can occasionally come out slightly bigger than qcf, so that we create
    ! very small negative qcf.  Explicitly avoid this
    ! (only need to fix this here if applying process rates to fluxes,
    !  since in that case the snowfall flux inherits the negative value;
    !  qcf itself has negative values removed in lsp_tidy).
    if ( l_proc_fluxes )  dqi_sub(i) = max( dqi_sub(i), -qcf(i) )

        !-----------------------------------------------
        ! Update cloud fractions
        !-----------------------------------------------
    deltacf = area_ice_2(i) * ( sqrt( max(                                     &
               one+dqi_sub(i)*cfice(i)/(qcft(i)*area_ice_2(i)),                &
                 zero) ) - one )

    ! The above calculation is ill-conditioned when the square
    ! root argument is slightly non-zero, which often occurs
    ! spuriously due to rounding error. The following test uses
    ! tolerances to detect when this situation has occured, and
    ! revises the calculation of deltacf accordingly:
    if (abs( dqi_sub(i)+qcft(i) )     < 1.0e-16_real_lsprec .and.              &
        abs( cfice(i)-area_ice_2(i) ) < 1.0e-12_real_lsprec) then
      deltacf = -area_ice_2(i)
    end if

    cff(i) = cff(i) + deltacf
    cf (i) = cf (i) + deltacf

    cftransfer(i) = cftransfer(i)  + deltacf * one_over_tsi
    cfftransfer(i)= cfftransfer(i) + deltacf * one_over_tsi

  end if  ! dqi_sub lt 0.0

end do

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do c = 1, npts

  i = ix(c)

    !-----------------------------------------------
    ! Calculate liquid water change
    !-----------------------------------------------
  if (cfliq(i) >  zero .and. area_mix(i) >  zero                               &
          .and. qcl(i) >  zero) then

        ! Deposition removes some liquid water content
        ! First estimate of the liquid water removed is explicit

    dqil = max (min ( dqi_dep(i)*area_mix(i)                                   &
           /(area_mix(i)+area_ice_1(i)),                                       &
           qcl(i)*areamix_over_cfliq(i)) ,zero)

      !-----------------------------------------------
      ! Update liquid cloud fraction (and new liquid water est.)
      !-----------------------------------------------
        ! First estimate of the liquid cloud fraction is based
        ! on this explicit estimate of liquid lost by deposition

    cfltemp=cfl(i)*sqrt(max(one-dqil/qcl(i),zero))

        ! Now form a half timestep estimate of the proportion
        ! of the depositing volume which contains liquid cloud

    cfltemp=half*max(area_mix(i)-(cfl(i)-cfltemp),zero)                        &
                       /(area_mix(i)+area_ice_1(i))                            &
               + half*area_mix(i)/ (area_mix(i)+area_ice_1(i))

        ! Recalculate an estimate of the liquid water removed
    dqil=max( min( dqi_dep(i)*cfltemp,qcl(i)*areamix_over_cfliq(i)),zero)

          ! Update liquid cloud fraction and transfer rate
    deltacfl = cfl(i)

    cfl(i) = cfl(i) * sqrt(max(one-dqil/qcl(i), zero))

    cfltransfer(i) = cfltransfer(i) + (cfl(i) - deltacfl) * one_over_tsi

  else
        ! Deposition does not remove any liquid water content
    dqil=zero
  end if ! cfliq gt 0.0 etc.

      !-----------------------------------------------
      ! Adjust ice content
      !-----------------------------------------------
  qcf(i) = qcf(i) + dqi_dep(i) + dqi_sub(i)


      !-----------------------------------------------
      ! Adjust liquid and vapour contents (liquid adjusts first)
      !-----------------------------------------------

  qcl(i) = qcl(i) - dqil  ! Bergeron Findeisen acts first
  t(i) = t(i) + lfrcp * dqil
  dqi = dqi_dep(i) + dqi_sub(i)- dqil

  q(i) = q(i) - dqi
  t(i) = t(i) + lsrcp * dqi

    !-----------------------------------------------
    ! Store depostion/sublimation rate
    !-----------------------------------------------
  ptransfer(i) = ptransfer(i) + (dqi_dep(i) + dqi_sub(i)) * one_over_tsi


    !-----------------------------------------------
    ! Store phase change amounts if using water tracers
    !-----------------------------------------------
  if (l_wtrac) then
    wtrac_mp_cpr_old%qchange(i) = - dqi_sub(i)
    wtrac_mp_cpr_old%depos_l(i) = dqil
    wtrac_mp_cpr_old%depos_v(i) = dqi_dep(i) - dqil
  end if

  if ( l_fix_gr_autoc ) then
    ! Store the tendency due to deposition only for use in graupel
    ! autoconversion (not net of deposition and sublimation, as the latter only
    ! occurs outside the liquid cloud and so shouldn't affect graupel production
    ! inside the mixed-phase cloud).
    ! dqi_dep should always be positive, but when running with single-precision
    ! it occasionally goes very slightly negative due to rounding error issues.
    ! Clip negative values to zero.
    graut_psdep(i) = max(dqi_dep(i), zero) * one_over_tsi
  else
    ! Use the buggy deposition for graupel autoconversion. This will produce
    ! graupel from snow in sublimating ice cloud (which is clearly wrong)
    ! due to the addition of the dqi_sub term in the equation below. However,
    ! it is being retained for now to preserve answers in the lightning scheme
    ! which is significantly impacted by this change.
    graut_psdep(i) = (dqi_dep(i) + dqi_sub(i)) * one_over_tsi
  end if ! l_fix_gr_autoc

end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_deposition
end module lsp_deposition_mod
