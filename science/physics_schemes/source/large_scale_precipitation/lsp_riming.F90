! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme.
!
!  Contains the following subroutines:
!
!  lsp_riming (Riming of ice particles including shape-dependence)
!  lsp_riming_sphere  (Riming of spherical ice particles. Used for Graupel)
!  lsp_riming_graupel (wrapper around lsp_riming_sphere, to setup
!                      graupel-specific input arrays)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation
!
! Subroutine Interface:
module lsp_riming_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_RIMING_MOD'

contains

! Subroutine lsp_riming only called if l_shape_rime is true
subroutine lsp_riming(                                                         &
  points, timestep,                                                            &
                                          ! Number of points and tstep
  qcl, qcf, t,                                                                 &
                                          ! Water contents and temp
  area_liq, area_mix, cfliq, cficei,                                           &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
  rho, m0, tcg, tcgi, corr, ice_nofall,                                        &
                                          ! Parametrization information
  lfrcp , ice_type,                                                            &
                                          ! Microphysical information
  l_psd,                                                                       &
                                          ! Code options
  ptransfer,                                                                   &
                                          ! Mass transfer diagnostic
  one_over_tsi,                                                                &
                                          ! 1/(timestep*iterations)
  l_use_agg_vt,                                                                &
  wtrac_mp_cpr_old                                                             &
                                          ! Store used for water tracers
  )

! Use in reals in lsprec precision, both microphysics related and general atmos
use lsprec_mod, only: cx, constp,qclmin_rime, area_ratio_prefac,               &
                               area_ratio_expn, zerodegc,                      &
                               one

  ! Microphysics modules
use mphys_constants_mod, only: ice_type_offset
use mphys_inputs_mod,    only: l_diff_icevt

use science_fixes_mod,   only: l_fix_riming

use um_types,            only: real_lsprec

! Water tracers
use free_tracers_inputs_mod, only: l_wtrac
use wtrac_mphys_mod,         only: mp_cpr_old_wtrac_type

! Dr Hook modules
use yomhook,             only: lhook, dr_hook
use parkind1,            only: jprb, jpim

! Large scale precipitation modules
use lsp_moments_mod, only: lsp_moments
implicit none

! Purpose:
!   Update cloud prognostics as a result of ice particle riming

! Method:
!   Solve the microphysical transfer equation for a specified
!   distribution of ice particles sweeping out liquid water.

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

! Riming is a source term for ice crystals
! and a sink term for supercooled water.
! Occurs if there are ice crystals and liquid water drops present
! in a mixed phase overlap region. The Numerical solution is implicit
! in liquid water.
! Riming does not, in this formulation, update cloud fractions
! which is why these variables have been commented out.
! This code allows for possible dependence of the the riming rate
! on the cross-sectional area of the ice/snow crystals.
! The parametrization of area in terms of particle size from Heymsfield and
! Miloshevich (J. Atmos. Sci., 60, 936 (2003)) is used.
! Also included is a minimum qcl, above which liquid water
! is available for riming. This follows the findings of
! Marimaya (J. Meteor. Japan, 53(6), 384-392) that
! there is a minimum droplet size for riming to occur.

!
! Subroutine Arguments

integer, intent(in) ::                                                         &
  points,                                                                      &
                        ! Number of points to calculate
  ice_type          ! Type of ice (0 - crystals, 1 - aggregates
                        !              3 - graupel)

real (kind=real_lsprec), intent(in) ::                                         &
  timestep,                                                                    &
                        ! Timestep / s
  area_liq(points),                                                            &
                        ! Fraction of gridbox with liquid-only cloud
  area_mix(points),                                                            &
                        ! Fraction of gridbox with mixed phase cloud
  cfliq(points),                                                               &
                        ! Fraction of gridbox with liquid cloud
  cficei(points),                                                              &
                        ! 1/Fraction of gridbox with ice cloud
    rho(points),                                                               &
                          ! Air density / kg m-3
    m0,                                                                        &
                          ! Seed ice water content / kg kg-1
    tcg(points),                                                               &
                          ! T dependent function in ice size distribution
    tcgi(points),                                                              &
                          ! 1/tcg (no units)
    corr(points),                                                              &
                          ! Fall velocity correction factor (no units)
    ice_nofall(points),                                                        &
                          ! Fraction of qcf that is not falling out
    lfrcp,                                                                     &
                          ! Latent heat of fusion
                          ! / heat capacity of air / K
    one_over_tsi
                          ! 1/(timestep*iterations)

real (kind=real_lsprec), intent(in out) ::                                     &
  qcl(points),                                                                 &
                        ! Liquid water content / kg kg-1
  qcf(points),                                                                 &
                        ! Ice water content    / kg kg-1
  t(points),                                                                   &
                        ! Temperature / K
  ptransfer(points) ! Mass rimed in this timestep / kg kg-1

logical, intent(in) ::                                                         &
  l_psd,                                                                       &
                        ! Use generic ice particle size distribution
  l_use_agg_vt(points)
                        ! Determines fallspeed branch to use

type(mp_cpr_old_wtrac_type), intent(in out) :: wtrac_mp_cpr_old

! Local Variables

integer ::                                                                     &
  i,                                                                           &
                        ! Loop counter for points
  cry_offset        ! Index offset for ice crystals

real (kind=real_lsprec) ::                                                     &
  qclnew,                                                                      &
                        ! For value of qcl after riming  / kg kg-1
  dqi,                                                                         &
                        ! Amount of ice rimed  / kg kg-1
  m_2_di(points),                                                              &
                        ! riming moment of particle size distribution
  m_2_dic(points),                                                             &
                        ! riming moment of particle size distribution
  ic_qcl,                                                                      &
                        ! In-cloud liquid water (qcl)
  qcl_rime,                                                                    &
                        ! LWC available for riming
  rime_eff_agg,                                                                &
                        ! Shape-dependent riming efficency of aggregates
  rime_eff_cry,                                                                &
                        ! Shape-dependent riming efficency of crystals
  rime_expn_agg,                                                               &
                        ! Order of riming moment for aggregate,
                        ! modified to account for cross-sectional area
  rime_expn_cry,                                                               &
                        ! Order of riming moment for crystals,
                        ! modified to account for cross-sectional area
  rime_fac
                        ! Factor in the implicit solution used for
                        ! qcl after riming

logical :: l_riming_conditions_met
real (kind=real_lsprec) :: qcl_min_threshold
  ! Minimum threshold of qcl for riming

! Amount of qcf that is not falling out
real (kind=real_lsprec) :: qcf_nofall(points)

! Threshold cloud-fraction below-which to skip riming
real (kind=real_lsprec), parameter :: cloud_tol_riming = 0.001_real_lsprec

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_RIMING'

    !-----------------------------------------------
    ! Select appropriate ice parametrization (see c_lspsiz)
    !-----------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set minimum threshold for qcl. This needs to be such that the in-cloud
! qcl ( i.e. qcl / cloud fraction) will be greater than qclmin_rime.
! To obtain this, multiply the minimum in-cloud value by the minimum
! acceptable cloud fraction
qcl_min_threshold = cloud_tol_riming * qclmin_rime

! Set other properties
cry_offset        = ice_type * ice_type_offset
rime_expn_agg     = cx(81) + area_ratio_expn
rime_expn_cry     = cx(181) + area_ratio_expn
rime_eff_agg      = area_ratio_prefac
rime_eff_cry      = area_ratio_prefac

! Set amount of qcf that is not falling out this timestep
do i = 1, points
  qcf_nofall(i) = qcf(i) * ice_nofall(i)
end do

if (l_psd) then
      ! Use the generic ice particle size distribution
      ! Calculate the 2+di (cx(81)) moment of the
      ! ice particle size distribution.

  if (.not. l_diff_icevt) then
        ! Only one set of vt-D parameters
    call lsp_moments(points,rho,t,qcf_nofall,cficei,rime_expn_agg,m_2_di)
  else

    ! The vt-D relation which gives the
    ! least mass weighted mean fallspeed will be used so
    ! calculate both required moments
    call lsp_moments(points,rho,t,qcf_nofall,cficei,rime_expn_cry,m_2_dic)
                      ! psd moment with crystal parameters
    call lsp_moments(points,rho,t,qcf_nofall,cficei,rime_expn_agg,m_2_di)
                      ! psd moment with aggregate parameters
  end if
end if

! May be required for 32-bit PROC comparability for some configs
! m_2_di(:)  = max(m_2_di,epsilon(one))
! m_2_dic(:) = max(m_2_dic,epsilon(one))

do i = 1, points

  ! Set riming conditions to False and then check if they are being met
  ! due to a different set of conditions. The first set should not produce
  ! very small values of qcl, the second set should.
  l_riming_conditions_met = .false.

  if ( l_fix_riming .and. qcf_nofall(i) > m0 .and.                             &
       qcl(i) > qcl_min_threshold .and.                                        &
       t(i) < zerodegc .and. area_mix(i) > cloud_tol_riming            .and.   &
       cfliq(i) > cloud_tol_riming                                    ) then
    ! Only rime if the in-cloud qcl is bigger than the threshold set
    ! and cloud fractions are bigger than the minimum tolerance set by the
    ! cloud scheme. By definition qcl(i) must be positive if ic_qcl is
    ! above qcl_rime and the cloud fractions are positive and above the
    ! tolerances.

    l_riming_conditions_met = .true.

  else if ( qcf_nofall(i) > m0 .and. qcl(i) > 0.0 .and. t(i) < zerodegc        &
           .and. area_mix(i) > 0.0 .and. cfliq(i) > 0.0 ) then

    ! l_fix_riming turned off: rime only if original conditions have been met
    ! which may lead to small numbers being passed through the model

    l_riming_conditions_met = .true.

  end if

  if ( l_riming_conditions_met ) then

    ! Determine what the in-cloud qcl value is - used in a few places
    ic_qcl = qcl(i) / cfliq(i)

        !-----------------------------------------------
        ! Calculate water content of mixed phase region
        !-----------------------------------------------
    if (l_psd) then

          ! Calculate the riming rate using the generic PSD
      if (.not. l_diff_icevt) then
         ! Only one set of vt-D parameters
         ! constp(81) = (pi/4) ci
        rime_fac = constp(81)                                                  &
                   *corr(i)*timestep*rime_eff_agg*m_2_di(i)

      else
        if (l_use_agg_vt(i)) then
           ! Use aggregate parameters
          rime_fac = constp(81)                                                &
                   *corr(i)*timestep*rime_eff_agg*m_2_di(i)

        else
           ! Use crystal parameters
          rime_fac = constp(181)                                               &
                   *corr(i)*timestep*rime_eff_cry*m_2_dic(i)

        end if
      end if
    else
          ! Use the defined gamma distribution

      rime_fac = rime_eff_agg*constp(9+cry_offset)*tcg(i)                      &
                *corr(i)*timestep*(rho(i)*qcf_nofall(i)*cficei(i)              &
                *constp(5+cry_offset)*tcgi(i))**cx(6+cry_offset)

    end if

    ! Only the fraction of qcl above qclmin is
    ! available for riming.
    ! Note that if in-liquid-cloud qcl < qclmin_rime
    ! then it is unchanged by riming, i.e., qclnew=qcl/cfliq
    !
    qcl_rime  = ic_qcl + rime_fac * qclmin_rime

    qclnew = qcl_rime / (one + rime_fac)

    qclnew = min(ic_qcl, qclnew)

      !-----------------------------------------------
      ! Convert to new grid box total water content
      !-----------------------------------------------
      ! The mixed-phase cloud now contains the result
      ! of the riming process, which may be the same
      ! qcl as before riming if qclmin is not reached
    qclnew = qcl(i)*area_liq(i)/cfliq(i)+qclnew*area_mix(i)
    dqi = ( qcl(i) - qclnew )


        !-----------------------------------------------
        ! Store process rate / kg kg-1 s-1 and cloud fraction changes
        !-----------------------------------------------
    ptransfer(i) = ptransfer(i) + dqi * one_over_tsi

        !-----------------------------------------------
        ! Update water contents
        !-----------------------------------------------
    qcf(i) = qcf(i) + dqi
    t(i)   = t(i) + lfrcp * dqi
    qcl(i) = qclnew
    if (l_wtrac) wtrac_mp_cpr_old%qchange(i) = dqi

  end if ! qcf_nofall(i) >  m0 etc.

end do  ! Points

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_riming


! Subroutine lsp_riming only called if l_shape_rime is false
subroutine lsp_riming_sphere(                                                  &
  points, timestep,                                                            &
                                          ! Number of points and tstep
  qcl, qcf, t,                                                                 &
                                          ! Water contents and temp
  area_liq, area_mix, cfliq, cficei,                                           &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
  rho, m0, tcg, tcgi, corr, ice_nofall,                                        &
                                          ! Parametrization information
  lfrcp , ice_type,                                                            &
                                          ! Microphysical information
  l_psd,                                                                       &
                                          ! Code options
  ptransfer,                                                                   &
                                          ! Mass transfer diagnostic
  one_over_tsi,                                                                &
                                          ! 1/(timestep*iterations)
  l_use_agg_vt                                                                 &
  )

use lsprec_mod,         only: cx, constp, zerodegc, zero

  ! Microphysics modules
use mphys_constants_mod, only: ice_type_offset
use mphys_inputs_mod,    only: l_diff_icevt

! General atmosphere modules

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,             only: real_lsprec

! Dr Hook modules
use yomhook,             only: lhook, dr_hook
use parkind1,            only: jprb, jpim

! Large scale precipitation modules
use lsp_moments_mod, only: lsp_moments
implicit none

! Purpose:
!   Update cloud prognostics as a result of ice particle riming
!   assuming particles are spherical

! Method:
!   Solve the microphysical transfer equation for a specified
!   distribution of ice particles sweeping out liquid water.

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

! Riming is a source term for ice crystals
! and a sink term for supercooled water.
! Occurs if there are ice crystals and liquid water drops present
! in a mixed phase overlap region. The Numerical solution is implicit
! in liquid water.
! Riming does not, in this formulation, update cloud fractions
! which is why these variables have been commented out.
! In this rountine the ice particles are assumed to be spherical,
! which is appropriate when accumulation of cloud droplets by
! graupel calculated
!

! Subroutine Arguments

integer, intent(in) ::                                                         &
  points,                                                                      &
                        ! Number of points to calculate
  ice_type          ! Type of ice (0 - crystals, 1 - aggregates
                        !              3 - graupel)

real (kind=real_lsprec), intent(in) ::                                         &
  timestep,                                                                    &
                        ! Timestep / s
  area_liq(points),                                                            &
                        ! Fraction of gridbox with liquid-only cloud
  area_mix(points),                                                            &
                        ! Fraction of gridbox with mixed phase cloud
  cfliq(points),                                                               &
                        ! Fraction of gridbox with liquid cloud
  cficei(points),                                                              &
                        ! 1/Fraction of gridbox with ice cloud
    rho(points),                                                               &
                          ! Air density / kg m-3
    m0,                                                                        &
                          ! Seed ice water content / kg kg-1
    tcg(points),                                                               &
                          ! T dependent function in ice size distribution
    tcgi(points),                                                              &
                          ! 1/tcg (no units)
    corr(points),                                                              &
                          ! Fall velocity correction factor (no units)
    ice_nofall(points),                                                        &
                          ! Fraction of ice mass that is not falling out
    lfrcp,                                                                     &
                          ! Latent heat of fusion
                          ! / heat capacity of air / K
    one_over_tsi
                          ! 1/(timestep*iterations)

real (kind=real_lsprec), intent(in out) ::                                     &
  qcl(points),                                                                 &
                        ! Liquid water content / kg kg-1
  qcf(points),                                                                 &
                        ! Ice water content    / kg kg-1
  t(points),                                                                   &
                        ! Temperature / K
  ptransfer(points) ! Mass rimed in this timestep / kg kg-1

logical, intent(in) ::                                                         &
  l_psd,                                                                       &
                        ! Use generic ice particle size distribution
  l_use_agg_vt(points)
                        ! Determines fallspeed branch to use

! Local Variables

integer ::                                                                     &
  i,                                                                           &
                        ! Loop counter for points
  cry_offset        ! Index offset for ice crystals

real (kind=real_lsprec) ::                                                     &
  qclnew,                                                                      &
                        ! For value of qcl after riming  / kg kg-1
  dqi,                                                                         &
                        ! Amount of ice rimed  / kg kg-1
  m_2_di(points),                                                              &
                        ! 2+DI moment of particle size distribution
  m_2_dic(points)
                        ! 2+DIC moment of particle size distribution

! Ice mass that is not falling out
real (kind=real_lsprec) :: qcf_nofall(points)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_RIMING_SPHERE'


    !-----------------------------------------------
    ! Select appropriate ice parametrization (see c_lspsiz)
    !-----------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

cry_offset=ice_type*ice_type_offset

! Set amount of qcf that is not falling out this timestep
do i = 1, points
  qcf_nofall(i) = qcf(i) * ice_nofall(i)
end do

if (l_psd) then
      ! Use the generic ice particle size distribution
      ! Calculate the 2+di (cx(81)) moment of the
      ! ice particle size distribution.

  if (.not. l_diff_icevt) then
        ! Only one set of vt-D parameters

    call lsp_moments(points,rho,t,qcf_nofall,cficei,                           &
                            cx(81),m_2_di)
  else
    ! The vt-D relation which gives the
    ! least mass weighted mean fallspeed will be used so
    ! calculate both required moments
    call lsp_moments(points,rho,t,qcf_nofall,cficei,                           &
                            cx(181),m_2_dic)
                      ! psd moment with crystal parameters
    call lsp_moments(points,rho,t,qcf_nofall,cficei,                           &
                            cx(81),m_2_di)
                      ! psd moment with aggregate parameters
  end if
end if

do i = 1, points

  if (qcf_nofall(i) >  m0 .and. qcl(i) >  zero .and. t(i)  <   zerodegc        &
      .and. area_mix(i) >  zero .and. cfliq(i) >  zero) then

        !-----------------------------------------------
        ! Calculate water content of mixed phase region
        !-----------------------------------------------
    if (l_psd) then

          ! Calculate the riming rate using the generic PSD
      if (.not. l_diff_icevt) then
         ! Only one set of vt-D parameters
         ! constp(81) = (pi/4) ci
        qclnew = qcl(i) /                                                      &
                  (cfliq(i)+cfliq(i)*constp(81)                                &
                   *corr(i)*timestep*m_2_di(i))
      else
        if (l_use_agg_vt(i)) then
           ! Use aggregate parameters
          qclnew = qcl(i) /                                                    &
                  (cfliq(i)+cfliq(i)*constp(81)                                &
                   *corr(i)*timestep*m_2_di(i))
        else
           ! Use crystal parameters
          qclnew = qcl(i) /                                                    &
                  (cfliq(i)+cfliq(i)*constp(181)                               &
                   *corr(i)*timestep*m_2_dic(i))
        end if
      end if
    else
          ! Use the defined gamma distribution

      qclnew = qcl(i) /                                                        &
                (cfliq(i)+cfliq(i)*constp(9+cry_offset)*tcg(i)                 &
                *corr(i)*timestep*(rho(i)*qcf_nofall(i)*cficei(i)              &
                *constp(5+cry_offset)*tcgi(i))**cx(6+cry_offset))
    end if

        !-----------------------------------------------
        ! Convert to new grid box total water content
        !-----------------------------------------------
    qclnew=qcl(i)*area_liq(i)/cfliq(i)+qclnew*area_mix(i)
    dqi=(qcl(i)-qclnew)


        !-----------------------------------------------
        ! Store process rate / kg kg-1 s-1 and cloud fraction changes
        !-----------------------------------------------
    ptransfer(i) = ptransfer(i) + dqi * one_over_tsi

        !-----------------------------------------------
        ! Update water contents
        !-----------------------------------------------
    qcf(i) = qcf(i) + dqi
    t(i)   = t(i) + lfrcp * dqi
    qcl(i) = qclnew

  end if ! qcf_nofall(i) >  m0 etc.

end do  ! Points

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_riming_sphere


subroutine lsp_riming_graupel( points, timestep_mp, one_over_tsi,              &
                               qrain, qcl, qgraup, t, cfliq,                   &
                               rho, tcgg, tcggi, corr, pgacw,                  &
                               rainfraci, rain_liq, rain_mix, graup_nofall,    &
                               dqprec_liq, dqprec_mix, precfrac_k )

use um_types,         only: real_lsprec
use lsprec_mod,       only: zero, one, m0, lfrcp, small_number
use mphys_inputs_mod, only: not_generic_size_dist, i_update_precfrac,          &
                            i_homog_areas, i_sg_correl
use lsp_combine_precfrac_mod, only: lsp_combine_precfrac

! Dr Hook Modules
use yomhook,          only: lhook, dr_hook
use parkind1,         only: jprb, jpim

implicit none

! Number of points
integer, intent(in) :: points
! Microphysics timestep
real (kind=real_lsprec), intent(in) :: timestep_mp
! 1.0/(timestep_mp*iterations)
real (kind=real_lsprec), intent(in) :: one_over_tsi

! Rain water content
real (kind=real_lsprec), intent(in) :: qrain(points)
! Grid-mean liquid cloud mass
real (kind=real_lsprec), intent(in out) :: qcl(points)
! Grid-mean graupel mass
real (kind=real_lsprec), intent(in out) :: qgraup(points)
! Temperature
real (kind=real_lsprec), intent(in out) :: t(points)
! Liquid cloud fraction
real (kind=real_lsprec), intent(in) :: cfliq(points)

! Air density
real (kind=real_lsprec), intent(in) :: rho(points)
! Temperature Factor for graupel size distribution (no units)
real (kind=real_lsprec), intent(in) :: tcgg(points)
! Inverse of TCGC (no units)
real (kind=real_lsprec), intent(in) :: tcggi(points)
! Fall velocity correction factor (no units)
real (kind=real_lsprec), intent(in) :: corr(points)
! Tendency diagnostic
real (kind=real_lsprec), intent(in out) :: pgacw(points)

! 1.0 / total precip fraction (contains both rain and graupel)
real (kind=real_lsprec), intent(in) :: rainfraci(points)
! Overlap of precip fraction with liquid-only and mixed-phase cloud
real (kind=real_lsprec), intent(in) :: rain_liq(points)
real (kind=real_lsprec), intent(in) :: rain_mix(points)

! Fraction of graupel mass that is not falling out
real (kind=real_lsprec), intent(in) :: graup_nofall(points)

! Increments to precipitation (rain+graupel) in liquid, mixed-phase partitions
real (kind=real_lsprec), intent(in out) :: dqprec_liq(points)
real (kind=real_lsprec), intent(in out) :: dqprec_mix(points)

! Prognostic precipitation fraction
real (kind=real_lsprec), intent(in out) :: precfrac_k(points)

! Area fraction in which riming occurs
real (kind=real_lsprec) :: cfliq_rim(points)
! Area of liquid cloud in-which no riming occurs (no overlap with graupel)
real (kind=real_lsprec) :: cfliq_norim(points)
! Graupel mass increment due to riming
real (kind=real_lsprec) :: dqgraup(points)
! Precip mass before increment
real (kind=real_lsprec) :: qprec(points)

! Flag for aggregate vs crystal branch of the fall-speed relation
logical :: l_use_agg_vt(points)

! Fraction of increment applying to the rain_liq area
real (kind=real_lsprec) :: frac

! Loop counter
integer :: i


! Stuff for DrHook
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_RIMING_GRAUPEL'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! Setup various inputs to lsp_riming_sphere...
do i = 1, points

  ! Area in-which riming of graupel occurs = area of overlap of graupel
  ! with both liquid-only and mixed-phase cloud
  cfliq_rim(i) = rain_mix(i) + rain_liq(i)

  ! Area of liquid-cloud in which no riming of graupel occurs
  cfliq_norim(i) = cfliq(i) - cfliq_rim(i)

  ! Assuming graupel should use a crystal-like fall-speed relation not
  ! aggregate??
  l_use_agg_vt(i) = .false.

  ! Store value before for computing increment
  dqgraup(i) = qgraup(i)

end do

! Call riming routine
call lsp_riming_sphere(points, timestep_mp,                                    &
                       qcl, qgraup, t,                                         &
                       cfliq_norim, cfliq_rim, cfliq, rainfraci,               &
                       rho, m0, tcgg, tcggi, corr, graup_nofall, lfrcp, 3,     &
                       not_generic_size_dist,                                  &
                       pgacw, one_over_tsi,                                    &
                       l_use_agg_vt                                            &
                      )

if ( i_update_precfrac == i_homog_areas ) then

  do i = 1, points
    ! Difference qgraup after with qgraup before to compute increment
    dqgraup(i) = qgraup(i) - dqgraup(i)

    ! Increment applies in both the rain_liq and rain_mix area partitions;
    ! divvy it up between the two in proportion to area
    frac = rain_liq(i) / max( rain_liq(i) + rain_mix(i), small_number )
    dqprec_liq(i) = dqprec_liq(i) +      frac  * dqgraup(i)
    dqprec_mix(i) = dqprec_mix(i) + (one-frac) * dqgraup(i)
  end do

else if ( i_update_precfrac == i_sg_correl ) then

  do i = 1, points
    ! Store precip mass before increment was added (ignoring negative values)
    qprec(i) = max( dqgraup(i), zero ) + max( qrain(i), zero )
    ! Difference qgraup after with qgraup before to compute increment
    dqgraup(i) = max( qgraup(i), zero ) - max( dqgraup(i), zero )
  end do

  ! Calculate combined fraction of the existing precip and increment
  ! (precfrac_k is updated)
  call lsp_combine_precfrac( points,                                           &
                             qprec, dqgraup, precfrac_k, cfliq_rim )

end if


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_riming_graupel


end module lsp_riming_mod
