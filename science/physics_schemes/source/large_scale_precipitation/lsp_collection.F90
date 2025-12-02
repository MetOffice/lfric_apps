! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Capture of one ice category by
!  another.
module lsp_collection_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_COLLECTION_MOD'

contains

subroutine lsp_collection(                                                     &
  points, timestep,                                                            &
                                          ! Number of points and tstep
  qcf1, qcf2, t,                                                               &
                                          ! Water contents and temp
  area_mix, area_ice,                                                          &
  cficei1, cficei2,                                                            &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
  rho, rhor, m0, tcg1, tcg1i, tcg2, tcg2i,                                     &
                                          ! Parametrization information
  corr, ice1_nofall, ice2_nofall, ice_type1, ice_type2,                        &
                                          ! Microphysical information
  l_psd_2 ,  l_use_area,                                                       &
  l_no_t_check,                                                                &
                                          ! Code options
  ptransfer,                                                                   &
                                          ! Mass transfer diagnostic
  one_over_tsi,                                                                &
                                          ! 1/(timestep*iterations)
  l_use_agg_vt,                                                                &
  dqprec_mix, dqprec_ice, qrain, precfrac_k                                    &
  )

!Use in reals in lsprec precision, both microphysics related and general atmos
use lsprec_mod, only: cx, constp, zerodegc,                                    &
                      zero, one, two, small_number

  ! Microphysics modules- logicals and integers
! use mphys_constants_mod, only: cx, constp, ice_type_offset
use mphys_constants_mod, only: ice_type_offset
use mphys_inputs_mod,    only: l_diff_icevt, l_mcr_precfrac,                   &
                               l_subgrid_graupel_frac, i_update_precfrac,      &
                               i_homog_areas, i_sg_correl

use lsp_combine_precfrac_mod, only: lsp_combine_precfrac

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,             only: real_lsprec

! Dr Hook Modules
use yomhook,             only: lhook, dr_hook
use parkind1,            only: jprb, jpim

! Large scale precipitation modules
use lsp_moments_mod,     only: lsp_moments

implicit none

! Purpose:
!   Update cloud prognostics as a result of collisions between different
!   categories of ice particles

!  Method:
!   Solve the microphysical transfer equation for a specified
!   distribution of ice particles (species 1) sweeping out a different
!   specified distribution of ice particles (species 2) and adding
!   combining the mass into species 1.
!   To clarify, qcf2 in this routine is not nessecesarily the same thing
!   as the UM variable called qcf2!  qcf1 and cficei1 are the ice water
!   content and 1/fraction of any collecting ice species, while
!   qcf2 and cficei2 are the content and 1/fraction of whichever ice
!   species is collected by it.
!   Note that the current formulation is limited to the collection of
!   snow by graupel or the collection of ice by snow. If we want to do
!   graupel to ice then we need to look again at the prefactors
!   constp(16 to 19 +ice_offset1) because they depend on both the
!   species, not just species 1.
!   The commented out code referring to cloud fraction updates is to
!   highlight where in the subroutine changes would need to be made
!   if it was later decided to update cloud fractions in some way.
!   We note that we are only allowing use of the generic ice particle
!   size distribution for second species, intending its use for the
!   collection of aggregates (generic psd) by graupel and *not* the
!   collection of crystals by aggregates (this would go against the
!   idea that the generic psd represents both the crystals and the
!   aggregates). The logic in the call to the routine should prevent
!   any inconsistencies since l_mcr_qcf2 and l_psd should not both be
!   true.
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
  ice_type1,                                                                   &
  ice_type2
                        ! Type of ice (0 - crystals, 1 - aggregates
                        !              3 - graupel)

real (kind=real_lsprec), intent(in) ::                                         &
  timestep,                                                                    &
                        ! Timestep / s
  area_mix(points),                                                            &
                        ! Fraction of gridbox with mixed phase cloud
  area_ice(points),                                                            &
                        ! Fraction of gridbox with ice-only cloud
  cficei1(points),                                                             &
                        ! 1/Fraction of gridbox with 1st ice species
  cficei2(points),                                                             &
                        ! 1/Fraction of gridbox with 2nd ice species
    rho(points),                                                               &
                          ! Air density / kg m-3
    rhor(points),                                                              &
                          ! 1/Air density / m3 kg-1
    m0,                                                                        &
                          ! Seed ice water content / kg kg-1
    tcg1(points),                                                              &
                          ! T dependent function in ice size distribution
    tcg1i(points),                                                             &
                          ! 1/tcg (no units)
    tcg2(points),                                                              &
                          ! T dependent function in ice size distribution
    tcg2i(points),                                                             &
                          ! 1/tcg (no units)
    corr(points),                                                              &
                          ! Fall velocity correction factor (no units)
    one_over_tsi,                                                              &
                          ! 1/(timestep*iterations)
    t(points),                                                                 &
                          ! Temperature (k)
    ice1_nofall(points),                                                       &
                          ! Fraction of species 1's mass that's not falling out
    ice2_nofall(points)
                          ! Fraction of species 2's mass that's not falling out

real (kind=real_lsprec), intent(in out) ::                                     &
  qcf1(points),                                                                &
                        ! Ice water content of 1st species / kg kg-1
  qcf2(points),                                                                &
                        ! Ice water content of 2nd species / kg kg-1
  ptransfer(points)
                        ! Mass rimed in this timestep / kg kg-1

logical, intent(in) ::                                                         &
  l_psd_2,                                                                     &
                        ! Use generic ice particle size distribution
  l_use_area,                                                                  &
                        ! Use ice area amount in calculating
                        ! gridbox mean transfer rate
  l_no_t_check,                                                                &
                        ! Do not check that temperature is less
                        ! than 0 deg C before proceeding
  l_use_agg_vt(points)
                        ! Determines which vt-D parameters to use

! Stuff for updating the prognostic precip fraction;
! increments to rain+graupel in the mixed-phase and ice-only cloud regions
real(kind=real_lsprec), optional, intent(in out) :: dqprec_mix(points)
real(kind=real_lsprec), optional, intent(in out) :: dqprec_ice(points)
real(kind=real_lsprec), optional, intent(in) :: qrain(points)
real(kind=real_lsprec), optional, intent(in out) :: precfrac_k(points)

! Local Variables

integer ::                                                                     &
  i,                                                                           &
                        ! Loop counter for points
  ice_offset1,                                                                 &
                        ! Index offset for first ice species
  ice_offset2
                        ! Index offset for second ice species

real (kind=real_lsprec) ::                                                     &
  dqi(points),                                                                 &
                        ! Transfer of mixing ratio  / kg kg-1
  vi1,                                                                         &
                        ! Average fall speed of first ice
                        ! particle species  / m s-1
  vi2,                                                                         &
                        ! Average fall speed of second ice
                        ! particle species  / m s-1
  fv1,                                                                         &
                        ! Average fall speed difference between
                        ! ice species  / m s-1
  lami1,                                                                       &
                        ! Reciprocal of slope parameter in
                        ! first ice particle size distribution  / m
  lami2,                                                                       &
                        ! Reciprocal of slope parameter in
                        ! second ice particle size distribution  / m
  lamfac1,                                                                     &
                        ! Combination of lamr1 and lamr2
  collision_eff,                                                               &
                        ! Collision efficiency / no units
  m_0(points), m_1(points), m_2(points),                                       &
                        ! zero, 1st and 2nd moments of the ice PSD
  m_bi_di(points),                                                             &
                        ! bi+di moment of the generic ice size distn
  m_bic_dic(points)
                        ! bic+dic moment of the generic ice size distn

! Area fraction where graupel created
real(kind=real_lsprec) :: area_inc(points)
! Precip content (rain+graupel) before the increment is applied
real(kind=real_lsprec) :: qprec(points)
! Fraction of increment occuring in the ice-only vs mixed-phase cloud
real(kind=real_lsprec) :: frac

! Ice species 1 and 2 mass that is not falling out
real (kind=real_lsprec) :: qcf1_nofall(points)
real (kind=real_lsprec) :: qcf2_nofall(points)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_COLLECTION'


    !-----------------------------------------------
    ! Select appropriate ice parametrization (see c_lspsiz)
    !-----------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ice_offset1=ice_type1*ice_type_offset
ice_offset2=ice_type2*ice_type_offset

! Set amount of qcf1 and qcf2 that is not falling out this timestep
do i = 1, points
  qcf1_nofall(i) = qcf1(i) * ice1_nofall(i)
  qcf2_nofall(i) = qcf2(i) * ice2_nofall(i)
end do

    !-----------------------------------------------
    ! Calculate moments of size distribution if appropriate
    !-----------------------------------------------
if (l_psd_2) then
      ! Calculate the 0th, 1st, 2nd and bi+di (cx(82)) moments of the
      ! ice particle size distribution
  call lsp_moments(points,rho,t,qcf2_nofall,cficei2,zero,m_0)
  call lsp_moments(points,rho,t,qcf2_nofall,cficei2,one,m_1)
  call lsp_moments(points,rho,t,qcf2_nofall,cficei2,two,m_2)
  if (.not. l_diff_icevt ) then
    ! Use only one vt-D relation
    call lsp_moments(points,rho,t,qcf2_nofall,cficei2,cx(82),m_bi_di)
  else
    ! The vt-D relation which gives the
    ! least mass weighted mean fallspeed will be used so
    ! calculate both required moments
    call lsp_moments(points,rho,t,qcf2_nofall,cficei2,cx(182),m_bic_dic)
            ! ice mass flux moment with crystal parameters
    call lsp_moments(points,rho,t,qcf2_nofall,cficei2,cx(82),m_bi_di)
            ! ice mass flux moment with aggregate parameters
  end if
end if  ! l_psd_2

do i = 1, points

  if (qcf1_nofall(i)  >   m0 .and. qcf2_nofall(i)  >   m0                      &
      .and. ( t(i) < zerodegc .or. l_no_t_check)                               &
      .and. ( (area_ice(i)+area_mix(i))  >   zero                              &
             .or. (.not. l_use_area) ) ) then

        !-----------------------------------------------
        ! Calculate first ice mass-weighted fallspeed
        !-----------------------------------------------
        ! Use size distribution based on intercepts
    vi1 = constp(4+ice_offset1) * corr(i)                                      &
        * ( rho(i) * qcf1_nofall(i) * cficei1(i)                               &
          * constp(5+ice_offset1) * tcg1i(i) )**cx(3+ice_offset1)

        !-----------------------------------------------
        ! Calculate second ice mass-weighted fallspeed
        !-----------------------------------------------
    if (l_psd_2) then
          ! Use generic ice size distribution
      if (.not. l_diff_icevt) then
         ! Only one set of vt-D parameters
        vi2 = constp(82) * corr(i) * m_bi_di(i)                                &
            / (rho(i) * qcf2_nofall(i) * cficei2(i))
      else
        if (l_use_agg_vt(i)) then
           ! Use aggregate parameters
          vi2 = constp(82) * corr(i) * m_bi_di(i)                              &
              / (rho(i) * qcf2_nofall(i) * cficei2(i))
        else
           ! Use crystal parameters
          vi2 = constp(182) * corr(i) * m_bic_dic(i)                           &
              / (rho(i) * qcf2_nofall(i) * cficei2(i))
        end if
      end if
    else
          ! Use size distribution based on intercepts
      vi2 = constp(4+ice_offset2) * corr(i)                                    &
          * ( rho(i) * qcf2_nofall(i) * cficei2(i)                             &
            * constp(5+ice_offset2) * tcg2i(i) )**cx(3+ice_offset2)
    end if  ! l_psd_2

        !-----------------------------------------------
        ! Estimate the mean absolute differences in velocities
        !-----------------------------------------------
    fv1 = max(abs(vi1-vi2),(vi1+vi2)/8.0_real_lsprec)

        !-----------------------------------------------
        ! Calculate reciprocal of lambda for first ice distribution
        !-----------------------------------------------
    lami1 = ( rho(i) * qcf1_nofall(i) * cficei1(i)                             &
            * constp(5+ice_offset1)*tcg1i(i) )**(-cx(7+ice_offset1))

        !------------------------------------------------
        ! Calculate transfer
        !------------------------------------------------
    collision_eff = 0.02_real_lsprec*exp(0.08_real_lsprec*(t(i)-zerodegc))

    if (l_psd_2) then
          ! Use the generic ice particle size distribution
          ! constp(80)=pi**2/24 x1g ag (80 = 20+ice_offset for graup)
          ! constp(91)=gamma(bi+3+x4i)
          ! constp(92)=2 gamma(bi+2+x4i)
          ! constp(93)=gamma(bi+1+x4i)
      dqi(i) = collision_eff * fv1 *  timestep * rhor(i) *                     &
               tcg1(i) * constp(20+ice_offset1) *                              &
               lami1**(-cx(11+ice_offset1))  *                                 &
            (constp(91) * m_0(i) * lami1**cx(13+ice_offset1) +                 &
             constp(92) * m_1(i) * lami1**cx(14+ice_offset1) +                 &
             constp(93) * m_2(i) * lami1**cx(15+ice_offset1))
    else
        !-----------------------------------------------
        ! Calculate reciprocal of lambda for second ice distribution
        !-----------------------------------------------
      lami2 = (rho(i) * qcf2_nofall(i) * cficei2(i)                            &
        *constp(5+ice_offset2)*tcg2i(i))**(-cx(7+ice_offset2))

          ! Use the distribution defined by intercepts
      lamfac1 =                                                                &
         constp(16+ice_offset1)*(lami2**cx(8+ice_offset2)                      &
                               * lami1**cx(13+ice_offset1)) +                  &
         constp(17+ice_offset1)*(lami2**cx(9+ice_offset2)                      &
                               * lami1**cx(14+ice_offset1)) +                  &
         constp(18+ice_offset1)*(lami2**cx(10+ice_offset2)                     &
                               * lami1**cx(15+ice_offset1))

      dqi(i) = collision_eff * tcg1(i) * tcg2(i)                               &
               * constp(19+ice_offset1)                                        &
               * lami1**(-cx(11+ice_offset1))                                  &
               * lami2**(-cx(11+ice_offset2))                                  &
               * fv1 * lamfac1 * timestep * rhor(i)

    end if  ! l_psd_2

    if (l_use_area) then
          ! The calculations above have used in-cloud quantities to
          ! obtain the transfer rates. Multiply the transfer by the
          ! ice cloud amount to get a gridbox mean.
      dqi(i) = dqi(i) * (area_mix(i) + area_ice(i))
    end if

        !------------------------------------------------
        ! Limit transfer to the mass of species 2 that is available
        !------------------------------------------------
    dqi(i) = min(dqi(i),qcf2(i))

    ptransfer(i) = ptransfer(i) + dqi(i) * one_over_tsi

        !------------------------------------------------
        ! Adjust ice species contents
        !------------------------------------------------
    qcf1(i)  = qcf1(i)   + dqi(i)
    qcf2(i)  = qcf2(i)   - dqi(i)

  else

    dqi(i) = zero

  end if ! qcf1_nofall(i) >  m0 etc.

end do  ! Points


if ( l_mcr_precfrac .and. l_subgrid_graupel_frac .and. ice_type1==3 ) then
  ! If the collecting species is graupel, and using sub-grid fraction
  ! for graupel, need to store the contribution of collection of ice-cloud
  ! by graupel to the total precip increment in the ice-cloud partitions.
  if ( i_update_precfrac == i_homog_areas ) then

    if ( present(dqprec_mix) .and. present(dqprec_ice) ) then
      ! Increment applies in both the area_ice and area_mix area partitions;
      ! divvy it up between the two in proportion to area
      do i = 1, points
        frac = area_ice(i) / max( area_ice(i) + area_mix(i), small_number )
        dqprec_ice(i) = dqprec_ice(i) +      frac  * dqi(i)
        dqprec_mix(i) = dqprec_mix(i) + (one-frac) * dqi(i)
      end do
    end if

  else if ( i_update_precfrac == i_sg_correl ) then

    if ( present(qrain) .and. present(precfrac_k) ) then
      ! Set precip mass before the increment was applied
      ! (ignoring negative values), and area fraction of the increment
      do i = 1, points
        qprec(i) = max( qcf1(i) - dqi(i), zero )  ! qcf1 is qgraup
        dqi(i) = max( qcf1(i), zero ) - qprec(i)
        qprec(i) = qprec(i) + max( qrain(i), zero )
        area_inc(i) = area_ice(i) + area_mix(i)
      end do
      ! Calculate combined fraction of the existing precip and increment
      ! (precfrac_k is updated)
      call lsp_combine_precfrac( points,                                       &
                                 qprec, dqi, precfrac_k, area_inc )
    end if

  end if
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_collection
end module lsp_collection_mod
