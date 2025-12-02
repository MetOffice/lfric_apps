! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Reflectivity calculation

module mphys_reflec_mod
! Description:
! Calculates radar reflectivity in dBZ for all available
! hydrometer species.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

use um_types, only: real_umphys

implicit none

! Number of reflectivity thresholds used for calculating echo-top height
integer, parameter :: n_dbz_thres = 5

! Reflectivity thresholds used for calculating echo-top height
real(kind=real_umphys)                                                         &
    , parameter  :: dbz_thres(n_dbz_thres) = [0.0, 10.0, 20.0, 30.0, 40.0]

character(len=*), parameter, private :: ModuleName='MPHYS_REFLEC_MOD'

contains

! Subroutine Interface:
subroutine mphys_reflec( points, rho, t, qgraup, qcf, qcf2, qrain, qcl, ndrop, &
                         cfice, cfliq, cfrain, cfgraup, tcg, tcgc,             &
                         dbz_tot, dbz_g, dbz_i, dbz_i2, dbz_l, dbz_r )

!Use in reals in lsprec precision, both microphysics related and general atmos
use lsprec_mod,           only: cx, constp, kliq, kice, mm6m3, ref_lim,        &
                                ref_lim_lin, mr_lim, cf_lim, nd_lim, ref_mom,  &
                                zero, one, ten

!- logicals and integers
use mphys_inputs_mod,     only: l_psd, l_mcr_qgraup, l_no_cf,                  &
                                l_subgrid_graupel_frac
use mphys_bypass_mod,     only: l_crystals

use lsp_moments_mod,      only: lsp_moments

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,             only: real_lsprec

! Dr Hook Modules
use yomhook,              only: lhook, dr_hook
use parkind1,             only: jprb, jpim

implicit none

!------------------------------------------------------------------------------
! Purpose:
!   Calculates any radar reflectivity required for diagnostics
!   Documentation: UMDP 26.

!------------------------------------------------------------------------------
! Subroutine Arguments

integer, intent(in) :: points      ! Number of points calculation is done over

real (kind=real_lsprec), intent(in) :: rho(points)
  ! Air density [kg m-3]
real (kind=real_lsprec), intent(in) :: t(points)
  ! Temperature [K]
real (kind=real_lsprec), intent(in) :: qgraup(points)
  ! Graupel mixing ratio [kg kg-1]
real (kind=real_lsprec), intent(in) :: qcf(points)
  ! Ice Agg mixing ratio [kg kg-1]
real (kind=real_lsprec), intent(in) :: qcf2(points)
  ! Ice Cry mixing ratio [kg kg-1]
real (kind=real_lsprec), intent(in) :: qrain(points)
  ! Rain mixing ratio [kg kg-1]
real (kind=real_lsprec), intent(in) :: qcl(points)
  ! Cloud liquid mixing ratio [kg kg-1]
real (kind=real_lsprec), intent(in) :: ndrop(points)
  ! Cloud droplet number [m-3]
real (kind=real_lsprec), intent(in) :: cfice(points)
  ! Ice cloud fraction
real (kind=real_lsprec), intent(in) :: cfrain(points)
  ! Rain fraction
real (kind=real_lsprec), intent(in) :: cfgraup(points)
  ! Graupel fraction
real (kind=real_lsprec), intent(in) :: cfliq(points)
! Liquid cloud fraction

real (kind=real_lsprec), intent(in) :: tcg(points)
! Agg temperature intercept function
real (kind=real_lsprec), intent(in) :: tcgc(points)
! Cry temperature intercept function

real (kind=real_lsprec), intent(out) :: dbz_tot(points)
  ! Total reflectivity [dBZ]
real (kind=real_lsprec), intent(out) :: dbz_g(points)
  ! Reflectivity due to graupel [dBZ]
real (kind=real_lsprec), intent(out) :: dbz_i(points)
  ! Reflectivity due to ice aggregates [dBZ]
real (kind=real_lsprec), intent(out) :: dbz_r(points)
  ! Reflectivity due to rain [dBZ]
real (kind=real_lsprec), intent(out) :: dbz_i2(points)
  ! Reflectivity due to ice crystals [dBZ]
real (kind=real_lsprec), intent(out) :: dbz_l(points)
  ! Reflectivity due to liquid cloud [dBZ]

!------------------------------------------------------------------------------
! Local Variables

integer :: i

real (kind=real_lsprec) :: one_over_cfice(points)
  ! 1.0 / ice cloud fraction

real (kind=real_lsprec) :: kgraup
  ! Reflectivity prefactor due to graupel
real (kind=real_lsprec) :: krain
  ! Reflectivity prefactor due to rain
real (kind=real_lsprec) :: kclw
  ! Reflectivity prefactor due to cloud liquid water
real (kind=real_lsprec) :: kice_a
  ! Reflectivity prefactor due to ice aggregates
real (kind=real_lsprec) :: kice_c
  ! Reflectivity prefactor due to ice crystals
real (kind=real_lsprec) :: gwc
  ! Graupel water content [kg m-3]
real (kind=real_lsprec) :: iwc
  ! Ice water content: Aggregates [kg m-3]
real (kind=real_lsprec) :: rwc
  ! Rain water content [kg m-3]
real (kind=real_lsprec) :: lwc
  ! Liquid water content [kg m-3]

real (kind=real_lsprec) :: lamr
  ! Lambda (slope parameter) for rain [m-1]
real (kind=real_lsprec) :: lamg
  ! Lambda (slope parameter) for graupel [m-1]
real (kind=real_lsprec) :: lamic
  ! Lambda (slope parameter) for ice cry [m-1]
real (kind=real_lsprec) :: lamia
  ! Lambda (slope parameter) for ice agg [m-1]
real (kind=real_lsprec) :: ze_g(points)
  ! Linear reflectivity due to graupel [mm6 m-3]
real (kind=real_lsprec) :: ze_i(points)
  ! Linear reflectivity due to ice agg [mm6 m-3]
real (kind=real_lsprec) :: ze_i2(points)
  ! Linear reflectivity due to ice cry [mm6 m-3]
real (kind=real_lsprec) :: ze_r(points)
  ! Linear reflectivity due to rain    [mm6 m-3]
real (kind=real_lsprec) :: ze_l(points)
  ! Linear reflectivity due to liq cld [mm6 m-3]
real (kind=real_lsprec) :: ze_tot
  ! Total linear reflectivity [mm6 m-3]

real (kind=real_lsprec) :: mom4(points)
  ! 4th moment of the ice particle size distribution

! Declarations for Dr Hook:
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='MPHYS_REFLEC'

!==============================================================================
! Start of calculations
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise all output dBZ variables. These are initialised to the
! minimum reflectivity as this is usually less than zero.
dbz_tot(:) = ref_lim
dbz_g(:)   = ref_lim
dbz_i(:)   = ref_lim
dbz_r(:)   = ref_lim
dbz_i2(:)  = ref_lim
dbz_l(:)   = ref_lim

! Initialise local variables:
ze_g(:)   = zero
ze_i(:)   = zero
ze_i2(:)  = zero
ze_r(:)   = zero
ze_l(:)   = zero

kgraup = kice / 0.93_real_lsprec
krain  = kliq / 0.93_real_lsprec
kclw   = kliq / 0.93_real_lsprec
kice_a = kice / 0.93_real_lsprec
kice_c = kice / 0.93_real_lsprec

! Generate one_over_cfice - used in a few places.

if (l_no_cf) then

  one_over_cfice(:) = one

else ! l_no_cf

  do i = 1, points
    if (cfice(i) > cf_lim) then
      one_over_cfice(i) = one / cfice(i)
    else
      one_over_cfice(i) = zero
    end if
  end do

end if ! l_no_cf

if (l_psd) then

  !Convert to correct precision

  call lsp_moments(points, rho, t, qcf, one_over_cfice,                        &
                          ref_mom, mom4)

  do i = 1, points

    if ( qcf(i) > mr_lim ) then

      ! This will already have an in-cloud value due to one_over_cfice
      ! being passed into lsp_moments above. However, need to ensure
      ! cloud fraction is dealt with appropriately below.

      if (l_no_cf) then
        ze_i(i) = mm6m3 * cx(90) * mom4(i)
      else
        ze_i(i) = cfice(i) * mm6m3 * cx(90) * mom4(i)
      end if
      ! Here cx(90) = 0.224 * (6.0 * ai/pi/900)**2

    end if ! qcf(i) > mr_lim

  end do ! points

else ! not l_psd

  ! Two ice categories (qcf and qcf2) need to be treated separately

  do i = 1, points

    ! Non-psd Crystals.
    if (l_crystals .and. qcf2(i) > mr_lim .and.                                &
        one_over_cfice(i) /= zero ) then

      ! Work out ice crystal water content
      if (l_no_cf) then
        iwc = qcf2(i) * rho(i)
      else
        iwc = qcf2(i) * rho(i) * one_over_cfice(i)
      end if ! l_no_cf

      ! Determine the slope parameter, Lambda
      ! Note that constp(101) does not include the temperature
      ! intercept from lsp_init. Therefore this must be added
      ! at this stage - see UMDP26 section on inferring PSD relations
      ! from mixing ratios and fluxes
      ! N.B. cx(94) =  1.0 / (bic + 1.0 + x4ic - x2ic)

      lamic = ( ( tcgc(i) * constp(101) ) / iwc ) ** cx(94)

      ! Calculate linear radar reflectivity
      ! N.B. cx(97)      = -1.0 * (1.0 + x4ic + (2.0*bic) - x2ic)
      !      constp(104) = x1ic * (aic**2) * gamma(1.0 + x4ic + ( 2.0 * bic))
      !                         * ( 6.0 /(pi * rho_i2))**2
      ! rho_i2 = 900.0 kg m-3

      if (l_no_cf) then
        ze_i(i) = mm6m3 * kice_c * constp(104) * lamic**cx(97)
      else
        ze_i(i) = mm6m3 * kice_c * cfice(i) * constp(104) * lamic**cx(97)
      end if

    end if ! l_crystals

    ! Non-psd Aggregates.
    if (qcf(i) > mr_lim .and. one_over_cfice(i) /= zero ) then

      ! Work out ice crystal water content
      if (l_no_cf) then
        iwc = qcf(i) * rho(i)
      else
        iwc = qcf(i) * rho(i) * one_over_cfice(i)
      end if ! l_no_cf

      ! Determine the slope parameter, Lambda
      ! Note that constp(102) does not include the temperature
      ! intercept from lsp_init. Therefore this must be added
      ! at this stage - see UMDP26 section on inferring PSD relations
      ! from mixing ratios and fluxes
      ! N.B. cx(95) =  1.0 / (bi + 1.0 + x4i - x2i)

      lamia = ( (tcg(i) * constp(102)) / iwc ) ** cx(95)

      ! Calculate linear radar reflectivity
      ! N.B. cx(98)      = -1.0 * (1.0 + x4i + (2.0*bi) - x2i)
      !      constp(105) = x1i  * (ai**2) * gamma(1.0 + x4i + ( 2.0 * bi))
      !                         * ( 6.0 /(pi * rho_i))**2
      ! rho_i = 900.0 kg m-3

      if (l_no_cf) then
        ze_i(i)  = mm6m3 * kice_a * constp(105) * lamia**cx(98)
      else
        ze_i(i)  = mm6m3 * kice_a * cfice(i) * constp(105) * lamia**cx(98)
      end if ! l_no_cf

    end if ! qcf(i) > mr_lim

  end do ! points

end if ! l_psd

do i = 1, points

  ! Graupel

  if (l_mcr_qgraup .and. qgraup(i) > mr_lim ) then

    ! Work out graupel water content
    if ( l_subgrid_graupel_frac .and. (.not. l_no_cf) ) then
      gwc  = qgraup(i) * rho(i)/cfgraup(i)
    else
      gwc  = qgraup(i) * rho(i)
    end if

    ! Determine the slope parameter, Lambda
    ! (  pi/6 rhog * x1r * gamma(4 + x4g) ) ** (1.0 / (bg  + 1.0 + x4g  - x2g))
    ! (---------------------------------- )
    ! (               GWC                 )

    lamg = (constp(100) / gwc ) ** cx(93)
    ! cx(93) = (1.0 / (bg  + 1.0 + x4g  - x2g))

    ! Calculate linear radar reflectivity (no cloud fraction assumed).
    if ( l_subgrid_graupel_frac .and. (.not. l_no_cf) ) then
      ze_g(i) = cfgraup(i) * mm6m3 * kgraup * constp(103) * lamg**cx(96)
    else
      ze_g(i) = mm6m3 * kgraup * constp(103) * lamg**cx(96)
    end if
    ! constp(103) = x1g  * (ag**2)  * gref1x4g  * ( 6.0 /(pi * rho_g))**2
    ! cx(96)      = -1.0 * (1.0 + x4g  + (2.0*bg)  - x2g)

  end if ! l_mcr_qgraup

  ! Rain
  if (qrain(i) > mr_lim .and. cfrain(i) > cf_lim ) then

    ! Work out rain water content
    if (l_no_cf) then
      rwc = qrain(i) * rho(i)
    else
      rwc = qrain(i) * rho(i)/cfrain(i)
    end if

    ! define a lamr, equivalent to the following equation

    ! (  pi/6 rhow * x1r * gamma(4 + x4r) ) ** (1.0 / (4 + x4r -x2r))
    ! ( ----------------------------------)
    ! (              RWC                  )

    lamr = (constp(106) / rwc ) ** cx(52)
    ! cx(52) = (1.0 / (4 + x4r -x2r))

    ! Then define the reflectivity in mm6 m-3.
    ! Reminder:
    ! constp(107) = gamma(1.0 + (2.0 * 3.0) +x4r) * x1r

    if (l_no_cf) then
      ze_r(i) = mm6m3 * krain * constp(107) * lamr**cx(92)
    else
      ze_r(i) = cfrain(i) * mm6m3 * krain * constp(107) * lamr**cx(92)
    end if

  end if !qrain > mr_lim

  ! Liquid Cloud
  if (qcl(i) > mr_lim .and. ndrop(i) > nd_lim .and. cfliq(i) > cf_lim) then

    ! Equation A12 of Stein et al (2014), Monthly Weather Review.
    ! The value of 201.6 has since been found to be an error and
    ! the correct value of this constant is 5.6, which is used here.
    if (l_no_cf) then
      lwc = qcl(i) * rho(i)
      ze_l(i) = mm6m3 * kclw * (5.6_real_lsprec/ndrop(i)) * (lwc*cx(91))**2
    else
      lwc = qcl(i) * rho(i)/cfliq(i)
      ze_l(i) = cfliq(i) * mm6m3 * kclw * (5.6_real_lsprec/ndrop(i))           &
                * (lwc*cx(91))**2
    end if

  end if ! qcl > mr_lim

  ! Compute total reflectivity
  ze_tot = ze_g(i) + ze_i(i) + ze_i2(i) + ze_r(i) + ze_l(i)

  ! Convert from linear (mm^6 m^-3) to dBZ
  if (ze_tot > ref_lim_lin) then
    dbz_tot(i) = max(ref_lim, ten * log10(ze_tot))
  end if

  if (ze_g(i) > ref_lim_lin) then
    dbz_g(i) = max(ref_lim, ten * log10(ze_g(i)))
  end if

  if (ze_i(i) > ref_lim_lin) then
    dbz_i(i) = max(ref_lim, ten * log10(ze_i(i)))
  end if

  if (ze_i2(i) > ref_lim_lin) then
    dbz_i2(i) = max(ref_lim, ten * log10(ze_i2(i)))
  end if

  if (ze_r(i) > ref_lim_lin) then
    dbz_r(i) = max(ref_lim, ten * log10(ze_r(i)))
  end if

  if (ze_l(i) > ref_lim_lin) then
    dbz_l(i)   = max(ref_lim, ten * log10(ze_l(i)))
  end if

end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine mphys_reflec
!==============================================================================
subroutine echo_top_height(dbz_tot, deltaz, dbz_thres, top)

! Grid bounds module
use atm_fields_bounds_mod, only: tdims

use missing_data_mod,      only: rmdi

implicit none

!------------------------------------------------------------------------------
! Purpose:
!   Calculates the echo-top height of a given reflectivity threshold dbz_thres
!
!------------------------------------------------------------------------------
! Subroutine Arguments
! Reflectivity
real(kind=real_umphys), intent(in) ::                                          &
                    dbz_tot( tdims%i_start : tdims%i_end,                      &
                             tdims%j_start : tdims%j_end,                      &
                                         1 : tdims%k_end )

! Vertical grid spacing
real(kind=real_umphys), intent(in) ::                                          &
                    deltaz(  tdims%i_start : tdims%i_end,                      &
                             tdims%j_start : tdims%j_end,                      &
                                         1 : tdims%k_end )

real(kind=real_umphys), intent(in) :: dbz_thres

! Height of top level with reflectivity 0dbz in cloud layer (m agl)
real(kind=real_umphys), intent(in out) ::                                      &
                        top( tdims%i_start : tdims%i_end,                      &
                             tdims%j_start : tdims%j_end )

!------------------------------------------------------------------------------
! Local Variables

! Loop counters
integer :: i, j, k

! Height of lowest level with reflectivity 0dbz in cloud layer (m agl)
real(kind=real_umphys) ::                                                      &
        base(  tdims%i_start : tdims%i_end,                                    &
               tdims%j_start : tdims%j_end )

! Height above ground level (metres)
real(kind=real_umphys) ::                                                      &
        z_agl( tdims%i_start : tdims%i_end,                                    &
               tdims%j_start : tdims%j_end )

! Minimum layer thickness for finding echo-top height
real(kind=real_umphys), parameter :: min_reflec_layer_depth = 500.0

! Determines whether echo-top height/base needs definining or not
logical :: l_base_found
logical :: l_top_found

! For each column work out the lowest and highest level at which
! reflectivity equals dbz_thres
! Ignore layers which are less than 500m thick so that (in theory)
! only convective clouds are considered. This should hopefully
! reject fog layers and thin cumulus layers if they exist.

!$OMP PARALLEL do                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(none)                                                            &
!$OMP private(k,j,i,l_base_found,l_top_found)                                  &
!$OMP SHARED(tdims,base,z_agl,dbz_tot,deltaz,dbz_thres,top)
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end

    ! For each new column, define that the base and top have
    ! not been found
    l_base_found = .false.
    l_top_found  = .false.

    base(i,j) = rmdi
    top(i,j) = rmdi

    ! Reset z_agl to zero to avoid contamination with any previous
    ! values it may have
    z_agl(i,j) = 0.0

    do k = 1, tdims%k_end

      ! Increment the height above ground level (z_agl) with the
      ! diagnosed microphysics layer thickness.
      ! This will give z_agl as the top of the model level
      z_agl(i,j) = z_agl(i,j) + deltaz(i,j,k)

      ! Check whether the base has already been found then find
      ! what level dbz_tot is dbz_thres

      if ( .not. l_base_found .and. dbz_tot(i,j,k) >= 0.0 ) then
        l_base_found = .true.
        base(i,j) = z_agl(i,j)
      end if

      ! If the base has already been found check where dbz_tot
      ! becomes dbz_thres again to find the top.

      if ( .not. l_top_found .and. l_base_found .and.                          &
           dbz_tot(i,j,k) < dbz_thres ) then
        l_top_found = .true.
        top(i,j) = z_agl(i,j)
      end if

      ! Check depth of layer and reject the layer if its
      ! thickness is less than 500m.

      if ( l_top_found .and. l_base_found .and. (top(i,j) -                    &
           base(i,j)) < min_reflec_layer_depth ) then
        base(i,j) = rmdi
        top(i,j)  = rmdi
        l_top_found = .false.
        l_base_found = .false.
      end if

    end do ! k

  end do ! i
end do ! j
!$OMP end PARALLEL do

end subroutine echo_top_height
!==============================================================================
end module mphys_reflec_mod
