! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Autoconversion of liquid to rain
! Subroutine Interface:
module lsp_autoc_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_AUTOC_MOD'

contains

subroutine lsp_autoc(                                                          &
  points, timestep,                                                            &
                                          ! Number of points and tstep
  qgraup, qcl, qrain, t, p,                                                    &
                                          ! Water contents, temp and p
  cfliq, rhcpt,                                                                &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
  area_liq, area_mix, area_ice,                                                &
                                          ! Cloud fraction overlaps
  rainfrac, rain_liq, rain_mix,                                                &
                                          ! Rain fractions for updating
  rain_ice, rain_clear, rain_new,                                              &
  rho, rhor, corr2,                                                            &
                                          ! Parametrization information
  lcrcp,                                                                       &
                                          ! Microphysical information
  one_over_tsi,                                                                &
                                          ! Number of iterations and
                                          ! 1/(timestep*iterations)
  ptransfer, rftransfer,                                                       &
                                          ! Mass and rainfrac transfer
  n_drop_tpr, n_drop_out,                                                      &
                                          ! Droplet concs.
  r_theta_levels_c, fv_cos_theta_latitude_c,                                   &
  f_arr1, f_arr2, f_arr3,                                                      &
  dqprec_liq, dqprec_mix, dqprec_new,                                          &
                                          ! Total precip incs within liq-cloud
  precfrac_k,                                                                  &
                                          ! Prognostic precip fraction
  wtrac_mp_cpr_old                                                             &
                                          ! Store used for water tracers
  )

! Use in reals in lsprec precision, both microphysics related and general atmos
use lsprec_mod,           only: r_auto, n_auto, power_droplet_auto,            &
                                power_qcl_auto, power_rho_auto,                &
                                consts_auto, aut_pref, aut_qc, aut_nc,         &
                                ec_auto, qcfmin, r, repsilon, lc, pi,          &
                                f_cons, fsd_eff_lam, fsd_eff_phi,              &
                                rad_mcica_sigma, two_d_fsd_factor,             &
                                zero, half, one, two, small_number

  ! Microphysics modules- logicals and integers
use mphys_inputs_mod,     only: l_warm_new, l_fsd_generator, l_auto_debias,    &
                                l_mcr_precfrac, l_subgrid_graupel_frac,        &
                                i_update_precfrac, i_homog_areas, i_sg_correl
use mphys_constants_mod,  only: l_inhomog

! General/atmosphere modules- logicals and integers
use gen_phys_inputs_mod,  only: l_mr_physics

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,             only: real_lsprec

use qsat_mod,             only: qsat_wat, qsat_wat_mix

! Water tracers
use free_tracers_inputs_mod, only: l_wtrac
use wtrac_mphys_mod,         only: mp_cpr_old_wtrac_type

! Dr Hook modules
use yomhook,              only: lhook, dr_hook
use parkind1,             only: jprb, jpim

! FSD parameters module- logicals and integers
use fsd_parameters_mod,   only: ip_fsd_constant

! Constant FSD value- logicals and integers
use rad_input_mod,        only: i_fsd

use lsp_combine_precfrac_mod, only: lsp_combine_precfrac

implicit none

! Purpose:
!   Update cloud and rain due to autoconversion of liquid to rain

! Method:
!   Use a rate based on a power law of liquid content and droplet
!   number with a minimum liquid content for autoconversion.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

!   There are many routes through this code depending upon the
!   selected options. Vector options have been removed since
!   they are very difficult to trace through.
!   1) Calculate number of droplets based upon either a
!      prescribed number or the amount of aerosol (murk or
!      sulphate/sea-salt/biomass/fossil-fuel organic carbon/
!      nitrate) in the model.
!   2) Calculate the autoconversion rate based upon a
!      Tripoli and Cotton power law formulation.
!   3) Calculate the autoconversion limit based either on
!      the concentration of droplets above 20 um diameter or
!      on a Tripoli and Cotton type argument

! Subroutine Arguments


integer, intent(in) ::                                                         &
  points
                        ! Number of points to calculate

real (kind=real_lsprec), intent(in) ::                                         &
  timestep,                                                                    &
                        ! Timestep / s
  qgraup(points),                                                              &
                        ! Graupel content / kg kg-1
  p(points),                                                                   &
                        ! Air pressure / N m-2
  t(points),                                                                   &
                        ! Temperature / K
  cfliq(points),                                                               &
                        ! Fraction of gridbox with liquid cloud
    area_liq(points),                                                          &
                          ! Fraction of gridbox with liquid but no ice
    area_mix(points),                                                          &
                          ! Fraction of gridbox with liquid and ice
    area_ice(points),                                                          &
                          ! Fraction of gridbox with ice but no liquid
    rhcpt(points),                                                             &
                          ! Rhcrit on each point (no units)
    rho(points),                                                               &
                          ! Air density / kg m-3
    rhor(points),                                                              &
                          ! 1 / air density / m3 kg-1
    corr2(points),                                                             &
                          ! Temperature correction factor (no units)
    lcrcp,                                                                     &
                          ! Latent heat of condensation/cP / K
    n_drop_tpr(points),                                                        &
                          ! Droplet number determined by
                          ! droplet taper curves
    r_theta_levels_c(points),                                                  &
                          ! Distance from centre of Earth and ...
    fv_cos_theta_latitude_c(points),                                           &
                          ! ... grid info for working out gridbox size.
    one_over_tsi,                                                              &
                          ! 1/(timestep*iterations)
    f_arr1(points),                                                            &
    f_arr2(points),                                                            &
    f_arr3(points)

real (kind=real_lsprec), intent(in out) ::                                     &
  qcl(points),                                                                 &
                        ! Liquid water content / kg kg-1
  qrain(points),                                                               &
                        ! Rain water content / kg kg-1
  ptransfer(points),                                                           &
                        ! Autoconversion rate / kg kg-1 s-1
    rftransfer(points),                                                        &
                          ! Rate of change of rain fraction / s-1
    rainfrac(points),                                                          &
                          ! Rain fraction (no units)
    rain_liq(points),                                                          &
                          ! Overlap of rain with liquid (no units)
    rain_mix(points),                                                          &
                          ! Overlap of rain with mixed phase region
    rain_ice(points),                                                          &
                          ! Overlap of rain with ice
    rain_clear(points),                                                        &
                          ! Overlap of rain with clear sky
    rain_new(points),                                                          &
                          ! Area where rain produced in air that had no rain
    n_drop_out(points)
                          ! Droplet number from autoconversion

! Stuff for updating the prognostic precip fraction;
! increments to rain+graupel in the mixed-phase and liquid-only cloud regions
real(kind=real_lsprec), intent(in out) :: dqprec_liq(points)
real(kind=real_lsprec), intent(in out) :: dqprec_mix(points)
! Increment in air which didn't have any rain already
real(kind=real_lsprec), intent(in out) :: dqprec_new(points)
! Prognostic precip fraction
real(kind=real_lsprec), intent(in out) :: precfrac_k(points)

type(mp_cpr_old_wtrac_type), intent(in out) :: wtrac_mp_cpr_old
                        ! Water tracers store

!-----------------------------
! Local Variables

integer ::                                                                     &
  i
                        ! Loop counter

real (kind=real_lsprec) ::                                                     &
  n_drop(points),                                                              &
                        ! Number of droplets / m-3
  dpr(points),                                                                 &
                        ! Amount of mass autoconverted / kg kg-1
  qsl_tl(points),                                                              &
                        ! Saturated humidity wrt liquid at
                        !  temperature TL / kg kg-1
  t_l(points),                                                                 &
                        ! Liquid temperature / K
  r_mean,                                                                      &
                        ! Mean radius of droplets / m
  r_mean0,                                                                     &
                        ! Factor in calculation of r_mean
  a_factor,                                                                    &
                        ! 3 r_mean / m
  b_factor,                                                                    &
                        ! 0.5 n_drop / m-3
  n_gt_20,                                                                     &
                        ! Concentration of droplets with
                        !  diameter gt 20 um / m-3
  autorate(points),                                                            &
                        ! Rate constants for autoc.
  autolim(points),                                                             &
                        ! Autoconversion threshold / kg kg-1
  qc,                                                                          &
                        ! Maximum amout of liquid to remove / kg kg-1
    rainfracnew,                                                               &
                          ! Updated rain fraction
    rho_1
                          ! rho^(power_rho_auto - power_qcl_auto + 1)

real (kind=real_lsprec) ::                                                     &
                        ! For debiasing code
  alpha_l,                                                                     &
                        ! dqsat/dT at T_L / kg kg-1 K-1
  a_l,                                                                         &
                        ! 1 / (1 + L/cp alpha)
  sigma_s,                                                                     &
                        ! Width of moisture PDF
  g_l,                                                                         &
                        ! Factor in debiasing calculation
  gacb,                                                                        &
                        ! Factor in debiasing calculation
  ac_factor,                                                                   &
                        ! Multiplying factor for debiasing
  ltiny,                                                                       &
                        ! Largest real number represented on platform
  min_float
                        ! Smallest allowed real value, for safety-check

real (kind=real_lsprec) ::                                                     &
  fsd(points),                                                                 &
          ! fractional standard dev of in-cloud liquid water content
  bias(points),                                                                &
          ! autoconversion rate bias
  x_in_km ! grid-box size in km

! Area fraction of rain produced by autoconversion where no pre-existing rain
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
                        ! Number of points to compute
  kk,                                                                          &
                        ! Compressed point index
  ix(points)
                        ! Original index

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_AUTOC'

    !-----------------------------------------------
    ! Part 1. Calculate droplet number concentration
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


ltiny = log(tiny(one))      ! Initialize ltiny
min_float = tiny(one)

npts = 0 ! Reset index counter for condensed points

! Use pre-determined droplet number concentration from
! lsp_taper_ndrop.F90 subroutine.

do i = 1, points

  if (qcl(i) > zero  .and. cfliq(i) > zero) then

    npts = npts + 1

    ix(npts) = i
    n_drop(npts)  = n_drop_tpr(i)
    n_drop_out(i) = n_drop(npts)

  else

    n_drop_out(i) = zero
    dpr(i) = zero

  end if

end do

! In all the following loops, the 'i' variable is used for the
! original index while the 'kk' variable is used for the compressed
! index.

    !-----------------------------------------------
    ! Part 2: Calculating the autoconversion rate and limit
    !-----------------------------------------------
if (l_warm_new) then

   !------------------------------------------------
   ! Use new autoconversion scheme
   !------------------------------------------------
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do kk = 1, npts

    i = ix(kk)

        ! calculate the in-cloud autoconversion rate
    autorate(i) = aut_pref *                                                   &
                  ((qcl(i)/cfliq(i)) ** aut_qc) *                              &
                  ((1.0e-6_real_lsprec*n_drop(kk)) ** aut_nc)

  end do

  if (l_inhomog) then
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do kk = 1, npts

      i = ix(kk)

      ! calculate bias factor E based on Boutle et al 2012 QJ
      x_in_km = 0.001_real_lsprec*sqrt (r_theta_levels_c(i) * fsd_eff_lam      &
                             * r_theta_levels_c(i) * fsd_eff_phi               &
                             * fv_cos_theta_latitude_c(i)     )

      if (l_fsd_generator) then ! Use same FSD param as cld generator
        if (i_fsd == ip_fsd_constant) then
          fsd(i) = rad_mcica_sigma
        else
          if ( cfliq(i) < one ) then
            fsd(i) = (f_arr2(i)-f_arr3(i)*cfliq(i))                            &
              *(((x_in_km*cfliq(i))**0.333_real_lsprec)                        &
              *((f_cons(1)*x_in_km*cfliq(i))**f_cons(2)+one)                   &
             **(f_cons(3)))
          else
            fsd(i) = f_arr1(i) * ((x_in_km**0.333_real_lsprec)                 &
              *((f_cons(1)*x_in_km)**f_cons(2)+one)**(f_cons(3)))
          end if
        end if
      else ! Use FSD param from Boutle et al 2012 QJ
        if ( cfliq(i) < one ) then
          fsd(i) = (0.45_real_lsprec-0.25_real_lsprec*cfliq(i))*               &
                (((x_in_km*cfliq(i))**0.333_real_lsprec)                       &
                *((0.06_real_lsprec*x_in_km*cfliq(i))**1.5_real_lsprec+one)    &
                **(-0.17_real_lsprec))
        else
          fsd(i) = 0.11_real_lsprec*(((x_in_km*cfliq(i))**0.333_real_lsprec)   &
                *((0.06_real_lsprec*x_in_km*cfliq(i))**1.5_real_lsprec+one)    &
                **(-0.17_real_lsprec))
        end if
      end if
      fsd(i)=fsd(i)*two_d_fsd_factor

    end do

  end if

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do kk = 1, npts

    i = ix(kk)

    if (l_inhomog) then

      bias(i) = ((one+fsd(i)**2)**(-half*aut_qc))*                             &
             ((one+fsd(i)**2)**(half*aut_qc**2))

    else ! no inhomog param

      bias(i) = one

    end if !l_inhomog

    dpr(i) = autorate(i) * bias(i) * timestep * cfliq(i)
    dpr(i) = max(min(dpr(i),qcl(i)-qcfmin),zero)

      !-----------------------------------------------
      ! Update liquid water content and rain
      !-----------------------------------------------
    qcl(i)   = qcl(i)   - dpr(i)
    qrain(i) = qrain(i) + dpr(i)

      !-----------------------------------------------
      ! Store process rate
      !-----------------------------------------------
    ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

    if (dpr(i) >  zero .and. (.not. l_mcr_precfrac) ) then
      ! (not doing this here if using prognostic precip fraction,
      !  as in that case rainfrac is updated later)
        !-----------------------------------------------
        ! Calculate change in rain fraction
        !-----------------------------------------------
      rainfracnew = max(rainfrac(i),cfliq(i))
      rftransfer(i) = rftransfer(i)                                            &
                      + (rainfracnew - rainfrac(i)) * one_over_tsi

        !-----------------------------------------------
        ! Update rain fractions
        !-----------------------------------------------
      rainfrac(i)   = rainfracnew
      rain_liq(i)   = min(area_liq(i),rainfrac(i))
      rain_mix(i)   = min(area_mix(i),rainfrac(i)-rain_liq(i))
      rain_ice(i)   =                                                          &
          min(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i))
      rain_clear(i) =                                                          &
          rainfrac(i)-rain_liq(i)-rain_mix(i)-rain_ice(i)

    end if  ! dpr(i) gt 0

  end do

else ! original autoconversion etc

      !-----------------------------------------------
      ! Use the later 3C/3D method for autoconversion
      !-----------------------------------------------

    ! Calculate multiplying factor for droplet size
  r_mean0=(27.0_real_lsprec/(80.0*pi*1000.0_real_lsprec))                      &
          **(-one/3.0_real_lsprec)

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do kk = 1, npts

    i = ix(kk)

      ! Calculate inverse of mean droplet size
    r_mean = r_mean0                                                           &
           * max( rho(i) * qcl(i) / (cfliq(i) * n_drop(kk)),                   &
                  min_float ) ** (-one/3.0_real_lsprec)
      ! Safety-check avoids getting NaN due to floating-point underflow
      ! in the term     rho qcl / cfliq n_drop    when qcl is extremely small
      ! (raising 0 to the power -1/3 is a div-by-zero).
      ! Therefore that term is limited to not fall below tiny.

      ! Calculate numerical factors
    b_factor = 3.0_real_lsprec * r_mean
    a_factor = n_drop(kk) * half

      !-----------------------------------------------
      ! Calculate droplet number concentration greater
      ! than threshold radius
      !-----------------------------------------------

    ! Only compute Exponent if it will not generate an inexact signal.
    ! Required for initial ENDGame runs on IBM.

    if ( (b_factor * r_auto) < -ltiny) then

      n_gt_20 = (b_factor**2 * a_factor) * r_auto**2                           &
                                * exp(-b_factor * r_auto)                      &
                 + (two * b_factor * a_factor) * r_auto                        &
                                * exp(-b_factor * r_auto)                      &
                 + (two * a_factor) * exp(-b_factor * r_auto)

    else

      n_gt_20 = zero

    end if

      !-----------------------------------------------
      ! Test to see if there is a sufficient concentration of droplets
      ! with diameters > threshold for autoconversion to proceed.
      !-----------------------------------------------
    if (n_gt_20  >=  n_auto) then
        ! Calculate autoconversion rate
      autorate(i) = consts_auto * ec_auto                                      &
                * n_drop(kk) ** power_droplet_auto
    else
        ! No autoconversion
      autorate(i)=zero
    end if ! n_gt_20 >= n_auto

  end do

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do kk = 1, npts

    i = ix(kk)

      !-----------------------------------------------
      ! Calculate the autoconversion limit
      !-----------------------------------------------
      ! Calculate value of local qcl at which the droplet concentration
      ! with radii greater than 20um will fall below a threshold (1000 m-3).
      ! This is a hardwired numerical approximation
      !-----------------------------------------------
    autolim(i)=(6.20e-31_real_lsprec*n_drop(kk)**3)-                           &
                (5.53e-22_real_lsprec*n_drop(kk)**2)                           &
           +(4.54e-13_real_lsprec*n_drop(kk))+(3.71e-6_real_lsprec)            &
           -(7.59_real_lsprec/n_drop(kk))

  end do

    !-----------------------------------------------
    ! Part 3: Optionally debias the autoconversion rate
    !-----------------------------------------------
  if (l_auto_debias) then

      !-----------------------------------------------
      ! Calculate qsat(TL) in order to allow subgrid
      ! debiasing of the autoconversion rate.
      !-----------------------------------------------
    do i = 1, points
      t_l(i) = t(i) - (lcrcp * qcl(i) )
    end do

    if (l_mr_physics) then
      call qsat_wat_mix(qsl_tl, t_l, p, points)
    else
      call qsat_wat(qsl_tl, t_l, p, points)
    end if

    do i = 1, points
      if (qcl(i) >  zero .and. cfliq(i) >  zero) then

          ! Autoconversion is active

          !-----------------------------------------------
          ! De-bias the autoconversion rate based on a
          ! triangular Smith PDF, following Wood et al
          ! (Atmos. Res., 65, 109-128, 2002).
          !-----------------------------------------------

        if (qcl(i)  <   1.0e-15_real_lsprec) then
            ! Water contents are small so set factor to 1
          ac_factor=one

        else  ! qcl lt 1e-15
          alpha_l = repsilon * lc * qsl_tl(i) / ( r * t_l(i)**2 )
          a_l = one / (one+(lcrcp*alpha_l))
          sigma_s = (one - rhcpt(i)) * a_l * qsl_tl(i) / sqrt(6.0_real_lsprec)
          g_l = 1.15_real_lsprec * (power_qcl_auto-one) * sigma_s
          gacb = exp(-one * qcl(i) / g_l)

            ! Calculate autoconversion rate multiplication factor
          ac_factor = max(                                                     &
                  (cfliq(i)**(power_qcl_auto-one))*one/(one-gacb),             &
                   one)

        end if  ! qcl lt 1e-15

          !-----------------------------------------------
          ! Apply the debiasing factor
          !-----------------------------------------------
        autorate(i) = ac_factor * autorate(i)

      end if ! qcl > 0

    end do  ! points

  end if  ! l_auto_debias

    !-----------------------------------------------
    ! Part 4. Calculate the autoconversion
    !-----------------------------------------------

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do kk = 1, npts

    i = ix(kk)

      !-----------------------------------------------
      ! Set the dependence of autoconversion on air density
      !-----------------------------------------------
      ! power_rho_auto and power_qcl_auto are set in c_lspmic.
    rho_1 = rho(i) ** (power_rho_auto-power_qcl_auto+one)

        ! Calculate maximum amount of liquid that can be removed
        ! from the grid box
    qc = min( autolim(i) * cfliq(i) * rhor(i) , qcl(i) )

        !-----------------------------------------------
        ! Calculate the autoconversion amount
        !-----------------------------------------------
    dpr(i) = min(autorate(i)                                                   &
              *(rho(i)*qcl(i)/cfliq(i))**(power_qcl_auto-one)                  &
              *timestep*qcl(i)*rho_1/corr2(i),qcl(i)-qc)

      !-----------------------------------------------
      ! Update liquid water content and rain
      !-----------------------------------------------
    qcl(i)   = qcl(i)   - dpr(i)
    qrain(i) = qrain(i) + dpr(i)

      !-----------------------------------------------
      ! Store process rate
      !-----------------------------------------------
    ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

    if (dpr(i) >  zero .and. (.not. l_mcr_precfrac) ) then
      ! (not doing this here if using prognostic precip fraction,
      !  as in that case rainfrac is updated later)
        !-----------------------------------------------
        ! Calculate change in rain fraction
        !-----------------------------------------------
      rainfracnew = max(rainfrac(i),cfliq(i))
      rftransfer(i) = rftransfer(i)                                            &
                      + (rainfracnew - rainfrac(i)) * one_over_tsi

        !-----------------------------------------------
        ! Update rain fractions
        !-----------------------------------------------
      rainfrac(i)   = rainfracnew
      rain_liq(i)   = min(area_liq(i),rainfrac(i))
      rain_mix(i)   = min(area_mix(i),rainfrac(i)-rain_liq(i))
      rain_ice(i)   =                                                          &
          min(area_ice(i),rainfrac(i)-rain_liq(i)-rain_mix(i))
      rain_clear(i) =                                                          &
          rainfrac(i)-rain_liq(i)-rain_mix(i)-rain_ice(i)

    end if  ! dpr(i) gt 0

  end do

end if ! l_warm_new

if (l_wtrac) then     ! Store phase change amount for water tracer use
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do kk = 1, npts
    i = ix(kk)
    wtrac_mp_cpr_old%qchange(i) = dpr(i)
  end do
end if

if ( l_mcr_precfrac ) then
  ! If using prognostic precipitation fraction...
  if ( i_update_precfrac == i_homog_areas ) then
    ! Need to store contribution of rain autoconversion to the
    ! total precip increment in the liquid-cloud partition.
    do kk = 1, npts
      i = ix(kk)
      ! If part of the liquid-cloud area doesn't already contain rain/graupel,
      ! need to assign part of the precip mass increment to "new" area...

      ! Compute area of newly-created rain
      area_new = max( zero, area_liq(i) - rain_liq(i) )                        &
               + max( zero, area_mix(i) - rain_mix(i) )

      ! Compute 1.0 / total area producing rain by autoconversion
      tmp = one / max( rain_liq(i) + rain_mix(i) + area_new, small_number )

      ! Increment applies in the rain_liq, rain_mix  and new area partitions;
      ! divvy it up between them in proportion to area
      dqprec_liq(i) = dqprec_liq(i) + dpr(i) * rain_liq(i) * tmp
      dqprec_mix(i) = dqprec_mix(i) + dpr(i) * rain_mix(i) * tmp
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
      area_inc(i) = area_liq(i) + area_mix(i)
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
  end if  ! ( i_update_precfrac )
end if  ! ( l_mcr_precfrac )

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_autoc
end module lsp_autoc_mod
