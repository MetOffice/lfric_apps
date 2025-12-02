! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Accretion of droplets by raindrops
! Subroutine Interface:
module lsp_accretion_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_ACCRETION_MOD'

contains

subroutine lsp_accretion(                                                      &
  points, timestep,                                                            &
                                          ! Number of points and tstep
  qgraup, qcl, qrain,                                                          &
                                          ! Water contents
  cfliq,                                                                       &
                                          ! Cloud fraction information
                                          ! at start of microphysics ts
  rainfrac, rain_liq, rain_mix,                                                &
                                          ! Rain fraction information
  rho, corr, dhir, rain_nofall,                                                &
                                          ! Parametrization information
  ptransfer,                                                                   &
                                          ! Mass transfer diagnostic
  one_over_tsi,                                                                &
                                          ! 1/(timestep*iterations)
  r_theta_levels_c, fv_cos_theta_latitude_c,                                   &
  f_arr1, f_arr2, f_arr3,                                                      &
  l_update_dqprec, dqprec_liq, dqprec_mix, precfrac_k,                         &
  wtrac_mp_cpr_old                                                             &
                                          ! Store used for water tracers
  )

!Use in reals in lsprec precision, both microphysics related and general atmos
use lsprec_mod, only: cx, constp, acc_pref, acc_qc, acc_qr, qcfmin,            &
                      f_cons, fsd_eff_lam, fsd_eff_phi, rad_mcica_sigma,       &
                      two_d_fsd_factor, zero, half, one, small_number

  ! Microphysics Modules- logicals and integers
use mphys_constants_mod, only: l_inhomog
use mphys_inputs_mod,    only: l_warm_new, c_r_correl,                         &
                               l_mcr_qrain, l_fsd_generator, l_mcr_precfrac,   &
                               l_subgrid_graupel_frac, i_update_precfrac,      &
                               i_homog_areas, i_sg_correl

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,            only: real_lsprec

! Water tracers
use free_tracers_inputs_mod, only: l_wtrac
use wtrac_mphys_mod,         only: mp_cpr_old_wtrac_type

! Dr Hook Modules
use yomhook,             only: lhook, dr_hook
use parkind1,            only: jprb, jpim

! FSD parameters module- logicals and integers
use fsd_parameters_mod,  only: ip_fsd_constant
! Constant FSD value and ratio between 1D and 2D FSD
use rad_input_mod,       only: i_fsd

use lsp_combine_precfrac_mod, only: lsp_combine_precfrac

implicit none

! Purpose:
!   Update cloud prognostics as a result of accretion of cloud
!   droplets by rain drops

! Method:
!   Solve implicitly the microphysical transfer equation for a
!   specified distribution of raindrops sweeping out a volume
!   of air uniformally inhabited by cloud water droplets.
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
  one_over_tsi,                                                                &
                        ! 1/(timestep*iterations)
  rain_liq(points),                                                            &
                        ! Overlap fraction of rain and liquid
  rain_mix(points),                                                            &
                        ! Overlap fraction of rain and mixed phase
  rainfrac(points),                                                            &
                        ! Rain fraction
  cfliq(points),                                                               &
                        ! Liquid cloud fraction at start of timestep
    rho(points),                                                               &
                          ! Air density / kg m-3
    corr(points),                                                              &
                          ! Fall velocity correction factor (no units)
    r_theta_levels_c(points),                                                  &
                          ! Distance from centre of Earth and ...
    fv_cos_theta_latitude_c(points),                                           &
                          ! ... grid info for working out gridbox size.
    dhir(points),                                                              &
                          ! Depth of layer / timestep  / m s-1
    rain_nofall(points),                                                       &
                          ! Fraction of the rain-mass that is not falling out
    f_arr1(points),                                                            &
    f_arr2(points),                                                            &
    f_arr3(points),                                                            &
    qgraup(points)
                        ! Graupel content / kg kg-1

real (kind=real_lsprec), intent(in out) ::                                     &
  qcl(points),                                                                 &
                        ! Liquid water content / kg kg-1
  qrain(points),                                                               &
                        ! Rain mixing ratio / kg kg-1
    ptransfer(points) ! Mass rimed in this timestep / kg kg-1


! Flag for whether this call actually updates qrain
! (used to disable updating of the dqprec stores below when this
!  routine is called in a non-updating way in the orographic
!  seeder-feeder calculation).
logical, intent(in) :: l_update_dqprec

! Precipitation mass increment within each sub-region of the precipitation
! fraction (bits that overlap liquid-cloud).
! Used to update the prognostic precipitation fraction.
real (kind=real_lsprec), intent(in out) :: dqprec_liq(points)
real (kind=real_lsprec), intent(in out) :: dqprec_mix(points)
! Prognostic precip fraction
real (kind=real_lsprec), intent(in out) :: precfrac_k(points)

type(mp_cpr_old_wtrac_type), intent(in out) :: wtrac_mp_cpr_old
                        ! Water tracers store

! Local Variables

integer ::                                                                     &
  i                 ! Loop counter for points

real (kind=real_lsprec) ::                                                     &
  dpr(points),                                                                 &
                        ! Transfer of mixing ratio  / kg kg-1
  temp1,                                                                       &
                        ! Temporary in rate calculations
  lambda,                                                                      &
                        ! Temporary lambda for Abel/Shipway
  lambda_h1,                                                                   &
                        ! Temporary lambda + h1r for Abel/Shipway
  lambda_h2,                                                                   &
                        ! Temporary lambda + h2r for Abel/Shipway
  temp7,                                                                       &
                        ! Rain free liquid cloud fraction (no units)
  qclnew            ! Updated value of liquid water / kg kg-1

real (kind=real_lsprec) ::                                                     &
  fsd_qc(points),                                                              &
          ! fractional standard dev of cloud water content
  fsd_qr(points),                                                              &
          ! fractional standard dev of rain water content
  bias,                                                                        &
          ! accretion rate bias
  x_in_km ! grid-box size in km

! Fraction of increment occuring in the liquid-only vs mixed-phase cloud
real(kind=real_lsprec) :: frac
! Area fraction where rain created
real(kind=real_lsprec) :: area_inc(points)
! Precip content (rain+graupel) before the increment is applied
real(kind=real_lsprec) :: qprec(points)

! Local compression variable
integer ::                                                                     &
  npts,                                                                        &
                        ! Number of points to compute
  c,                                                                           &
                        ! Compressed point index
  ix(points)
                        ! Original index

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_ACCRETION'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! When working in 32-bit, a Cray compiler bug on power breaks PROC
! comparability.  Conditionally using NOVECTOR makes this go
! away. This directive will be valid until the end of the program unit
! or when a VECTOR directive (without any extra string) is used.
! Review upon upgrade.
#if defined(LSPREC_32B) && defined (CRAYFTN_VERSION) && (CRAYFTN_VERSION <8004000)
!DIR$ NOVECTOR
#endif

! In all the following loops, the 'i' variable is used for the
! original index while the 'c' variable is used for the compressed
! index.

! Identify the points that need to be calculated.
npts = 0
do i = 1, points

  if (rainfrac(i) >  zero .and. qrain(i) >  zero                               &
      .and. cfliq(i) >  zero) then

    npts = npts + 1
    ix(npts) = i

  else

    dpr(i) = zero

  end if

end do

if (l_warm_new) then

  !---------------------------------------------
  ! Use new accretion scheme
  !---------------------------------------------
  if (l_inhomog) then

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do c = 1, npts

      i = ix(c)

      ! calculate bias E factor based on Boutle et al 2012 QJ
      x_in_km = 0.001_real_lsprec*sqrt (r_theta_levels_c(i) * fsd_eff_lam      &
                             * r_theta_levels_c(i) * fsd_eff_phi               &
                             * fv_cos_theta_latitude_c(i)     )

      if (l_fsd_generator) then ! Use same FSD param as cld generator
        if (i_fsd == ip_fsd_constant) then
          fsd_qc(i) = rad_mcica_sigma
        else
          if ( cfliq(i) < one ) then
            fsd_qc(i) = (f_arr2(i)-f_arr3(i)*cfliq(i))                         &
              *(((x_in_km*cfliq(i))**0.333_real_lsprec)                        &
              *((f_cons(1)*x_in_km*cfliq(i))**f_cons(2)+one)                   &
             **(f_cons(3)))
          else
            fsd_qc(i) = f_arr1(i) * ((x_in_km**0.333_real_lsprec)              &
              *((f_cons(1)*x_in_km)**f_cons(2)+one)**(f_cons(3)))
          end if
        end if
      else ! Use FSD param from Boutle et al 2012 QJ
        if ( cfliq(i) < one ) then
          fsd_qc(i)=(0.45_real_lsprec-0.25_real_lsprec*cfliq(i))               &
                  *(((x_in_km*cfliq(i))**0.333_real_lsprec)                    &
                 *((0.06_real_lsprec*x_in_km*cfliq(i))**1.5_real_lsprec+one)   &
                 **(-0.17_real_lsprec))
        else
          fsd_qc(i)=0.11_real_lsprec*(((x_in_km*cfliq(i))**0.333_real_lsprec)  &
                   *((0.06_real_lsprec*x_in_km*cfliq(i))                       &
                   **1.5_real_lsprec+one)**(-0.17_real_lsprec))
        end if
      end if

      fsd_qr(i) = (1.1_real_lsprec-0.8_real_lsprec*rainfrac(i))                &
               *(((x_in_km*rainfrac(i))**0.333_real_lsprec)                    &
               *((0.11_real_lsprec*x_in_km*rainfrac(i))                        &
               **1.14_real_lsprec+one)**(-0.22_real_lsprec))

      fsd_qc(i) = fsd_qc(i)*two_d_fsd_factor
      fsd_qr(i) = fsd_qr(i)*two_d_fsd_factor

    end do
  else ! no inhomog param

    bias = one

  end if

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do c = 1, npts

    i = ix(c)

    if (l_inhomog) then

      bias = ((one+fsd_qc(i)**2)**(-half*acc_qc))*                             &
             ((one+fsd_qc(i)**2)**(half*acc_qc**2))*                           &
             ((one+fsd_qr(i)**2)**(-half*acc_qr))*                             &
             ((one+fsd_qr(i)**2)**(half*acc_qr**2))*                           &
             exp(c_r_correl*acc_qc*acc_qr*                                     &
             sqrt(log(one+fsd_qc(i)**2)*                                       &
             log(one+fsd_qr(i)**2)))

    end if


    dpr(i) = acc_pref * bias *                                                 &
             ( (qcl(i)/cfliq(i))**acc_qc) *                                    &
             ( (qrain(i)*rain_nofall(i)/rainfrac(i))**acc_qr )                 &
             * timestep * (rain_liq(i)+rain_mix(i))
    dpr(i) = max(min(dpr(i),                                                   &
             qcl(i)*(rain_liq(i)+rain_mix(i))/cfliq(i)-qcfmin),zero)

    ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

    !------------------------------------------------
    ! Adjust liquid and rain contents
    !------------------------------------------------
    qrain(i) = qrain(i) + dpr(i)
    qcl(i)   = qcl(i)   - dpr(i)

  end do

else ! original accretion

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do c = 1, npts

    i = ix(c)


    if (l_mcr_qrain) then
            ! Use the prognostic rain formulation
      temp1 = qrain(i)*rain_nofall(i)*rho(i)*constp(50)/(rainfrac(i))
    else
            ! Use the diagnostic rain formulation
      temp1 = qrain(i)*rho(i)*dhir(i)                                          &
               /(rainfrac(i)*constp(42)*corr(i))

    end if  ! l_mcr_qrain

    !-----------------------------------------------
    ! Calculate the new local value of qcl
    !-----------------------------------------------

    if (l_mcr_qrain) then

      !Prognostic Abel and Shipway version

      !First need to calculate lambda:
      ! lambda = (1 / (rain fraction temp1)) to power of cx(52)

      lambda = ( one / (rainfrac(i)*temp1) )**cx(52)

      lambda_h1 = lambda+cx(56)
      lambda_h2 = lambda+cx(57)

      qclnew = qcl(i) / (                                                      &
                 cfliq(i)+(cfliq(i)*timestep*corr(i)*  (                       &
                 (constp(51) * (lambda**cx(46))) /                             &
                  (lambda_h1**cx(61))     +                                    &
                 (constp(52) * (lambda**cx(46))) /                             &
                  (lambda_h2**cx(62))     )  )                                 &
                         )


    else

         !Diagnostic Abel and Shipway version
         !First need lambda, which for this case is
         !lambda = (1 /(rain fraction temp1)) to power of cx(48)

      lambda = ( one / (rainfrac(i)*temp1) )**cx(48)

      lambda_h1 = lambda+cx(56)
      lambda_h2 = lambda+cx(57)

         !Now we have lambda, should be the same as above.
         !Keeping prognostic and diagnostic separate for now so
         !I can change one without affecting the other.

      qclnew = qcl(i) / (                                                      &
                cfliq(i)+(cfliq(i)*timestep*corr(i)*  (                        &
                (constp(51) * (lambda**cx(46))) /                              &
                 (lambda_h1**cx(61))     +                                     &
                (constp(52) * (lambda**cx(46))) /                              &
                 (lambda_h2**cx(62))     )  )                                  &
                        )

    end if ! l_mcr_qrain

        !-----------------------------------------------
        ! Convert qclnew to a gridbox mean
        !-----------------------------------------------
        ! temp7 is the rain free region of liquid cloud
    temp7 = max(cfliq(i)-rain_liq(i)-rain_mix(i),zero)
    qclnew = qclnew*(rain_liq(i)+rain_mix(i))                                  &
               +qcl(i)/cfliq(i)*temp7

        !-----------------------------------------------
        ! Calculate transfer
        !-----------------------------------------------
    dpr(i) = qcl(i) - qclnew

    ptransfer(i) = ptransfer(i) + dpr(i) * one_over_tsi

      !------------------------------------------------
      ! Adjust liquid and rain contents
      !------------------------------------------------
    qrain(i) = qrain(i) + dpr(i)
    qcl(i)   = qcl(i)   - dpr(i)

  end do

end if ! l_warm_new

if (l_wtrac .and. l_update_dqprec) then
  ! Store phase change amount for water tracer use
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
  do c = 1, npts
    i = ix(c)
    wtrac_mp_cpr_old%qchange(i) = dpr(i)
  end do
end if

if ( l_mcr_precfrac .and. l_update_dqprec ) then
  ! If using prognostic precipitation fraction...
  if ( i_update_precfrac == i_homog_areas ) then
    ! Need to store contribution of rain accretion to the
    ! total precip increment in the liquid-cloud partition.
!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
    do c = 1, npts
      i = ix(c)
      ! Increment applies in both the rain_liq and rain_mix area partitions;
      ! divvy it up between the two in proportion to area
      frac = rain_liq(i) / max( rain_liq(i) + rain_mix(i), small_number )
      dqprec_liq(i) = dqprec_liq(i) +      frac  * dpr(i)
      dqprec_mix(i) = dqprec_mix(i) + (one-frac) * dpr(i)
    end do
  else if ( i_update_precfrac == i_sg_correl ) then
    ! Set precip mass before the increment was applied
    ! (ignoring negative values), and area fraction of the increment
    do i = 1, points
      qprec(i) = max( qrain(i) - dpr(i), zero )
      dpr(i) = max( qrain(i), zero ) - qprec(i)
      area_inc(i) = rain_liq(i) + rain_mix(i)
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
end subroutine lsp_accretion
end module lsp_accretion_mod
