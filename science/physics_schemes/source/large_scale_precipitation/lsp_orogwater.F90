! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Large-scale precipitation scheme. Subgrid orographic water calculation

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

! Subroutine Interface:
module lsp_orogwater_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_OROGWATER_MOD'

contains

subroutine lsp_orogwater(                                                      &
  points,                                                                      &
                     ! Number of points
  hmteff,                                                                      &
                     ! Effective peak-to-trough mountain height (m)
  zb,                                                                          &
                     ! Blocked layer depth (m)
  r_theta_levels_c,r_theta_surf_c, fv_cos_theta_latitude_c,                    &
                     ! For obtaining altitude
  p, t, q, qsl, esw, qcl,                                                      &
                     ! Grid-box thermodynamics
  ql_orog, cf_orog                                                             &
                     ! Orographic cloud water MR (kg/kg) and fraction
  )

use um_types,             only: real_lsprec
use lsprec_mod,           only: zerodegc, pi, delta_lambda, delta_phi, lc,     &
                                qcfmin, repsilon, r, cp, g, rv,                &
                                zero, one, two, half,                          &
! Horizontal scale of sub-grid sinusoidal ridge is (nscalesf * grid-length)
! in the vertical decay equation - longer scales decay more slowly
                                 nscalesf

! Dr Hook Modules
use yomhook,              only: lhook, dr_hook
use parkind1,             only: jprb, jpim

implicit none

!------------------------------------------------------------------------------
! Purpose:
!   Calculates sub-grid orographic water MR and cloud fraction
!   No standard model variables are modified here.
!   Assumes moist neutral flow over a 2D ridge of height determined by sd_orog
!   Its horizontal wavelength is assumed equal to the grid-spacing, which
!   determines the rate at which orographic displacements decay with altitude.
!   If air is already saturated and some cloud water is present then the LCL
!   is assumed to lie below the model level at the distance required to
!   evaporate qcl,  Otherwise the LCL is above the model level by an amount
!   determined by the difference in the T and Td lapse rates.

!   Documentation: UMDP 026

!------------------------------------------------------------------------------

! Subroutine Arguments

integer, intent(in) :: points
!   Number of points calculation is done over

real (kind=real_lsprec), intent(in) ::                                         &
  hmteff(points),                                                              &
!         Effective mountain peak-to-trough height (metres)
  zb(points),                                                                  &
!         Sub-grid orographic blocked layer depth (metres)
  p(points),                                                                   &
!         Air pressure at this level (Pa).
  qcl(points),                                                                 &
!         Cloud liquid water (kg water per kg air).
  t(points),                                                                   &
!         Temperature at this level (K).
  q(points),                                                                   &
!         Specific humidity at this level (kg water per kg air).
  qsl(points),                                                                 &
!         Saturation specific humidity (kg water per kg air).
  esw(points),                                                                 &
!         Saturation vapour pressure wrt water at all T (Pa).
  r_theta_levels_c(points),                                                    &
                     ! Distance from centre of the Earth
  r_theta_surf_c(points),                                                      &
                     ! ...and near surface (k=0)
  fv_cos_theta_latitude_c(points)
                     ! Finite volume cosine of latitude.

real (kind=real_lsprec), intent(out) ::                                        &
  ql_orog(points),                                                             &
!         Grid-box mean orographic cloud water MR (kg/kg)
  cf_orog(points)
!        Orographic cloud fraction

!------------------------------------------------------------------------------
! Local Variables

integer :: i

real (kind=real_lsprec) :: altitude

real (kind=real_lsprec) :: gridspacing,waveno,tdew,dzcond,salr,                &
        xlim,ydis,dqldz,maxdisp,gridspacing_surf,                              &
        tmidcld,pmidcld,esatmidcld,qsatmidcld,                                 &
        dewlr,satdz,evapdist,satdzmin,hmin,tc,tdewc,tmidcldc

real (kind=real_lsprec),parameter ::                                           &
! Constant in Tv equation used within tmidcld equation
    d1 = 0.61,                                                                 &
! Dry Adiabatic Lapse Rate (K/km)
    dalr = 9.8e-3,                                                             &
! Constants in dewpoint and saturation vapour pressure equations as
! suggested by Alduchov and Eskridge (1996)
    a1 = 17.625,                                                               &
    b1 = 243.04,                                                               &
! Multiple in saturation vapour pressure estimate (Pa)
    c1 = 610.94

! Declarations for Dr Hook:
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_OROGWATER'

!==============================================================================
! Start of calculations

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise arrays passed back to lsp_ice
do i = 1, points
  ql_orog(i) = zero
  cf_orog(i) = zero
end do  ! For each cloudy datapoint

! Minimum subgrid orography able to enhance rain
hmin = one
satdzmin = 0.1_real_lsprec

! Calculate sub-grid orographic water MR when needed
! (Need rain or snow to be able to enhance accretion)
do i = 1, points

  !   Only create extra water if above blocked layer depth zb
  altitude = max( zero, r_theta_levels_c(i) - r_theta_surf_c(i))

  !   Only give non-zero ql_orog if above the blocking level
  !   Do not proceed if negative Q has been passed in as this will
  !    produce NaNs and cause the model to crash
  if ( altitude > zb(i) .and. (half*hmteff(i)) > hmin  .and. q(i) > zero ) then

    !         Horizontal gridspacing at this altitude for averaging
    gridspacing = sqrt ( r_theta_levels_c(i) * delta_lambda                    &
                       * r_theta_levels_c(i) * delta_phi                       &
                       * fv_cos_theta_latitude_c(i)  )

    !         Horizontal surface gridspacing to characterise orography
    gridspacing_surf = sqrt ( r_theta_surf_c(i) * delta_lambda                 &
                            * r_theta_surf_c(i) * delta_phi                    &
                            * fv_cos_theta_latitude_c(i)  )
    waveno = two * pi / gridspacing_surf

    !         Max vertical orographic displacement accounting for exp decay
    !         Effective hill height is peak-to-trough value so halve it
    !         Horizontal wavelength used by decay function is scaled by
    !         nscalesf (longer scales decay more slowly aloft)

    maxdisp = half*hmteff(i) * exp(-altitude*waveno/nscalesf)

    !         Cloudy grid-box:  saturated with water present
    if ( qcl(i) > qcfmin .and. q(i)/qsl(i) > one ) then

      !           Rate of change of Qsat per metre of saturated descent
      dqldz = max( zero,  g*(one + lc*q(i)/r/t(i))/                            &
         (cp + repsilon*lc*lc*q(i)/r/t(i)/t(i))                                &
         *(repsilon + qsl(i))*qsl(i)*lc/r/t(i)/t(i)                            &
         -qsl(i)*p(i)*g/(p(i) - esw(i))/r/t(i)  )

      evapdist = qcl(i) / dqldz

      if (maxdisp > evapdist) then

        xlim = (one/waveno)*acos(-evapdist/maxdisp)
        ydis = maxdisp*two*sin(waveno*xlim)/waveno                             &
                 - two*evapdist*((half*gridspacing)-xlim)
        ql_orog(i) = max(zero, dqldz*ydis/gridspacing)
        !             Make sure orog water is no larger than available vapour
        ql_orog(i) = min( ql_orog(i), q(i) )

        !             Cloud fraction reduced by descent, assuming cfliq=1
        cf_orog(i) = max(zero, two*xlim/gridspacing)

      end if

    else     ! Clear gridbox - no resolved cloud

      !           Convert temperature to C for use in tdew
      tc = t(i)-zerodegc

      !           Dewpoint temperature (C)
      tdewc = b1*( log(q(i)/qsl(i)) + (a1*tc/(b1 + tc)) ) /                    &
          ( a1 - log(q(i)/qsl(i)) - (a1*tc/(b1 + tc)) )

      !           Convert dewpoint temperature to K
      tdew = tdewc + zerodegc

      !           Dewpoint temperature adiabatic lapse rate
      dewlr = (tdew * tdew * g * rv)/(lc * r * t(i))

      !           Condensation level (ascent required to reach saturation)
      dzcond = max( (t(i) - tdew)/(dalr - dewlr), zero)

      !           Only carry on if there is saturated ascent
      !           Will need distance to mid-cloud level satdz
      satdz = half * (maxdisp - dzcond)

      if (satdz > satdzmin) then

        !             Estimate T and SALR at cloud base - use q which remains
        !             constant during dry ascent and = (qsat at LCL)
        tmidcld = t(i) - (dzcond * dalr)

        salr = g * (one + (lc*q(i))/(r*tmidcld) ) /                            &
               (cp + (lc*lc*q(i)*repsilon)/(r*tmidcld*tmidcld))

        !             Find variables at mid-cloud level - further decrease T
        tmidcld = tmidcld - (satdz*salr)

        !             New pressure from total dry+moist ascent (use mean Tv)
        pmidcld = p(i) * exp( -( (dzcond + satdz) * g)/                        &
                    (r*half*(t(i) + tmidcld)*(one + d1*q(i)) ) )

        !             Convert mid-cloud temperature to C
        tmidcldc = tmidcld - zerodegc

        !             New esat and qat to get dqldz
        esatmidcld = c1 * exp( a1*tmidcldc/(tmidcldc + b1) )

        qsatmidcld = repsilon * esatmidcld / pmidcld

        salr = g * (one + (lc*qsatmidcld)/(r*tmidcld) ) /                      &
               (cp + (lc*lc*qsatmidcld*repsilon)/                              &
               (r*tmidcld*tmidcld))

        dqldz = max( zero, salr*(repsilon + qsatmidcld)                        &
                  *qsatmidcld*lc/r/tmidcld/tmidcld                             &
                  -qsatmidcld*pmidcld*g/                                       &
                  (pmidcld - esatmidcld)/r/tmidcld )

        xlim = (one/waveno)*acos(dzcond/maxdisp)
        ydis = maxdisp*two*sin(waveno*xlim)/waveno                             &
                  - (two*xlim*dzcond)

        ql_orog(i) = max(zero, dqldz*ydis/gridspacing)
        ql_orog(i) = min( ql_orog(i), q(i) )

        cf_orog(i) = max(zero, two*xlim/gridspacing)

      end if ! Water streamline above corrected dzcond

    end if  ! Clear or cloudy gridbox

  end if

end do  ! For each cloudy datapoint

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_orogwater
end module lsp_orogwater_mod
