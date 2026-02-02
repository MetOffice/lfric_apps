! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

#if !defined(LFRIC)
module set_qsat_mod

implicit none

contains

! UM qsat functions for mixing ratios

!----------------------------------------------------------------
! Returns saturation mixing ratio with respect to liquid at all
! temperatures
!----------------------------------------------------------------
subroutine set_qsat_liq( npnts, t, p, qs )

use comorph_constants_mod, only: real_cvprec, R_dry, R_vap,                    &
                     zerodegc => melt_temp
use qsat_data, only: t_low, t_high, delta_t, es_wat

implicit none

! Subroutine Arguments:
integer, intent(in) :: npnts
real(kind=real_cvprec), intent(in)  :: t(npnts), p(npnts)
real(kind=real_cvprec), intent(out) :: qs(npnts)

! Local scalars
integer                 :: itable, i
real(kind=real_cvprec) :: atable, fsubw, tt

! Local parameters allows simple templating of KINDs
real(kind=real_cvprec), parameter :: one   = 1.0_real_cvprec,                  &
                                     pconv = 1.0e-8_real_cvprec,               &
                                     term1 = 4.5_real_cvprec,                  &
                                     term2 = 6.0e-4_real_cvprec,               &
                                     term3 = 1.1_real_cvprec

! Ratio of gas constants for dry air and water vapour
real(kind=real_cvprec) :: repsilon

repsilon = R_dry / R_vap

do i=1, npnts
  ! Compute the factor that converts from sat vapour pressure in a
  ! pure water system to sat vapour pressure in air, fsubw.
  ! This formula is taken from equation A4.7 of Adrian Gill's book:
  ! atmosphere-ocean dynamics. Note that his formula works in terms
  ! of pressure in mb and temperature in celsius, so conversion of
  ! units leads to the slightly different equation used here.
  fsubw = one + pconv * p(i) *                                                 &
  (term1 + term2 * (t(i) - zerodegc) * (t(i) - zerodegc))

  ! Use the lookup table to find saturated vapour pressure. Store it in qs.
  tt     = max(t_low,t(i))
  tt     = min(t_high,tt)
  atable = (tt - t_low + delta_t) / delta_t
  itable = int( atable )
  atable = atable - real( itable, real_cvprec )
  qs(i)  = (one - atable) * es_wat(itable) + atable * es_wat(itable+1)

  ! Multiply by fsubw to convert to saturated vapour pressure in air
  ! (equation A4.6 OF Adrian Gill's book).
  qs(i)  = qs(i) * fsubw

  ! Now form the accurate expression for qs, which is a rearranged
  ! version of equation A4.3 of Gill's book.
  ! Note that at very low pressures we apply a fix, to prevent a
  ! singularity (qsat tends to 1. kg/kg).
  qs(i)  = (repsilon * qs(i)) / (max(p(i), term3 * qs(i)) - qs(i))
end do

return
end subroutine set_qsat_liq


!----------------------------------------------------------------
! Returns saturation mixing ratio with respect to ice when
! below freezing; defaults to liquid where above freezing
!----------------------------------------------------------------
subroutine set_qsat_ice( npnts, t, p, qs )

use comorph_constants_mod, only: real_cvprec, R_dry, R_vap,                    &
                     zerodegc => melt_temp
use qsat_data, only: t_low, t_high, delta_t, es

implicit none

! Subroutine Arguments:
integer, intent(in) :: npnts
real(kind=real_cvprec), intent(in)  :: t(npnts), p(npnts)
real(kind=real_cvprec), intent(out) :: qs(npnts)

! Local scalars
integer                :: itable, i
real(kind=real_cvprec) :: atable, fsubw, tt

! Local parameters allows simple templating of KINDs
real(kind=real_cvprec), parameter :: one   = 1.0_real_cvprec,                  &
                                     pconv = 1.0e-8_real_cvprec,               &
                                     term1 = 4.5_real_cvprec,                  &
                                     term2 = 6.0e-4_real_cvprec,               &
                                     term3 = 1.1_real_cvprec

! Ratio of gas constants for dry air and water vapour
real(kind=real_cvprec) :: repsilon

repsilon = R_dry / R_vap

do i=1, npnts
  ! Compute the factor that converts from sat vapour pressure in a
  ! pure water system to sat vapour pressure in air, fsubw.
  ! This formula is taken from equation A4.7 of Adrian Gill's book:
  ! atmosphere-ocean dynamics. Note that his formula works in terms
  ! of pressure in mb and temperature in celsius, so conversion of
  ! units leads to the slightly different equation used here.
  fsubw = one + pconv * p(i) *                                                 &
  (term1 + term2 * (t(i) - zerodegc) * (t(i) - zerodegc))

  ! Use the lookup table to find saturated vapour pressure. Store it in qs.
  tt     = max(t_low,t(i))
  tt     = min(t_high,tt)
  atable = (tt - t_low + delta_t) / delta_t
  itable = int( atable )
  atable = atable - real( itable, real_cvprec )
  qs(i)  = (one - atable) * es(itable) + atable * es(itable+1)

  ! Multiply by fsubw to convert to saturated vapour pressure in air
  ! (equation A4.6 OF Adrian Gill's book).
  qs(i)  = qs(i) * fsubw

  ! Now form the accurate expression for qs, which is a rearranged
  ! version of equation A4.3 of Gill's book.
  ! Note that at very low pressures we apply a fix, to prevent a
  ! singularity (qsat tends to 1. kg/kg).
  qs(i)  = (repsilon * qs(i)) / (max(p(i), term3 * qs(i)) - qs(i))
end do

return
end subroutine set_qsat_ice


end module set_qsat_mod
#endif
