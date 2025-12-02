! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_qsat
!
! Included via qsat_mod_qsat.h

! Subroutine Arguments:
integer, intent(in)                 :: npnts
real (kind=prec),    intent(in)     :: t(npnts), p(npnts)
real (kind=prec),    intent(in out) :: qs(npnts)
#if defined(QSAT_WITH_IDX)
! Indirect indexing
integer,             intent(in)    :: idx(npnts)
integer,             intent(in)    :: ni
#endif


! Local scalars
integer                          :: itable, i
real (kind=prec)                 :: atable, fsubw, tt

! Local parameters allows simple templating of KINDs
real (kind=prec), parameter      :: one   = 1.0_prec,                          &
                                    pconv = 1.0e-8_prec,                       &
                                    term1 = 4.5_prec,                          &
                                    term2 = 6.0e-4_prec


#if defined(QSAT_WITH_IDX)
integer                          :: m

!DIR$ IVDEP
!DIR$ VECTOR ALWAYS
do m=1, ni
  i = idx(m)
#else
do i=1, npnts
#endif
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
  itable = int(atable)
  atable = atable - real(itable, kind=prec)
  qs(i)  = (one - atable) * es(itable) + atable * es(itable+1)

  ! Multiply by fsubw to convert to saturated vapour pressure in air
  ! (equation A4.6 OF Adrian Gill's book).
  qs(i)  = qs(i) * fsubw

  ! Now form the accurate expression for qs, which is a rearranged
  ! version of equation A4.3 of Gill's book.
  ! Note that at very low pressures we apply a fix, to prevent a
  ! singularity (qsat tends to 1. kg/kg).
  qs(i)  = (repsilon * qs(i)) / (max(p(i), qs(i)) - one_minus_epsilon * qs(i))
end do
return
