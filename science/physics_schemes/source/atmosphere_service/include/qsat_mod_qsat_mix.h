! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_qsat

! Subroutine Arguments:
integer, intent(in)              :: npnts
real (kind=prec),    intent(in)  :: t(npnts), p(npnts)
real (kind=prec),    intent(out) :: qs(npnts)

! Local scalars
integer                          :: itable, i
real (kind=prec)                 :: atable, fsubw, tt

! Local parameters allows simple templating of KINDs
real (kind=prec), parameter      :: one   = 1.0_prec,                          &
                                    pconv = 1.0e-8_prec,                       &
                                    term1 = 4.5_prec,                          &
                                    term2 = 6.0e-4_prec,                       &
                                    term3 = 1.1_prec

do i=1, npnts
  ! Compute the factor that converts from sat vapour pressure in a
  ! pure water system to sat vapour pressure in air, fsubw.
  ! This formula is taken from equation A4.7 of Adrian Gill's book:
  ! atmosphere-ocean dynamics. Note that his formula works in terms
  ! of pressure in mb and temperature in celsius, so conversion of
  ! units leads to the slightly different equation used here.
  fsubw = one + pconv * p(i) *                                                 &
  (term1 + term2 * (t(i) - zerodegc) * (t(i) - zerodegc))

  ! Required or crashes happen
  tt = t(i)
  if (t(i) < t_low) then
    tt = t_low
  else if (t(i) > t_high) then
    tt = t_high
  end if

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

  !Refactored to prevent decompostion errors
  if (p(i) > term3 * qs(i)) then
    qs(i)  = (repsilon * qs(i)) / (p(i) - qs(i))
  else
    qs(i)  = (repsilon * qs(i)) / ((term3 * qs(i)) - qs(i))
  end if
end do
return