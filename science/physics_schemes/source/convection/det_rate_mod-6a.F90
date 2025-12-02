! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Calculates forced detrainment rate in layer K

module det_rate_6a_mod

use um_types, only: real_umphys

implicit none

!
! Description:
!   Calculates forced detrainment rate in layer K
!
! Method:
!   See UM Documentation paper No 27
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 v8.3 programming standards.


character(len=*), parameter, private :: ModuleName='DET_RATE_6A_MOD'

contains

subroutine det_rate_6a (npnts, exkp1,                                          &
                        thek, thekp1, thpk, thpkp1, thrk,                      &
                        qek, qekp1, qpk, qpkp1, qrk,                           &
                        Qlkp1, Qfkp1, Frezkp1,                                 &
                        ekp14, ekp34,                                          &
                        deltak, idx, ni)

use water_constants_mod, only: lc, lf
use planet_constants_mod, only: cp
use cv_derived_constants_mod, only: ls
use cv_run_mod, only: fdet_opt
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use ereport_mod, only: ereport
use errormessagelength_mod, only: errormessagelength

implicit none

! Subroutine arguments

!----------------------------------------------------------------------
! Variables which are input
!----------------------------------------------------------------------

integer, intent(in) :: npnts     ! Number of points

real(kind=real_umphys),intent(in) :: exkp1(npnts)
                                 ! Exner ratio at mid-point of layer k+1
real(kind=real_umphys),intent(in) :: thek(npnts)
                                 ! Env. potential temperature in layer k (K)
real(kind=real_umphys),intent(in) :: thekp1(npnts)
                                 ! Env. potential temperature in layer k+1
real(kind=real_umphys),intent(in) :: thpk(npnts)
                                 ! Par. potential temperature in layer k (K)
real(kind=real_umphys),intent(in) :: thpkp1(npnts)
                                 ! Par. potential temperature in layer k+1 (K)
real(kind=real_umphys),intent(in) :: thrk(npnts)
                                 ! p. temperature of forced detrained
                                 ! parcel in layer k (K)
real(kind=real_umphys),intent(in) :: qek(npnts)
                                 ! Env. specific humidity in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qekp1(npnts)
                                 ! Env. spec. humidity in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: qpk(npnts)
                                 ! Par. specific humidity in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qpkp1(npnts)
                                 ! Par. specific humidity in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: qrk(npnts)
                                 ! Specific humidity of forced detrained
                                 ! parcel in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: Qlkp1(npnts)
                                 ! Amount of condensation to liquid water
                                 ! in the parcel (kg/kg)
real(kind=real_umphys),intent(in) :: Qfkp1(npnts)
                                 ! Amount of deposition to ice water
                                 ! in the parcel (kg/kg)
real(kind=real_umphys),intent(in) :: Frezkp1(npnts)
                                 ! Amount of freezing from liquid
                                 ! to ice water in the parcel (kg/kg)
real(kind=real_umphys),intent(in) :: ekp14(npnts)
                                 ! Entrainment coefficient at level k+1/4
                                 ! multiplied by appropriate layer thickness
real(kind=real_umphys),intent(in) :: ekp34(npnts)
                                 ! Entrainment coefficient at level k+3/4
                                 ! multiplied by appropriate layer thickness

!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------

real(kind=real_umphys),intent(in out) :: deltak(npnts)
                                    ! Parcel forced detrainment rate in
                                    ! layer k multiplied by layer thickness


! Indices of work points
integer, intent(in)  :: idx(npnts)
integer, intent(in)  :: ni

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------
real(kind=real_umphys) :: deltaqk(npnts)
                                 ! Parcel forced detrainment rate based on
                                 ! specific humidity eqn in
                                 ! layer k multiplied by layer thickness
real(kind=real_umphys) :: deltathk(npnts)
                                 ! Parcel forced detrainment rate based on
                                 ! potential temperature eqn in
                                 ! layer k multiplied by layer thickness

integer :: i,m                   ! loop counter
integer :: errorstatus           ! error status

real(kind=real_umphys) :: Factor
                    ! factor used to calculate the detrainment rate.
real(kind=real_umphys) :: Denom
                    ! denominator used to calculate the detrainment rate.
real(kind=real_umphys) :: Numer
                    ! numerator used to calculate the detrainment rate.
real(kind=real_umphys), parameter :: eps=epsilon(Denom)
                                      !Smallest allowable denominator

character(len=errormessagelength) :: cmessage    ! error message
character(len=*), parameter ::  RoutineName = 'DET_RATE_6A'

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle


!---------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!---------------------------------------------------------------------
! Calculate the forced detrainment rates for humidity and theta
! Both versions are required because checks are always made on
! whether either is outside of their 0.0 to 1.0 limits.
!---------------------------------------------------------------------

!---------------------------------------------------------------------
! Calculate the forced detrainment rate for humidity
!---------------------------------------------------------------------
!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  Factor = ekp14(i)*qek(i) + (1.0+ekp14(i))*ekp34(i)*qekp1(i)                  &
           - (1.0+ekp14(i))*(1.0+ekp34(i))                                     &
           * (qpkp1(i) + Qlkp1(i) + Qfkp1(i))

  Numer  = qpk(i) + Factor
  Denom  = qrk(i) + Factor

  if (abs(Denom) <= eps*abs(Numer)) then
    !If the denominator is close to zero then set
    !the detrainment rate to one (plus a bit to ensure
    !that this is taken account of later).
    deltaqk(i) = 1.1
  else
    deltaqk(i) = Numer/Denom
  end if

end do

!---------------------------------------------------------------------
! Calculate the forced detrainment rate for theta
!---------------------------------------------------------------------
!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  Factor = ekp14(i)*thek(i) + (1.0+ekp14(i))*ekp34(i)*thekp1(i)                &
           - (1.0+ekp14(i))*(1.0+ekp34(i))                                     &
           * (thpkp1(i) - (lc*Qlkp1(i) + ls*Qfkp1(i) + lf*Frezkp1(i))          &
           / (cp*exkp1(i)))

  Numer  = thpk(i) + Factor
  Denom  = thrk(i) + Factor

  if (abs(Denom) <= eps*abs(Numer)) then
    !If the denominator is close to zero then set
    !the detrainment rate to one (plus a bit to ensure
    !that this is taken account of later).
    deltathk(i) = 1.1
  else
    deltathk(i) = Numer/Denom
  end if

end do

!---------------------------------------------------------------------
! Copy the either the humidity or theta based detrainment rate to the
! output detrainment rate
!---------------------------------------------------------------------
if (fdet_opt == 0) then
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    deltak(i) = deltaqk(i)
  end do
else if (fdet_opt == 1 .or. fdet_opt == 2 .or. fdet_opt == 3) then
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    deltak(i) = deltathk(i)
  end do
else
  errorstatus = 1   ! will cause model to fail
  write (cmessage,'(a48)') 'Invalid value for fdet_opt. Valid values=0,1,2,3'
  call ereport(routinename, errorstatus, cmessage)
end if

if (fdet_opt == 0 .or. fdet_opt == 1) then
  !---------------------------------------------------------------------
  ! But if either deltaqk or deltathk are outside of the acceptable
  ! range then reset to either 0.0 or 1.0 as appropriate.
  !---------------------------------------------------------------------
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    if (deltaqk(i) <= 0.0 .or. deltathk(i) <= 0.0) then
      deltak(i) = 0.0
    else if (deltaqk(i) >= 1.0 .or. deltathk(i) >= 1.0) then
      deltak(i) = 1.0
    end if
  end do
else if (fdet_opt == 2 .or. fdet_opt == 3) then
  !---------------------------------------------------------------------
  ! If deltathk are outside of the acceptable
  ! range then reset to either 0.0 or 1.0 as appropriate.
  !---------------------------------------------------------------------
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    if (deltathk(i) <= 0.0) then
      deltak(i) = 0.0
    else if (deltaqk(i) >= 1.0) then
      deltak(i) = 1.0
    end if
  end do
end if


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine det_rate_6a
end module det_rate_6a_mod
