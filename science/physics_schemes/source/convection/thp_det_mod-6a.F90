! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Calculates theta of the parcel in layer k+1 after forced detrainment
!
module thp_det_6a_mod

use um_types, only: real_umphys

implicit none

!
! Description:
!   Calculates potential temperature of the parcel in layer k+1
!   after forced detrainment in layer k.
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


character(len=*), parameter, private :: ModuleName='THP_DET_6A_MOD'

contains

subroutine thp_det_6a(npnts, exkp1, pkp1, thekp1, qekp1, watldekp1, watldpkp1, &
                      xsbmin, bwkp1, bgmkp1, qpkp1, thpkp1, idx, ni)

use water_constants_mod, only: lc
use cv_derived_constants_mod, only: ls
use planet_constants_mod, only: c_virtual, rv

use qsat_mod, only: qsat

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments

!----------------------------------------------------------------------
! Variables which are input
!----------------------------------------------------------------------

integer,intent(in) :: npnts         ! Number of points

real(kind=real_umphys),intent(in) :: exkp1(npnts)
                                    ! Exner ratio at mid-point of layer k+1
real(kind=real_umphys),intent(in) :: pkp1(npnts)
                                    ! pressure at mid-point of layer k+1 (Pa)
real(kind=real_umphys),intent(in) :: thekp1(npnts)
                                    ! Env. p. temperature in layer k+1 (K)
real(kind=real_umphys),intent(in) :: qekp1(npnts)
                                    ! Env. spec. humidity in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: watldekp1(npnts)
                                    ! Env. water loading in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: watldpkp1(npnts)
                                    ! Par. water loading in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: xsbmin(npnts)
                                    ! Threshold buoyancy for forced
                                    ! detrainment (K)

logical,intent(in) :: bwkp1(npnts)  ! Mask for whether condensate is
                                    ! liquid in layer k+1
logical,intent(in) :: bgmkp1(npnts) ! Mask for parcels which are
                                    ! saturated in layer k+1

!----------------------------------------------------------------------
! Variables which are input and output
!----------------------------------------------------------------------

real(kind=real_umphys),intent(in out) :: qpkp1(npnts)
                                     ! Par. spec. humidity in layer k+1 (kg/kg)
                                     ! after forced detrainment

!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------

real(kind=real_umphys),intent(in out) :: thpkp1(npnts)
                                       ! Par. p. temperature in layer k+1 (K)
                                       ! after forced detrainment


! Indices of work points
integer, intent(in)  :: idx(npnts)
integer, intent(in)  :: ni

!-------------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer :: i, j, m ! loop counters

real(kind=real_umphys) :: tpkp1(npnts)
                     ! estimate of the detrained parcel's temperature (K)
real(kind=real_umphys) :: qspkp1(npnts)! qsat at tpkp1
real(kind=real_umphys) :: dqsdth
                ! Rate of change of qsat with potential temperature
real(kind=real_umphys) :: el
                ! Latent heat of gas-to-whatever-condenses PC2 defn (J/kg)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='THP_DET_6A'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  if ( bgmkp1(i) ) then
    ! If saturated then an iterative calculation is required
    ! Set the first estimate of the parcel p.temp to that of the
    ! environment
    thpkp1(i) = thekp1(i)
    tpkp1(i)  = thpkp1(i) * exkp1(i)
  else
    ! If not saturated the detrained p.temp and humidity can be directly
    ! calculated
    thpkp1(i) = ( thekp1(i)*(1.0 + c_virtual*qekp1(i) - watldekp1(i))          &
                + xsbmin(i) ) / (1.0 + c_virtual*qpkp1(i) - watldpkp1(i))
    !tpkp1 is set to prevent qsat operating on uninitialised data.
    tpkp1(i)  = thpkp1(i) * exkp1(i)
    ! The humidity is unchanged.
  end if
end do

call qsat(qspkp1, tpkp1, pkp1, npnts, idx, ni)

do j=1,3  !Three iterations should be sufficient
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    if ( bgmkp1(i) ) then !if saturated
      if ( bwkp1(i) ) then
        el=lc
      else
        el=ls
      end if

      !Calculate the gradient of qsat. nb qspkp1 is saturated humidity at tpkp1
      dqsdth = el * qspkp1(i) / ( rv * exkp1(i) * thpkp1(i) * thpkp1(i) )

      !Calculate the next estimate of the parcel's p.temp
      thpkp1(i) = ( thekp1(i)*(1.0 + c_virtual*qekp1(i) - watldekp1(i))        &
                - thpkp1(i)*c_virtual*(qspkp1(i) - thpkp1(i)*dqsdth)           &
                + xsbmin(i) )                                                  &
                / ( 1.0 + c_virtual*thpkp1(i)*dqsdth - watldpkp1(i) )

      !update the parcel's temperature
      tpkp1(i)  = thpkp1(i) * exkp1(i)
    end if    !if saturated
  end do      !i

  ! Assume that the detrained parcel is saturated and update its humidity.
  call qsat(qspkp1, tpkp1, pkp1, npnts, idx, ni)

end do  !j

!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  if ( bgmkp1(i) ) then
    ! If saturated update qpkp1 to saturated humidity at tpkp1
    qpkp1(i) = qspkp1(i)
  end if
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine thp_det_6a
end module thp_det_6a_mod
