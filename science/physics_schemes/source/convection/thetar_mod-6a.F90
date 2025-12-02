! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Calculates theta and q of detraining air in layer k
!
module thetar_6a_mod

use um_types, only: real_umphys

implicit none

!
! Description:
!   Calculates the potential temperature and the humidity of the parcel
!   undergoing forced detrainment in layer k
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


character(len=*), parameter, private :: ModuleName='THETAR_6A_MOD'

contains

subroutine thetar_6a(npnts, exk, exkp1, pk, thek, thekp1, qek, qekp1,          &
                     qpk, watldek, watldekp1, watldpk, watldpkp1,              &
                     bwk, bgmk, thrk, qrk, idx, ni)

use water_constants_mod, only: lc
use cv_derived_constants_mod, only: ls
use planet_constants_mod, only: c_virtual, rv, recip_kappa, pref
use cv_run_mod, only: fdet_opt

use qsat_mod, only: qsat

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments

!----------------------------------------------------------------------
! Variables which are input
!----------------------------------------------------------------------

integer,intent(in) :: npnts        ! Number of points

real(kind=real_umphys),intent(in) :: exk(npnts)
                                   ! Exner ratio at mid-point of layer k
real(kind=real_umphys),intent(in) :: exkp1(npnts)
                                   ! Exner ratio at mid-point of layer k+1
real(kind=real_umphys),intent(in) :: pk(npnts)
                                   ! pressure at mid-point of layer k (Pa)
real(kind=real_umphys),intent(in) :: thek(npnts)
                                   ! Env. p. temperature in layer k (K)
real(kind=real_umphys),intent(in) :: thekp1(npnts)
                                   ! Env. pot. temperature in layer k+1 (K)
real(kind=real_umphys),intent(in) :: qek(npnts)
                                   ! Env. specific humidity in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qekp1(npnts)
                                   ! Env. spec. humidity in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: qpk(npnts)
                                   ! Par. specific humidity in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: watldek(npnts)
                                   ! Env. water loading in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: watldekp1(npnts)
                                   ! Env. water loading in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: watldpk(npnts)
                                   ! Par. water loading in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: watldpkp1(npnts)
                                   ! Par. water loading in layer k+1 (kg/kg)

logical,intent(in) :: bwk(npnts)   ! Mask for whether condensate is
                                   ! liquid in layer k+1
logical,intent(in) :: bgmk(npnts)  ! Mask for parcels which are
                                   ! saturated in layer k

!----------------------------------------------------------------------
! Variables which are output
!----------------------------------------------------------------------

real(kind=real_umphys),intent(in out) :: thrk(npnts)
                                      ! Pot. temperature of forced detrained
                                      ! parcel in layer k (K)
real(kind=real_umphys),intent(in out) :: qrk(npnts)
                                      ! Specific humidity of forced detrained
                                      ! parcel in layer k (kg/kg)


! Indices of work points
integer, intent(in)  :: idx(npnts)
integer, intent(in)  :: ni

!-------------------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

integer :: i, j, m ! loop counters

real(kind=real_umphys) :: trk(npnts)
                        ! estimate of the detrained parcel's temp (K)
real(kind=real_umphys) :: thve(npnts)
                        ! env. thetav that sets the detrained parcel's thv (K)
real(kind=real_umphys) :: qsrk(npnts)     ! qsat at trk
real(kind=real_umphys) :: watldr(npnts)
                        ! Par. water loading of detrained parcel (kg/kg)
real(kind=real_umphys) :: exr(npnts)      ! Exner ratio for the detrained parcel
real(kind=real_umphys) :: pr(npnts)
                        ! Pressure of the detrained parcel (Pa)
real(kind=real_umphys) :: dqsdth
                        ! Rate of change of qsat with potential temperature
real(kind=real_umphys) :: el
                        ! Latent heat of gas-to-whatever-condenses
                        ! PC2 defn (J/kg)

real(kind=real_umphys), parameter :: a_smth = 0.5
                                ! The weighting between at k and k-1 for
                                ! fdet_opt=2

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='THETAR_6A'


!----------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (fdet_opt == 3) then
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    ! Set the env. properties to the weighted average of k and k+1
    thve(i)   = a_smth *(thek(i)  *(1.0+c_virtual*qek(i)  - watldek(i)))       &
            +(1-a_smth)*(thekp1(i)*(1.0+c_virtual*qekp1(i)- watldekp1(i)))
    exr(i)    = a_smth*exk(i)     +(1.0-a_smth)*exkp1(i)
    ! Set the water loading of the detrained parcel to the weighted
    ! average of k and k+1
    watldr(i) = a_smth*watldpk(i) +(1.0-a_smth)*watldpkp1(i)

    if ( bgmk(i) ) then
      ! If saturated then an iterative calculation is required
      ! Set the first estimate of the detrained p.temp to that of the
      ! environment
      thrk(i)   = a_smth*thek(i)  +(1.0-a_smth)*thekp1(i)
    else
      ! If not saturated the detrained p.temp and humidity can be directly
      ! calculated
      qrk(i)  = 0.5*(qpk(i) + a_smth * qek(i) +(1.0-a_smth)*qekp1(i))
      thrk(i) = thve(i) / (1.0 + c_virtual*qrk(i) - watldr(i))
    end if
    !update the parcel's temperature
    trk(i)  = thrk(i) * exk(i)
  end do

!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    ! calculate pressure at forced detrainment to be consistent with exner
    pr(i)     = pref*exr(i)**recip_kappa
  end do

  call qsat(qsrk, trk, pr, npnts, idx, ni)

  do j=1,3  !Three iterations should be sufficient
!DIR$ IVDEP
    do m=1, ni
      i = idx(m)
      if ( bgmk(i) ) then !if saturated
        if ( bwk(i) ) then
          el=lc
        else
          el=ls
        end if

        !Calculate the gradient of qsat. nb qsrk is saturated humidity at trk
        dqsdth = el * qsrk(i) / ( rv * exr(i) * thrk(i) * thrk(i) )

        !Calculate the next estimate of the parcel's p.temp
        thrk(i) = ( thve(i)                                                    &
                  - thrk(i)*c_virtual*(qsrk(i) - thrk(i)*dqsdth) )             &
                  / ( 1.0 + c_virtual*thrk(i)*dqsdth - watldr(i) )

        !update the parcel's temperature
        trk(i)  = thrk(i) * exr(i)
      end if    !if saturated
    end do      !i

    ! Assume that the detrained parcel is saturated and update its humidity.
    call qsat(qsrk, trk, pr, npnts, idx, ni)

  end do  !j
else ! fdet_opt != 3
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    if ( bgmk(i) ) then
      ! If saturated then an iterative calculation is required
      ! Set the first estimate of the detrained p.temp to that of the
      ! environment
      thrk(i) = thek(i)
    else
      ! If not saturated the detrained p.temp and humidity can be directly
      ! calculated
      if (fdet_opt == 2) then
        thrk(i) = thek(i)*(1.0 + c_virtual*qek(i) - watldek(i))                &
                  /(1.0 + c_virtual*0.5*(qpk(i) + qek(i)) - watldpk(i))
        qrk(i)  = 0.5*(qpk(i) + qek(i))
      else
        thrk(i) = thek(i)*(1.0 + c_virtual*qek(i) - watldek(i))                &
                  /(1.0 + c_virtual*qpk(i) - watldpk(i))
        qrk(i)  = qpk(i)
      end if
    end if
    !update the parcel's temperature
    trk(i)  = thrk(i) * exk(i)
  end do

  call qsat(qsrk, trk, pk, npnts, idx, ni)

  do j=1,3  !Three iterations should be sufficient
!DIR$ IVDEP
    do m=1, ni
      i = idx(m)
      if ( bgmk(i) ) then !if saturated
        if ( bwk(i) ) then
          el=lc
        else
          el=ls
        end if

        !Calculate the gradient of qsat. nb qsrk is saturated humidity at trk
        dqsdth = el * qsrk(i) / ( rv * exk(i) * thrk(i) * thrk(i) )

        !Calculate the next estimate of the parcel's p.temp
        thrk(i) = ( thek(i)*(1.0 + c_virtual*qek(i) - watldek(i))              &
                  - thrk(i)*c_virtual*(qsrk(i) - thrk(i)*dqsdth) )             &
                  / ( 1.0 + c_virtual*thrk(i)*dqsdth - watldpk(i) )

        !update the parcel's temperature
        trk(i)  = thrk(i) * exk(i)
      end if    !if saturated
    end do      !i

    call qsat(qsrk, trk, pk, npnts, idx, ni)

  end do  !j
end if

!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  if ( bgmk(i) ) then
    ! If saturated update qrk to saturated humidity at trk
    qrk(i) = qsrk(i)
  end if
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine thetar_6a
end module thetar_6a_mod
