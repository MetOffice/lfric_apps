! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Lifts Undilute Parcel from layer k to k+1

module lift_undil_par_6a_mod

use um_types, only: real_umphys

implicit none

!
! Description:
!   Lifts Parcel from layer k to k+1 taking account of
!   condensation but assuming zero parcel condensate
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


character(len=*), parameter, private :: ModuleName='LIFT_UNDIL_PAR_6A_MOD'

contains

subroutine lift_undil_par_6a (npnts,                                           &
                        pkp1, exkp1,                                           &
                        thpk, qpk,                                             &
                        Frezkp1, bwkp1,                                        &
                        !Out
                        thpkp1, qpkp1,                                         &
                        !Indirect indexing
                        idx,ni)

use water_constants_mod, only: lc
use cv_derived_constants_mod, only: ls, lfrcp
use planet_constants_mod, only: cp, rv

use qsat_mod, only: qsat

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

! Subroutine arguments
! Arguments with intent(in):
integer, intent(in) :: npnts    ! Number of points

! Properties of the cloud environment
real(kind=real_umphys),intent(in) :: pkp1(npnts)
                                    ! Pressure at level k+1 (Pa)
real(kind=real_umphys),intent(in) :: exkp1(npnts)
                                    ! Exner ratio at mid-point of layer k+1

! Properties of the undilute parcel at layer k
real(kind=real_umphys),intent(in) :: thpk(npnts)
                                    ! Par. pot. temperature in layer k (K)
real(kind=real_umphys),intent(in) :: qpk(npnts)
                                    ! Par. mixing ratio of in layer k (kg/kg)

! Process rates at layer k+1 (from dilute ascent)
real(kind=real_umphys),intent(in) :: Frezkp1(npnts)
                                    ! Amount of freezing from liquid
                                    ! to ice water in the parcel (kg/kg)
! Logical Array arguments with intent(in):
logical,intent(in) :: bwkp1(npnts)  ! Mask for whether condensate is liquid
                                    ! in layer k+1

! Array  arguments with intent(INOUT):
! Properties of the parcel at layer k+1 after condensation
real(kind=real_umphys),intent(in out) :: thpkp1(npnts)
                                     ! Par. pot. temperature in layer k+1 (K)
real(kind=real_umphys),intent(in out) :: qpkp1(npnts)
                                     ! Par. mixing ratio of in layer k+1 (kg/kg)

! Indices of work points
integer, intent(in)  :: idx(npnts)
integer, intent(in)  :: ni


!----------------------------------------------------------------------
! Local variables
!----------------------------------------------------------------------
integer :: i, j, m       ! Loop counters

real(kind=real_umphys) :: tpkp1(npnts)     ! Par. temperature in layer k+1 (K)
real(kind=real_umphys) :: thpkp1dry(npnts) ! Par. pot. temperature in layer k+1
                         ! after dry ascent (K)
real(kind=real_umphys) :: qpkp1dry(npnts)  ! Par. mixing ratio of in layer k+1
                         ! after dry ascent (kg/kg)
real(kind=real_umphys) :: qspkp1(npnts)    ! Saturation Mixing ratio of parcel
                         ! after dry ascent (kg/kg)
real(kind=real_umphys) :: dqsdth
                         ! Rate of change of qsat with potential temperature
real(kind=real_umphys) :: el
                         ! Latent heat of gas-to-whatever-condenses PC2 defn
                         ! (J/kg)
                         ! Amount of moist processes
real(kind=real_umphys) :: Qlkp1(npnts)
                         ! Amount of condensation to liquid water
                         ! multiplied by the layer thickness (kg/kg)
real(kind=real_umphys) :: Qfkp1(npnts)     ! Amount of deposition to ice water
                         ! multiplied by the layer thickness (kg/kg)

real(kind=real_umphys) :: lbycpexner    !L/(cp*exner)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LIFT_UNDIL_PAR_6A'

!----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
!   Initial 'dry' ascent
! ----------------------------------------------------------------------
!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  thpkp1(i)  = thpk(i)
  qpkp1(i)   = qpk(i)
end do   ! i


!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  tpkp1(i)     = thpkp1(i) * exkp1(i)
  ! Save dry ascent values
  thpkp1dry(i) = thpkp1(i)
  qpkp1dry(i)  = qpkp1(i)
end do

! ----------------------------------------------------------------------
!  Adjust parcel temperature for freezing (or melting).
!  Assume the same freezing/melting rate as the dilute ascent
! ----------------------------------------------------------------------
!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  thpkp1(i)  = thpkp1(i) + Frezkp1(i) * lfrcp / exkp1(i)
end do   ! i

call qsat(qspkp1,tpkp1,pkp1,npnts,idx,ni)

! ----------------------------------------------------------------------
!       Calculate theta and q if the parcel was brought to saturation
! ----------------------------------------------------------------------

do j=1,3  !Three iterations should be sufficient

!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    if ( bwkp1(i) ) then
      el=lc
    else
      el=ls
    end if

    dqsdth = el * qspkp1(i) / ( rv * exkp1(i) * thpkp1(i) * thpkp1(i) )
    lbycpexner=el/(cp*exkp1(i))

    !Calculate the next estimate of the parcel's p.temp after condensation
    thpkp1(i) = ( thpkp1dry(i) + lbycpexner*(qpkp1dry(i) - qspkp1(i)           &
                 + thpkp1(i)*dqsdth) ) /                                       &
                 (1.0 + lbycpexner*dqsdth)

    !Calculate the next estimate of the parcel's temp after condensation
    tpkp1(i) = thpkp1(i) * exkp1(i)

  end do

  ! Calculate qsat at the next estimate of the parcel's temp after condensation
  call qsat(qspkp1,tpkp1,pkp1,npnts,idx,ni)

end do !End iterations


! ----------------------------------------------------------------------
! Update parcel properties to take account of condensation
! or evaporation
! ----------------------------------------------------------------------

!DIR$ IVDEP
do m=1, ni
  i = idx(m)

  if ( bwkp1(i) ) then
    Qlkp1(i) = qpkp1dry(i) - qspkp1(i)
    Qfkp1(i) = 0.0
  else
    Qfkp1(i) = qpkp1dry(i) - qspkp1(i)
    Qlkp1(i) = 0.0
  end if

  ! Limit to allow condensation only
  Qlkp1(i)   = max(Qlkp1(i), 0.0)
  Qfkp1(i)   = max(Qfkp1(i), 0.0)

  thpkp1(i)  = thpkp1dry(i) + (lc*Qlkp1(i) + ls*Qfkp1(i))                      &
               / (cp * exkp1(i))
  qpkp1(i)   = qpkp1dry(i) - Qlkp1(i) - Qfkp1(i)

end do  !points

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lift_undil_par_6a
end module lift_undil_par_6a_mod
