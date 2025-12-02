! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Lifts Parcel from layer k to k+1

module lift_par_6a_mod

use um_types, only: real_umphys

implicit none

!
! Description:
!   Lifts Parcel from layer k to k+1 taking account of
!   entrainment, detrainment, phase changes and condensation
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


character(len=*), parameter, private :: ModuleName='LIFT_PAR_6A_MOD'

contains

subroutine lift_par_6a (k, npnts, n_wtrac, thek, thekp1,                       &
                        qek, qekp1, qclek, qcfek,                              &
                        qclekp1, qcfekp1,                                      &
                        pkp1, exkp1,                                           &
                        thpk, qpk, qclpk, qcfpk,                               &
                        ekp14, ekp34,                                          &
                        l_q_interact, bwkp1, wtrac_e,                          &
                        !Inout
                        wtrac_p,                                               &
                        !Out
                        bgmkp1, thpkp1, qpkp1,                                 &
                        qclpkp1, qcfpkp1,                                      &
                        Qlkp1, Qfkp1, Frezkp1,                                 &
                        !Indirect indexing
                        idx,ni)

use water_constants_mod, only: lc
use cv_derived_constants_mod, only: ls, lfrcp

use planet_constants_mod, only: cp, rv
use qsat_mod, only: qsat

use wtrac_conv_mod,               only: l_wtrac_conv, conv_e_wtrac_type,       &
                                        conv_p_wtrac_type
use wtrac_conv_store_mod,         only: conv_old_wtrac_type,                   &
                                        wtrac_alloc_conv_store1,               &
                                        wtrac_dealloc_conv_store1
use lift_par_phase_chg_wtrac_mod, only: lift_par_phase_chg_wtrac

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
! Arguments with intent(in):
integer, intent(in) :: k        ! Current vertical level
integer, intent(in) :: npnts    ! Number of points
integer, intent(in) :: n_wtrac  ! Number of water tracers

! Properties of the cloud environment
real(kind=real_umphys),intent(in) :: thek(npnts)
                                    ! Env. pot. temperature in layer k (K)
real(kind=real_umphys),intent(in) :: thekp1(npnts)
                                    ! Env. pot. temperature in layer k+1 (K)
real(kind=real_umphys),intent(in) :: qek(npnts)
                                    ! Env. mixing ratio of in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qekp1(npnts)
                                    ! Env. mixing ratio of in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: qclek(npnts)
                                    ! Env. liquid condensate mixing ratio
                                    ! in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qcfek(npnts)
                                    ! Env. frozen condensate mixing ratio
                                    ! in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qclekp1(npnts)
                                    ! Env. liquid condensate mixing ratio
                                    ! in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: qcfekp1(npnts)
                                    ! Env. frozen condensate mixing ratio
                                    ! in layer k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: pkp1(npnts)
                                    ! Pressure at level k+1 (Pa)
real(kind=real_umphys),intent(in) :: exkp1(npnts)
                                    ! Exner ratio at mid-point of layer k+1

! Properties of the parcel at layer k
real(kind=real_umphys),intent(in) :: thpk(npnts)
                                    ! Par. pot. temperature in layer k (K)
real(kind=real_umphys),intent(in) :: qpk(npnts)
                                    ! Par. mixing ratio of in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qclpk(npnts)
                                    ! Par. liquid condensate mixing ratio
                                    ! in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qcfpk(npnts)
                                    ! Par. frozen condensate mixing ratio
                                    ! in layer k (kg/kg)

!Entrainment and detrainment rates
real(kind=real_umphys),intent(in) :: ekp14(npnts)
                                    ! entrainment coefficient at level k+1/4
                                    ! multiplied by appropriate layer thickness
real(kind=real_umphys),intent(in) :: ekp34(npnts)
                                    ! entrainment coefficient at level k+3/4
                                    ! multiplied by appropriate layer thickness

logical,intent(in) :: l_q_interact  ! True if PC2 is switched on

! Array  arguments with intent(in):
logical,intent(in) :: bwkp1(npnts)  ! Mask for whether condensate is liquid
                                    ! in layer k+1

type(conv_e_wtrac_type), intent(in) :: wtrac_e(n_wtrac)
                                     ! Structure containing environment water
                                     ! tracer arrays

type(conv_p_wtrac_type),intent(in out) :: wtrac_p(n_wtrac)
                                    ! Structure containing parcel
                                    ! water tracer fields

! Array  arguments with intent(out):
! Properties of the parcel at layer k+1 after entrainment, phase changes
! and condensation
logical,intent(out) :: bgmkp1(npnts)! Is Parcel saturated in layer k+1?

real(kind=real_umphys),intent(out) :: thpkp1(npnts)
                                    ! Par. pot. temperature in layer k+1 (K)
real(kind=real_umphys),intent(out) :: qpkp1(npnts)
                                    ! Par. mixing ratio of in layer k+1 (kg/kg)
real(kind=real_umphys),intent(out) :: qclpkp1(npnts)
                                    ! Par. liquid condensate mixing ratio
                                    ! in layer k+1 (kg/kg)
real(kind=real_umphys),intent(out) :: qcfpkp1(npnts)
                                    ! Par. frozen condensate mixing ratio
                                    ! in layer k+1 (kg/kg)
! Amount of moist processes
real(kind=real_umphys),intent(out) :: Qlkp1(npnts)
                                    ! Amount of condensation to liquid water
                                    ! multiplied by the layer thickness (kg/kg)
real(kind=real_umphys),intent(out) :: Qfkp1(npnts)
                                    ! Amount of deposition to ice water
                                    ! multiplied by the layer thickness (kg/kg)
real(kind=real_umphys),intent(out) :: Frezkp1(npnts)
                                    ! Amount of freezing multiplied
                                    ! by the layer thickness (kg/kg)

! Indices of work points
integer, intent(in)  :: idx(npnts)
integer, intent(in)  :: ni


!----------------------------------------------------------------------
! Local variables
!----------------------------------------------------------------------
integer :: i, j, m, i_wt ! Loop counters

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
real(kind=real_umphys) :: Factor(npnts)    ! Factor used in update calculation

real(kind=real_umphys) :: lbycpexner    !L/(cp*exner)
! Water tracer structure to store water values prior to phase change
type(conv_old_wtrac_type) :: wtrac_conv_old

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LIFT_PAR_6A'

!----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
!   Initial 'dry' ascent
! ----------------------------------------------------------------------
!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  Factor(i)  = 1.0/((1.0+ekp14(i))*(1.0+ekp34(i)))
  thpkp1(i)  = ( thpk(i)                                                       &
               + ekp14(i)*thek(i)                                              &
               + (1.0+ekp14(i))*ekp34(i)*thekp1(i) )                           &
               * Factor(i)

  qpkp1(i)   = ( qpk(i)                                                        &
               + ekp14(i)*qek(i)                                               &
               + (1.0+ekp14(i))*ekp34(i)*qekp1(i) )                            &
               * Factor(i)
end do   ! i

if (l_wtrac_conv) then
  ! Initial dry ascent for water tracer
  do i_wt= 1, n_wtrac
!DIR$ IVDEP
    do m=1, ni
      i = idx(m)
      wtrac_p(i_wt)%q(i,k+1) = ( wtrac_p(i_wt)%q(i,k)                          &
                          + ekp14(i)*wtrac_e(i_wt)%q(i,k)                      &
                          + (1.0+ekp14(i))* ekp34(i)*wtrac_e(i_wt)%q(i,k+1) )  &
                          * Factor(i)
    end do   ! m
  end do     ! i_wt
end if ! l_wtrac_conv

if (l_q_interact) then
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    !PC2 so entrainment from the environment
    qclpkp1(i) = ( qclpk(i)                                                    &
               + ekp14(i)*qclek(i)                                             &
               + (1.0+ekp14(i))*ekp34(i)*qclekp1(i) )                          &
               * Factor(i)

    qcfpkp1(i) = ( qcfpk(i)                                                    &
               + ekp14(i)*qcfek(i)                                             &
               + (1.0+ekp14(i))*ekp34(i)*qcfekp1(i) )                          &
               * Factor(i)
  end do  ! i

  if (l_wtrac_conv) then
    do i_wt = 1, n_wtrac
!DIR$ IVDEP
      do m=1, ni
        i = idx(m)
        !PC2 so entrainment from the environment
        wtrac_p(i_wt)%qcl(i,k+1) = ( wtrac_p(i_wt)%qcl(i,k)                    &
                         + ekp14(i)*wtrac_e(i_wt)%qcl(i,k)                     &
                         + (1.0+ekp14(i))*ekp34(i)*wtrac_e(i_wt)%qcl(i,k+1) )  &
                         * Factor(i)

        wtrac_p(i_wt)%qcf(i,k+1) = ( wtrac_p(i_wt)%qcf(i,k)                    &
                         + ekp14(i)*wtrac_e(i_wt)%qcf(i,k)                     &
                         + (1.0+ekp14(i))*ekp34(i)*wtrac_e(i_wt)%qcf(i,k+1) )  &
                         * Factor(i)
      end do  ! m
    end do    ! i_wt
  end if ! l_wtrac_conv

else
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)
    !Not PC2 so no entrainment from the environment
    qclpkp1(i) =  qclpk(i) * Factor(i)

    qcfpkp1(i) =  qcfpk(i) * Factor(i)
  end do  ! i

  ! No water tracer code here as only available for PC2 runs

end if ! l_q_interact

if (l_wtrac_conv) then
  ! Set up arrays and store current values of water for water tracer use
  call wtrac_alloc_conv_store1(npnts, ni, idx, qpkp1, qclpkp1, qcfpkp1,        &
                               wtrac_conv_old)
end if

! ----------------------------------------------------------------------
!       Currently mixed phase parcel is forbidden. Melt or freeze the
!       entrained layer cloud and adjust parcel temperature accordingly.
! ----------------------------------------------------------------------
!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  if (bwkp1(i)) then
    Frezkp1(i) = -qcfpkp1(i)
    qcfpkp1(i) = 0.0
    qclpkp1(i) = qclpkp1(i) - Frezkp1(i)
  else
    Frezkp1(i) = qclpkp1(i)
    qclpkp1(i) = 0.0
    qcfpkp1(i) = qcfpkp1(i) + Frezkp1(i)
  end if
  thpkp1(i)  = thpkp1(i) + Frezkp1(i) * lfrcp / exkp1(i)

end do   ! i

if (l_wtrac_conv) then
  ! Update water tracers for melting/freezing and calculate water tracer
  ! equivalent of Frezkp1
  call lift_par_phase_chg_wtrac(k, npnts, n_wtrac, ni, idx, qclpkp1, qcfpkp1,  &
                                Frezkp1, 'Frezkp1',  wtrac_conv_old%qclpkp1,   &
                                wtrac_conv_old%qcfpkp1, wtrac_p)
end if ! l_wtrac_conv

!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  tpkp1(i)     = thpkp1(i) * exkp1(i)
  ! Save dry ascent values
  thpkp1dry(i) = thpkp1(i)
  qpkp1dry(i)  = qpkp1(i)
end do

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

  end do  !npnts

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

  ! Ql and Qf can be negative if evaporation rather than condensation.
  ! In this case, limit evaporation of condensate to the amount of
  ! available condensate.
  Qlkp1(i)   = max(Qlkp1(i), -qclpkp1(i))
  Qfkp1(i)   = max(Qfkp1(i), -qcfpkp1(i))

  ! To prevent evaporation comment out previous two lines and uncomment the
  ! following two lines
  !  Qlkp1(i)   = max(Qlkp1(i), 0.0)
  !  Qfkp1(i)   = max(Qfkp1(i), 0.0)

  thpkp1(i)  = thpkp1dry(i) + (lc*Qlkp1(i) + ls*Qfkp1(i))                      &
               / (cp * exkp1(i))
  tpkp1(i)   = thpkp1(i) * exkp1(i)

  qclpkp1(i) = qclpkp1(i) + Qlkp1(i)
  qcfpkp1(i) = qcfpkp1(i) + Qfkp1(i)

  qpkp1(i)   = qpkp1dry(i) - Qlkp1(i) - Qfkp1(i)

end do  !points

if (l_wtrac_conv) then

  ! Update water tracers for phase changes and calculate the water tracer
  ! equivalents to Qlkp1 and Qfkpq

  ! Condensation/evaporation
  call lift_par_phase_chg_wtrac(k, npnts, n_wtrac, ni, idx, qpkp1, qclpkp1,    &
                                Qlkp1, 'Qlkp1', wtrac_conv_old%qpkp1,          &
                                wtrac_conv_old%qclpkp1, wtrac_p)

  ! Deposition/sublimation
  call lift_par_phase_chg_wtrac(k, npnts, n_wtrac, ni, idx, qpkp1, qcfpkp1,    &
                                Qfkp1,  'Qfkp1', wtrac_conv_old%qpkp1,         &
                                wtrac_conv_old%qcfpkp1, wtrac_p)

  ! Deallocate working arrays
  call wtrac_dealloc_conv_store1(wtrac_conv_old)

end if

! Calculate qsat at the parcel's temp after condensation
! or evaporation
call qsat(qspkp1,tpkp1,pkp1,npnts,idx,ni)

!DIR$ IVDEP
do m=1, ni
  i = idx(m)
  bgmkp1(i) = ( qpkp1(i) > 0.9999 * qspkp1(i) )
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lift_par_6a
end module lift_par_6a_mod
