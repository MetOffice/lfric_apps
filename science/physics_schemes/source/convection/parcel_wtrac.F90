! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module parcel_wtrac_mod

use um_types, only: real_umphys

implicit none

! Description:
!   Calculate the water tracer content of the parcel at level k+1.
!
! Method:
!   Using the same method as used to update the prognostic water variables in
!   subroutine parcel.  (Assuming the use of PC2 here so l_q_interact = T.)
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private :: ModuleName = 'PARCEL_WTRAC_MOD'

contains

! Subroutine Interface:
subroutine parcel_wtrac(k, npnts, ni, n_wtrac, idx, tiny_value_condensate,     &
                        deltak, ekp14, ekp34, qpk, qrk, qclpkp1, qcfpkp1,      &
                        wtrac_e, bterm, wtrac_p, qrk_wtrac)

use wtrac_conv_mod,       only: conv_e_wtrac_type, conv_p_wtrac_type
use wtrac_calc_ratio_mod, only: wtrac_calc_ratio_fn

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
integer,intent(in) :: k          ! Present model layer
integer,intent(in) :: npnts      ! No. of points
integer,intent(in) :: ni         ! No. of work points
integer,intent(in) :: n_wtrac    ! Number of water tracers

integer,intent(in) :: idx(npnts) ! Index of work points

real(kind=real_umphys),intent(in) :: tiny_value_condensate
                                 ! Condensate has been removed if below this
                                 ! value
real(kind=real_umphys),intent(in) :: deltak(npnts)
                                 ! Parcel forced detrainment rate in
                                 ! layer k multiplied by layer thickness
real(kind=real_umphys),intent(in) :: ekp14(npnts)
                                 ! Entrainment coefficient at level k+1/4
                                 ! multiplied by appropriate layer thickness
real(kind=real_umphys),intent(in) :: ekp34(npnts)
                                 ! Entrainment coefficient at level k+3/4
                                 ! multiplied by appropriate layer thickness

real(kind=real_umphys),intent(in) :: qpk(npnts)
                                 ! Parcel specific humidity in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qrk(npnts)
                                 ! Specific humidity of forced detrained
                                 ! parcel in layer k (kg/kg)
real(kind=real_umphys),intent(in) :: qclpkp1(npnts)
                                 ! Parcel liquid content at k+1 (kg/kg)
real(kind=real_umphys),intent(in) :: qcfpkp1(npnts)
                                 ! Parcel ice content at k+1 (kg/kg)

type(conv_e_wtrac_type), intent(in) :: wtrac_e(n_wtrac)
                                 ! Structure containing environment
                                 ! water tracer fields

logical,intent(in) :: bterm(npnts) ! True if parcel is terminating

type(conv_p_wtrac_type), intent(in out) :: wtrac_p(n_wtrac)
                                 ! Structure containing parcel
                                 ! water tracer fields

real(kind=real_umphys),intent(in out) :: qrk_wtrac(npnts,n_wtrac)
                                 ! Water tracer content of detraining vapour

! Local variables
integer :: i, m, i_wt    ! Loop counters

real(kind=real_umphys) :: ratio_qrk     ! Ratio used in calculating qrk_wtrac
real(kind=real_umphys) :: Factor1       ! Factor used in update calculation
real(kind=real_umphys) :: Factor2       ! Factor used in update calculation

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='PARCEL_WTRAC'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do m = 1, ni
  i = idx(m)
  if (bterm(i)) then
    ! Parcel is terminating

    do i_wt = 1, n_wtrac
      wtrac_p(i_wt)%q(i,k+1)   = 0.0
      wtrac_p(i_wt)%qcl(i,k+1) = 0.0
      wtrac_p(i_wt)%qcf(i,k+1) = 0.0
      qrk_wtrac(i,i_wt)        = wtrac_p(i_wt)%q(i,k)
    end do

  else
    ! Parcel is not terminating. Calculate parcel properties at level k+1

    do i_wt = 1, n_wtrac
      ! Calculate water tracer content of detraining vapour (based on ratio
      ! of water tracer vapour to water vapour in parcel)

      ratio_qrk         = wtrac_calc_ratio_fn                                  &
                         (i_wt,wtrac_p(i_wt)%q(i,k),qpk(i))
      qrk_wtrac(i,i_wt) = ratio_qrk * qrk(i)

      ! Calculate parcel content at k+1

      Factor1 = 1.0/((1.0+ekp14(i))*(1.0+ekp34(i)))
      Factor2 = Factor1/(1.0-deltak(i))

      wtrac_p(i_wt)%q(i,k+1) = ( wtrac_p(i_wt)%q(i,k)                          &
             - deltak(i)*qrk_wtrac(i,i_wt)                                     &
             + (1.0-deltak(i))*ekp14(i)*wtrac_e(i_wt)%q(i,k)                   &
             + (1.0-deltak(i))                                                 &
             * (1.0+ekp14(i))*ekp34(i)*wtrac_e(i_wt)%q(i,k+1) )                &
             * Factor2                                                         &
             - wtrac_p(i_wt)%Qlkp1(i) - wtrac_p(i_wt)%Qfkp1(i)

      !PC2 so entrainment from the environment
      wtrac_p(i_wt)%qcl(i,k+1) = ( wtrac_p(i_wt)%qcl(i,k)                      &
                          + ekp14(i)*wtrac_e(i_wt)%qcl(i,k)                    &
                          + (1.0+ekp14(i))                                     &
                          * ekp34(i)*wtrac_e(i_wt)%qcl(i,k+1) )                &
                          * Factor1                                            &
                          + wtrac_p(i_wt)%Qlkp1(i) - wtrac_p(i_wt)%Frezkp1(i)

      wtrac_p(i_wt)%qcf(i,k+1) = ( wtrac_p(i_wt)%qcf(i,k)                      &
                          + ekp14(i)*wtrac_e(i_wt)%qcf(i,k)                    &
                          + (1.0+ekp14(i))                                     &
                          * ekp34(i)*wtrac_e(i_wt)%qcf(i,k+1) )                &
                          * Factor1                                            &
                          + wtrac_p(i_wt)%Qfkp1(i) + wtrac_p(i_wt)%Frezkp1(i)

      ! If any tiny amounts of qcl and qcf have been removed, ensure
      ! water tracers are also removed
      if (abs(qclpkp1(i)) < tiny_value_condensate) then
        wtrac_p(i_wt)%qcl(i,k+1) = 0.0
      end if
      if (abs(qcfpkp1(i)) < tiny_value_condensate) then
        wtrac_p(i_wt)%qcf(i,k+1) = 0.0
      end if

    end do ! i_wt
  end if   ! bterm
end do     ! i


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine parcel_wtrac

end module parcel_wtrac_mod

