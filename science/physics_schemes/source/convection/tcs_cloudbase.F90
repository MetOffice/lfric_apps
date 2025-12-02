! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module to calculate the fluxes through cloud base.
!
module tcs_cloudbase

  ! Working note: need to sort out level and jump calculations done
  ! here with those done in the main routine


use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use um_types, only: real_umphys

implicit none
!
! Description:
! This module defines the tcs warm rain "cloud_input" derived type and
! subroutines for allocating and deallocating an instance of this class.
!
!
! Method:
!   <reference to documentation to go here, once available>
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.2 programming standards.


character(len=*), parameter, private :: ModuleName='TCS_CLOUDBASE'

contains

subroutine calc_cloudbase( n_xx, ntra, l_tracer,                               &
     wthetav_surf, dtracer_cb,                                                 &
     dthetav_cb,dthetal_cb, drt_cb,                                            &
     wtheta_plus_cb, wthetal_cb,wq_plus_cb,                                    &
     wqt_cb,wql_cb, wthetav_cb, wthetav_cb_m2, wtracer_cb)
!-----------------------------------------------------------------------
!
! Description:
!   Calculate the fluxes at cloud base
!-----------------------------------------------------------------------

use tcs_parameters_warm,      only:                                            &
   wql_cb_a1, jump_cb_a1, wthv_cb_a1, sat_a1
use tcs_constants,            only:                                            &
   lc_o_cp, c_virtual, lc, rv
use tcs_common_warm ,         only:                                            &
   scales, cb, cb_m1

implicit none
!------------------------------------------------------------------
! Subroutine Arguments
!------------------------------------------------------------------
!
! Arguments with intent in:
!

integer, intent(in) ::                                                         &
   n_xx                                                                        &
                   ! No. of congestus convection points
,  ntra            ! No. of tracers

logical, intent(in) ::                                                         &
   l_tracer        ! true for tracers

real(kind=real_umphys), intent(in) ::                                          &
  wthetav_surf(n_xx)                                                           &
                         ! w'thetav' at surface
, dtracer_cb(n_xx,ntra)                                                        &
                         ! change in tracers across cloud base
                         ! (kg/kg)
, dthetav_cb(n_xx)
                         ! dthetav across cloud base
!
! Arguments with intent INOUT:
!
!          none

!
! Arguments with intent out:
!

real(kind=real_umphys), intent(out) ::                                         &
  wthetal_cb(n_xx)                                                             &
                            !  w'thetal' at cloud base
, wqt_cb(n_xx)                                                                 &
                            !  w'qt'  at cloud base
, wql_cb(n_xx)                                                                 &
                            !  w'ql' at cloud base
, wq_plus_cb(n_xx)                                                             &
                            !  w'q' at cloud base
, wtheta_plus_cb(n_xx)                                                         &
                            !  w'theta' at cloud base
, wthetav_cb(n_xx)                                                             &
                            !  w'thetav' at cloud base
, wthetav_cb_m2(n_xx)                                                          &
                            !  w'thetav' at cloud base using
                            !  method 2
, dthetal_cb(n_xx)                                                             &
                         ! dthetal across cloud base
, drt_cb(n_xx)                                                                 &
                         ! drt across cloud base
, wtracer_cb(n_xx,ntra)  ! w'tracer' at cloud base

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
real(kind=real_umphys) ::                                                      &
  T_plus(n_xx)
                         ! T at + side of cloud base

! temporary variables
real(kind=real_umphys) ::                                                      &
  term_a(n_xx)                                                                 &
, term_c(n_xx)                                                                 &
, term_d(n_xx)                                                                 &
, r_rsat_term(n_xx)                                                            &
                   ! r-rsat   term
, drsatdt_topt(n_xx)
                   ! drsat/dT at the top of the transition region

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CALC_CLOUDBASE'


!-------------------------
! Loop counters
!-------------------------
integer :: ktra

!-----------------------------------------------------------------
! Calculate Some derived fields
!-----------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

t_plus(:)  = cb%theta(:)*cb%exner_theta(:)

!-----------------------------------------------------------------
! The convective parametrization scheme assumes that the
! enviromental air has qcl=0 and qcf=0. i.e thetal=theta and rt=r
! (this may not be the case in the model).
!-----------------------------------------------------------------
dthetal_cb(:) = cb%theta(:) - cb_m1%theta(:)

drt_cb(:)     = cb%q_mix(:) - cb_m1%q_mix(:)

!-----------------------------------------------------------------
! assume w'ql'+ is value at cloud base
!-----------------------------------------------------------------

wql_cb(:) = wql_cb_a1*wthetav_surf(:)                                          &
     /(lc_o_cp/cb%exner_theta(:)-(1.0+c_virtual)*cb%theta(:))

!-----------------------------------------------------------------
! Note that this is the value of wthetav on the lower
! interface of cloudbase and should be negative.  Here the
! value should be the same as that for wthetavl since there is
! no condensed water.
!-----------------------------------------------------------------
wthetav_cb(:) = -jump_cb_a1*scales%mb(:)*dthetav_cb(:) ! = wthv_minus

wthetav_cb_m2(:) = -wthv_cb_a1*wthetav_surf(:)  ! = wthv_minus

!-----------------------------------------------------------------
! expression for dqsat/dT
!-----------------------------------------------------------------
drsatdt_topt(:) = Lc *cb%qse(:) /(Rv * t_plus(:)*t_plus(:))

r_rsat_term(:) = cb%q_mix(:) - cb%qse(:)

!-----------------------------------------------------------------
! Invert saturation condition to get wqt and wthetal
!-----------------------------------------------------------------
term_a(:) = wthetav_cb(:)

term_c(:) = wql_cb(:)*(1.0+lc_o_cp*drsatdt_topt)/sat_a1

term_d(:) = scales%mb(:)*r_rsat_term

wqt_cb(:) = (term_c(:) - term_d(:)                                             &
     + cb%exner_theta(:)*drsatdt_topt(:)*term_a(:) )/                          &
     (1.0+c_virtual*cb%exner_theta(:)*cb_m1%theta(:)*drsatdt_topt)

wthetal_cb(:) = term_a(:) - c_virtual*cb_m1%theta(:)*wqt_cb(:)

wq_plus_cb(:) = wqt_cb(:) - wql_cb(:)

wtheta_plus_cb(:) = wthetal_cb(:)                                              &
                     + lc_o_cp*wql_cb(:)/cb_m1%exner_rho(:)


!      ! Working note: The moist static energy flux through cloudbase
!      ! needs some further investigation from LEM results.
!      ! Working note: not used in this version - fqt_temp was
!      ! passed through the tcs_temp_module

      ! Note expect wthetal & wtheta to be negative at cloud base and wqt
      ! wq and wql to be positive.

      !-----------------------------------------------------------------------
      ! tracers - assumed to behave as other conserved variables

if (l_tracer) then
  do ktra = 1,ntra
    wtracer_cb(:,ktra) = -jump_cb_a1*scales%mb(:)*dtracer_cb(:,ktra)
  end do
end if
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine calc_cloudbase

end module tcs_cloudbase
