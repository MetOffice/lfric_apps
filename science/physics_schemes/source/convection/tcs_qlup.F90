! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Module to calculate cumulus ensemble liquid water mixing ratio, ql_up
!
module tcs_qlup


use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use um_types, only: real_umphys

implicit none
!
! Description:
! Module to calculate cumulus ensemble liquid water mixing ratio, ql_up
!
! Method:
!   <reference to documentation to go here, once available>
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!

character(len=*), parameter, private :: ModuleName='TCS_QLUP'

contains

subroutine calc_qlup_warm( n_xx, max_cldlev, wthetal_cb, wqt_cb,               &
   dthetal_cb, dqt_cb, dthetal_cld, dqt_cld, cld_in, sim, ql_up)

use tcs_constants,        only:                                                &
   g, repsilon, lc_o_cp, c_virtual
use tcs_class_similarity, only:                                                &
   similarity
use tcs_class_cloud,      only:                                                &
   cloud_input
use tcs_common_warm ,     only:                                                &
   scales

implicit none
!----------------------------------------------------------------
! Subroutine Arguments
!----------------------------------------------------------------
integer, intent(in) ::                                                         &
   n_xx                                                                        &
                       ! No. of congestus convection points
   ,  max_cldlev
                       ! Maximum number of convective cloud levels


real(kind=real_umphys), intent(in) ::                                          &
   dthetal_cb(n_xx)                                                            &
                            ! dthetal across cloud base
   , dqt_cb(n_xx)                                                              &
                            ! dqt across cloud base
   , dthetal_cld(n_xx)                                                         &
                            ! dthetal across cloud
   , dqt_cld(n_xx)                                                             &
                            ! dqt across cloud
   , wthetal_cb(n_xx)                                                          &
                            ! w'thetal' at cloud base
   , wqt_cb(n_xx)
                            ! w'qt'  at cloud base

type(cloud_input), intent(in) :: cld_in ! fields on cloud levels

type(similarity), intent(in) :: sim !Similarity functions

real(kind=real_umphys), intent(out) ::                                         &
   ql_up(n_xx,max_cldlev)    ! cumulus ensemble ql up    (kg/kg)


!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

real(kind=real_umphys) ::                                                      &
    thetal_up(n_xx,max_cldlev)                                                 &
                            ! thetal  up
   , thetav_up(n_xx,max_cldlev)                                                &
                            ! thetav  up
   , qt_up(n_xx,max_cldlev)
                            ! qt  up

real(kind=real_umphys) ::                                                      &
   thetal_star2(n_xx)                                                          &
                            ! thetal star (squared)
   , qt_star2(n_xx)                                                            &
                            ! qt star (squared)
   , thetal_cb_up2(n_xx)                                                       &
                            ! thetal at cloud base (squared)
   , qt_cb_up2(n_xx)        ! qt at cloud base (squared)

real(kind=real_umphys) ::                                                      &
   thetav_cld(n_xx,max_cldlev)    ! thetav on cloud levels


! temporary variables

real(kind=real_umphys) ::                                                      &
   qt_up2                                                                      &
   , thetal_up2                                                                &
   , virt_term(n_xx)                                                           &
                            ! virtual_term
   , d_term
                            ! denominator term

!-------------------------
! Loop counters
!-------------------------
integer :: i,k

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CALC_QLUP_WARM'

!-----------------------------------------------------------------------
! 1.0 Initialise arrays
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do k = 1,max_cldlev
  do i = 1,n_xx

    ql_up(i,k)    = 0.0

  end do
end do

!-----------------------------------------------------------------------
! 2.0 calculate terms independent of level
!-----------------------------------------------------------------------

do i = 1,n_xx
  !
  ! thetal_star & qt_star  squared
  !
  thetal_star2(i) = scales%mb_new_o_wsc(i)*dthetal_cld(i)*dthetal_cld(i)

  qt_star2(i)     = scales%mb_new_o_wsc(i)*dqt_cld(i)*dqt_cld(i)


  ! upward cloud base fluxes

  qt_cb_up2(i)     = -0.175 * wqt_cb(i) * dqt_cb(i) / scales%mb_new(i)

  thetal_cb_up2(i) = -0.175*wthetal_cb(i)*dthetal_cb(i)/ scales%mb_new(i)

  !
  ! virtual term
  !
  virt_term(i) = sqrt(scales%wstar_up(i)/scales%mb_new(i))                     &
     *scales%wstar_up(i)*scales%wstar_up(i)/(g*scales%zcld(i))
end do


!-----------------------------------------------------------------------
!   Loop over levels - calculating various terms requried for ql
!-----------------------------------------------------------------------
! Note at present doing caluations over a lot of points not in cloud
!-----------------------------------------------------------------------

do k = 1,max_cldlev
  do i = 1,n_xx
    !
    ! Calculate virtual potential temperature in cloud
    !
    thetav_cld(i,k) = cld_in%theta(i,k)*                                       &
       (1.0+cld_in%q_mix(i,k)/repsilon)/(1.0+cld_in%q_mix(i,k))

    !
    ! Ensemble thetal up
    !
    thetal_up2     = thetal_cb_up2(i)*sim%f0_func(i,k)                         &
       + thetal_star2(i)*sim%f1_func(i,k)

    ! Take negative root as cumulus updraughts have lower thetal than the
    ! environment as cloud liquid water forms

    ! thetal_up2 should be positive but as it is the result of combining
    ! 2 terms one of which can be negative the result can be negative.

    if (thetal_up2 <  0.0) then
      thetal_up(i,k) = 0.0
    else
      thetal_up(i,k) = -1.0*sqrt(thetal_up2)
    end if
    !
    ! Ensemble qt up
    !
    qt_up2     = qt_cb_up2(i)*sim%f0_func(i,k)                                 &
       + qt_star2(i)*sim%f1_func(i,k)

    ! Take positive root as cumulus updraughts are greater than mean
    ! environment
    ! qt_up2 should be positive but as it is the result of combining
    ! 2 terms one of which can be negative the result can be negative.

    if (qt_up2 <  0.0) then
      qt_up(i,k) = 0.0
    else
      qt_up(i,k) = sqrt(qt_up2)
    end if

    !
    ! Ensemble virtual potential temperature up
    !
    thetav_up(i,k) = thetav_cld(i,k)*virt_term(i)                              &
       *sim%ftheta_func(i,k)

    !
    ! original calculation of ql_up using thetav_up & qt_up
    !

    d_term = lc_o_cp/cld_in%exner_theta(i,k) -                                 &
                          cld_in%theta(i,k)/repsilon

    ql_up(i,k)=(thetav_up(i,k)-thetal_up(i,k)                                  &
       -c_virtual*cld_in%theta(i,k)*qt_up(i,k))                                &
       /d_term

  end do
end do
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine calc_qlup_warm

end module tcs_qlup
