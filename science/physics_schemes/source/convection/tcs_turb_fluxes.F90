! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Module to calclulate the turbulent fluxes w'theta' and w'q'
! on in-cloud levels
!
module tcs_turb_fluxes


use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use um_types, only: real_umphys

implicit none
!
! Description:
! Module to calclulate the turbulent fluxes w'theta' and w'q'
! on in-cloud levels
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


character(len=*), parameter, private :: ModuleName='TCS_TURB_FLUXES'

contains

subroutine calc_turb_fluxes(n_xx,ntra, max_cldlev, nlev                        &
   ,                      max_cldtrlev, trlev                                  &
   ,                      l_tracer                                             &
   ,                      kterm_thv, kterm_tracer                              &
   ,                      sim                                                  &
   ,                      wql_cld, mf_cld, cld_in                              &
   ,                      dqsatdt_cld                                          &
   ,                      precip_product_cld,wup_h_cld                         &
   ,                      scales,wthetav_cb,wtracer_cb                         &
   ,                      wtheta, wq, wthetavl, wthetav                        &
   ,                      wthetal, wqt, wh, wtracer)

use tcs_parameters_warm,    only:                                              &
   icong_options, wthvl_factor
use tcs_constants,          only:                                              &
    lc_o_cp, c_virtual, g
use tcs_class_similarity,   only:                                              &
    similarity
use tcs_class_scales,       only:                                              &
    scales_conv
use tcs_class_cloud,        only:                                              &
    cloud_input

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
   , ntra                                                                      &
                            ! No. of tracers
   , max_cldlev                                                                &
                            ! Maximum number of convective cloud levels
   , nlev                                                                      &
                            ! Maximum number of convective cloud levels
   , max_cldtrlev                                                              &
                            ! Maximum number of convective cloud levels
   , trlev
! Maximum number of tracer levels

logical, intent(in) ::                                                         &
   l_tracer      ! true - tracers present

real(kind=real_umphys), intent(in) ::                                          &
   kterm_thv(n_xx,nlev)
                            ! thetav gradient (kterm)

type(similarity) :: sim   ! similarity functions

real(kind=real_umphys), intent(in) ::                                          &
   wql_cld(n_xx,nlev)                                                          &
                            ! wql on uv cloud lev (similarity func)
   , mf_cld(n_xx,nlev)
! mass_flux on theta levels in cloud

type(cloud_input), intent(in) :: cld_in
! input fields on cloud levels

real(kind=real_umphys), intent(in) ::                                          &
   dqsatdt_cld(n_xx,nlev)                                                      &
                            !dqsat/dT on theta levels in cloud
   , precip_product_cld(n_xx,nlev)                                             &
                            ! precip production rho levels
   , wup_h_cld (n_xx,nlev)                                                     &
                            ! wup on rho levels
   , kterm_tracer(n_xx,nlev,ntra)
! tracer gradient (kterm)
!                                   dimensioned on nlev not trlev

type(scales_conv), intent(in) :: scales

real(kind=real_umphys), intent(in) ::                                          &
   wthetav_cb(n_xx)                                                            &
                            ! w'thetav' at cloud base (K m/s)
   , wtracer_cb(n_xx,ntra)  ! w'tracer' at cloud base (kg/kg m/s)
!
! Arguments with intent INOUT:
!
!     None
!
! Arguments with intent out:
!
real(kind=real_umphys), intent(out) ::                                         &
   wtheta(n_xx,nlev)                                                           &
                            ! wtheta flux on cloud levels  (K m/s)
   ,wq(n_xx,nlev)                                                              &
                            ! wq flux on cloud levels      (kg/kg m/s)
   ,wthetavl(n_xx,nlev)                                                        &
                            ! wthetavl flux on cloud levels  (K m/s)
   ,wthetav(n_xx,nlev)                                                         &
                            ! wthetav flux on cloud levels  (K m/s)
   ,wthetal(n_xx,nlev)                                                         &
                            ! wthetal flux on cloud levels  (K m/s)
   ,wqt(n_xx,nlev)                                                             &
                            ! wqt flux on cloud levels      (kg/kg m/s)
   ,wh(n_xx,nlev)                                                              &
                            ! wh flux on cloud levels      (K m/s)
   ,wtracer(n_xx,trlev,ntra) ! wtracer flux on cloud levels(kg/kg m/s)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------
real(kind=real_umphys) ::                                                      &
   wtheta_vl_cb(n_xx)
! wtheta_vl flux at cloud base.

! temporary stores
real(kind=real_umphys) ::                                                      &
   fngterm2,  A_term                                                           &
   , div_term

!-------------------------
! Loop counters
!-------------------------
integer :: i,k,ktra

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CALC_TURB_FLUXES'

!-----------------------------------------------------------------------
! 1.0 Initialise arrays - set all functions to zero
!-----------------------------------------------------------------------
! Initialise dhdz

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do k = 1,max_cldlev
  do i = 1,n_xx
    wq(i,k) = 0.0
    wqt(i,k) = 0.0
    wh(i,k) = 0.0
    wtheta(i,k) = 0.0
    wthetavl(i,k) = 0.0
    wthetav(i,k) = 0.0
    wthetal(i,k) = 0.0
  end do
end do


!-----------------------------------------------------------------------
! 2.0 Level loop to calculate in cloud values of fluxes on uv levels
!-----------------------------------------------------------------------
!  Numbering of levels
!  level 1 is cloud base for uv - in cloud levels
!
!-----------------------------------------------------------------------
do i = 1,n_xx

  wtheta_vl_cb(i) = wthetav_cb(i)

end do


do k = 1,max_cldlev
  do i = 1,n_xx

    ! thetav flux
    fngterm2 = wtheta_vl_cb(i)*sim%fng_func_rho(i,k)


    ! additional term for precip added to original w'thetav'
    wthetavl(i,k) = -kterm_thv(i,k) + fngterm2                                 &
       + (lc_o_cp/cld_in%exner_theta(i,k)-c_virtual*cld_in%theta(i,k))         &
       *scales%zcld(i)*(wup_h_cld(i,k)/scales%wstar_up(i))                     &
       *precip_product_cld(i,k)/cld_in%rho(i,k)

  end do
end do
! Here we've split up the loop so that wthetavl can be renormalized if
! required


if (icong_options == 7) then
  ! Rescale wthetavl using TKE scaling
  do i = 1,n_xx
    wthetavl(i,1:max_cldlev)=wthetavl(i,1:max_cldlev)                          &
       *max_cldlev/sum(wthetavl(i,1:max_cldlev))                               &
       *sign(1.0,wthetavl(i,1:max_cldlev))                                     &
       *wthvl_factor*(cld_in%theta(i,1:max_cldlev)/g)                          &
       *(scales%mb(i)/scales%wstar_up(i))**(0.5)                               &
       *scales%wstar_up(i)**3/scales%zcld(i)
  end do
end if

do k = 1,max_cldlev
  do i = 1,n_xx

    ! calculations depend on interpolation for q and qsat
    ! preferred method
    if (sim%gql_func_rho(i,k) /= 0.0) then
      A_term = (1.0+lc_o_cp*dqsatdt_cld(i,k))*wql_cld(i,k)                     &
         /sim%gql_func_rho(i,k) - mf_cld(i,k)                                  &
         *(cld_in%q_mix(i,k) - cld_in%qse(i,k))
    else  ! not in cloud
      a_term=0.0
    end if

    ! Should I use mf on th or uv levels?

    div_term=1.0+c_virtual*cld_in%exner_theta(i,k)*cld_in%theta(i,k)           &
       *dqsatdt_cld(i,k)

    wthetav(i,k)=wthetavl(i,k)                                                 &
       +(lc_o_cp/cld_in%exner_theta(i,k)-c_virtual*cld_in%theta(i,k))          &
       *wql_cld(i,k)


    wthetal(i,k)=(wthetavl(i,k)-c_virtual*cld_in%theta(i,k)*a_term)            &
       /div_term

    wqt(i,k) =( A_term + wthetavl(i,k)*cld_in%exner_theta(i,k)                 &
       *dqsatdt_cld(i,k))/div_term


    wqt(i,k)=max(0.0,wqt(i,k))

  end do
end do


! ----------------------------------------------------------------------
! 3.0 Tracers
! ----------------------------------------------------------------------

if (l_tracer) then

  ! initialise
  do ktra = 1,ntra
    do k = 1,max_cldtrlev
      do i = 1,n_xx
        wtracer(i,k,ktra) = 0.0
      end do
    end do
  end do

  do ktra = 1,ntra
    do k = 1,max_cldtrlev
      do i = 1,n_xx

        ! non-gradient term
        fngterm2 = wtracer_cb(i,ktra)*sim%fng_func_rho(i,k)

        wtracer(i,k,ktra) = -kterm_tracer(i,k,ktra) + fngterm2

      end do
    end do
  end do
end if
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine calc_turb_fluxes

end module tcs_turb_fluxes
