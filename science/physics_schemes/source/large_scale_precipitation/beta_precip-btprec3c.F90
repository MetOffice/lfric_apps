! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  subroutine BETA_PRECIP --------------------------------------------
!     PURPOSE:
! Process fields of precipitation intensity to give scattering coefft
! in 1/metres.
! Calculated at model level (eg bottom eta level 25m)
! or level within surface layer eg screen ht ( 1.5M )
!  Programming standard: U M Doc. Paper No. 4
!
!  External documentation
!    Forecasting Research Scientific Paper NO.4
!    Diagnosis of visibility in the UK Met Office Mesoscale Model
!    and the use of a visibility analysis to constrain initial
!    conditions.  SP Ballard, BJ Wright, BW Golding    1992
!      NIMROD diagnostic:
!      Wright, B. J., 1997: Improvements to the Nimrod Visibility
!         Analysis/Forecast System. FR-Div. Tech. Rep., No. 217.
!      Wright, B. J., 1997: A New Visibility Analysis/Forecast System
!         for Nimrod. Met. Office FR Tech Rep., No. 222.
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation
module beta_precip_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='BETA_PRECIP_MOD'

contains

subroutine beta_precip                                                         &
           (ls_rain, ls_snow, c_rain, c_snow, qcf, qrain,                      &
                                                      !INPUT
            rho, t, pressure, snownumber, rainnumber,                          &
                                                      !INPUT
            lca,cca,pct,avg,                                                   &
                                                      !INPUT
            p_field,points,k1stpt,                                             &
                                                      !INPUT
            beta_ls_rain, beta_ls_snow,                                        &
                                                      !OUTPUT
            beta_c_rain, beta_c_snow )                !OUTPUT

  ! General atmosphere modules
use planet_constants_mod, only: pref
use water_constants_mod, only: rho_water
use conversions_mod,     only: pi

! Microphysics modules
use mphys_psd_mod,       only: ri, si, cr, dr, x2i, x4i, ci0, di0
use mphys_inputs_mod,    only: x1r, x2r, ai,  bi, l_psd
use mphys_constants_mod, only: x1i, x4r

! Error reporting modules and print manager
use errormessagelength_mod, only: errormessagelength
use ereport_mod,            only: ereport
use umprintmgr,             only: newline

! Dr Hook Modules
use yomhook,             only: lhook, dr_hook
use parkind1,            only: jprb, jpim

! Large scale precipitation modules
use gammaf_mod,          only: gammaf
use lsp_moments_mod,     only: lsp_moments

!For vis calculation using casim
use mphys_inputs_mod,           only: l_casim
use thresholds,                 only: qs_small, qr_small, nr_small, ns_small
use atm_fields_bounds_mod,      only:  tdims, tdims_l
use mphys_parameters,           only: snow_params, rain_params
implicit none
!---------------------------------------------------------------------
! input variables-----------------------------------------------------
!---------------------------------------------------------------------
integer ::                                                                     &
 p_field,                                                                      &
                                      ! in NO. points in field.
 points,                                                                       &
              ! in Number of gridpoints being processed.
 k1stpt
              ! in First gridpoint processed within complete field.

real, intent(in) ::  rainnumber(tdims_l%i_start:tdims_l%i_end,                 &
                                tdims_l%j_start:tdims_l%j_end,                 &
                                tdims_l%k_start:tdims_l%k_end)
! Rain number concentration from Casim
real, intent(in) ::  snownumber(tdims_l%i_start:tdims_l%i_end,                 &
                                tdims_l%j_start:tdims_l%j_end,                 &
                                tdims_l%k_start:tdims_l%k_end)
! Snow number concentration from Casim
real(kind=real_umphys) ::                                                      &
 ls_rain(p_field),                                                             &
                                      ! in Large scale Rain
 ls_snow(p_field),                                                             &
                                      ! in Large scale Snow
 c_rain(p_field),                                                              &
                                      ! in Convective Rain
 c_snow(p_field),                                                              &
                                      ! in Convective Snow
 qcf(p_field),                                                                 &
                                      ! in large-scale snow / kg kg-1
 qrain(p_field),                                                               &
                                      ! in large-scale rain on levels / kg kg-1
 rho(p_field),                                                                 &
                                      ! in Air density / kg m-3
 t(p_field),                                                                   &
                                      ! in Temperature / K
 pressure(p_field),                                                            &
                                      ! in Pressure
 lca(p_field),                                                                 &
                                      ! in Total Layer Cloud.
 cca(p_field)                         ! in Convective Cloud.


logical ::                                                                     &
 pct,                                                                          &
                                      ! in T:Cloud amounts are in %
 avg
                                      ! in T:Precip =local*prob

!---------------------------------------------------------------------
! output variables----------------------------------------------------
!---------------------------------------------------------------------
real(kind=real_umphys) ::                                                      &
 beta_ls_rain(p_field),                                                        &
                                      ! out Scattering in LS Rain.
 beta_ls_snow(p_field),                                                        &
                                      ! out Scattering in LS Snow.
 beta_c_rain(p_field),                                                         &
                                      ! out Scattering in Conv Rain
 beta_c_snow(p_field)                 ! out Scattering in Conv Snow
!---------------------------------------------------------------------

real(kind=real_umphys)  ::                                                     &
 powerr,                                                                       &
 poweri,                                                                       &
 factorr1,                                                                     &
 factorr2,                                                                     &
 factori1,                                                                     &
 factori2
real(kind=real_umphys) ::                                                      &
 inst_ls_rain(p_field),                                                        &
                                        ! Local Large scale Rain
 inst_ls_snow(p_field),                                                        &
                                        ! Local Large scale Snow
 inst_c_rain(p_field),                                                         &
                                        ! Local Convective Rain
 inst_c_snow(p_field),                                                         &
                                        ! Local Convective Snow
 lcai(p_field),                                                                &
                                        ! 1 / large-scale cloud (lca)
 m_si(p_field)
                        ! si moment of ice particle size distribution

! Local variables:-----------------------------------------------------
real(kind=real_umphys) ::                                                      &
 pfactor(p_field),                                                             &
 gammar,                                                                       &
 gammar1,                                                                      &
 gammai,                                                                       &
 gammai1,                                                                      &
 smallvalue
parameter(smallvalue=1.0e-7)

!-----------------------------------------------------------------------
!  Define local variables ----------------------------------------------
integer :: i       ! Loop counters: I - horizontal field index;
integer :: error_status ! Error status code

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='BETA_PRECIP'
character(len=errormessagelength) :: err_msg ! Error message

!for casim calculations
real :: mui, mur, ar, br, xi, yi, xr, yr, aic, bic
real :: lam_casim_r, lam_casim_s
real :: gam_rain_cas1, gam_rain_cas2, gam_rain_cas3
real :: gam_snow_cas1, gam_snow_cas2, gam_snow_cas3
real :: rain_num(p_field), snow_num(p_field)
integer :: icnt, j,  k

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

error_status = 0

if ((k1stpt+points-1) >  p_field) then

  error_status = 1
  err_msg = 'Attempting to process over a field size which is '//newline//     &
            'larger than the input arrays specified.          '//newline//     &
            'Please check input code. '

  call ereport(RoutineName, error_status, err_msg)

end if

do i = k1stpt, k1stpt+points-1
  pfactor(i)=(pref/pressure(i))**0.4
end do


if (l_casim) then

  aic = snow_params%c_x
  bic = snow_params%d_x
  mui = snow_params%fix_mu
  xi=pi/4.0 !area=xi*D^yi
  yi=2.0
  ! Need ai and bi for use with convection related vis for global.
  ! Set to them to casim values
  ai=aic
  bi=bic

  ar = rain_params%c_x
  br = rain_params%d_x
  xr = pi/4.0 *0.5   !factor
  yr=2.0
  mur = rain_params%fix_mu

  call gammaf(1.0+mur+br, gam_rain_cas1 )
  call gammaf(1.0+mur, gam_rain_cas2)
  call gammaf(1.0+mur+yr, gam_rain_cas3)

  call gammaf(1.0+mui+bic, gam_snow_cas1 )
  call gammaf(1.0+mui, gam_snow_cas2)
  call gammaf(1.0+mui+yi, gam_snow_cas3)

  icnt=0
  do k = 1, 1
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        icnt = icnt+1
        rain_num(icnt) = rainnumber(i,j,k)
        snow_num(icnt) = snownumber(i,j,k)
      end do
    end do
  end do

end if !l_casim

powerr = (3.0-x2r+x4r)/(dr+4.0-x2r+x4r)
poweri = (1.0+si-x2i+x4i)/(bi+di0+1.0-x2i+x4i)
factorr1 = 0.5*pi*x1r
factorr2 = pi/6.0*rho_water*cr*x1r
factori1 = 2.0*ri*x1i
factori2 = ai*ci0*x1i

call gammaf(x4r+dr+4.0,gammar)
call gammaf(x4r+3.0,gammar1)
call gammaf(x4i+bi+di0+1.0,gammai)
call gammaf(x4i+si+1.0,gammai1)

if (l_psd) then
      ! Find inverse of layer cloud amount
  do i = 1, p_field
    lcai(i)=1.0/max(lca(i),0.01)
  end do

      ! Use the generic ice particle size distribution to
      ! calculate the si moment of the ice particle size distribution
  call lsp_moments(p_field,rho,t,qcf,lcai,si,m_si )

end if  ! l_psd

if (avg) then

  if (pct) then

    do i = k1stpt, k1stpt+points-1

      if (lca(i)  >   smallvalue .and. cca(i)  <   100.0) then
        inst_ls_rain(i) = 10000.0*ls_rain(i)/(100.0-cca(i))/lca(i)
        inst_ls_snow(i) = 10000.0*ls_snow(i)/(100.0-cca(i))/lca(i)
      else
        inst_ls_rain(i) = 0.0
        inst_ls_snow(i) = 0.0
      end if
      if (cca(i)  >   smallvalue) then
        inst_c_rain(i) = 100.0*c_rain(i)/cca(i)
        inst_c_snow(i) = 100.0*c_snow(i)/cca(i)
      else
        inst_c_rain(i) = 0.0
        inst_c_snow(i) = 0.0
      end if

    end do

  else

    do i = k1stpt, k1stpt+points-1

      if (lca(i)  >   smallvalue .and. cca(i)  <   1.0) then
        inst_ls_rain(i) = ls_rain(i)/(1.0-cca(i))/lca(i)
        inst_ls_snow(i) = ls_snow(i)/(1.0-cca(i))/lca(i)
      else
        inst_ls_rain(i) = 0.0
        inst_ls_snow(i) = 0.0
      end if
      if (cca(i)  >   smallvalue) then
        inst_c_rain(i) = c_rain(i)/cca(i)
        inst_c_snow(i) = c_snow(i)/cca(i)
      else
        inst_c_rain(i) = 0.0
        inst_c_snow(i) = 0.0
      end if

    end do

  end if

else

  do i = k1stpt, k1stpt+points-1
    inst_ls_rain(i) = ls_rain(i)
    inst_ls_snow(i) = ls_snow(i)
    inst_c_rain(i)  = c_rain(i)
    inst_c_snow(i)  = c_snow(i)
  end do

end if

do i = k1stpt, k1stpt+points-1

  if (.not. l_casim) then

    if (inst_ls_rain(i) > smallvalue) then
      beta_ls_rain(i) = factorr1*gammar1*(inst_ls_rain(i)/                     &
                      (pfactor(i)*factorr2*gammar))**powerr
    else
      beta_ls_rain(i) = 0.0
    end if

    if (inst_ls_snow(i) > smallvalue) then
      if (l_psd) then
            ! Use the generic ice particle size distribution
        beta_ls_snow(i) = 2.0*ri*m_si(i)
      else
        beta_ls_snow(i) = factori1*gammai1*(inst_ls_snow(i)/                   &
                        (pfactor(i)*factori2*gammai))**poweri
      end if  ! l_psd
    else
      beta_ls_snow(i) = 0.0
    end if


  else !casim
    !beta is simply 2*integrated cross sectional area for geometric optics
    if ((qrain(i) > qr_small) .and. (rain_num(i) > nr_small)) then
      lam_casim_r = (ar*rain_num(i)/qrain(i)*                                  &
                   gam_rain_cas1/gam_rain_cas2)**(1.0/br)
      beta_ls_rain(i) =  2.0  *xr*rain_num(i)/(lam_casim_r**yr)*               &
                   gam_rain_cas3/gam_rain_cas2
    else
      beta_ls_rain(i) = 0.0
    end if

    if ((qcf(i) > qs_small) .and. (snow_num(i) > ns_small)) then
      lam_casim_s = (aic*snow_num(i)/qcf(i)*                                   &
                  gam_snow_cas1/gam_snow_cas2)**(1.0/bic)
      beta_ls_snow(i) = 2.0*xi*snow_num(i)/(lam_casim_s**yi)*                  &
                  gam_snow_cas3/gam_snow_cas2
    else
      beta_ls_snow(i) = 0.0
    end if


  end if





  if (inst_c_rain(i)  >   smallvalue) then
    beta_c_rain(i) =factorr1*gammar1*(inst_c_rain(i) /                         &
                    (pfactor(i)*factorr2*gammar))**powerr
  else
    beta_c_rain(i) =0.0
  end if

  if (inst_c_snow(i)  >   smallvalue) then
        ! Use the 3C size distribution of snow for the convective
        ! contribution since we do not have the equivalent of qcf
        ! easily available.
    beta_c_snow(i) =factori1*gammai1*(inst_c_snow(i) /                         &
                    (pfactor(i)*factori2*gammai))**poweri
  else
    beta_c_snow(i) =0.0
  end if

end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine beta_precip
end module beta_precip_mod
