! ------------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
module deep_turb_cmt_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'DEEP_TURB_CMT_MOD'
contains

subroutine deep_turb_cmt (n_dp, nterm, nlev, deep_cmt_opt,                     &
                          ntml, kterm, cu_term, freeze_lev,                    &
                          timestep,                                            &
                          uw0, vw0, mb, wcld, wstar, zlcl,                     &
                          mass_flux,                                           &
                          z_rho, z_theta,                                      &
                          rho, rho_theta,                                      &
                          r2rho, r2rho_th, dr_across_th, dr_across_rh,         &
                          u, v,                                                &
                          dubydt, dvbydt, uw_diag, vw_diag)

! Modules

! numbers indicating type of convection
use conv_type_defs,    only: deep_conv
! parameters from CRM fits for deep CMT
use tcs_cmt_params_dp, only: a_kcmt_deep, a_wup_deep, b_wup_deep

use planet_constants_mod, only: g

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

use deep_turb_grad_stress_mod, only: deep_turb_grad_stress
use tcs_cb_stress_mod,         only: tcs_cb_stress
use tcs_cmt_incr_mod,          only: tcs_cmt_incr
implicit none

! ------------------------------------------------------------------------------
! Description:
!   This routine calculates convective momentum transport for deep convection
!   using turbulence ideas. This version is designed for use with the mass flux
!   convection scheme.
!
!Code Owner: Please refer to the UM file CodeOwners.txt
!This file belongs in section: Convection
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.
!
! ------------------------------------------------------------------------------
! Subroutine arguments

integer, intent(in) ::                                                         &
  n_dp                 & ! number of deep columns
 ,nterm                & ! Number of deep columns which actually convected
 ,nlev                 & ! Number of model levels for calculations
 ,deep_cmt_opt           ! 3 or 4 to indicate treatment of non-local part

integer, intent(in) ::                                                         &
  ntml(n_dp)           & ! Cloud base level information
 ,kterm(n_dp)          & ! cloud top level information
 ,cu_term(n_dp)        & ! Location of deep points which actually convected
 ,freeze_lev(n_dp)       ! First theta level with T < 273.15K

real(kind=real_umphys), intent(in) ::                                          &
  timestep               ! Convection timestep (s)

real(kind=real_umphys), intent(in) ::                                          &
  uw0(n_dp)            & ! surface uw stress (N/m2)
 ,vw0(n_dp)            & ! surface vw stress (N/m2)
 ,mb(n_dp)             & ! Cloud base mass flux (Pa/s)
 ,wcld(n_dp)           & ! convective velocity scale (m/s)
 ,wstar(n_dp)          & ! convective sub-cloud velocity scale (m/s)
 ,zlcl(n_dp)           & ! Exact height of LCL (m) ? (or model level height)
 ,mass_flux(n_dp,nlev) & ! Mass flux profile (on theta levels)
 ,z_rho(n_dp,nlev)     & ! Height above surface of rho levels (m)
 ,z_theta(n_dp,nlev)   & ! Height above surface of theta levels (m)
 ,rho(n_dp,nlev)       & ! rho levels (kg/m3)
 ,rho_theta(n_dp,nlev) & ! rho  theta levels (kg/m3)
 ,r2rho(n_dp,nlev)     & ! r*r*rho  rho levels (kg/m)
 ,r2rho_th(n_dp,nlev)  & ! r*r*rho  theta levels (kg/m)
 ,dr_across_th(n_dp,nlev) & ! thickness of theta levels (m)
 ,dr_across_rh(n_dp,nlev) & ! thickness of rho levels (m)
 ,u(n_dp,nlev)            & ! U component of wind (m/s)
 ,v(n_dp,nlev)              ! V component of wind (m/s)

real(kind=real_umphys), intent(out) ::                                         &
  dubydt(n_dp,nlev)    & ! dU/dt due to deep CMT (m/s/s)
 ,dvbydt(n_dp,nlev)    & ! dV/dt due to deep CMT (m/s/s)
 ,uw_diag(n_dp,nlev)   & ! uw stress profile on theta levels (N/m2)
 ,vw_diag(n_dp,nlev)     ! vw stress profile on theta levels (N/m2)


! Local declarations:
integer  ::                                                                    &
   i,k,ii             &  ! loop counter
  ,klev               &  ! level in cloud
  ,ilev               &  ! level in cloud
  ,itop               &  ! top level
  ,ntop_max           &  ! maximum cloud top level
  ,max_cldlev         &  ! maximum number of cloud levels
  ,max_ntml              ! maximum number of below cloud levels

integer  ::                                                                    &
  ntml_uv(n_dp)         ! lcl for UV calculations


real(kind=real_umphys) ::                                                      &
  zp                  & !
 ,eta_val             & !
 , factor

real(kind=real_umphys) ::                                                      &
  uw(n_dp,nlev)       & ! uw stress profile on theta levels (N/m2)
 ,vw(n_dp,nlev)         ! vw stress profile on theta levels (N/m2)

! arrays compressed to deep points which actually convected

integer ::                                                                     &
  ncld_thlev(nterm)  & ! number of cloud levels
 ,nstart(nterm)        ! level for ust

real(kind=real_umphys) ::                                                      &
  uw_cb(nterm)       & ! uw at cloud base (N/m2)
 ,vw_cb(nterm)       & ! vw at cloud base (N/m2)
 ,uw0_term(nterm)    & ! uw surface stress compressed to nterm (N/m2)
 ,vw0_term(nterm)    & ! vw surface stress compressed to nterm (N/m2)
 ,mb_term(nterm)     & ! cloud base mass flux in (m/s)
 ,zcld(nterm)        & ! cloud depth (m)
 ,zlcl_term(nterm)   & ! lifting condensation level (m)
 ,dz_cb(nterm)       & ! depth of cloud base layer (m)
 ,wup_cb(nterm)      & ! wup/wcld at cloud base
 ,fng_nlclp1(nterm)  & ! value of non-gradient function at lcl+1
 ,du_start(nterm)    & ! du across cloud base
 ,dv_start(nterm)    & ! dv across cloud base
 ,wup_peak(nterm)    & ! peak value of wup/wcld
 ,zcld_freeze(nterm) & ! depth of cloud from frezing level to top
 ,zfreeze(nterm)     & ! height of freezing level
 ,zsurf(nterm)         ! depth of surface layer

real(kind=real_umphys) ::                                                      &
  uw_cld(nterm,nlev)           & ! uw on cloud levels (N/m2)
 ,vw_cld(nterm,nlev)           & ! vw on cloud levels (N/m2)
 ,u_cld(nterm,nlev)            & ! u on cloud levels
 ,v_cld(nterm,nlev)            & ! v on cloud levels
 ,uth_cld(nterm,nlev)          & ! u on theta cloud levels
 ,vth_cld(nterm,nlev)          & ! v on theta cloud levels
 ,dr_across_rh_cld(nterm,nlev) & ! thickness of rho levels
 ,dr_across_th_cld(nterm,nlev) & ! thickness of theta levels
 ,r2rho_cld(nterm,nlev)        & ! r2*rho on rho levels
 ,r2rho_theta_cld(nterm,nlev)  & ! r2*rho on theta levels
 ,eta_rho(nterm,nlev)          & ! eta (non-dimensional in cloud depth)
                                 ! on rho levels
 ,eta_theta(nterm,nlev)        & ! eta (non-dimensional in cloud depth)
                                 ! on theta levels
 ,mf_cld(nterm,nlev)           & ! mass flux on cloud levels
 ,wup_cld(nterm,nlev)          & ! wup on cloud levels
 ,k_func(nterm,nlev)           & ! k function
 ,fng_func(nterm,nlev)           ! non-gradient function

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DEEP_TURB_CMT'


!-------------------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-------------------------------------------------------------------------------
! Initialise output arrays
!-------------------------------------------------------------------------------

do k=1,nlev
  do i=1,n_dp
    dubydt(i,k) = 0.0
    dvbydt(i,k) = 0.0
    uw(i,k)     = 0.0
    vw(i,k)     = 0.0
    uw_diag(i,k)  = 0.0
    vw_diag(i,k)  = 0.0
  end do
end do

!-------------------------------------------------------------------------------
! Picture of levels used for calculation
!-------------------------------------------------------------------------------

!   ---------------------  T,q

!   + + + + + + + + + + +  u,v   Cloud top  (kterm+1)          ^
!                                                              |
!   ---------------------  T,q   uw,vw     kterm               |
!                                                              |
!   - - - - - - - - - - -  u,v                                 |
!                                                              |
!   ---------------------  T,q  uw,vw                          |
!                                                              |
!   - - - - - - - - - - -  u,v                                zcld
!                                                              |
!   ---------------------  T,q   uw,vw                         |
!                                                              |
!   - - - - - - - - - - -  u,v                                 |
!                                                              |
!   ---------------------  T,q  uw,vw                          |
!                                                              |
!   - - - - - - - - - - -  u,v                                 |
!                                                              |
!   ---------------------  T,q  uw,vw                          |
!                                                              |
!   + + + + + + + + + + +  u,v   Cloud base  (ntml+1)  zlcl    v uw_cb
!
!   ---------------------  T,q   ntml
!
!   - - - - - - - - - - -  u,v
!
!   ---------------------  T,q level 1

!   - - - - - - - - - - -  u,v level 1

!   ---------------------  Surface
!   /////////////////////

! There are Kterm -ntml + 1 wind levels used for in cloud calculations
! (includes cloud base and cloud top winds) to produce uw, vw fluxes
! on  kterm - ntml in cloud levels.
! Below cloud base the flux is assumed to reduce linearly with height from
! the cloud base flux value.

!-------------------------------------------------------------------------------
! From this point the code operates on those points which actually convected
! which for the mass flux scheme is not equal to those first diagnosed as deep
! on this time step.
!-------------------------------------------------------------------------------

! work out maximum level values - used later to reduce loops

ntop_max = 0        ! maximum cloud top level
max_cldlev = 0      ! maximum number of cloud levels
max_ntml = 0        ! maximum below cloud

do i = 1,nterm
  ii = cu_term(i)
  if (kterm(ii) >  ntop_max) then
    ntop_max=kterm(ii)+1
  end if
  if (ntml(ii) >  max_ntml) then
    max_ntml = ntml(ii)
  end if

  ! Number of cloud levels - i.e. number of u,v levels required for in cloud
  ! calculation.

  ncld_thlev(i) = kterm(ii) - ntml(ii)+1

  if (ncld_thlev(i)>  max_cldlev) then
    max_cldlev = ncld_thlev(i)
  end if

end do    ! nterm

! Note this should never happen unless problem of convection going to the
! model top.
if (ntop_max > nlev-1) then
  ntop_max = nlev-1
end if

! Convert cloud base mass flux from Pa/s to m/s (divide by rho*g)

do i=1,nterm

  ii = cu_term(i)        ! location in full input arrays

  ! Convert cloud base mass flux from Pa/s to m/s (divide by rho*g)
  mb_term(i) = mb(ii)/(g*rho(ii,ntml(ii)+1))

  ntml_uv(i) = ntml(ii) + 1

  ! cloud depth
  zcld(i)  = z_rho(ii,kterm(ii)+1) - zlcl(ii)

  ! compress uw0 and vw0 to just those point which did deep convection
  uw0_term(i) = uw0(ii)
  vw0_term(i) = vw0(ii)
  zlcl_term(i) = zlcl(ii)

end do

!-------------------------------------------------------------------------------
! Compress to cloud levels and points which actually did deep convection
!-------------------------------------------------------------------------------
!CDIR NOUNROLL
do k=1,max_cldlev+1
  do i = 1,nterm
    ii = cu_term(i)
    klev = ntml(ii)+k

    if (klev >= nlev-1) then   ! problem will go outside array
      klev = nlev-1
    end if

    u_cld(i,k)   = u(ii,klev)
    v_cld(i,k)   = v(ii,klev)
    dr_across_rh_cld(i,k) = dr_across_rh(ii,klev)
    dr_across_th_cld(i,k) = dr_across_th(ii,klev)
    r2rho_cld(i,k)        = r2rho(ii,klev)
    r2rho_theta_cld(i,k)  = r2rho_th(ii,klev)

    ! mass flux converted from Pa/s to m/s

    mf_cld(i,k)  = mass_flux(ii,klev)/(g*rho_theta(ii,klev))

    if (k <= kterm(ii)+1) then

      eta_rho(i,k)   = (z_rho(ii,klev)  - zlcl_term(i))/zcld(i)
      eta_theta(i,k) = (z_theta(ii,klev)- zlcl_term(i))/zcld(i)

      ! visocity function  on theta levels (i.e. uw stress levels)
      !           k(eta) = a*eta

      k_func(i,k) = a_kcmt_deep*eta_theta(i,k)

      ! non-gradient function on theta levels

      fng_func(i,k) = exp(-eta_theta(i,k))

    else       ! above cloud top

      fng_func(i,k) = 0.0
      k_func(i,k) = 0.0
      eta_rho(i,k)   = 0.0
      eta_theta(i,k) = 0.0

    end if

  end do
end do


! Calculate u & v on theta cloud levels using linear interpolation

if (deep_cmt_opt == 4) then
  do k=1,max_cldlev+1
    do i = 1,nterm
      ii = cu_term(i)
      klev = ntml(ii)+k

      if (klev >= nlev-1) then
        klev = nlev-1
      end if
      factor = (z_theta(ii,klev) - z_rho(ii,klev))/dr_across_th(ii,klev)
      uth_cld(i,k)   = u(ii,klev+1)*factor+(1.0-factor)*u(ii,klev)
      vth_cld(i,k)   = v(ii,klev+1)*factor+(1.0-factor)*v(ii,klev)
    end do
  end do
end if


! wup = wup/wcld
! -------------------------------------------------------------------------
! Scheme assumes a shape for wup
! Peak value at freezing level wup = wcld  from CRM fit (provided freezing
! level is above cloud base).
! Value at cloud base ~ a+b* wstar  (sub-cloud velocity scale)
! Note encountered problems with fit. At many points in a full atmospheric run
! wcld is small or wstar bigger so that wup_cb > wcld and the shape wrong.
! Decided to enforce a condition on wup_cb/wcld <= 0.8.

do i=1,nterm
  ii = cu_term(i)
  wup_cb(i)   = (a_wup_deep+b_wup_deep*wstar(ii))/wcld(ii)

  if (wup_cb(i) > 0.8) then
    ! Don't use wstar and wcld information just assume a shape.
    wup_cb(i) = 0.8
  end if

  if (freeze_lev(ii) > ntml(ii)+1) then  ! freezing level above cloud base

    ! Value at freezing level wup_peak = wcld/wcld
    wup_peak(i)    = 1.0
    zcld_freeze(i) = z_rho(ii,kterm(ii)+1) - z_theta(ii,freeze_lev(ii))
    zfreeze(i)     = z_theta(ii,freeze_lev(ii))

  else                                   ! set to cloud base value

    ! peak value assume to be cloud value of wup
    wup_peak(i)    = wup_cb(i)
    zcld_freeze(i) = zcld(i)
    zfreeze(i)     = zlcl_term(i)

  end if

end do

!CDIR NOUNROLL
do k=1,max_cldlev+1
  do i = 1,nterm
    ii = cu_term(i)     ! location in full array
    klev = ntml(ii)+k

    if (klev < kterm(ii)+1) then    ! work out in cloud wup

      if (klev < freeze_lev(ii)) then  ! below freezing level

        eta_val = (z_theta(ii,klev) - zlcl_term(i))                            &
                                      /(zfreeze(i)-zlcl_term(i))
        factor = wup_cb(i)*wup_cb(i)
        wup_cld(i,k) = (factor + (1.0-factor)*eta_val)**0.5

      else        ! above freezing level

        ! straight line from peak to zero at top

        zp = z_theta(ii,klev) - zfreeze(i)
        wup_cld(i,k) = wup_peak(i)* (1.0-zp/zcld_freeze(i))

      end if       ! test on freezing level

    else                            ! above cloud top

      wup_cld(i,k) = 0.0

    end if                          ! test on klev in cloud

  end do
end do

!-------------------------------------------------------------------------------
! Gradient part of stress
!-------------------------------------------------------------------------------

call deep_turb_grad_stress(nterm, nlev, max_cldlev, ncld_thlev                 &
                      , timestep                                               &
                      , zcld                                                   &
                      , mf_cld, wup_cld, k_func                                &
                      , u_cld, v_cld                                           &
                      , r2rho_cld, r2rho_theta_cld                             &
                      , dr_across_rh_cld, dr_across_th_cld                     &
                      ,uw_cld,vw_cld)


!-------------------------------------------------------------------------------
! Non-gradient part of stress
! Note this code does not calculate the non-gradient stress in exactly the
! same way as the old turbluence based code documented by A Grant.
!-------------------------------------------------------------------------------
if (deep_cmt_opt == 4) then      ! version using m(ust-u)*(1-eta)
                                 ! This version is still experimental.

  ! need level near surface  = 0.1*zlcl
  do i=1,nterm
    zsurf(i) = 0.1*zlcl_term(i)
  end do
  k=1
  do i=1,nterm
    ii = cu_term(i)
    if (zsurf(i) <= z_theta(ii,k)) then
      nstart(i) = k
    end if
  end do
  !CDIR NOUNROLL
  do k=2,max_ntml
    do i=1,nterm
      ii = cu_term(i)
      if (zsurf(i) <= z_theta(ii,k) .and. zsurf(i) > z_theta(ii,k-1)) then
        nstart(i) = k
      end if
    end do
  end do

  !CDIR NOUNROLL
  do k=1,max_cldlev+1
    do i = 1,nterm
      ii = cu_term(i)
      klev = ntml(ii)+k
      if (klev < kterm(ii)+1) then    ! only in cloud
        uw_cld(i,k) = uw_cld(i,k)                                              &
                       + mf_cld(i,k)*(u(ii,nstart(i))-uth_cld(i,k))*           &
                          (1.0-eta_theta(i,k))
        vw_cld(i,k) = vw_cld(i,k)                                              &
                       + mf_cld(i,k)*(v(ii,nstart(i))-vth_cld(i,k))*           &
                          (1.0-eta_theta(i,k))
      end if
    end do
  end do
  do i=1,nterm
    ii = cu_term(i)
    uw_cb(i) = mb_term(i)*(u(ii,nstart(i))-u(ii,ntml(ii)+1))
    vw_cb(i) = mb_term(i)*(v(ii,nstart(i))-v(ii,ntml(ii)+1))
  end do


else              ! option 3 non-local term uses cloud base stress

  ! Values across cloud base - at present uses adjacent model levels
  ! but could be altered.

  do i=1,nterm
    ii = cu_term(i)
    dz_cb(i) = z_rho(ii,ntml_uv(i)) - z_rho(ii,ntml_uv(i)-1)
    du_start(i) = u(ii,ntml_uv(i)) - u(ii,ntml_uv(i)-1)
    dv_start(i) = v(ii,ntml_uv(i)) - v(ii,ntml_uv(i)-1)

    !  Non gradient function at first in cloud level

    fng_nlclp1(i) = exp(-eta_rho(i,1))  ! Required at a rho level

  end do


  ! Calculate cloud base stresses

  call tcs_cb_stress (deep_conv, nterm                                         &
                     ,timestep                                                 &
                     ,uw0_term, vw0_term, mb_term, zlcl_term                   &
                     ,du_start, dv_start, dz_cb, fng_nlclp1                    &
                     ,uw_cb, vw_cb )

  ! Add non-gradient part of stress to gradient stress

  !CDIR NOUNROLL
  do k=1,max_cldlev+1
    do i = 1,nterm

      uw_cld(i,k) = uw_cld(i,k) + uw_cb(i)*fng_func(i,k)
      vw_cld(i,k) = vw_cld(i,k) + vw_cb(i)*fng_func(i,k)

    end do
  end do

end if        ! test on non-local option

!-------------------------------------------------------------------------------
! Expand in cloud values on uv levels

!CDIR NOUNROLL
do k = 1,max_cldlev+1
  do i = 1,nterm
    ii = cu_term(i)
    itop = kterm(ii)+1
    ilev = ntml(ii)+k
    if (ilev < itop) then
      uw(ii,ilev)  = uw_cld(i,k)
      vw(ii,ilev)  = vw_cld(i,k)
    end if
  end do
end do

! Below cloud base assume uwcb linearly to zero

!CDIR NOUNROLL
do k = 1,max_ntml
  do i = 1,nterm
    ii = cu_term(i)
    ilev = ntml(ii)
    if (k <= ilev) then
      uw(ii,k)  = uw_cb(i)*z_theta(ii,k)/zlcl_term(i)
      vw(ii,k)  = vw_cb(i)*z_theta(ii,k)/zlcl_term(i)
    end if
  end do
end do

!-------------------------------------------------------------------------------
! Calculate increments to U and V - full arrays


call tcs_cmt_incr(n_dp,nlev, ntop_max, kterm                                   &
                 ,r2rho, r2rho_th                                              &
                 ,dr_across_rh                                                 &
                 ,uw ,vw                                                       &
                 !out
                 ,dubydt,dvbydt)

!-------------------------------------------------------------------------------
! Copy stress to output arrays multiplying by density on theta levels
! so that diagnostic output has the correct units for stress rather than
! m2/s2.

do k=1,ntop_max+1
  do i=1,n_dp
    uw_diag(i,k)=uw(i,k)*rho_theta(i,k)
    vw_diag(i,k)=vw(i,k)*rho_theta(i,k)
  end do
end do
!-------------------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

!-------------------------------------------------------------------------------
end subroutine deep_turb_cmt
end module deep_turb_cmt_mod
