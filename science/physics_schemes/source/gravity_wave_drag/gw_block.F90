! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
module gw_block_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'GW_BLOCK_MOD'
contains

subroutine gw_block(levels,points,gw_seg_size,dt,u,v,rho,nsq,ulow,vlow,        &
                    modulow,rho_levels,theta_levels,mt_high,sd_orog,           &
                    slope,anis,                                                &
                    mtdir,zb,banis,canis,dudt,dvdt,dtdt,                       &
                    fbcd,fcrit,l_drag,l_fb_heating,                            &
!diagnostics
                    du_dt_wake,points_du_dt_wake,du_dt_wake_on,                &
                    dv_dt_wake,points_dv_dt_wake,dv_dt_wake_on,                &
                    stress_ud,stress_ud_on,points_stress_ud,                   &
                    stress_ud_p_on,                                            &
                    stress_vd,stress_vd_on,points_stress_vd,                   &
                    stress_ud_wake,points_stress_ud_wake,stress_ud_wake_on,    &
                    stress_vd_wake,points_stress_vd_wake,stress_vd_wake_on,    &
                    fr_d, fr_d_on,points_fr_d,                                 &
                    bld_d, bld_d_on,points_bld_d,                              &
                    bldt_d, bldt_d_on, points_bldt_d,                          &
                    tausx_d,tausx_d_on,points_tausx_d,                         &
                    tausy_d,tausy_d_on,points_tausy_d)

! subroutine gw_block to calculate the vertical profile of stress and
!            associated wind increments due to flow blocking

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use c_gwave_mod, only: lambdaz_min, lambdaz_max,                               &
                       nsq_neutral, zav_converge, zav_iterate
use planet_constants_mod, only: cp
implicit none

! Description:
!     calculates the flow blocking stress profiles and wind increments.
!     1. calculate stress profile and wind increments for
!        the blocked flow. Calculations based on Lott and Miller (1997)
!        with modification to blocking layer calculation
!        from Vosper et al (2009).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: gravity_wave_drag
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.3 programming standards.

! Local constants

!----------------------
! Intent(in) variables
!----------------------
integer, intent(in) :: &!block for intent(in)
  levels,              &!num of vertical levels
  points,              &!num of land points
  gw_seg_size

integer, intent(in)  :: &! intent(in) for diags
  points_stress_ud_wake,                                                       &
  points_stress_vd_wake,                                                       &
  points_stress_ud,                                                            &
  points_stress_vd,                                                            &
  points_du_dt_wake,                                                           &
  points_dv_dt_wake,                                                           &
  points_fr_d,                                                                 &
  points_bld_d,                                                                &
  points_bldt_d,                                                               &
  points_tausx_d,                                                              &
  points_tausy_d

real(kind=real_umphys), intent(in) ::           &!block for intent(in)
  u(gw_seg_size,levels),           &!zonal wind
  v(gw_seg_size,levels),           &!meridional wind
  rho(gw_seg_size,levels),         &!density
  nsq(gw_seg_size,levels),         &!brunt-vaisala freq squared
  rho_levels(gw_seg_size,levels),  &!height (rho_levels)
  theta_levels(gw_seg_size,levels),&!height (theta_levels)
  slope(points),                   &!sso params - slope
  anis(points),                    &!sso params - anisotropy
  mtdir(points),                   &!sso params - angle of major axis
  mt_high(points),                 &!sso params - height
  sd_orog(points),                 &!standard deviation of orography (m)
  ulow(points),                    &!u averaged from z=0.5mt_high to z=mt_high
  vlow(points),                    &!v averaged from z=0.5mt_high to z=mt_high
  modulow(points),                 &!modulus of hrz wind (sqrt(ulow^2+vlow^2)
  banis(points),                   &!const (function of anisotropy)
  canis(points),                   &!const (function of anisotropy)
  fcrit,                           &!critical Froude number
  fbcd,                            &!flow blocking drag coefficient
  dt                                !time-step

logical, intent(in) ::   &!
  l_drag(points),        &!whether point has a non-zero stress or not
  l_fb_heating            !calculate heating tendency if true

! Below are the stash flags for calculating diagnostics:
logical, intent(in) ::  &!intent(in) for diags
  stress_ud_wake_on,    &!u wake stress
  stress_vd_wake_on,    &!v wake stress
  stress_ud_on,         &!total u stress
  stress_ud_p_on,       &!total u stress on p-points
  stress_vd_on,         &!total v stress
  du_dt_wake_on,        &!u accel blocked flow
  dv_dt_wake_on,        &!v accel blocked flow
  fr_d_on,              &!fr_d switch
  bld_d_on,             &!bld_d switch
  bldt_d_on,            &!bldt_d switch
  tausx_d_on,           &!tausx_d switch
  tausy_d_on             !tausy_d switch
!----------------------
! Intent(out) variables
!----------------------
real(kind=real_umphys), intent(out) ::     &!block for intent(out)
  dudt(gw_seg_size,levels),                                                    &
  dvdt(gw_seg_size,levels),   &!profiles of acceleration
  dtdt(gw_seg_size,levels),   &!profiles of heating
  zb(points)                   !depth of flow blocking layer

real(kind=real_umphys), intent(out)  ::                                        &
  stress_ud      ( points_stress_ud      , 0:levels ),                         &
  stress_vd      ( points_stress_vd      , 0:levels ),                         &
  stress_ud_wake ( points_stress_ud_wake , 0:levels ),                         &
  stress_vd_wake ( points_stress_vd_wake , 0:levels ),                         &
  du_dt_wake ( points_du_dt_wake , levels ),                                   &
  dv_dt_wake ( points_dv_dt_wake , levels ),                                   &
  fr_d(points_fr_d),                                                           &
  bld_d(points_bld_d),                                                         &
  bldt_d(points_bldt_d),                                                       &
  tausx_d(points_tausx_d),                                                     &
  tausy_d(points_tausy_d)

!----------------------
! Local variables
!----------------------

real(kind=real_umphys) ::       &!block for local variables
 psi(points,levels),            &!wind direction rel. to major axis of SSO
 local_xstress(points,0:levels),&!used in stress calculation for diags
 local_ystress(points,0:levels)

real(kind=real_umphys) ::          &!block for local variables
 zneu(points),   &!depth of near surface neutral layer
 fav(points),    &!Froude number for calc zb (zb=max(0,mt_high(fcrit-fav))
 zav(points),    &!depth used to calc fav
 zav1(points),   &!used to calculate zav1
 zav_new(points),&!used in calculation of zav
 u_n(points),    &!wind speed div by buoyancy freq
 nav(points),    &!bulk averaged n^2 from z=0 to z=zav
 uav(points),    &!u averaged from z=0 to z=zav
 vav(points)      !v averaged from z=0 to z=zav

real(kind=real_umphys) ::          &!block for local variables
 cpsi,        &!(cos(psi))**2
 spsi,        &!(sin(psi))**2
 ratio,       &!aspect ratio of SSO as seen by incident flow
 modu,        &!modulus of horizontal wind (i.e. sqrt(u^2+v^2)
 width,       &!mountain width seen by incident flow (i.e.sqrt((zb-z)/(z+h))
 zzd1,        &!direction of blocking drag (i.e.bcos(psi)^2+csin(psi)^2)
 dblk,        &!flow blocking drag = dblk*u
 wind,        &!wind speed resolved in direction of low-level flow.
 dzt,         &!layer thickness used to calc average u etc.
 dzb,         &!layer thickness used to calc average u etc.
 rdt,         &!reciprocal of timestep (ie. 1/dt)
 delta_z,     &!for stress diag calc
 dzz,         &!layer thickness used to calculate heating increment
 ududt,       &!Rate of change of kinetic energy after wind increments
 vdvdt,       &!Rate of change of kinetic energy after wind increments
 uhat,        &!u on theta level
 vhat          !v on theta level
integer :: i, ii, k
integer :: ktop_fb(points)

logical :: l_cont(points),                                                     &
           l_cont2(points),                                                    &
           l_cont3(points)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='GW_BLOCK'

!-------------------------------------------------------------------
!   1.0 start  preliminaries
! initialise increment and increment diagnostics
!------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
if (tausx_d_on) then
  do i=1,points
    tausx_d(i) = 0.0
  end do
end if
if (tausy_d_on) then
  do i=1,points
    tausy_d(i) = 0.0
  end do
end if
if (fr_d_on) then
  do i = 1, points
    fr_d(i)   = 0.0
  end do
end if
if (bld_d_on) then
  do i = 1, points
    bld_d(i)  = 0.0
  end do
end if
if (bldt_d_on) then
  do i = 1, points
    bldt_d(i) = 0.0
  end do
end if
if (du_dt_wake_on) then
  do k = 1, levels
    do i = 1, points
      du_dt_wake(i,k)=0.0
    end do
  end do
end if
if (dv_dt_wake_on) then
  do k = 1, levels
    do i = 1, points
      dv_dt_wake(i,k)=0.0
    end do
  end do
end if
if (stress_ud_wake_on) then
  do k = 0, levels
    do i = 1, points
      stress_ud_wake(i,k)=0.0
    end do
  end do
end if
if (stress_vd_wake_on) then
  do k = 0, levels
    do i = 1, points
      stress_vd_wake(i,k)=0.0
    end do
  end do
end if
if (stress_ud_on .or. stress_ud_p_on) then
  do k = 0, levels
    do i = 1, points
      stress_ud(i,k)=0.0
    end do
  end do
end if
if (stress_vd_on) then
  do k = 0, levels
    do i = 1, points
      stress_vd(i,k)=0.0
    end do
  end do
end if
if (stress_ud_on .or. stress_ud_p_on .or. stress_ud_wake_on .or.               &
    tausx_d_on) then
  do k = 0, levels
    do i = 1, points
      local_xstress(i,k)=0.0
    end do
  end do
end if
if (stress_vd_on .or. stress_vd_wake_on .or. tausy_d_on) then
  do k = 0, levels
    do i = 1, points
      local_ystress(i,k)=0.0
    end do
  end do
end if
do k = 1, levels
  do i = 1, points
    dudt(i,k) = 0.0
    dvdt(i,k) = 0.0
    dtdt(i,k) = 0.0
  end do
end do

do i = 1, points
  zav_new(i) = 0.0
  zneu(i)    = 0.0
  uav(i)     = 0.0
  vav(i)     = 0.0
  nav(i)     = 0.0
  fav(i)     = -1.0
  zb(i)      = 0.0
  ktop_fb(i) = 2
  l_cont(i)  = .true.
  l_cont2(i) = .true.
  l_cont3(i) = .true.
end do
rdt = 1.0 / dt

! ------------------------------------------------
!  calculate zav following Vosper et al (2009)
! ------------------------------------------------

! calculate depth of neutral layer - zneu
do k = 1, levels
  do i = 1, points
    if (l_drag(i)) then
      if (l_cont(i)) then
        if (nsq(i,k) < nsq_neutral) then
          zneu(i)    = rho_levels(i,k)
        else
          l_cont(i)  = .false.
        end if !nsq(i,k) < nsq_neutral
      end if !l_cont
    end if !l_drag(i)
  end do!i = 1, points
end do !k = 1, levels

! find max(mt_high,zneu) for zav
do i = 1, points
  if (l_drag(i)) then
    if (zneu(i) > mt_high(i)) then
      zav1(i) = zneu(i)
    else
      zav1(i) = mt_high(i)
    end if
    zav(i)     = zav1(i)
  end if!l_drag(i)
end do !i = 1, points

do ii = 1, zav_iterate
  !Reset l_cont3 (logical to test if rho_levels(k) lt zav)
  do i = 1, points
    l_cont3(i) = .true.
  end do !i = 1, points
  do k = 2, levels-1
    do i = 1, points
      if (l_drag(i)) then
        if (l_cont2(i)) then !l_cont2(i) tests if zav is converged
          !-----------------------------------------------------
          ! need an if test for case where zav doesn't
          ! converge in zav_iterate iterations?
          !-----------------------------------------------------
          if (l_cont3(i)) then !l_cont3(i) tests if rho_levels(k) lt zav
            if (theta_levels(i,k) < zav(i)) then
              dzt = theta_levels(i,k)
              if (k == 2) then
                dzb    = dzt
              else
                dzb    = theta_levels(i,k) -  theta_levels(i,k-1)
              end if !(k==2)
            else
              dzt    = zav(i)
              if (k == 2) then
                dzb    = dzt
              else
                dzb     = zav(i) - theta_levels(i,k-1)
                ktop_fb(i) = k
                l_cont3(i)  = .false.
              end if !(k==2)
            end if !theta_levels(i,k)<zav(i)
            !-----------------------------------------------------
            ! average u,v and n from z = 0 to current level
            ! which when k = ktop_fb become z = 0 tO z = zav(i) averages
            ! nav(i) calc is equivalent to bulk average n
            ! i.e. n^2 = sqrt(g/theta0*thetaav-theta0/zav(i))
            !-----------------------------------------------------
            uav(i) = (uav(i)*theta_levels(i,k-1) + u(i,k)*dzb)  / dzt
            vav(i) = (vav(i)*theta_levels(i,k-1) + v(i,k)*dzb)  / dzt
            nav(i) = (nav(i)*theta_levels(i,k-1) + nsq(i,k)*dzb)/ dzt
          end if    ! l_cont3
        end if    ! l_cont2
      end if !l_drag
    end do ! i=1, points
  end do   ! loop over k levels
  do i = 1, points
    if (l_drag(i)) then
      if (l_cont2(i)) then !l_cont2(i) tests if zav is converged
        ! resolve wind in the direction of the low-level flow
        wind    = (ulow(i)*uav(i) + vlow(i)*vav(i)) /modulow(i)
        wind    = abs(wind)
        if (nav(i) > nsq_neutral) then
          ! first consider stable cases
          u_n(i) = wind/sqrt(nav(i))
        else
          ! now consider neutral (and near neutral) cases
          u_n(i) = wind/sqrt(nsq_neutral)
        end if
        ! limit u_n to sensible values
        u_n(i) = max(lambdaz_min,u_n(i))
        u_n(i) = min(lambdaz_max,u_n(i))
        zav_new(i) = zav1(i) + u_n(i)
        ! currently set zav_converge to 0.05
        if ((zav(i) < zav_new(i)*(1.0+zav_converge)) .and.                     &
           ( zav(i) > zav_new(i)*(1.0-zav_converge))) then
          l_cont2(i)  =  .false.
        end if ! test if zav(i) is converged
        zav(i)    = zav_new(i)
      end if  ! l_cont2(i)
    end if!l_drag(i)
  end do!i = 1, points
end do  !ii= 1, zav_iterate

! -----------------------------------------
!  calculate fav and blocked layer depth
! -----------------------------------------
do i = 1, points
  l_cont(i)  =  .true.
  if (l_drag(i)) then
    ! calculate froude number
    ! prevent div by zero if nav(i) =0.
    if (nav(i) > 0.0) then
      !limit fav with wavelength too
      fav(i) = u_n(i)/mt_high(i)
    else
      fav(i) = -1.0
    end if
    ! find zb and variable drag coefficient
    if ((fav(i) > 0.0) .and. ((fcrit - fav(i)) > 0.0)) then
      zb(i) = mt_high(i)*(1 - fav(i)/fcrit)
    else
      zb(i) = 0.0
    end if
    ! find level at top of blocked layer
    ktop_fb(i) = 0
  end if!l_drag(i)
end do !i = 1, points
do k = 1, levels
  do i = 1, points
    if ((l_drag(i)) .and. (zb(i) > 0.0)) then
      if (l_cont(i)) then
        if (rho_levels(i,k) >= zb(i)) then
          ktop_fb(i)   = k
          l_cont(i) = .false.
        end if  !rho_levels(i,k)>=zb(i)
        ! calculate psi at each level
        psi(i,k) = atan2 (v(i,k),u(i,k))
        psi(i,k) = mtdir(i) - psi(i,k)
      end if  ! l_cont
    end if!l_drag(i)
  end do!i = 1, points
end do !k=1, levels
! --------------------------------
! calculate flow blocking drag
! --------------------------------
do k = 1, levels
  do i = 1, points
    if ((l_drag(i)) .and. (zb(i) > 0.0)) then
      if (k <  ktop_fb(i)) then
        cpsi  = (cos(psi(i,k)))**2
        spsi  = 1.0-cpsi
        ratio = 0.0
        !Note ratio expression revised from LM97 due to mistake in derivation in LM97
        if ((anis(i)**2*cpsi + spsi) /=0) then
          ratio = sqrt((cpsi + anis(i)**2*spsi) / (anis(i)**2*cpsi + spsi))
        end if
        if (ratio /=0) then
          ratio =  2.0 - 1.0/ratio
        end if
        if (ratio < 0) then
          ratio = 0.0
        end if
        ! calculate modu using uav and vav at previous time-step
        modu  = sqrt(uav(i)**2 + vav(i)**2)
        width = sqrt((zb(i) - rho_levels(i,k)) /                               &
                    (rho_levels(i,k) + sd_orog(i)))
        zzd1  = banis(i)*cpsi + canis(i)*spsi
        !------------------------------------------------------------------
        ! calculate tendencies using partially implicit formulation.
        !------------------------------------------------------------------
        dblk         = -fbcd*ratio*slope(i)/(2.0*sd_orog(i))                   &
                       *width*zzd1*0.5
        dblk         = ((1.0/(1.0-dblk*modu*dt))-1.0)*rdt
        dudt(i,k)    =  u(i,k) * dblk
        dvdt(i,k)    =  v(i,k) * dblk
      end if !(k <  ktop_fb(i))
    end if !(l_drag(i))
  end do !Loop over points
end do !Loop over levels

!-----------------------------------------------------------------
! Calculate heating due to dissipation
!-----------------------------------------------------------------
if ( l_fb_heating ) then

  do k = 1, levels-1
    do i = 1, points
      if ((l_drag(i)) .and. (zb(i) > 0.0)) then

        dzb          = theta_levels(i,k)   -  rho_levels(i,k)
        dzt          = rho_levels(i,k+1)   -  theta_levels(i,k)
        dzz          = rho_levels(i,k+1)   -   rho_levels(i,k)

        !          u and v on theta_level(k)
        uhat       =  dzt * u(i,k) + dzb * u(i,k+1)
        vhat       =  dzt * v(i,k) + dzb * v(i,k+1)

        !          u*du/dt abd v*dv/dt on theta_level(k)
        ududt      = uhat *( dzt * dudt(i,k) + dzb * dudt(i,k+1) )
        vdvdt      = vhat *( dzt * dvdt(i,k) + dzb * dvdt(i,k+1) )

        !          dT/dt on theta_level(k)
        dtdt(i,k)  = - (ududt + vdvdt) / ( cp * dzz * dzz )

      end if !(l_drag(i))
    end do !Loop over points
  end do !Loop over levels

end if !l_fb_heating


!-----------------------------------------------------------------
! 4 diagnostics
!-----------------------------------------------------------------

if ( bld_d_on ) then
  do i=1,points
    if (l_drag(i)) then
      bld_d(i) = zb(i)
    end if
  end do
end if
if ( bldt_d_on ) then
  do i=1,points
    if ((dudt(i,1) /=  0.0) .or. (dvdt(i,1) /= 0.0 )) then
      bldt_d(i) = 100.0
    else
      bldt_d(i) = 0.0
    end if
  end do
end if
if ( fr_d_on ) then
  do i=1,points
    if (l_drag(i)) then
      fr_d(i) = fav(i)
    end if
  end do
end if
if ( du_dt_wake_on ) then
  do k=1, levels
    do i=1,points
      if ( l_drag(i) .and. k <= ktop_fb(i) ) then
        du_dt_wake(i,k) = dudt(i,k)
      end if
    end do
  end do
end if
if ( dv_dt_wake_on ) then
  do k=1, levels
    do i=1,points
      if ( l_drag(i) .and. k <= ktop_fb(i) ) then
        dv_dt_wake(i,k) = dvdt(i,k)
      end if
    end do
  end do
end if

if (stress_ud_on .or. stress_ud_p_on .or. stress_ud_wake_on .or.               &
    tausx_d_on) then
  do k=levels-1, 0, -1
    do i=1,points
      if ( l_drag(i) .and. k <= ktop_fb(i) ) then
        if ( k  ==  0 ) then
          delta_z = theta_levels(i,k+1)
        else
          delta_z = theta_levels(i,k+1) - theta_levels(i,k)
        end if
        local_xstress(i,k) = local_xstress(i,k+1) -                            &
                              (dudt(i,k+1)*rho(i,k+1)*delta_z)
      end if
    end do
  end do
  if ( stress_ud_wake_on) then
    do k=0, levels-1
      do i=1,points
        stress_ud_wake(i,k) = local_xstress(i,k)
      end do
    end do
    do i=1,points
      stress_ud_wake(i,levels) = stress_ud_wake(i,levels-1)
    end do
  end if
  if ( stress_ud_on .or. stress_ud_p_on) then
    do k=0, levels-1
      do i=1,points
        stress_ud(i,k) = local_xstress(i,k)
      end do
    end do
    do i=1,points
      stress_ud(i,levels) = stress_ud(i,levels-1)
    end do
  end if
  if (tausx_d_on) then
    do i=1,points
      if ( l_drag(i)) then
        tausx_d(i) = local_xstress(i,0)
      end if
    end do
  end if
end if

if (stress_vd_on .or. stress_vd_wake_on .or. tausy_d_on) then
  do k=levels-1, 0, -1
    do i=1,points
      if ( l_drag(i) .and. k <= ktop_fb(i) ) then
        if ( k  ==  0 ) then
          delta_z = theta_levels(i,k+1)
        else
          delta_z = theta_levels(i,k+1) - theta_levels(i,k)
        end if
        local_ystress(i,k) = local_ystress(i,k+1) -                            &
                              (dvdt(i,k+1)*rho(i,k+1)*delta_z)
      end if
    end do
  end do
  if ( stress_vd_wake_on) then
    do k=0, levels-1
      do i=1,points
        stress_vd_wake(i,k) = local_ystress(i,k)
      end do
    end do
    do i=1,points
      stress_vd_wake(i,levels) = stress_vd_wake(i,levels-1)
    end do
  end if
  if ( stress_vd_on ) then
    do k=0, levels-1
      do i=1,points
        stress_vd(i,k) = local_ystress(i,k)
      end do
    end do
    do i=1,points
      stress_vd(i,levels) = stress_vd(i,levels-1)
    end do
  end if
  if (tausy_d_on) then
    do i=1,points
      if ( l_drag(i)) then
        tausy_d(i) = local_ystress(i,0)
      end if
    end do
  end if
end if


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine gw_block
end module gw_block_mod
