! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Large-scale precipitation scheme. Moist Froude calculation

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

! Subroutine Interface:
module lsp_froude_moist_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='LSP_FROUDE_MOIST_MOD'

contains

subroutine lsp_froude_moist(levels, points,                                    &
                           u, v, theta,                                        &
                           rho_levels, theta_levels,                           &
                           hmteff, zb)

! Code to calculate moist Brunt-Vaisala frequency squared
! which is then used to calculate the moist Froude number
! and the amount of blocking
! Copied gw_setup.F90  gw_block.F90 from gravity_wave_drag
! and modified

use planet_constants_mod, only: g
use mphys_inputs_mod,     only: fcrit

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use c_gwave_mod, only: lambdaz_min, lambdaz_max,                               &
                       nsq_neutral, zav_converge, zav_iterate

implicit none

!----------------------
! Intent(in) variables
!----------------------
integer, intent(in) :: &! INPUT array dimensions
  levels,              &!num of vertical levels
  points                !num of land points

real(kind=real_umphys), intent(in) ::           &! ARRAYS
  u(points,levels),           &!zonal wind (on P levels)
  v(points,levels),           &!meridional wind (on P levels)
  theta(points,levels),       &!potential temperature
  theta_levels(points,levels), &!height(theta_levels)
  rho_levels(points,levels)    !height (rho_levels)

!----------------------
! Intent(inout) variables
!----------------------
real(kind=real_umphys), intent(in out) ::                                      &
  zb(points),            & !Depth of flow blocking layer
  hmteff(points)           !Effective mountain height producing ascent

! Local variables
real(kind=real_umphys) ::                                                      &
  mt_high(points),    &!sso height (n_sigma*sd_orog)
  nsq(points,levels), &!brunt-vaisala freq squared
  ulow(points),       &!u averaged from z=0.5mt_high to z=mt_high
  vlow(points),       &!v averaged from z=0.5mt_high to z=mt_high
  modu(points),       &!modulus of horizontal wind (i.e. sqrt(u^2+v^2)
  nlow(points),       &!N bulk averaged from z=0.5mt_high to z=mt_high
  zneu(points),       &!depth of near surface neutral layer
  zav(points),        &!depth used to calc fav
  zav1(points),       &!used to calculate zav1
  zav_new,            &!used in calculation of zav
  u_n(points),        &!wind speed div by buoyancy freq
  nav(points),        &!bulk averaged n^2 from z=0 to z=zav
  uav(points),        &!u averaged from z=0 to z=zav
  vav(points),        &!v averaged from z=0 to z=zav
  fav(points)          !Froude number for calc zb
                       ! (zb=max(0,mt_high(fcrit-fav))

real(kind=real_umphys) ::                                                      &
 wind,        &!wind speed resolved in direction of low-level flow.
 dzt,         &!layer thickness used to calc average u etc.
 dzb,         &!layer thickness used to calc average u etc.
 dzu

integer :: i, k, ii
integer :: kbot(points)
integer :: ktop(points)    ! From setup

logical :: l_sfblock(points),  &   !whether point has a non-zero stress or not
           l_cont(points),                                                     &
           l_cont2(points),                                                    &
           l_cont3(points)


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_FROUDE_MOIST'
! ----------------------------------------
! initialise arrays
! ----------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise BV freq profile
!$OMP PARALLEL DEFAULT(NONE) PRIVATE(i,k) SHARED(levels,points,nsq,hmteff,     &
!$OMP mt_high,ulow,vlow,nlow,ktop,kbot,zneu,uav,vav,nav,l_sfblock,u, l_cont,   &
!$OMP l_cont2,l_cont3,theta,theta_levels,g)
!$OMP DO SCHEDULE(STATIC)
do k = 1, levels
  do i =1, points
    nsq(i,k)=0.0
  end do
end do
!$OMP END DO

! Initialise low-level parameters
!$OMP DO SCHEDULE(STATIC)
do i=1,points
  ! Mountain height
  mt_high(i) = hmteff(i)
  ! Other local variables
  ulow(i)    = 0.0
  vlow(i)    = 0.0
  nlow(i)    = 0.0
  ktop(i)    = 2
  kbot(i)    = 1
  zneu(i)    = 0.0
  uav(i)     = 0.0
  vav(i)     = 0.0
  nav(i)     = 0.0
  l_sfblock(i)  = .true.
  !u(i,1)=0 implies at north/south pole in global model
  if ((mt_high(i) <= 0.0) .or. (u(i,1) == 0.0)) then
    l_sfblock(i) = .false.
  end if
  l_cont(i)  = .true.
  l_cont2(i) = .true.
  l_cont3(i) = .true.
end do
!$OMP END DO

! Code copied from gw_setup
! --------------------------
! Estimate N squared at each level - DRY
!$OMP DO SCHEDULE(STATIC)
do k = 2, levels-1
  do i = 1, points
    if (l_sfblock(i)) then
      nsq(i,k) = 2.0*g*  ( theta(i,k) - theta(i,k-1) )                         &
                   /( ( theta(i,k) + theta(i,k-1) )                            &
                     *( theta_levels(i,k) - theta_levels(i,k-1) ) )
    end if!(l_sfblock(i)).
  end do !i = 1, points
end do !k= 2, levels-1
!$OMP END DO

!$OMP DO SCHEDULE(STATIC)
do i = 1, points
  !Set nsq(1)=nsq(2) as nsq is undefined on level 1
  nsq(i,1)       = nsq(i,2)
  !Set nsq(levels)=nsq(levels-1) as nsq is undefined on top level
  nsq(i,levels)  = nsq(i,levels-1)
end do !i = 1, points
!$OMP END DO
!$OMP END PARALLEL

! Get bottom and top of averaging layer (top half mountain)
do k = 1,levels
  do i = 1, points
    if (l_sfblock(i)) then
      if (theta_levels(i,k) <= 0.5*mt_high(i)) then
        kbot(i) = k
      end if
      if (theta_levels(i,k) < mt_high(i)) then
        ktop(i) = k+1
      end if
    end if !(l_sfblock(i)).
  end do !i = 1, points
end do !k = 1,levels

! Find low-level averages (top half mountain)
do k = 2, levels-1
  do i = 1, points
    if (l_sfblock(i)) then
      if ((k > kbot(i)) .and. (k <= ktop(i))) then
        dzt = theta_levels(i,k)   -  0.5 * mt_high(i)
        dzb = theta_levels(i,k)   -  theta_levels(i,k-1)
        if (k  ==  kbot(i)+1) then
          dzu = 0
          dzb = dzt
        else
          dzu = theta_levels(i,k-1) - 0.5 * mt_high(i)
        end if ! (k == bot(i)+1)
        if (k == ktop(i)) then
          if (k  ==  kbot(i)+1) then
            dzb = dzt
          else
            dzt = mt_high(i) -  0.5 * mt_high(i)
            dzb = mt_high(i) -  theta_levels(i,k-1)
          end if !(k  ==  kbot(i)+1)
        end if ! (k == ktop(i))
        !------------------------------------------------------
        ! average u,v and rho from z = 0.5h to current level
        !-----------------------------------------------------
        ulow(i)  = (ulow(i)   * dzu  + u(i,k)  *dzb)/dzt
        vlow(i)  = (vlow(i)   * dzu  + v(i,k)  *dzb)/dzt
        ! bulk average N from z = 0.5mt_high tO z = mt_high
        nlow(i)  = (nlow(i)   * dzu  + nsq(i,k)*dzb)/dzt
      end if ! (k > kbot(i) .and k <= ktop(i))
    end if !(l_sfblock(i))
  end do !i = 1, points
end do  ! k = 2, levels-1

do i = 1, points
  if (l_sfblock(i)) then

    modu(i)    = sqrt(ulow(i)**2 + vlow(i)**2)

    if ( modu(i) <= 0.0 ) then
      l_sfblock(i) = .false.
      zb(i) = mt_high(i)
      hmteff(i) = 0.0
      fav(i) = -2
    end if

    if (nlow(i) > 0.0) then
      nlow(i) = sqrt(nlow(i))
    else
      nlow(i)   = 0.0
      l_sfblock(i) = .false.
    end if !(nlow(i) > 0.)

  end if !(l_sfblock(i)).
end do !i=1,points


! Code from gw_block to get Froude and blocking
! ---------------------------------------------

! ------------------------------------------------
!  calculate zav following Vosper et al (2009)
! ------------------------------------------------
! calculate depth of neutral layer - zneu
do k = 1, levels
  do i = 1, points
    if (l_sfblock(i)) then
      if (l_cont(i)) then
        if (nsq(i,k) < nsq_neutral) then
          zneu(i)    = rho_levels(i,k)
        else
          l_cont(i)  = .false.
        end if !nsq(i,k) < nsq_neutral
      end if !l_cont
    end if !l_sfblock(i)
  end do!i = 1, points
end do !k = 1, levels

! find max(mt_high,zneu) for zav
do i = 1, points
  if (l_sfblock(i)) then
    if (zneu(i) > mt_high(i)) then
      zav1(i) = zneu(i)
    else
      zav1(i) = mt_high(i)
    end if
    zav(i)     = zav1(i)
  end if!l_sfblock(i)
end do !i = 1, points

do ii = 1, zav_iterate
  !Reset l_cont3 (logical to test if rho_levels(k) lt zav)
  do i = 1, points
    l_cont3(i) = .true.
  end do !i = 1, points
  do k = 2, levels-1
    do i = 1, points
      if (l_sfblock(i)) then
        if (l_cont2(i)) then !l_cont2(i) tests if zav is converged
          !-----------------------------------------------------
          ! need an if test for case where zav does not
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
      end if !l_sfblock
    end do ! i=1, points
  end do   ! loop over k levels

  do i = 1, points
    if (l_sfblock(i)) then
      if (l_cont2(i)) then !l_cont2(i) tests if zav is converged
        ! resolve wind in the direction of the low-level flow
        wind    = (ulow(i)*uav(i) + vlow(i)*vav(i)) /modu(i)
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
        zav_new = zav1(i) + u_n(i)
        ! currently set zav_converge to 0.05
        if ((zav(i) < zav_new*(1.0+zav_converge)) .and.                        &
           ( zav(i) > zav_new*(1.0-zav_converge))) then
          l_cont2(i)  =  .false.
        end if ! test if zav(i) is converged
        zav(i)    = zav_new
      end if  ! l_cont2(i)
    end if!l_sfblock(i)
  end do!i = 1, points
end do  !ii= 1, zav_iterate

! -----------------------------------------
!  calculate fav and blocked layer depth
! -----------------------------------------
do i = 1, points
  l_cont(i)  =  .true.
  if (l_sfblock(i)) then
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
      hmteff(i) = mt_high(i) - zb(i)
    else
      zb(i) = 0.0
      hmteff(i) = mt_high(i)
    end if
  end if!l_sfblock(i)

end do !i = 1, points

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_froude_moist
end module lsp_froude_moist_mod
