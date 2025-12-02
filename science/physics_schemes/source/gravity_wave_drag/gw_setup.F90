! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
module gw_setup_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'GW_SETUP_MOD'
contains

subroutine gw_setup(levels,points,gw_seg_size,u,v,rho,theta,                   &
                    ! Inputs to calculate moist buoyancy frequency
                    temp,q,qcl,qcf,press,                                      &
                    ! Outputs from moist buoyancy frequency calculations
                    nsq,nsq_dry,nsq_unsat,nsq_sat,                             &
                    dzcond,l_lapse,kbot,                                       &
                    ulow,vlow,rholow,nlow,psilow,psi1,modu,                    &
                    theta_levels,                                              &
                    sd_orog,grad_xx,grad_xy,grad_yy,mt_high,                   &
                    slope,anis,banis,canis,mtdir,ktop,l_drag,                  &
                    nsigma,                                                    &
!diagnostics
                    u_s_d,u_s_d_on,points_u_s_d,                               &
                    v_s_d,v_s_d_on,points_v_s_d)

! subroutine gw_setup to set up variables for gravity wave and
!            flow blocking schemes

!Code to set up variables for gravity wave and flow blocking schemes
!Calculates
!(i)   Buoyancy frequency squared
!(ii)  Low-level average variables (for FB and GW schemes)
!(iii) SSO variables
!(iv)  ktop
use planet_constants_mod, only: g

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use um_types, only: real_32

! Constants needed for calculation of moist buoyancy frequency
use planet_constants_mod, only: g, repsilon, r, cp, rv
use lsprec_mod,           only: zerodegc, qcfmin, lc
use g_wave_input_mod,     only: i_moist

implicit none
! Description:
!     code to set up variables for gravity wave and flow blocking schemes
!     calculates:
!     1. Buoyancy frequency squared
!     2. Low-level average variables (for FB and GW schemes)
!     3. Sub-gridscale orography (SSO) variables
!     4. First model level above mountain top
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

integer, intent(in) ::  &!intent(in)
  points_u_s_d,                                                                &
  points_v_s_d

real(kind=real_umphys) ::                                                      &
                     !,intent(in):
 nsigma              ! Scaling factor for sd_orog_land to
                     ! give estimate of sub-grid mountain tops

real(kind=real_umphys), intent(in) ::           &!block for intent(in)
  u(gw_seg_size,levels),           &!zonal wind
  v(gw_seg_size,levels),           &!meridional wind
  rho(gw_seg_size,levels),         &!density
  theta(gw_seg_size,levels),       &!potential temperature
  theta_levels(gw_seg_size,levels),&!height(theta_levels)
  sd_orog(points),                 &!standard deviation of orography (m)
  grad_xx(points),                 &!dh/dx squared gradient orography
  grad_xy(points),                 &!(dh/dx)(dh/dy) gradient orography
  grad_yy(points)                   !dh/dy squared gradient orography

real(kind=real_umphys), intent(in) ::                                          &
  ! Inputs for calculating moist buoyancy frequency
  temp(gw_seg_size,levels),        &! temperature
  press(gw_seg_size,levels),       &! pressure
  q(gw_seg_size,levels),           &! water vapour mixing ratio
  qcl(gw_seg_size,levels),         &! resolved cloud water mixing ratio
  qcf(gw_seg_size,levels)           ! resolved cloud ice mixing ratio

logical, intent(in) ::        &!intent(in)
  u_s_d_on,                                                                    &
  v_s_d_on

!-------------------------------------------------------------
! Intent(out) variables (kbot for recalculating nlow in gw_wave)
!------------------------------------------------------------
integer, intent(out) ::     &!block for intent(out)
  ktop(points),             &!First model level above mountain top
  kbot(points)               !Model level at half mountain height

real(kind=real_umphys), intent(out) ::  &!block for intent(out)
  mt_high(points),               &!sso height (n_sigma*sd_orog)
  slope(points),                 &!sso params - slope
  anis(points),                  &!sso params - anisotropy
  mtdir(points),                 &!sso params - angle of major axis
  banis(points),                 &!const (function of anisotropy)
  canis(points),                 &!const (function of anisotropy)
  nsq(gw_seg_size,levels),       &!moist buoyancy freqency squared
  nsq_dry(gw_seg_size,levels),   &!dry buoyancy frequency squared
  nsq_unsat(gw_seg_size,levels), &!unsaturated buoyancy frequency squared
  nsq_sat(gw_seg_size,levels),   &!saturated buoyancy frequency squared
  dzcond(gw_seg_size,levels),    &!Lifting condensation level
                                  ! or descent to evaporate qcl
  ulow(points),                  &!u averaged from z=0.5mt_high to z=mt_high
  vlow(points),                  &!v averaged from z=0.5mt_high to z=mt_high
  modu(points),                & !modulus of horizontal wind (i.e. sqrt(u^2+v^2)
  nlow(points),                & !N bulk averaged from z=0.5mt_high to z=mt_high
  psilow(points),                &!psi averaged from z=0.5mt_high to z=mt_high
  rholow(points),                &!rho averaged from z=0.5mt_high to z=mt_high
  psi1(points)                    !atan(vlow/ulow)

real(kind=real_umphys), intent(out)  ::         &!intent(out) diagnostics
  u_s_d(points_u_s_d),      &!u averaged between 0.5*mt_high and mt_high
  v_s_d(points_v_s_d)        !v averaged between 0.5*mt_high and mt_high

logical, intent(out)  ::                                                       &
   l_drag(points)            !whether point has a non-zero stress or not

! Logical for moist NSQ calculation
logical, intent(out)  ::                                                       &
   l_lapse(gw_seg_size,levels)    !whether latent heating has any effect on the
                                  !buoyancy frequency

!------------------------
! Diagnostic variables
!- ----------------------

!----------------------
! Local variables
!----------------------
integer :: i,k
real(kind=real_umphys)    :: smallp !Small positive number
real(kind=real_32)  :: real4  !dummy variable real*4
real(kind=real_umphys)    :: lmk,lml,lmm,&!Used to calculate the sso params
           plm,mlm      !Used to calculate the sso params
real(kind=real_umphys)    ::                                                   &
     dzt,                                                                      &
     dzb,                                                                      &
     dzu

! Local variables for moist buoyancy frequency calculation
real  ::                                                                       &
  qtotal(points,levels),    & ! Total water mixing ratio
  drydz(points,levels),     & ! Subsaturated subgrid vertical displacement
  satdz(points,levels),     & ! Saturated sub-grid vertical displacement
  thetav(points,levels),    & ! Virtual potential temp
  salr(points,levels)         ! Saturated adiabatic lapse rate

real  ::                                                                       &
   esat,                    & ! Saturation vapour pressure
   qsat,                    & ! Saturation vapour mixing ratio
   tc,                      & ! Temperature in degrees C
   tcloud,                  & ! Temperature at cloud base
   tdew,                    & ! Dewpoint temperature
   tdewc,                   & ! Dewpoint temperature in degrees C
   dewlr,                   & ! Dewpoint vertical lapse rate
   dqldz,                   & ! Adiabatic rate of liquid water formation
   dalr                       ! Dry adiabatic lapse rate

real,parameter ::                                                              &
!  Constants in dewpoint and saturation vapour pressure equations as
!  suggested by Alduchov and Eskridge (1996)
   a1 = 17.625,                                                                &
   b1 = 243.04,                                                                &
!  Multiple in saturation vapour pressure estimate (Pa)
   c1 = 610.94,                                                                &
!  Constant used in virtual potential temperature calculation
   conthetav = 0.61

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='GW_SETUP'
! ----------------------------------------
! initialise arrays
! ----------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

smallp = epsilon(real4) / 100.0
do k = 1, levels
  do i =1, points
    nsq(i,k)=0.0
    nsq_dry(i,k)=0.0
  end do
end do
if (u_s_d_on) then
  do i=1,points
    u_s_d(i)   = 0.0
  end do
end if
if (v_s_d_on) then
  do i=1,points
    v_s_d(i)   = 0.0
  end do
end if
do i=1,points
  ulow(i)    = 0.0
  vlow(i)    = 0.0
  rholow(i)  = 0.0
  psilow(i)  = 0.0
  psi1(i)    = 0.0
  nlow(i)    = 0.0
  anis(i)    = 1.0
  banis(i)   = 0.0
  canis(i)   = 0.0
  slope(i)   = 0.0
  mtdir(i)   = 0.0
  ktop(i)    = 2
  kbot(i)    = 1
  mt_high(i) = nsigma*sd_orog(i)
  l_drag(i)  = .true.
  !u(i,1)=0 implies at north/south pole in global model
  if ((mt_high(i) <= 0.0) .or. (u(i,1) == 0.0)) then
    l_drag(i) = .false.
  end if
end do

! Initialise variables for dry and moist buoyancy frequency calculations
if (i_moist == 1 .or. i_moist==2) then
  do k = 1, levels
    do i =1, points
      nsq_unsat(i,k)=0.0
      nsq_sat(i,k)=0.0
      ! Ascent to reach saturation (or descent to evaporate qcl)
      dzcond(i,k)=0.0
      satdz(i,k)=0.0
      drydz(i,k)=0.0
      ! Total water content
      qtotal(i,k) = q(i,k) + qcl(i,k) + qcf(i,k)
      ! Virtual potential temp
      thetav(i,k) = theta(i,k) * (1.0 + (conthetav*q(i,k))                     &
                    - qcl(i,k) - qcf(i,k))
      ! Saturated adiabatic lapse rate
      salr(i,k) = 0.0
      ! Logical indicating that saturated adiabatic lapse rate is smaller than
      ! the dry adiabatic lapse rate
      l_lapse(i,k)  = .false.
    end do
  end do
  esat = 0.0
  qsat = 0.0
  tc = 0.0
  tdew = 0.0
  tdewc = 0.0
  dewlr = 0.0
  tcloud = 0.0
  dqldz = 0.0
  ! Set Dry adiabatic lapse rate
  dalr = g / cp
end if ! if i_moist


! ----------------------------------------
!  calculate sso parameters and constants
! ----------------------------------------
do i=1,points
  if (l_drag(i)) then
    lmk          = 0.5*(grad_xx(i)+grad_yy(i))
    lml          = 0.5*(grad_xx(i)-grad_yy(i))
    lmm          = grad_xy(i)
    !Denominator of equation (1) in technical documentation
    plm          = lmk + sqrt(lml**2+lmm**2)
    !Numerator of equation (1) in technical documentation
    mlm          = abs(lmk - sqrt(lml**2+lmm**2))
    !Note do not take atan(slope) as slope is always used as tan(slope)
    !in scheme (i.e. keep slope as non-dim number rather than convert
    !to angle
    !Equation (3) in technical documentation
    slope(i)     = sqrt(plm)
    if ((slope(i) == 0.0) ) then
      l_drag(i) = .false.
    end if
    if ((abs(lmm) <= smallp) .and. (abs(lml) <= smallp)) then
      mtdir(i) = 0.0
    else
      !Use atan2 to give correct range, as need mtdir to be in range
      ![-pi/2,pi/2] (using atan would give range [-pi/4,pi/4]
      !Equation (2) in technical documentation
      mtdir(i) = 0.5*atan2(lmm,lml)
    end if

    if (plm <= smallp) then
      anis(i) = 1.0
    else
      !Equation (1) in technical documentation
      anis(i) = sqrt(mlm/plm)
    end if

    banis(i)    = 1-0.18*anis(i) - 0.04*anis(i)**2
    canis(i)    =   0.48*anis(i) + 0.3 *anis(i)**2

  end if !(l_drag(i))
end do ! i=1 ,points

do k = 2, levels-1
  do i = 1, points
    if (l_drag(i)) then
      nsq_dry(i,k) = 2.0*g*( theta(i,k) - theta(i,k-1) )                       &
                     /( ( theta(i,k) + theta(i,k-1) )                          &
                     *( theta_levels(i,k) - theta_levels(i,k-1) ) )
    end if !(l_drag(i))
  end do !i = 1, points
end do !k= 2, levels-1

do i = 1, points
  !Set nsq(1)=nsq(2) as nsq is undefined on level 1
  nsq_dry(i,1)       = nsq_dry(i,2)
  !Set nsq(levels)=nsq(levels-1) as nsq is undefined on top level
  nsq_dry(i,levels)  = nsq_dry(i,levels-1)
end do !i = 1, points

! Moist NSQ calculations
if (i_moist == 1 .or. i_moist==2) then

  !Calculate subsaturated and saturated NSQ profiles

  ! l_lapse indicates that the saturated adiabatic lapse rate is smaller than
  ! the dry adiabatic lapse rate so that latent heat is important
  do k = 1, levels
    do i = 1, points
      if (l_drag(i)) then

        ! Saturated Adiabatic Lapse Rate at this level
        salr(i,k) = g * ( (r*temp(i,k)**2) + (lc*q(i,k)*temp(i,k)) ) /         &
               ( (cp*r*temp(i,k)**2) + (lc**2*q(i,k)*repsilon) )

        if ( (100.0*(dalr - salr(i,k))/salr(i,k)) > 0.3 ) then
          l_lapse(i,k) = .true.
        end if

      end if  ! if l_drag
    end do ! i
  end do ! k

  ! Calculate ascent required to reach saturation and thus
  ! the amount of dry/saturated subgrid orographic ascent to mt_high
  ! Test for q>0 avoids occasional occurence producing NaNs
  do k = 1, levels
    do i = 1, points
      if (l_drag(i) .and. l_lapse(i,k) .and. q(i,k) > 0.0) then

        ! Convert temperature to C for use in tdew
        tc = temp(i,k) - zerodegc

        ! Get saturation vapour mixing ratio for RH
        esat = c1 * exp( a1*tc/(tc + b1) )
        qsat = repsilon * esat / press(i,k)

        if ( (q(i,k)/qsat) < 1.0 ) then  ! Sub-saturated air

          ! Dewpoint temperature (C) then convert to K
          tdewc = b1*( log(q(i,k)/qsat ) + (a1*tc/(b1 + tc)) ) /               &
            ( a1 - log(q(i,k)/qsat ) - (a1*tc/(b1 + tc)) )
          tdew = tdewc + zerodegc

          ! Dewpoint temperature adiabatic lapse rate
          dewlr = (tdew * tdew * g * rv)/(lc * r * temp(i,k))

          ! Ascent to condensation level (ensure positive)
          dzcond(i,k) = max( (temp(i,k) - tdew)/(dalr - dewlr), 0.0)

          ! Assume subgrid orography rises mt_high above model surface
          ! as in the Froude calculation in gw_block
          ! ...Saturated ascent
          satdz(i,k) = max( (mt_high(i) - dzcond(i,k)) , 0.0 )
          ! ...Sub-saturated ascent
          drydz(i,k) = min( dzcond(i,k), mt_high(i) )

        else    ! Air already saturated

          ! All subgrid ascent is saturated
          satdz(i,k) = mt_high(i)
          drydz(i,k) = 0.0

          ! dzcond is now the descent required to evaporate qcl
          ! ...will be needed to recalculate moist NSQ for GWs
          if ( qcl(i,k) > qcfmin ) then

            ! Adiabatic rate of water evaporation per metre of descent
            dqldz = max( 0.0, ( (salr(i,k)*(repsilon + qsat)*qsat*lc)/         &
                      (r*temp(i,k)*2)  )  -                                    &
              ( (qsat*press(i,k)*g)/((press(i,k) - esat)*(r*temp(i,k)) ) )  )

            ! Descent (negative) required to evaporate resolved cloud
            dzcond(i,k) = -qcl(i,k) / dqldz

          else    ! Just saturated with very little water
            dzcond(i,k) = 0.0
          end if

        end if  !  Saturation test

      end if ! if l_drag, l_lapse, q > 0
    end do  ! i
  end do ! k

  ! Estimate various moist buoyancy frequency profiles
  do k = 2, levels-1
    do i = 1, points
      if (l_drag(i)) then

        ! Unsaturated moist buoyancy frequency
        nsq_unsat(i,k) = 2.0*g*  ( thetav(i,k) - thetav(i,k-1) )               &
                     /( ( thetav(i,k) + thetav(i,k-1) )                        &
                     *( theta_levels(i,k) - theta_levels(i,k-1) ) )
        ! Dont allow this to be larger than the dry value (unphysical)
        if ( nsq_unsat(i,k) > nsq_dry(i,k) ) then
          nsq_unsat(i,k) = nsq_dry(i,k)
        end if

        ! Moisture only important if
        ! 1) the lapse rate is modified by latent heat release and
        ! 2) some of the sub-grid orographic ascent is saturated, and
        ! 3) moisture is present (occasionally it is not, producing errors)
        if (l_lapse(i,k) .and. satdz(i,k) > 0.0 .and. q(i,k) > 0.0) then

          ! Find cloud base temperature (q unchanged by dry adiabatic ascent)
          if ( dzcond(i,k) > 0.0 ) then
            tcloud  = temp(i,k)   - ( dalr * dzcond(i,k) )
          else
            tcloud  = temp(i,k)
          end if

          ! Saturated Adiabatic Lapse Rate of this parcel when it
          ! has ascended (dry adiabatically) to cloud base
          salr(i,k) = g * ( (r*tcloud**2) + (lc*q(i,k)*tcloud) ) /             &
                ( (cp*r*tcloud**2) + (lc**2*q(i,k)*repsilon) )

          ! Calculate saturated NSQ profile valid above dzcond
          nsq_sat(i,k) = salr(i,k) + ( (temp(i,k) - temp(i,k-1)) /             &
                         (theta_levels(i,k) - theta_levels(i,k-1)) )
          nsq_sat(i,k) = nsq_sat(i,k)*( (g*2.0) / (temp(i,k) + temp(i,k-1)) )

          nsq_sat(i,k) = nsq_sat(i,k)*( 1.0 + (lc*q(i,k))/(r*tcloud) )

          nsq_sat(i,k) = nsq_sat(i,k) - (g/(1.0 + qtotal(i,k))) *              &
                           (qtotal(i,k) - qtotal(i,k-1)) /                     &
                           (theta_levels(i,k) - theta_levels(i,k-1))

          ! Dont allow this to be larger than the dry value (unphysical)
          if ( nsq_sat(i,k) > nsq_unsat(i,k) ) then
            nsq_sat(i,k) = nsq_unsat(i,k)
          end if

        else  ! No need to account for latent heating

          nsq_sat(i,k) = nsq_unsat(i,k)

        end if  ! if l_lapse, satdz>0, q>0

      end if ! if l_drag
    end do  ! i
  end do ! k = 2 to levels-1

  do i = 1, points
    !Set nsq(1)=nsq(2) as nsq is undefined on level 1
    nsq_unsat(i,1)       = nsq_unsat(i,2)
    nsq_sat(i,1)       = nsq_sat(i,2)
    !Set nsq(levels)=nsq(levels-1) as nsq is undefined on top level
    nsq_unsat(i,levels)  = nsq_unsat(i,levels-1)
    nsq_sat(i,levels)  = nsq_sat(i,levels-1)
  end do !i = 1, points

  do k = 1, levels
    do i = 1, points
      if (l_drag(i)) then
        if ( l_lapse(i,k) .and. satdz(i,k) > 0.0 .and. q(i,k) > 0.0 ) then
          nsq(i,k) = ((nsq_unsat(i,k)*drydz(i,k)) +                            &
                    (nsq_sat(i,k)*satdz(i,k)))                                 &
                        / mt_high(i)
        else
          nsq(i,k) = nsq_unsat(i,k)
        end if  ! If l_lapse
      end if ! if l_drag
    end do  ! i
  end do ! k = 1 to levels


else    ! Dry NSQ

  do k = 1, levels
    do i = 1, points
      nsq(i,k) = nsq_dry(i,k)
    end do  ! i
  end do ! k = 1 to levels

end if  ! i_moist

! Start low-level averaging
do k = 1,levels
  do i = 1, points
    if (l_drag(i)) then
      if (theta_levels(i,k) <= 0.5*mt_high(i)) then
        kbot(i) = k
      end if
      if (theta_levels(i,k) < mt_high(i)) then
        ktop(i) = k+1
      end if
    end if !(l_drag(i))
  end do !i = 1, points
end do !k = 1,levels

do k = 2, levels-1
  do i = 1, points
    if (l_drag(i)) then
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
        rholow(i)= (rholow(i) * dzu  + rho(i,k)*dzb)/dzt

        ! bulk average N from z = 0.5mt_high tO z = mt_high
        nlow(i)  = (nlow(i)   * dzu  + nsq(i,k)*dzb)/dzt

      end if ! (k > kbot(i) .and k <= ktop(i))
    end if !(l_drag(i))
  end do !i = 1, points
end do  ! k = 2, levels-1

do i = 1, points
  if (l_drag(i)) then

    modu(i)    = sqrt(ulow(i)**2 + vlow(i)**2)

    if ( modu(i) <= 0.0 ) then
      l_drag(i) = .false.
    end if

    psi1(i)    = atan2(vlow(i),ulow(i))
    psilow(i)  = mtdir(i) - psi1(i)
    if (nlow(i) > 0.0) then
      nlow(i) = sqrt(nlow(i))
    else
      nlow(i)   = 0.0
      l_drag(i) = .false.
    end if !(nlow(i) > 0.)
  end if !(l_drag(i))
end do !i=1,points
!-----------------------------------------------------------------
! 4 diagnostics
!-----------------------------------------------------------------
if ( u_s_d_on ) then
  do i=1,points
    u_s_d(i) = ulow(i)
  end do
end if

if ( v_s_d_on ) then
  do i=1,points
    v_s_d(i) = vlow(i)
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine gw_setup
end module gw_setup_mod
