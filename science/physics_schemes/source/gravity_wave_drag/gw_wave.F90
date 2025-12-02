! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module gw_wave_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'GW_WAVE_MOD'
contains

subroutine gw_wave(levels,points,gw_seg_size,u,v,rho,                          &
                   ! Inputs for recalculating moist buoyancy frequency profile
                   nsq_dry,nsq_unsat,nsq_sat,dzcond,l_lapse,                   &
                   kbot,                                                       &
                   ulow,vlow,rholow,psi1,                                      &
                   psilow,nlow,modu,ktop,rho_levels,theta_levels,              &
                   delta_lambda,delta_phi,latitude,mt_high,sd_orog,            &
                   slope,zb,banis,canis,orog_f1,orog_f2,orog_f3,orog_amp,      &
                   dudt,dvdt,dtdt,                                             &
                   l_dynbeta,l_nonhydro,l_smooth,fsat,Gsharp,                  &
                   l_drag,l_gw_heating,                                        &
                   ! diagnostics
                   du_dt_satn,points_du_dt_satn,du_dt_satn_on,                 &
                   du_dt_satn_p_on,                                            &
                   dv_dt_satn,points_dv_dt_satn,dv_dt_satn_on,                 &
                   stress_ud,points_stress_ud ,stress_ud_on,                   &
                   stress_ud_p_on,                                             &
                   stress_vd,points_stress_vd ,stress_vd_on,                   &
                   stress_ud_satn,points_stress_ud_satn,                       &
                   stress_ud_satn_on,                                          &
                   stress_vd_satn,points_stress_vd_satn,                       &
                   stress_vd_satn_on,                                          &
                   tausx_d, tausx_d_on   , points_tausx_d  ,                   &
                   tausy_d, tausy_d_on   , points_tausy_d,                     &
                   nsq_s_d,nsq_s_d_on,points_nsq_s_d)

! subroutine gw_wave to calculate the vertical profile of gravity wave
!            stress and associated wind increments

! Constants needed for calculation of moist NSQ
use planet_constants_mod, only: planet_radius, cp
use conversions_mod, only: pi
use g_wave_input_mod,     only: i_moist, scale_aware, middle, var

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use c_gwave_mod, only: beta_fix, frac_wl, lambdaz_min, lambdaz_max,            &
                       nsq_neutral

implicit none

! Description:
!     calculates the gwd stress profiles and wind increments.
!     1. calculate stress profile and wind increments for
!        linear hydrostatic waves.
!     2. stress may be deposited either over a single model
!        level or over a vertical gravity wave wavelength.
!     3. non-hydrostatic waves may be allowed to propagate outside
!        of the column.
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
integer, intent(in) ::        &!block for intent(in)
  levels,                     &!num of vertical levels
  points,                     &!num of land points
  gw_seg_size

integer, intent(in) ::  &!intent(in)
  points_nsq_s_d

integer, intent(in) ::        &!block for intent(in)
  ktop(points),               &!half mountain top height
  kbot(points)                 !first model level above mountain top

integer, intent(in) ::        &! intent(in) for diags
  points_stress_ud,                                                            &
  points_stress_vd,                                                            &
  points_stress_ud_satn,                                                       &
  points_stress_vd_satn,                                                       &
  points_du_dt_satn,                                                           &
  points_dv_dt_satn,                                                           &
  points_tausx_d,                                                              &
  points_tausy_d

real(kind=real_umphys), intent(in)    ::        &!block for intent(in)
  u(gw_seg_size,levels),           &!zonal wind
  v(gw_seg_size,levels),           &!meridional wind
  rho(gw_seg_size,levels),         &!density
  rho_levels(gw_seg_size,levels),  &!height(rho_levels)
  theta_levels(gw_seg_size,levels),&!height (theta_levels)
  slope(points),                   &!sso params - slope
  mt_high(points),                 &!sso params - height
  sd_orog(points),                 &!standard deviation of orography (m)
  banis(points),                   &!const (function of anisotropy)
  canis(points),                   &!const (function of anisotropy)
  orog_f1(points),                 &!orographic F1 coefficient
  orog_f2(points),                 &!orographic F2 coefficient
  orog_f3(points),                 &!orographic F3 coefficient
  orog_amp(points),                &!orographic amplitude (max - min)
  ulow(points),                    &!u averaged from z=0.5mt_high to z=mt_high
  vlow(points),                    &!v averaged from z=0.5mt_high to z=mt_high
  modu(points),                    &!modulus of horizontal wind ( sqrt(u^2+v^2))
  psilow(points),                  &!psi averaged from z=0.5mt_high to z=mt_high
  rholow(points),                  &!rho averaged from z=0.5mt_high to z=mt_high
  psi1(points),                    &!atan(vlow/ulow)
  zb(points),                      &!depth of flow blocking layer
  latitude(points),                &!latitude
  delta_lambda,                    &!spacing between points in the i direction.
  delta_phi,                       &!spacing between points in the j direction.
  fsat,                            &!saturation Froude number
  Gsharp                            !function of mtn sharpness

real(kind=real_umphys), intent(in)    ::& !intent(in) for moist NSQ calculations
  nsq_dry(gw_seg_size,levels),     &!Dry moist Brunt-Vaisala freq squared
  nsq_unsat(gw_seg_size,levels), & !Unsaturated moist Brunt-Vaisala freq squared
  nsq_sat(gw_seg_size,levels),     &!Saturated moist Brunt-Vaisala freq squared
  dzcond(gw_seg_size,levels)        !Displacement to LCL

logical, intent(in) ::  &!intent(in)
  nsq_s_d_on

logical, intent(in) ::      &!intent(in) for moist NSQ calculations
  l_lapse(gw_seg_size,levels) !Indicates whether Latent Heating has effect on N

logical, intent(in) ::                                                         &
  l_dynbeta,                &!dynamically adjusting beta(angle of group vel)
  l_nonhydro,               &!nonhydro scheme
  l_smooth,                 &!lambda_z smoothing of acc
  l_gw_heating               !calculate heating tendency if true

logical, intent(in)  ::      &!intent(in) for diags
  stress_ud_on,              &!u stress
  stress_ud_p_on,            &!u stress on pressure points
  stress_vd_on,              &!v stress
  stress_ud_satn_on,         &!u satn stress
  stress_vd_satn_on,         &!v satn stress
  du_dt_satn_on,             &!u accel (saturation)
  du_dt_satn_p_on,           &!u accel (saturation) on pressure points
  dv_dt_satn_on,             &!v accel (saturation)
  tausx_d_on,                &!tausx_d switch
  tausy_d_on                  !tausy_d switch
!----------------------
! Intent(in out) variables
!----------------------
real(kind=real_umphys), intent(in out) ::       &!block for intent(in out)
  dudt(gw_seg_size,levels),                                                    &
  dvdt(gw_seg_size,levels),       &!profiles of acceleration
  dtdt(gw_seg_size,levels)         !profiles of heating

real(kind=real_umphys),intent(in out)  ::           &!intent(inout) diags
  stress_ud (points_stress_ud,0:levels),                                       &
  stress_vd (points_stress_vd,0:levels),                                       &
  tausx_d(points_tausx_d),                                                     &
  tausy_d(points_tausy_d)


! Variables for moist buoyancy frequency calculation
real(kind=real_umphys), intent(in out) ::                                      &
  nlow(points)               !N bulk averaged from 0.5mt_high to mt_high

logical, intent(in out) ::                                                     &
  l_drag(points)             !whether point has a non-zero stress or not

!----------------------
! Intent(out) variables
!----------------------
real(kind=real_umphys),intent(out)    ::           &!intent(out) diags
  stress_ud_satn (points_stress_ud_satn, 0:levels),                            &
  stress_vd_satn (points_stress_vd_satn, 0:levels),                            &
  du_dt_satn (points_du_dt_satn,levels),                                       &
  dv_dt_satn (points_dv_dt_satn,levels)

real(kind=real_umphys), intent(out)  ::         &!intent(out) diags
  nsq_s_d(points_nsq_s_d)      !0.5mt_high-mt_high av n (nlow) diag

!----------------------
! Local variables
!----------------------
real(kind=real_umphys) ::&!block for local variables
 wave_amp(points,levels),&!Profile of wave amplitude
 lambdaz(points,levels), &!Vertical wavelength
 tau(points,levels),     &!stress profile
 beta(points,levels),    &!ratio of vert group vel to hrz group vel
 nonhyd(points,levels),  &!used to calc frac of stress leaving grid-column
 k_wave(points),         &!horizontal wavenumber
 heff(points),           &!cut-off mtn height (h-zb)
 pd1(points),            &!Directional term
 pd2(points),            &!Directional term
 pdmod(points),          &!Directional term
 tau_sfc(points),        &!used in calc of sfc stress
 taux_sfc(points),       &!used in calc of sfc stress
 tauy_sfc(points),       &!used in calc of sfc stress
 spd(points,levels),     &!wind resolved in direction of low-level stress
 n(points,levels),       &!buoyancy freq sqrt(nsq)
 deltaz(points,levels),  &!z vertical grid-length
 u_n(points,levels),     &!wind speed div by n
 maxgwdlev,              &!max height for gwd
 deltax(points,levels),  &!x horizontal grid-length
 deltay(points,levels),  &!y horizontal grid-length
 uacc,                   &!dudt at level k (prior to spreading over lambdaz)
 vacc,                   &!dvdt at level k (prior to spreading over lambdaz)
 rhoav,                  &!Average density over lambdaz
 rhob,                   &!rho on theta grid at level k
 rhob_l,                 &!rho on theta grid at level k-1
 ub,                     &!u on theta grid at level k
 vb,                     &!v on theta grid at level k
 spsi,                   &!(sin(psi))**2
 scpsi,                  &!(sin(psi))*(cos(psi))
 cpsi1_dzdy,             &!cos(psi1)*deltaz(k)*deltay
 spsi1_dzdx,             &!sin(psi1)*deltaz(k)*deltax
 ztau,                   &!sfc stress(before directional terms)
 zvt1,zvt2,                                                                    &
 ztemp,                                                                        &
 amplow,                 &!Wave amplitude at level k-1
 acrit,                  &!Critical wave amplitude
 dzb,dzu,dzl,            &!layer thickness used for various averages
 uhat,                   &!u after gravity wave drag increment
 vhat,                   &!v after gravity wave drag increment
 ududt,                  &!Rate of change of kinetic energy after wind incs
 vdvdt,                  &!Rate of change of kinetic energy after wind incs
 zb_scaled,              &!zb / hamp_orog, normalised blocking depth
 alpha                    !Function of zb_scaled to reduce stress



real(kind=real_umphys) ::&!block for local variables used to calc diags
 dudt_gw(points,levels),      &!u acc due to gravity waves
 dvdt_gw(points,levels),      &!v acc due to gravity waves
 local_xstress(points,0:levels),&!for calculation of stress diagnostics
 local_ystress(points,0:levels),&!for calculation of stress diagnostics
 dz                            !delta z on theta grid for stress diags

integer :: i, k, kl, kk
integer :: khigh(points,levels),klow(points,levels)

logical :: l_cont(points), l_cont2(points)


! Local variables required for moist NSQ calculation
real(kind=real_umphys) ::                                                      &
 nsq(points,levels),     &!Moist buoyancy frequency squared
 nsq_new(points,levels), &!..and temporary value for GW propagation iteration
 n_u(points,levels),     &!unlimited n div by windspeed
 n_u_av(points),         &!Running average of n_u from top of subgrid hill
 avlambdaz,              &!Convert n_u_av to mean GW vertical wavelength
 maxup,                  &!Maximum GW ascent for final moist NSQ
 maxdown,                &!Maximum GW descent for final moist NSQ
 satdz,                  &!Saturated subgrid orographic ascent
 drydz,                  &!Sub-saturated subgrid orographic ascent
 fracwave,               &!Fraction through current vertical wavelength
 dzt,                    &!theta layer height above ht or ht/2
 dzw,                    &!Altitude above previous lambdaz/2
 nsq_converge             !Minimum NSQ change for convergence

integer :: n_iterate, ii, nwave

logical :: l_cont3(points,levels)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='GW_WAVE'

!---------------------------------------------------------------------
! 1.0 initialisation
!---------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
if (stress_ud_satn_on) then
  do k=0,levels
    do i=1,points
      stress_ud_satn(i,k) = 0.0
    end do
  end do
end if
if (stress_vd_satn_on) then
  do k=0,levels
    do i=1,points
      stress_vd_satn(i,k) = 0.0
    end do
  end do
end if
if (du_dt_satn_on .or. du_dt_satn_p_on) then
  do k=1,levels
    do i=1,points
      du_dt_satn(i,k) = 0.0
    end do
  end do
end if
if (dv_dt_satn_on) then
  do k=1,levels
    do i=1,points
      dv_dt_satn(i,k) = 0.0
    end do
  end do
end if
if ( stress_ud_on .or. stress_ud_p_on .or. stress_ud_satn_on                   &
                 .or. tausx_d_on ) then
  do k=0,levels
    do i=1,points
      local_xstress(i,k) =0.0
    end do
  end do
end if
if ( stress_vd_on .or. stress_vd_satn_on .or. tausy_d_on ) then
  do k=0,levels
    do i=1,points
      local_ystress(i,k) =0.0
    end do
  end do
end if
do k=1,levels
  do i=1,points
    tau(i,k)     = 0.0
    nonhyd(i,k)  = 0.0
    dudt_gw(i,k) = 0.0
    dvdt_gw(i,k) = 0.0
  end do
end do
do i = 1, points
  tau_sfc(i) = 0.0
  taux_sfc(i) = 0.0
  tauy_sfc(i) = 0.0
end do

! Initialise NLOW diagnostic
if (nsq_s_d_on) then
  do i=1,points
    nsq_s_d(i) = 0.0
  end do
end if
!
! Set maxgwdlev
!
maxgwdlev = 1.0e10

! Initialise NSQ with dry value
! Reset later if moisture requested
do k = 1, levels
  do i = 1, points
    if (l_drag(i)) then
      nsq(i,k) = nsq_dry(i,k)
    end if
  end do
end do

if ( L_nonhydro ) then
  do k = 1, levels
    do i = 1, points
      deltay(i,k) = (planet_radius+theta_levels(i,k))
      deltax(i,k) = cos(latitude(i))*delta_lambda*deltay(i,k)
      deltay(i,k) = deltay(i,k)*delta_phi
    end do
  end do
end if

do i = 1, points
  if (l_drag(i)) then
    Heff(i)   = mt_high(i)-zb(i)
    scpsi     = sin(psilow(i))*cos(psilow(i))
    spsi      = (sin(psilow(i)))**2
    pd1(i)    = banis(i)-(banis(i)-canis(i))*spsi
    pd2(i)    = (banis(i)-canis(i))*scpsi
    pdmod(i)  = sqrt(pd1(i)**2+pd2(i)**2)
    k_wave(i)  = slope(i)/mt_high(i)  !horizontal wave number = 1/L
  end if
end do


! Recalculate moist Nsq for low-level air rising over Heff
! Only levels below sub-grid mountain top are utilised for NLOW
! (NSQ is recalculated later to levels above mountain top)
! -------------------------------------------------------
if (i_moist == 1 .or. i_moist==2) then

  ! Recalculate moist Nsq for low-level air rising over Heff
  ! for use in new moist NLOW
  ! All levels - gives first estimate of NSQ for GW propagation later
  do k = 1, levels
    do i = 1, points
      if (l_drag(i)) then
        if ( l_lapse(i,k) ) then  ! Only necessary if SALR < DALR
          if (dzcond(i,k) >= 0.0) then
            drydz = min( dzcond(i,k), Heff(i) )
            satdz = max( Heff(i) - drydz, 0.0 )
            nsq(i,k) = ( (nsq_unsat(i,k)*drydz) +                              &
                       (nsq_sat(i,k)*satdz) ) / ( Heff(i) )
          else
            nsq(i,k) = nsq_sat(i,k)
          end if
        else
          nsq(i,k) = nsq_unsat(i,k)
        end if  !  If l_lapse
      end if ! if l_drag
    end do ! i
  end do ! k

  !Re-initialise nlow
  do i = 1, points
    nlow(i) = 0.0
  end do ! i

  ! Recalculate low-level average of new NSQ
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
          ! bulk average N from z = 0.5mt_high tO z = mt_high
          nlow(i)  = (nlow(i) * dzu  + nsq(i,k)*dzb)/dzt
        end if ! (k > kbot(i) .and k <= ktop(i))
      end if !(l_drag(i))
    end do !i = 1, points
  end do  ! k = 2, levels-1

  ! New nlow and switch off drag if <=0
  ! Used in tau_sfc which would be zero or negative
  do i = 1, points
    if (l_drag(i)) then
      if (nlow(i) > 0.0) then
        nlow(i) = sqrt(nlow(i))
      else
        nlow(i)   = 0.0
        l_drag(i) = .false.
      end if !(nlow(i) > 0.)
    end if !(l_drag(i))
  end do !i=1,points

end if ! if i_moist


! Total surface drag calculation
do i = 1, points
  if (l_drag(i)) then
    !-----------------------------------
    ! Calculate surface gravity wave stress
    ! following LM97 (but with cut-off mtn)
    !-----------------------------------
    ! here assume that maximum depth of zb is mt_high
    if (scale_aware) then
      if ((zb(i) > 0.0) .and. (zb(i) < orog_amp(i))                            &
         .and. (zb(i) < mt_high(i))) then
        zb_scaled = zb(i)/orog_amp(i)
        alpha = 0.5*( 1.0 - tanh( (zb_scaled - middle) / var ) )
      else if (zb(i) == 0.0) then
        alpha = 1.0
      else
        alpha = 0.0
      end if !zb(i) > 0.0
      taux_sfc(i)  = rholow(i)*nlow(i)*alpha*Gsharp*                           &
                    (ulow(i)*orog_f1(i) + vlow(i)*orog_f2(i))
      tauy_sfc(i)  = rholow(i)*nlow(i)*alpha*Gsharp*                           &
                    (ulow(i)*orog_f2(i) + vlow(i)*orog_f3(i))
      tau_sfc(i)  = sqrt(taux_sfc(i)**2 + tauy_sfc(i)**2)
      if (tau_sfc(i) == 0.0) then
        l_drag(i) = .false.
      end if
    else
      ztau        = rholow(i)*modu(i)*nlow(i)*(0.25*Heff(i)**2)                &
                  *(slope(i)/sd_orog(i))*Gsharp
      tau_sfc(i)  = ztau*pdmod(i)
    end if ! scale_aware
  end if ! l_drag
end do


! Start vertical GW propagation

!Resolve wind into direction of sfc stress (ULOW, VLOW)
!on theta (theta_levels) grid
do k=1, levels-1
  do i = 1, points
    if (l_drag(i)) then
      dzl          = theta_levels(i,k)    -  rho_levels(i,k)
      dzu          = rho_levels(i,k+1)   -  theta_levels(i,k)
      dzb          = rho_levels(i,k+1)   -   rho_levels(i,k)
      ub           = (dzu*u(i,k)       + dzl*u(i,k+1))/dzb
      vb           = (dzu*v(i,k)       + dzl*v(i,k+1))/dzb
      if (scale_aware) then
        spd(i,k)     = (ub*taux_sfc(i) + vb*tauy_sfc(i))/tau_sfc(i)
      else
        zvt1         =  ulow(i)*ub + vlow(i)*vb
        zvt2         = -vlow(i)*ub + ulow(i)*vb
        spd(i,k)     = (zvt1*pd1(i) + zvt2*pd2(i))/(modu(i)*pdmod(i))
      end if ! scale_aware
    end if ! l_drag
  end do
end do

! GW launch of magnitude Heff below mountain top (k<ktop)
! Uses moist NSQ valid for Heff (only NSQ above this recalculated later)
! or dry NSQ if i_moist=0
do k=1, levels-1
  do i = 1, points
    if (l_drag(i)) then
      if (k < ktop(i)) then
        !set stress and wave amp to zero if encounter a neutral layer
        !or a critical layer below the mountain top so that wave
        !propagation loop doesn't calc for cases where
        !there is a neutral layer or critical layer at k=ktop-1
        if ((nsq(i,k) <= nsq_neutral) .or.                                     &
                 (spd(i,k+1)*spd(i,k) <= 0.0)) then
          tau(i,k)      = 0.0
          wave_amp(i,k) = 0.0
        else
          !otherwise set wave amp and stress as follows
          wave_amp(i,k) = max(0.0,Heff(i))
          tau(i,k) = tau_sfc(i)
        end if
      else if (k == ktop(i)) then
        wave_amp(i,k) = max(0.0,Heff(i))
      else if (k > ktop(i)) then
        wave_amp(i,k) = 0.0
      end if
    end if
  end do
end do

! If moist N is to be used for vertical wave propagation/saturation (which
! determines the levels on which the total wave stress is deposited), then
! iteratively calculate GW amplitude, vertical wavelength and therefore
! moist NSQ at each altitude.
! -----------------------------------------------------------
if (i_moist == 2) then

  n_iterate = 10.0     ! Number of iterations
  nsq_converge = 0.05  ! NSQ convergence limit (5%)

  ! At each successive model level, require an average vertical
  ! GW wavelength avlambdaz over all previous levels
  ! Unlike the calculation of u_n and lambdaz used for beta and vertical
  ! spreading of GW drag, this is not forced to remain within strict limits
  ! Used to determine the wave phase at each level

  ! Initialise convergence identifier
  do k=1, levels-1
    do i = 1, points
      l_cont3(i,k) = .true.  !NSQ not yet converged
    end do
  end do

  ! Start iteration loop for vertical GW proagation and NSQ
  do ii = 1, n_iterate

    ! Find original U_N to give BETA and lambdaz
    do k=1, levels-1
      do i = 1, points
        if (l_drag(i)) then
          if (l_cont3(i,k)) then  !If moist NSQ not yet converged
            if (nsq(i,k) > nsq_neutral) then
              ! first consider stable cases
              n(i,k)   = sqrt(nsq(i,k))
              u_n(i,k) = spd(i,k)/n(i,k)
            else
              ! now consider neutral (and near neutral) cases
              u_n(i,k) = spd(i,k)/sqrt(nsq_neutral)
            end if
            ! limit u_n to sensible values
            ! only used for beta and for vertically spreading GW drag
            u_n(i,k)     = max(lambdaz_min,u_n(i,k))
            u_n(i,k)     = min(lambdaz_max,u_n(i,k))
            lambdaz(i,k) = frac_wl*(2.0*pi*u_n(i,k))
            deltaz(i,k)  = rho_levels(i,k+1)-rho_levels(i,k)
            ! set beta as beta_fix
            beta(i,k)    = beta_fix
          end if
        end if
      end do
    end do

    if (l_dynbeta) then
      do k=1, levels-1
        do i = 1, points
          if (l_drag(i)) then
            if (l_cont3(i,k)) then
              !Calculate dynamically adjusting group velocity
              !angle at each height (Uses limited U_N)
              if (spd(i,k) /= 0.0) then
                if ((k_wave(i)**2*u_n(i,k)**2) < 1.0) then
                  beta(i,k) = sqrt(1.0-k_wave(i)**2*u_n(i,k)**2)               &
                             /(k_wave(i)*u_n(i,k))
                else
                  beta(i,k) = 1.0
                end if
              else
                beta(i,k) = 1.0e16
              end if
            end if !l_cont
          end if !l_drag
        end do ! i
      end do !k
    end if ! if l_dynbeta

    !-----------------------------------------------------------------
    ! Calculate stress lost out of grid-box due to non-hydrostatic prop.
    ! and deposit over height that it is lost
    !-----------------------------------------------------------------
    if (l_nonhydro) then
      do k = 1, levels-1
        do i = 1, points
          if (l_drag(i)) then
            if (l_cont3(i,k)) then  !If moist NSQ not yet converged
              cpsi1_dzdy  = abs(cos(psi1(i)))*deltaz(i,k)*deltay(i,k)
              spsi1_dzdx  = abs(sin(psi1(i)))*deltaz(i,k)*deltax(i,k)
              nonhyd(i,k) = (cpsi1_dzdy + spsi1_dzdx)                          &
                        /(beta(i,k)*deltax(i,k)*deltay(i,k)                    &
                         + cpsi1_dzdy + spsi1_dzdx)
            end if
          end if
        end do
      end do
    end if !l_nonhydro


    !-----------------------------------
    ! Wave propagation following McFarlane
    !-----------------------------------
    ! Start loop at mountain top
    ! Just sfc stress as launch stress from mtn top
    do k=2, levels-1
      do i = 1, points
        if (l_drag(i)) then
          if (l_cont3(i,k)) then  !If moist NSQ not yet converged
            if (k >= ktop(i) .and. theta_levels(i,k) <=maxgwdlev) then
              acrit       = 0.0
              kl          = k-1
              if ((wave_amp(i,kl)  /=  0.0) .and. (tau(i,kl) /= 0.0)) then
                ! All wave stress deposited if dth/dz < 0 or wind more than
                ! pi/2 to direction of surface stress
                if ((nsq(i,k) <= nsq_neutral) .or.                             &
                    (spd(i,k) <= 0.0)) then
                  tau(i,k) = 0.0
                else
                  ! Test whether wave amplitude exceeds critical amplitude,
                  ! acrit, for wave breaking
                  dzl           = theta_levels(i,k)   -  rho_levels(i,k)
                  dzu           = rho_levels(i,k+1)   -  theta_levels(i,k)
                  dzb           = rho_levels(i,k+1)   -   rho_levels(i,k)
                  rhob          = (dzu*rho(i,k)       + dzl*rho(i,k+1))/dzb

                  dzl           = theta_levels(i,k-1) -  rho_levels(i,k-1)
                  dzu           = rho_levels(i,k)     -  theta_levels(i,k-1)
                  dzb           = rho_levels(i,k)     -   rho_levels(i,k-1)
                  rhob_l        = (dzu*rho(i,k-1)     + dzl*rho(i,k))  /dzb
                  amplow        = rhob_l*n(i,kl)*spd(i,kl)
                  wave_amp(i,k) = wave_amp(i,kl)*                              &
                              sqrt(amplow/(rhob*spd(i,k)*n(i,k)))
                  acrit         = fsat*(spd(i,k)/n(i,k))
                  if (wave_amp(i,k) > acrit) then
                    wave_amp(i,k) = acrit
                  end if
                  tau(i,k) = tau(i,kl)*                                        &
                    (wave_amp(i,k) / wave_amp(i,kl))**2*                       &
                    (rhob*n(i,k)*spd(i,k)) / (rhob_l*n(i,kl)*spd(i,kl)) *      &
                    (1.0-nonhyd(i,k))
                end if  !  n(i,k)<= 0 and spd(i,k)<= 0
              end if !((wave_amp(i,kl)  /=  0.0) .and. (tau(i,kl) /= 0.0))
            end if !  (k >= ktop(i) .and. theta_levels(i,k) <=maxgwdlev)
          end if !(l_cont3(i,k))
        end if !(l_drag(i))
      end do !  i = 1, points
    end do !k=1, levels-1

    !Must initialise running average for each iteration
    do i = 1, points
      if (l_drag(i)) then
        n_u_av(i) = 0.0 ! Running average
      end if
    end do
    do k=1, levels-1
      do i = 1, points
        n_u(i,k) = 0.0  ! N_U at each level
      end do
    end do

    ! Recalculate new moist NSQ using GW amplitude and wavelength
    ! Running U_N average must be done ragardless of convergence at that
    ! level as may be needed at higher levels
    do k = 2, levels
      do i = 1, points
        if (l_drag(i)) then
          if (k >= ktop(i) .and. theta_levels(i,k) <=maxgwdlev) then
            if ( (wave_amp(i,k) /=  0.0) .and. l_lapse(i,k) ) then

              ! Unlimited N/wind (vertical wavelength) at this level
              if (nsq(i,k) > nsq_neutral) then
                ! first consider stable cases
                n(i,k)   = sqrt(nsq(i,k))
                n_u(i,k) = n(i,k)/spd(i,k)
              else
                ! now consider neutral (and near neutral) cases
                n_u(i,k) = sqrt(nsq_neutral)/spd(i,k)
              end if

              ! Height above sub-grid hill peak
              dzt = theta_levels(i,k) - mt_high(i)
              ! Vertically average all levels from ktop to present level
              if (k == ktop(i)) then    ! Depth of this level
                dzb = dzt  ! Depth of this level
                dzl = 0.0  ! No previous level used
              else
                ! Level depth
                dzb = theta_levels(i,k) - theta_levels(i,k-1)
                ! Height above SG mountain peak of level below
                dzl = theta_levels(i,k-1) - mt_high(i)
              end if
              n_u_av(i) = (n_u_av(i)*dzl + n_u(i,k)*dzb) / dzt

              ! Average vertical wavelength to this level
              avlambdaz = 2.0 * pi / n_u_av(i)

              ! Calculate new NSQ if it has not already converged
              if (l_cont3(i,k)) then

                fracwave = dzt / avlambdaz
                nwave = int( fracwave )      ! Number of full waves below
                fracwave = fracwave - nwave  ! Fraction of current wave

                ! Maximum ascent and descent at this level
                if (fracwave < 0.5) then   !  Lower half of the wave
                  dzw = dzt - (nwave * avlambdaz)
                  maxdown =  wave_amp(i,k) * dzw*2.0/avlambdaz
                  maxup   =  wave_amp(i,k) * (1.0 - dzw*2.0/avlambdaz)
                else   !  Upper half of the wave
                  dzw = dzt - 0.5 * avlambdaz - (nwave * avlambdaz)
                  maxup   = wave_amp(i,k) * dzw*2.0/avlambdaz
                  maxdown  = wave_amp(i,k) * (1.0 - dzw*2.0/avlambdaz)
                end if

                ! Saturated and unsaturated regions of vertical motion
                if ( dzcond(i,k) >= 0.0 ) then
                  satdz = max(0.0, maxup - dzcond(i,k))
                  drydz = maxdown + min( dzcond(i,k), maxup )
                else
                  satdz = maxup + min( abs(dzcond(i,k)), maxdown )
                  drydz = max(0.0, maxdown - abs(dzcond(i,k)) )
                end if

                ! New NSQ for current iteration
                nsq_new(i,k) = ( (nsq_unsat(i,k)*drydz) +                      &
                      (nsq_sat(i,k)*satdz) ) / ( wave_amp(i,k) )

                ! Test for convergence of the buoyancy frequency
                if ((nsq(i,k) < nsq_new(i,k)*(1.0+nsq_converge)) .and.         &
                    (nsq(i,k) > nsq_new(i,k)*(1.0-nsq_converge))) then
                  l_cont3(i,k)  =  .false.
                end if ! test if NSQ(i,k) is converged

                ! Put new NSQ value into actual NSQ array
                nsq(i,k)  = nsq_new(i,k)

              end if  !  If l_cont3 (NSQ not converged)

            end if ! Wave amplitude and l_lapse not zero

          end if  !  (k >= ktop(i) .and. theta_levels(i,k) <=maxgwdlev)

        end if ! Sub-grid mountain present  l_drag
      end do ! i = 1, points
    end do ! k

  end do   ! iterate GW vertical propagation calculations until NSQ converges

else if (i_moist == 1) then

  ! Used Nmoist elsewhere so now use dry for propagation
  ! NB if i_moist=0 this is already set
  do k = 1, levels
    do i = 1, points
      if (l_drag(i)) then
        nsq(i,k) = nsq_dry(i,k)
      end if
    end do
  end do

end if   ! i_moist


! Actual wave propagation for final NSQ value (moist or dry)
! ----------------------------------------------------------
do k=1, levels-1
  do i = 1, points
    if (l_drag(i)) then
      if (nsq(i,k) > nsq_neutral) then
        ! first consider stable cases
        n(i,k)   = sqrt(nsq(i,k))
        u_n(i,k) = spd(i,k)/n(i,k)
      else
        ! now consider neutral (and near neutral) cases
        u_n(i,k) = spd(i,k)/sqrt(nsq_neutral)
      end if
      ! limit u_n to sensible values
      u_n(i,k)     = max(lambdaz_min,u_n(i,k))
      u_n(i,k)     = min(lambdaz_max,u_n(i,k))
      lambdaz(i,k) = frac_wl*(2.0*pi*u_n(i,k))
      deltaz(i,k)  = rho_levels(i,k+1)-rho_levels(i,k)
      ! set beta as beta_fix
      beta(i,k)    = beta_fix
    end if
  end do
end do

if (l_dynbeta) then
  do k=1, levels-1
    do i = 1, points
      if (l_drag(i)) then
        !Calculate dynamically adjusting group velocity angle at each height
        if (spd(i,k) /= 0.0) then
          if ((k_wave(i)**2*u_n(i,k)**2) < 1.0) then
            beta(i,k) = sqrt(1.0-k_wave(i)**2*u_n(i,k)**2)                     &
                             /(k_wave(i)*u_n(i,k))
          else
            beta(i,k) = 1.0
          end if
        else
          beta(i,k) = 1.0e16
        end if
      end if
    end do
  end do
end if

!-----------------------------------------------------------------
! Calculate stress lost out of grid-box due to non-hydrostatic prop.
! and deposit over height that it is lost
!-----------------------------------------------------------------
if (l_nonhydro) then
  do k = 1, levels-1
    do i = 1, points
      if (l_drag(i)) then
        cpsi1_dzdy  = abs(cos(psi1(i)))*deltaz(i,k)*deltay(i,k)
        spsi1_dzdx  = abs(sin(psi1(i)))*deltaz(i,k)*deltax(i,k)
        nonhyd(i,k) = (cpsi1_dzdy + spsi1_dzdx)                                &
                      /(beta(i,k)*deltax(i,k)*deltay(i,k)                      &
                       + cpsi1_dzdy + spsi1_dzdx)
      end if
    end do
  end do
end if !l_nonhydro

!-----------------------------------
! Wave propagation following McFarlane
!-----------------------------------

! Start loop at mountain top
! Just sfc stress as launch stress from mtn top

do k=2, levels-1
  do i = 1, points
    if (l_drag(i)) then
      if (k >= ktop(i) .and. theta_levels(i,k) <=maxgwdlev) then
        acrit       = 0.0
        kl          = k-1
        if ((wave_amp(i,kl)  /=  0.0) .and. (tau(i,kl) /= 0.0)) then
          !   All wave stress deposited if dth/dz < 0 or wind more than
          !   pi/2 to direction of surface stress
          if ((nsq(i,k) <= nsq_neutral) .or.                                   &
              (spd(i,k) <= 0.0)) then
            tau(i,k) = 0.0
          else
            ! Test whether wave amplitude exceeds critical amplitude,
            ! acrit, for wave breaking
            dzl           = theta_levels(i,k)   -  rho_levels(i,k)
            dzu           = rho_levels(i,k+1)   -  theta_levels(i,k)
            dzb           = rho_levels(i,k+1)   -   rho_levels(i,k)
            rhob          = (dzu*rho(i,k)       + dzl*rho(i,k+1))/dzb

            dzl           = theta_levels(i,k-1) -  rho_levels(i,k-1)
            dzu           = rho_levels(i,k)     -  theta_levels(i,k-1)
            dzb           = rho_levels(i,k)     -   rho_levels(i,k-1)
            rhob_l        = (dzu*rho(i,k-1)     + dzl*rho(i,k))  /dzb
            amplow        = rhob_l*n(i,kl)*spd(i,kl)
            wave_amp(i,k) = wave_amp(i,kl)*                                    &
                            sqrt(amplow/(rhob*spd(i,k)*n(i,k)))
            acrit         = fsat*(spd(i,k)/n(i,k))
            if (wave_amp(i,k) > acrit) then
              wave_amp(i,k) = acrit
            end if
            tau(i,k) = tau(i,kl)*                                              &
              (wave_amp(i,k) / wave_amp(i,kl))**2*                             &
              (rhob*n(i,k)*spd(i,k)) / (rhob_l*n(i,kl)*spd(i,kl)) *            &
              (1.0-nonhyd(i,k))
          end if  !  n(i,k)<= 0 and spd(i,k)<= 0
        end if !((wave_amp(i,kl)  /=  0.0) .and. (tau(i,kl) /= 0.0))
      end if !  (k >= ktop(i) .and. theta_levels(i,k) <=maxgwdlev)
    end if !(l_drag(i))
  end do !  i = 1, points
end do !k=1, levels-1
! End of repeat of wave propagation for final moist NSQ (or using Ndry)


! For this if test should i have some epsilon so we only calc stress
! where neccessary and never diagnose accs due to rounding error?
if (l_smooth) then
  do k=1, levels-1
    do i = 1, points
      if (l_drag(i)) then
        if (k >= ktop(i) .and. theta_levels(i,k) <=maxgwdlev) then
          kl = k-1
          if ((wave_amp(i,kl)  /=  0.0) .and. (tau(i,kl) /= 0.0)) then
            if (tau(i,k) /= tau(i,kl)) then
              ! Find region to apply stress over (vertical wavelength)
              l_cont(i)  = .true.
              l_cont2(i) = .true.
              khigh(i,k) = levels-1
              klow(i,k)  = 2
              do kk = k, 1, -1
                if ((theta_levels(i,kk) <=                                     &
                    (theta_levels(i,k)-0.5*lambdaz(i,k))) .and.                &
                    (l_cont(i))) then
                  klow(i,k) =  kk
                  l_cont(i) = .false.
                end if
              end do !kk, level k down to surface
              do kk = k, levels
                if ((theta_levels(i,kk) >=                                     &
                    (theta_levels(i,k)+0.5*lambdaz(i,k))) .and.                &
                    (l_cont2(i))) then
                  khigh(i,k) =  kk
                  l_cont2(i) = .false.
                end if
              end do !kk, level k up to model top
              if (klow(i,k)  <=  ktop(i)) then
                klow(i,k)    = 1
                ! lambdaz is hydrostatic vertical wavelength on theta levels
                lambdaz(i,k) = theta_levels(i,khigh(i,k))
              else
                lambdaz(i,k) = theta_levels(i,khigh(i,k)) -                    &
                               theta_levels(i,klow(i,k))
              end if !(klow(i,k) <=ktop(i))
            end if ! (tau(i,k) /= tau(i,kl))
          end if !((wave_amp(i,kl)  /=  0.0) .and. (tau(i,kl) /= 0.0))
        end if !  (k >= ktop(i) .and. theta_levels(i,k) <=maxgwdlev)
      end if !(l_drag(i))
    end do !  i = 1, points
  end do !k=1, levels-1
else !if l_smooth is false
  do k=1, levels-1
    do i = 1, points
      if (l_drag(i)) then
        if (k >= ktop(i) .and. theta_levels(i,k) <=maxgwdlev) then
          kl          = k-1
          if ((wave_amp(i,kl)  /=  0.0) .and. (tau(i,kl) /= 0.0)) then
            if (tau(i,k) /= tau(i,kl)) then
              khigh(i,k)  = k
              klow(i,k)   = k
              lambdaz(i,k)= deltaz(i,k)
              if (k  ==  ktop(i)) then
                klow(i,k)    = 1
                lambdaz(i,k) = theta_levels(i,k)
              end if !(k  ==  ktop(i))
            end if !(tau(i,k) /= tau(i,kl))
          end if !((wave_amp(i,kl)  /=  0.0) .and. (tau(i,kl) /= 0.0))
        end if !  (k >= ktop(i) .and. theta_levels(i,k) <=maxgwdlev)
      end if !(l_drag(i))
    end do !  i = 1, points
  end do !k=1, levels-1
end if !(l_smooth)
!------------------------------------------------------------------
! Calculate drag from vertical stress divergence
! Stress in theta (theta_levels) levels so stress divergence
! and hence acceleration are on rho (z) levels
!------------------------------------------------------------------

! need to resolve tau into x and y direction
do k=1, levels-1
  do i = 1, points
    if (l_drag(i)) then
      if (k >= ktop(i) .and. theta_levels(i,k) <=maxgwdlev) then
        kl          = k-1
        if ((wave_amp(i,kl)  /=  0.0) .and. (tau(i,kl) /= 0.0)) then
          if (tau(i,k) /= tau(i,kl)) then
            ztemp=(tau(i,k) - tau(i,kl))  /  lambdaz(i,k)
            uacc =((ulow(i)*pd1(i)-vlow(i)*pd2(i))* ztemp/pdmod(i))/modu(i)
            vacc =((vlow(i)*pd1(i)+ulow(i)*pd2(i))* ztemp/pdmod(i))/modu(i)
            if (scale_aware) then
              uacc = taux_sfc(i)*ztemp/tau_sfc(i)
              vacc = tauy_sfc(i)*ztemp/tau_sfc(i)
            end if ! scale_aware
            rhoav = 0.0
            if (klow(i,k) == 1) then
              rhoav = rho(i,1) * theta_levels(i,1)
              do kk = 2, khigh(i,k)
                rhoav = rhoav + rho(i,kk)*                                     &
                       (theta_levels(i,kk)-theta_levels(i,kk-1))
              end do
            else
              do kk = klow(i,k), khigh(i,k)
                rhoav = rhoav + rho(i,kk)*                                     &
                       (theta_levels(i,kk)-theta_levels(i,kk-1))
              end do
            end if !klow(i,k) == 1
            rhoav = rhoav / lambdaz(i,k)
            do kk = klow(i,k), khigh(i,k)
              dudt_gw(i,kk) = uacc/rhoav  + dudt_gw(i,kk)
              dvdt_gw(i,kk) = vacc/rhoav  + dvdt_gw(i,kk)
            end do
          end if !tau(i,k) /= tau(i,kl)
        end if    ! wave_amp(kl)  /=  0.0
      end if !k >= ktop(i) and (theta_levels(i,k) <=maxgwdlev)
    end if !(l_drag(i)
  end do ! Loop over points
end do ! Loop over levels

!update arrays with total acceleration (gravity wave + flow blocking)
do k=1, levels-1
  do i = 1, points
    dudt(i,k) = dudt(i,k)  + dudt_gw(i,k)
    dvdt(i,k) = dvdt(i,k)  + dvdt_gw(i,k)
  end do ! Loop over points
end do ! Loop over levels

!-----------------------------------------------------------------
! Calculate heating due to gravity wave dissipation
!-----------------------------------------------------------------
if ( l_gw_heating ) then

  do k = 1, levels-1
    do i = 1, points

      dzb    = theta_levels(i,k)   -  rho_levels(i,k)
      dzu    = rho_levels(i,k+1)   -  theta_levels(i,k)
      dzl    = rho_levels(i,k+1)   -   rho_levels(i,k)

      ! u and v on theta_level(k)
      uhat   =  dzu * u(i,k) + dzb * u(i,k+1)
      vhat   =  dzu * v(i,k) + dzb * v(i,k+1)

      ! u*du/dt abd v*dv/dt on theta_level(k)
      ududt  = uhat *( dzu * dudt_gw(i,k) + dzb * dudt_gw(i,k+1) )
      vdvdt  = vhat *( dzu * dvdt_gw(i,k) + dzb * dvdt_gw(i,k+1) )

      ! dT/dt on theta_level(k)
      dtdt(i,k)  = dtdt(i,k) - (ududt + vdvdt) / ( cp * dzl * dzl )

    end do !Loop over points
  end do !Loop over levels

end if !l_gw_heating





!------------------------------------------------------------------
! diagnostics
!------------------------------------------------------------------

if ( du_dt_satn_on .or. du_dt_satn_p_on ) then
  do k=1, levels-1
    do i = 1, points
      if (l_drag(i)) then
        du_dt_satn(i,k) = dudt_gw(i,k)
      end if !(l_drag(i)
    end do ! Loop over points
  end do ! Loop over levels
end if

if ( dv_dt_satn_on ) then
  do k=1, levels-1
    do i = 1, points
      if (l_drag(i)) then
        dv_dt_satn(i,k) = dvdt_gw(i,k)
      end if !(l_drag(i)
    end do ! Loop over points
  end do ! Loop over levels
end if

if ( stress_ud_on .or. stress_ud_p_on .or. stress_ud_satn_on                   &
                 .or. tausx_d_on ) then
  do k =levels-1,0,-1
    do i=1,points
      if (l_drag(i)) then
        if ( k  ==  0 ) then
          dz = theta_levels(i,k+1)
        else
          dz = theta_levels(i,k+1) - theta_levels(i,k)
        end if
        local_xstress(i,k) = local_xstress(i,k+1) -                            &
                            (dudt_gw(i,k+1)*rho(i,k+1)*dz)
      end if
    end do
  end do
  if ( stress_ud_on .or. stress_ud_p_on) then
    do k =0, levels-1
      do i=1,points
        if (l_drag(i)) then
          stress_ud(i,k) = stress_ud(i,k) + local_xstress(i,k)
        end if
      end do
    end do
    do i=1,points
      stress_ud(i,levels) = stress_ud(i,levels-1)
    end do
  end if
  if ( stress_ud_satn_on ) then
    do k =0, levels-1
      do i=1,points
        if (l_drag(i)) then
          stress_ud_satn(i,k) = local_xstress(i,k)
        end if
      end do
    end do
    do i=1,points
      stress_ud_satn(i,levels) = stress_ud_satn(i,levels-1)
    end do
  end if
  if ( tausx_d_on ) then
    do i=1,points
      if (l_drag(i)) then
        tausx_d(i) = tausx_d(i) + local_xstress(i,0)
      end if
    end do
  end if
end if

if ( stress_vd_on .or. stress_vd_satn_on .or. tausy_d_on ) then
  do k =levels-1,0,-1
    do i=1,points
      if (l_drag(i)) then
        if ( k  ==  0 ) then
          dz = theta_levels(i,k+1)
        else
          dz = theta_levels(i,k+1) - theta_levels(i,k)
        end if
        local_ystress(i,k) = local_ystress(i,k+1) -                            &
                             (dvdt_gw(i,k+1)*rho(i,k+1)*dz)
      end if
    end do
  end do
  if ( stress_vd_on ) then
    do k =0, levels-1
      do i=1,points
        if (l_drag(i)) then
          stress_vd(i,k) =  stress_vd(i,k) + local_ystress(i,k)
        end if
      end do
    end do
    do i=1,points
      stress_vd(i,levels) = stress_vd(i,levels-1)
    end do
  end if
  if ( stress_vd_satn_on ) then
    do k =0, levels-1
      do i=1,points
        if (l_drag(i)) then
          stress_vd_satn(i,k) = local_ystress(i,k)
        end if
      end do
    end do
    do i=1,points
      stress_vd_satn(i,levels) = stress_vd_satn(i,levels-1)
    end do
  end if
  if ( tausy_d_on ) then
    do i=1,points
      if (l_drag(i)) then
        tausy_d(i) = tausy_d(i) + local_ystress(i,0)
      end if
    end do
  end if
end if

if ( nsq_s_d_on ) then
  do i=1,points
    nsq_s_d(i) = nlow(i)
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine gw_wave
end module gw_wave_mod
