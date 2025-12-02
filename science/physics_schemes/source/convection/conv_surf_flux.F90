! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module conv_surf_flux_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='CONV_SURF_FLUX_MOD'

contains

! Calculate surface fluxes for convective diagnosis

subroutine conv_surf_flux(                                                     &
  row_length, rows                                                             &
, l_spec_z0 , land_mask                                                        &
, pstar, tstar_land, tstar_sea, tstar_sice, zh, flandg                         &
, ice_fract, u_p, v_p, u_0_p, v_0_p                                            &
, flux_e, flux_h,  ustar_in, z0msea, z0m_scm, z0h_scm                          &
, z_full, q, theta, exner_theta_levels, ccp_strength                           &
! INOUT
, tstar , fb_surf, tv1_sd, bl_vscale2                                          &
 )

!-----------------------------------------------------------------------
! Purpose:
!  Calculate the surface buoyancy flux.
!  Also calculates an approximate standard deviation of virtual
!  temperature at level 1
!
!  Part of convective diagnosis routines
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!  Language: FORTRAN90
!  This code is written to UMDP3  programming standards v8.3
!
!-----------------------------------------------------------------------

! Definitions of prognostic variable array sizes
use atm_fields_bounds_mod, only:                                               &
  pdims_s, tdims

use cv_run_mod, only:                                                          &
  tv1_sd_opt, cnv_cold_pools, l_ccp_trig

use cv_param_mod, only:                                                        &
  ccp_off, ccp_fb_coeff, ccp_fb_offset

use planet_constants_mod, only:                                                &
  vkman, cp, kappa, r, repsilon, g

use water_constants_mod, only: lc

use cv_derived_constants_mod, only:   gamma_dry

use bl_option_mod, only: flux_bc_opt, interactive_fluxes,                      &
                         specified_fluxes_only, specified_fluxes_cd

use jules_sea_seaice_mod, only: z0sice, z0h_z0m_sice, z0hsea

use qsat_mod, only: qsat, qsat_mix

use gen_phys_inputs_mod, only: l_mr_physics

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

!-----------------------------------------------------------------------
! Subroutine Arguments
!-----------------------------------------------------------------------

integer, intent(in) ::                                                         &
  row_length            & ! Local number of points on a row
 ,rows                    ! Local number of rows in a theta field

logical,intent(in) ::                                                          &
  l_spec_z0               ! true if roughness length has been specIFied

logical,intent(in) ::                                                          &
  land_mask(row_length, rows)    ! T if land, F elsewhere.

real(kind=real_umphys), intent(in) ::                                          &
  pstar(row_length, rows)      & ! Surface pressure (Pa)
 ,tstar_land(row_length, rows) & ! Surface T on land
 ,tstar_sea(row_length, rows)  & ! Surface T on sea
 ,tstar_sice(row_length, rows) & ! Surface T on sea-ice
 ,zh(row_length,rows)          & ! Height above surface of top
                                 !  of boundary layer (metres).
 ,ice_fract(row_length,rows)     ! fraction of sea that has ice

real(kind=real_umphys), intent(in) ::                                          &
  flandg(pdims_s%i_start:pdims_s%i_end,         & ! Land fraction of gridbox
         pdims_s%j_start:pdims_s%j_end)         & ! on all points
 ,z_full(row_length,rows,1:tdims%k_end)         & ! height th lev (m)
 ,q(tdims%i_start:tdims%i_end,                  & ! water vapour (kg/kg)
    tdims%j_start:tdims%j_end,                                                 &
                1:tdims%k_end)                                                 &
 ,theta(tdims%i_start:tdims%i_end,              & ! Theta (Kelvin)
        tdims%j_start:tdims%j_end,                                             &
                    1:tdims%k_end)                                             &
 ,exner_theta_levels(tdims%i_start:tdims%i_end, & ! exner pressure theta lev
                     tdims%j_start:tdims%j_end, & !  (Pa)
                                 1:tdims%k_end)

real(kind=real_umphys), intent(in) ::                                          &
  u_p(row_length, rows)        & ! U(1) on P-grid.
 ,v_p(row_length, rows)        & ! V(1) on P-grid.
 ,u_0_p(row_length,rows)       & ! W'ly component of surface current
                                 !    (metres per second) on P-grid.
 ,v_0_p(row_length,rows)       & ! S'ly component of surface current
                                 !    (metres per second) on P-grid.
 ,flux_e(row_length,rows)      & ! Specified surface
                                 !    latent heat flux (W/m^2)
 ,flux_h(row_length,rows)      & ! Specified surface
                                 !    sensible heat fluxes (in W/m2)
 ,ustar_in(row_length,rows)    & ! Specified surface friction velocity (m/s)
 ,z0msea(row_length,rows)      & ! Sea roughness length for momentum (m)
 ,z0m_scm(row_length,rows)     & ! Namelist input z0m (if >0)
 ,z0h_scm(row_length,rows)       ! Namelist input z0h (if >0)

real(kind=real_umphys), intent(in) ::                                          &
  ccp_strength(row_length,rows) ! measure of cold-pool strength (dimensionless)

real(kind=real_umphys), intent(in out) ::                                      &
  tstar(row_length,rows)        & ! Surface temperature
                                  ! (= top soil layer temperature) (K).
 ,fb_surf(row_length, rows)     & ! Change in theta_v from surface
                                  ! to layer 1 (note diff from BL)
 ,tv1_sd( row_length*rows)      & ! Approx to standard dev of level 1
                                  ! virtual temperature (K).(unstable points)
 ,bl_vscale2(row_length*rows)     ! Velocity scale squared for
                                  ! boundary layer eddies (m2/s2)

!-----------------------------------------------------------------------
! Variables defined locally
!-----------------------------------------------------------------------

character(len=*), parameter ::  RoutineName = 'CONV_SURF_FLUX'

integer ::                                                                     &
  i,j,ii       ! Local Loop counter (horizontal field index).

real(kind=real_umphys) ::                                                      &
  qs_star(row_length, rows)     & ! Saturated sp humidity at surface
 ,qs_star_sice(row_length, rows)& ! Saturated sp humidity at sea-ice surface
 ,dqsdt(row_length,rows)        & ! d(qsat)/dT
 ,z0(row_length,rows)           & ! roughness length (m)
 ,z0m_land(row_length,rows)     & ! roughness length for momentum over land (m)
 ,z0h_land(row_length,rows)     & ! roughness length for heat over land (m)
 ,z0m_sea(row_length,rows)      & ! roughness length for momentum over sea(m)
 ,z0h_sea(row_length,rows)        ! roughness length for heat over sea(m)


! Used in calculation to decide on unstable points

real(kind=real_umphys) ::                                                      &
  theta1          &  ! Potential temperature in layer 1
 ,ushear          &  ! U wind shear from level 1 to surface
 ,vshear          &  ! V wind shear from level 1 to surface
 ,wshr1           &  ! magnitude of surface wind shear
 ,wshr2           &  ! (wshr1)**2
 ,rhostar         &  ! surface air density
 ,theta_star      &  ! theta at surface
 ,wthvbar         &  ! surface buoyancy flux
 ,cd              &  ! bulk transfer coefficient for momentum
 ,ch              &  ! bulk transfer coefficient for heat
 ,rib             &  ! Ri for surface exchange
 ,ustar           &  ! surface friction velocity
 ,w_m             &  !
 ,w_s_cubed       &  !
 ,ustar2_sea      &  ! ustar^2 over sea
 ,ustar2_land     &  ! ustar^2 over land
 ,wthvbar_sea     &  ! surface buoyancy flux over sea
 ,wthvbar_land       ! surface buoyancy flux

real(kind=real_umphys)  :: tv1_sd_temp(row_length, rows),                      &
         bl_vscale2_temp(row_length, rows)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Set up roughness lengths
!-----------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(SHARED) private(i, j)

if (tv1_sd_opt == 2) then

  ! using grid-box mean in coastal points, hence need land and sea
  ! separately
!$OMP do SCHEDULE(STATIC)
  do j=1, rows
    do i=1, row_length
      ! Land: assume z0m=0.1m and z0h=z0m/10.
      z0m_land(i,j) = 0.1
      z0h_land(i,j) = 0.01
      ! Sea: use parametrized values
      !      z0h_sea updated later for low wind speed limit
      z0m_sea(i,j) = z0msea(i,j)
      z0h_sea(i,j) = max( 2.56e-9/z0msea(i,j), 7.0e-08 )
    end do ! i
  end do ! j
!$OMP end do

  if ( L_spec_z0 ) then
    ! Code to use z0mh_scm if namelist specifies it
!$OMP do SCHEDULE(STATIC)
    do j=1, rows
      do i=1, row_length
        if ( z0m_scm(i,j)  >   0.0 ) then
          z0m_sea(i,j)  = z0m_scm(i,j)
          z0m_land(i,j) = z0m_scm(i,j)
        end if
        if ( z0h_SCM(i,j)  >   0.0 ) then
          z0h_sea(i,j)  = z0h_scm(i,j)
          z0h_land(i,j) = z0h_scm(i,j)
        end if
      end do ! i
    end do ! j
!$OMP end do
  end if
end if

!$OMP do SCHEDULE(STATIC)
do j=1, rows
  do i=1, row_length
    if (land_mask(i,j)) then
      ! Approximate z0 for land as 0.1.
      z0(i,j) = 0.1
    else
      z0(i,j) = z0hsea
    end if
  end do ! i
end do ! j
!$OMP end do

! Code to use z0h_scm if Namelist specIFies it

if ( L_spec_z0 ) then
!$OMP do SCHEDULE(STATIC)
  do j=1, rows
    do i=1, row_length
      if ( z0h_SCM(i,j)  >   0.0 ) then
        z0(i,j) = z0h_SCM(i,j)
      end if! z0h_scm
    end do ! i
  end do ! j
!$OMP end do
end if

!$OMP end PARALLEL

!-----------------------------------------------------------------------
! 1.5 Surface buoyancy flux and calculation of unstable points
!-----------------------------------------------------------------------

if ( flux_bc_opt == interactive_fluxes ) then    ! used by most UM runs

  if (tv1_sd_opt == 2) then
    !-----------------------------------------------------------------------
    ! Calculate the surface buoyancy flux
    ! new method includes stability dependence and area mean for coastal
    ! and sea-ice points
    !-----------------------------------------------------------------------

    if ( l_mr_physics ) then
      call qsat_mix(qs_star,tstar_sea,pstar,row_length,rows)
      call qsat_mix(qs_star_sice,tstar_sice,pstar,row_length,rows)
    else
      call qsat(qs_star,tstar_sea,pstar,row_length,rows)
      call qsat(qs_star_sice,tstar_sice,pstar,row_length,rows)
    end if

!$OMP  PARALLEL DEFAULT(SHARED) private(theta1, rhostar, Ushear, wshr1,        &
!$OMP  Vshear, wthvbar_sea, rib, cd, ustar2_sea, ch, w_m, wshr2,               &
!$OMP  wthvbar_land, ustar2_land, w_s_cubed, i, j, ustar, wthvbar)

!$OMP do SCHEDULE(STATIC)
    do j=1,rows
      do i=1,row_length

        theta1 = theta(i,j,1)
        rhostar = pstar(i,j) / ( r*tstar(i,j) )

        ushear = u_p(i,j) - u_0_p(i,j)
        vshear = v_p(i,j) - v_0_p(i,j)
        wshr2 = max (1.0e-6 , ushear*ushear + vshear*vshear)
        wshr1 = sqrt(wshr2)
        !-------------------------------------------------------------
        ! Sea
        !-------------------------------------------------------------
        wthvbar_sea  = 0.0
        ustar2_sea   = 0.0
        if ( flandg(i,j) < 0.99 ) then
          ! Include a crude stability dependence
          rib = - ( g / tstar_sea(i,j) ) *                                     &
            (tstar_sea(i,j) - ( theta1*exner_theta_levels(i,j,1)               &
                    +gamma_dry*z_full(i,j,1) ) ) * z_full(i,j,1) / wshr2
          cd = ( vkman / log(z_full(i,j,1)/z0m_sea(i,j)) )**2
          cd = cd * (1.0 + 0.7 * (max(0.0, -rib))**0.33 )
          ustar2_sea = cd * wshr2
          if ( .not. l_spec_z0 ) then
            ! include low wind speed limit (larger z0)
            z0h_sea(i,j) = max( 2.52e-6/(sqrt(ustar2_sea)+1.0e-05),            &
                                 z0h_sea(i,j) )
          end if
          ch = vkman**2 / ( log(z_full(i,j,1)/z0m_sea(i,j)) *                  &
                             log(z_full(i,j,1)/z0h_sea(i,j) ) )
          ch = ch * (1.0 + (max(0.0, -rib))**0.33 )
          wthvbar_sea = ch * wshr1 * ( tstar_sea(i,j) -                        &
             ( theta1*exner_theta_levels(i,j,1) + gamma_dry*z_full(i,j,1) )    &
               + 0.61*theta1*(qs_star(i,j)-q(i,j,1)) )
          !-------------------------------------------------------------
          ! Sea-ice
          !-------------------------------------------------------------
          if ( ice_fract(i,j) > 0.01 ) then
            ! Include a crude stability dependence
            rib = - ( g / tstar_sice(i,j) ) *                                  &
              (tstar_sice(i,j) - ( theta1 * exner_theta_levels(i,j,1)          &
                    +gamma_dry*z_full(i,j,1) ) ) * z_full(i,j,1) / wshr2
            cd = ( vkman / log(z_full(i,j,1)/z0sice) )**2
            cd = cd * (1.0 + 0.7 * (max(0.0, -rib))**0.33 )
            ustar2_sea = (1.0-ice_fract(i,j)) * ustar2_sea +                   &
                               ice_fract(i,j) * cd*wshr2
            ch = vkman**2 / ( log(z_full(i,j,1)/z0sice) *                      &
                              log(z_full(i,j,1)/(z0sice*z0h_z0m_sice)) )
            ch = ch * (1.0 + (max(0.0, -rib))**0.33 )
            wthvbar_sea = (1.0-ice_fract(i,j)) * wthvbar_sea +                 &
                        ice_fract(i,j) * ch*wshr1*( tstar_sice(i,j) -          &
              ( theta1*exner_theta_levels(i,j,1) + gamma_dry*z_full(i,j,1) )   &
                + 0.61*theta1*(qs_star_sice(i,j)-q(i,j,1)) )
          end if
        end if

        !-------------------------------------------------------------
        ! Land
        !-------------------------------------------------------------
        wthvbar_land = 0.0
        ustar2_land   = 0.0
        if ( flandg(i,j) > 0.01 ) then
          ! Include a crude stability dependence
          rib = - ( g / tstar_land(i,j) ) *                                    &
            (tstar_land(i,j) - ( theta1 * exner_theta_levels(i,j,1)            &
             +gamma_dry*z_full(i,j,1) ) ) * z_full(i,j,1) / wshr2
          cd = ( vkman / log(z_full(i,j,1)/z0m_land(i,j)) ) ** 2
          cd = cd * (1.0 + 0.7 * (max(0.0, -rib))**0.33 )
          ustar2_land = cd * wshr2
          ch = vkman**2 / ( log(z_full(i,j,1)/z0m_land(i,j)) *                 &
                             log(z_full(i,j,1)/z0h_land(i,j) ) )
          ch = ch * (1.0 + (max(0.0, -rib))**0.33 )
          wthvbar_land = ch * wshr1 *                                          &
            ( tstar_land(i,j) - ( theta1 * exner_theta_levels(i,j,1)           &
             +gamma_dry*z_full(i,j,1) ) )
        end if
        !-------------------------------------------------------------
        ! Combine to cell average values and then take sqrt for ustar
        !-------------------------------------------------------------
        if ( flandg(i,j) < 0.01 ) then
          ustar = ustar2_sea
          wthvbar = wthvbar_sea
        else if ( flandg(i,j) > 0.99 ) then
          ustar = ustar2_land
          wthvbar = wthvbar_land
        else
          ! Take area-weighted mean
          ustar = (1.0 - flandg(i,j) ) * ustar2_sea +                          &
                     flandg(i,j) * ustar2_land
          wthvbar = (1.0 - flandg(i,j) ) * wthvbar_sea +                       &
                     flandg(i,j) * wthvbar_land
        end if

        ustar = sqrt(ustar)
        fb_surf(i,j) = g * wthvbar /( rhostar * theta1*(1.0+0.61*q(i,j,1)) )
        if ((cnv_cold_pools > ccp_off) .and. l_ccp_trig) then
          fb_surf(i,j) = fb_surf(i,j)                                          &
                        + ccp_fb_coeff * ( ccp_strength(i,j) - ccp_fb_offset )
        end if

        if (fb_surf(i,j)  >   0.0) then
          w_s_cubed  = 0.25 * zh(i,j) * fb_surf(i,j)
          w_m        = (w_s_cubed + ustar*ustar*ustar)**(1.0/3.0)
          tv1_sd_temp(i, j) = 1.93 * wthvbar/( rhostar * w_m )
          bl_vscale2_temp(i,j) = 2.5 * 2.52 * w_m * w_m
                  ! 2.5 from Beare (2008), 2.52 to undo the 0.25 factor
          bl_vscale2_temp(i,j) = max( 0.01, bl_vscale2_temp(i,j) ) ! for safety
        end if! (fb_surf > 0.0)

      end do ! i
    end do ! j
!$OMP end do

!$OMP end PARALLEL

    ii = 0

    do j=1, rows
      do i=1, row_length
        if (fb_surf(i,j)  >   0.0) then
          ii= ii+1
          tv1_sd(ii) = tv1_sd_temp(i,j)
          bl_vscale2(ii) = bl_vscale2_temp(i,j)
        end if
      end do
    end do

  else  ! tv1_sd_opt /= 2
    ! old (neutral stability) method

    ! qsat at surface
    if ( l_mr_physics ) then
      call qsat_mix(qs_star,tstar,pstar,row_length,rows)
    else
      call qsat(qs_star,tstar,pstar,row_length,rows)
    end if

    !----------------------------------------------------------------
    ! Standard UM code for surface T boundary condition:
    ! Calculate the surface buoyancy flux using
    ! approximation for unstable cd as 1.5*neutral value (defined
    ! Garratt p54) and that ch=cd.
    !----------------------------------------------------------------
    ii=0
    do j=1,rows
      do i=1,row_length

        theta1 = theta(i,j,1)
        theta_star = tstar(i,j)*((100000.0/pstar(i,j))**kappa)
        rhostar = pstar(i,j) / ( r*tstar(i,j) )

        ushear = u_p(i,j) - u_0_p(i,j)
        vshear = v_p(i,j) - v_0_p(i,j)
        wshr2 = max (1.0e-6 , ushear*ushear + vshear*vshear)
        wshr1 = sqrt(wshr2)
        cd = 1.5 * ( vkman/log(z_full(i,j,1)/z0(i,j)) )**2

        if (land_mask(i,j)) then         ! land
          wthvbar = wshr1 * cd * ( theta_star - theta1 )
        else                             ! sea
          wthvbar = wshr1 * cd *                                               &
                  ( theta_star - theta1 + 0.61*theta1*(qs_star(i,j)-q(i,j,1)) )
        end if

        ustar = sqrt(cd * wshr2)
        fb_surf(i,j) = g * wthvbar / ( rhostar * theta1*(1.0+0.61*q(i,j,1)) )

        if (fb_surf(i,j)  >   0.0) then
          ii= ii+1
          ! improved method uses BL depth from the previous timestep
          w_s_cubed  = 0.25 * zh(i,j) * fb_surf(i,j)
          w_m        = (w_s_cubed + ustar*ustar*ustar)**(1.0/3.0)
          tv1_sd(ii) = 1.93 * wthvbar/( rhostar * w_m )
          bl_vscale2(ii) = 2.5 * 2.52 * w_m * w_m
                 ! 2.5 from Beare (2008), 2.52 to undo the 0.25 factor
          bl_vscale2(ii) = max( 0.01, bl_vscale2(ii) ) ! for safety
        end if

      end do ! i
    end do ! j
  end if ! test on tv1_sd_opt /= 2

else ! if flux_bc_opt/=interactive      (used by some SCM runs)

  !-----------------------------------------------------------
  ! Code for specified surface flux boundary condition.
  ! Assumes a saturated surface (ie. only appropriate for sea)
  !-----------------------------------------------------------

  ii=0
  if (flux_bc_opt == specified_fluxes_only) then
    ! tstar not specified so need to estimate from fluxes
    do j=1, rows
      do i=1, row_length
        ! For taylor expansion about T0=SL(K=1)
        tstar(i,j) = theta(i,j,1)*exner_theta_levels(i,j,1)+                   &
                                 gamma_dry*z_full(i,j,1)
      end do
    end do
  end if

  if ( l_mr_physics ) then
    call qsat_mix(qs_star,tstar,pstar,row_length,rows)
  else
    call qsat(qs_star,tstar,pstar,row_length,rows)
  end if

  do j=1, rows
    do i=1, row_length
      dqsdt(i,j) = (repsilon * lc * qs_star(i,j))                              &
                      / ( r * tstar(i,j) * tstar(i,j) )
    end do
  end do

  ii=0
  do j=1,rows
    do i=1,row_length

      ushear = u_p(i,j) - u_0_p(i,j)
      vshear = v_p(i,j) - v_0_p(i,j)
      ! Need to have a higher minimum wind speed limit with
      ! specified fluxes in order not to generate huge tstar
      wshr2 = max (0.1, ushear*ushear + vshear*vshear)
      wshr1 = sqrt(wshr2)

      ! Calculate wthv from namelist flux_h and flux_e (in W/m2)
      wthvbar = ((flux_h(i,j)/cp)+0.61*(flux_e(i,j)/lc))                       &
               * ((100000.0/pstar(i,j))**kappa)

      cd = 1.5 * ( vkman/log(z_full(i,j,1)/z0(i,j)) )**2

      theta1 = theta(i,j,1)
      ! Taylor expansion for qsat(T*) about SL(k=1)

      if (flux_bc_opt == specified_fluxes_only) then
        ! tstar not specified so need to estimate from fluxes
        tstar(i,j) = ( theta1 + (wthvbar/(wshr1*cd))                           &
               -   0.61*theta1                                                 &
               *   (qs_star(i,j)-q(i,j,1)-dqsdt(i,j)*tstar(i,j)))              &
               /   ( (100000.0/pstar(i,j))**kappa +                            &
               0.61*theta1*dqsdt(i,j) )
      end if
      rhostar = pstar(i,j) / ( r*tstar(i,j) )

      if (flux_bc_opt == specified_fluxes_cd) then
        ustar = ustar_in(i,j)
      else
        ustar = sqrt(cd * wshr2)
      end if
      fb_surf(i,j) = g * wthvbar /( rhostar * theta1*(1.0+0.61*q(i,j,1)) )

      if (fb_surf(i,j)  >   0.0) then
        ii= ii+1
        if (tv1_sd_opt == 0) then
          ! old method assumes the BL depth zh=300.
          tv1_sd(ii) = 1.93 * wthvbar / ( rhostar * ( 75.0 *                   &
                            fb_surf(i,j) + ustar*ustar*ustar)**(1.0/3.0) )
          w_s_cubed  = 0.25 * zh(i,j) * fb_surf(i,j)
          w_m        = (w_s_cubed + ustar*ustar*ustar)**(1.0/3.0)
          bl_vscale2(ii) = 2.5 * 2.52 * w_m * w_m
                  ! 2.5 from Beare (2008), 2.52 to undo the 0.25 factor
        else
          ! improved method uses BL depth from the previous timestep
          w_s_cubed  = 0.25 * zh(i,j) * fb_surf(i,j)
          w_m        = (w_s_cubed + ustar*ustar*ustar)**(1.0/3.0)
          tv1_sd(ii) = 1.93 * wthvbar/( rhostar * w_m )
          bl_vscale2(ii) = 2.5 * 2.52 * w_m * w_m
                  ! 2.5 from Beare (2008), 2.52 to undo the 0.25 factor
        end if
        bl_vscale2(ii) = max( 0.01, bl_vscale2(ii) ) ! for safety
      end if! (fb_surf > 0.0)

    end do ! i
  end do ! j

end if!  flux_bc_opt

!----------------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine conv_surf_flux

end module conv_surf_flux_mod
