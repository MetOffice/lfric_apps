! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! Purpose: Calculate form drag profiles for distributed
!           drag parametrization

!           Based on work by Wood, Brown and Hewer (2001),
!           Quart. J. Roy. Met. Soc., 127, 759--777.

! Programming standard: UMDP 3

! Description:
!  Turbulent form drag due to sub-grid orography

! Method:
!  Calculates the form drag and stress profiles due to sub-grid scale
!  orography. The stress is later added as an additional explicit
!  stress to the boundary-layer equations in ex_flux_uv. The orographic
!  stress is added to the surface stress in bdy_expl2.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer

! Code Description:
!  Language: FORTRAN 95
!  This code is written to UMDP3 programming standards.
!---------------------------------------------------------------------
module fm_drag_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'FM_DRAG_MOD'
contains
subroutine fm_drag (                                                           &
! in levels
  land_pts, land_index, bl_levels,                                             &
! in fields
  u_p, v_p, tl, qw, bt_gb, bq_gb, rho_wet_tq,                                  &
  z_uv, z_tq, z0m, zh, rib_surf, sil_orog_land,                                &
! out fields
  tau_fd_x, tau_fd_y                                                           &
  )

use atm_fields_bounds_mod, only: pdims, tdims
use bl_option_mod, only: on, zero, one, one_half
use c_surf, only: ri_crit
use conversions_mod, only: pi => pi_bl
use jules_surface_mod, only: fd_stability_dep, orog_drag_param, use_bulk_ri,   &
                             fd_hill_option, steep_hill, low_hill,             &
                             capped_lowhill
use planet_constants_mod, only: vkman => vkman_bl, grcp => grcp_bl, g => g_bl
use stochastic_physics_run_mod, only: l_rp2, orog_drag_param_rp, rp_idx,       &
                                      i_rp_scheme, i_rp2b
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

integer, intent(in):: land_pts,                                                &
                                           ! Number of land points
                      bl_levels        ! No. of levels for which
                                           ! boundary layer fluxes
                                           ! are calculated

integer, intent(in):: land_index(land_pts) ! Index for compressed
                                               ! land point array; ith
                                               ! element holds
                                               ! position in the FULL
                                               ! field of the ith land
                                               ! point to be processed

real(kind=r_bl), intent(in)::                                                  &
 u_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),           &
                                    ! Wind component in x direction
                                    ! horizontally interpolated to
                                    ! P-grid (m/s)
 v_p(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),           &
                                    ! Wind component in y direction
                                    ! horizontally interpolated to
                                    ! P-grid (m/s)
 tl(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),           &
                                    ! Ice/liquid water temperature
 qw(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end, bl_levels),           &
                                    ! Total water content
 bt_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                                    ! grid-box mean buoyancy param for
                                    ! tl on th-levels
 bq_gb(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),         &
                                    ! grid-box mean buoyancy param for
                                    ! qw on th-levels
 rho_wet_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,               &
            bl_levels),                                                        &
                                    ! For a vertically staggered grid
                                    ! with a u,v-level first above the
                                    ! surface, RHO_WET_TQ(*,K) is the
                                    ! density of the k-th T,q-level
                                    ! above the surface;
                                    ! for an unstaggered grid the
                                    ! densities at the layer interface
                                    ! (half-levels) 1.5 to BL_LEVELS+0
                                    ! should be input to elements 1 to
                                    ! BL_LEVELS.
                                    ! (Value for BL_LEVELS not used
                                    ! in either case.)
 z_uv(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels+1),        &
                                    ! For a vertically staggered grid
                                    ! with a u,v-level first above the
                                    ! surface, Z_UV(*,K) is the height
                                    ! of the k-th u,v-level(half level
                                    ! k-1/2) above the surface;
                                    ! for an unstaggered grid the
                                    ! heights of the half-levels
                                    ! 0.5 to BL_LEVELS-0.5 should be
                                    ! input to elements 1 to BL_LEVELS
                                    ! (1st value not used in either
                                    !  case.)
 z_tq(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,bl_levels),          &
                                    ! For a vertically staggered grid
                                    ! with a u,v-level first above the
                                    ! surface, Z_TQ(*,K) is the height
                                    ! of the k-th T,q-level (full
                                    ! level k) above the surface;
                                    ! for an unstaggered grid the
                                    ! heights of the half levels
                                    ! 1.5 to BL_LEVELS+0.5 should be
                                    ! input to elements 1 to BL_LEVELS
                                    ! (Value for BL_LEVELS not used
                                    ! in either case.)
 z0m(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                     &
                                    ! Roughness length for momentum (m
 zh(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                      &
                                    ! Boundary layer height (actually
                                    ! ZH_PREV)
 rib_surf(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end),                &
                                    ! Bulk Richardson number for
                                    ! lowest layer
 sil_orog_land(land_pts)
                                    ! Silhouette area of unresolved
                                    ! orography per unit hoz. area

real(kind=r_bl), intent(out)::                                                 &
 tau_fd_x(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          bl_levels),                                                          &
                                         ! X-comp of orographic stress
 tau_fd_y(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end,                 &
          bl_levels) ! Y-comp of orographic stress
                                         ! (N/m^2)

!     Local and other symbolic constants:

!     Define local storage

!     Local arrays:

real(kind=r_bl)::                                                              &
  h_m(land_pts),                                                               &
                                   ! height at which form drag is
                                   ! calculated (m)
  u_hm(land_pts),                                                              &
                                   ! X-component of wind at height h_m
  v_hm(land_pts),                                                              &
                                   ! Y-component of wind at height h_m
  tl_hm(land_pts),                                                             &
                                   ! TL at height h_m
  qw_hm(land_pts),                                                             &
                                   ! QW at height h_m
  rib(land_pts),                                                               &
                                   ! Bulk Richardson number for surface
                                   ! to scale height, h_m
  db_surf(land_pts),                                                           &
                                   ! Buoyancy difference between the
                                   ! surface and level 1
  fp_x(land_pts),                                                              &
                                   ! X-component of pressure force
  fp_y(land_pts)                   ! Y-component of pressure force

integer::                                                                      &
  k_for_buoy(land_pts)             ! index marking level below h_m,
                                   ! used for buoyancy coefficients

integer::                                                                      &
  land_index_i(land_pts), land_index_j(land_pts)
                                   ! arrays to store horizontal grid indices
                                   ! for land points

!     Local scalars:

real(kind=r_bl)::                                                              &
  zeta,                                                                        &
                    ! Log(h_m/z0m)
  height_fac,                                                                  &
                    ! Height dependency factor for form drag
  wta, wtb,                                                                    &
                    ! weights for interpolation between levels
  tausx,tausy,                                                                 &
                    ! Surface stress
  fp_x_low,fp_y_low,fp_x_steep,fp_y_steep,                                     &
  rib_fn            ! Richardson number function for stability
                    ! correction to drag

integer::                                                                      &
  i,                                                                           &
                    ! Loop counter (horizontal field index)
  j,                                                                           &
                    ! Loop counter (offset within I loop)
  k,                                                                           &
                    ! Loop counter (vertical level index)
  l             ! Loop counter (horizontal land field index)

! Local parameters
      ! Tunable parameters in calculation of explicit orographic stresss
real(kind=r_bl),parameter:: alpha    = 12.0_r_bl,                              &
                                           ! Tunable parameter for form
                 beta     = one,                                               &
                                           ! drag calculation.
                 fd_decay = 0.3333_r_bl,                                       &
                                           ! Height scale factors for
                 max_ht_scale  = 300.0_r_bl,                                   &
                                           ! stress profiles
                 min_ht_scale  = 100.0_r_bl,                                   &
                 pi_squared = pi * pi

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='FM_DRAG'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if ( l_rp2 .and. i_rp_scheme == i_rp2b ) then
  orog_drag_param = orog_drag_param_rp(rp_idx)
end if

!------------------------------------------------------------------
! 0. Set stresses to zero
!------------------------------------------------------------------
!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP private(i,j,l,height_fac,wta,wtb,rib_fn,zeta,tausx,tausy,k,              &
!$OMP         fp_x_low, fp_y_low,fp_x_steep,fp_y_steep)                        &
!$OMP SHARED(bl_levels,tdims,tau_fd_x, tau_fd_y,land_pts,land_index,h_m,       &
!$OMP        k_for_buoy,u_hm,v_hm,tl_hm,qw_hm,rib,z_uv,v_p,fd_stability_dep,   &
!$OMP        z_tq,tl,qw,db_surf,rib_surf,bt_gb,grcp,g,z0m,fd_hill_option,      &
!$OMP        fp_x, rho_wet_tq, sil_orog_land,fp_y,pdims,zh,u_p,bq_gb,          &
!$OMP        orog_drag_param,land_index_i,land_index_j)

!$OMP do SCHEDULE(STATIC)
do k = 1, bl_levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      tau_fd_x(i,j,k) = zero
      tau_fd_y(i,j,k) = zero
    end do
  end do ! land_pts
end do ! bl_levels
!$OMP end do

! The rest of the routine is only interested in land points.
if (land_pts > 0) then

  !----------------------------------------------------------------
  ! 1. Calculate the height scale h_m and interpolate the wind and
  !    density to this height.
  !----------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
  do l = 1, land_pts
    j = (land_index(l)-1)/pdims%i_end + 1
    i = land_index(l) - (j-1)*pdims%i_end
    land_index_i(l) = i
    land_index_j(l) = j

    h_m(l)=min(max_ht_scale,fd_decay*zh(i,j))
    h_m(l)=max(min_ht_scale,h_m(l))
    k_for_buoy(l) = 0
    u_hm(l) =zero
    v_hm(l) =zero
    tl_hm(l)=zero
    qw_hm(l)=zero
    rib(l)  =zero
  end do ! land_pts
!$OMP end do

  ! Interpolate to get U and V at z=h_m

  do k = 2, bl_levels
!$OMP do SCHEDULE(STATIC)
    do l = 1, land_pts

      i = land_index_i(l)
      j = land_index_j(l)

      if (h_m(l) <= z_uv(i,j,k) .and. h_m(l) >= z_uv(i,j,k-1)) then

        k_for_buoy(l) = k-1   ! theta-level below that containing h_m
        wta = ( h_m(l) - z_uv(i,j,k-1) )                                       &
                 /( z_uv(i,j,k) - z_uv(i,j,k-1) )
        wtb = ( z_uv(i,j,k) - h_m(l) )                                         &
                 /( z_uv(i,j,k) - z_uv(i,j,k-1) )
        u_hm(l) = wta*u_p(i,j,k) + wtb*u_p(i,j,k-1)
        v_hm(l) = wta*v_p(i,j,k) + wtb*v_p(i,j,k-1)

      end if

    end do ! land_pts
!$OMP end do
  end do ! bl_levels

  if ( fd_stability_dep == use_bulk_ri ) then
    do k = 2, bl_levels
!$OMP do SCHEDULE(STATIC)
      do l = 1, land_pts
        i = land_index_i(l)
        j = land_index_j(l)
        if ( h_m(l)<=z_tq(i,j,k) .and. h_m(l)>=z_tq(i,j,k-1) ) then
          wta = ( h_m(l) - z_tq(i,j,k-1) )                                     &
                      /( z_tq(i,j,k) - z_tq(i,j,k-1) )
          wtb = ( z_tq(i,j,k) - h_m(l)   )                                     &
                      /( z_tq(i,j,k) - z_tq(i,j,k-1) )
          tl_hm(l) = wta*tl(i,j,k) + wtb*tl(i,j,k-1)
          qw_hm(l) = wta*qw(i,j,k) + wtb*qw(i,j,k-1)
        end if
      end do ! land_pts
!$OMP end do
    end do ! bl_levels

!$OMP do SCHEDULE(STATIC)
    do l = 1, land_pts
      if (k_for_buoy(l) > 1) then
        k = k_for_buoy(l)
        i = land_index_i(l)
        j = land_index_j(l)
        db_surf(l) = rib_surf(i,j)*z_tq(i,j,1)*                                &
                     (u_p(i,j,1)*u_p(i,j,1)+v_p(i,j,1)*v_p(i,j,1))/            &
                     (z_uv(i,j,1)*z_uv(i,j,1))
        rib(l) = h_m(l)*(                                                      &
            g*( bt_gb(i,j,k)*( tl_hm(l)-tl(i,j,1)+                             &
                               grcp*(h_m(l)-z_tq(i,j,1)) )                     &
              + bq_gb(i,j,k)*(qw_hm(l)-qw(i,j,1)) ) + db_surf(l) )/            &
                      max(1.0e-12_r_bl, (u_hm(l)*u_hm(l)+v_hm(l)*v_hm(l)) )
      end if
    end do ! land_pts
!$OMP end do
  end if

  !-----------------------------------------------------------------------
  ! 2. Calculate the pressure force from the wind components and density
  !    at height h_m, the frontal silhouette area and surface roughness
  !    length for momentum, z_0m.
  !-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
  do l = 1, land_pts

    i = land_index_i(l)
    j = land_index_j(l)

    if ( fd_stability_dep >= on ) then
      if ( fd_stability_dep == use_bulk_ri ) then
        rib_fn=one-rib(l)/ri_crit
      else
        rib_fn=one-rib_surf(i,j)/ri_crit
      end if
      if (rib_fn >  one)rib_fn=one
      if (rib_fn <  zero)rib_fn=zero
    else
      rib_fn=one
    end if

    zeta = log( h_m(l)/z0m(i,j) )

    if (fd_hill_option == low_hill) then

      ! Compute Wood and Mason (1993) low-hill drag expression

      tausx=(vkman/zeta)*(vkman/zeta)*u_hm(l)*                                 &
            sqrt(u_hm(l)*u_hm(l)+v_hm(l)*v_hm(l))
      tausy=(vkman/zeta)*(vkman/zeta)*v_hm(l)*                                 &
            sqrt(u_hm(l)*u_hm(l)+v_hm(l)*v_hm(l))
      fp_x(l) = rho_wet_tq(i,j,1)*alpha*beta*pi_squared                        &
              *sil_orog_land(l)*sil_orog_land(l)                               &
              *rib_fn*tausx
      fp_y(l) = rho_wet_tq(i,j,1)*alpha*beta*pi_squared                        &
              *sil_orog_land(l)*sil_orog_land(l)                               &
              *rib_fn*tausy
    else if (fd_hill_option == steep_hill) then

      ! Compute steep hill drag expression

      tausx=u_hm(l)*sqrt(u_hm(l)*u_hm(l)+v_hm(l)*v_hm(l))
      tausy=v_hm(l)*sqrt(u_hm(l)*u_hm(l)+v_hm(l)*v_hm(l))
      fp_x(l)=one_half*rho_wet_tq(i,j,1)*orog_drag_param*                      &
              sil_orog_land(l)*rib_fn*tausx
      fp_y(l)=one_half*rho_wet_tq(i,j,1)*orog_drag_param*                      &
              sil_orog_land(l)*rib_fn*tausy

    else if (fd_hill_option == capped_lowhill) then

      ! Compute steep hill drag expression

      tausx=u_hm(l)*sqrt(u_hm(l)*u_hm(l)+v_hm(l)*v_hm(l))
      tausy=v_hm(l)*sqrt(u_hm(l)*u_hm(l)+v_hm(l)*v_hm(l))
      fp_x_steep=one_half*rho_wet_tq(i,j,1)*orog_drag_param*                   &
                 sil_orog_land(l)*rib_fn*tausx
      fp_y_steep=one_half*rho_wet_tq(i,j,1)*orog_drag_param*                   &
                 sil_orog_land(l)*rib_fn*tausy

      ! Compute Wood and Mason (1993) low-hill drag expression

      tausx=(vkman/zeta)*(vkman/zeta)*tausx
      tausy=(vkman/zeta)*(vkman/zeta)*tausy
      fp_x_low = rho_wet_tq(i,j,1)*alpha*beta*pi_squared                       &
              *sil_orog_land(l)*sil_orog_land(l)                               &
              *rib_fn*tausx
      fp_y_low = rho_wet_tq(i,j,1)*alpha*beta*pi_squared                       &
              *sil_orog_land(l)*sil_orog_land(l)                               &
              *rib_fn*tausy

      ! Take the minimum - effectively caps WM93 for large A/S or large z0

      fp_x(l) = min(fp_x_low, fp_x_steep)
      fp_y(l) = min(fp_y_low, fp_y_steep)

    end if

    tau_fd_x(i,j,1) = fp_x(l)
    tau_fd_y(i,j,1) = fp_y(l)

  end do ! land_pts
!$OMP end do

  !-----------------------------------------------------------------------
  ! 3. Calculate the vertical profiles of the explicit orographic stress
  !-----------------------------------------------------------------------
!$OMP do SCHEDULE(STATIC)
  do k = 2, bl_levels
    do l = 1, land_pts

      i = land_index_i(l)
      j = land_index_j(l)

      ! Limit the exponent so that height_fac doesn't get too big for
      ! single precision
      height_fac = exp(min(z_tq(i,j,k-1)/h_m(l),80.0_r_bl))
      tau_fd_x(i,j,k) = tau_fd_x(i,j,1)/height_fac
      tau_fd_y(i,j,k) = tau_fd_y(i,j,1)/height_fac

    end do ! land_pts
  end do ! bl_levels
!$OMP end do

end if ! land_pts > 0

!$OMP end PARALLEL

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine fm_drag
end module fm_drag_mod
