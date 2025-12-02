! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module conv_diag_6a_mod

use planet_constants_mod, only: vkman => vkman_bl, cp => cp_bl,                &
     kappa => kappa_bl, r => rd_bl, repsilon => repsilon_bl

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName='CONV_DIAG_6A_MOD'

contains
  !
  !   to diagnose convective occurrence and type
  !
subroutine conv_diag_6a(                                                       &

   ! in values defining field dimensions and subset to be processed :
    row_length, rows                                                           &

   ! in values defining vertical grid of model atmosphere :
   , bl_levels                                                                 &
   , p, P_theta_lev, exner_rho                                                 &
   , rho_only, rho_theta, z_full, z_half, r_theta_levels, r_rho_levels         &

   ! in Model switches
   , l_extra_call                                                              &
   , no_cumulus                                                                &

   ! in Cloud data :
   , qcf, qcl, cloud_fraction                                                  &

   ! in everything not covered so far :
   , pstar, q, theta, exner_theta_levels, u_p, v_p, u_0_p, v_0_p               &
   , tstar_land, tstar_sea, tstar_sice, z0msea                                 &
   , flux_e, flux_h, ustar_in, L_spec_z0, z0m_scm, z0h_scm                     &
   , tstar, land_mask, frac_land, ice_fract                                    &
   , w_copy, w_max                                                             &
   , conv_prog_precip                                                          &
   , g_ccp, h_ccp, ccp_strength                                                &

   ! in surface fluxes
   , fb_surf, ustar                                                            &

   ! SCM Diagnostics (dummy values in full UM)
   , nscmdpkgs, l_scmdiags                                                     &

   ! out data required elsewhere in UM system :

   , zh,zhpar,dzh,qcl_inv_top,z_lcl,z_lcl_uv,delthvu,ql_ad                     &
   , ntml,ntpar,nlcl                                                           &
   , cumulus, l_shallow,l_congestus, l_congestus2, conv_type, cin              &
   , cape, wstar, wthvs, entrain_coef, qsat_lcl                                &
   , error_code                                                                &
   ! tnuc_dust fields
   , tnuc_new, tnuc_nlcl)

  !
  ! Purpose:
  !    To diagnose convection occurrence and type
  !
  ! Code Owner: Please refer to the UM file CodeOwners.txt
  ! This file belongs in section: Convection
  !
  ! Code Description:
  !  Language: FORTRAN90
  !  This code is written to UMDP3 v6 programming standards
  !
  !-----------------------------------------------------------------------

! Definitions of prognostic variable array sizes
use atm_fields_bounds_mod, only:                                               &
  pdims, pdims_s, tdims_s, tdims, wdims, tdims_l


use bl_option_mod, only: one_third, zero, one
use cv_run_mod, only:                                                          &
   icvdiag, l_jules_flux, cnv_cold_pools, phi_ccp, l_ccp_blv
use cv_param_mod, only: ccp_off

use planet_constants_mod, only: g => g_bl, c_virtual => c_virtual_bl
use missing_data_mod,    only: rmdi
use model_domain_mod,    only: model_type, mt_single_column
use s_scmop_mod,         only: default_streams,                                &
                               t_inst, d_wet, d_point,                         &
                               scmdiag_bl, scmdiag_conv

! subroutines
use conv_surf_flux_mod, only: conv_surf_flux
use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook

use conv_diag_comp_6a_mod, only: conv_diag_comp_6a
use mphys_inputs_mod,    only: l_progn_tnuc
implicit none

!-------------------------------------------------------------------
! Subroutine Arguments
!-------------------------------------------------------------------
!
! Arguments with intent in:
!
! (a) Defining horizontal grid and subset thereof to be processed.

integer, intent(in) ::                                                         &
   row_length                                                                  &
                            ! Local number of points on a row
   , rows
                            ! Local number of rows in a v field

! (b) Defining vertical grid of model atmosphere.


integer, intent(in) ::                                                         &
     bl_levels
                            ! Max. no. of "boundary" levels allowed.

real(kind=r_bl), intent(in) ::                                                 &
    p(pdims_s%i_start:pdims_s%i_end,    & ! pressure  on rho levels (Pa)
      pdims_s%j_start:pdims_s%j_end,                                           &
      pdims_s%k_start:pdims_s%k_end)                                           &
  , P_theta_lev(tdims%i_start:tdims%i_end,     & ! P on theta lev (Pa)
                tdims%j_start:tdims%j_end,                                     &
                            1:tdims%k_end)                                     &
  , exner_rho(pdims_s%i_start:pdims_s%i_end,   & ! Exner on rho level
              pdims_s%j_start:pdims_s%j_end,   & !
              pdims_s%k_start:pdims_s%k_end)                                   &
  , rho_only(row_length,rows,1:tdims%k_end)    & ! density (kg/m3)
  , rho_theta(row_length,rows,1:tdims%k_end-1) & ! rho th lev (kg/m3)
  , r_theta_levels(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,&
                   0:tdims%k_end) & ! dist of theta lev from Earth centre (m)
  , r_rho_levels(tdims_l%i_start:tdims_l%i_end,tdims_l%j_start:tdims_l%j_end,  &
                 tdims%k_end) & ! dist of rho lev from Earth centre (m)
  , z_full(row_length,rows,1:tdims%k_end)      & ! height th lev (m)
  , z_half(row_length,rows,1:tdims%k_end)        ! height rho lev (m)


logical,  intent(in) ::                                                        &
    L_spec_z0             & ! true if roughness length has been specified
  , L_extra_call            ! true this is an additional call to conv_diag
                            ! within a timestep

logical,  intent(in) ::                                                        &
   no_cumulus(row_length,rows)   ! Points overruled by BL

! (c) Cloud data.

real(kind=r_bl), intent(in) ::                                                 &
  qcf(tdims%i_start:tdims%i_end,             & ! Cloud ice (kg per kg air)
      tdims%j_start:tdims%j_end,                                               &
                  1:tdims%k_end)                                               &
 ,qcl(tdims%i_start:tdims%i_end,             & ! Cloud liquid water (kg/kg air)
      tdims%j_start:tdims%j_end,                                               &
                  1:tdims%k_end)                                               &
 ,cloud_fraction(tdims%i_start:tdims%i_end,  & ! Cloud fraction
                 tdims%j_start:tdims%j_end,                                    &
                             1:tdims%k_end)

    ! (d) Atmospheric + any other data not covered so far, incl control.

real(kind=r_bl), intent(in) ::                                                 &
  pstar(row_length, rows)                       & ! Surface pressure (Pa)
 ,q(tdims%i_start:tdims%i_end,                  & ! water vapour (kg/kg)
    tdims%j_start:tdims%j_end,                                                 &
                1:tdims%k_end)                                                 &
 ,theta(tdims%i_start:tdims%i_end,              & ! Theta (Kelvin)
        tdims%j_start:tdims%j_end,                                             &
                    1:tdims%k_end)                                             &
 ,exner_theta_levels(tdims%i_start:tdims%i_end, & ! exner pressure theta lev
                     tdims%j_start:tdims%j_end, & !  (Pa)
                                 1:tdims%k_end)
real(kind=r_bl), intent(in) ::                                                 &
    u_p(row_length, rows)                                                      &
                           ! U(1) on P-grid.
  , v_p(row_length, rows)                                                      &
                           ! V(1) on P-grid.
  , u_0_p(row_length,rows)                                                     &
                           ! W'ly component of surface current
                           !    (metres per second) on P-grid.
  , v_0_p(row_length,rows)                                                     &
                           ! S'ly component of surface current
                           !    (metres per second) on P-grid.
  ,flux_e(row_length,rows)                                                     &
                           ! Specified surface
                           !    latent heat flux (W/m^2)
  ,flux_h(row_length,rows)                                                     &
                           ! Specified surface
                           !    sensible heat fluxes (in W/m2)
  ,ustar_in(row_length,rows)                                                   &
                           ! Specified surface friction velocity (m/s)
  ,z0msea(row_length,rows)                                                     &
                           ! Sea roughness length for momentum (m)
  ,z0m_scm(row_length,rows)                                                    &
                           ! Namelist input z0m (if >0)
  ,z0h_scm(row_length,rows)! Namelist input z0h (if >0)

real(kind=r_bl), intent(in) ::                                                 &
 tstar_land(row_length, rows) & ! Surface T on land
,tstar_sea(row_length, rows)  & ! Surface T on sea
,tstar_sice(row_length, rows)   ! Surface T on sea-ice

real(kind=r_bl), intent(in out) ::                                             &
   tstar(row_length,rows)    ! Surface temperature (K).
                            ! NOTE only inout for SCM
real(kind=r_bl), intent(in) ::                                                 &
    tnuc_new( tdims%i_start:tdims%i_end,                                       &
             tdims%j_start:tdims%j_end,                                        &
             1:tdims%k_end )
                            ! ice nucln. temp. as function of dust (deg cel)
real(kind=r_bl), intent(in out) ::                                             &
    tnuc_nlcl(row_length, rows)
                            ! ice nucln. temp. as function of dust (deg cel)
                            ! indexed using nlcl
logical, intent(in) ::                                                         &
   land_mask(row_length, rows)  ! T If land, F Elsewhere.

real(kind=r_bl), intent(in) ::                                                 &
frac_land(pdims_s%i_start:pdims_s%i_end, & ! fraction of gridbox that is
        pdims_s%j_start:pdims_s%j_end)  & ! land, on all points
,ice_fract(row_length,rows)                 ! fraction of sea that has ice

real(kind=r_bl), intent(in) ::                                                 &
 w_copy(wdims%i_start:wdims%i_end,  & ! vertical velocity W (m/s)
     wdims%j_start:wdims%j_end,                                                &
                 0:wdims%k_end)  & ! Not exact match to module values
,w_max(row_length,rows)               ! Column max vertical velocity (m/s)

! Surface precipitation based 3d convective prognostic in kg/m2/s
real(kind=r_bl), intent(in) ::                                                 &
                    conv_prog_precip( tdims_s%i_start:tdims_s%i_end,           &
                                      tdims_s%j_start:tdims_s%j_end,           &
                                      tdims_s%k_start:tdims_s%k_end )

! convective cold-pool prognostics
real(kind=r_bl), intent (in) ::                                                &
                     g_ccp( pdims%i_start:pdims%i_end,                         &
                            pdims%j_start:pdims%j_end)                         &
!            gridbox c.c.p. reduced gravity (m/s^2)
!
                   , h_ccp( pdims%i_start:pdims%i_end,                         &
                            pdims%j_start:pdims%j_end)
!            gridbox c.c.p. height (m)

! diagnosed cold-pool strength
real(kind=r_bl), intent (in) ::                                                &
                     ccp_strength( pdims%i_start:pdims%i_end,                  &
                                   pdims%j_start:pdims%j_end)

real(kind=r_bl), intent(in) ::                                                 &
   ustar(tdims%i_start:tdims%i_end,tdims%j_start:tdims%j_end)
                                ! GBM surface friction velocity

! Additional variables for SCM diagnostics which are dummy in full UM
integer ::                                                                     &
   nscmdpkgs          ! No of diagnostics packages
!
logical ::                                                                     &
   l_scmdiags(nscmdpkgs) ! Logicals for diagnostics packages
!
real(kind=r_bl), intent(in out) ::                                             &
   zh(row_length,rows)    ! Height above surface of top
!  of boundary layer (metres).
real(kind=r_bl), intent(out) ::                                                &
   zhpar(row_length,rows)                                                      &
                            ! Height of max parcel ascent (m)
   ,dzh(row_length,rows)                                                       &
                            ! Height of inversion top (m)
   ,qcl_inv_top(row_length,rows)                                               &
                            ! Parcel water content at inversion top
   ,z_lcl(row_length,rows)                                                     &
                            ! Height of lifting condensation level (m)
   ,z_lcl_uv(row_length,rows)                                                  &
                            ! Height of lIfting condensation
                            !     level on uv grid (m)
   ,delthvu(row_length,rows)                                                   &
                            ! Integral of undilute parcel buoyancy
                            ! over convective cloud layer
                            ! (for convection scheme)
   ,ql_ad(row_length,rows)                                                     &
                            ! adiabatic liquid water content at
                            ! inversion or cloud top (kg/kg)
   ,entrain_coef(row_length,rows)                                              &
                            ! Entrainment coefficient
   ,qsat_lcl(row_length,rows)
                            ! qsat at cloud base (kg/kg)

real(kind=r_bl), intent(out) ::                                                &
     cape(row_length, rows)                                                    &
                            ! CAPE from parcel ascent (m2/s2)
   , cin(row_length, rows)
                            ! CIN from parcel ascent (m2/s2)

integer, intent(in out)  ::                                                    &
   ntml(row_length,rows)
                            ! Number of model levels in the
                            ! turbulently mixed layer.
integer, intent(out)  ::                                                       &
   ntpar(row_length,rows)                                                      &
                            ! Max levels for parcel ascent
   ,nlcl(row_length,rows)
                            ! No. of model layers below the
                            ! lifting condensation level.

! Convective type array ::
integer, intent(out) :: conv_type(row_length, rows)
                            ! Integer index describing convective type:
                            !    0=no convection
                            !    1=non-precipitating shallow
                            !    2=drizzling shallow
                            !    3=warm congestus
                            !    ...
                            !    8=deep convection

logical, intent(out) ::                                                        &
   cumulus(row_length,rows)                                                    &
                            ! Logical indicator for convection
   ,l_shallow(row_length,rows)                                                 &
                            ! Logical indicator for shallow Cu
   ,l_congestus(row_length,rows)                                               &
                            ! Logical indicator for congestus Cu
   ,l_congestus2(row_length,rows) ! Logical ind 2 for congestus Cu

! Required for extra call to conv_diag as values from original BL call
! may return zero values based on initial timestep profiles.
real(kind=r_bl), intent(in out) ::                                             &
  wstar(row_length,rows)  & ! Convective sub-cloud velocity scale (m/s)
 ,wthvs(row_length,rows)    ! surface flux of wthv (K m/s)

integer, intent(out)  ::                                                       &
   error_code     ! 0 - no error in this routine

!-----------------------------------------------------------------------
! Local variables
!-----------------------------------------------------------------------

character (len=*), parameter ::  RoutineName = 'CONV_DIAG_6A'

integer ::                                                                     &
   i,j           & ! LOCAL Loop counter (horizontal field index).
   , ii          & ! Local compressed array counter.
   , k           & ! LOCAL Loop counter (vertical level index).
   ,nunstable      ! total number of unstable points


! Declare arrays for scm diagnostics

! uncompressed arrays - all points

real(kind=r_bl) ::                                                             &
    fb_surf(row_length, rows)     ! Change in theta_v from surface
                                  ! to layer 1 (note diff from BL)


! Arrays added for extra conv_diag calls
integer ::                                                                     &
   ntml_copy(row_length,rows)

! compressed arrays store only values for unstable columns

integer ::                                                                     &
   index_i(row_length*rows)   & ! column number of unstable points
   , index_j(row_length*rows)     ! row number of unstable points

real(kind=r_bl)  ::                                                            &
   tv1_sd( row_length*rows)    & ! Approx to standard dev of level
                                 ! 1 virtual temperature (K).
   ,bl_vscale2(row_length*rows)& ! Velocity scale squared for
                                 ! boundary layer eddies (m2/s2)
   ,frac_land_c(row_length*rows) ! fraction of land compressed

real(kind=r_bl) :: w_m ! velocity scale
!
! Model constants:
!

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
! Mixing ratio, r,  versus specific humidity, q
!
! In most cases the expression to first order are the same
!
!  Tl = T - (lc/cp)qcl - [(lc+lf)/cp]qcf
!  Tl = T - (lc/cp)rcl - [(lc+lf)/cp]rcf  - equally correct definition
!
! thetav = theta(1+cvq)         accurate
!        = theta(1+r/repsilon)/(1+r) ~ theta(1+cvr) approximate
!
! svl = (Tl+gz/cp)*(1+(1/repsilon-1)qt)
!     ~ (Tl+gz/cp)*(1+(1/repsilon-1)rt)
!
! dqsat/dT = repsilon*Lc*qsat/(R*T*T)
! drsat/dT = repsilon*Lc*rsat/(R*T*T)  equally approximate
!
! Only altering the expression for vapour pressure
!
!  e = qp/repsilon       - approximation
!  e = rp/(repsilon+r)   - accurate
!-----------------------------------------------------------------------
! 1.0 Initialisation
!-----------------------------------------------------------------------

error_code = 0

!-----------------------------------------------------------------------
! 1.1 Verify grid/subset definitions.
!-----------------------------------------------------------------------

if ( bl_levels <  1 .or. rows <  1 .or. tdims%k_end <  1 ) then
  error_code = 1
  return

end if
!-----------------------------------------------------------------------
! 1.1a copy ntml values where do not want to overwrite
!-----------------------------------------------------------------------
if (l_extra_call) then
  do j=1, rows
    do i=1,row_length
      if (no_cumulus(i,j)) then
        ntml_copy(i,j) = ntml(i,j)
      end if
    end do
  end do
end if

!-----------------------------------------------------------------------
! 1.2 initialisation of output arrays
!-----------------------------------------------------------------------

do j=1,rows
  do i=1,row_length

    cumulus(i,j)      = .false.
    l_shallow(i,j)    = .false.
    l_congestus(i,j)  = .false.
    l_congestus2(i,j) = .false.
    conv_type(i,j) = 0
    ntml(i,j) = 1
    nlcl(i,j) = 1
    ntpar(i,j) = 1
    delthvu(i,j)  = zero
    ql_ad(i,j)    = zero
    cape(i,j)     = zero
    cin(i,j)      = zero
    z_lcl(i,j)    = zero       ! set to zero
    zhpar(i,j)    = zero
    dzh(i,j)      = rmdi      ! initialise to missing data
    qcl_inv_top(i,j)= zero

    ! Set LCL for UV grid to level below z_lcl (note that (at vn5.1) UV
    ! levels are half levels (below) relative to P-grid. Consistent with
    ! BL scheme's treatment of staggered vertical grid.)

    z_lcl_uv(i,j)= z_full(i,j,nlcl(i,j))

    entrain_coef(i,j) = -99.0_r_bl    ! indicates not set

  end do
end do

!-----------------------------------------------------------------------
! 1.4 Calculation of surface buoyancy flux and level one standard
!  deviation of virtual temperature.
!  Also modifies tstar for some SCM runs.
!-----------------------------------------------------------------------

if (l_jules_flux) then

  ii = 0

  do j=1, rows
    do i=1, row_length
      if (fb_surf(i,j) > zero) then
        ii= ii+1
        w_m = ( 0.25_r_bl * zh(i,j) * fb_surf(i,j) +                           &
                ustar(i,j) * ustar(i,j) * ustar(i,j) ) ** one_third
        tv1_sd(ii) = theta(i,j,1) *                                            &
                    (one+c_virtual*q(i,j,1)-qcl(i,j,1)-qcf(i,j,1)) *           &
                    ( 1.93_r_bl * fb_surf(i,j) / ( w_m * g) )
                     ! tv1_sd = 1.93 * wthvbar / w_m
                     ! fb_surf = g * wthvbar / thv
        bl_vscale2(ii) = max(0.01_r_bl, 2.5_r_bl * 2.52_r_bl * w_m * w_m)
                       ! 2.5 from Beare (2008), 2.52 to undo the 0.25 factor
      end if
    end do
  end do

else
end if
!-----------------------------------------------------------------------
! 2.0 Decide on unstable points  ( fb_surf > 0.0 )
!     Only work on these points for the rest of the calculations.
!-----------------------------------------------------------------------

nunstable = 0           ! total number of unstable points
do j=1,rows
  do i=1,row_length
    if ( fb_surf(i,j)  >   zero ) then
      nunstable = nunstable + 1
      index_i(nunstable) = i
      index_j(nunstable) = j
    end if
  end do
end do

! land fraction on just unstable points
if (nunstable > 0) then
  do ii=1,nunstable
    i = index_i(ii)
    j = index_j(ii)
    frac_land_c(ii) = frac_land(i,j)
  end do

  ! add cold-pool contribution to the TKE velocity scale for BL eddies
  if ((cnv_cold_pools > ccp_off) .and. l_ccp_blv) then
    do ii=1,nunstable
      i = index_i(ii)
      j = index_j(ii)
      bl_vscale2(ii) = bl_vscale2(ii) + phi_ccp*g_ccp(i,j)*h_ccp(i,j)
    end do
  end if

  !-----------------------------------------------------------------------
  ! Work on just unstable points to diagnose whether convective.
  !-----------------------------------------------------------------------
  call conv_diag_comp_6a(                                                      &
    row_length, rows                                                           &
  , bl_levels, nunstable                                                       &
  , index_i, index_j                                                           &
  , p, P_theta_lev, exner_rho                                                  &
  , rho_only, rho_theta, z_full, z_half, r_theta_levels, r_rho_levels          &
  , qcf, qcl, cloud_fraction                                                   &
  , pstar, q, theta, exner_theta_levels                                        &
  , frac_land_c                                                                &
  , w_copy ,w_max, tv1_sd, bl_vscale2                                          &
  , conv_prog_precip                                                           &
  , nSCMDpkgs, L_SCMDiags                                                      &
  , ntml, ntpar, nlcl                                                          &
  , cumulus, L_shallow, l_congestus, l_congestus2, conv_type                   &
  , zh,zhpar,dzh,qcl_inv_top,z_lcl,z_lcl_uv,delthvu,ql_ad                      &
  , cin, cape,entrain_coef                                                     &
  , qsat_lcl                                                                   &
    )

else

  ! Reset zh as done in conv_diag_comp for all points
  do j=1,rows
    do i=1,row_length
      zh(i,j) = z_half(i,j,2)
    end do
  end do

end if        ! unstable points only
!=======================================================================
! using nlcl to index tnuc_new !
if (l_progn_tnuc) then
  do j=1,rows
    do i=1,row_length
      tnuc_nlcl(i,j) = tnuc_new(i,j,nlcl(i,j))
    end do
  end do
end if
 !=======================================================================
 ! Extra calculation required if conv_diag called more than once per
 ! model physics timestep
 !=======================================================================
if (l_extra_call) then

  do j=1, rows
    do i=1,row_length
      if (fb_surf(i,j) > zero) then
        wstar(i,j) = (zh(i,j)*fb_surf(i,j))**one_third
        wthvs(i,j) = fb_surf(i,j)*theta(i,j,1)                                 &
           *exner_theta_levels(i,j,1)/g
      else
        wstar(i,j) = zero
        wthvs(i,j) = zero
      end if
      ! BL overruled first sweep diagnosis of cumulus for a good reason
      ! So prevent subsequent sweeps from diagnosing convection
      if (no_cumulus(i,j)) then
        cumulus(i,j)   = .false.
        L_shallow(i,j) = .false.
        conv_type(i,j) = 0

        ! Copy back orignal ntml value
        ntml(i,j) = ntml_copy(i,j)

      end if
    end do
  end do

end if        ! test on l_extra_call



!-----------------------------------------------------------------------
!  SCM Diagnostics for no unstable points
!  Unstable diagnostics are output from conv_diag_comp but this routine
!  is not called for stable points so null (zero values) are output from
!  this routine for stable points.
!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine conv_diag_6a

end module conv_diag_6a_mod
