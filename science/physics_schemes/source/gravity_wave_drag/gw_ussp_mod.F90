! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Gravity Wave (Ultra-Simple Spectral Parametrization) Scheme.
! Subroutine Interface:

module gw_ussp_mod

use um_types, only: real_umphys, real_usprec

implicit none

character(len=*), parameter, private :: ModuleName='GW_USSP_MOD'

contains

subroutine gw_ussp(levels, rows,                                               &
  row_length,                                                                  &
  global_row_length,                                                           &
  r_rho_levels, r_theta_levels, p_layer_boundaries,                            &
  sin_theta_latitude,                                                          &
  theta, rho, u, v, totalppn, timestep,                                        &
  r_u, r_v, T_inc,                                                             &
  L_ussp_heating,                                                              &
  gwspec_eflux,gwspec_sflux,gwspec_wflux,gwspec_nflux,                         &
  gwspec_ewacc,gwspec_nsacc,                                                   &
  gwspec_eflux_on, gwspec_eflux_p_on, gwspec_sflux_on,                         &
  gwspec_wflux_on, gwspec_wflux_p_on, gwspec_nflux_on,                         &
  gwspec_ewacc_on, gwspec_ewacc_p_on, gwspec_nsacc_on)
!
! purpose:    This subroutine calculates the vertical momentum flux
!             divergence due to gravity waves as parametrised by the
!             Warner and McIntyre Ultra Simple Spectral gravity wave
!             Parametrization adapted for use in the UM.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: gravity_wave_drag
!
! code description:
!   language: fortran 90
!   this code is written to umdp3 programming standards.
!   documentation: Unified Model Documentation Paper 34 (Non-Orog GW)

! Processor addresses (N, S)

use g_wave_input_mod, only: L_add_cgw

use gw_ussp_prec_mod, only: ussp_launch_factor, wavelstar, cgw_scale_factor,   &
                            two_omega, gw_ussp_prec_set_reals

use planet_constants_mod, only: g, r, cp,                                      &
                                planet_radius, recip_a2
use gw_ussp_params_mod, only: idir

use tuning_segments_mod, only: l_autotune_segments, ussp_seg_size

! Definitions of prognostic variable array sizes
use atm_fields_bounds_mod, only:                                               &
   udims, vdims, pdims, tdims, pdims,                                          &
   udims_s, vdims_s, pdims_s, tdims_s, pdims_l, tdims_l

! Model level values
use level_heights_mod, only:                                                   &
   eta_theta_levels           ! Eta values of theta levels

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
!$    use omp_lib, only: omp_get_num_threads
use model_domain_mod, only: model_type, mt_global

use segments_mod, only:                                                        &
  segment_type, meta_segment_type,                                             &
  segments_mod_seg_meta, segments_mod_segments

! ----------------------------------------------------------------------+-------
! Subroutines defined from modules
! ----------------------------------------------------------------------+-------
use gw_ussp_core_mod, only: gw_ussp_core

use p_to_t_mod, only: p_to_t
use u_to_p_mod, only: u_to_p
use v_to_p_mod, only: v_to_p

! ----------------------------------------------------------------------+-------
! Subroutines that are not used when compiling in LFRic
! ----------------------------------------------------------------------+-------


! ----------------------------------------------------------------------+-------

implicit none

! ----------------------------------------------------------------------+-------
! Subroutine arguments of GW_USSP
! ----------------------------------------------------------------------+-------
!     Fixed starting theta-level (e.g. for P_layer_boundaries)
integer, parameter :: tkfix0start    = 0
!
!     Fixed starting theta-level (e.g. for STASH diagnostic arrays)
integer, parameter :: tkfix1start    = 1
!
! Number of model levels
integer, intent(in) :: levels
!
! Number of rows for u field
integer, intent(in) :: rows
!
! Number of points on a row of a full global field
integer, intent(in) :: global_row_length
!
! Number of grid points in row
integer, intent(in) :: row_length
! ----------------------------------------------------------------------+-------
!
! Timestep
real(kind=real_umphys), intent(in)    :: timestep
!
! P-GRID Latitudes
real(kind=real_umphys), intent(in)    ::                                       &
                       sin_theta_latitude(tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end)
!
! Centre of earth distance to theta-level points
real(kind=real_umphys), intent(in)    ::                                       &
                       r_theta_levels(tdims_l%i_start:tdims_l%i_end,           &
                                      tdims_l%j_start:tdims_l%j_end,           &
                                                    0:tdims_l%k_end)
!
! Centre of earth distance to rho-level points
real(kind=real_umphys), intent(in)    ::                                       &
                       r_rho_levels(pdims_l%i_start:pdims_l%i_end,             &
                                    pdims_l%j_start:pdims_l%j_end,             &
                                                    pdims_l%k_end)
!
! Pressure on layer boundaries
real(kind=real_umphys), intent(in)    ::                                       &
                       p_layer_boundaries(tdims%i_start:tdims%i_end,           &
                                          tdims%j_start:tdims%j_end,           &
                                            tkfix0start:tdims%k_end)
!
! Total Precipitation for CGW source
real(kind=real_umphys), intent(in)    ::                                       &
                       totalppn(tdims%i_start:tdims%i_end,                     &
                                     tdims%j_start:tdims%j_end)
!
! Primary model array for theta
real(kind=real_umphys), intent(in)    ::                                       &
                       theta(tdims_s%i_start:tdims_s%i_end,                    &
                             tdims_s%j_start:tdims_s%j_end,                    &
                                 tkfix1start:tdims_s%k_end)
!
! Primary model array for Density x (radius earth)^2.
real(kind=real_umphys), intent(in)    ::                                       &
                       rho(pdims_s%i_start:pdims_s%i_end,                      &
                           pdims_s%j_start:pdims_s%j_end,                      &
                           pdims_s%k_start:pdims_s%k_end)
!
! Primary model array for U-wind field
real(kind=real_umphys), intent(in)    ::                                       &
                       u(udims_s%i_start:udims_s%i_end,                        &
                         udims_s%j_start:udims_s%j_end,                        &
                         udims_s%k_start:udims_s%k_end)
!
! Primary model array for V-wind field
real(kind=real_umphys), intent(in)    ::                                       &
                       v(vdims_s%i_start:vdims_s%i_end,                        &
                         vdims_s%j_start:vdims_s%j_end,                        &
                         vdims_s%k_start:vdims_s%k_end)
!
! FLux Fp in E azimuth
real(kind=real_umphys), intent(in out) ::                                      &
                       gwspec_eflux(udims%i_start:udims%i_end,                 &
                                    udims%j_start:udims%j_end,                 &
                                      tkfix1start:tdims%k_end)
!
! FLux Fp in S azimuth
real(kind=real_umphys), intent(in out) ::                                      &
                       gwspec_sflux(vdims%i_start:vdims%i_end,                 &
                                    vdims%j_start:vdims%j_end,                 &
                                      tkfix1start:tdims%k_end)
!
! FLux Fp in W azimuth
real(kind=real_umphys), intent(in out) ::                                      &
                       gwspec_wflux(udims%i_start:udims%i_end,                 &
                                    udims%j_start:udims%j_end,                 &
                                      tkfix1start:tdims%k_end)
!
! Fp in N azimuth
real(kind=real_umphys), intent(in out) ::                                      &
                       gwspec_nflux(vdims%i_start:vdims%i_end,                 &
                                    vdims%j_start:vdims%j_end,                 &
                                      tkfix1start:tdims%k_end)
!
! Accel of U wind
real(kind=real_umphys), intent(in out) ::                                      &
                       gwspec_ewacc(udims%i_start:udims%i_end,                 &
                                    udims%j_start:udims%j_end,                 &
                                    udims%k_start:udims%k_end)
!
! Accel of V wind
real(kind=real_umphys), intent(in out) ::                                      &
                       gwspec_nsacc(vdims%i_start:vdims%i_end,                 &
                                    vdims%j_start:vdims%j_end,                 &
                                    vdims%k_start:vdims%k_end)
!
! Temperature increment
real(kind=real_umphys), intent(in out) ::                                      &
                       T_inc(tdims%i_start:tdims%i_end,                        &
                             tdims%j_start:tdims%j_end,                        &
                               tkfix1start:tdims%k_end)
!
! U-wind increment diagnostic
real(kind=real_umphys), intent(in out) ::                                      &
                       r_u(udims_s%i_start:udims_s%i_end,                      &
                           udims_s%j_start:udims_s%j_end,                      &
                           udims_s%k_start:udims_s%k_end)
!
! V-wind increment diagnostic
real(kind=real_umphys), intent(in out) ::                                      &
                       r_v(vdims_s%i_start:vdims_s%i_end,                      &
                           vdims_s%j_start:vdims_s%j_end,                      &
                           vdims_s%k_start:vdims_s%k_end)
! ----------------------------------------------------------------------+-------
!
! Switch to calculate heating tendency
logical, intent(in) :: L_ussp_heating
!
! Switch for diagnostic of Fp in azimuthal direction East
logical, intent(in) :: gwspec_eflux_on
!
! Switch for diagnostic of Fp in azimuthal direction East interpolated to p-levs
logical, intent(in) :: gwspec_eflux_p_on
!
! Switch for diagnostic of Fp in azimuthal direction South
logical, intent(in) :: gwspec_sflux_on
!
! Switch for diagnostic of Fp in azimuthal direction West
logical, intent(in) :: gwspec_wflux_on
!
! Switch for diagnostic of Fp in azimuthal direction West interpolated to p-levs
logical, intent(in) :: gwspec_wflux_p_on
!
! Switch for diagnostic of Fp in azimuthal direction North
logical, intent(in) :: gwspec_nflux_on
!
!Switch for net Eastward acceleration diagnostic
logical, intent(in) :: gwspec_ewacc_on
!
!Switch for net Eastward acceleration diagnostic interpolated to p-levels
logical, intent(in) :: gwspec_ewacc_p_on
!
!Switch for net Northward acceleration diagnostic
logical, intent(in) :: gwspec_nsacc_on
!
!
! ----------------------------------------------------------------------+-------
! Local parameters
! ----------------------------------------------------------------------+-------
!
! Eta (hybrid height) level for model launch
real(kind=real_umphys), parameter :: etalaunch         = 0.045

! Constant radius height level for model launch (m):
! In future, etalauch should be calculated by determining the
! nearest model level to zlaunch.
! real(kind=real_umphys), parameter :: zlaunch           = 3825.0
!
! ----------------------------------------------------------------------+-------
! Security parameters
! ----------------------------------------------------------------------+-------
!
! Minimum allowed value of buoyancy frequency squared
real(kind=real_umphys), parameter ::  sqnmin           = 1.0e-4
!
! ----------------------------------------------------------------------+-------
! Local Constants (a) Physical
! ----------------------------------------------------------------------+-------
!
! Cos(phi_j) - azimuthal direction
real(kind=real_umphys)            :: cosphi(4)
!
! Sin(phi_j) - azimuthal direction
real(kind=real_umphys)            :: sinphi(4)
!
! Reciprocal of middle atmosphere mean scale height for pressure,
! conventionally assumed to be around 7km: 1 / H = g / (R * 239.145)
real(kind=real_umphys) :: rscale_h
!
! ----------------------------------------------------------------------+-------
! Local variables (scalars) used in GW_USSP
! Some effectively expanded to workspace (using vector registers)
! ----------------------------------------------------------------------+-------
!
! Longitude index
integer         :: i
!
! Latitude index
integer         :: j
!
! Level index
integer         :: k
!
! Azimuthal direction index
integer         :: jdir
!
! omp block iterator
integer         :: jj
!
!
! Level for gw launch - in principle this could be launchlev(i,j)
integer         :: launchlev
!
! Minimum level in launchlev(i,j)
integer         :: minlaunchlev
!
! Maximum level in launchlev(i,j)
integer         :: maxlaunchlev
integer         :: ii             ! Looper
integer         :: ss             ! Addressing in segment
integer         :: num_parallel   ! Segmentation variable
integer         :: ipar           ! Segmentation variable
integer         :: num_segments   ! Number of segments
!

! Wave-induced force per unit mass due to azimuthal sectors (m s^-2)
real(kind=real_umphys)            :: g_g
!
! Layer depths for heating calculations
real(kind=real_umphys)            :: dzb, dzu, dzl
!
! Velocities on theta levels
real(kind=real_umphys)            :: uhat, vhat
!
! Kinetic energy tendencies on theta levels
real(kind=real_umphys)            :: ududt, vdvdt
!
! ----------------------------------------------------------------------+-------
! Local variables (dynamic arrays) used in GW_USSP
! ----------------------------------------------------------------------+-------
!
! Pseudomomentum flux integrated over azimuthal sector (kg m^-1 s^-2)
real(kind=real_umphys)            ::                                           &
                   fptot(tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end,                            &
                           tkfix1start:tdims%k_end,idir)
!
! Zonal component of wave-induced force on density level (m s^-2)
! Note that in our notation G_X equates to DU_DT, the zonal wind tendency
real(kind=real_umphys)            ::                                           &
                   g_x(pdims%i_start:pdims%i_end,                              &
                       pdims%j_start:pdims%j_end,                              &
                       pdims%k_start:pdims%k_end)
!
! Meridional component of wave-induced force on density level (m s^-2)
! Note that in our notation G_Y equates to DV_DT, meridional wind tendency
real(kind=real_umphys)            ::                                           &
                   g_y(pdims%i_start:pdims%i_end,                              &
                       pdims%j_start:pdims%j_end,                              &
                       pdims%k_start:pdims%k_end)
!
! Buoyancy [Brunt Vaisala] frequency on half-levels (rad s^-1)
real(kind=real_umphys)            ::                                           &
                   nbv(tdims%i_start:tdims%i_end,                              &
                       tdims%j_start:tdims%j_end,                              &
                         tkfix1start:tdims%k_end)
!
! Rho on theta levels
real(kind=real_umphys)            ::                                           &
                   rho_th(tdims%i_start:tdims%i_end,                           &
                          tdims%j_start:tdims%j_end,                           &
                            tkfix1start:tdims%k_end)
!
! Component of wind in phi_jdir azimuth direction (m s^-1)
real(kind=real_umphys)            ::                                           &
                   udotk(tdims%i_start:tdims%i_end,                            &
                         tdims%j_start:tdims%j_end,                            &
                           tkfix1start:tdims%k_end,idir)
!
! Zonal wind U interpolated onto Rho grid
real(kind=real_umphys)            ::                                           &
                   uonp(pdims%i_start:pdims%i_end,                             &
                        pdims%j_start:pdims%j_end,                             &
                        pdims%k_start:pdims%k_end)
!
! Meridional wind V interpolated onto Rho grid
real(kind=real_umphys)            ::                                           &
                   vonp(pdims%i_start:pdims%i_end,                             &
                        pdims%j_start:pdims%j_end,                             &
                        pdims%k_start:pdims%k_end)
!
! G_X with halo
real(kind=real_umphys)            ::                                           &
                   g_xp_smallhalo(pdims_s%i_start:pdims_s%i_end,               &
                                  pdims_s%j_start:pdims_s%j_end,               &
                                  pdims_s%k_start:pdims_s%k_end)
!
! G_Y with halo
real(kind=real_umphys)            ::                                           &
                   g_yp_smallhalo(pdims_s%i_start:pdims_s%i_end,               &
                                  pdims_s%j_start:pdims_s%j_end,               &
                                  pdims_s%k_start:pdims_s%k_end)
!
! Rho*(radius_earth)^2 on Theta grid with halo
real(kind=real_umphys)            ::                                           &
                   rhont_smallhalo(tdims_s%i_start:tdims_s%i_end,              &
                                   tdims_s%j_start:tdims_s%j_end,              &
                                       tkfix1start:tdims_s%k_end)
!
real(kind=real_usprec), allocatable :: s_sin_theta_lat(:)
real(kind=real_usprec), allocatable :: s_totalppn(:)
real(kind=real_usprec), allocatable :: s_rho_th(:,:)
real(kind=real_usprec), allocatable :: s_nbv(:,:)
real(kind=real_usprec), allocatable :: s_udotk(:,:,:)
real(kind=real_usprec), allocatable :: s_fptot(:,:,:)

! Variables for segmentation
type(segment_type),      allocatable  :: segments(:)
type(meta_segment_type)               :: meta_segments

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle


character(len=*), parameter :: RoutineName='GW_USSP'

!
!  End of Header
!
! ==Main Block==--------------------------------------------------------+-------
data cosphi/0,-1,0,1/
data sinphi/1,0,-1,0/
!
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call gw_ussp_prec_set_reals()
! ----------------------------
! Local Constants (a) Physical
! ----------------------------
rscale_h = g / (r * 239.145)
!
! ----------------------------------------------------------------------+-------
! Find model level for launch
! ----------------------------------------------------------------------+-------
! Geographically fixed value - rethink if variable height sources are used
! WARNING: Not defined relative to tdims%k_start as endgame alters it
launchlev=tkfix1start
do k=(tkfix1start + 1),tdims%k_end
  if (eta_theta_levels(k) >= etalaunch .and.                                   &
      eta_theta_levels(k-1) <  etalaunch) then
    if ((eta_theta_levels(k)-etalaunch) <                                      &
        (etalaunch-eta_theta_levels(k-1))) then
      launchlev=k
    else
      launchlev=k-1
    end if
  end if
end do
!
! ----------------------------------------------------------------------+-------
!     Interpolate : [Vertical]   RHO onto T grid
!             and : [Horizontal] U,V onto P grid
! ----------------------------------------------------------------------+-------
!$OMP PARALLEL DEFAULT(SHARED)                                                 &
!$OMP          private( i, j, k, jdir)
call p_to_t (row_length, rows, pdims_l%halo_i, pdims_l%halo_j,                 &
             pdims_s%halo_i, pdims_s%halo_j,levels-1,r_theta_levels,           &
             r_rho_levels,rho,rhont_smallhalo)
! Need to barrier here to make sure previous routine has fully completed
! before we use the results
!$OMP BARRIER
! Extrapolate topmost level to be a scale height from level below
!$OMP do SCHEDULE(STATIC)
Rows_do_init: do j=tdims%j_start,tdims%j_end
  Row_length_do_init: do i=tdims%i_start,tdims%i_end
    rhont_smallhalo(i,j,tdims%k_end) =                                         &
      rhont_smallhalo(i,j,tdims%k_end-1) *                                     &
      exp(-(r_theta_levels(i,j,tdims%k_end)-                                   &
            r_theta_levels(i,j,tdims%k_end-1)) * rscale_h)
  end do  Row_length_do_init
end do  Rows_do_init
!$OMP end do NOWAIT
!
call u_to_p(u,                                                                 &
                  udims_s%i_start,udims_s%i_end,                               &
                  udims_s%j_start,udims_s%j_end,                               &
                  pdims%i_start,pdims%i_end,                                   &
                  pdims%j_start,pdims%j_end,                                   &
                  levels, uonp )

call v_to_p(v,                                                                 &
                  vdims_s%i_start,vdims_s%i_end,                               &
                  vdims_s%j_start,vdims_s%j_end,                               &
                  pdims%i_start,pdims%i_end,                                   &
                  pdims%j_start,pdims%j_end,                                   &
                  levels, vonp )


! Need to make sure all parallelised calculations in these routines are
! complete before we use the computed values
!$OMP BARRIER

!
! ----------------------------------------------------------------------+-------
!     Initialize local arrays to zero
! ----------------------------------------------------------------------+-------
! ---------------------------------
!     Set winds with halos to zero
! ---------------------------------

!$OMP do SCHEDULE(STATIC)
Levels_do1: do k=pdims_s%k_start,pdims_s%k_end
  Rows_do1: do j=pdims_s%j_start,pdims_s%j_end
    Row_length_do1: do i=pdims_s%i_start,pdims_s%i_end
      g_xp_smallhalo(i,j,k) = 0.0
      g_yp_smallhalo(i,j,k) = 0.0
    end do  Row_length_do1
  end do  Rows_do1
end do  Levels_do1
!$OMP end do NOWAIT
!
!$OMP do SCHEDULE(STATIC)
Levels_do1a: do k=pdims%k_start,pdims%k_end
  Rows_do1a: do j=pdims%j_start,pdims%j_end
    Row_length_do1a: do i=pdims%i_start,pdims%i_end
      ! ----------------------------------------------------------------------+-------
      !     Zero vertical divergence of pseudomomentum flux
      ! ----------------------------------------------------------------------+-------
      g_x(i,j,k) = 0.0
      g_y(i,j,k) = 0.0
    end do  Row_length_do1a
  end do  Rows_do1a
end do  Levels_do1a
!$OMP end do NOWAIT

!
! ----------------------------------------------------------------------+-------
! 1.0   Set variables that are to be defined on all model levels
! ----------------------------------------------------------------------+-------
!

!$OMP do SCHEDULE(STATIC)
Levels_do2: do k=(tkfix1start + 1),(tdims%k_end - 1)
  ! ----------------------------------------------------------------------+-------
  ! 1.1 Density, buoyancy frequency and altitude for middle levels
  ! ----------------------------------------------------------------------+-------
  Rows_do2: do j=tdims%j_start,tdims%j_end
    Row_length_do2: do i=tdims%i_start,tdims%i_end
      rho_th(i,j,k) = rhont_smallhalo(i,j,k) * recip_a2
      !     Buoyancy (Brunt-Vaisala) frequency calculation
      nbv(i,j,k) = ( g*(theta(i,j,k+1) - theta(i,j,k-1)) /                     &
                          (theta(i,j,k) *                                      &
        (r_theta_levels(i,j,k+1) - r_theta_levels(i,j,k-1)) ) )
      nbv(i,j,k) = max(nbv(i,j,k), sqnmin)
      nbv(i,j,k) = sqrt(nbv(i,j,k))
    end do  Row_length_do2
  end do  Rows_do2
end do  Levels_do2
!$OMP end do

!
! ----------------------------------------------------------------------+-------
!     Density, Buoyancy (Brunt-Vaisala) frequency at top and bottom
! ----------------------------------------------------------------------+-------

!$OMP do SCHEDULE(STATIC)
Row_length_do3: do i=tdims%i_start,tdims%i_end
  j = 1
  !   Set rho and theta level value of Z_TH.
  rho_th(i,j,tkfix1start) = rhont_smallhalo(i,j,tkfix1start) *               &
                             recip_a2
  nbv(i,j,tkfix1start)    = nbv(i,j,(tkfix1start + 1))
  rho_th(i,j,tdims%k_end) = rhont_smallhalo(i,j,tdims%k_end) *               &
                             recip_a2
  nbv(i,j,tdims%k_end)    = nbv(i,j,tdims%k_end-1)
end do  Row_length_do3
!$OMP end do


!
Levels_do4: do k=tdims%k_end,(tkfix1start + 1),-1
  ! ----------------------------------------------------------------------+-------
  ! 1.2 Set buoyancy frequency constant up to 1km altitude
  ! ----------------------------------------------------------------------+-------
  !$OMP do SCHEDULE(STATIC)
  Row_length_do4: do i=tdims%i_start,tdims%i_end
    j = 1
    if ( (r_theta_levels(i,j,k) - planet_radius) <  1.0e3)                   &
      nbv(i,j,k-1) = nbv(i,j,k)
  end do  Row_length_do4
  !$OMP end do
end do  Levels_do4

!
! ----------------------------------------------------------------------+-------
! 1.3 Compute component of wind U in each wave-propagation direction.
!     U is the dot product of (u,v) with k_0 but n.b. UDOTK is half-way
!     between rho levels.
!     Interpolate : [Vertical]   RHO onto T grid
! ----------------------------------------------------------------------+-------

!$OMP do SCHEDULE(STATIC)
IDir_do1: do jdir=1,idir
  !
  Levels_do5: do k=tkfix1start,(tdims%k_end - 1)
    Rows_do5: do j=tdims%j_start,tdims%j_end
      Row_length_do5: do i=tdims%i_start,tdims%i_end
        udotk(i,j,k,jdir) =                                                    &
        0.5*(uonp(i,j,k) + uonp(i,j,k+1))*cosphi(jdir) +                       &
        0.5*(vonp(i,j,k) + vonp(i,j,k+1))*sinphi(jdir)
        !       Assumes theta levels are half way between rho levels.
      end do  Row_length_do5
    end do  Rows_do5
  end do  Levels_do5
  !
  ! ----------------------------------------------------------------------+-------
  ! Set wind component for top level, to be equal to that on the top Rho level
  ! ----------------------------------------------------------------------+-------
  Row_length_do5a: do i=tdims%i_start,tdims%i_end
    j = 1
    udotk(i,j,tdims%k_end,jdir) = uonp(i,j,pdims%k_end)*cosphi(jdir)         &
                                + vonp(i,j,pdims%k_end)*sinphi(jdir)
  end do  Row_length_do5a
end do  IDir_do1
!$OMP end do

!$OMP end PARALLEL

!
! ----------------------------------------------------------------------+-------
! Core calculations on Physics Grid now within separate subroutine
! ----------------------------------------------------------------------+-------
!
! As presently coded, a single launch level is set but in principle it could be
! variable. Optimization is easier if we define level limits such that:
! 1 <= k < minlaunchlev            : initial total flux = 0.
! minlaunchlev <= k < maxlaunchlev : initial total flux = 0. / launch flux(i,j)
! maxlaunchlev <= k                : initial total flux = launch flux(i,j)
! For now they are simply assigned (also replicated in the core code)
minlaunchlev = launchlev
maxlaunchlev = launchlev
! ----------------------------------------------------------------------+-------
! 3.0 Initialize gravity wave spectrum variables
! ----------------------------------------------------------------------+-------
! 4.0 Calculations carried out for each azimuthal direction
! ----------------------------------------------------------------------+-------
!
! L_add_cgw = F : calculate standard USSP isotropic GW launch flux
! L_add_cgw = T : calculate variable CGW launch flux

! Number of threads is 1 at this point: we are not load-balancing between
! threads, then segmenting here.  Rather, set up segmentation on one thread,
! then distribute the segments among threads.
! Variables needed for segmentation
num_parallel = 1
ipar         = 1
num_segments = -99

!Set up segment meta-information.
call segments_mod_seg_meta(meta_segments, ipar, num_parallel,                  &
                           rows*row_length, ussp_seg_size, num_segments)

! Allocate space for segmentation arrays
allocate(  segments( meta_segments%num_segments ) )

! Work out starting points and sizes of each segment individually
call segments_mod_segments(segments, meta_segments, row_length)
! Main parallelisation and segmentation loop
! For this we allocate space for each segment and then copy necessary data
! in to these variables before calling the ussp_core with them. The result
! is then copied back out. Looping is a little awkward as segments may well
! span partial rows.
!$OMP PARALLEL DEFAULT(SHARED)                                                 &
!$OMP private(ii, jj, ss, i, k, jdir, s_sin_theta_lat, s_totalppn,             &
!$OMP         s_rho_th, s_nbv, s_udotk, s_fptot)

!$OMP do SCHEDULE(DYNAMIC)
do i = 1, meta_segments%num_segments

  allocate(s_sin_theta_lat(segments(i)%seg_points))
  allocate(s_totalppn(segments(i)%seg_points))
  allocate(s_rho_th(segments(i)%seg_points,tkfix1start:tdims%k_end))
  allocate(s_nbv(segments(i)%seg_points,tkfix1start:tdims%k_end))
  allocate(s_udotk(segments(i)%seg_points, tkfix1start:tdims%k_end, idir))
  allocate(s_fptot(segments(i)%seg_points, tkfix1start:tdims%k_end, idir))

  ii = segments(i)%first_x
  jj = segments(i)%first_y
  ss = 1
  do while (ss <= segments(i)%seg_points)
    s_sin_theta_lat(ss) = real(sin_theta_latitude(ii,jj), kind=real_usprec)
    s_totalppn(ss)      = real(totalppn(ii,jj), kind=real_usprec)
    ss = ss + 1
    if (ii == tdims%i_end) then
      ii = tdims%i_start
      jj = jj + 1
    else
      ii = ii + 1
    end if
  end do

  do k = tkfix1start, tdims%k_end
    ii = segments(i)%first_x
    jj = segments(i)%first_y
    ss = 1
    do while (ss <= segments(i)%seg_points)
      s_nbv(ss,k) = real(nbv(ii,jj,k), kind=real_usprec)
      s_rho_th(ss,k) = real(rho_th(ii,jj,k), kind=real_usprec)
      ss = ss + 1
      if (ii == tdims%i_end) then
        ii = tdims%i_start
        jj = jj + 1
      else
        ii = ii + 1
      end if
    end do
  end do

  do jdir = 1, idir
    do k = tkfix1start, tdims%k_end
      ii = segments(i)%first_x
      jj = segments(i)%first_y
      ss = 1
      do while (ss <= segments(i)%seg_points)
        s_udotk(ss,k,jdir)    = real(udotk(ii,jj,k,jdir), kind=real_usprec)
        ss = ss + 1
        if (ii == tdims%i_end) then
          ii = tdims%i_start
          jj = jj + 1
        else
          ii = ii + 1
        end if
      end do
    end do
  end do

  call gw_ussp_core(segments(i)%seg_points, launchlev,                         &
         tkfix1start, tdims%k_end, ussp_launch_factor, wavelstar,              &
         cgw_scale_factor, two_omega, s_sin_theta_lat, s_rho_th, s_nbv,        &
         s_udotk, s_totalppn, s_fptot, L_add_cgw)

  do jdir = 1, idir
    do k = tkfix1start, tdims%k_end
      ii = segments(i)%first_x
      jj = segments(i)%first_y
      ss = 1
      do while (ss <= segments(i)%seg_points)
        fptot(ii,jj,k,jdir) = real(s_fptot(ss,k,jdir), kind=real_umphys)
        ss = ss + 1
        if (ii == tdims%i_end) then
          ii = tdims%i_start
          jj = jj + 1
        else
          ii = ii + 1
        end if
      end do
    end do
  end do

  deallocate(s_fptot)
  deallocate(s_udotk)
  deallocate(s_nbv)
  deallocate(s_rho_th)
  deallocate(s_totalppn)
  deallocate(s_sin_theta_lat)
end do
!$OMP end do

!$OMP end PARALLEL
deallocate(segments)

! ----------------------------------------------------------------------+-------
!

!$OMP  PARALLEL DEFAULT(none) SHARED(fptot, g_x, g_y, cosphi,                  &
!$OMP  p_layer_boundaries, g_xp_smallhalo, g_yp_smallhalo, timestep,           &
!$OMP  sinphi, pdims, minlaunchlev, L_ussp_heating, r_theta_levels,            &
!$OMP  r_rho_levels, uonp, vonp, T_inc, tdims, g, cp)                          &
!$OMP  private(g_g, i, j, k, jdir, jj, dzb, dzu, dzl, uhat, vhat,              &
!$OMP  ududt, vdvdt )

!
! ----------------------------------------------------------------------+-------
! 5.0   Compute vertical divergence of pseudomomentum flux.
!       Column Integral Converts : [Vertical]   T onto RHO grid
! ----------------------------------------------------------------------+-------

! gives each thread the largest block possible to execute


!$OMP do SCHEDULE(STATIC)
Levels_do14: do k=minlaunchlev+1, pdims%k_end
    Rows_do14: do j=pdims%j_start,pdims%j_end
      Row_length_do14: do i=pdims%i_start,pdims%i_end
        IDir_do5: do jdir=1,idir
          !       Pseudomomentum flux
          g_g = g * (fptot(i,j,k,jdir) - fptot(i,j,k-1,jdir)) /                &
           (p_layer_boundaries(i,j,k) - p_layer_boundaries(i,j,k-1))
          g_x(i,j,k) = g_x(i,j,k) + g_g* cosphi(jdir)
          g_y(i,j,k) = g_y(i,j,k) + g_g* sinphi(jdir)
        end do  IDir_do5
      end do  Row_length_do14
    end do  Rows_do14
end do  Levels_do14
!$OMP end do

!
! ----------------------------------------------------------------------+-------
! 5.1   Wind and temperature increments from wave dissipation
! ----------------------------------------------------------------------+-------

!$OMP do SCHEDULE(STATIC)
Levels_do15: do k=minlaunchlev,pdims%k_end
  Rows_do15: do j=pdims%j_start,pdims%j_end
    Row_length_do15: do i=pdims%i_start,pdims%i_end
      !     g_xp_smallhalo(i,j,k) = g_x(i,j,k)
      !     g_yp_smallhalo(i,j,k) = g_y(i,j,k)
      !! ABOVE is preferred but following should reduce bit comparison differences
      !! for purposes of initial testing. REMOVE when scheme goes LIVE!
      g_xp_smallhalo(i,j,k) = timestep * g_x(i,j,k)
      g_yp_smallhalo(i,j,k) = timestep * g_y(i,j,k)
    end do  Row_length_do15
  end do  Rows_do15
end do  Levels_do15
!$OMP end do

!-----------------------------------------------------------------------+-------
! Calculate heating due to gravity wave dissipation
!-----------------------------------------------------------------------+-------
GW_heating: if ( L_ussp_heating ) then

!$OMP do SCHEDULE(STATIC)
  Levels_do16: do k = tkfix1start,tdims%k_end-1
    Rows_do16: do j = tdims%j_start,tdims%j_end
      Row_length_do16: do i = tdims%i_start,tdims%i_end

        dzb = r_theta_levels(i,j,k) -  r_rho_levels(i,j,k)
        dzu = r_rho_levels(i,j,k+1) -  r_theta_levels(i,j,k)
        dzl = r_rho_levels(i,j,k+1) -  r_rho_levels(i,j,k)

        !       u and v on theta_level(k)
        uhat   = dzu * uonp(i,j,k) + dzb * uonp(i,j,k+1)
        vhat   = dzu * vonp(i,j,k) + dzb * vonp(i,j,k+1)

        !       u*du and v*dv on theta_level(k)
        ududt  = uhat *( dzu * g_x(i,j,k) + dzb * g_x(i,j,k+1) )
        vdvdt  = vhat *( dzu * g_y(i,j,k) + dzb * g_y(i,j,k+1) )

        !       dT/dt on theta_level(k)
        T_inc(i,j,k)=T_inc(i,j,k) - timestep * ( ududt + vdvdt ) /             &
                    ( cp * dzl * dzl )

      end do  Row_length_do16
    end do  Rows_do16
  end do  Levels_do16
!$OMP end do

end if GW_heating

!$OMP end PARALLEL

!
! ----------------------------------------------------------------------+-------
! Put U,V increments onto U,V grid
! Interpolate : [Horizontal] P onto U,V grid
! ----------------------------------------------------------------------+-------

!-------------------------------------------------------------------------------
! Block of code below that performs swap_bounds, interpolation to u or v points
! and computes diagnostics is not compiled in LFRic
! An #else statement below sets r_u and r_v increments without swap_bounds
! or interpolation when compiling in LFRic
!-------------------------------------------------------------------------------
do k=minlaunchlev,udims%k_end
  do j=udims%j_start,udims%j_end
    do i=udims%i_start,udims%i_end
      r_u(i,j,k) = r_u(i,j,k) + (g_xp_smallhalo(i,j,k))
    end do
  end do
  do j=vdims%j_start,vdims%j_end
    do i=vdims%i_start,vdims%i_end
      r_v(i,j,k) = r_v(i,j,k) + (g_yp_smallhalo(i,j,k))
    end do
  end do
end do
! ----------------------------------------------------------------------+-------
! Diagnostic output : Northward Flux of Horiz Pseudomom (on theta levels)
! ----------------------------------------------------------------------+-------
if (gwspec_nflux_on) then
  do k=tkfix1start,tdims%k_end
    do j=tdims%j_start,tdims%j_end
      do i=tdims%i_start,tdims%i_end
        gwspec_nflux(i,j,k) = fptot(i,j,k,1)
      end do
    end do
  end do
end if
! ----------------------------------------------------------------------+-------
! Diagnostic output : Westward Flux of Horiz Pseudomom (on theta levels)
! ----------------------------------------------------------------------+-------
if (gwspec_wflux_on .or. gwspec_wflux_p_on) then
  do k=tkfix1start,tdims%k_end
    do j=tdims%j_start,tdims%j_end
      do i=tdims%i_start,tdims%i_end
        gwspec_wflux(i,j,k) = fptot(i,j,k,2)
      end do
    end do
  end do
end if
! ----------------------------------------------------------------------+-------
! Diagnostic output : Southward Flux of Horiz Pseudomom (on theta levels)
! ----------------------------------------------------------------------+-------
if (gwspec_sflux_on) then
  do k=tkfix1start,tdims%k_end
    do j=tdims%j_start,tdims%j_end
      do i=tdims%i_start,tdims%i_end
        gwspec_sflux(i,j,k) = fptot(i,j,k,3)
      end do
    end do
  end do
end if
! ----------------------------------------------------------------------+-------
! Diagnostic output : Eastward Flux of Horiz Pseudomom (on theta levels)
! ----------------------------------------------------------------------+-------
if (gwspec_eflux_on .or. gwspec_eflux_p_on) then
  do k=tkfix1start,tdims%k_end
    do j=tdims%j_start,tdims%j_end
      do i=tdims%i_start,tdims%i_end
        gwspec_eflux(i,j,k) = fptot(i,j,k,4)
      end do
    end do
  end do
end if

!
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
!
end subroutine gw_ussp

end module gw_ussp_mod
