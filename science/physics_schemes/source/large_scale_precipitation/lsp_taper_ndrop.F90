! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Cloud droplet number calculator
! Subroutine Interface:
module lsp_taper_ndrop_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='LSP_TAPER_NDROP_MOD'

contains


subroutine lsp_taper_ndrop(                                                    &
                   ! (Full) Aerosol tracers
                            so4_acc, so4_dis, sea_salt_film,                   &
                            biogenic, sea_salt_jet, bmass_agd,                 &
                            bmass_cld, ocff_agd, ocff_cld,                     &
                            nitr_acc, nitr_diss,                               &
                            n_arcl_compnts,                                    &
                            i_arcl_compnts, arcl,                              &
                   ! Murk aerosol
                            aerosol,                                           &
                   ! CDNC from UKCA
                            ukca_cdnc,                                         &
                            cdnc_dim1, cdnc_dim2, cdnc_dim3,                   &
                   ! CDNC from EasyAerosol
                            easyaerosol_cdnc,                                  &
                   ! Other parameters
                            rhodz_dry, rhodz_moist,                            &
                            deltaz,                                            &
                            snow_depth, land_fract,                            &
                   ! Output parameter of n_drop_tpr
                            n_drop_tpr                                         &
                                 )
! Microphysics modules

use mphys_inputs_mod,      only: ndrop_surf,                                   &
                                 l_autoconv_murk, l_droplet_tpr,               &
                                 l_mcr_arcl, arcl_inhom_sc,                    &
                                 l_use_sulphate_autoconv,                      &
                                 l_use_bmass_autoconv,                         &
                                 l_use_ocff_autoconv,                          &
                                 l_use_nitrate_autoconv,                       &
                                 l_use_seasalt_autoconv

use mphys_constants_mod,   only: ntot_land, ntot_sea, max_drop,                &
                                 n0_murk, m0_murk,                             &
                                 mu_g_to_kg, min_n_aer

use lsp_autoc_consts_mod,  only: power_murk, min_drop_alt, eta_peak,           &
                                 eta_low_nd, level_peak, level_surf,           &
                                 vala_fac1, vala_fac2,                         &
                                 half_range

! Dynamics and grid bounds modules

use level_heights_mod,     only: eta_theta_levels
use atm_fields_bounds_mod, only: tdims

! General constants modules
use arcl_mod,              only: npd_arcl_compnts,                             &
                                 ip_arcl_sulp_ac,  ip_arcl_sulp_di,            &
                                 ip_arcl_sslt_fi,  ip_arcl_sslt_jt,            &
                                 ip_arcl_biom_ag, ip_arcl_biom_ic,             &
                                 ip_arcl_ocff_ag, ip_arcl_ocff_ic

use conversions_mod,       only: pi
use rad_input_mod,         only: l_use_biogenic, l_use_arclbiom,               &
                                 l_use_arclsslt, l_use_arclocff
use gen_phys_inputs_mod,   only: l_mr_physics
use ukca_option_mod,       only: l_ukca_aie2
use glomap_clim_option_mod, only: l_glomap_clim_aie2
use murk_inputs_mod,       only: l_murk
use number_droplet_mod,    only: number_droplet, min_cdnc_sea_ice,             &
                                 min_cdnc_land, snow_depth_thresh,             &
                                 land_frac_thresh

! Stochastic physics module
use stochastic_physics_run_mod, only: l_rp2, i_rp_scheme, i_rp2b,              &
                                      ndrop_surf_rp, rp_idx

! EasyAerosol
use easyaerosol_option_mod,only: l_easyaerosol_autoconv
use def_easyaerosol, only: t_easyaerosol_cdnc

! Dr Hook modules
!--------------------------------------
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
!--------------------------------------

implicit none

!---------------------------------------------------------------------
! Purpose:
!   Calculates cloud drop number concentration for the microphysics
!   scheme.

!   By using the tapering method, we can also reduce droplet number
!   in the atmospheric boundary layer, from a peak at a given altitude
!   (z_peak_nd) to a user defined value (ndrop_surf) or a variable
!   value dependent on aerosol amounts.

! Method:
!  1) Determine the height (specifically eta value) below which droplet
!     number is to taper. When the taper curve is inactive, this value
!     is not used.

!  2) Above this height, calculate droplet number using the
!     Clark-Jones formulae for Murk aerosol
!     (as in the autoconversion routine).

!     Clark et al (2008): (aerosol number) :
!              n_aer = n0_murk * ( Aerosol / m0_murk ) ** power_murk

!     Jones  et al (1994): Aerosol number to droplet number:
!              n_d   = 3.75e8 * ( 1 - exp ( -2.5e-9 * n_aer ) )

!     If full prognostic aerosols are used or aerosol climatologies
!     are used, then the routine calls the number_droplet function
!     to generate cloud droplet number concentration.

!     If no aerosol species are available and droplet tapering is
!     requested, then the routine uses a simple profile of
!     pre-determined values for droplet number.

!  3) Using the first height above z_peak_nd, determine the droplet
!     profile at each level towards the surface

!    n_drop = ndrop_surf + vala * log ( z / z_surf )

!    where vala = ( nd(z_peak_nd) - ndrop_surf ) / log ( z_peak_nd / z_surf)

!    z_surf is the height at which the cloud drop number reaches ndrop_surf
!    and it is kept constant below this

!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.

!------------------------------------------------
! Subroutine arguments
!------------------------------------------------

integer :: cdnc_dim1, cdnc_dim2, cdnc_dim3

integer, intent(in) ::                                                         &
  n_arcl_compnts,                                                              &
                    ! Number of requested arcl components
  i_arcl_compnts(npd_arcl_compnts)
                    ! Array index of each aerosol clim component

real(kind=real_umphys), intent(in) :: ukca_cdnc(cdnc_dim1, cdnc_dim2, cdnc_dim3)
! CDNC from UKCA

type (t_easyaerosol_cdnc), intent(in) :: easyaerosol_cdnc
! CDNC from EasyAerosol

real(kind=real_umphys), intent(in) ::                                          &
!-----------
! Aerosols
!-----------
  so4_acc      ( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end,                                  &
                             1 : tdims%k_end ),                                &
! SO4
! acc
  so4_dis      ( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end,                                  &
                             1 : tdims%k_end ),                                &
! SO4
! dis
  sea_salt_film( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end,                                  &
                             1 : tdims%k_end ),                                &
! Sea salt
! film
  biogenic     ( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end,                                  &
                             1 : tdims%k_end ),                                &
! Biogenic
! aerosol
  sea_salt_jet ( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end,                                  &
                             1 : tdims%k_end ),                                &
! Sea salt
! Jet
  bmass_agd    ( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end,                                  &
                             1 : tdims%k_end ),                                &
! Aged
! Biomass
  bmass_cld    ( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end,                                  &
                             1 : tdims%k_end ),                                &
! Cloudy
! Biomass
  ocff_agd     ( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end,                                  &
                             1 : tdims%k_end ),                                &
! ocff
! aged
  ocff_cld     ( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end,                                  &
                             1 : tdims%k_end ),                                &
! ocff
! cloud
  nitr_acc     ( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end,                                  &
                             1 : tdims%k_end ),                                &
! Nitr
! aged
  nitr_diss    ( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end,                                  &
                             1 : tdims%k_end ),                                &
!nitrate
!diss
  arcl         ( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end,                                  &
                             1 : tdims%k_end,                                  &
                             1 : n_arcl_compnts ),                             &
! Aerosol climatologies for droplet number calculations

  rhodz_dry    ( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end,                                  &
                             1 : tdims%k_end ),                                &
  rhodz_moist  ( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end,                                  &
                             1 : tdims%k_end ),                                &
  deltaz       ( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end,                                  &
                             1 : tdims%k_end ),                                &
!Grid box information

  land_fract   ( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end ),                                &
! Land fraction
!(for Jones et al,
! 1994 method)
  snow_depth   ( tdims%i_start : tdims%i_end ,                                 &
                 tdims%j_start : tdims%j_end )
! Depth of snow
!(for Jones et al,
! 1994 method)

real(kind=real_umphys), intent(in out) ::                                      &
  aerosol      ( tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end,                                  &
                             1 : tdims%k_end )
! Murk aerosol for
! droplet number
! calculations

real(kind=real_umphys), intent(out) ::                                         &
! intended output
    n_drop_tpr( tdims%i_start : tdims%i_end,                                   &
                tdims%j_start : tdims%j_end, 1:tdims%k_end )
! Droplet number after tapering has taken place

!------------------------------------------------
! Local variables
!------------------------------------------------

real(kind=real_umphys) ::                                                      &
      n_aer,                                                                   &
         ! Aerosol number from
         ! Clark et al (2008)
      vala,n_aer2,                                                             &
         ! Droplet number determined constants
      ndrop_surf2( tdims%i_start : tdims%i_end,                                &
                   tdims%j_start : tdims%j_end )
         ! Variable surface droplet number

! 3D version of Air density in kg/m3
real(kind=real_umphys) ::                                                      &
        rho3d(  tdims%i_start : tdims%i_end,                                   &
                tdims%j_start : tdims%j_end,                                   &
                1 : tdims%k_end )

integer :: i, j, k              ! Loop counters
integer :: vec_length

! Logical to set nitrate climatology. Currently hardwired to .false.
! as a nitrate climatology is not yet available.
logical, parameter :: l_use_arclnitr = .false.

! Is the input "sulphate" aerosol in the form of
! ammonium sulphate (T) or just sulphur (F)? Used in the
! same way for the nitrate aerosol, in form of ammonium
! nitrate (T) or just nitrogen (F).
logical, parameter :: l_nh42so4 = .false.

real(kind=real_umphys) :: temp1, temp2
! temp for plane invariant reciprocals
real(kind=real_umphys) ::  tempv
! temp  variable for intermediate computations

!------------------------------------------------
! Dr Hook subroutine timer details
!------------------------------------------------
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_TAPER_NDROP'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
!---------------------------------------------------------------------
! Start of physics
!---------------------------------------------------------------------

! First update value of ndrop_surf if random parameters 2b
! scheme is switched on
if ( l_rp2 .and. i_rp_scheme == i_rp2b ) then
  ndrop_surf = ndrop_surf_rp(rp_idx)
end if

vec_length  = tdims%i_len*tdims%j_len

if (l_droplet_tpr) then

  !------------------------------------------------
  ! Droplet tapering is on, so we need to calculate
  ! droplet number concentration and add in a taper
  !------------------------------------------------

  if (l_murk .and. l_autoconv_murk) then

    !----------------------------------------------
    ! Use murk aerosol to calculate droplet number
    ! and taper this profile
    !----------------------------------------------

    j = 1 ! Arrays passed in by LFRic have a size of 1 in the j dimension.
!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC) private(i)                    &
!$OMP SHARED(tdims,ndrop_surf2,ndrop_surf,j)
    do i = tdims%i_start, tdims%i_end
      ndrop_surf2(i,j) = ndrop_surf
    end do ! tdims%i
!$OMP end PARALLEL do

!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC) private(i,j,k,n_aer2,tempv)   &
!$OMP SHARED(tdims,level_peak,aerosol,land_fract,snow_depth,n_drop_tpr,        &
!$OMP  n0_murk,m0_murk)
    do k = level_peak,tdims%k_end

      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end

          ! Above the taper level, so use Clark-Jones formulae:

          n_aer2 = max( aerosol(i,j, k) / m0_murk * mu_g_to_kg, min_n_aer)
          ! 1.0E-9 converts from ug/kg to kg/kg

          !-----------------------------------------------
          ! Calculation of the aerosol number
          !-----------------------------------------------

          tempv=n_aer2**power_murk

          !-----------------------------------------------
          ! Convert to CCN using the Jones et al (1994)
          ! relationship (follows number_droplet routine
          ! but reproduced here to avoid compiler issues).
          !-----------------------------------------------

          n_aer2 = -2.5e-9 *  n0_murk  * tempv

          n_drop_tpr(i,j,k) = max_drop * ( 1.0e+00 - exp(n_aer2)  )

        end do
      end do

      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end

          if ( land_fract(i,j) >= land_frac_thresh .and.                       &
               snow_depth(i,j) < snow_depth_thresh    ) then

            n_drop_tpr(i,j,k) = max( n_drop_tpr(i,j,k), min_cdnc_land )

          else

            n_drop_tpr(i,j,k) = max( n_drop_tpr(i,j,k), min_cdnc_sea_ice  )

          end if ! land fract > 0.2 etc

        end do
      end do
    end do
!$OMP end PARALLEL do

    ! Start tapering

    do k=1, level_peak-1

      temp1 =  log( eta_theta_levels(k) / eta_theta_levels(1) )
      temp2 =  log( eta_theta_levels(k) / eta_theta_levels(level_surf) )

      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end

          vala = (n_drop_tpr(i,j,level_peak) - ndrop_surf2(i,j))

          if (n_drop_tpr(i,j,level_peak) > ndrop_surf2(i,j)) then

            ! if drop number is increasing with height, keep it at
            ! its surface value below z_surf before tapering up to
            ! the z_peak value

            n_drop_tpr( i, j, k ) =  max( ndrop_surf2(i,j),                    &
                 ndrop_surf2(i,j) + ( vala * vala_fac2 ) * temp2 )

          else

            ! if drop number is decreasing with height, make sure it
            ! only reaches its surface value at level 1

            n_drop_tpr( i, j, k ) = ndrop_surf2(i,j) +                         &
                                    ( vala * vala_fac1 ) * temp1

          end if

        end do ! tdims%i
      end do   ! tdims%j

    end do ! eta values below peak

  else if (l_mcr_arcl) then

    !------------------------------------------------------------
    ! Full droplet number calculation with aerosol climatologies
    !------------------------------------------------------------
    ! Step 1: Calculate air density, rho
!$OMP  PARALLEL DEFAULT(none)                                                  &
!$OMP  SHARED( level_peak, tdims, l_mr_physics, rho3d, rhodz_dry,              &
!$OMP          deltaz, rhodz_moist )                                           &
!$OMP  private( i, j, k )
    if (l_mr_physics) then
!$OMP do SCHEDULE(STATIC)
      do k = 1, tdims%k_end
        do j = tdims%j_start, tdims%j_end
          do i = tdims%i_start, tdims%i_end
            ! rho is the dry density
            rho3d(i,j,k) = rhodz_dry(i,j,k) / deltaz(i,j,k)
          end do
        end do
      end do
!$OMP end do

    else ! l_mr_physics
!$OMP do SCHEDULE(STATIC)
      do k = 1, tdims%k_end
        do j = tdims%j_start, tdims%j_end
          do i = tdims%i_start, tdims%i_end
            ! rho is the moist density
            rho3d(i,j,k) = rhodz_moist(i,j,k) / deltaz(i,j,k)
          end do
        end do
      end do
!$OMP end do
    end if  ! l_mr_physics
!$OMP end PARALLEL

    ! Step 2: Call number_droplet routine
    call number_droplet(                                                       &
                         tdims%i_start, tdims%i_end,                           &
                         tdims%j_start, tdims%j_end,                           &
                         1, tdims%k_end,                                       &
                         1, tdims%k_end,                                       &
                         l_mcr_arcl,                                           &
                         l_nh42so4,                                            &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_sulp_ac)),          &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_sulp_di)),          &
                         l_use_arclsslt,                                       &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_sslt_fi)),          &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_sslt_jt)),          &
                         l_use_biogenic, biogenic,                             &
                         l_use_arclbiom,                                       &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_biom_ag)),          &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_biom_ic)),          &
                         l_use_arclocff,                                       &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_ocff_ag)),          &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_ocff_ic)),          &
                         l_use_arclnitr,                                       &
                         rho3d,                                                &
                         snow_depth, land_fract,                               &
                         n_drop_tpr, aerosol=aerosol)

!$OMP  PARALLEL DEFAULT(none)                                                  &
!$OMP  SHARED( level_peak, tdims, n_drop_tpr, arcl_inhom_sc,                   &
!$OMP          ndrop_surf, eta_theta_levels, level_surf, vala_fac1,            &
!$OMP          vala_fac2 )                                                     &
!$OMP  private( i, j, k, temp1, temp2, vala )
!$OMP  do SCHEDULE(STATIC)
    do k = level_peak, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          ! scaling the n_drop by the inhomogeneity scaling:
          n_drop_tpr(i,j,k) = n_drop_tpr(i,j,k) * arcl_inhom_sc

        end do
      end do
    end do
!$OMP end do

!$OMP do SCHEDULE(STATIC)
    do k = 1, level_peak-1
      temp1 =  log( eta_theta_levels(k) / eta_theta_levels(1) )
      temp2 =  log( eta_theta_levels(k) / eta_theta_levels(level_surf) )
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end

          vala = (n_drop_tpr(i,j,level_peak) - ndrop_surf)

          if (n_drop_tpr(i,j,level_peak) > ndrop_surf) then

            ! if drop number is increasing with height, keep it at
            ! its surface value below z_surf before tapering up to
            ! the z_peak value

            n_drop_tpr( i, j, k ) =  max( ndrop_surf,                          &
                 ndrop_surf + ( vala * vala_fac2 ) * temp2 )

          else

            ! if drop number is decreasing with height, make sure it
            ! only reaches its surface value at level 1

            n_drop_tpr( i, j, k ) =  ndrop_surf + ( vala * vala_fac1 ) * temp1

          end if
        end do ! tdims%i
      end do ! tdims%j
    end do ! tdims%k
!$OMP end do
!$OMP end PARALLEL

  else if (l_use_sulphate_autoconv .or. l_ukca_aie2 .or.                       &
           l_glomap_clim_aie2 .or. l_easyaerosol_autoconv) then

    !------------------------------------------------------------
    ! Full droplet number calculation based on prognostic aerosol
    !------------------------------------------------------------

!$OMP  PARALLEL DEFAULT(none)                                                  &
!$OMP  SHARED( l_ukca_aie2, l_glomap_clim_aie2, l_mr_physics, tdims,           &
!$OMP          n_drop_tpr, rho3d,                                              &
!$OMP          level_peak, ukca_cdnc, rhodz_dry, deltaz, rhodz_moist,          &
!$OMP          l_easyaerosol_autoconv, easyaerosol_cdnc )                      &
!$OMP  private( i, j, k )
    if (l_easyaerosol_autoconv) then
!$OMP do SCHEDULE(STATIC)
      do k = level_peak, tdims%k_end
        do j = tdims%j_start, tdims%j_end
          do i = tdims%i_start, tdims%i_end
            n_drop_tpr(i,j,k) = easyaerosol_cdnc%cdnc(i,j,k)
          end do
        end do
      end do
!$OMP end do
    else
      if (l_ukca_aie2 .or. l_glomap_clim_aie2) then
!$OMP do SCHEDULE(STATIC)
        do k = level_peak, tdims%k_end
          do j = tdims%j_start, tdims%j_end
            do i = tdims%i_start, tdims%i_end
              n_drop_tpr(i,j,k) = ukca_cdnc(i,j,k)
            end do
          end do
        end do
!$OMP end do
      else  ! (not l_ukca_aie2 and not l_glomap_clim_aie2)

            ! Step 1: Calculate air density, rho

        if (l_mr_physics) then
!$OMP do SCHEDULE(STATIC)
          do k = level_peak, tdims%k_end
            do j = tdims%j_start, tdims%j_end
              do i = tdims%i_start, tdims%i_end

                ! rho is the dry density

                rho3d(i,j,k) = rhodz_dry(i,j,k) / deltaz(i,j,k)
              end do
            end do
          end do
!$OMP end do
        else ! l_mr_physics
!$OMP do SCHEDULE(STATIC)
          do k = level_peak, tdims%k_end
            do j = tdims%j_start, tdims%j_end
              do i = tdims%i_start, tdims%i_end

                ! rho is the moist density

                rho3d(i,j,k) = rhodz_moist(i,j,k) / deltaz(i,j,k)
              end do
            end do
          end do
!$OMP end do
        end if  ! l_mr_physics

      end if  ! l_ukca_aie2 .or. l_glomap_clim_aie2

    end if ! l_easyaerosol_autoconv
!$OMP end PARALLEL

    ! Step 2: Call number_droplet routine
    if (.not. l_easyaerosol_autoconv .and. .not.                               &
                                    (l_ukca_aie2 .or. l_glomap_clim_aie2)) then

      call number_droplet(                                                     &
                       tdims%i_start, tdims%i_end,                             &
                       tdims%j_start, tdims%j_end,                             &
                       1, tdims%k_end,                                         &
                       level_peak, tdims%k_end,                                &
                       l_use_sulphate_autoconv,                                &
                       l_nh42so4,                                              &
                       so4_acc,                                                &
                       so4_dis,                                                &
                       l_use_seasalt_autoconv,                                 &
                       sea_salt_film,                                          &
                       sea_salt_jet,                                           &
                       l_use_biogenic, biogenic,                               &
                       l_use_bmass_autoconv,                                   &
                       bmass_agd,                                              &
                       bmass_cld,                                              &
                       l_use_ocff_autoconv,                                    &
                       ocff_agd, ocff_cld,                                     &
                       l_use_nitrate_autoconv,                                 &
                       rho3d,                                                  &
                       snow_depth, land_fract,                                 &
                       n_drop_tpr,                                             &
                       nitr_acc,                                               &
                       nitr_diss)
    end if


!$OMP  PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                              &
!$OMP  SHARED( level_peak, tdims, n_drop_tpr, ndrop_surf, level_surf,          &
!$OMP          eta_theta_levels, vala_fac1, vala_fac2 )                        &
!$OMP  private( i, j, k, vala, temp1, temp2 )
    do k = 1, level_peak-1
      temp1 =  log( eta_theta_levels(k) / eta_theta_levels(1) )
      temp2 =  log( eta_theta_levels(k) / eta_theta_levels(level_surf) )
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end

          vala = (n_drop_tpr(i,j,level_peak) - ndrop_surf)

          if (n_drop_tpr(i,j,level_peak) > ndrop_surf) then

            ! if drop number is increasing with height, keep it at
            ! its surface value below z_surf before tapering up to
            ! the z_peak value

            n_drop_tpr( i, j, k ) =  max( ndrop_surf,                          &
                 ndrop_surf + ( vala * vala_fac2 ) * temp2 )

          else

            ! if drop number is decreasing with height, make sure it
            ! only reaches its surface value at level 1

            n_drop_tpr( i, j, k ) =  ndrop_surf + ( vala * vala_fac1 ) * temp1

          end if
        end do ! tdims%i
      end do ! tdims%j
    end do ! tdims%k
!$OMP end PARALLEL do

  else   ! not l_murk, l_mcr_arcl, l_use_sulphate_autoconv or l_ukca_aie2
         ! or l_glomap_clim_aie2 or l_easyaerosol_autoconv

    !-----------------------------------------
    ! Use a simple profile for tapering
    !-----------------------------------------

    ! Do not allow a peak value of eta to exceed that of the
    ! pre-determined low droplet number (usually around 2km)

    if (eta_peak > eta_low_nd) eta_peak = eta_low_nd

!$OMP  PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                              &
!$OMP  SHARED( tdims, eta_theta_levels, eta_low_nd, eta_peak,                  &
!$OMP          n_drop_tpr, ndrop_surf, vala_fac1, vala_fac2,                   &
!$OMP          level_surf )                                                    &
!$OMP  private( i, j, k, vala, temp1, temp2 )
    do k = 1, tdims%k_end

      if ( eta_theta_levels(k) >= eta_low_nd .and.                             &
           eta_theta_levels(k) >= eta_peak         ) then

        ! Above taper level and level set for minimum droplet
        ! number, hence set n_drop at this level to low value

        do j = tdims%j_start, tdims%j_end
          do i = tdims%i_start, tdims%i_end
            n_drop_tpr(i,j,k) = min_drop_alt
          end do
        end do

      else if ( eta_theta_levels(k) >= eta_peak .and.                          &
                eta_theta_levels(k) < eta_low_nd     ) then

        ! Above taper level yet below level set for minimum
        ! droplet number, so use a simple function of eta to
        ! determine the droplet number. This is intended to
        ! be independent of model levels

        do j = tdims%j_start, tdims%j_end
          do i = tdims%i_start, tdims%i_end
            n_drop_tpr(i,j,k) = half_range * (- cos (pi + pi *                 &
                          ( ( eta_theta_levels (k) - eta_peak ) /              &
                          ( eta_low_nd - eta_peak ) ) ) )                      &
                          + half_range + min_drop_alt
          end do
        end do

        ! No need to set up peak droplet as this should be max_drop

      else ! eta_theta_levels

        vala = (max_drop - ndrop_surf)

        if (max_drop > ndrop_surf) then

          ! if drop number is increasing with height, keep it at
          ! its surface value below z_surf before tapering up to
          ! the z_peak value
          temp2 =  log( eta_theta_levels(k) / eta_theta_levels(level_surf) )

          do j = tdims%j_start, tdims%j_end
            do i = tdims%i_start, tdims%i_end
              n_drop_tpr( i, j, k ) =  max( ndrop_surf,                        &
                   ndrop_surf + ( vala / vala_fac2 ) * temp2 )
            end do ! tdims%i
          end do ! tdims%j

        else

          ! if drop number is decreasing with height, make sure it
          ! only reaches its surface value at level 1
          temp1 =  log( eta_theta_levels(k) / eta_theta_levels(1) )

          do j = tdims%j_start, tdims%j_end
            do i = tdims%i_start, tdims%i_end
              n_drop_tpr( i, j, k ) = ndrop_surf + ( vala * vala_fac1 )       &
                                      * temp1
            end do ! tdims%i
          end do ! tdims%j

        end if

      end if !eta_levels above or below peak values

    end do  ! tdims%k
!$OMP end PARALLEL do

  end if ! l_murk / l_use_sulphate_autoconv

else ! l_droplet_tpr

  !-------------------------------------------------
  ! Droplet tapering is not active, but we still
  ! need to calculate potential cloud drop number
  ! concentration for use in autoconversion
  !-------------------------------------------------

  if (l_murk .and. l_autoconv_murk) then

    ! Calculate using MURK aerosol
!$OMP  PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                              &
!$OMP  SHARED( tdims, aerosol, n_drop_tpr, land_fract, snow_depth, n0_murk,    &
!$OMP          m0_murk  )                                                      &
!$OMP  private( i, j, k, n_aer)
    do k = 1,tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end

          ! Use Jones-Clark Formulae
          n_aer = max( aerosol(i,j, k) / m0_murk * 1.0e-9, 0.0001)
          ! 1.0E-9 converts from ug/kg to kg/kg

          !-----------------------------------------------
          ! Calculation of the aerosol number
          !-----------------------------------------------
          n_aer = n0_murk * n_aer ** power_murk
          n_drop_tpr(i,j,k) = max_drop * ( 1.0e+00-exp( -2.5e-9 * n_aer ) )

        end do
      end do

      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end

          if ( land_fract(i,j) >= land_frac_thresh .and.                       &
               snow_depth(i,j) < snow_depth_thresh    ) then

            n_drop_tpr(i,j,k) = max( n_drop_tpr(i,j,k), min_cdnc_land )

          else

            n_drop_tpr(i,j,k) = max( n_drop_tpr(i,j,k), min_cdnc_sea_ice  )

          end if ! land fract > 0.2 etc

        end do ! i (tdims%i)
      end do ! j (tdims%j)
    end do ! k (tdims%k)
!$OMP end PARALLEL do

  else if (l_mcr_arcl) then

    ! Calculate using aerosol climatologies

!$OMP  PARALLEL DEFAULT(none)                                                  &
!$OMP  SHARED( l_mr_physics, tdims, rho3d, rhodz_dry, rhodz_moist,             &
!$OMP          deltaz )                                                        &
!$OMP  private( i, j, k )
    if (l_mr_physics) then
!$OMP do SCHEDULE(STATIC)
      do k = 1,tdims%k_end
        do j = tdims%j_start, tdims%j_end
          do i = tdims%i_start, tdims%i_end

            ! Step 1: Calculate air density, rho

              ! rho is the dry density
            rho3d(i,j,k) = rhodz_dry(i,j,k) / deltaz(i,j,k)
          end do
        end do
      end do
!$OMP end do

    else ! l_mr_physics

!$OMP do SCHEDULE(STATIC)
      do k = 1,tdims%k_end
        do j = tdims%j_start, tdims%j_end
          do i = tdims%i_start, tdims%i_end

            ! rho is the moist density
            rho3d(i,j,k) = rhodz_moist(i,j,k) / deltaz(i,j,k)

          end do
        end do
      end do
!$OMP end do
    end if  ! l_mr_physics
!$OMP end PARALLEL

    ! Step 2: Call number_droplet routine
    call number_droplet(                                                       &
                         tdims%i_start, tdims%i_end,                           &
                         tdims%j_start, tdims%j_end,                           &
                         1, tdims%k_end,                                       &
                         1, tdims%k_end,                                       &
                         l_mcr_arcl,                                           &
                         l_nh42so4,                                            &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_sulp_ac)),          &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_sulp_di)),          &
                         l_use_arclsslt,                                       &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_sslt_fi)),          &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_sslt_jt)),          &
                         l_use_biogenic, biogenic,                             &
                         l_use_arclbiom,                                       &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_biom_ag)),          &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_biom_ic)),          &
                         l_use_arclocff,                                       &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_ocff_ag)),          &
                         arcl(:,:,:,i_arcl_compnts(ip_arcl_ocff_ic)),          &
                         l_use_arclnitr,                                       &
                         rho3d,                                                &
                         snow_depth, land_fract,                               &
                         n_drop_tpr, aerosol=aerosol)

!$OMP  PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                              &
!$OMP  SHARED( tdims, n_drop_tpr, arcl_inhom_sc )                              &
!$OMP  private( i, j, k )
    do k = 1,tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          ! scaling the n_drop by the inhomogeneity scaling:
          n_drop_tpr(i,j,k) = n_drop_tpr(i,j,k) * arcl_inhom_sc

        end do
      end do
    end do
!$OMP end PARALLEL do

  else if (l_use_sulphate_autoconv .or. l_ukca_aie2 .or.                       &
           l_glomap_clim_aie2 .or. l_easyaerosol_autoconv) then

    !------------------------------------------------------------
    ! Full droplet number calculation based on prognostic aerosol
    !------------------------------------------------------------

!$OMP  PARALLEL DEFAULT(none)                                                  &
!$OMP  SHARED( l_ukca_aie2, l_glomap_clim_aie2, l_mr_physics, tdims,           &
!$OMP          n_drop_tpr,                                                     &
!$OMP          ukca_cdnc, rho3d, rhodz_dry, deltaz, rhodz_moist,               &
!$OMP          l_easyaerosol_autoconv, easyaerosol_cdnc )                      &
!$OMP  private( i, j, k )
    if (l_easyaerosol_autoconv) then
!$OMP do SCHEDULE(STATIC)
      do k = 1, tdims%k_end
        do j = tdims%j_start, tdims%j_end
          do i = tdims%i_start, tdims%i_end
            n_drop_tpr(i,j,k) = easyaerosol_cdnc%cdnc(i,j,k)
          end do
        end do
      end do
!$OMP end do
    else
      if (l_ukca_aie2 .or. l_glomap_clim_aie2) then
!$OMP do SCHEDULE(STATIC)
        do k = 1, tdims%k_end
          do j = tdims%j_start, tdims%j_end
            do i = tdims%i_start, tdims%i_end
              n_drop_tpr(i,j,k) = ukca_cdnc(i,j,k)
            end do
          end do
        end do
!$OMP end do
      else  ! (not l_easyaerosol_autoconv and not
            !                               (l_ukca_aie2 or l_glomap_clim_aie2)

            ! Step 1: Calculate air density, rho

        if (l_mr_physics) then
!$OMP do SCHEDULE(STATIC)
          do k = 1, tdims%k_end
            do j = tdims%j_start, tdims%j_end
              do i = tdims%i_start, tdims%i_end

                ! rho is the dry density

                rho3d(i,j,k) = rhodz_dry(i,j,k) / deltaz(i,j,k)
              end do
            end do
          end do
!$OMP end do
        else ! l_mr_physics
!$OMP do SCHEDULE(STATIC)
          do k = 1, tdims%k_end
            do j = tdims%j_start, tdims%j_end
              do i = tdims%i_start, tdims%i_end

                ! rho is the moist density

                rho3d(i,j,k) = rhodz_moist(i,j,k) / deltaz(i,j,k)
              end do
            end do
          end do
!$OMP end do
        end if  ! l_mr_physics

      end if  ! l_ukca_aie2 or l_glomap_clim_aie2

    end if ! l_easyaerosol_autoconv
!$OMP end PARALLEL

    if (.not. (l_ukca_aie2 .or. l_glomap_clim_aie2) .and. .not.                &
                                                   l_easyaerosol_autoconv) then

      ! Step 2: Call number_droplet routine
      call number_droplet(                                                     &
                         tdims % i_start, tdims % i_end,                       &
                         tdims % j_start, tdims % j_end,                       &
                         1, tdims % k_end,                                     &
                         1, tdims % k_end,                                     &
                         l_use_sulphate_autoconv,                              &
                         l_nh42so4,                                            &
                         so4_acc,                                              &
                         so4_dis,                                              &
                         l_use_seasalt_autoconv,                               &
                         sea_salt_film,                                        &
                         sea_salt_jet,                                         &
                         l_use_biogenic, biogenic,                             &
                         l_use_bmass_autoconv,                                 &
                         bmass_agd,                                            &
                         bmass_cld,                                            &
                         l_use_ocff_autoconv,                                  &
                         ocff_agd, ocff_cld,                                   &
                         l_use_nitrate_autoconv,                               &
                         rho3d,                                                &
                         snow_depth, land_fract,                               &
                         n_drop_tpr,                                           &
                         nitr_acc,                                             &
                         nitr_diss)
    end if ! l_ukca_aie2 or l_glomap_clim_aie2

  else ! (not l_murk, l_mcr_arcl, l_use_sulphate_autoconv or l_ukca_aie2 or
       !  l_glomap_clim_aie2 or l_easyaerosol_autoconv)

    !-----------------------------------------------------------
    ! In this case, STASH 4/211 will be a simple land-sea split
    ! dependant on ntot_land and ntot_sea
    !-----------------------------------------------------------

    ! Set the first level
!$OMP  PARALLEL DEFAULT(none)                                                  &
!$OMP  SHARED( tdims, land_fract, n_drop_tpr )                                 &
!$OMP  private( i, k, j )
    j = 1 ! Arrays passed in by LFRic have a size of 1 in the j dimension.
!$OMP do SCHEDULE(STATIC)
    do i = tdims%i_start, tdims%i_end

      if (land_fract(i,j) >= 0.5) then
        n_drop_tpr(i,j,1) = ntot_land
      else ! land_fract
        n_drop_tpr(i,j,1) = ntot_sea
      end if ! land fract

    end do ! i (tdims%i)
!$OMP end do

    ! copy up to higher levels - thus avoiding too many branches
    ! in loops
!$OMP do SCHEDULE(STATIC)
    do k = 2, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          n_drop_tpr(i,j,k) = n_drop_tpr(i,j,1)
        end do
      end do
    end do
!$OMP end do
!$OMP end PARALLEL

  end if ! l_autoconv_murk, l_use_sulphate_autoconv

end if ! l_droplet_tpr

!---------------------------------------------------------------------
! End of physics
!---------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_taper_ndrop
end module lsp_taper_ndrop_mod
