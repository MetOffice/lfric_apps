! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module bm_ez_diagnosis_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='BM_EZ_DIAGNOSIS_MOD'

contains
! Entrainment zone diagnosis.
! Subroutine Interface:
subroutine bm_ez_diagnosis( p_theta_levels, tgrad_bm, z_theta,                 &
                            ri_bm, zh, zhsc, dzh, bl_type_7, levels, t, q,     &
                            l_mixing_ratio, kez_inv, kez_bottom, kez_top)

use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim
use atm_fields_bounds_mod, only: pdims, tdims
use planet_constants_mod,  only: r,lcrcp, kappa, repsilon, grcp
use water_constants_mod,   only: lc
use pc2_constants_mod,     only: bm_negative_init

use qsat_mod, only: qsat_wat, qsat_wat_mix

use cloud_inputs_mod, only: i_bm_ez_opt, i_bm_ez_orig, i_bm_ez_subcrit,        &
                            ez_max_bm

implicit none

! Purpose:
!   This subroutine is called from either the PC2 initiation scheme or the
!   bimodal cloud scheme and detects entrainment zones near temperature
!   inversions

! Method:
!   Find discontinuities based on local maxima in the theta_liq gradient
!   (i.e. tgrad_bm). If an inversion is detected, search down to find the
!   bottom of the entrainment zone. Keep looking down as long as the
!   gradient monotonically decreases downward, but remains larger than
!   1.1*g/cp. Then search upward from inversion level to find properties
!   of air from above the EZ, identified as at least the level directly
!   above the inversion level, or the driest level anywhere within 250m
!   above the inversion level, whichever level is farthest from the
!   inversion level.


! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_cloud
!
! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP No. 39


!  Global Variables:----------------------------------------------------

!  Subroutine Arguments:------------------------------------------------
integer, intent(in) :: levels
!       No. of levels being processed.

real(kind=real_umphys), intent(in) ::                                          &
   p_theta_levels(pdims%i_start:pdims%i_end,                                   &
                  pdims%j_start:pdims%j_end,levels),                           &
!       pressure at all points (Pa).
   z_theta(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,levels),                           &
!       Height of levels
   tgrad_bm(      tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,levels),                           &
!       Gradient of liquid potential temperature interpolated on theta levels
   ri_bm(         tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,levels),                           &
!       Richardson Number
   zh(            tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Boundary-layer height for bimodal cloud scheme
   zhsc(          tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Decoupled layer height for bimodal cloud scheme
   dzh(           tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Inversion Thickness for bimodal cloud scheme
   bl_type_7(     tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end),                                  &
!       Shear-driven boundary layer indicator
   q(             tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,levels),                           &
!       Total water content (QW) (kg per kg air).
   t(             tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,levels)
!       Liquid/frozen water temperature (TL) (K).

logical, intent(in) :: l_mixing_ratio
!       True if using mixing ratios

integer, intent(out) ::                                                        &
   kez_inv(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,levels),                           &
!       if greater than 1, this is a flag indicating that this level belongs to
!       an entrainment zone (EZ) and its value is set to the k-level of the
!       inversion above this EZ.
   kez_bottom(    tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,levels),                           &
!       if greater than 1, this is the k-level identified as the mode from the
!       bottom of the EZ for the EZ that this level belongs to.
   kez_top(       tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end,levels)
!       if greater than 1, this is the k-level identified as the mode
!       representative for the air above the inversion for the EZ that this
!       level belongs to.

!  Local scalars--------------------------------------------------------

real(kind=real_umphys) ::                                                      &
   zh_eff(        tdims%i_start:tdims%i_end,                                   &
                  tdims%j_start:tdims%j_end)
!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
real(kind=real_umphys) ::                                                      &
 alphal,                                                                       &
                      ! Local gradient of clausius-clapeyron
 alphl,                                                                        &
                      ! repsilon*lc/r
 mux,                                                                          &
                      ! Local first moment of the s-distribution
 mukp1,                                                                        &
                      ! First moment of the s-distribution at level k+1
 alx,                                                                          &
                      ! Local latent-heat correction term
 tlx,                                                                          &
                      ! Local liquid potential temperature
 qs
                      ! Local saturation specific humidity

! Flag for turbulent entrainment zone
logical :: l_turb

!  (b) Others.
integer :: k,i,j, kk  ! Loop counters: K, KK - vertical level index.
                      !                I, J - horizontal field indices.

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='BM_EZ_DIAGNOSIS'

!- End of Header

! ==Main Block==--------------------------------------------------------
! Subroutine structure :
! Loop round levels to be processed.
! ----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

  ! ----------------------------------------------------------------------------
  ! --  Section 1 - initialisations                                           --
  ! ----------------------------------------------------------------------------

alphl=repsilon*lc/r

!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED( tdims, levels, kez_inv, kez_bottom, kez_top, zh_eff, bl_type_7,  &
!$OMP         zh, zhsc, dzh, i_bm_ez_opt )                                     &
!$OMP private( i, j, k )
!$OMP do SCHEDULE(STATIC)
do k = 1, levels
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      kez_inv    (i,j,k) = 1
      kez_bottom (i,j,k) = 1
      kez_top    (i,j,k) = 1
    end do
  end do
end do
!$OMP end do NOWAIT

if (i_bm_ez_opt == i_bm_ez_subcrit) then
!$OMP  do SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      ! Calculate an effective PBL depth
      ! include resolved inversion thickness, dzh
      ! (which is set to rmdi when subgrid)
      zh_eff(i,j) = zh(i,j) + max( dzh(i,j), 0.0 )
      if ( bl_type_7(i,j) > 0.5_real_umphys ) then
        ! use cloud top (if higher) in decoupled PBLs that are
        ! shear-dominated
        zh_eff(i,j) = max( zh_eff(i,j), zhsc(i,j) )
      end if
    end do
  end do
!$OMP end do
else
!$OMP do SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end
      zh_eff(i,j)   = 0.0
    end do
  end do
!$OMP end do
end if
!$OMP end PARALLEL

! ----------------------------------------------------------------------------
! --  Section 2 - Detection of entrainment zones                            --
! -- Find discontinuities based on local maxima in the theta_liq gradient   --
! -- (i.e. tgrad_bm). If an inversion is detected, search down to find the  --
! -- bottom of the entrainment zone. Keep looking down as long as the       --
! -- gradient monotonically decreases downward, but remains larger than     --
! -- 1.1*g/cp. Then search upward from inversion level to find properties   --
! -- of air from above the EZ, identified as at least the level directly    --
! -- above the inversion level, or the driest level anywhere within 250m    --
! -- above the inversion level, whichever level is farthest from the        --
! -- inversion level.                                                       --
! -- Avoid spurious entrainment zones by requiring levels above and below   --
! -- to all have initialised values of tgrad_bm (k and k+2 values >         --
! -- bm_negative_init)                                                      --
! -- Inversion are only searched for in the lowest 6 km above the surface   --
! -- (in case l_bm_ez_subcrit_only is true, otherwise, inversions are only  --
! -- searched for in the lowest 3 km above the surface, by not initialising --
! -- tgrad_bm above this level in bdy_expl2).                               --
! -- The depth of the entrainment zone is restricted to be ez_max_bm, as    --
! -- set in the namelist, typically between 300 and 500m                    --
! ----------------------------------------------------------------------------

!$OMP  PARALLEL                                                                &
!$OMP  DEFAULT(none)                                                           &
!$OMP  SHARED(tdims,t,q,p_theta_levels,l_mixing_ratio,grcp,                    &
!$OMP  tgrad_bm,lcrcp,kappa,repsilon,r,z_theta,alphl,levels,ri_bm,             &
!$OMP  zh_eff,i_bm_ez_opt,kez_top,kez_bottom,kez_inv,ez_max_bm)                &
!$OMP  private(j,i,k,kk,qs,alphal,alx,tlx,mux,mukp1,l_turb)
do k = 2, levels-3
!$OMP do SCHEDULE(DYNAMIC)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end

      if (tgrad_bm(i,j,k+1) > tgrad_bm(i,j,k)   .and.                          &
          tgrad_bm(i,j,k+1) > tgrad_bm(i,j,k+2) .and.                          &
          tgrad_bm(i,j,k+1) > 0.1*grcp          .and.                          &
          tgrad_bm(i,j,k+2) > bm_negative_init  .and.                          &
          tgrad_bm(i,j,k)   > 0.0               .and.                          &
          ((z_theta(i,j,k)    < 6.0e3 .and. i_bm_ez_opt==i_bm_ez_subcrit) .or. &
          (i_bm_ez_opt==i_bm_ez_orig))) then

        ! inversion found, so preset the level of the bottom of the EZ and
        ! the level representative of the air above the EZ.
        kez_inv(i,j,k)    = k
        kez_bottom(i,j,k) = k-1
        kez_top(i,j,k)    = k+1

        ! Look down to find bottom of EZ and set kez_inv flag for each level
        ! being part of the EZ to value k. Look down as long as the liquid
        ! potential temperature gradient is monotonically decreasing, but
        ! remains larger than 0.1*grcp. Never look further down than ez_max_bm
        ! below the inversion.
        kk=k-1
        kez_inv(i,j,kk) = k

        do while ( tgrad_bm(i,j,kk) < tgrad_bm(i,j,kk+1)        .and.          &
                   tgrad_bm(i,j,kk) > 0.1*grcp                  .and.          &
                   kk > 2                                       .and.          &
                   z_theta(i,j,k)-z_theta(i,j,kk) < ez_max_bm )

          kez_inv(i,j,kk)   = k
          kez_bottom(i,j,k) = kk
          kk=kk-1

        end do

        ! Look up to find level representative of air above the EZ. Start
        ! with presetting the mean saturation departure for level k and k+1,
        ! then keep looking upward as long as the free troposphere is
        ! monotonically drying, but don't look further than 250m above the
        ! discontinuity. Use at least one level above inversion
        ! (k+2, since inversion is at k+1).

        tlx = t(i,j,k)*(p_theta_levels(i,j,k)/                                 &
              p_theta_levels(i,j,k))**kappa
        ! Calculate the saturation specific humidity for layers
        if ( l_mixing_ratio ) then
          call qsat_wat_mix(qs,tlx,p_theta_levels(i,j,k))
        else
          call qsat_wat(qs,tlx,p_theta_levels(i,j,k))
        end if
        alphal = alphl * qs / (tlx * tlx)
        alx = 1.0 / (1.0 + (lcrcp * alphal))
        mux   = alx*(q(i,j,k) - qs)

        tlx = t(i,j,k+1)*(p_theta_levels(i,j,k)/                               &
              p_theta_levels(i,j,k+1))**kappa
        ! Calculate the saturation specific humidity for layers
        if ( l_mixing_ratio ) then
          call qsat_wat_mix(qs,tlx,p_theta_levels(i,j,k))
        else
          call qsat_wat(qs,tlx,p_theta_levels(i,j,k))
        end if
        alphal = alphl * qs / (tlx * tlx)
        alx = 1.0 / (1.0 + (lcrcp * alphal))
        mukp1 = alx*(q(i,j,k+1) - qs)

        kk=k+1
        do while ( (kk < levels-1) .and.                                       &
                   ((z_theta(i,j,kk) < 1.25*z_theta(i,j,k+1)         .and.     &
                     z_theta(i,j,kk)-z_theta(i,j,k+1) < 250.0        .and.     &
                     mux > mukp1)                                              &
                    .or. (kk == k+2) .or. (kk == k+1)) )


          tlx = t(i,j,kk)*(p_theta_levels(i,j,k)/                              &
                p_theta_levels(i,j,kk))**kappa
          ! Calculate the saturation specific humidity for layers
          if ( l_mixing_ratio ) then
            call qsat_wat_mix(qs,tlx,p_theta_levels(i,j,k))
          else
            call qsat_wat(qs,tlx,p_theta_levels(i,j,k))
          end if
          alphal = alphl * qs / (tlx * tlx)
          alx = 1.0 / (1.0 + (lcrcp * alphal))
          mux   = alx*(q(i,j,kk) - qs)

          tlx = t(i,j,kk+1)*(p_theta_levels(i,j,k)/                            &
                p_theta_levels(i,j,kk+1))**kappa
          ! Calculate the saturation specific humidity for layers
          if ( l_mixing_ratio ) then
            call qsat_wat_mix(qs,tlx,p_theta_levels(i,j,k))
          else
            call qsat_wat(qs,tlx,p_theta_levels(i,j,k))
          end if
          alphal = alphl * qs / (tlx * tlx)
          alx = 1.0 / (1.0 + (lcrcp * alphal))
          mukp1 = alx*(q(i,j,kk+1) - qs)

          kez_top(i,j,k) = kk

          kk=kk+1

        end do

        if ( i_bm_ez_opt==i_bm_ez_subcrit ) then
          ! Option to only diagnose entrainment zones if turbulent

          ! Diagnose whether this entrainment zone is turbulent
          l_turb = .false.
          do kk = kez_bottom(i,j,k), kez_top(i,j,k)
            ! Turbulent if any part of the zone has a sub-critical Ri
            ! (using long-tailed critical Richardson number of 1.0)
            if ( ri_bm(i,j,kk) < 1.0 ) then
              l_turb = .true.
            end if
          end do

          if ( .not. ( l_turb .or.                                             &
                       z_theta(i,j,kez_bottom(i,j,k)) <= zh_eff(i,j) ) ) then
            ! Only allow entrainment zones to be diagnosed if they
            ! overlap with a sub-critical layer (and so might be turbulent),
            ! or if they overlap with boundary-layer non-local fluxes.
            ! Reset entrainment zone variables to initialised value 1 if not.
            do kk = kez_bottom(i,j,k), kez_top(i,j,k)
              kez_inv(i,j,kk) = 1
            end do
            kez_top(i,j,k) = 1
            kez_bottom(i,j,k) = 1
          end if

        end if  ! ( i_bm_ez_opt==i_bm_ez_subcrit )

      end if
    end do
  end do
!$OMP end do
end do
!$OMP end PARALLEL

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine bm_ez_diagnosis
end module bm_ez_diagnosis_mod
