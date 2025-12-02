! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Microphysics precipitation scheme. Air density calculation

module mphys_air_density_mod
! Description:
! Calculates air density in one place for use within the microphysics

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='MPHYS_AIR_DENSITY_MOD'

contains

! Subroutine Interface:
subroutine mphys_air_density( r_theta_levels, r_rho_levels,                    &
                              rho_dry, rho_r2, qdims,                          &
                              q, qcl, qcf, qcf2, qrain, qgraup,                &
                              rhodz_dry, rhodz_moist, deltaz )

use atm_fields_bounds_mod, only: array_dims, tdims, pdims_s, tdims_l

use mphys_inputs_mod,      only: l_mcr_qrain, l_mcr_qcf2, l_mcr_qgraup,        &
                                 l_mphys_nonshallow

use gen_phys_inputs_mod,   only: l_mr_physics, l_vol_interp_rho

! Dr Hook Modules
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim

implicit none

real(kind=real_umphys), intent(in)  ::                                         &
                    r_theta_levels( tdims_l%i_start:tdims_l%i_end,             &
                                    tdims_l%j_start:tdims_l%j_end,             &
                                    tdims_l%k_start:tdims_l%k_end ),           &
                                    ! in height of theta levels
                    r_rho_levels( tdims_l%i_start:tdims_l%i_end,               &
                                  tdims_l%j_start:tdims_l%j_end,               &
                                  1:tdims_l%k_end ),                           &
                                  ! in height of rho levels
                    rho_dry( pdims_s%i_start : pdims_s%i_end,                  &
                             pdims_s%j_start : pdims_s%j_end,                  &
                             pdims_s%k_start : pdims_s%k_end )
                             ! in Unscaled dry-density / kg m-3

real(kind=real_umphys), intent(in)  ::                                         &
                     rho_r2( pdims_s%i_start : pdims_s%i_end,                  &
                             pdims_s%j_start : pdims_s%j_end,                  &
                             pdims_s%k_start : pdims_s%k_end )
                             ! in Air density * earth radius**2

! Structure storing dimensions of the q arrays, to allow flexibility
! to call this routine using arrays with or without halos.
type(array_dims), intent(in) :: qdims

real(kind=real_umphys), intent(in) :: q    (qdims%i_start:qdims%i_end,         &
                                            qdims%j_start:qdims%j_end,         &
                                            qdims%k_start:qdims%k_end)
                                       ! Water vapour (kg per kg air).
real(kind=real_umphys), intent(in) :: qcf  (qdims%i_start:qdims%i_end,         &
                                            qdims%j_start:qdims%j_end,         &
                                            qdims%k_start:qdims%k_end)
                                       ! Cloud ice (kg per kg air).
real(kind=real_umphys), intent(in) :: qcl  (qdims%i_start:qdims%i_end,         &
                                            qdims%j_start:qdims%j_end,         &
                                            qdims%k_start:qdims%k_end)
                                       ! Cloud liquid water (kg per kg air)
real(kind=real_umphys), intent(in) ::  qcf2  (qdims%i_start:qdims%i_end,       &
                                              qdims%j_start:qdims%j_end,       &
                                              qdims%k_start:qdims%k_end)
                                       ! Ice (kg per kg air)
real(kind=real_umphys), intent(in) ::  qrain (qdims%i_start:qdims%i_end,       &
                                              qdims%j_start:qdims%j_end,       &
                                              qdims%k_start:qdims%k_end)
                                       ! Rain (kg per kg air)
real(kind=real_umphys), intent(in) ::  qgraup(qdims%i_start:qdims%i_end,       &
                                              qdims%j_start:qdims%j_end,       &
                                              qdims%k_start:qdims%k_end)

real(kind=real_umphys), intent(out) ::                                         &
                      rhodz_dry( tdims%i_start : tdims%i_end,                  &
                                 tdims%j_start : tdims%j_end,                  &
                                             1 : tdims%k_end )
                                       ! Dry density

real(kind=real_umphys), intent(out) ::                                         &
                      rhodz_moist( tdims%i_start : tdims%i_end,                &
                                   tdims%j_start : tdims%j_end,                &
                                               1 : tdims%k_end )
                                       ! Moist density

real(kind=real_umphys), intent(out) ::                                         &
                      deltaz( tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end,                     &
                                          1 : tdims%k_end )

                                       ! Layer thickness

!------------------------------------------------------------------------------
! Purpose:
!   Calculates air density in a single place, for use in either the 3D or
!   CASIM microphysics
!   Documentation: UMDP 26.

!------------------------------------------------------------------------------
! Subroutine Arguments


!------------------------------------------------------------------------------
! Local Variables

real(kind=real_umphys) :: rho1
real(kind=real_umphys) :: rho2

real(kind=real_umphys) :: q_total

! Factors used in layer-mass calculations for spherical coordinates
real :: r_sq_factor  ! Geometrical factor ( R(k) / R(surf) )^2
real :: interp       ! Interpolation weight for theta-levels from rho-levels
real :: rho_th       ! Density interpolated to theta-levels

integer :: i, j, k ! Loop counters

! Declarations for Dr Hook:
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle
character(len=*), parameter   :: RoutineName='MPHYS_AIR_DENSITY'

!==============================================================================
! Start of calculations
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! ----------------------------------------------------------------------
! Calculate the (non-hydrostatic) layer thicknesses (deltaz) and air
! densities multiplied by deltaz (rhodz_moist and rhodz_dry).
! ----------------------------------------------------------------------
! We should note that this formulation, although better than the
! hydrostatic formulation, is still not entirely conservative. To ensure
! conservation we would need to rewrite the large-scale precipitation
! scheme to consider masses in terms of rho<q>, and
! not the current <rho>q formulation.

! We only need to calculate averages for the moist levels

if ( l_mphys_nonshallow ) then
  ! Calculations accounting for the curvature of the Earth, which means
  ! the atmospheric mass on a model-level increases with the width of the
  ! grid-box as we rise further away from the centre of the Earth.
  ! To get the precip fall calculations to conserve moisture in the
  ! global model's spherical coordinate framework, we need to multiply
  ! the rho deltaz factors by ( R(k) / R(surf) )^2, where R is the
  ! height above the centre of the Earth.
  ! Also, under this switch:
  ! a) The vertical interpolation of density and the
  !    treatment of the upper boundary are slightly modified to be
  !    consistent with the densities and layer thicknesses used in the
  !    boundary-layer scheme.
  ! b) We use the start-of-timestep dry-density passed in from top-level
  !    instead of recomputing it from the water species fields.
  ! So this switch changes answers even if the model is in
  ! Cartesian coordinates.


  ! Compute deltaz (thickness of the model theta-levels)
  ! Note special treatment at the top (next rho-level doesn't exist)
  ! and bottom (lower boundary is the surface, not the intervening rho-level).
!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP SHARED  ( tdims, r_rho_levels, r_theta_levels, deltaz,                   &
!$OMP           rho_dry, rhodz_dry, rho_r2, rhodz_moist, l_vol_interp_rho )    &
!$OMP private ( i, j, k, interp, rho_th )
  k = 1
  j = 1 ! Arrays passed in by LFRic have a size of 1 in the j dimension.
!$OMP do SCHEDULE(STATIC)
  do i = tdims%i_start, tdims%i_end
    deltaz(i,j,k) = r_rho_levels(i,j,k+1) - r_theta_levels(i,j,0)
  end do
!$OMP end do NOWAIT
!$OMP do SCHEDULE(STATIC)
  do k = 2, tdims%k_end - 1
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        deltaz(i,j,k) = r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k)
      end do
    end do
  end do
!$OMP end do NOWAIT
  k = tdims%k_end
  j = 1 ! Arrays passed in by LFRic have a size of 1 in the j dimension.
!$OMP do SCHEDULE(STATIC)
  do i = tdims%i_start, tdims%i_end
    deltaz(i,j,k) = r_theta_levels(i,j,k) - r_rho_levels(i,j,k)
  end do
!$OMP end do

  ! Interpolate densities onto theta-levels and scale by deltaz
  if (l_vol_interp_rho) then
!$OMP do SCHEDULE(STATIC)
    do k = 1, tdims%k_end - 1
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end

          interp = ( r_theta_levels(i,j,k) - r_rho_levels(i,j,k) )             &
                 / ( r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k) )

          rho_th =        interp * rho_dry(i,j,k)                              &
                 + (1.0-interp)  * rho_dry(i,j,k+1)
          rhodz_dry(i,j,k) = rho_th * deltaz(i,j,k)

          rho_th = interp * rho_r2(i,j,k) / (r_rho_levels(i,j,k)               &
                                            *r_rho_levels(i,j,k))              &
                 + (1.0-interp) * rho_r2(i,j,k+1) / (r_rho_levels(i,j,k+1)     &
                                                    *r_rho_levels(i,j,k+1))
          rhodz_moist(i,j,k) = rho_th * deltaz(i,j,k)

        end do
      end do
    end do
!$OMP end do NOWAIT
  else ! l_vol_interp_rho
!$OMP do SCHEDULE(STATIC)
    do k = 1, tdims%k_end - 1
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end

          interp = ( r_theta_levels(i,j,k) - r_rho_levels(i,j,k) )             &
                 / ( r_rho_levels(i,j,k+1) - r_rho_levels(i,j,k) )

          rho_th = (1.0-interp) * rho_dry(i,j,k)                               &
                 +      interp  * rho_dry(i,j,k+1)
          rhodz_dry(i,j,k) = rho_th * deltaz(i,j,k)

          rho_th = (1.0-interp) * rho_r2(i,j,k) / (r_rho_levels(i,j,k)         &
                                                  *r_rho_levels(i,j,k))        &
                 +      interp  * rho_r2(i,j,k+1) / (r_rho_levels(i,j,k+1)     &
                                                    *r_rho_levels(i,j,k+1))
          rhodz_moist(i,j,k) = rho_th * deltaz(i,j,k)

        end do
      end do
    end do
!$OMP end do NOWAIT
  end if ! l_vol_interp_rho

  ! Special case of top model-level (rho(k+1) doesn't exist)
  k = tdims%k_end
  j = 1 ! Arrays passed in by LFRic have a size of 1 in the j dimension.
!$OMP do SCHEDULE(STATIC)
  do i = tdims%i_start, tdims%i_end

    rhodz_dry(i,j,k) = rho_dry(i,j,k) * deltaz(i,j,k)

    rhodz_moist(i,j,k) = ( rho_r2(i,j,k) / (r_rho_levels(i,j,k)              &
                          * r_rho_levels(i,j,k)) )                           &
                          * deltaz(i,j,k)

  end do
!$OMP end do

!$OMP end PARALLEL

  ! If using spherical coordinates, include factor (r/r_surf)^2
  ! in the layer mass, to account for the grid-columns getting
  ! wider with height.  This means the precip rates are all
  ! defined per unit area at Earth's surface, not per unit area at
  ! the current level.

!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none)                               &
!$OMP SHARED  ( tdims, r_theta_levels, rhodz_dry, rhodz_moist )                &
!$OMP private ( i, j, k, r_sq_factor )
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end

        r_sq_factor = ( r_theta_levels(i,j,k) * r_theta_levels(i,j,k) )      &
                    / ( r_theta_levels(i,j,0) * r_theta_levels(i,j,0) )

        rhodz_dry(i,j,k)   = rhodz_dry(i,j,k)   * r_sq_factor
        rhodz_moist(i,j,k) = rhodz_moist(i,j,k) * r_sq_factor

      end do
    end do
  end do
!$OMP end PARALLEL do


else  ! ( l_mphys_nonshallow )
  ! Original calculation of layer-masses; assumes Cartesian geometry
  ! (which is not quite right except in idealised runs), and the
  ! vertical interpolation of density is different.


!$OMP  PARALLEL do SCHEDULE(STATIC) DEFAULT(SHARED)                            &
!$OMP  private(i,j,k,q_total,rho1,rho2)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end

              ! Calculate densities at the boundaries of the layer
              ! by removing the r**2 term from rho_r2.
              ! Rho1 is the density at the lower boundary.
        rho1 = rho_r2(i,j,k)/( r_rho_levels(i,j,k) * r_rho_levels(i,j,k) )

              ! Check whether there is a rho level above the current
              ! moist level.
        if ( k  <   tdims%k_end ) then
                ! Rho2 is the density at the upper boundary.
          rho2 = rho_r2(i,j,k+1)/( r_rho_levels(i,j,k+1)                       &
                                 * r_rho_levels(i,j,k+1) )

                ! Calculate the average value of rho across the layer
                ! multiplied by the layer thickness and the layer
                ! thickness.
          if (l_vol_interp_rho) then
            rhodz_moist(i,j,k) =                                               &
                        rho1 * ( r_theta_levels(i,j,k) -                       &
                                           r_rho_levels(i,j,k) )               &
                     +  rho2 * ( r_rho_levels(i,j,k+1) -                       &
                                      r_theta_levels(i,j,k) )
          else
            rhodz_moist(i,j,k) =                                               &
                        rho2 * ( r_theta_levels(i,j,k) -                       &
                                           r_rho_levels(i,j,k) )               &
                     +  rho1 * ( r_rho_levels(i,j,k+1) -                       &
                                      r_theta_levels(i,j,k) )
          end if
          deltaz(i,j,k) = r_rho_levels(i,j,k+1)                                &
                            - r_rho_levels(i,j,k)

          if (k  ==  1) then
            ! For the lowest layer we need to extend the lower
            ! boundary from the first rho level to the surface.
            ! The surface is the 0'th theta level.
            deltaz(i,j,1) = r_rho_levels(i,j,2)                                &
                              - r_theta_levels(i,j,0)
            rhodz_moist(i,j,1) = rhodz_moist(i,j,1)*deltaz(i,j,1)              &
                     / (r_rho_levels(i,j,2)-r_rho_levels(i,j,1))
          end if  ! k  ==  1

        else
          ! For a top layer higher than the highest rho level
          ! we can calculate a pseudo rho level. We will assume
          ! it has a similar density to the rho level below
          ! and that the intervening theta level is in the centre
          ! of the layer.
          deltaz(i,j,k) = 2.0*(r_theta_levels(i,j,k)                           &
                                -r_rho_levels(i,j,k))
          rhodz_moist(i,j,k) = rho1 * deltaz(i,j,k)

        end if  ! k  <   tdims%k_end

        ! Calculate total moisture
        q_total = q(i,j,k) + qcl(i,j,k) + qcf(i,j,k)

        if (l_mcr_qcf2) then
          q_total = q_total + qcf2(i,j,k)
        end if  ! l_mcr_qcf2

        if (l_mcr_qrain) then
          q_total = q_total + qrain(i,j,k)
        end if  ! l_mcr_qrain

        if (l_mcr_qgraup) then
          q_total = q_total + qgraup(i,j,k)
        end if  ! l_mcr_qgraup


        ! Rho_r2 uses the moist density of air. If the mixing
        ! ratio framework is in place then we need to also know
        ! the dry density of air.
        if (l_mr_physics) then
          rhodz_dry(i,j,k) = rhodz_moist(i,j,k) / (1.0 + q_total)
        else
          rhodz_dry(i,j,k) = rhodz_moist(i,j,k) * (1.0 - q_total)
        end if  ! l_mr_physics

      end do  ! i
    end do  ! j
  end do  ! k
!$OMP end PARALLEL do


end if  ! ( l_mphys_nonshallow )

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine mphys_air_density
end module mphys_air_density_mod
