! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!     subroutine IMP_MIX -----------------------------------------------

!    Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer

!    Purpose: Calculate turbulent mixing increments for a passive tracer
!             using an implicit numerical scheme.  The tridiagonal
!             matices are inverted using simple Gaussian elimination.

!    Programming standard: UM Documentation Paper No. 3

!    Documentation: UM Documentation Paper No 24.

!----------------------------------------------------------------------
module imp_mix_mod
use um_types, only: real_umphys

implicit none
character(len=*), parameter, private :: ModuleName='IMP_MIX_MOD'

contains
subroutine imp_mix (                                                           &
bl_levels,dtrdz,r_theta_levels,r_rho_levels,r_dims,                            &
gamma_rhokh_rdz, gamma_rhok_dep,f_field,surf_dep_flux,field                    &
)

use atm_fields_bounds_mod, only:  pdims, array_dims
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
!$ use omp_lib, only: omp_get_max_threads
implicit none

!  Inputs :-

type(array_dims), intent(in) :: r_dims ! in size of r_X_levels fields

integer, intent(in) ::                                                         &
 bl_levels                   ! in No. of atmospheric levels for
                               !    which boundary layer fluxes are
                               !    calculated.
real(kind=real_umphys), intent(in) ::                                          &
 r_theta_levels(r_dims%i_start:r_dims%i_end,r_dims%j_start:r_dims%j_end,       &
                 0:bl_levels),                                                 &
 r_rho_levels(r_dims%i_start:r_dims%i_end,r_dims%j_start:r_dims%j_end,         &
               bl_levels),                                                     &
                               ! in height of model rho and theta levels
 dtrdz(pdims%i_start:,pdims%j_start:,:)
                               ! in  dt/(rho*r*r*dz) for scalar
                               !     flux divergence

!  Next 2 arrays are in as "explicit" fluxes and out as "implicit"
!  fluxes.

real(kind=real_umphys), intent(in out) ::                                      &
 gamma_rhokh_rdz(pdims%i_start:,pdims%j_start:,2:),                            &
                               ! INOUT Turbulent mixing coefs. above
                               !    surface, =gamma(K)*RHOKH(,K)
                               !    *RDZ(K) for K>=2 (from KMKH).
 gamma_rhok_dep(pdims%i_start:,pdims%j_start:),                                &
                               ! INOUT Surface exchange coefficient
                               !    for surface deposition*gamma(1)
 f_field(pdims%i_start:,pdims%j_start:,:),                                     &
                                       ! INOUT Flux of tracer
 surf_dep_flux(pdims%i_start:,pdims%j_start:),                                 &
                                       ! INOUT surface deposition
                                       !       flux
 field(pdims%i_start:,pdims%j_start:,:)
                                       ! INOUT Amount of tracer

!    Local and other symbolic constants :-
!   Workspace :-
!   4*BL_LEVELS + 4 blocks of real workspace are required.
real(kind=real_umphys) ::                                                      &
 af(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,bl_levels),            &
                                       ! Elements in rows in matrix
                               ! equation (modified during
                               ! Gaussian elimination calculations).
 d_field(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end,                  &
         bl_levels)
                                       ! Delta FIELD (tracer field)
                               ! elements of vector on RHS, then
                               ! LHS, of eqn P244.79.

!  Local scalars :-
real(kind=real_umphys) ::                                                      &
 cf,                                                                           &
            ! Matrix element for local increments
 rbf,                                                                          &
            ! Reciprocal of B for local increments
 r_sq,                                                                         &
            ! r*r
 rr_sq    ! 1/(r*r)

integer ::                                                                     &
 blm1,                                                                         &
            ! BL_LEVELS minus 1.
 i,                                                                            &
            ! Loop counter (horizontal field index).
 j,                                                                            &
            ! Offset version of I.
 k,                                                                            &
            ! Loop counter (vertical index).
 km1,                                                                          &
            ! K minus 1.
 kp1,                                                                          &
            ! K plus 1.
 ii,                                                                           &
            ! omp blocking counter
 pdims_omp_block,                                                              &
            ! omp block size
 pdims_seg_block,                                                              &
            ! omp segment length
 max_threads

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='IMP_MIX'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

blm1 = bl_levels-1

max_threads = 1
!$ max_threads = omp_get_max_threads()
pdims_omp_block = ceiling(real(pdims%i_end)/max_threads)
pdims_seg_block = min(pdims_omp_block, pdims%i_len)

!$OMP PARALLEL DEFAULT(none) SHARED(pdims_seg_block, f_field, dtrdz, field, pdims,   &
!$OMP  gamma_rhokh_rdz, r_rho_levels, surf_dep_flux, d_field, af,              &
!$OMP  r_theta_levels, gamma_rhok_dep, bl_levels, blm1)                        &
!$OMP  private(r_sq, cf, rbf, kp1, km1, ii, k, j, i, rr_sq)

! ----------------------------------------------------------------------
!  (A) Calculations on P-grid.
!-----------------------------------------------------------------------
! 1.0 For simulations on a sphere use spherical geometry for vertical
!     derivatives so multiply fluxes and rho*K_h by r*r
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do k = 2, bl_levels
  do j = pdims%j_start, pdims%j_end
    do i = pdims%i_start, pdims%i_end
      r_sq = r_rho_levels(i,j,k)*r_rho_levels(i,j,k)
      gamma_rhokh_rdz(i,j,k) = r_sq * gamma_rhokh_rdz(i,j,k)
      f_field(i,j,k)         = r_sq * f_field(i,j,k)
    end do
  end do
end do
!$OMP end do

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    r_sq = r_theta_levels(i,j,0)*r_theta_levels(i,j,0)
    gamma_rhok_dep(i,j) = r_sq * gamma_rhok_dep(i,j)
    f_field(i,j,1)      = r_sq * f_field(i,j,1)
    surf_dep_flux(i,j)  = r_sq * surf_dep_flux(i,j)

    ! ----------------------------------------------------------------------
    !  4.  Calculate those matrix and vector elements on the LHS of eqn
    !      which are to do with implicit solution of the tracer
    !      transport problem at the surface and between all levels.
    !      Begin with "upward sweep" through lower half of matrix).
    !-----------------------------------------------------------------------

        !-----------------------------------------------------------------------
        !  4.2 Lowest atmospheric layer FIELD row of matrix.
        !-----------------------------------------------------------------------
                  ! "Explicit" increment to FIELD(1)
    d_field(i,j,1) = -dtrdz(i,j,1) *                                           &
                   ( f_field(i,j,2) - f_field(i,j,1)                           &
                                  - surf_dep_flux(i,j) )

    cf = -dtrdz(i,j,1) * gamma_rhok_dep(i,j)
    af(i,j,1) = -dtrdz(i,j,1) * gamma_rhokh_rdz(i,j,2)
    rbf = 1.0 / ( 1.0 - af(i,j,1) - cf )
    d_field(i,j,1) = rbf * d_field(i,j,1)
    af(i,j,1) = rbf * af(i,j,1)
  end do
end do
!$OMP end do

!-----------------------------------------------------------------------
!  4.3 Rows of matrix applying to FIELD transport into model layers in
!      the "middle" of the "boundary" layer, i.e. all but the bottom
!      layer and the top "boundary" layer.
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do ii = pdims%i_start, pdims%i_end, pdims_seg_block
  do k = 2, blm1
    kp1 = k+1
    km1 = k-1
    do j = pdims%j_start, pdims%j_end
      do i = ii, min(ii+pdims_seg_block-1, pdims%i_end)

        !   "Explicit" flux divergence across layer giving explicit FIELD
        !   increment due to mixing

        d_field(i,j,k) = -dtrdz(i,j,k) *                                       &
                       (f_field(i,j,kp1) - f_field(i,j,k))
        af(i,j,k) = -dtrdz(i,j,k) * gamma_rhokh_rdz(i,j,kp1)
        cf = -dtrdz(i,j,k) * gamma_rhokh_rdz(i,j,k)
        rbf = 1.0 / ( 1.0 - af(i,j,k)                                          &
                          - cf * ( 1.0 + af(i,j,km1) ) )
        d_field(i,j,k) = rbf * ( d_field(i,j,k)                                &
                                  - cf*d_field(i,j,km1) )
        af(i,j,k) = rbf * af(i,j,k)
      end do
    end do
  end do
end do
!$OMP end do

!-----------------------------------------------------------------------
!  4.4 Top "boundary" layer FIELD row of matrix. FIELD for this layer
!      can then be, and is, updated.
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do j = pdims%j_start, pdims%j_end
  do i = pdims%i_start, pdims%i_end
    d_field(i,j,bl_levels) = dtrdz(i,j,bl_levels) *                            &
                             f_field(i,j,bl_levels)

    cf = -dtrdz(i,j,bl_levels) * gamma_rhokh_rdz(i,j,bl_levels)
    rbf = 1.0 / ( 1.0 - cf*( 1.0 + af(i,j,blm1) ) )
    d_field(i,j,bl_levels) = rbf * ( d_field(i,j,bl_levels)                    &
                                        - cf*d_field(i,j,blm1) )
    field(i,j,bl_levels) = field(i,j,bl_levels) +                              &
                           d_field(i,j,bl_levels)
  end do
end do
!$OMP end do

!-----------------------------------------------------------------------
!  5.  "Downward sweep" through whole matrix.  FIELD is updated when
!      the final implicit increments have been calculated.
!-----------------------------------------------------------------------
!  5.1 Remaining FIELD rows of matrix and add implicit increments
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do ii = pdims%i_start, pdims%i_end, pdims_seg_block
  do k = blm1, 1, -1
    do j = pdims%j_start, pdims%j_end
      do i = ii, min(ii+pdims_seg_block-1,pdims%i_end)
        d_field(i,j,k) = d_field(i,j,k) -                                      &
                       af(i,j,k)*d_field(i,j,k+1)
        field(i,j,k) = field(i,j,k) + d_field(i,j,k)
      end do
    end do
  end do
end do
!$OMP end do

!-----------------------------------------------------------------------
!  6.  Calculate final implicit flux of tracer.
!-----------------------------------------------------------------------
! 6.0 First divide by r*r for diagnostics.
!-----------------------------------------------------------------------

!$OMP do SCHEDULE(STATIC)
do k = 1, bl_levels
  if ( k > 1 ) then
    !-----------------------------------------------------------------------
    !  6.1 Fluxes at layer interfaces above the surface.
    !-----------------------------------------------------------------------
    km1 = k-1
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        rr_sq = 1.0 / ( r_rho_levels(i,j,k)*r_rho_levels(i,j,k) )
        gamma_rhokh_rdz(i,j,k) = rr_sq * gamma_rhokh_rdz(i,j,k)
        f_field(i,j,k)         = rr_sq * f_field(i,j,k)
        !  Calculate and store implicit fluxes due to local mixing.

        f_field(i,j,k) = f_field(i,j,k) - gamma_rhokh_rdz(i,j,k)               &
                         * ( d_field(i,j,k) - d_field(i,j,km1) )
      end do
    end do
  else
    !-----------------------------------------------------------------------
    !  6.2 Surface fluxes.
    !-----------------------------------------------------------------------
    do j = pdims%j_start, pdims%j_end
      do i = pdims%i_start, pdims%i_end
        rr_sq = 1.0 / ( r_theta_levels(i,j,0)*r_theta_levels(i,j,0) )
        gamma_rhok_dep(i,j) = rr_sq * gamma_rhok_dep(i,j)
        f_field(i,j,1)      = rr_sq * f_field(i,j,1)
        surf_dep_flux(i,j)  = rr_sq * surf_dep_flux(i,j)

        surf_dep_flux(i,j) = surf_dep_flux(i,j)                                &
                            - gamma_rhok_dep(i,j) * d_field(i,j,1)
      end do
    end do
  end if
end do
!$OMP end do
!$OMP end PARALLEL

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine imp_mix
end module imp_mix_mod
