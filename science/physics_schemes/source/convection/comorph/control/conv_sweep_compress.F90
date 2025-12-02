! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module conv_sweep_compress_mod

implicit none

contains


! Subroutine to compress all the fields required by conv_level_step
! onto the current compression list of convecting points
subroutine conv_sweep_compress( k, k_half, dk, max_points, ij_first, ij_last,  &
                                n_fields_tot, l_last_level, cmpr,              &
                                grid, turb, fields,                            &
                                virt_temp, layer_mass, sum_massflux,           &
                                l_within_bl, grid_k_super, grid_half_super,    &
                                env_k_fields, env_k_super, env_half_super,     &
                                layer_mass_step, frac_level_step,              &
                                sum_massflux_cmpr, index_ij )

use comorph_constants_mod, only: real_cvprec, real_hmprec, zero,               &
                                 nx_full, ny_full, k_bot_conv, k_top_conv,     &
                                 l_cv_cloudfrac
use cmpr_type_mod, only: cmpr_type
use grid_type_mod, only: grid_type, i_height, n_grid, grid_compress
use turb_type_mod, only: turb_type
use fields_type_mod, only: fields_type, field_positive, i_cf_liq, i_cf_bulk,   &
                           i_wind_w
use env_half_mod, only: n_env_half, i_virt_temp, i_wind_w_half,                &
                        env_half_interp

use set_l_within_bl_mod, only: set_l_within_bl
use compress_mod, only: compress
use force_cloudfrac_consistency_mod, only: force_cloudfrac_consistency

implicit none

! Index of current full model-level
integer, intent(in) :: k
! Index of model-level interface
integer, intent(in) :: k_half
! k-increment to next full level beyond k_half (+ or - 1)
integer, intent(in) :: dk

! Size of compression arrays
integer, intent(in) :: max_points
! i,j indices of the first and last point on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Number of primary fields, accounting for whether doing tracer transport
integer, intent(in) :: n_fields_tot

! Flag for reached the last model-level
logical, intent(in) :: l_last_level

! Structure storing compression indices
type(cmpr_type), intent(in) :: cmpr

! Structure containing pointers to model grid fields
! (full 3-D arrays, possibly with halos);
! contains model-level heights, pressures and dry-density
type(grid_type), intent(in) :: grid

! Structure containing pointers to turbulence fields passed
! in from the boundary-layer scheme
type(turb_type), intent(in) :: turb

! Structure containing pointers to the _np1 fields;
! these are the primary fields already updated with any other
! increments computed before convection.
type(fields_type), intent(in) :: fields

! Full 3-D array of environment virtual temperature
real(kind=real_hmprec), intent(in) :: virt_temp                                &
       ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Full 3-D array of dry-mass per unit surface area
! contained in each grid-cell
real(kind=real_hmprec), intent(in) :: layer_mass                               &
       ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Sum of parcel mass-fluxes over types / layers
real(kind=real_cvprec), intent(in) :: sum_massflux(ij_first:ij_last)

! OUTPUT COMPRESSED ARRARYS...

! Flag for whether level k is within the boundary-layer
! (as defined by the BL-scheme's BL-top height passed in)
logical, intent(out) :: l_within_bl(max_points)

! Height and pressure
real(kind=real_cvprec), intent(out) :: grid_k_super(max_points,n_grid)
real(kind=real_cvprec), intent(out) :: grid_half_super(max_points,n_grid)

! Environment primary fields at the current level
real(kind=real_cvprec), intent(out) :: env_k_fields(max_points,n_fields_tot)

! Arrays containing certain environment fields which are needed
! at both half-levels and full-levels
real(kind=real_cvprec), intent(out) :: env_k_super(max_points,n_env_half)
real(kind=real_cvprec), intent(out) :: env_half_super(max_points,n_env_half)

! Layer-mass on current half-level-step
real(kind=real_cvprec), intent(out) :: layer_mass_step(max_points)
! Fraction of the layer-mass on current half of the level-step
real(kind=real_cvprec), intent(out) :: frac_level_step(max_points)
! Sum of mass-fluxes, compressed onto points where current type is active
real(kind=real_cvprec), intent(out) :: sum_massflux_cmpr(max_points)

! Indices used to reference a collapsed horizontal coordinate from the parcel
integer, intent(out) :: index_ij(max_points)

! Array lower and upper bounds
integer :: lb(3), ub(3), lb2(2), ub2(2)

! Loop counters
integer :: ic, i_field


! Set flag for whether level k is within the
! boundary-layer at each point.
lb = lbound( grid % height_half )
ub = ubound( grid % height_half )
lb2 = lbound( turb % z_bl_top )
ub2 = ubound( turb % z_bl_top )
call set_l_within_bl( cmpr, k, lb, ub, grid % height_half,                     &
                      lb2, ub2, turb % z_bl_top,                               &
                      l_within_bl )

! Compress grid fields into super-arrays for the
! current level k and the model-level interface
call grid_compress( grid, cmpr,                                                &
                    k     =k,      grid_k_super   =grid_k_super,               &
                    k_half=k_half, grid_half_super=grid_half_super )

! Compress primary fields (optionally includes tracers)
do i_field = 1, n_fields_tot
  lb = lbound(fields%list(i_field)%pt)
  ub = ubound(fields%list(i_field)%pt)
  call compress( cmpr, lb(1:2), ub(1:2), fields%list(i_field)%pt(:,:,k),       &
                 env_k_fields(:,i_field) )
  ! Remove any spurious negative values from input data
  if ( field_positive(i_field) ) then
    do ic = 1, cmpr%n_points
      env_k_fields(ic,i_field) = max( env_k_fields(ic,i_field), zero )
    end do
  end if
end do
! If cloud-fractions included...
if ( l_cv_cloudfrac ) then
  ! Rounding errors when converting the cloud-fractions
  ! to 32-bit in compress can cause them to become
  ! slightly inconsistent; correct if needed:
  call force_cloudfrac_consistency( cmpr%n_points, max_points,                 &
                                    env_k_fields(:,i_cf_liq:i_cf_bulk) )
end if

! Set environment virtual temperature and w wind at level k
lb = [1,1,k_bot_conv]
ub = [nx_full,ny_full,k_top_conv]
call compress( cmpr, lb(1:2), ub(1:2), virt_temp(:,:,k),                       &
               env_k_super(:,i_virt_temp) )
do ic = 1, cmpr%n_points
  env_k_super(ic,i_wind_w_half) = env_k_fields(ic,i_wind_w)
end do

! Compress / interpolate required environment fields onto the model-level
! interface (k+1/2 or k-1/2, depending on whether going up or down)
call env_half_interp( k+dk, max_points, l_last_level, cmpr,                    &
                      grid % height_full,                                      &
                      grid_k_super(:,i_height), grid_half_super(:,i_height),   &
                      fields % wind_w, virt_temp,                              &
                      env_k_super(:,i_wind_w_half), env_k_super(:,i_virt_temp),&
                      env_half_super )

! Compress layer-mass on level k
call compress( cmpr, lb(1:2), ub(1:2), layer_mass(:,:,k),                      &
               layer_mass_step )

! Compress height of half-level the other side of k from k_half
lb = lbound(grid%height_half)
ub = ubound(grid%height_half)
call compress( cmpr, lb(1:2), ub(1:2), grid%height_half(:,:,k_half-dk),        &
               frac_level_step )
do ic = 1, cmpr%n_points
  ! Set fraction of layer k on curent half-level-step
  frac_level_step(ic)                                                          &
    = ( grid_k_super(ic,i_height) - grid_half_super(ic,i_height) )             &
    / ( frac_level_step(ic)       - grid_half_super(ic,i_height) )
  ! Scale layer_mass_step to get mass on current half-level-step
  layer_mass_step(ic) = layer_mass_step(ic) * frac_level_step(ic)
end do

! Set indices for referencing collapsed horizontal coordinate
do ic = 1, cmpr%n_points
  index_ij(ic) = nx_full*(cmpr%index_j(ic)-1) + cmpr%index_i(ic)
end do

! Compress sums of mass-fluxes over types and layers
do ic = 1, cmpr%n_points
  sum_massflux_cmpr(ic) = sum_massflux(index_ij(ic))
end do


return
end subroutine conv_sweep_compress

end module conv_sweep_compress_mod
