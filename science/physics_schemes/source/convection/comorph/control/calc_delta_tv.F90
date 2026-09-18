! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_delta_tv_mod

implicit none

contains

! Subroutine to pre-estimate the subsidence warming per unit mass-flux,
! for use in the implicit detrainment calculations.
subroutine calc_delta_tv( l_kpdk, l_last_level,                                &
                          cmpr, k, k_next, dk, ij_first, ij_last,              &
                          virt_temp, layer_mass, grid, fields,                 &
                          delta_tv )

use comorph_constants_mod, only: real_cvprec, real_hmprec,                     &
                                 zero, one,                                    &
                                 nx_full, ny_full, k_bot_conv, k_top_conv,     &
                                 comorph_timestep, alpha_detrain
use cmpr_type_mod, only: cmpr_type
use grid_type_mod, only: grid_type, n_grid, i_height, i_pressure,              &
                         grid_compress
use fields_type_mod, only: fields_type, i_temperature, i_q_vap,                &
                           i_qc_first, i_qc_last, field_positive

use compress_mod, only: compress
use dry_adiabat_mod, only: dry_adiabat

implicit none

! Flag for 2nd call when estimating subsidence increment at next
logical, intent(in) :: l_kpdk
! Flag for last model-level (in-which case k+dk doesn't exist)
logical, intent(in) :: l_last_level

! Compression indices for points where this needs to be calculated
type(cmpr_type), intent(in) :: cmpr

! Index of current full-level
integer, intent(in) :: k
! Index of next model-level interface
integer, intent(in) :: k_next
! Vertical level increment (+1 for updrafts, -1 for downdrafts)
integer, intent(in) :: dk

! i,j indices of the first and last point on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Full 3-D array of environment virtual temperature
real(kind=real_hmprec), intent(in) :: virt_temp                                &
       ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Full 3-D array of dry-mass per unit surface area contained in each grid-cell
real(kind=real_hmprec), intent(in) :: layer_mass                               &
       ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Structure containing pointers to model grid fields
! (full 3-D arrays, possibly with halos);
! contains model-level heights, pressures and dry-density
type(grid_type), intent(in) :: grid

! Structure containing pointers to the _np1 fields;
! these are the primary fields already updated with any other
! increments computed before convection.
type(fields_type), intent(in) :: fields

! Estimated subsidence warming per unit mass-flux
real(kind=real_cvprec), intent(in out) :: delta_tv ( ij_first:ij_last )

! Compressed copies of required inputs
real(kind=real_cvprec) :: layer_mass_k(cmpr%n_points)
real(kind=real_cvprec) :: height_work(cmpr%n_points)
real(kind=real_cvprec) :: grid_1(cmpr%n_points,n_grid)
real(kind=real_cvprec) :: grid_2(cmpr%n_points,n_grid)
real(kind=real_cvprec) :: virt_temp_1(cmpr%n_points)
real(kind=real_cvprec) :: virt_temp_2(cmpr%n_points)
real(kind=real_cvprec) :: fields_k(cmpr%n_points,i_temperature:i_qc_last)

! Fraction of level-step contained on current half-level
real(kind=real_cvprec) :: frac_level_step(cmpr%n_points)

! Exner ratio for subsidence of Tv over half-level step
real(kind=real_cvprec) :: exner_ratio(cmpr%n_points)

! Compressed subsidence warming per unit mass-flux
real(kind=real_cvprec) :: delta_tv_cmpr(cmpr%n_points)

! Vertical interpolation weight
real(kind=real_cvprec) :: interp

! Array lower and upper bounds
integer :: lb(3)
integer :: ub(3)

! Loop counters
integer :: ic, i_field, ij


!------------------------------------------------------------------------------
! 1) Compress required fields
!------------------------------------------------------------------------------

! Compress layer-mass
lb = [1,1,1]
ub = [nx_full,ny_full,1]
call compress( cmpr, lb(1:2), ub(1:2), layer_mass(:,:,k),                      &
               layer_mass_k )

! Compress required fields from levels k
do i_field = i_temperature, i_qc_last
  lb = lbound( fields % list(i_field)%pt )
  ub = ubound( fields % list(i_field)%pt )
  call compress( cmpr, lb(1:2), ub(1:2), fields % list(i_field)%pt(:,:,k),     &
                 fields_k(:,i_field) )
  ! Remove any spurious negative values from input data
  if ( field_positive(i_field) ) then
    do ic = 1, cmpr%n_points
      fields_k(ic,i_field)   = max( fields_k(ic,i_field),   zero )
    end do
  end if
end do

! Compress heights, pressures and virtual temperatures
if ( l_kpdk ) then
  ! 2nd half-level step from k to next...

  ! Get height and pressure from k (_1) and next (_2)
  call grid_compress( grid, cmpr,                                              &
                      k = k,         grid_k_super = grid_1,                    &
                      k_half=k_next, grid_half_super=grid_2 )
  lb = lbound(grid%height_half)
  ub = ubound(grid%height_half)
  call compress( cmpr, lb(1:2), ub(1:2), grid%height_half(:,:,k_next-dk),      &
                 height_work )
  ! Set fraction of layer k on 2nd half-level-step
  do ic = 1, cmpr%n_points
    frac_level_step(ic)                                                        &
      = ( grid_2(ic,i_height) - grid_1(ic,i_height) )                          &
      / ( grid_2(ic,i_height) - height_work(ic) )
  end do
  ! Compress env Tv at k, k+dk and interpolate to next (_2)
  lb = [1,1,k_bot_conv]
  ub = [nx_full,ny_full,k_top_conv]
  call compress( cmpr, lb(1:2), ub(1:2), virt_temp(:,:,k),                     &
                 virt_temp_1 )
  if ( l_last_level ) then
    ! Can't interpolate at last level as fields don't exist at k+dk,
    ! so just assume constant Tv
    do ic = 1, cmpr%n_points
      virt_temp_2(ic) = virt_temp_1(ic)
    end do
  else
    call compress( cmpr, lb(1:2), ub(1:2), virt_temp(:,:,k+dk),                &
                   virt_temp_2 )
    lb = lbound(grid%height_full)
    ub = ubound(grid%height_full)
    call compress( cmpr, lb(1:2), ub(1:2), grid%height_full(:,:,k+dk),         &
                   height_work )
    do ic = 1, cmpr%n_points
      interp = ( grid_2(ic,i_height) - grid_1(ic,i_height) )                   &
             / ( height_work(ic)     - grid_1(ic,i_height) )
      virt_temp_2(ic) = (one-interp) * virt_temp_1(ic)                         &
                      +      interp  * virt_temp_2(ic)
    end do
  end if

else
  ! 1st half-level step from prev to k...

  ! Get height and pressure from prev (_1) and k (_2)
  call grid_compress( grid, cmpr,                                              &
                      k = k,            grid_k_super = grid_2,                 &
                      k_half=k_next-dk, grid_half_super=grid_1 )
  lb = lbound(grid%height_half)
  ub = ubound(grid%height_half)
  call compress( cmpr, lb(1:2), ub(1:2), grid%height_half(:,:,k_next),         &
                 height_work )
  ! Set fraction of layer k on 1st half-level-step
  do ic = 1, cmpr%n_points
    frac_level_step(ic)                                                        &
      = ( grid_2(ic,i_height) - grid_1(ic,i_height) )                          &
      / ( height_work(ic)     - grid_1(ic,i_height) )
  end do
  ! Compress env Tv at k, k-dk and interpolate to prev (_1)
  lb = [1,1,k_bot_conv]
  ub = [nx_full,ny_full,k_top_conv]
  call compress( cmpr, lb(1:2), ub(1:2), virt_temp(:,:,k),                     &
                 virt_temp_2 )
  if ( l_last_level ) then
    ! Setting Tv(prev) = Tv(k) at last level to preserve KGO,
    ! but this is wrong; we should just interpolate from k and k-dk as usual.
    do ic = 1, cmpr%n_points
      virt_temp_1(ic) = virt_temp_2(ic)
    end do
  else
    call compress( cmpr, lb(1:2), ub(1:2), virt_temp(:,:,k-dk),                &
                   virt_temp_1 )
    lb = lbound(grid%height_full)
    ub = ubound(grid%height_full)
    call compress( cmpr, lb(1:2), ub(1:2), grid%height_full(:,:,k-dk),         &
                   height_work )
    do ic = 1, cmpr%n_points
      interp = ( grid_1(ic,i_height) - grid_2(ic,i_height) )                   &
             / ( height_work(ic)     - grid_2(ic,i_height) )
      virt_temp_1(ic) = (one-interp) * virt_temp_2(ic)                         &
                      +      interp  * virt_temp_1(ic)
    end do
  end if

end if

! Compute mass on the current half-level step
do ic = 1, cmpr%n_points
  layer_mass_k(ic) = layer_mass_k(ic) * frac_level_step(ic)
  exner_ratio(ic) = one
end do

!------------------------------------------------------------------------------
! 2) Compute change in Tv when subsiding air from 2 to 1
!------------------------------------------------------------------------------

! Adiabatically adjust temperature from k+dk to k
call dry_adiabat( cmpr%n_points, cmpr%n_points,                                &
                  grid_2(:,i_pressure), grid_1(:,i_pressure),                  &
                  fields_k(:,i_q_vap),                                         &
                  fields_k(:,i_qc_first:i_qc_last),                            &
                  exner_ratio )

! Compute Tv difference over layer-mass ( 1/rho dTv/dz )
do ic = 1, cmpr%n_points
  delta_tv_cmpr(ic) = ( virt_temp_2(ic) * exner_ratio(ic) - virt_temp_1(ic) )  &
                    * comorph_timestep * alpha_detrain                         &
                    / layer_mass_k(ic)
end do


!------------------------------------------------------------------------------
! 3) Scatter calculated subsidence Tv increment back to full array
!------------------------------------------------------------------------------

do ic = 1, cmpr%n_points
  ij = cmpr%index_i(ic) + nx_full * ( cmpr%index_j(ic) - 1 )
  delta_tv(ij) = delta_tv_cmpr(ic)
end do


return

end subroutine calc_delta_tv

end module calc_delta_tv_mod
