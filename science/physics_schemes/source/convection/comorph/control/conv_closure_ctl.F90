! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module conv_closure_ctl_mod

implicit none

contains


! Subroutine to calculate and apply closure rescaling to the
! resolved-scale source terms (and diagnostics for consistency).
subroutine conv_closure_ctl( n_updraft_layers, n_dndraft_layers,               &
                             ij_first, ij_last,                                &
                             n_fields_tot, cmpr_any, layer_mass,               &
                             fields_np1,                                       &
                             updraft_fallback_par_gen,                         &
                             dndraft_fallback_par_gen,                         &
                             updraft_res_source,                               &
                             updraft_fallback_res_source,                      &
                             dndraft_res_source,                               &
                             dndraft_fallback_res_source,                      &
                             updraft_fields_2d,                                &
                             dndraft_fields_2d,                                &
                             comorph_diags,                                    &
                             updraft_diags_super,                              &
                             updraft_fallback_diags_super,                     &
                             dndraft_diags_super,                              &
                             dndraft_fallback_diags_super )

use comorph_constants_mod, only: real_cvprec, real_hmprec, zero, one,          &
                     nx_full, ny_full, k_bot_conv, k_top_conv,                 &
                     n_updraft_types, n_dndraft_types,                         &
                     l_updraft_fallback, l_dndraft_fallback,                   &
                     comorph_timestep, max_cfl, sqrt_min_float,                &
                     l_calc_mfw_cape

use cmpr_type_mod, only: cmpr_type
use parcel_type_mod, only: parcel_type
use res_source_type_mod, only: res_source_type
use comorph_diags_type_mod, only: comorph_diags_type
use draft_diags_type_mod, only: draft_diags_super_type,                        &
                                draft_diags_scaling
use fields_type_mod, only: fields_type, ragged_array_type,                     &
                           i_q_vap, i_qc_last
use fields_2d_mod, only: n_fields_2d
use cfl_limit_indep_mod, only: cfl_limit_indep,                                &
                               cfl_limit_indep_fallback
use cfl_limit_sum_ent_mod, only: cfl_limit_sum_ent,                            &
                                 cfl_limit_sum_ent_fallback
use cfl_limit_scaling_mod, only: cfl_limit_scaling
use apply_scaling_mod, only: apply_scaling
use closure_scale_integrals_mod, only: closure_scale_integrals

implicit none

! Highest number of distinct updraft and downdraft layers
! occuring on the current segment
integer, intent(in) :: n_updraft_layers
integer, intent(in) :: n_dndraft_layers

! Position indices ( nx*(j-1) + i ) of the first and last
! convecting point on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Total number of fields to loop over (depends on whether or
! not this call is updating the passive tracers)
integer, intent(in) :: n_fields_tot

! Structure storing compression list indices for the points where
! any kind of convection is occuring
type(cmpr_type), intent(in) :: cmpr_any( k_bot_conv:k_top_conv )

! Full 3-D field of layer-masses (used for finding CFL limit)
real(kind=real_hmprec), intent(in) :: layer_mass                               &
                     ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Structure containing pointers to the primary fields just before
! being updated by convection increments
type(fields_type), intent(in) :: fields_np1

! Initiation mass-sources for fallback flows.
! These are used for applying the CFL limit for the fall-back
! flows; the fall-back initiating mass-sources don't count
! towards the entrainment out of the layer, since they come
! from detrainment by primary drafts, not from the environment.
type(parcel_type), intent(in) :: updraft_fallback_par_gen(:,:,:)
type(parcel_type), intent(in) :: dndraft_fallback_par_gen(:,:,:)
! Will be dimensioned either
! ( n_updraft_types, n_updraft_layers, k_bot_conv:k_top_conv )
! or
! ( 1, 1, k_bot_conv:k_top_conv )
! depending on whether they are used.

! Arrays of structures containing resolved-scale source terms
! in separate compression lists for updrafts & downdrafts
! and their fall-back flows
type(res_source_type), intent(in out) ::                                       &
                       updraft_res_source(:,:,:)
type(res_source_type), intent(in out) ::                                       &
                       updraft_fallback_res_source(:,:,:)
type(res_source_type), intent(in out) ::                                       &
                       dndraft_res_source (:,:,:)
type(res_source_type), intent(in out) ::                                       &
                       dndraft_fallback_res_source(:,:,:)
! These are allocatable arrays, and we don't know for certain
! what shape they will have been allocated to.
! If n_updraft_layers > 0, shape is
! ( n_updraft_types, n_updraft_layers, k_bot_conv:k_top_conv )
! But if n_updraft_layers = 0, shape is
! ( 1, 1, k_bot_conv:k_top_conv )

! Super-arrays storing 2D work arrays
real(kind=real_cvprec), intent(in out) :: updraft_fields_2d                    &
                                          ( ij_first:ij_last, n_fields_2d,     &
                                            n_updraft_types, n_updraft_layers )
real(kind=real_cvprec), intent(in out) :: dndraft_fields_2d                    &
                                          ( ij_first:ij_last, n_fields_2d,     &
                                            n_dndraft_types, n_dndraft_layers )

! Structure containing switches and meta-data for diagnostics
type(comorph_diags_type), intent(in) :: comorph_diags

! Compressed super-arrays for diagnostics from each draft.
! Some of the contained diag fields may need to be rescaled
! to be consistent with the rescaling of the mass-flux
type(draft_diags_super_type), intent(in out) ::                                &
                              updraft_diags_super
type(draft_diags_super_type), intent(in out) ::                                &
                              updraft_fallback_diags_super
type(draft_diags_super_type), intent(in out) ::                                &
                              dndraft_diags_super
type(draft_diags_super_type), intent(in out) ::                                &
                              dndraft_fallback_diags_super


! Closure scalings for each draft / type / layer on ij-grid
real(kind=real_cvprec) :: updraft_scaling                                      &
         ( ij_first:ij_last, n_updraft_types, n_updraft_layers )
real(kind=real_cvprec) :: dndraft_scaling                                      &
         ( ij_first:ij_last, n_dndraft_types, n_dndraft_layers )
! Note: the fall-back flows must always use the same closure
! scaling as their parent primary updraft or downdraft

! Rescalings applied for CFL restriction
real(kind=real_cvprec) :: updraft_cfl_scaling                                  &
         ( ij_first:ij_last, n_updraft_types, n_updraft_layers )
real(kind=real_cvprec) :: dndraft_cfl_scaling                                  &
         ( ij_first:ij_last, n_dndraft_types, n_dndraft_layers )

! CFL rescaling on each level, on the cmpr_any compression list;
! This is the mass on the model-level divided by the
! sum of entrainment over all drafts / types / layers.
! used in the CFL limit calculation.
type(ragged_array_type), allocatable :: cfl_scaling(:)

! Sum of entrained masses on current level
! (used to calculate the above)
real(kind=real_cvprec) :: sum_ent_k(ij_first:ij_last)
! Sum of resolved-scale forcing of water species on current level
! (used in the check to avoid creating negative water)
real(kind=real_cvprec) :: sum_res_source_q_k                                   &
                          ( ij_first:ij_last, i_q_vap:i_qc_last )

! Compression arrays for layer-mass and water-species
real(kind=real_cvprec) :: layer_mass_cmpr(ij_first:ij_last)
real(kind=real_cvprec) :: q_cmpr

! Res source and layer terms compared against CFL limit
real(kind=real_cvprec) :: q_removed
real(kind=real_cvprec) :: q_available

! Grid to store indices of points in the cmpr_any compression
! list, for the current level
integer :: index_ic_k( ij_first:ij_last )

! Loop counters
integer :: ic, ij, i, j, k, i_type, i_layr, i_field


! Initialise all closure scalings to 1
if ( n_updraft_types > 0 .and. n_updraft_layers > 0 ) then
  do i_layr = 1, n_updraft_layers
    do i_type = 1, n_updraft_types
      do ij = ij_first, ij_last
        updraft_scaling(ij,i_type,i_layr) = one
      end do
    end do
  end do
end if
if ( n_dndraft_types > 0 .and. n_dndraft_layers > 0 ) then
  do i_layr = 1, n_dndraft_layers
    do i_type = 1, n_dndraft_types
      do ij = ij_first, ij_last
        dndraft_scaling(ij,i_type,i_layr) = one
      end do
    end do
  end do
end if


!----------------------------------------------------------------
! 1) Calculate some physically-based closure scalings
!----------------------------------------------------------------

! (not yet implemented)



!----------------------------------------------------------------
! 2) Apply CFL limits to closure scalings.
!----------------------------------------------------------------

! This is essentially just enforcing the restriction that
! convection cannot remove more mass from a model-level in one
! timestep than exists on that level.

! Note: sometimes we may have multiple vertically-separated
! layers of convection.  If we reduce the closure scaling to
! meet the CFL limit for one layer, we don't want to needlessly
! scale down other layers in the same column.  So a simple
! "scale down the whole column" approach is not desirable.
! We proceed by first scaling down each draft / layer / type
! of convection independently such that it would meet the CFL
! limit if it were the only convection present.  We then sum
! the total-mass entrained on each model-level over all
! drafts / layers / types.  Each draft / layer / type is
! then further rescaled so as to satisfy the summed CFL
! restriction on all model-levels where it is active.


! First, impose CFL limit on closure scaling for each convective
! draft of each type / layer independently:

! For updrafts (and their fall-backs if used)
if ( n_updraft_types > 0 .and. n_updraft_layers > 0 ) then
  if ( l_updraft_fallback ) then
    call cfl_limit_indep_fallback( n_updraft_types,                            &
                                   n_updraft_layers,                           &
                                   ij_first, ij_last, cmpr_any,                &
                                   layer_mass, fields_np1,                     &
                                   updraft_res_source,                         &
                                   updraft_fallback_res_source,                &
                                   updraft_fallback_par_gen,                   &
                                   updraft_scaling )
  else
    call cfl_limit_indep( n_updraft_types, n_updraft_layers,                   &
                          ij_first, ij_last, cmpr_any,                         &
                          layer_mass, fields_np1,                              &
                          updraft_res_source,                                  &
                          updraft_scaling )
  end if
end if

! For downdrafts (and their fall-backs if used)
if ( n_dndraft_types > 0 .and. n_dndraft_layers > 0 ) then
  if ( l_dndraft_fallback ) then
    call cfl_limit_indep_fallback( n_dndraft_types,                            &
                                   n_dndraft_layers,                           &
                                   ij_first, ij_last, cmpr_any,                &
                                   layer_mass, fields_np1,                     &
                                   dndraft_res_source,                         &
                                   dndraft_fallback_res_source,                &
                                   dndraft_fallback_par_gen,                   &
                                   dndraft_scaling )
  else
    call cfl_limit_indep( n_dndraft_types, n_dndraft_layers,                   &
                          ij_first, ij_last, cmpr_any,                         &
                          layer_mass, fields_np1,                              &
                          dndraft_res_source,                                  &
                          dndraft_scaling )
  end if
end if


! Next, impose CFL limit on closure scaling due to sum of
! entrainment over all drafts / types / layers

! Allocate space for CFL-limit scaling on each level
allocate( cfl_scaling( k_bot_conv:k_top_conv ) )
do k = k_bot_conv, k_top_conv
  if ( cmpr_any(k) % n_points > 0 )                                            &
     allocate( cfl_scaling(k) % f( cmpr_any(k) % n_points ) )
end do

! Initialise sum of all entrainment to zero
do ij = ij_first, ij_last
  sum_ent_k(ij) = zero
end do
! Initialise sum of resolved-scale source terms for water species
do i_field = i_q_vap, i_qc_last
  do ij = ij_first, ij_last
    sum_res_source_q_k(ij,i_field) = zero
  end do
end do


! Loop over levels to sum total entrainment on each level:
do k = k_bot_conv, k_top_conv
  if ( cmpr_any(k) % n_points > 0 ) then

    if ( n_updraft_types > 0 .and. n_updraft_layers > 0 ) then
      ! Add contributions from primary updrafts
      call cfl_limit_sum_ent( n_updraft_types, n_updraft_layers,               &
                              ij_first, ij_last, updraft_scaling,              &
                              updraft_res_source(:,:,k),                       &
                              sum_ent_k, sum_res_source_q_k )
      ! Add contributions from updraft fall-backs
      if ( l_updraft_fallback ) then
        call cfl_limit_sum_ent_fallback(                                       &
                              n_updraft_types, n_updraft_layers,               &
                              ij_first, ij_last, updraft_scaling,              &
                              updraft_fallback_res_source(:,:,k),              &
                              updraft_fallback_par_gen(:,:,k),                 &
                              sum_ent_k, sum_res_source_q_k )
      end if
    end if
    if ( n_dndraft_types > 0 .and. n_dndraft_layers > 0 ) then
      ! Add contributions from primary downdrafts
      call cfl_limit_sum_ent( n_dndraft_types, n_dndraft_layers,               &
                              ij_first, ij_last, dndraft_scaling,              &
                              dndraft_res_source(:,:,k),                       &
                              sum_ent_k, sum_res_source_q_k )
      ! Add contributions from downdraft fall-backs
      if ( l_dndraft_fallback ) then
        call cfl_limit_sum_ent_fallback(                                       &
                              n_dndraft_types, n_dndraft_layers,               &
                              ij_first, ij_last, dndraft_scaling,              &
                              dndraft_fallback_res_source(:,:,k),              &
                              dndraft_fallback_par_gen(:,:,k),                 &
                              sum_ent_k, sum_res_source_q_k )
      end if
    end if

    ! Compress layer-mass onto convecting points
    do ic = 1, cmpr_any(k) % n_points
      i = cmpr_any(k) % index_i(ic)
      j = cmpr_any(k) % index_j(ic)
      ij = nx_full*(j-1)+i
      layer_mass_cmpr(ij) = real( layer_mass(i,j,k), real_cvprec)
    end do

    ! Calculate a maximum allowable closure rescaling
    ! (relative to the scalings already calculated)
    ! for each point
    do ic = 1, cmpr_any(k) % n_points
      i = cmpr_any(k) % index_i(ic)
      j = cmpr_any(k) % index_j(ic)
      ij = nx_full*(j-1)+i
      cfl_scaling(k) % f(ic) = min( one,                                       &
         max_cfl * layer_mass_cmpr(ij)                                         &
         / ( comorph_timestep * max( sum_ent_k(ij),                            &
                                     sqrt_min_float ) )  )
          ! Check avoids div-by-zero where no entrainment
      ! Reset to zero ready for next level
      sum_ent_k(ij) = zero
    end do
    ! Also apply limit to avoid creating negative water species
    do i_field = i_q_vap, i_qc_last
      do ic = 1, cmpr_any(k) % n_points
        i = cmpr_any(k) % index_i(ic)
        j = cmpr_any(k) % index_j(ic)
        ij = nx_full*(j-1)+i
        ! If there is a net sink of this water species
        if ( sum_res_source_q_k(ij,i_field) < zero ) then
          q_cmpr = max( real( fields_np1%list(i_field)%pt(i,j,k),              &
                              real_cvprec ), zero )

          ! Calculate water-mass removed and water-mass available
          q_removed = -comorph_timestep * sum_res_source_q_k(ij,i_field)
          q_available = max_cfl * q_cmpr * layer_mass_cmpr(ij)

          ! Apply CFL limit if removal term exceeds available term
          if ( q_removed > q_available ) then
            cfl_scaling(k) % f(ic) = min( cfl_scaling(k) % f(ic),              &
                                          q_available / q_removed )
          end if

        end if
        ! Reset to zero ready for next level
        sum_res_source_q_k(ij,i_field) = zero
      end do
    end do

  end if  ! ( cmpr_any(k) % n_points > 0 )
end do  ! k = k_bot_conv, k_top_conv


! Initialise CFL-limit rescalings
if ( n_updraft_types > 0 .and. n_updraft_layers > 0 ) then
  do i_layr = 1, n_updraft_layers
    do i_type = 1, n_updraft_types
      do ij = ij_first, ij_last
        updraft_cfl_scaling(ij,i_type,i_layr) = one
      end do
    end do
  end do
end if
if ( n_dndraft_types > 0 .and. n_dndraft_layers > 0 ) then
  do i_layr = 1, n_dndraft_layers
    do i_type = 1, n_dndraft_types
      do ij = ij_first, ij_last
        dndraft_cfl_scaling(ij,i_type,i_layr) = one
      end do
    end do
  end do
end if

! Loop over levels to find minimum of cfl_scaling over all
! model-levels where each convective draft / type / layer is
! active.  This sets the CFL rescaling to apply to each.
do k = k_bot_conv, k_top_conv
  if ( cmpr_any(k) % n_points > 0 ) then

    ! Scatter cmpr_any compression indices for current level
    ! into the grid
    do ic = 1, cmpr_any(k) % n_points
      i = cmpr_any(k) % index_i(ic)
      j = cmpr_any(k) % index_j(ic)
      index_ic_k( nx_full*(j-1)+i ) = ic
    end do

    if ( n_updraft_types > 0 .and. n_updraft_layers > 0 ) then
      ! Minimum cfl_scaling over each primary updraft
      call cfl_limit_scaling( n_updraft_types, n_updraft_layers,               &
                              cmpr_any(k) % n_points,                          &
                              cfl_scaling(k) % f,                              &
                              updraft_res_source(:,:,k),                       &
                              ij_first, ij_last, index_ic_k,                   &
                              updraft_cfl_scaling )
      ! Minimum cfl_scaling over each updraft fall-back
      if ( l_updraft_fallback ) then
        call cfl_limit_scaling(                                                &
                              n_updraft_types, n_updraft_layers,               &
                              cmpr_any(k) % n_points,                          &
                              cfl_scaling(k) % f,                              &
                              updraft_fallback_res_source(:,:,k),              &
                              ij_first, ij_last, index_ic_k,                   &
                              updraft_cfl_scaling )
      end if
    end if
    if ( n_dndraft_types > 0 .and. n_dndraft_layers > 0 ) then
      ! Minimum cfl_scaling over each primary downdraft
      call cfl_limit_scaling( n_dndraft_types, n_dndraft_layers,               &
                              cmpr_any(k) % n_points,                          &
                              cfl_scaling(k) % f,                              &
                              dndraft_res_source(:,:,k),                       &
                              ij_first, ij_last, index_ic_k,                   &
                              dndraft_cfl_scaling )
      ! Minimum cfl_scaling over each downdraft fall-back
      if ( l_dndraft_fallback ) then
        call cfl_limit_scaling(                                                &
                              n_dndraft_types, n_dndraft_layers,               &
                              cmpr_any(k) % n_points,                          &
                              cfl_scaling(k) % f,                              &
                              dndraft_fallback_res_source(:,:,k),              &
                              ij_first, ij_last, index_ic_k,                   &
                              dndraft_cfl_scaling )
      end if
    end if

  end if  ! ( cmpr_any(k) % n_points > 0 )
end do  ! k = k_bot_conv, k_top_conv

! Deallocate cfl_scaling ragged array
deallocate( cfl_scaling )

! Apply the CFL rescalings to the closure scalings
if ( n_updraft_types > 0 .and. n_updraft_layers > 0 ) then
  do i_layr = 1, n_updraft_layers
    do i_type = 1, n_updraft_types
      do ij = ij_first, ij_last
        updraft_scaling(ij,i_type,i_layr)                                      &
          = updraft_scaling(ij,i_type,i_layr)                                  &
          * updraft_cfl_scaling(ij,i_type,i_layr)
      end do
    end do
  end do
end if
if ( n_dndraft_types > 0 .and. n_dndraft_layers > 0 ) then
  do i_layr = 1, n_dndraft_layers
    do i_type = 1, n_dndraft_types
      do ij = ij_first, ij_last
        dndraft_scaling(ij,i_type,i_layr)                                      &
          = dndraft_scaling(ij,i_type,i_layr)                                  &
          * dndraft_cfl_scaling(ij,i_type,i_layr)
      end do
    end do
  end do
end if



!----------------------------------------------------------------
! 3) Apply the final closure scalings to the resolved-scale
!    source terms and diagnostics
!----------------------------------------------------------------

! Loop over levels to apply closure scaling to the
! resolved-scale source terms.
do k = k_bot_conv, k_top_conv
  if ( cmpr_any(k) % n_points > 0 ) then

    if ( n_updraft_types > 0 .and. n_updraft_layers > 0 ) then
      ! Primary updrafts
      call apply_scaling( n_updraft_types, n_updraft_layers,                   &
                          n_fields_tot, cmpr_any(k) % n_points,                &
                          ij_first, ij_last, updraft_scaling,                  &
                          updraft_res_source(:,:,k) )
      ! Updraft fall-backs
      if ( l_updraft_fallback ) then
        call apply_scaling( n_updraft_types, n_updraft_layers,                 &
                            n_fields_tot, cmpr_any(k) % n_points,              &
                            ij_first, ij_last, updraft_scaling,                &
                            updraft_fallback_res_source(:,:,k) )
      end if
    end if
    if ( n_dndraft_types > 0 .and. n_dndraft_layers > 0 ) then
      ! Primary downdrafts
      call apply_scaling( n_dndraft_types, n_dndraft_layers,                   &
                          n_fields_tot, cmpr_any(k) % n_points,                &
                          ij_first, ij_last, dndraft_scaling,                  &
                          dndraft_res_source(:,:,k) )
      ! Downdraft fall-backs
      if ( l_dndraft_fallback ) then
        call apply_scaling( n_dndraft_types, n_dndraft_layers,                 &
                            n_fields_tot, cmpr_any(k) % n_points,              &
                            ij_first, ij_last, dndraft_scaling,                &
                            dndraft_fallback_res_source(:,:,k) )
      end if
    end if

  end if  ! ( cmpr_any(k) % n_points > 0 )
end do  ! k = k_bot_conv, k_top_conv


! Apply closure scalings to vertical integrals used to calculate
! mass-flux-weighted CAPE
if ( l_calc_mfw_cape ) then

  if ( n_updraft_types > 0 .and. n_updraft_layers > 0 ) then
    call closure_scale_integrals( ij_first, ij_last,                           &
                                  n_updraft_types, n_updraft_layers,           &
                                  updraft_scaling, updraft_fields_2d )
  end if

  if ( n_dndraft_types > 0 .and. n_dndraft_layers > 0 ) then
    call closure_scale_integrals( ij_first, ij_last,                           &
                                  n_dndraft_types, n_dndraft_layers,           &
                                  dndraft_scaling, dndraft_fields_2d )
  end if

end if


! Apply closure scalings to diagnostics
if ( n_updraft_types > 0 .and. n_updraft_layers > 0 ) then
  ! Primary updrafts
  if ( comorph_diags % updraft % n_diags > 0 ) then
    call draft_diags_scaling( n_updraft_types, n_updraft_layers,               &
                              ij_first, ij_last, updraft_scaling,              &
                              comorph_diags % updraft,                         &
                              updraft_diags_super )
  end if
  ! Updraft fall-backs
  if ( l_updraft_fallback ) then
    if ( comorph_diags % updraft_fallback % n_diags > 0 ) then
      call draft_diags_scaling( n_updraft_types,n_updraft_layers,              &
                                ij_first,ij_last,updraft_scaling,              &
                                comorph_diags % updraft_fallback,              &
                                updraft_fallback_diags_super )
    end if
  end if
end if
if ( n_dndraft_types > 0 .and. n_dndraft_layers > 0 ) then
  ! Primary downdrafts
  if ( comorph_diags % dndraft % n_diags > 0 ) then
    call draft_diags_scaling( n_dndraft_types, n_dndraft_layers,               &
                              ij_first, ij_last, dndraft_scaling,              &
                              comorph_diags % dndraft,                         &
                              dndraft_diags_super )
  end if
  ! Downdraft fall-backs
  if ( l_dndraft_fallback ) then
    if ( comorph_diags % dndraft_fallback % n_diags > 0 ) then
      call draft_diags_scaling( n_dndraft_types,n_dndraft_layers,              &
                                ij_first,ij_last,dndraft_scaling,              &
                                comorph_diags % dndraft_fallback,              &
                                dndraft_fallback_diags_super )
    end if
  end if
end if


return
end subroutine conv_closure_ctl


end module conv_closure_ctl_mod
