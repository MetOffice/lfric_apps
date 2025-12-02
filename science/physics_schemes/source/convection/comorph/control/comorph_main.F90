! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module comorph_main_mod

implicit none

contains

! main comorph routine; calls various other routines to diagnose
! initiation mass-sources, perform updraft and downdraft parcel sweeps,
! and update primary fields with convective increments
! (all only at points on the current segment).
subroutine comorph_main( max_points, ij_first, ij_last,                        &
                         n_fields_tot, l_tracer,                               &
                         grid, turb, cloudfracs, fields,                       &
                         layer_mass, virt_temp_n, virt_temp_np1,               &
                         l_init_poss, comorph_diags )

use comorph_constants_mod, only: real_hmprec, real_cvprec,                     &
                                 nx_full, ny_full,                             &
                                 k_bot_conv, k_top_conv, k_top_init,           &
                                 n_updraft_types, n_dndraft_types,             &
                                 l_updraft_fallback, l_dndraft_fallback,       &
                                 l_calc_mfw_cape,                              &
                                 i_cfl_check, i_cfl_check_closure,             &
                                 i_cfl_check_hybrid

use cmpr_type_mod, only: cmpr_type
use grid_type_mod, only: grid_type
use turb_type_mod, only: turb_type
use cloudfracs_type_mod, only: cloudfracs_type
use fields_type_mod, only: fields_type
use parcel_type_mod, only: parcel_type, parcel_dealloc
use res_source_type_mod, only: res_source_type
use comorph_diags_type_mod, only: comorph_diags_type
use draft_diags_type_mod, only: draft_diags_super_type,                        &
                                draft_diags_compute_means
use diags_2d_type_mod, only: diags_2d_compute_means

use conv_genesis_ctl_mod, only: conv_genesis_ctl
use conv_sweep_ctl_mod, only: conv_sweep_ctl
use normalise_integrals_mod, only: normalise_integrals
use find_cmpr_any_mod, only: find_cmpr_any
use conv_closure_ctl_mod, only: conv_closure_ctl
use conv_incr_ctl_mod, only: conv_incr_ctl
use calc_diag_conv_cloud_mod, only: zero_diag_conv_cloud

implicit none

! Max size the compression lists can possibly need
! (i.e. the number of columns in the current segment)
integer, intent(in) :: max_points

! i,j indices of the first and last point on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Total number of fields, including tracers if they are used
integer, intent(in) :: n_fields_tot

! Flag for whether tracers are used
logical, intent(in) :: l_tracer

! Structure containing pointers to model grid and related
! fields; contains model-level heights, pressures and
! dry-density
type(grid_type), intent(in) :: grid

! Structure containing pointers to turbulence fields passed
! in from the boundary-layer scheme
type(turb_type), intent(in) :: turb

! Structure containing cloud fractions, in the event that they
! are diagnostic rather than prognostic and so aren't included
! in the primary fields structure.
type(cloudfracs_type), intent(in out) :: cloudfracs

! Structure containing pointers to the primary fields
type(fields_type), intent(in out) :: fields

! Dry-mass on each model-level, per unit surface area
real(kind=real_hmprec), intent(in) :: layer_mass                               &
                ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Environment virtual temperature profile (full 3-D array):
! Start-of-timestep
real(kind=real_hmprec), intent(in) :: virt_temp_n                              &
                ( nx_full, ny_full, k_bot_conv:k_top_conv )
! Latest fields
real(kind=real_hmprec), intent(in) :: virt_temp_np1                            &
                ( nx_full, ny_full, k_bot_conv:k_top_conv )

! 3-D mask of points where convective initiation mass-sources
! for updrafts or downdrafts might be possible
logical, intent(in) :: l_init_poss                                             &
                ( nx_full, ny_full, k_bot_conv:k_top_init )

! Structure containing pointers to output diagnostics and
! miscellaneous outputs required elsewhere in the model
type(comorph_diags_type), intent(in out) :: comorph_diags
! (holds outputs, but intent inout to retain the assigned status
!  of the contained pointers on input)


! Array of structures storing compression list indices
! for the combined list of points where any kind of convection
! is occuring on each level
! (used by the resolved-scale increment calculation in
!  conv_incr_ctl).
type(cmpr_type), allocatable :: cmpr_any(:)


! The following data are stored as arrays of derived-type
! structures. The structures contain separate compression
! lists on each model-level. i.e. the compression lists maybe
! of different lengths on each level, so that they form
! ragged arrays.
! Storing the data in this way has 2 advantages:
! a) It consumes much less memory.  No data has to be stored
!    for levels/points where there is no convection!
! b) The data don't have to be recompressed onto convecting
!    points when used.

! Arrays of structures for updraft and downdraft initiation
! mass sources, for each convection type,
! each vertically distinct convecting layer, and each model-level
! (hence 3-D)
type(parcel_type), allocatable :: updraft_par_gen(:,:,:)
type(parcel_type), allocatable :: dndraft_par_gen(:,:,:)
! These arrays store initiating parcel properties on
! thermodynamic-levels (the levels they are entrained on).
! To convert them to initial parcel properties on flux-levels,
! the pressure (and temperature) need to be adjusted
! consistently, and resulting water phase-changes performed.
! This parcel initialisation is done in init_par_conv,
! which is called from conv_sweep_ctl.

! Initiating mass source properties for the fall-back flows
! (overshooting detrained air that needs to fall back to its
!  neutral buoyancy level in a separate downdraft calculation)
type(parcel_type), allocatable :: updraft_fallback_par_gen(:,:,:)
type(parcel_type), allocatable :: dndraft_fallback_par_gen(:,:,:)

! Arrays of structures to store resolved-scale source terms
! from convection
! (rates of entrainment and detrainment of dry-mass,
!  tendencies of the primary fields due to removal / injection
!  of convective parcels with different properties to the
!  grid-mean, and tendencies due to non-mass-exchanging forces
!  and precipitation)...

! ...for primary updrafts and downdrafts
type(res_source_type), allocatable ::                                          &
                       updraft_res_source(:,:,:)
type(res_source_type), allocatable ::                                          &
                       dndraft_res_source(:,:,:)

! ...for fall-back flows from the updrafts and downdrafts
type(res_source_type), allocatable ::                                          &
                       updraft_fallback_res_source(:,:,:)
type(res_source_type), allocatable ::                                          &
                       dndraft_fallback_res_source(:,:,:)

! Super-arrays to hold 2D work arrays (i.e. anything that isn't needed
! separately on each model-level, e.g. vertical integrals)
real(kind=real_cvprec), allocatable :: updraft_fields_2d(:,:,:,:)
real(kind=real_cvprec), allocatable :: dndraft_fields_2d(:,:,:,:)

! Max number of vertically distinct convection layers found
! within any one column in the current segment
! (set to zero if no convecting mass sources found at all)
integer :: n_updraft_layers
integer :: n_dndraft_layers


! Structures to store diagnostics from the updraft
! and downdraft calculations in compressed super-arrays
type(draft_diags_super_type) :: updraft_diags_super
type(draft_diags_super_type) :: dndraft_diags_super
type(draft_diags_super_type) :: updraft_fallback_diags_super
type(draft_diags_super_type) :: dndraft_fallback_diags_super
! Note: these need to be stored outside conv_sweep_ctl, as the values
! actually output for the diags may want to be the mass-flux
! weighted means over convection types / layers, but we won't
! know the mass-fluxes to use for the weighting until after the
! closure scaling has been applied.

! Flags input to the updraft / downdraft sweep routine
logical :: l_down
logical :: l_fallback
logical :: l_output_fallback

! Lower and upper bounds of pressure array, passed into conv_incr_ctl
integer :: lb_p(3), ub_p(3)

! Loop counters
integer :: k, i_type, i_layr


!--------------------------------------------------------------
! 1) Calculate initiation mass sources from each model-level
!--------------------------------------------------------------

! Call routine for initiation mass source calculation,
! for both updrafts and downdrafts
! Note: this routine calculates convective initiation based
! on the instability of the profile as measure by the
! gradient of the start-of-timestep virtual temperature,
! since that should be a well-balanced profile.
! But the initiating parcel properties are derived from
! the latest fields with other increments added on before
! convection, for consistency / better numerical behaviour
! when resolved-scale increments are computed from the
! difference of the parcel fields from the latest environment.
call conv_genesis_ctl( max_points, ij_first, ij_last,                          &
                       n_fields_tot, l_tracer,                                 &
                       grid, turb, cloudfracs, fields,                         &
                       layer_mass, virt_temp_n, l_init_poss,                   &
                       n_updraft_layers, updraft_par_gen,                      &
                       n_dndraft_layers, dndraft_par_gen )
! Note: the _par_gen arrays get allocated inside
! conv_genesis_ctl.


! Only anything else to do on the current segment if at least
! one convective mass-source layer has been found there...
if ( n_updraft_layers > 0 .or. n_dndraft_layers > 0 ) then


  !------------------------------------------------------------
  ! 2) Perform the updraft parcel ascent calculations
  !------------------------------------------------------------

  ! If updrafts are on,
  ! and at least one updraft mass-source layer has been found:
  if ( n_updraft_types > 0 .and. n_updraft_layers > 0 ) then

    ! Call routine to perform the primary updraft sweep.
    ! Computes convective mass-flux / transport and
    ! entrainment/detrainment sources for the primary fields
    l_down = .false.
    l_fallback = .false.
    l_output_fallback = l_updraft_fallback
    call conv_sweep_ctl( n_fields_tot,                                         &
                         n_updraft_types, n_updraft_layers,                    &
                         max_points, ij_first, ij_last,                        &
                         l_tracer, l_down, l_fallback,                         &
                         l_output_fallback,                                    &
                         grid, layer_mass, turb,                               &
                         fields, virt_temp_np1,                                &
                         updraft_par_gen,                                      &
                         updraft_res_source,                                   &
                         updraft_fields_2d,                                    &
                         comorph_diags % updraft,                              &
                         updraft_diags_super,                                  &
      fallback_par_gen = updraft_fallback_par_gen )

    ! If updraft fall-backs are on:
    if ( l_updraft_fallback ) then

      ! Call routine to perform the updraft fall-back
      ! (downwards) sweep.
      l_down = .true.
      l_fallback = .true.
      l_output_fallback = .false.
      call conv_sweep_ctl( n_fields_tot,                                       &
                           n_updraft_types, n_updraft_layers,                  &
                           max_points, ij_first, ij_last,                      &
                           l_tracer, l_down, l_fallback,                       &
                           l_output_fallback,                                  &
                           grid, layer_mass, turb,                             &
                           fields, virt_temp_np1,                              &
                           updraft_fallback_par_gen,                           &
                           updraft_fallback_res_source,                        &
                           updraft_fields_2d,                                  &
                           comorph_diags % updraft_fallback,                   &
                           updraft_fallback_diags_super )

    end if  ! ( l_updraft_fallback )

    if ( l_calc_mfw_cape ) then
      ! Apply any final normalisations needed to 2D output quantities
      call normalise_integrals( ij_first, ij_last,                             &
                                n_updraft_types, n_updraft_layers,             &
                                updraft_fields_2d )
    end if

  end if  ! ( n_updraft_types > 0 .and. n_updraft_layers > 0 )

  ! Allocate resolved-scale source term and fall-back
  ! mass-source arrays to minimal size if not used
  if ( .not. allocated(updraft_res_source) ) then
    allocate( updraft_res_source(1,1,k_bot_conv:k_bot_conv) )
  end if
  if ( .not. allocated(updraft_fallback_res_source) ) then
    allocate( updraft_fallback_res_source(1,1,k_bot_conv:k_bot_conv) )
  end if
  if ( .not. allocated(updraft_fallback_par_gen) ) then
    allocate( updraft_fallback_par_gen(1,1,k_bot_conv:k_bot_conv) )
  end if
  if ( .not. allocated(updraft_fields_2d) ) then
    allocate( updraft_fields_2d(1,1,1,1) )
  end if


  !------------------------------------------------------------
  ! 3) Perform the downdraft parcel descent calculations
  !------------------------------------------------------------

  ! If downdrafts are on,
  ! and at least one downdraft mass-source layer has been found:
  if ( n_dndraft_types > 0 .and. n_dndraft_layers > 0 ) then

    ! Call routine to perform the primary downdraft sweep.
    ! Computes convective mass-flux / transport and
    ! entrainment/detrainment sources for the primary fields
    l_down = .true.
    l_fallback = .false.
    l_output_fallback = l_dndraft_fallback
    call conv_sweep_ctl( n_fields_tot,                                         &
                         n_dndraft_types, n_dndraft_layers,                    &
                         max_points, ij_first, ij_last,                        &
                         l_tracer, l_down, l_fallback,                         &
                         l_output_fallback,                                    &
                         grid, layer_mass, turb,                               &
                         fields, virt_temp_np1,                                &
                         dndraft_par_gen,                                      &
                         dndraft_res_source,                                   &
                         dndraft_fields_2d,                                    &
                         comorph_diags % dndraft,                              &
                         dndraft_diags_super,                                  &
      fallback_par_gen = dndraft_fallback_par_gen )

    ! If downdraft fall-backs are on:
    if ( l_dndraft_fallback ) then

      ! Call routine to perform the downdraft fall-back
      ! (upwards) sweep.
      l_down = .false.
      l_fallback = .true.
      l_output_fallback = .false.
      call conv_sweep_ctl( n_fields_tot,                                       &
                           n_dndraft_types, n_dndraft_layers,                  &
                           max_points, ij_first, ij_last,                      &
                           l_tracer, l_down, l_fallback,                       &
                           l_output_fallback,                                  &
                           grid, layer_mass, turb,                             &
                           fields, virt_temp_np1,                              &
                           dndraft_fallback_par_gen,                           &
                           dndraft_fallback_res_source,                        &
                           dndraft_fields_2d,                                  &
                           comorph_diags % dndraft_fallback,                   &
                           dndraft_fallback_diags_super )

    end if  ! ( l_dndraft_fallback )

    if ( l_calc_mfw_cape ) then
      ! Apply any final normalisations needed to 2D output quantities
      call normalise_integrals( ij_first, ij_last,                             &
                                n_dndraft_types, n_dndraft_layers,             &
                                dndraft_fields_2d )
    end if

  end if  ! ( n_dndraft_types > 0 .and. n_dndraft_layers > 0 )

  ! Allocate resolved-scale source term and fall-back
  ! mass-source arrays to minimal size if not used
  if ( .not. allocated(dndraft_res_source) ) then
    allocate( dndraft_res_source(1,1,k_bot_conv:k_bot_conv) )
  end if
  if ( .not. allocated(dndraft_fallback_res_source) ) then
    allocate( dndraft_fallback_res_source(1,1,k_bot_conv:k_bot_conv) )
  end if
  if ( .not. allocated(dndraft_fallback_par_gen) ) then
    allocate( dndraft_fallback_par_gen(1,1,k_bot_conv:k_bot_conv) )
  end if
  if ( .not. allocated(dndraft_fields_2d) ) then
    allocate( dndraft_fields_2d(1,1,1,1) )
  end if


  !------------------------------------------------------------
  ! 4) Find the compression list of points on each level
  !    where any type of convection is active...
  !------------------------------------------------------------

  allocate( cmpr_any (k_bot_conv:k_top_conv) )
  ! This routine allocates arrays on each level to store
  ! compression indices for points where any kind of
  ! convection is active, and saves the i,j indices of
  ! those points.
  call find_cmpr_any( n_updraft_layers, n_dndraft_layers, ij_first, ij_last,   &
                      updraft_res_source,                                      &
                      updraft_fallback_res_source,                             &
                      dndraft_res_source,                                      &
                      dndraft_fallback_res_source,                             &
                      cmpr_any )


  !------------------------------------------------------------
  ! 5) Convective closure
  !    (rescales the fields in the res_source structures)
  !------------------------------------------------------------

  ! Currently the only closure available is a CFL-limiter
  if ( i_cfl_check == i_cfl_check_closure .or.                                 &
       i_cfl_check == i_cfl_check_hybrid ) then
    call conv_closure_ctl( n_updraft_layers, n_dndraft_layers,                 &
                           ij_first, ij_last,                                  &
                           n_fields_tot, cmpr_any, layer_mass,                 &
                           fields,                                             &
                           updraft_fallback_par_gen,                           &
                           dndraft_fallback_par_gen,                           &
                           updraft_res_source,                                 &
                           updraft_fallback_res_source,                        &
                           dndraft_res_source,                                 &
                           dndraft_fallback_res_source,                        &
                           updraft_fields_2d,                                  &
                           dndraft_fields_2d,                                  &
                           comorph_diags,                                      &
                           updraft_diags_super,                                &
                           updraft_fallback_diags_super,                       &
                           dndraft_diags_super,                                &
                           dndraft_fallback_diags_super )
  end if


  !------------------------------------------------------------
  ! 6) Calculate the updated primary model fields after
  !    convection ("large-scale increment" calculation)
  !------------------------------------------------------------

  ! Call routine to update the fields with increments due to:
  ! - entrainment and detrainment
  ! - non-mass-exchanging forces represented in the plume-model
  ! - optionally, compensating subsidence
  lb_p = lbound(grid%pressure_full)
  ub_p = ubound(grid%pressure_full)
  call conv_incr_ctl( n_updraft_layers, n_dndraft_layers,                      &
                      max_points, ij_first, ij_last,                           &
                      n_fields_tot, cmpr_any,                                  &
                      updraft_res_source,                                      &
                      updraft_fallback_res_source,                             &
                      dndraft_res_source,                                      &
                      dndraft_fallback_res_source,                             &
                      lb_p, ub_p, grid%pressure_full,                          &
                      layer_mass, fields, cloudfracs,                          &
                      comorph_diags )


  !------------------------------------------------------------
  ! 7) Compute means of draft diagnostics over types / layers
  !    and decompress into output arrays
  !------------------------------------------------------------

  ! The calls to draft_diags_compute_means below also deallocate the
  ! compressed diagnostic arrays, so are called in reverse
  ! order compared to the conv_sweep_ctl calls which allocated
  ! the arrays.

  if ( n_dndraft_types > 0 .and. n_dndraft_layers > 0 ) then
    if ( l_dndraft_fallback ) then
      ! Downdraft fall-back diagnostics:
      call draft_diags_compute_means( n_dndraft_types,n_dndraft_layers,        &
                                      max_points, ij_first, ij_last,           &
                                      comorph_diags%dndraft_fallback,          &
                                      dndraft_fallback_diags_super )
    end if
    ! Primary downdraft diagnostics
    call draft_diags_compute_means( n_dndraft_types, n_dndraft_layers,         &
                                    max_points, ij_first, ij_last,             &
                                    comorph_diags%dndraft,                     &
                                    dndraft_diags_super )
    ! 2D diagnostics
    if ( comorph_diags % dndraft_diags_2d % n_diags > 0 ) then
      call diags_2d_compute_means( n_dndraft_types, n_dndraft_layers,          &
                                   ij_first, ij_last,                          &
                                   comorph_diags % dndraft_diags_2d,           &
                                   dndraft_fields_2d )
    end if
  end if

  if ( n_updraft_types > 0 .and. n_updraft_layers > 0 ) then
    if ( l_updraft_fallback ) then
      ! Updraft fall-back diagnostics:
      call draft_diags_compute_means( n_updraft_types,n_updraft_layers,        &
                                      max_points, ij_first, ij_last,           &
                                      comorph_diags%updraft_fallback,          &
                                      updraft_fallback_diags_super )
    end if
    ! Primary updraft diagnostics
    call draft_diags_compute_means( n_updraft_types, n_updraft_layers,         &
                                    max_points, ij_first, ij_last,             &
                                    comorph_diags%updraft,                     &
                                    updraft_diags_super )
    ! 2D diagnostics
    if ( comorph_diags % updraft_diags_2d % n_diags > 0 ) then
      call diags_2d_compute_means( n_updraft_types, n_updraft_layers,          &
                                   ij_first, ij_last,                          &
                                   comorph_diags % updraft_diags_2d,           &
                                   updraft_fields_2d )
    end if
  end if


  ! Deallocate cmpr_any compression indices
  deallocate( cmpr_any )

  ! Deallocate resolved-scale source term and fall-back
  ! mass-source arrays, and 2D fields super-arrays
  deallocate( dndraft_fields_2d )
  deallocate( dndraft_fallback_par_gen )
  deallocate( dndraft_fallback_res_source )
  deallocate( dndraft_res_source )
  deallocate( updraft_fields_2d )
  deallocate( updraft_fallback_par_gen )
  deallocate( updraft_fallback_res_source )
  deallocate( updraft_res_source )
  ! Note: we could have got away with deallocating the
  ! fall-back mass sources *_fallback_par_gen sooner.
  ! However, doing so may create memory fragmentation
  ! problems, since the sub-arrays within the res_source and
  ! fall_back_par_gen arrays are all allocated together
  ! level-by-level.  Safest to deallocate them all together.


else  ! ( n_updraft_layers > 0 .or. n_dndraft_layers > 0 )


  ! Any calculations that need to be done even when there
  ! are no convecting layers:

  ! Convective cloud calculations at non-convecting points
  allocate( cmpr_any (k_bot_conv:k_top_conv) )
  do k = k_bot_conv, k_top_conv
    cmpr_any(k) % n_points = 0
    call zero_diag_conv_cloud( ij_first, ij_last,                              &
                               k, cmpr_any, cloudfracs )
  end do
  deallocate( cmpr_any )


end if  ! ( n_updraft_layers > 0 .or. n_dndraft_layers > 0 )


! Deallocate the initiating parcel property arrays
do k = k_bot_conv, k_top_conv
  do i_layr = 1, n_updraft_layers
    do i_type = 1, n_updraft_types
      if ( updraft_par_gen(i_type,i_layr,k) % cmpr % n_points > 0 ) then
        call parcel_dealloc( updraft_par_gen(i_type,i_layr,k) )
      end if
    end do
  end do
end do
deallocate( dndraft_par_gen )
deallocate( updraft_par_gen )


return
end subroutine comorph_main

end module comorph_main_mod
