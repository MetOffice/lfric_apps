! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module conv_sweep_ctl_mod

implicit none

contains

! Top-level subroutine for either an upwards or downwards
! sweep of the convection scheme's diagnostic parcel
! ascent/descent model.
! This generic routine handles updrafts or downdrafts,
! or overshoot fall-back flows from either.

! This routine includes:
! 1) The main loop over model-levels
!    (going up or down depending on the input switch l_down)
! 2) Call to calc_sum_massflux, to calculate the sum of existing
!    mass-fluxes at the current level over all convection types
!    (used for numerical safety checks within parcel calculation)
! 3) The main loop over convection types
!    (and the loop over distinct convecting layers)
! 4) At each model-level, for each convection type /
!    convecting layer:
!    a) Initialisation of new convecting parcels from the
!       initiating mass source properties calculated earlier.
!    b) Compression of environment fields onto points where
!       the current type / layer is convecting on the
!       current level.
!    c) Call to conv_level_step, which does one level-step
!       of the actual parcel ascent / descent model,
!       including calculation of resolved-scale source terms.
!    d) Consolidating and merging of the parcel compression
!       lists to account for points which terminate or newly
!       initiate.
!
subroutine conv_sweep_ctl( n_fields_tot,                                       &
                           n_conv_types, n_conv_layers,                        &
                           max_points, ij_first, ij_last,                      &
                           l_tracer, l_down, l_fallback,                       &
                           l_output_fallback,                                  &
                           grid, layer_mass, turb,                             &
                           fields, virt_temp,                                  &
                           par_gen, res_source, fields_2d,                     &
                           draft_diags, draft_diags_super,                     &
                           fallback_par_gen )

use comorph_constants_mod, only: nx_full, ny_full, k_bot_conv, k_top_conv,     &
                                 real_cvprec, real_hmprec, zero,               &
                                 max_ent_frac_up, max_ent_frac_up_fb,          &
                                 max_ent_frac_dn, max_ent_frac_dn_fb,          &
                                 i_check_bad_values_cmpr, i_check_bad_none,    &
                                 name_length, l_homog_conv_bl

use cmpr_type_mod, only: cmpr_type, cmpr_copy, cmpr_merge,                     &
                         cmpr_alloc, cmpr_dealloc
use grid_type_mod, only: grid_type, grid_compress, n_grid, i_pressure
use fields_type_mod, only: fields_type, fields_k_conserved_vars
use turb_type_mod, only: turb_type
use env_half_mod, only: n_env_half
use parcel_type_mod, only: parcel_type, i_massflux_d,                          &
                           parcel_alloc, parcel_compress, parcel_expand,       &
                           parcel_copy, parcel_combine,                        &
                           parcel_check_bad_values
use res_source_type_mod, only: res_source_type,                                &
                               res_source_alloc,                               &
                               res_source_dealloc,                             &
                               res_source_check_bad_values,                    &
                               res_source_init_zero,                           &
                               res_source_combine,                             &
                               res_source_expand
use fields_2d_mod, only: n_fields_2d
use entdet_res_source_mod, only: entdet_res_source
use draft_diags_type_mod, only: draft_diags_type,                              &
                                draft_diags_super_type
use diags_super_type_mod, only: diags_super_type,                              &
                                diags_super_alloc, diags_super_init_zero,      &
                                diags_super_expand, diags_super_combine
use parcel_diags_type_mod, only: parcel_diags_copy

use conv_sweep_compress_mod, only: conv_sweep_compress
use calc_sum_massflux_mod, only: calc_sum_massflux
use calc_delta_tv_mod, only: calc_delta_tv
use conv_level_step_mod, only: conv_level_step
use save_parcel_bl_top_mod, only: save_parcel_bl_top
use homog_conv_bl_ctl_mod, only: homog_conv_bl_ctl

implicit none

! Number of fields, accounting for whether tracer transport
! is required
integer, intent(in) :: n_fields_tot

! Number of convection types
integer, intent(in) :: n_conv_types

! Largest number of convective mass-source layers found in any
! column on the current segment
integer, intent(in) :: n_conv_layers

! Max size the compression lists can possibly need
! (i.e. the number of columns in the current segment)
integer, intent(in) :: max_points

! i,j indices of the first and last point on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Flag for whether to include passive tracers in the calculations
logical, intent(in) :: l_tracer

! Flag for whether this routine is being called for a downwards
! sweep (primary downdraft or updraft fall-back) as opposed to an
! upwards sweep (primary updraft or downdraft fall-back).
logical, intent(in) :: l_down

! Flag for whether this is the call for the fall-back flows;
! passed into conv_level_step, where it maybe used to
! specify a different detrainment calculation for
! fall-backs vs primary updrafts
logical, intent(in) :: l_fallback

! Flag for whether to ouput fall-back mass sources based
! on detrained mass which is sufficiently negatively buoyant
! to subsequently fall down below its detrainment level.
logical, intent(in) :: l_output_fallback

! Structure containing pointers to model grid fields
! (full 3-D arrays, possibly with haloes);
! contains model-level heights, pressures and dry-density
type(grid_type), intent(in) :: grid

! Full 3-D array of dry-mass per unit surface area
! contained in each grid-cell
real(kind=real_hmprec), intent(in) :: layer_mass                               &
       ( nx_full, ny_full, k_bot_conv:k_top_conv )
! Needed for calculating the max allowed entrainment
! by each convection type on each level;
! for numerical stability, the total entrainment
! summed over all types over the timestep must not exceed
! the mass of the layer

! Structure containing pointers to turbulence fields passed
! in from the boundary-layer scheme
type(turb_type), intent(in) :: turb

! Structure containing pointers to the _np1 fields;
! these are the primary fields already updated with any other
! increments computed before convection.
type(fields_type), intent(in) :: fields

! Full 3-D array of environment virtual temperature on full-levels
real(kind=real_hmprec), intent(in) :: virt_temp                                &
       ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Array of input structures containing the initiating parcel
! properties on each level for the current draft.
! Note: if this is the call for fall-back flows, then par_gen
! contains the fall-back mass-sources, and the output
! argument fallback_par_gen below is not used.
type(parcel_type), intent(in out) :: par_gen                                   &
       ( n_conv_types, n_conv_layers, k_bot_conv:k_top_conv )
! intent inout because we convert the initiating parcel properties
! to conserved variable form to combine with any existing parcel
! on each level.

! The following outputs need intent(inout) to preserve the
! assigned status of the contained allocatable arrays on input:

! Array of structures to store entrainment / detrainment and
! source terms used for calculation of large-scale increments
type(res_source_type), allocatable, intent(in out) :: res_source(:,:,:)

! Super-array storing 2D work arrays
real(kind=real_cvprec), allocatable, intent(in out) :: fields_2d(:,:,:,:)

! Structure storing flags and meta-data for diagnostics
type(draft_diags_type), intent(in) :: draft_diags

! Structure containing arrays storing diagnostics in compressed
! super-array form
type(draft_diags_super_type), intent(in out) :: draft_diags_super

! Array of structures containing the
! initiation mass source properties from fall-back flows
! originating in overshooting air detrained by the primary
! updraft or downdraft.
! If this call is for the primary updraft / downdraft
! (l_fallback==.FALSE.), but fall-backs are in use
! (l_output_fallback==.TRUE.), then this is an output.
! Otherwise it is not used.
! Beware: if you set both l_fallback and l_output_fallback to
! true (i.e. try to initiate fall-back flows from a fall-back
! flow), a horrid hell of memory-related bugs may be unleashed...
type(parcel_type), allocatable, intent(in out):: fallback_par_gen(:,:,:)
! This is an output, but needs intent inout to preserve its
! allocation status on input.


! LOCAL VARIABLES

! Parcel properties that are carried during the ascent/descent
type(parcel_type), allocatable :: par_conv(:,:)

! Sum of parcel mass-fluxes over types/layers at previous model-level interface
real(kind=real_cvprec) :: sum_massflux(ij_first:ij_last)
! Sum of initiating mass-sources over types/layers at current full-level
real(kind=real_cvprec) :: sum_massinit(ij_first:ij_last)

! Copy of the above compressed onto points where the current
! convection type/layer is active.
real(kind=real_cvprec) :: sum_massflux_cmpr(max_points)

! Difference in virtual temperature between level k and next full-level,
! divided by layer-mass (used to pre-estimate subsidence increment to Tv)
real(kind=real_cvprec) :: delta_tv(ij_first:ij_last)

! Value compressed onto current type / layer points
real(kind=real_cvprec) :: delta_tv_cmpr(max_points)

! Flags for whether any convection types / layers are
! active at the current level
logical :: l_any_conv
logical :: l_any_source

! Loop bounds and stride to use in the main level-loop
integer :: k_first, k_last, dk

! Model-level index for the next model-level interface
! (k+1/2 for updrafts, k-1/2 for downdrafts)
integer :: k_next

! Note: we have a convention throughout that "_next" refers to the
! next model-level interface, NOT the next full-level (the latter
! is denoted by k+dk).  This is because the convection scheme
! integrates over level-steps from k-1/2 to k+1/2 (or the other
! way round for downdrafts).  Thus "prev" refers to variables at
! the previous model-level interface where each level-step begins,
! while "next" refers to the updated variables at the end of
! each level-step.
! Therefore, k_next is set to subscript the fields defined on the
! next model-level interface / half-level (rho-level in UM-speak).

! Flag for when reached the last model-level
! (any remaining parcels must fully detrain when this is true)
logical :: l_last_level

! Flag for first half of the level-step, from prev to level k
logical :: l_to_full_level


! Local arrays used to store various fields compressed onto
! convecting points at the current level,
! for the current convection type only:

! (note that these arrays are dimensioned by max_points,
!  which is the largest size the compression lists will ever
!  need to be over all levels and types; this avoids ever
!  having to deallocate / reallocate them during the loop
!  over levels and convection types).

! Model grid arrays
real(kind=real_cvprec) :: grid_k_super(max_points,n_grid)
real(kind=real_cvprec) :: grid_prev_super(max_points,n_grid)
real(kind=real_cvprec) :: grid_next_super(max_points,n_grid)

! Environment primary fields at the current level
real(kind=real_cvprec) :: env_k_fields(max_points,n_fields_tot)

! Layer mass on current level-step
real(kind=real_cvprec) :: layer_mass_step(max_points)
! Fraction of layer on current half-level step
real(kind=real_cvprec) :: frac_level_step(max_points)

! Arrays containing certain environment fields which are needed
! at both half-levels and full-levels
real(kind=real_cvprec) :: env_k_super(max_points,n_env_half)
real(kind=real_cvprec) :: env_prev_super(max_points,n_env_half)
real(kind=real_cvprec) :: env_next_super(max_points,n_env_half)

! Flag for whether level k is within the boundary-layer
! (as defined by the BL-scheme's BL-top height passed in)
logical :: l_within_bl(max_points)

! Maximum allowed entrained fraction of layer mass for current draft
real(kind=real_cvprec) :: max_ent_frac

! Space for plume-model diagostics from the 2nd conv_level_step call
type(diags_super_type) :: plume_model_diags_step

! Resolved-scale source terms due to the initiating
! mass source only
type(res_source_type) :: res_source_gen

! Parcel properties at the diagnosed boundary-layer top
! (used if homogenising the increments below the BL-top)
type(parcel_type), allocatable :: par_bl_top(:,:,:)

! Compression indices used for combining the existing and initiating parcels
type(cmpr_type) :: cmpr_combined ( n_conv_types, n_conv_layers )
integer :: index_ic_conv ( max_points, n_conv_types, n_conv_layers )
integer :: index_ic_gen  ( max_points, n_conv_types, n_conv_layers )

! Indices used to reference a collapsed horizontal coordinate from the parcel
integer :: index_ij(max_points)

! Points where any type / layer active at current height
type(cmpr_type) :: cmpr_any
! Points for calculating delta_tv
type(cmpr_type) :: cmpr_tmp

! Flag passed into fields_k_conserved_vars
logical :: l_reverse
! Flag input to entdet_res_source to indicate entrainment vs detrainment
logical :: l_ent

! Character strings used in error messages:
! String identifying what sort of draft this call is for
character(len=name_length) :: draft_string
! Description of where we are in the code
character(len=name_length) :: where_string

! Dimensions of the diagnostics super-array passed into conv_level_step
integer :: dgdims(2)

! Dimension of the initiating parcel super-arrays
integer :: n_points_gen

! Loop counters
integer :: ic, ij, k, i_type, i_layr, i_field


! Set character string indicating updraft or downdraft etc
if ( l_down ) then
  if ( l_fallback ) then
    draft_string = "updraft fall-back"
  else
    draft_string = "primary downdraft"
  end if
else
  if ( l_fallback ) then
    draft_string = "downdraft fall-back"
  else
    draft_string = "primary updraft"
  end if
end if


!----------------------------------------------------------------
! 1) Allocate and initialise parcel arrays etc...
!----------------------------------------------------------------

! Setup arrays for the output resolved-scale source terms
allocate( res_source( n_conv_types, n_conv_layers,                             &
                      k_bot_conv:k_top_conv ) )
! Initialise number of points in the compression lists to zero
do k = k_bot_conv, k_top_conv
  do i_layr = 1, n_conv_layers
    do i_type = 1, n_conv_types
      res_source(i_type,i_layr,k) % cmpr % n_points = 0
    end do
  end do
end do
! Note: the compression arrays contained in the res_source
! structures are not allocated yet,
! since we don't yet know how many convecting points there will
! be on each model-level.  We allocate them level-by-level
! during the parcel ascent/descent, as and when we have data
! to populate them with.

if ( .not. l_fallback ) then
  ! The following arrays are already initialised in the primary sweep call
  ! if this is a fall-back sweep call...
  if ( n_fields_2d > 0 ) then
    ! If any 2D work arrays are needed...
    ! Allocate super-array for 2D work fields
    allocate( fields_2d( ij_first:ij_last, n_fields_2d,                        &
                         n_conv_types, n_conv_layers ) )
    ! Initialise to zero
    do i_layr = 1, n_conv_layers
      do i_type = 1, n_conv_types
        do i_field = 1, n_fields_2d
          do ij = ij_first, ij_last
            fields_2d(ij,i_field,i_type,i_layr) = zero
          end do
        end do
      end do
    end do
  else
    ! Minimal allocation when not used
    allocate( fields_2d( 1, 1, n_conv_types, n_conv_layers ) )
  end if
end if  ! ( .NOT. l_fallback )

! If fall-back flows are on, but this is not the fall-back call,
! allocate arrays for outputting the fall-back mass-sources to
! pass to the fall-back call
if ( .not. l_fallback ) then
  ! Needs allocating even if not used, to avoid error when subscripting
  ! in the argument list to conv_level_step
  allocate( fallback_par_gen( n_conv_types, n_conv_layers,                     &
                              k_bot_conv:k_top_conv ) )
  if ( l_output_fallback ) then
    ! Initialise number of points to zero
    do k = k_bot_conv, k_top_conv
      do i_layr = 1, n_conv_layers
        do i_type = 1, n_conv_types
          fallback_par_gen(i_type,i_layr,k) % cmpr % n_points = 0
        end do
      end do
    end do
  end if
end if

! Setup arrays for output diagnostics
! Note: these need to be allocated even if not used, to avoid
! error due to subscripting the array and passing it
! into other routines when not allocated.
allocate( draft_diags_super % par                                              &
  ( n_conv_types, n_conv_layers, k_bot_conv:k_top_conv ) )
allocate( draft_diags_super % gen                                              &
  ( n_conv_types, n_conv_layers, k_bot_conv:k_top_conv ) )
allocate( draft_diags_super % plume_model                                      &
  ( n_conv_types, n_conv_layers, k_bot_conv:k_top_conv ) )
! Initialise number of points in the compression lists to zero
do k = k_bot_conv, k_top_conv
  do i_layr = 1, n_conv_layers
    do i_type = 1, n_conv_types
      draft_diags_super % par(i_type,i_layr,k) % cmpr % n_points = 0
      draft_diags_super % gen(i_type,i_layr,k) % cmpr % n_points = 0
      draft_diags_super % plume_model(i_type,i_layr,k) % cmpr % n_points = 0
    end do
  end do
end do

! Setup the parcel property structures
allocate( par_conv ( n_conv_types, n_conv_layers ) )
do i_layr = 1, n_conv_layers
  do i_type = 1, n_conv_types

    ! Allocate parcel arrays to max required size
    call parcel_alloc( l_tracer, max_points,                                   &
                       par_conv(i_type,i_layr) )

    ! Initialise number of points to 0 ready for start of level loop
    par_conv(i_type,i_layr) % cmpr % n_points = 0

  end do
end do

! Allocate compression arrays for the resolved-scale source terms
! from only the initiation mass-sources
call res_source_alloc( l_tracer, max_points, res_source_gen )
res_source_gen % cmpr % n_points = 0

! Allocate compression indices used for merging existing and initiating parcels
do i_layr = 1, n_conv_layers
  do i_type = 1, n_conv_types
    call cmpr_alloc( cmpr_combined(i_type,i_layr), max_points )
  end do
end do

! Allocate additional compression indices
call cmpr_alloc( cmpr_any, max_points )
call cmpr_alloc( cmpr_tmp, max_points )

! Initialise structures to store parcel properties at the boundary-layer
! top if needed
if ( l_homog_conv_bl ) then
  allocate( par_bl_top ( n_conv_types, n_conv_layers, k_bot_conv:k_top_conv ) )
  do k = k_bot_conv, k_top_conv
    do i_layr = 1, n_conv_layers
      do i_type = 1, n_conv_types
        par_bl_top(i_type,i_layr,k) % cmpr % n_points = 0
      end do
    end do
  end do
end if  ! ( l_homog_conv_bl )

if ( draft_diags % plume_model % n_diags > 0 ) then
  ! Allocate space for plume-model diagnostics from 2nd conv_level_step call
  call diags_super_alloc( max_points, draft_diags%plume_model%n_diags_super,   &
                          plume_model_diags_step )
else
  allocate( plume_model_diags_step % super(1,1) )
end if

! Select the correct max entrainment fraction for this draft
if ( l_down ) then
  if ( l_fallback ) then
    max_ent_frac = max_ent_frac_dn_fb
  else
    max_ent_frac = max_ent_frac_dn
  end if
else
  if ( l_fallback ) then
    max_ent_frac = max_ent_frac_up_fb
  else
    max_ent_frac = max_ent_frac_up
  end if
end if

! Set vertical level loop indices based on whether going
! up or down
if ( l_down ) then
  k_first = k_top_conv
  k_last = k_bot_conv
  dk = -1
else
  k_first = k_bot_conv
  k_last = k_top_conv
  dk = 1
end if


!----------------------------------------------------------------
! 2) Loop over model-levels
!----------------------------------------------------------------
do k = k_first, k_last, dk

  ! Set flags for whether there is any existing convection
  ! rising from level k-1/2, or any initiation mass-sources
  ! at level k, for any convection type/layer
  l_any_conv = .false.
  l_any_source = .false.
  do i_layr = 1, n_conv_layers
    do i_type = 1, n_conv_types
      if ( par_conv(i_type,i_layr) % cmpr                                      &
           % n_points > 0 )  l_any_conv = .true.
      if ( par_gen(i_type,i_layr,k) % cmpr                                     &
           % n_points > 0 )  l_any_source = .true.
    end do
  end do

  ! Only anything else to do if existing convection or new mass-source found
  if ( l_any_conv .or. l_any_source ) then

    ! Set flag for whether we've reached the final level
    l_last_level = ( k == k_last )

    ! Set k-index of the next model-level interface
    ! (k+1/2 for updrafts, k-1/2 for downdrafts)
    if ( l_down ) then
      k_next = k
    else
      k_next = k + 1
    end if
    ! k is vertical subscript for fields defined on the current
    !   theta-level.
    ! k_next is vertical subscript for fields defined on the next
    !   rho-level adjacent to k.
    ! k+dk is vertical subscript for fields defined on the next
    !   theta-level after the rho-level at k_next.
    !
    !   Model-level  Updrafts:  Downdrafts:
    !
    !  +~~  k+1  ~~+~~  k+dk ~~+~~       ~~+
    !  |           |           |           |
    !  |           |           |           |
    !  +-- k+1/2 --+-- k_next--+-- prev  --+
    !  |           |     ^     |     v     |
    !  |           |     ^     |     v     |
    !  +~~   k   ~~+~~   k   ~~+~~   k   ~~+
    !  |           |     ^     |     v     |
    !  |           |     ^     |     v     |
    !  +-- k-1/2 --+-- prev  --+-- k_next--+
    !  |           |           |           |
    !  |           |           |           |
    !  +~~  k-1  ~~+~~       ~~+~~  k+dk ~~+
    !
    ! Note that the model-level step moves the parcel
    ! from prev to next.  The fields at k+dk are only used
    ! for interpolating onto the next model-level interface
    ! at k_next.

    ! First sum the mass-fluxes and initiating masses over all layers / types.
    ! These are used in various checks...
    call calc_sum_massflux( n_conv_types, n_conv_layers, ij_first, ij_last,    &
                            par_conv, sum_massflux )
    call calc_sum_massflux( n_conv_types, n_conv_layers, ij_first, ij_last,    &
                            par_gen(:,:,k), sum_massinit )

    !------------------------------------------------------
    ! 3) Calculate numerical corrections
    !------------------------------------------------------

    ! Find points where any type / layer has any mass-flux or initiating
    ! mass-sources
    cmpr_any%n_points = 0
    do ij = ij_first, ij_last
      if ( sum_massflux(ij) > zero .or. sum_massinit(ij) > zero ) then
        cmpr_any%n_points = cmpr_any%n_points + 1
        cmpr_any%index_j(cmpr_any%n_points) = 1 + (ij-1) / nx_full
        cmpr_any%index_i(cmpr_any%n_points) =                                  &
                      ij - nx_full * (cmpr_any%index_j(cmpr_any%n_points)-1)
      end if
    end do

    ! Initialise delta_tv to zero at those points
    do ic = 1, cmpr_any%n_points
      ij = cmpr_any%index_i(ic) + nx_full * ( cmpr_any%index_j(ic) - 1 )
      delta_tv(ij) = zero
    end do

    ! TEMPORARILY COMMENTED-OUT TO PRESERVE KGO:
    ! When in the last model-level we can't really have any compensating
    ! subsidence to treat implicitly (since there is no next level to subside).
    ! For now, keep calculation of delta_tv from subsidence during the
    ! last level but remove this soon...
    !IF ( .NOT. l_last_level ) THEN
    if ( .true. ) then
       ! Calculation needs data from k+1; leave as zero at last level

      ! Find points where mass-flux is non-zero
      cmpr_tmp%n_points = 0
      do ic = 1, cmpr_any%n_points
        ij = cmpr_any%index_i(ic) + nx_full * ( cmpr_any%index_j(ic) - 1 )
        if ( sum_massflux(ij) > zero ) then
          cmpr_tmp%n_points = cmpr_tmp%n_points + 1
          cmpr_tmp%index_i(cmpr_tmp%n_points) = cmpr_any%index_i(ic)
          cmpr_tmp%index_j(cmpr_tmp%n_points) = cmpr_any%index_j(ic)
        end if
      end do

      ! Pre-estimate subsidence virtual temperature increment per unit
      ! mass-flux, used for the implicit detrainment
      if ( cmpr_tmp%n_points > 0 ) then
        call calc_delta_tv( .false., l_last_level,                             &
                            cmpr_tmp, k, k_next, dk, ij_first, ij_last,        &
                            virt_temp, layer_mass, grid, fields,               &
                            delta_tv )
      end if

    end if  ! ( .NOT. l_last_level )


    !------------------------------------------------------
    ! 4) Allocate arrays required for calculations on this level
    !------------------------------------------------------

    do i_layr = 1, n_conv_layers
      do i_type = 1, n_conv_types

        ! Find the combined compression list which includes all points
        ! from both the existing parcel and the initiating parcel
        call cmpr_merge( par_conv(i_type,i_layr) % cmpr,                       &
                         par_gen(i_type,i_layr,k) % cmpr,                      &
                         cmpr_combined(i_type,i_layr),                         &
                         index_ic_conv(:,i_type,i_layr),                       &
                         index_ic_gen(:,i_type,i_layr) )

        if ( cmpr_combined(i_type,i_layr) % n_points > 0 ) then
          ! If any points have either existing convection or initiation...

          ! Set size of resolved-scale source-term arrays and allocate
          call res_source_alloc( l_tracer,                                     &
                                 cmpr_combined(i_type,i_layr) % n_points,      &
                                 res_source(i_type,i_layr,k) )
          ! Initialise to zero
          res_source(i_type,i_layr,k) % cmpr % n_points                        &
            = cmpr_combined(i_type,i_layr) % n_points
          call res_source_init_zero( l_tracer, res_source(i_type,i_layr,k) )

          ! Allocate diagnostics super-array and initialise to zero
          if ( draft_diags % plume_model % n_diags > 0 ) then
            call diags_super_alloc( cmpr_combined(i_type,i_layr) % n_points,   &
                                    draft_diags % plume_model % n_diags_super, &
                      draft_diags_super % plume_model(i_type,i_layr,k) )
            ! Initialise to zero
            draft_diags_super % plume_model(i_type,i_layr,k)%cmpr%n_points     &
              = cmpr_combined(i_type,i_layr) % n_points
            call diags_super_init_zero( draft_diags%plume_model%n_diags_super, &
                      draft_diags_super % plume_model(i_type,i_layr,k) )
          else
            allocate( draft_diags_super % plume_model(i_type,i_layr,k)         &
                      % super(1,1) )
          end if

        end if  ! ( cmpr_combined(i_type,i_layr) % n_points > 0 )

      end do  ! i_type = 1, n_conv_types
    end do  ! i_layr = 1, n_conv_layers


    !------------------------------------------------------
    ! 5) Calculations for existing convection from previous level
    !------------------------------------------------------

    ! Set flag for 1st half of the level-step, from prev to level k
    l_to_full_level = .true.



    do i_layr = 1, n_conv_layers
      do i_type = 1, n_conv_types

        if ( par_conv(i_type,i_layr) % cmpr % n_points > 0 ) then
          ! If any existing convection of this type / layer...

          ! Set source-term and diagnostics compression indices
          ! the same as the parcel
          call cmpr_copy( par_conv(i_type,i_layr) % cmpr,                      &
                          res_source(i_type,i_layr,k) % cmpr )
          if ( draft_diags % plume_model % n_diags > 0 ) then
            call cmpr_copy( par_conv(i_type,i_layr) % cmpr,                    &
                      draft_diags_super % plume_model(i_type,i_layr,k) % cmpr )
          end if
          dgdims = ubound( draft_diags_super % plume_model(i_type,i_layr,k)    &
                           % super )

          ! Compress environment fields onto convecting points
          call conv_sweep_compress(                                            &
                 k, k_next, dk, max_points, ij_first, ij_last, n_fields_tot,   &
                 l_to_full_level, l_last_level,                                &
                 par_conv(i_type,i_layr) % cmpr,                               &
                 grid, turb, fields,                                           &
                 virt_temp, layer_mass, sum_massflux, delta_tv,                &
                 l_within_bl, grid_k_super, grid_prev_super,                   &
                 env_k_fields, env_k_super, env_prev_super,                    &
                 layer_mass_step, frac_level_step,                             &
                 delta_tv_cmpr, sum_massflux_cmpr, index_ij )

          ! Check input parcel properties for bad values (NaN, Inf, etc)
          if ( i_check_bad_values_cmpr > i_check_bad_none ) then
            where_string = "On input to 1st conv_level_step call for "      // &
                           trim(adjustl(draft_string)) // "; par_conv"
            call parcel_check_bad_values( par_conv(i_type,i_layr),             &
                                          n_fields_tot, k, where_string )
            where_string = "On input to 1st conv_level_step call for "      // &
                           trim(adjustl(draft_string)) // "; res_source"
            call res_source_check_bad_values( res_source(i_type,i_layr,k),     &
                                              n_fields_tot, k, where_string )
          end if

          ! Compute the parcel ascent / descent from the previous
          ! model-level interface to level k
          call conv_level_step(                                                &
                 par_conv(i_type,i_layr)%cmpr%n_points, max_points,            &
                 cmpr_combined(i_type,i_layr)%n_points,                        &
                 dgdims(1), dgdims(2), n_fields_tot,                           &
                 l_down, l_tracer, l_last_level, l_to_full_level,              &
                 l_within_bl,                                                  &
                 sum_massflux_cmpr, max_ent_frac, index_ij, ij_first, ij_last, &
                 par_conv(i_type,i_layr)%cmpr, k, draft_string,                &
                 grid_prev_super, grid_k_super,                                &
                 env_prev_super, env_k_super,                                  &
                 env_k_fields, layer_mass_step, frac_level_step,               &
                 delta_tv_cmpr,                                                &
                 par_conv(i_type,i_layr)%par_super,                            &
                 par_conv(i_type,i_layr)%mean_super,                           &
                 par_conv(i_type,i_layr)%core_super,                           &
                 res_source(i_type,i_layr,k)%res_super,                        &
                 res_source(i_type,i_layr,k)%fields_super,                     &
                 res_source(i_type,i_layr,k)%convcloud_super,                  &
                 fields_2d(:,:,i_type,i_layr),                                 &
                 draft_diags % plume_model,                                    &
                 draft_diags_super % plume_model(i_type,i_layr,k) % super )

          ! Check output parcel properties for bad values (NaN, Inf, etc)
          if ( i_check_bad_values_cmpr > i_check_bad_none ) then
            where_string = "On output from 1st conv_level_step call for "   // &
                           trim(adjustl(draft_string)) // "; par_conv"
            call parcel_check_bad_values( par_conv(i_type,i_layr),             &
                                          n_fields_tot, k, where_string )
            where_string = "On output from 1st conv_level_step call for "   // &
                           trim(adjustl(draft_string)) // "; res_source"
            call res_source_check_bad_values( res_source(i_type,i_layr,k),     &
                                              n_fields_tot, k, where_string )
          end if

        end if  ! ( par_conv(i_type,i_layr) % cmpr % n_points > 0 )

      end do  ! i_type = 1, n_conv_types
    end do  ! i_layr = 1, n_conv_layers


    !------------------------------------------------------
    ! 6) Add on new initiating mass-sources from level k
    !------------------------------------------------------

    do i_layr = 1, n_conv_layers
      do i_type = 1, n_conv_types

        ! Convert the updated parcel properties to conserved
        ! variable form before combining them with initiating
        ! mass-sources
        if ( par_conv(i_type,i_layr) % cmpr % n_points > 0 ) then
          l_reverse = .false.
          call fields_k_conserved_vars(                                        &
                 par_conv(i_type,i_layr) % cmpr % n_points,                    &
                 max_points, n_fields_tot, l_reverse,                          &
                 par_conv(i_type,i_layr) % mean_super )
        end if
        ! Note: this is done even if there are no initiating
        ! mass-sources, to ensure bit-reproducibility if we
        ! change the processor decomposition such that
        ! whether or not there are any mass-sources on the
        ! local proc domain changes.

        if ( par_gen(i_type,i_layr,k) % cmpr % n_points > 0 ) then
          ! If any initiating mass-source of this type / layer...

          ! Find super-array size for initiating parcel arrays
          n_points_gen = ubound( par_gen(i_type,i_layr,k) % mean_super, 1 )

          ! Convert initiating parcel properties to conserved-variable form
          l_reverse = .false.
          call fields_k_conserved_vars( par_gen(i_type,i_layr,k)%cmpr%n_points,&
                                        n_points_gen, n_fields_tot, l_reverse, &
                                        par_gen(i_type,i_layr,k) % mean_super )

          if ( par_conv(i_type,i_layr) % cmpr % n_points > 0 ) then
            ! If there is existing convection on the same level as
            ! initiating mass-sources...

            ! Initialise resolved-scale source terms to zero
            res_source_gen % cmpr % n_points                                   &
              = par_gen(i_type,i_layr,k) % cmpr % n_points
            call res_source_init_zero( l_tracer, res_source_gen )

            ! Calculate resolved-scale sources due to the initiating parcel
            ! (treated as entrainment)
            l_ent = .true.
            call entdet_res_source( par_gen(i_type,i_layr,k)%cmpr%n_points,    &
                                    n_points_gen,                              &
                                    max_points, n_fields_tot, l_ent,           &
                                    par_gen(i_type,i_layr,k)                   &
                                       % par_super(:,i_massflux_d),            &
                                    par_gen(i_type,i_layr,k) % mean_super,     &
                                    res_source_gen % res_super,                &
                                    res_source_gen % fields_super )

            if ( cmpr_combined(i_type,i_layr) % n_points >                     &
                 par_conv(i_type,i_layr) % cmpr % n_points ) then
              ! If initiating mass-sources have added new points to the
              ! resulting combined compression list, expand all the compression
              ! arrays onto the new list
              call parcel_expand( l_tracer,                                    &
                                  cmpr_combined(i_type,i_layr) % n_points,     &
                                  index_ic_conv(:,i_type,i_layr),              &
                                  par_conv(i_type,i_layr) )
              call cmpr_copy( cmpr_combined(i_type,i_layr),                    &
                              par_conv(i_type,i_layr) % cmpr )
              call res_source_expand( l_tracer,                                &
                                      cmpr_combined(i_type,i_layr) % n_points, &
                                      index_ic_conv(:,i_type,i_layr),          &
                                      res_source(i_type,i_layr,k) )
              if ( draft_diags % plume_model % n_diags > 0 ) then
                call diags_super_expand(                                       &
                      draft_diags % plume_model % n_diags_super,               &
                      cmpr_combined(i_type,i_layr) % n_points,                 &
                      index_ic_conv(:,i_type,i_layr),                          &
                      draft_diags_super % plume_model(i_type,i_layr,k) )
              end if
            end if

            ! Combine parcel properties and resolved-scale source-terms due
            ! to initiation with those from the existing convection
            call parcel_combine( l_tracer, l_down,                             &
                                 index_ic_gen(:,i_type,i_layr),                &
                                 par_gen(i_type,i_layr,k),                     &
                                 par_conv(i_type,i_layr) )
            call res_source_combine( l_tracer, index_ic_gen(:,i_type,i_layr),  &
                                     res_source_gen,                           &
                                     res_source(i_type,i_layr,k) )

          else  ! ( par_conv(i_type,i_layr) % cmpr % n_points == 0 )
            ! Initiating mass-source but no existing convection to combine with

            ! Calculate resolved-scale sources due to the initiating parcel
            ! (treated as entrainment)
            l_ent = .true.
            call entdet_res_source( par_gen(i_type,i_layr,k)%cmpr%n_points,    &
                                    n_points_gen,                              &
                                    par_gen(i_type,i_layr,k)%cmpr%n_points,    &
                                    n_fields_tot, l_ent,                       &
                                    par_gen(i_type,i_layr,k)                   &
                                         % par_super(:,i_massflux_d),          &
                                    par_gen(i_type,i_layr,k) % mean_super,     &
                                    res_source(i_type,i_layr,k) % res_super,   &
                                    res_source(i_type,i_layr,k) % fields_super)

            ! Copy the initiating parcel properties into par_conv
            call parcel_copy( l_tracer, par_gen(i_type,i_layr,k),              &
                              par_conv(i_type,i_layr) )

          end if  ! ( par_conv(i_type,i_layr) % cmpr % n_points == 0 )

          ! Update compression indices consistent with the new list
          call cmpr_copy( cmpr_combined(i_type,i_layr),                        &
                          res_source(i_type,i_layr,k) % cmpr )
          if ( draft_diags % plume_model % n_diags > 0 ) then
            call cmpr_copy( cmpr_combined(i_type,i_layr),                      &
                      draft_diags_super % plume_model(i_type,i_layr,k) % cmpr )
          end if

          ! Save any requested diagnostics of initiating parcel properties
          if ( draft_diags % gen % n_diags > 0 ) then
            ! Note: in the calls below, the diagnostics are valid
            ! at the full-level k
            ! Allocate super-array for diagnostics
            call diags_super_alloc(                                            &
                   par_gen(i_type,i_layr,k) % cmpr % n_points,                 &
                   draft_diags % gen % n_diags_super,                          &
                   draft_diags_super % gen(i_type,i_layr,k) )
            ! Copy parcel properties for diags
            call cmpr_copy( par_gen(i_type,i_layr,k) % cmpr,                   &
                            draft_diags_super % gen(i_type,i_layr,k) % cmpr )
            call grid_compress( grid, par_gen(i_type,i_layr,k) % cmpr,         &
                                k=k, grid_k_super=grid_k_super )
            call parcel_diags_copy( n_fields_tot,                              &
                   draft_diags % gen,                                          &
                   par_gen(i_type,i_layr,k), grid_k_super(:,i_pressure),       &
                   draft_diags_super % gen(i_type,i_layr,k) % super )
          end if

        end if  ! ( par_gen(i_type,i_layr,k) % cmpr % n_points > 0 )

        ! Convert parcel properties back from conserved variable form
        if ( par_conv(i_type,i_layr) % cmpr % n_points > 0 ) then
          l_reverse = .true.
          call fields_k_conserved_vars(                                        &
                 par_conv(i_type,i_layr) % cmpr % n_points,                    &
                 max_points, n_fields_tot, l_reverse,                          &
                 par_conv(i_type,i_layr) % mean_super )
        end if

      end do  ! i_type = 1, n_conv_types
    end do  ! i_layr = 1, n_conv_layers


    !------------------------------------------------------
    ! 7) Calculations for convection from level k to next
    !------------------------------------------------------

    ! Set flag for 2nd half of the level-step, from k to next
    l_to_full_level = .false.

    ! Calculate updated sum of mass-fluxes
    call calc_sum_massflux( n_conv_types, n_conv_layers, ij_first, ij_last,    &
                            par_conv, sum_massflux )

    ! TEMPORARILY COMMENTED-OUT TO PRESERVE KGO:
    ! When in the last model-level we can't really have any compensating
    ! subsidence to treat implicitly (since there is no next level to subside).
    ! For now, keep calculation of delta_tv from subsidence during the
    ! last level but remove this soon...
    !IF ( .NOT. l_last_level ) THEN
    if ( .true. ) then
      ! Calculation needs data from k+1; leave as zero at last level

      ! Find points where mass-flux is non-zero
      cmpr_tmp%n_points = 0
      do ic = 1, cmpr_any%n_points
        ij = cmpr_any%index_i(ic) + nx_full * ( cmpr_any%index_j(ic) - 1 )
        if ( sum_massflux(ij) > zero ) then
          cmpr_tmp%n_points = cmpr_tmp%n_points + 1
          cmpr_tmp%index_i(cmpr_tmp%n_points) = cmpr_any%index_i(ic)
          cmpr_tmp%index_j(cmpr_tmp%n_points) = cmpr_any%index_j(ic)
        end if
      end do

      ! Pre-estimate subsidence virtual temperature increment per unit
      ! mass-flux, used for the implicit detrainment
      if ( cmpr_tmp%n_points > 0 ) then
        call calc_delta_tv( .true., l_last_level,                              &
                            cmpr_tmp, k, k_next, dk, ij_first, ij_last,        &
                            virt_temp, layer_mass, grid, fields,               &
                            delta_tv )
      end if

    end if  ! ( .NOT. l_last_level )

    ! Loop over convection layers and types
    do i_layr = 1, n_conv_layers
      do i_type = 1, n_conv_types
        if ( par_conv(i_type,i_layr) % cmpr % n_points > 0 ) then
          ! If any convection of this type / layer...

          dgdims = ubound( plume_model_diags_step % super )

          ! Compress environment fields onto convecting points
          call conv_sweep_compress(                                            &
                 k, k_next, dk, max_points, ij_first, ij_last, n_fields_tot,   &
                 l_to_full_level, l_last_level,                                &
                 par_conv(i_type,i_layr) % cmpr,                               &
                 grid, turb, fields,                                           &
                 virt_temp, layer_mass, sum_massflux, delta_tv,                &
                 l_within_bl, grid_k_super, grid_next_super,                   &
                 env_k_fields, env_k_super, env_next_super,                    &
                 layer_mass_step, frac_level_step,                             &
                 delta_tv_cmpr, sum_massflux_cmpr, index_ij )

          ! Check input parcel properties for bad values (NaN, Inf, etc)
          if ( i_check_bad_values_cmpr > i_check_bad_none ) then
            where_string = "On input to 2nd conv_level_step call for "      // &
                           trim(adjustl(draft_string)) // "; par_conv"
            call parcel_check_bad_values( par_conv(i_type,i_layr),             &
                                          n_fields_tot, k, where_string )
            where_string = "On input to 2nd conv_level_step call for "      // &
                           trim(adjustl(draft_string)) // "; res_source"
            call res_source_check_bad_values( res_source(i_type,i_layr,k),     &
                                              n_fields_tot, k, where_string )
          end if

          ! Compute the parcel ascent / descent from level k to the
          ! next model-level interface
          call conv_level_step(                                                &
                 par_conv(i_type,i_layr)%cmpr%n_points, max_points,            &
                 res_source(i_type,i_layr,k)%cmpr%n_points,                    &
                 dgdims(1), dgdims(2), n_fields_tot,                           &
                 l_down, l_tracer, l_last_level, l_to_full_level,              &
                 l_within_bl,                                                  &
                 sum_massflux_cmpr, max_ent_frac, index_ij, ij_first, ij_last, &
                 par_conv(i_type,i_layr)%cmpr, k, draft_string,                &
                 grid_k_super, grid_next_super,                                &
                 env_k_super, env_next_super,                                  &
                 env_k_fields, layer_mass_step, frac_level_step,               &
                 delta_tv_cmpr,                                                &
                 par_conv(i_type,i_layr)%par_super,                            &
                 par_conv(i_type,i_layr)%mean_super,                           &
                 par_conv(i_type,i_layr)%core_super,                           &
                 res_source(i_type,i_layr,k)%res_super,                        &
                 res_source(i_type,i_layr,k)%fields_super,                     &
                 res_source(i_type,i_layr,k)%convcloud_super,                  &
                 fields_2d(:,:,i_type,i_layr),                                 &
                 draft_diags % plume_model,                                    &
                 plume_model_diags_step % super )

          ! Check output parcel properties for bad values (NaN, Inf, etc)
          if ( i_check_bad_values_cmpr > i_check_bad_none ) then
            where_string = "On output from 2nd conv_level_step call for "   // &
                           trim(adjustl(draft_string)) // "; par_conv"
            call parcel_check_bad_values( par_conv(i_type,i_layr),             &
                                          n_fields_tot, k, where_string )
            where_string = "On output from 2nd conv_level_step call for "   // &
                           trim(adjustl(draft_string)) // "; res_source"
            call res_source_check_bad_values( res_source(i_type,i_layr,k),     &
                                              n_fields_tot, k, where_string )
          end if

          do ic = 1, par_conv(i_type,i_layr) % cmpr % n_points
            index_ic_conv(ic,i_type,i_layr) = ic
          end do

          if ( draft_diags % plume_model % n_diags > 0 ) then
            ! Combine diagnostics from 1st half-level-step with those
            ! from the 2nd
            plume_model_diags_step % cmpr % n_points                           &
              = par_conv(i_type,i_layr) % cmpr % n_points
            call diags_super_combine(                                          &
                   draft_diags % plume_model % n_diags_super,                  &
                   draft_diags % plume_model % l_weight,                       &
                   index_ic_conv(:,i_type,i_layr),                             &
                   plume_model_diags_step,                                     &
                   draft_diags_super % plume_model(i_type,i_layr,k) )
          end if

        end if  ! ( par_conv(i_type,i_layr) % cmpr % n_points > 0 )
      end do  ! i_type = 1, n_conv_types
    end do  ! i_layr = 1, n_conv_layers


    !------------------------------------------------------
    ! 8) Save final parcel properties at next model-level interface
    !------------------------------------------------------
    do i_layr = 1, n_conv_layers
      do i_type = 1, n_conv_types
        if ( par_conv(i_type,i_layr) % cmpr % n_points > 0 ) then

          ! Convert parcel properties to conserved form
          l_reverse = .false.
          call fields_k_conserved_vars(                                        &
                 par_conv(i_type,i_layr) % cmpr % n_points,                    &
                 max_points, n_fields_tot, l_reverse,                          &
                 par_conv(i_type,i_layr) % mean_super )

          if ( l_homog_conv_bl ) then
            ! Save parcel properties at the first layer interface
            ! beyond the boundary-layer top
            call save_parcel_bl_top( n_fields_tot, k, dk,                      &
                                     grid, turb,                               &
                                     par_conv(i_type,i_layr),                  &
                                     par_bl_top(i_type,i_layr,k) )
            ! Note: the BL-top parcel is saved at the current
            ! full-level k.  However, it is actually now defined at the
            ! next model-level interface, which is just above for updrafts,
            ! just below for downdrafts.  Care needs to be taken to account
            ! for this difference correctly for both updrafts and downdrafts!
          end if

          ! Compress to remove points where convection has terminated
          ! (mass-flux gone to zero).
          ! NOTE: it is important that the call to save_parcel_bl_top
          ! is done BEFORE this compression, in order to not miss
          ! events where the parcel terminated within the model-level
          ! straddling the BL-top.  But the saving of parcel diagnostics
          ! should be done AFTER this compression, in order to not include
          ! data where the mass-flux is zero (the computing of
          ! mass-flux-weighted means of the diagnostics over layers / types
          ! for output can give silly answers if the mass-fluxes used as
          ! the averaging weights are all zero).
          call parcel_compress( l_tracer, par_conv(i_type,i_layr) )

          if ( par_conv(i_type,i_layr) % cmpr % n_points > 0 ) then
            ! If still any points left after compressing-out terminated points

            ! Save parcel property diagnostics in conserved form,
            ! ready for computing means over types
            if ( draft_diags % par % n_diags > 0 ) then
              ! Note: in the calls below, the diagnostics are valid
              ! at half-level k_next, not level k
              ! Allocate super-array for diagnostics
              call diags_super_alloc(                                          &
                     par_conv(i_type,i_layr) % cmpr % n_points,                &
                     draft_diags % par % n_diags_super,                        &
                     draft_diags_super % par(i_type,i_layr,k_next) )
              ! Copy parcel properties for diags
              call cmpr_copy( par_conv(i_type,i_layr) % cmpr,                  &
                              draft_diags_super                                &
                              % par(i_type,i_layr,k_next) % cmpr )
              call grid_compress( grid, par_conv(i_type,i_layr) % cmpr,        &
                                  k_half=k_next,                               &
                                  grid_half_super=grid_next_super )
              call parcel_diags_copy( n_fields_tot,                            &
                     draft_diags % par,                                        &
                     par_conv(i_type,i_layr), grid_next_super(:,i_pressure),   &
                     draft_diags_super % par(i_type,i_layr,k_next) % super )
            end if

            ! Convert final parcel properties back from conserved form
            l_reverse = .true.
            call fields_k_conserved_vars(                                      &
                   par_conv(i_type,i_layr) % cmpr % n_points,                  &
                   max_points, n_fields_tot, l_reverse,                        &
                   par_conv(i_type,i_layr) % mean_super )

          end if  ! ( par_conv(i_type,i_layr) % cmpr % n_points > 0 )

        end if  ! ( par_conv(i_type,i_layr) % cmpr % n_points > 0 )
      end do  ! i_type = 1, n_conv_types
    end do  ! i_layr = 1, n_conv_layers

  end if  ! ( l_any_conv .OR. l_any_source )

end do  ! k = k_first, k_last, dk
! END MAIN LOOP OVER MODEL-LEVELS


!----------------------------------------------------------------
! 9) Optionally vertically homogenize the resolved-scale
!    source terms within the boundary-layer
!----------------------------------------------------------------

if ( l_homog_conv_bl ) then
  call homog_conv_bl_ctl( n_conv_types, n_conv_layers,                         &
                          max_points, ij_first, ij_last,                       &
                          n_fields_tot, l_down,                                &
                          grid, fields, layer_mass,                            &
                          par_bl_top, turb, res_source )
  deallocate( par_bl_top )
end if


!----------------------------------------------------------------
! 10) Check output resolved-scale source terms for NaNs etc
!----------------------------------------------------------------

! Check for bad values in the output resolved-scale source terms
if ( i_check_bad_values_cmpr > i_check_bad_none ) then
  where_string = "End of conv_sweep_ctl call for "          //                 &
                 trim(adjustl(draft_string)) // "; "        //                 &
                 "res_source"
  do k = k_bot_conv, k_top_conv
    do i_layr = 1, n_conv_layers
      do i_type = 1, n_conv_types
        if ( res_source(i_type,i_layr,k)%cmpr%n_points > 0 ) then
          call res_source_check_bad_values(                                    &
                 res_source(i_type,i_layr,k), n_fields_tot, k,                 &
                 where_string )
        end if
      end do
    end do
  end do  ! k = k_bot_conv, k_top_conv
end if  ! ( i_check_bad_values_cmpr > i_check_bad_none )


!----------------------------------------------------------------
! 12) Deallocate work arrays
!----------------------------------------------------------------

do i_layr = 1, n_conv_layers
  do i_type = 1, n_conv_types
    call cmpr_dealloc( cmpr_combined(i_type,i_layr) )
  end do
end do
call res_source_dealloc( res_source_gen )
deallocate( par_conv )


return
end subroutine conv_sweep_ctl

end module conv_sweep_ctl_mod
