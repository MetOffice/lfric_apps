! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module conv_genesis_ctl_mod

implicit none

contains

! Control subroutine which calls code to compute the initiating
! mass-sources and inititating parcel properties from each
! model-level, for both updrafts and downdrafts.
! This routine mainly handles the compression of data input
! and output by the convection initiation calculation
! (init_mass_moist_frac) to save on memory and computations.
subroutine conv_genesis_ctl( max_points, ij_first, ij_last,                    &
                             n_fields_tot, l_tracer,                           &
                             grid, turb, cloudfracs, fields,                   &
                             layer_mass, virt_temp, l_init_poss,               &
                             n_updraft_layers, updraft_par_gen,                &
                             n_dndraft_layers, dndraft_par_gen )

use comorph_constants_mod, only: real_hmprec, real_cvprec, zero,               &
                     nx_full, ny_full, k_bot_conv, k_top_conv, k_top_init,     &
                     n_updraft_types, n_dndraft_types,                         &
                     l_turb_par_gen, l_cv_cloudfrac,                           &
                     i_check_bad_values_cmpr, i_check_bad_none, name_length

use cmpr_type_mod, only: cmpr_type, cmpr_alloc, cmpr_copy
use grid_type_mod, only: grid_type, n_grid, grid_compress
use turb_type_mod, only: turb_type, n_turb
use cloudfracs_type_mod, only: cloudfracs_type, n_cloudfracs,                  &
                               i_frac_liq, i_frac_ice,                         &
                               i_frac_bulk, i_frac_precip
use fields_type_mod, only: fields_type, field_positive, i_cf_liq !, i_cf_bulk
use parcel_type_mod, only: parcel_type, parcel_alloc,                          &
                           parcel_init_zero, parcel_compress,                  &
                           parcel_check_bad_values

use compress_mod, only: compress
use force_cloudfrac_consistency_mod, only:                                     &
                                 force_cloudfrac_consistency
use set_l_within_bl_mod, only: set_l_within_bl
use init_mass_moist_frac_mod, only: init_mass_moist_frac
use par_gen_distinct_layers_mod, only: par_gen_distinct_layers

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
type(cloudfracs_type), intent(in) :: cloudfracs

! Structure containing pointers to the primary fields
type(fields_type), intent(in) :: fields

! Dry-mass on each model-level, per unit surface area
real(kind=real_hmprec), intent(in) :: layer_mass                               &
                ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Virtual temperature profile
real(kind=real_hmprec), intent(in) :: virt_temp                                &
                ( nx_full, ny_full, k_bot_conv:k_top_conv )

! 3-D mask of points where initiation mass-sources might be
! possible, for either updrafts or downdrafts
logical, intent(in) :: l_init_poss                                             &
                ( nx_full, ny_full, k_bot_conv:k_top_init )

! Counter for the max number of convection mass-source
! layers found in a column
integer, intent(out) :: n_updraft_layers
integer, intent(out) :: n_dndraft_layers

! Array of structures of output initiation mass source properties
type(parcel_type), allocatable, intent(in out) ::                              &
                                updraft_par_gen(:,:,:)
type(parcel_type), allocatable, intent(in out) ::                              &
                                dndraft_par_gen(:,:,:)
! These are out, but need intent inout so that their
! non-allocated status is known on input.


! Arrays used to initially calculate the initiating mass source
! properties, before divying up into distinct layers and copying
! into the output par_gen arrays
type(parcel_type), allocatable :: updraft_par_gen_tmp(:,:)
type(parcel_type), allocatable :: dndraft_par_gen_tmp(:,:)

! Compressed copy of layer_mass
real(kind=real_cvprec) :: layer_mass_k(max_points)
! Compressed heights and pressures (full-levels, half-levels)
real(kind=real_cvprec) :: grid_full(max_points,n_grid,-1:1)
real(kind=real_cvprec) :: grid_half(max_points,n_grid,0:1)
! Compressed virtual temperatures
real(kind=real_cvprec) :: virt_temp_cmpr(max_points,-1:1)

! Compressed primary fields from level k and the
! full levels above and below
real(kind=real_cvprec) :: fields_cmpr                                          &
                          (max_points,n_fields_tot,-1:1)

! Compressed super-arrays containing the 3 cloud-fraction fields
! cf_liq, cf_ice, cf_bulk and rain fraction, at current level
real(kind=real_cvprec) :: cloudfracs_k(max_points,n_cloudfracs)

! Compressed copy of turbulence fields at upper and lower
! model-level interface
real(kind=real_cvprec) :: turb_cmpr(max_points,n_turb,0:1)

! Compressed turbulence length-scale at current full-level
real(kind=real_cvprec) :: turb_len_k(max_points)

! Compressed par_radius_amp factor
real(kind=real_cvprec) :: par_radius_amp_cmpr(max_points)

! Flag for whether each point is below the boundary-layer top
logical :: l_within_bl(max_points)

! Flag for updraft vs downdraft calculations
logical :: l_down

! Model-level indices for neighbouring full levels
integer :: kp1, km1

! Model-level indices for neighbouring half-levels (interfaces)
integer :: kph, kmh

! Offset versions of the above for referencing multi-levelled
! compression list arrays
integer :: k_c, kp1_c, km1_c, kph_c, kmh_c

! Structure to store compression indices of points where
! initiation mass-sources are possible
type(cmpr_type) :: cmpr_init

! Array lower and upper bounds to pass into other routines
integer :: lb(3),  ub(3)
integer :: lb2(2), ub2(2)

! Description of where we are in the code, for error messages
character(len=name_length) :: where_string


! Loop counter
integer :: i, j, ij, k, k2, k2_c, ic,                                          &
           i_field, i_frac, i_type, i_layr


! Note: some calculations in the updraft / downdraft model
! (init_par_conv) assume that next-level fields exist when
! adding in the initiation mass-sources.  If mass-sources
! are allowed to initiate in the last level, these calculations
! will reference the primary fields arrays out-of-bounds,
! creating a stinking heap of garbage, and there will be
! nowhere for such mass-sources to rise/fall to anyway!

! Therefore, updrafts cannot be allowed to initiate in
! the top model-level, and downdrafts cannot be allowed
! to initiate in the bottom model-level.

! Initialise number of distinct layers of convection to zero
n_updraft_layers = 0
n_dndraft_layers = 0

! Allocate arrays of structures to store initiating parcel properties
allocate( updraft_par_gen_tmp ( n_updraft_types, k_bot_conv:k_top_init ) )
allocate( dndraft_par_gen_tmp ( n_dndraft_types, k_bot_conv:k_top_init ) )

! Initialise number of mass-source points to zero for all levels
! and types
do k = k_bot_conv, k_top_init
  do i_type = 1, n_updraft_types
    updraft_par_gen_tmp(i_type,k) % cmpr % n_points = 0
  end do
  do i_type = 1, n_dndraft_types
    dndraft_par_gen_tmp(i_type,k) % cmpr % n_points = 0
  end do
end do

! Allocate arrays for i,j indices of points where initiation
! might occur
call cmpr_alloc( cmpr_init, max_points )


! Loop over convection initiation levels...
do k = k_bot_conv, k_top_init

  ! Set level indices
  kph = k + 1
  kmh = k
  kp1 = k + 1
  km1 = k - 1

  ! Avoid subscripting k out of bounds at top and bottom
  if ( k == k_bot_conv ) then
    km1 = k
  end if
  if ( k == k_top_conv ) then
    kp1 = k
  end if

  ! Set k-indices for picking out the data in the
  ! multi-level compression arrays
  km1_c = km1 - k
  kmh_c = kmh - k
  k_c   = 0
  kph_c = kph - k
  kp1_c = kp1 - k

  ! Count number of points on current segment where initiation
  ! might occur, and store the i,j indices of the points
  cmpr_init % n_points = 0
  ! Loop over only points on current segment!
  do ij = ij_first, ij_last
    ! Reverse ij = nx*(j-1)+i to get back i, j
    j = 1 + (ij-1) / nx_full
    i = ij - nx_full * (j-1)
    if ( l_init_poss(i,j,k) ) then
      cmpr_init % n_points = cmpr_init % n_points + 1
      cmpr_init % index_i( cmpr_init % n_points ) = i
      cmpr_init % index_j( cmpr_init % n_points ) = j
    end if
  end do

  ! If any initiating points
  if ( cmpr_init % n_points > 0 ) then

    ! Store number of points in the parcel structures,
    ! allocate space for initiating parcel properties,
    ! and initialise parcel fields to zero
    ! (but can't initiate updrafts from model-top or
    !  downdrafts from model-bottom)

    if ( n_updraft_types > 0 .and. k < k_top_conv ) then
      call parcel_alloc( l_tracer, cmpr_init % n_points,                       &
                         updraft_par_gen_tmp(1,k) )
      call cmpr_copy( cmpr_init, updraft_par_gen_tmp(1,k) % cmpr )
      call parcel_init_zero( l_tracer, updraft_par_gen_tmp(1,k) )
    end if

    if ( n_dndraft_types > 0 .and. k > k_bot_conv ) then
      call parcel_alloc( l_tracer, cmpr_init % n_points,                       &
                         dndraft_par_gen_tmp(1,k) )
      call cmpr_copy( cmpr_init, dndraft_par_gen_tmp(1,k) % cmpr )
      call parcel_init_zero( l_tracer, dndraft_par_gen_tmp(1,k) )
    end if


    ! COMPRESS VARIOUS FIELDS ONTO POINTS ON THE CURRENT LEVEL
    ! where CONVECTIVE INITIATION MIGHT OCCUR, FOR COMPRESSED
    ! call to INIT_MASS_MOIST_FRAC...
    ! (these compressions also convert everything to the
    !  convection scheme's native precision).

    ! Compress layer_mass
    lb2 = [1,1]
    ub2 = [nx_full,ny_full]
    call compress( cmpr_init, lb2, ub2, layer_mass(:,:,k), layer_mass_k )

    ! Compress height and pressure from levels k-1, k, k+1
    do k2 = km1, kp1
      k2_c = k2 - k
      call grid_compress( grid, cmpr_init,                                     &
                    k=k2, grid_k_super=grid_full(:,:,k2_c) )
    end do
    ! Compress height and pressure from the half-levels
    ! (next and previous model-level interfaces)
    do k2 = kmh, kph
      k2_c = k2 - k
      call grid_compress( grid, cmpr_init,                                     &
               k_half=k2, grid_half_super=grid_half(:,:,k2_c) )
    end do


    ! Compress grid-mean primary fields from current level
    ! and the levels above and below
    do i_field = 1, n_fields_tot
      lb = lbound(fields%list(i_field)%pt)
      ub = ubound(fields%list(i_field)%pt)
      do k2 = km1, kp1
        k2_c = k2 - k
        call compress( cmpr_init, lb(1:2), ub(1:2),                            &
                       fields%list(i_field)%pt(:,:,k2),                        &
                       fields_cmpr(:,i_field,k2_c) )
      end do
      ! Remove any spurious negative values from input data
      if ( field_positive(i_field) ) then
        do k2 = km1, kp1
          k2_c = k2 - k
          do ic = 1, cmpr_init % n_points
            fields_cmpr(ic,i_field,k2_c)                                       &
              = max( fields_cmpr(ic,i_field,k2_c), zero )
          end do
        end do
      end if
    end do  ! i_field = 1, n_fields_tot

    if ( l_cv_cloudfrac ) then
      ! Cloud fractions included as primary fields
      ! Rounding errors when converting the cloud-fractions to
      ! 32-bit in compress can cause them to become slightly
      ! inconsistent; correct if needed:
      do k2 = km1, kp1
        k2_c = k2 - k
        call force_cloudfrac_consistency( cmpr_init % n_points,                &
                                          max_points,                          &
                                   fields_cmpr(:,i_cf_liq,k2_c) )
        !                          fields_cmpr(:,i_cf_liq:i_cf_bulk,k2_c) )
                ! For some reason, if the array fields_cmpr is
                ! subscripted in the argument list using subsetting,
                ! ifort insists on making an array temporary,
                ! even though the array section is contiguous!
                ! So using old-school syntax of just specifying the
                ! first field index but implicitly passing in all 3
                ! cloud-fraction fields because they are adjacent
                ! in memory and 3 fields are required to fill the
                ! dummy argument inside.
      end do
    end if

    ! Set flag for whether level k is within the
    ! boundary-layer at each point.
    lb = lbound( grid % height_half )
    ub = ubound( grid % height_half )
    lb2 = lbound( turb % z_bl_top )
    ub2 = ubound( turb % z_bl_top )
    call set_l_within_bl( cmpr_init,                                           &
                          k, lb, ub, grid % height_half,                       &
                          lb2, ub2, turb % z_bl_top,                           &
                          l_within_bl )

    ! If using turbulence-based parcel perturbations
    if ( l_turb_par_gen ) then
      ! Compress turbulence fields.
      ! Note: these are defined on rho-levels, and we extract
      ! them from the rho-levels kph and kmh
      do i_field = 1, n_turb
        lb = lbound(turb%list(i_field)%pt)
        ub = ubound(turb%list(i_field)%pt)
        do k2 = kmh, kph
          k2_c = k2 - k
          call compress( cmpr_init, lb(1:2), ub(1:2),                          &
                         turb%list(i_field)%pt(:,:,k2),                        &
                         turb_cmpr(:,i_field,k2_c) )
        end do
      end do
      ! Also compress turbulence length-scale on full-levels
      lb = lbound(turb%lengthscale)
      ub = ubound(turb%lengthscale)
      call compress( cmpr_init, lb(1:2), ub(1:2),                              &
                     turb%lengthscale(:,:,k),                                  &
                     turb_len_k )
      ! And par_radius_amp
      lb2 = lbound(turb%par_radius_amp)
      ub2 = ubound(turb%par_radius_amp)
      call compress( cmpr_init, lb2, ub2,                                      &
                     turb%par_radius_amp(:,:),                                 &
                     par_radius_amp_cmpr )
    end if  ! ( l_turb_par_gen )

    ! Compress virtual temperatures
    lb2 = [1,1]
    ub2 = [nx_full,ny_full]
    do k2 = km1, kp1
      k2_c = k2 - k
      call compress( cmpr_init, lb2, ub2,                                      &
                     virt_temp(:,:,k2), virt_temp_cmpr(:,k2_c) )
    end do

    if ( l_cv_cloudfrac ) then
      ! Cloud fractions included as primary fields

      ! Copy cloud fractions into their own arrays
      do i_frac = i_frac_liq, i_frac_bulk
        i_field = i_cf_liq-1 + i_frac
        do ic = 1, cmpr_init % n_points
          cloudfracs_k(ic,i_frac) = fields_cmpr(ic,i_field,k_c)
        end do
      end do

    else  ! ( l_cv_cloudfrac )
      ! Cloud fractions only passed in as diagnostic fields

      ! Compress cloud fractions from diagnostic input fields
      lb = lbound(cloudfracs % frac_liq)
      ub = ubound(cloudfracs % frac_liq)
      call compress( cmpr_init, lb(1:2), ub(1:2),                              &
                     cloudfracs % frac_liq(:,:,k),                             &
                     cloudfracs_k(:,i_frac_liq) )
      lb = lbound(cloudfracs % frac_ice)
      ub = ubound(cloudfracs % frac_ice)
      call compress( cmpr_init, lb(1:2), ub(1:2),                              &
                     cloudfracs % frac_ice(:,:,k),                             &
                     cloudfracs_k(:,i_frac_ice) )
      lb = lbound(cloudfracs % frac_bulk)
      ub = ubound(cloudfracs % frac_bulk)
      call compress( cmpr_init, lb(1:2), ub(1:2),                              &
                     cloudfracs % frac_bulk(:,:,k),                            &
                     cloudfracs_k(:,i_frac_bulk) )

      do i_frac = i_frac_liq, i_frac_bulk
        ! Remove spurious negative values
        do ic = 1, cmpr_init % n_points
          cloudfracs_k(ic,i_frac)                                              &
            = max( cloudfracs_k(ic,i_frac), zero )
        end do
      end do

      ! Rounding errors when converting the cloud-fractions to
      ! 32-bit in compress can cause them to become slightly
      ! inconsistent; correct if needed:
      call force_cloudfrac_consistency( cmpr_init % n_points,                  &
                                        max_points,                            &
               cloudfracs_k(:,i_frac_liq:i_frac_bulk) )

    end if  ! ( l_cv_cloudfrac )

    lb = lbound(cloudfracs % frac_precip)
    ub = ubound(cloudfracs % frac_precip)
    ! Compress rain fractions
    call compress( cmpr_init, lb(1:2), ub(1:2),                                &
                   cloudfracs % frac_precip(:,:,k),                            &
                   cloudfracs_k(:,i_frac_precip) )
    ! Remove spurious negative values
    do ic = 1, cmpr_init % n_points
      cloudfracs_k(ic,i_frac_precip)                                           &
        = max( cloudfracs_k(ic,i_frac_precip), zero )
    end do

    ! Compute initiating parcel mass-flux and properties
    call init_mass_moist_frac(                                                 &
           cmpr_init % n_points, max_points,                                   &
           l_tracer, n_fields_tot, cmpr_init, k,                               &
           l_within_bl, layer_mass_k, turb_len_k, par_radius_amp_cmpr,         &
           turb_cmpr(:,:,kmh_c), turb_cmpr(:,:,kph_c),                         &
           grid_full(:,:,km1_c), grid_half(:,:,kmh_c),                         &
           grid_full(:,:,k_c),                                                 &
           grid_half(:,:,kph_c), grid_full(:,:,kp1_c),                         &
           fields_cmpr(:,:,km1_c), fields_cmpr(:,:,k_c),                       &
           fields_cmpr(:,:,kp1_c), cloudfracs_k,                               &
           virt_temp_cmpr(:,km1_c), virt_temp_cmpr(:,k_c),                     &
           virt_temp_cmpr(:,kp1_c),                                            &
           updraft_par_gen_tmp(1,k), dndraft_par_gen_tmp(1,k) )

    ! Check whether any points in the compression list
    ! have come out with zero mass source;
    ! if they have, recompress to remove them
    call parcel_compress( l_tracer, updraft_par_gen_tmp(1,k) )
    call parcel_compress( l_tracer, dndraft_par_gen_tmp(1,k) )


  end if  ! ( cmpr_init % n_points > 0 )

end do  ! k = k_bot_conv, k_top_init


! Divvy-up the initiating parcel data among different
! convection types

! Not yet implemented.


! Divvy-up the initiating parcel data into distinct
! convecting layers, copied into the output arrays in par_gen,
! and deallocate the par_gen_tmp structures
if ( n_updraft_types > 0 ) then
  l_down = .false.
  call par_gen_distinct_layers( l_down, l_tracer,                              &
                                ij_first, ij_last,                             &
                                n_updraft_types,n_updraft_layers,              &
                                updraft_par_gen_tmp,                           &
                                updraft_par_gen )
end if
deallocate( updraft_par_gen_tmp )
if ( n_dndraft_types > 0 ) then
  l_down = .true.
  call par_gen_distinct_layers( l_down, l_tracer,                              &
                                ij_first, ij_last,                             &
                                n_dndraft_types,n_dndraft_layers,              &
                                dndraft_par_gen_tmp,                           &
                                dndraft_par_gen )
end if
deallocate( dndraft_par_gen_tmp )


! If no convecting layers found, allocate the output par_gen
! arrays to minimal size (avoids error from passing an
! unallocated array into various routines from comorph_main).
if ( n_updraft_layers == 0 ) allocate( updraft_par_gen(1,1,1) )
if ( n_dndraft_layers == 0 ) allocate( dndraft_par_gen(1,1,1) )


! Check output initiating mass-source properties for bad values
if ( i_check_bad_values_cmpr > i_check_bad_none ) then

  if ( n_updraft_layers > 0 ) then
    where_string = "End of conv_genesis_ctl call; "           //               &
                   "updraft_par_gen"
    do k = k_bot_conv, k_top_conv
      do i_layr = 1, n_updraft_layers
        do i_type = 1, n_updraft_types
          if ( updraft_par_gen(i_type,i_layr,k)                                &
               % cmpr % n_points > 0 ) then
            call parcel_check_bad_values(                                      &
                   updraft_par_gen(i_type,i_layr,k), n_fields_tot,             &
                   k, where_string )
          end if
        end do
      end do
    end do  ! k = k_bot_conv, k_top_conv
  end if  ! ( n_updraft_layers > 0 )

  if ( n_dndraft_layers > 0 ) then
    where_string = "End of conv_genesis_ctl call; "           //               &
                   "dndraft_par_gen"
    do k = k_bot_conv, k_top_conv
      do i_layr = 1, n_dndraft_layers
        do i_type = 1, n_dndraft_types
          if ( dndraft_par_gen(i_type,i_layr,k)                                &
               % cmpr % n_points > 0 ) then
            call parcel_check_bad_values(                                      &
                   dndraft_par_gen(i_type,i_layr,k), n_fields_tot,             &
                   k, where_string )
          end if
        end do
      end do
    end do  ! k = k_bot_conv, k_top_conv
  end if  ! ( n_dndraft_layers > 0 )

end if  ! ( i_check_bad_values_cmpr > i_check_bad_none )


return
end subroutine conv_genesis_ctl


end module conv_genesis_ctl_mod
