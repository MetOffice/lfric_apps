! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module conv_incr_ctl_mod

implicit none

contains

! Subroutine to calculate the updated primary model fields
! after convection.  This involves adding contributions from:
!
! a) Entrainment and detrainment of air with different properties
!    to the grid-mean, transfer of hydrometeors by precipitation
!    fall, and transfer of momentum by pressure forces
!    (these terms are all stored in the res_source structures)
!
! b) Vertical advection due to compensating rearrangement of
!    mass within the column, in response to the imposed net
!    sub-grid vertical mass transport by the convective updrafts
!    and downdrafts.
!
! Note that Kuell & Bott (2008) proposed that the latter step (b)
! should be handled by the host model's dynamical core
! rather than by the convection scheme.
! This option is catered for by separating the computations
! for (b) from those for (a), so that (b) can easily be omitted
! here if desired, by setting the switch l_column_mass_rearrange
! to false.

subroutine conv_incr_ctl( n_updraft_layers, n_dndraft_layers,                  &
                          max_points, ij_first, ij_last,                       &
                          n_fields_tot, cmpr_any,                              &
                          updraft_res_source,                                  &
                          updraft_fallback_res_source,                         &
                          dndraft_res_source,                                  &
                          dndraft_fallback_res_source,                         &
                          lb_p, ub_p, pressure,                                &
                          layer_mass, fields, cloudfracs,                      &
                          comorph_diags )

use comorph_constants_mod, only: real_cvprec, real_hmprec, zero,               &
                     nx_full, ny_full, k_bot_conv, k_top_conv,                 &
                     n_updraft_types, n_dndraft_types,                         &
                     l_updraft_fallback, l_dndraft_fallback,                   &
                     l_column_mass_rearrange, l_cv_cloudfrac,                  &
                     i_convcloud, i_convcloud_none

use cmpr_type_mod, only: cmpr_type
use res_source_type_mod, only: res_source_type
use fields_type_mod, only: fields_type, ragged_super_type,                     &
                           fields_k_conserved_vars,                            &
                           i_cf_liq, i_cf_bulk
use cloudfracs_type_mod, only: cloudfracs_type
use comorph_diags_type_mod, only: comorph_diags_type

use copy_field_mod, only: copy_field_cmpr
use diff_field_mod, only: diff_field_cmpr
use compress_mod, only: compress
use decompress_mod, only: decompress
use add_res_source_mod, only: add_res_source
use mass_rearrange_mod, only: mass_rearrange
use calc_diag_conv_cloud_mod, only: calc_diag_conv_cloud,                      &
                                    zero_diag_conv_cloud

implicit none

! Highest number of distinct updraft and downdraft layers
! occuring on the current segment
integer, intent(in) :: n_updraft_layers
integer, intent(in) :: n_dndraft_layers

! Max size the compression lists can possibly need
! (i.e. the number of columns in the current segment)
integer, intent(in) :: max_points

! Position indices ( nx*(j-1) + i ) of the first and last
! point on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Total number of fields to loop over (depends on whether or
! not this call is updating the passive tracers)
integer, intent(in) :: n_fields_tot

! Structure storing compression list indices for the points where
! any kind of convection is occuring
type(cmpr_type), intent(in) :: cmpr_any( k_bot_conv:k_top_conv )

! Arrays of structures containing resolved-scale source terms
! in separate compression lists for updrafts & downdrafts
! and their fall-back flows
type(res_source_type), intent(in) :: updraft_res_source (:,:,:)
type(res_source_type), intent(in) :: updraft_fallback_res_source               &
                                                        (:,:,:)
type(res_source_type), intent(in) :: dndraft_res_source (:,:,:)
type(res_source_type), intent(in) :: dndraft_fallback_res_source               &
                                                        (:,:,:)
! These are allocatable arrays, and we don't know for certain
! what shape they will have been allocated to.
! If n_updraft_layers > 0, shape is
! ( n_updraft_types, n_updraft_layers, k_bot_conv:k_top_conv )
! But if n_updraft_layers = 0, shape is
! ( 1, 1, k_bot_conv:k_top_conv )

! Full 3-D array of pressure
! (and lower and upper bounds of the array, in case it has halos)
integer, intent(in) :: lb_p(3), ub_p(3)
real(kind=real_hmprec), intent(in) :: pressure                                 &
                          ( lb_p(1):ub_p(1), lb_p(2):ub_p(2), lb_p(3):ub_p(3) )

! Full 3-D field of layer-masses (not modified by convection)
real(kind=real_hmprec), intent(in) :: layer_mass                               &
                     ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Structure containing pointers to the 3-D primary fields
type(fields_type), intent(in out) :: fields

! Structure containing pointers to the output full 3-D
! convective cloud fields.
type(cloudfracs_type), intent(in out) :: cloudfracs


! Structure containing pointers to output diagnostics and
! miscellaneous outputs required elsewhere in the model
type(comorph_diags_type), intent(in out) :: comorph_diags


! Super-array storing updated grid-mean fields on current
! model-level on the cmpr_any compression list
real(kind=real_cvprec) :: fields_k_super                                       &
                          ( max_points, n_fields_tot )

! Layer-mass on current level on the cmpr_any compression list
real(kind=real_cvprec) :: layer_mass_k ( max_points )

! Separate compressed fields on each level and layer-mass,
! needed for mass-rearrangement calculation
type(ragged_super_type), allocatable :: fields_cmpr(:)
real(kind=real_cvprec), allocatable :: layer_mass_k2_added(:)

! k-index of other layer that fields from level k are being
! scattered into, used by mass-rearrangement calculation
integer, allocatable :: k2(:)

! Flag input to fileds_k_conserved_vars to reverse the conversion
logical :: l_reverse

! Array lower and upper bounds
integer :: lb(3),   ub(3)
integer :: lb_f(3), ub_f(3)  ! For primary field
integer :: lb_d(3), ub_d(3)  ! For diagnostic

! cmpr_any compression indices scattered into a common
! array-space for referencing from the individual convection
! type compression lists
integer :: index_ic( ij_first : ij_last,                                       &
                     k_bot_conv-1 : k_top_conv+1 )

! Loop counters
integer :: i, j, ij, ic, k, i_field, i_diag


!----------------------------------------------------------------
! 1) Initialisations
!----------------------------------------------------------------

! Initialise index_ic to zero
do k = k_bot_conv-1, k_top_conv+1
  do ij = ij_first, ij_last
    index_ic(ij,k) = 0
  end do
end do
! Scatter cmpr_any indices into index_ic
do k = k_bot_conv, k_top_conv
  do ic = 1, cmpr_any(k) % n_points
    i = cmpr_any(k) % index_i(ic)
    j = cmpr_any(k) % index_j(ic)
    ij = nx_full*(j-1)+i
    index_ic(ij,k) = ic
    ! Now we can find the correct address in the cmpr_any
    ! compression list from any i,j,k point, by computing ij
    ! and looking up ic in index_ic
  end do
end do


! If doing mass rearrangement, we need additional work variables:
if ( l_column_mass_rearrange ) then

  ! Compressed copies of fields on all levels
  allocate( fields_cmpr(k_bot_conv:k_top_conv) )
  do k = k_bot_conv, k_top_conv
    if ( cmpr_any(k) % n_points > 0 ) then
      ! Allocate work arrays
      allocate( fields_cmpr(k) % super( cmpr_any(k) % n_points,                &
                                        n_fields_tot ) )
      ! Initialise to zero
      do i_field = 1, n_fields_tot
        do ic = 1, cmpr_any(k) % n_points
          fields_cmpr(k) % super(ic,i_field) = zero
        end do
      end do
    end if
  end do

  ! Model-level that fields from level k are being scattered into,
  ! and amount of mass so-far added to that level
  allocate( k2(ij_first:ij_last) )
  allocate( layer_mass_k2_added(ij_first:ij_last) )
  do ij = ij_first, ij_last
    k2(ij) = 0
    layer_mass_k2_added(ij) = zero
  end do

end if  ! ( l_column_mass_rearrange )


!----------------------------------------------------------------
! 2) Save fields before updating, for convection increment diags
!----------------------------------------------------------------

if ( comorph_diags % fields_incr % n_diags > 0 ) then

  ! Copy the fields before convection into the field
  ! increment diagnostic arrays, before we modify them.
  ! This way, we can work out the total convective increments
  ! afterwards as:
  ! fields_incr = fields - fields_incr

  ! Loop over diagnostics list
  do i_diag = 1, comorph_diags % fields_incr % n_diags
    i_field = comorph_diags % fields_incr % list(i_diag)%pt                    &
              % i_field
    lb_f = lbound( fields % list(i_field)%pt )
    ub_f = ubound( fields % list(i_field)%pt )
    lb_d = lbound( comorph_diags % fields_incr % list(i_diag)%pt % field_3d )
    ub_d = ubound( comorph_diags % fields_incr % list(i_diag)%pt % field_3d )
    do k = k_bot_conv, k_top_conv
      if ( cmpr_any(k) % n_points > 0 ) then
        call copy_field_cmpr( cmpr_any(k),                                     &
               lb_f(1:2), ub_f(1:2), fields % list(i_field)%pt(:,:,k),         &
               lb_d(1:2), ub_d(1:2), comorph_diags % fields_incr               &
                                     % list(i_diag)%pt%field_3d(:,:,k) )
      end if
    end do
  end do

end if  ! ( comorph_diags % fields_incr % n_diags > 0 )


! MAIN LOOP OVER LEVELS
! To add resolved-scale source terms onto the fields...
do k = k_bot_conv, k_top_conv
  if ( cmpr_any(k) % n_points > 0 ) then


    !------------------------------------------------------------
    ! 3) Add resolved-scale source terms onto compressed copies
    !    of the fields on the current level k
    !------------------------------------------------------------

    ! Compress the fields into the super-array on each level
    do i_field = 1, n_fields_tot
      lb = lbound(fields%list(i_field)%pt)
      ub = ubound(fields%list(i_field)%pt)
      call compress( cmpr_any(k), lb(1:2), ub(1:2),                            &
                     fields%list(i_field)%pt(:,:,k),                           &
                     fields_k_super(:,i_field) )
    end do

    ! Compress layer-masses
    lb = [1,1,k_bot_conv]
    ub = [nx_full,ny_full,k_top_conv]
    call compress( cmpr_any(k), lb(1:2), ub(1:2),                              &
                   layer_mass(:,:,k), layer_mass_k )


    ! Convert compressed fields to variables which are
    ! conserved following dry-mass
    l_reverse = .false.
    call fields_k_conserved_vars( cmpr_any(k)%n_points,                        &
                                  max_points, n_fields_tot,                    &
                                  l_reverse, fields_k_super )


    ! Scale the fields by layer-masses
    do i_field = 1, n_fields_tot
      do ic = 1, cmpr_any(k) % n_points
        fields_k_super(ic,i_field) = fields_k_super(ic,i_field)                &
                                   * layer_mass_k(ic)
      end do
    end do


    ! Add on the contributions from resolved-scale source
    ! terms for all convective drafts, types and layers...

    if ( n_updraft_types > 0 .and. n_updraft_layers > 0 ) then
      ! Contribution from primary updrafts
      call add_res_source( max_points, n_fields_tot,                           &
                           n_updraft_types, n_updraft_layers,                  &
                           ij_first, ij_last, index_ic(:,k),                   &
                           updraft_res_source(:,:,k),                          &
                           fields_k_super, layer_mass_k )
      ! Contribution from updraft fall-backs
      if ( l_updraft_fallback ) then
        call add_res_source( max_points, n_fields_tot,                         &
                             n_updraft_types, n_updraft_layers,                &
                             ij_first, ij_last, index_ic(:,k),                 &
                             updraft_fallback_res_source(:,:,k),               &
                             fields_k_super, layer_mass_k )
      end if
    end if

    if ( n_dndraft_types > 0 .and. n_dndraft_layers > 0 ) then
      ! Contribution from primary downdrafts
      call add_res_source( max_points, n_fields_tot,                           &
                           n_dndraft_types, n_dndraft_layers,                  &
                           ij_first, ij_last, index_ic(:,k),                   &
                           dndraft_res_source(:,:,k),                          &
                           fields_k_super, layer_mass_k )
      ! Contribution from downdraft fall-backs
      if ( l_dndraft_fallback ) then
        call add_res_source( max_points, n_fields_tot,                         &
                             n_dndraft_types, n_dndraft_layers,                &
                             ij_first, ij_last, index_ic(:,k),                 &
                             dndraft_fallback_res_source(:,:,k),               &
                             fields_k_super, layer_mass_k )
      end if
    end if

    ! Divide the fields by the modified layer-mass
    ! (i.e. normalise to get the mean of each field, combining
    !  the pre-existing and entrained/detrained properties)
    do i_field = 1, n_fields_tot
      do ic = 1, cmpr_any(k) % n_points
        fields_k_super(ic,i_field) = fields_k_super(ic,i_field)                &
                                   / layer_mass_k(ic)
      end do
    end do

    ! Occasionally, the resolved-scale source-terms can create
    ! small negative values for cloud-fractions, since the CFL
    ! check in the closure doesn't explicitly limit the mass-flux
    ! to avoid this.  Remove negative values here.
    if ( l_cv_cloudfrac ) then
      do i_field = i_cf_liq, i_cf_bulk
        do ic = 1, cmpr_any(k) % n_points
          fields_k_super(ic,i_field) = max( zero,                              &
                       fields_k_super(ic,i_field) )
        end do
      end do
    end if


    !------------------------------------------------------------
    ! 4) Scatter updated fields at k back to grid
    !------------------------------------------------------------

    ! If doing mass-rearrangement
    if ( l_column_mass_rearrange ) then

      ! Call routine to scatter the
      ! updated fields onto the level(s) where they end up,
      ! not nessecarily the curent level k
      call mass_rearrange( max_points, n_fields_tot,                           &
                           ij_first, ij_last, index_ic,                        &
                           cmpr_any(k), lb_p, ub_p, pressure,                  &
                           fields_k_super, layer_mass, k, k2,                  &
                           layer_mass_k, layer_mass_k2_added,                  &
                           fields_cmpr, comorph_diags )

      ! Not doing mass rearrangement
    else  ! ( l_column_mass_rearrange )

      ! Convert conserved variables back to normal variables
      l_reverse = .true.
      call fields_k_conserved_vars( cmpr_any(k)%n_points,                      &
                                    max_points, n_fields_tot,                  &
                                    l_reverse, fields_k_super )

      if ( i_convcloud > i_convcloud_none ) then
        ! Update diagnosed convective cloud-fractions
        call calc_diag_conv_cloud( n_updraft_layers,                           &
                                   n_dndraft_layers,                           &
                                   max_points, n_fields_tot,                   &
                                   cmpr_any, k,                                &
                                   ij_first, ij_last, index_ic,                &
                                   updraft_res_source,                         &
                                   updraft_fallback_res_source,                &
                                   dndraft_res_source,                         &
                                   dndraft_fallback_res_source,                &
                                   lb_p, ub_p, pressure,                       &
                                   fields_k_super, cloudfracs )
      end if

      ! Scatter updated values back to the 3-D arrays
      do i_field = 1, n_fields_tot
        lb = lbound(fields%list(i_field)%pt)
        ub = ubound(fields%list(i_field)%pt)
        call decompress( cmpr_any(k),                                          &
                         fields_k_super(:,i_field),                            &
                         lb(1:2), ub(1:2), fields%list(i_field)%pt(:,:,k) )
      end do

    end if  ! ( l_column_mass_rearrange )

  end if  ! ( cmpr_any(k) % n_points > 0 )
end do  ! k = k_bot_conv, k_top_conv
! end OF MAIN LOOP OVER LEVELS


!----------------------------------------------------------------
! 5) Final decompression after mass rearrangement calculation
!----------------------------------------------------------------

! If doing mass-rearrangement, the updated 3-D fields are now
! stored in compressed form in fields_cmpr.
if ( l_column_mass_rearrange ) then

  do k = k_bot_conv, k_top_conv
    if ( cmpr_any(k) % n_points > 0 ) then

      ! Convert back to normal variables
      l_reverse = .true.
      call fields_k_conserved_vars( cmpr_any(k)%n_points,                      &
                                    cmpr_any(k)%n_points,                      &
                                    n_fields_tot, l_reverse,                   &
                                    fields_cmpr(k) % super )

      if ( i_convcloud > i_convcloud_none ) then
        ! Update diagnosed convective cloud-fractions
        call calc_diag_conv_cloud( n_updraft_layers,                           &
                                   n_dndraft_layers,                           &
                                   cmpr_any(k) % n_points,                     &
                                   n_fields_tot, cmpr_any, k,                  &
                                   ij_first, ij_last, index_ic,                &
                                   updraft_res_source,                         &
                                   updraft_fallback_res_source,                &
                                   dndraft_res_source,                         &
                                   dndraft_fallback_res_source,                &
                                   lb_p, ub_p, pressure,                       &
                                   fields_cmpr(k) % super,                     &
                                   cloudfracs )
      end if

      ! Scatter updated values back to the 3-D arrays
      do i_field = 1, n_fields_tot
        lb = lbound(fields%list(i_field)%pt)
        ub = ubound(fields%list(i_field)%pt)
        call decompress( cmpr_any(k),                                          &
                         fields_cmpr(k) % super(:,i_field),                    &
                         lb(1:2), ub(1:2), fields%list(i_field)%pt(:,:,k) )
      end do

    end if
  end do

end if


!----------------------------------------------------------------
! 6) Convective cloud calculations at non-convecting points
!    (currently just sets them to zero)
!----------------------------------------------------------------
if ( i_convcloud > i_convcloud_none ) then
  do k = k_bot_conv, k_top_conv
    call zero_diag_conv_cloud( ij_first, ij_last, k, cmpr_any,                 &
                               cloudfracs )
  end do
end if


!----------------------------------------------------------------
! 7) Compute fields after-before, for convection increment diags
!----------------------------------------------------------------

if ( comorph_diags % fields_incr % n_diags > 0 ) then

  ! fields_incr currently stores saved values from before
  ! being updated by convection, so
  ! fields_incr = fields - fields_incr

  ! Loop over diagnostics list
  do i_diag = 1, comorph_diags % fields_incr % n_diags
    i_field = comorph_diags % fields_incr % list(i_diag)%pt                    &
              % i_field
    lb_f = lbound( fields % list(i_field)%pt )
    ub_f = ubound( fields % list(i_field)%pt )
    lb_d = lbound( comorph_diags % fields_incr % list(i_diag)%pt % field_3d )
    ub_d = ubound( comorph_diags % fields_incr % list(i_diag)%pt % field_3d )
    do k = k_bot_conv, k_top_conv
      if ( cmpr_any(k) % n_points > 0 ) then
        call diff_field_cmpr( cmpr_any(k),                                     &
               lb_f(1:2), ub_f(1:2), fields % list(i_field)%pt(:,:,k),         &
               lb_d(1:2), ub_d(1:2), comorph_diags % fields_incr               &
                                     % list(i_diag)%pt%field_3d(:,:,k) )

      end if
    end do
  end do

end if  ! ( comorph_diags % fields_incr % n_diags > 0 )

! Environment pressure increment diagnostic;
! Lagrangian pressure change
! dp = pressure at k - source-layer pressure (stored in diag)
if ( comorph_diags % pressure_incr_env % flag .and.                            &
     l_column_mass_rearrange ) then
  lb_d = lbound( comorph_diags % pressure_incr_env % field_3d )
  ub_d = ubound( comorph_diags % pressure_incr_env % field_3d )
  do k = k_bot_conv, k_top_conv
    if ( cmpr_any(k) % n_points > 0 ) then
      call diff_field_cmpr( cmpr_any(k),                                       &
               lb_p(1:2), ub_p(1:2), pressure(:,:,k),                          &
               lb_d(1:2), ub_d(1:2), comorph_diags % pressure_incr_env         &
                                     % field_3d(:,:,k) )
    end if
  end do
end if


! Deallocate work arrays
if ( l_column_mass_rearrange ) then
  deallocate( layer_mass_k2_added )
  deallocate( k2 )
  deallocate( fields_cmpr )
end if


return
end subroutine conv_incr_ctl


end module conv_incr_ctl_mod
