! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module diags_super_type_mod

use comorph_constants_mod, only: real_cvprec
use cmpr_type_mod, only: cmpr_type

implicit none


!----------------------------------------------------------------
! Type definition containing a super-array to store
! diagnostics in compressed form on each model-level
!----------------------------------------------------------------
type :: diags_super_type
  type(cmpr_type) :: cmpr
  real(kind=real_cvprec), allocatable :: super(:,:)
end type diags_super_type


contains


!----------------------------------------------------------------
! Subroutine to allocate a diags super-array and setup its
! compression list info
!----------------------------------------------------------------
subroutine diags_super_alloc( n_points, n_diags, diags_super )

use cmpr_type_mod, only: cmpr_alloc

implicit none

! Number of points in the super-arrays
integer, intent(in) :: n_points
! Number of diags to be stored in the super-array
integer, intent(in) :: n_diags
! Structure containing diags super-array and compression
! indices to be setup
type(diags_super_type), intent(in out) :: diags_super

! Allocate the compression index arrays
call cmpr_alloc( diags_super % cmpr, n_points )

! Allocate the super-array
allocate( diags_super % super( n_points, n_diags ) )

return
end subroutine diags_super_alloc


!----------------------------------------------------------------
! Subroutine to initialise a diags super-array to zero
!----------------------------------------------------------------
subroutine diags_super_init_zero( n_diags, diags_super )

use comorph_constants_mod, only: zero

implicit none

! Number of diags to be stored in the super-array
integer, intent(in) :: n_diags
! Structure containing diags super-array
type(diags_super_type), intent(in out) :: diags_super

! Loop counters
integer :: ic, i_field

do i_field = 1, n_diags
  do ic = 1, diags_super % cmpr % n_points
    diags_super % super(ic,i_field) = zero
  end do
end do

return
end subroutine diags_super_init_zero


!----------------------------------------------------------------
! Subroutine to apply closure scaling to compressed diags
!----------------------------------------------------------------
subroutine diags_super_scaling( n_conv_types, n_conv_layers,                   &
                                ij_first, ij_last, draft_scaling,              &
                                n_diags_super, l_weight, diags_super )

use comorph_constants_mod, only: real_cvprec, k_bot_conv, k_top_conv, nx_full

implicit none

! Number of convection types
integer, intent(in) :: n_conv_types

! Number of distinct layers of convection
integer, intent(in) :: n_conv_layers

! First and last ij-indices on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Closure scaling on ij-grid, for each type / layer
real(kind=real_cvprec), intent(in) :: draft_scaling                            &
               ( ij_first:ij_last, n_conv_types, n_conv_layers )

! Number of diags in the super-array
integer, intent(in) :: n_diags_super

! List of flags indicating which diagnostic(s) in the super-array
! correspond to masses which should be scaled
logical, intent(in) :: l_weight(n_diags_super)

! Array of structures containing compressed diags super-arrays
! for each model-level, type and layer
type(diags_super_type), intent(in out) :: diags_super                          &
                        ( n_conv_types, n_conv_layers, k_bot_conv:k_top_conv )

! Loop counters
integer :: i, j, k, ic, i_type, i_layr, i_super

! Loop over diags in the super-array
do i_super = 1, n_diags_super
  if ( l_weight(i_super) ) then
    ! If this diag is a mass-weight...

    ! Loop over levels
    do k = k_bot_conv, k_top_conv
      ! Loop over convecting layers and convection types
      do i_layr = 1, n_conv_layers
        do i_type = 1, n_conv_types

          ! If any points in the diag compression list
          if ( diags_super(i_type,i_layr,k) % cmpr % n_points > 0 ) then

            do ic = 1, diags_super(i_type,i_layr,k) % cmpr % n_points
              ! Extract i,j coordinates for the current point
              i = diags_super(i_type,i_layr,k) % cmpr % index_i(ic)
              j = diags_super(i_type,i_layr,k) % cmpr % index_j(ic)

              ! Multiply mass-flux by the closure scaling value
              ! at the current ij index
              diags_super(i_type,i_layr,k) % super(ic,i_super)                 &
                = diags_super(i_type,i_layr,k) % super(ic,i_super)             &
                * draft_scaling( nx_full*(j-1)+i, i_type, i_layr)
            end do

          end if

        end do
      end do
    end do

  end if  ! ( l_weight(i_super) )
end do  ! i_super = 1, n_diags_super

return
end subroutine diags_super_scaling


!----------------------------------------------------------------
! Subroutine to compute means of compressed diagnostics over
! layers / types, and copy them into the output arrays
!----------------------------------------------------------------
subroutine diags_super_compute_means( n_conv_types, n_conv_layers,             &
                                      max_points, ij_first, ij_last,           &
                                      n_diags, n_diags_super, l_weight,        &
                                      diag_list, diags_super,                  &
                                      fields_diags1, fields_diags2 )

use comorph_constants_mod, only: zero, n_conv_layers_diag,                     &
                                 nx_full, k_bot_conv, k_top_conv
use cmpr_type_mod, only: cmpr_alloc, cmpr_copy
use diag_type_mod, only: diag_list_type
use fields_diags_type_mod, only: fields_diags_type, fields_diags_conserved_vars
use decompress_mod, only: decompress

implicit none

! Number of convection types
integer, intent(in) :: n_conv_types
! Number of vertically distinct convecting layers
integer, intent(in) :: n_conv_layers

! Max number of convecting points on current segment
integer, intent(in) :: max_points

! i,j indices of the first and last point on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Number of diags in the output list (maybe smaller than n_diags_super,
! as some extra fields maybe required in the super-array for computing
! other diags, even if not requested themselves)
integer, intent(in) :: n_diags

! Number of diags in the super-array
integer, intent(in) :: n_diags_super

! Flags for which fields in the super-array to use as weights
! for computing means of other fields
logical, intent(in) :: l_weight( n_diags_super )

! List of diags to be output (each element is a pointer to meta-data)
type(diag_list_type), intent(in out) :: diag_list(n_diags)
! Note: this is an input, but needs intent in out because this routine
! writes to the arrays referenced by the contained pointers, and some
! compilers flag this as an error because the structure is intent in.

! Array of structures containing compressed diags super-arrays
! for each type, layer and level
type(diags_super_type), intent(in out) :: diags_super                          &
     ( n_conv_types, n_conv_layers, k_bot_conv:k_top_conv )

! Structures containing meta-data for any primary fields diags
! contained in the list
type(fields_diags_type), optional, intent(in) :: fields_diags1
type(fields_diags_type), optional, intent(in) :: fields_diags2

! Array of structures storing means of diags over all
! distinct layers, for each type (on the current level)
type(diags_super_type) :: diags_super_mean1(n_conv_types)

! Means of diags over all layers and types (on the current level)
type(diags_super_type) :: diags_super_mean2

! Array bounds
integer :: lb(5), ub(5)

! Store for ij indices of active points
integer :: index_ij( max_points )
! Index lists for referencing compressed mean over all types from
! compressed data for one type only
integer :: index_ic2( max_points,                                              &
                      max( n_conv_types, n_conv_layers ) )

! Work array set to combined index at points where any
! types / layers are active
integer :: index_ic(ij_first:ij_last)

! Flag indicating whether any convecting points on current level
logical :: l_any_layr, l_any_type
! Flag indicating more than one overlapping convecting layer
logical :: l_multi_layr, l_multi_type

! Flag to pass into fields_diags_conserved_vars
logical, parameter :: l_reverse = .true.

! Loop counters
integer :: i, j, ij, ic, nc, k, i_type, i_layr, i_diag, i_super


! Allocate work arrays to max possible required size

! Arrays for means over layers for each type
do i_type = 1, n_conv_types
  allocate( diags_super_mean1(i_type) % super( max_points, n_diags_super ) )
  call cmpr_alloc( diags_super_mean1(i_type) % cmpr, max_points )
  diags_super_mean1(i_type) % cmpr % n_points = 0
end do

! Array for mean over the above for all types
allocate( diags_super_mean2 % super( max_points, n_diags_super ) )
call cmpr_alloc( diags_super_mean2 % cmpr, max_points )
diags_super_mean2 % cmpr % n_points = 0


! Loop over levels
do k = k_bot_conv, k_top_conv


  ! Initialise combined compression list sizes to zero
  do i_type = 1, n_conv_types
    diags_super_mean1(i_type) % cmpr % n_points = 0
  end do
  diags_super_mean2 % cmpr % n_points = 0

  ! Initialise flags for any / multiple convective types active
  l_any_type = .false.
  l_multi_type = .false.


  !--------------------------------------------------------------
  ! 1) For each type, compute means of diags over all layers
  !--------------------------------------------------------------

  ! For each type...
  do i_type = 1, n_conv_types

    ! Check whether any convecting points of this type are
    ! present, and whether there are multiple overlapping layers
    l_any_layr = .false.
    l_multi_layr = .false.
    do i_layr = 1, n_conv_layers
      if ( diags_super(i_type,i_layr,k) % cmpr % n_points > 0 ) then
        ! l_multi flag only set if another layer already found
        l_multi_layr = l_any_layr
        ! l_any set if any points found
        l_any_layr = .true.
      end if
    end do

    ! If any points found for current type
    if ( l_any_layr ) then
      ! Set flags for found something in current type
      l_multi_type = l_any_type
      l_any_type = .true.

      ! If multiple overlapping layers found
      if ( l_multi_layr ) then

        ! Find combined compression list of points where any
        ! layer is active for the current type and level...

        ! Initialise work array to zero
        do ij = ij_first, ij_last
          index_ic(ij) = 0
        end do

        ! Set index_ic array to 1 where any layer present
        do i_layr = 1, n_conv_layers
          if ( diags_super(i_type,i_layr,k) % cmpr % n_points > 0) then
            do ic = 1, diags_super(i_type,i_layr,k) % cmpr % n_points
              i = diags_super(i_type,i_layr,k) % cmpr % index_i(ic)
              j = diags_super(i_type,i_layr,k) % cmpr % index_j(ic)
              index_ic( nx_full*(j-1)+i ) = 1
            end do
          end if
        end do
        ! Store indices of points
        nc = 0
        do ij = ij_first, ij_last
          if ( index_ic(ij) > 0 ) then
            nc = nc + 1
            index_ij(nc) = ij
          end if
        end do
        ! Set compression indices of all found points
        diags_super_mean1(i_type) % cmpr % n_points = nc
        do ic = 1, diags_super_mean1(i_type) % cmpr % n_points
          ! Reverse ij = nx*(j-1)+i to get back i, j
          j = 1 + (index_ij(ic)-1) / nx_full
          i = index_ij(ic) - nx_full * (j-1)
          diags_super_mean1(i_type) % cmpr % index_i(ic) = i
          diags_super_mean1(i_type) % cmpr % index_j(ic) = j
          ! Store the compression list index in index_ic
          index_ic( index_ij(ic) ) = ic
        end do
        ! Store indices for referencing the combined compression
        ! list from each of the individual lists
        do i_layr = 1, n_conv_layers
          if ( diags_super(i_type,i_layr,k) % cmpr % n_points > 0) then
            do ic = 1, diags_super(i_type,i_layr,k) % cmpr % n_points
              i = diags_super(i_type,i_layr,k) % cmpr % index_i(ic)
              j = diags_super(i_type,i_layr,k) % cmpr % index_j(ic)
              index_ic2(ic,i_layr) = index_ic( nx_full*(j-1)+i )
            end do
          end if
        end do

        ! Initialise combined super-array to zero
        do i_super = 1, n_diags_super
          do ic = 1, diags_super_mean1(i_type) % cmpr % n_points
            diags_super_mean1(i_type) % super(ic,i_super) = zero
          end do
        end do

        ! Sum contributions from each layer that has any active points
        do i_layr = 1, n_conv_layers
          if ( diags_super(i_type,i_layr,k) % cmpr % n_points > 0 ) then
            call diags_super_combine( n_diags_super, l_weight,                 &
                                      index_ic2(:,i_layr),                     &
                                      diags_super(i_type,i_layr,k),            &
                                      diags_super_mean1(i_type) )
          end if
        end do

        ! Just one layer active at this level / type
      else  ! ( l_multi_layr )

        do i_layr = 1, n_conv_layers
          if ( diags_super(i_type,i_layr,k) % cmpr % n_points > 0) then

            ! Combined compression list is just a copy of the
            ! list from the one active layer
            call cmpr_copy( diags_super(i_type,i_layr,k) % cmpr,               &
                            diags_super_mean1(i_type) % cmpr )

            ! Just copy the diagnostics to form mean over layers
            do i_super = 1, n_diags_super
              do ic = 1, diags_super(i_type,i_layr,k) % cmpr % n_points
                diags_super_mean1(i_type) % super(ic,i_super)                  &
                  = diags_super(i_type,i_layr,k) % super(ic,i_super)
              end do
            end do

          end if
        end do

      end if  ! ( l_multi_layr )

    end if  ! ( l_any_layr )

  end do  ! i_type = 1, n_conv_types


  !--------------------------------------------------------------
  ! 2) Compute means of diags over all the type-means
  !--------------------------------------------------------------

  ! If any convective types are active at this level
  if ( l_any_type ) then

    ! If multiple convective types are active
    if ( l_multi_type ) then

      ! Find combined compression list of points where any
      ! type is active at the current level...

      ! Initialise work array to zero
      do ij = ij_first, ij_last
        index_ic(ij) = 0
      end do

      ! Set index_ic array to 1 where any layer present
      do i_type = 1, n_conv_types
        if ( diags_super_mean1(i_type)%cmpr % n_points > 0 ) then
          do ic = 1, diags_super_mean1(i_type) % cmpr % n_points
            i = diags_super_mean1(i_type) % cmpr % index_i(ic)
            j = diags_super_mean1(i_type) % cmpr % index_j(ic)
            index_ic( nx_full*(j-1)+i ) = 1
          end do
        end if
      end do
      ! Store indices of points
      nc = 0
      do ij = ij_first, ij_last
        if ( index_ic(ij) > 0 ) then
          nc = nc + 1
          index_ij(nc) = ij
        end if
      end do
      ! Set compression indices of all found points
      diags_super_mean2 % cmpr % n_points = nc
      do ic = 1, diags_super_mean2 % cmpr % n_points
        ! Reverse ij = nx*(j-1)+i to get back i, j
        j = 1 + (index_ij(ic)-1) / nx_full
        i = index_ij(ic) - nx_full * (j-1)
        diags_super_mean2 % cmpr % index_i(ic) = i
        diags_super_mean2 % cmpr % index_j(ic) = j
        ! Store the compression list index in index_ic
        index_ic( index_ij(ic) ) = ic
      end do
      ! Store indices for referencing the combined compression
      ! list from each of the individual lists
      do i_type = 1, n_conv_types
        if ( diags_super_mean1(i_type) % cmpr % n_points > 0) then
          do ic = 1, diags_super_mean1(i_type) % cmpr % n_points
            i = diags_super_mean1(i_type) % cmpr % index_i(ic)
            j = diags_super_mean1(i_type) % cmpr % index_j(ic)
            index_ic2(ic,i_type) = index_ic( nx_full*(j-1)+i )
          end do
        end if
      end do

      ! Initialise combined super-array to zero
      do i_super = 1, n_diags_super
        do ic = 1, diags_super_mean2 % cmpr % n_points
          diags_super_mean2 % super(ic,i_super) = zero
        end do
      end do

      ! Sum contributions from each type that has any active points
      do i_type = 1, n_conv_types
        if ( diags_super_mean1(i_type)%cmpr % n_points > 0 ) then
          call diags_super_combine( n_diags_super, l_weight,                   &
                                    index_ic2(:,i_type),                       &
                                    diags_super_mean1(i_type),                 &
                                    diags_super_mean2 )
        end if
      end do

      ! Just one type active at this level
    else  ! ( l_multi_type )

      do i_type = 1, n_conv_types
        if ( diags_super_mean1(i_type)%cmpr % n_points > 0 ) then

          ! Combined compression list is just a copy of the
          ! list from the one active type
          diags_super_mean2 % cmpr % n_points                                  &
              = diags_super_mean1(i_type) % cmpr % n_points
          call cmpr_copy( diags_super_mean1(i_type) % cmpr,                    &
                          diags_super_mean2 % cmpr )

          ! Just copy the diagnostics to form mean over layers
          do i_super = 1, n_diags_super
            do ic = 1, diags_super_mean1(i_type)%cmpr % n_points
              diags_super_mean2 % super(ic,i_super)                            &
                = diags_super_mean1(i_type) % super(ic,i_super)
            end do
          end do

        end if
      end do

    end if  ! ( l_multi_type )

  end if  ! ( l_any_type )


  ! Only anything else to do if combined compression list
  ! has at least one active point in it...
  if ( diags_super_mean2 % cmpr % n_points > 0 ) then

    !--------------------------------------------------------------
    ! 3) Convert primary fields diagnostics back from conserved
    ! variable form to normal
    !--------------------------------------------------------------

    if ( present( fields_diags1 ) ) then
      if ( fields_diags1 % n_diags > 0 ) then
        ! Conversion for fields from each type and layer
        do i_type = 1, n_conv_types
          do i_layr = 1, n_conv_layers
            if ( diags_super(i_type,i_layr,k) % cmpr % n_points > 0 ) then
              call fields_diags_conserved_vars(                                &
                     diags_super(i_type,i_layr,k) % cmpr % n_points,           &
                     diags_super(i_type,i_layr,k) % cmpr % n_points,           &
                     n_diags_super, l_reverse, fields_diags1,                  &
                     diags_super(i_type,i_layr,k) % super )
            end if
          end do
        end do
        ! Conversion for fields from each type
        do i_type = 1, n_conv_types
          if ( diags_super_mean1(i_type) % cmpr % n_points > 0 ) then
            call fields_diags_conserved_vars(                                  &
                   diags_super_mean1(i_type) % cmpr % n_points,                &
                   max_points,                                                 &
                   n_diags_super, l_reverse, fields_diags1,                    &
                   diags_super_mean1(i_type) % super )
          end if
        end do
        ! Conversion for means over all types
        call fields_diags_conserved_vars(                                      &
               diags_super_mean2 % cmpr % n_points,                            &
               max_points,                                                     &
               n_diags_super, l_reverse, fields_diags1,                        &
               diags_super_mean2 % super )
      end if  ! ( fields_diags1 % n_diags > 0 )
    end if  ! ( present( fields_diags1 ) )

    if ( present( fields_diags2 ) ) then
      if ( fields_diags2 % n_diags > 0 ) then
        ! Conversion for fields from each type and layer
        do i_type = 1, n_conv_types
          do i_layr = 1, n_conv_layers
            if ( diags_super(i_type,i_layr,k) % cmpr % n_points > 0 ) then
              call fields_diags_conserved_vars(                                &
                     diags_super(i_type,i_layr,k) % cmpr % n_points,           &
                     diags_super(i_type,i_layr,k) % cmpr % n_points,           &
                     n_diags_super, l_reverse, fields_diags2,                  &
                     diags_super(i_type,i_layr,k) % super )
            end if
          end do
        end do
        ! Conversion for fields from each type
        do i_type = 1, n_conv_types
          if ( diags_super_mean1(i_type) % cmpr % n_points > 0 ) then
            call fields_diags_conserved_vars(                                  &
                   diags_super_mean1(i_type) % cmpr % n_points,                &
                   max_points,                                                 &
                   n_diags_super, l_reverse, fields_diags2,                    &
                   diags_super_mean1(i_type) % super )
          end if
        end do
        ! Conversion for means over all types
        call fields_diags_conserved_vars(                                      &
               diags_super_mean2 % cmpr % n_points,                            &
               max_points,                                                     &
               n_diags_super, l_reverse, fields_diags2,                        &
               diags_super_mean2 % super )
      end if  ! ( fields_diags2 % n_diags > 0 )
    end if  ! ( present( fields_diags2 ) )


    !--------------------------------------------------------------
    ! 4) Copy requested diagnostics into output arrays
    !--------------------------------------------------------------

    ! For each diag in the list of those requested
    do i_diag = 1, n_diags

      ! Get super-array address
      i_super = diag_list(i_diag)%pt % i_super

      ! If diag requested as separate field for all types and layers
      if ( diag_list(i_diag)%pt % request % x_y_z_lay_typ ) then
        ! Find bounds of the output array
        lb(1:5) = lbound( diag_list(i_diag)%pt % field_5d )
        ub(1:5) = ubound( diag_list(i_diag)%pt % field_5d )
        ! Scatter data for each convection type and layer
        do i_type = 1, n_conv_types
          do i_layr = 1, min( n_conv_layers, n_conv_layers_diag )
            if ( diags_super(i_type,i_layr,k) % cmpr % n_points > 0 ) then
              call decompress( diags_super(i_type,i_layr,k) % cmpr,            &
                               diags_super(i_type,i_layr,k) % super(:,i_super),&
                               lb(1:2), ub(1:2),                               &
                               diag_list(i_diag)%pt                            &
                               % field_5d(:,:,k,i_layr,i_type) )
            end if
          end do
        end do
      end if

      ! If diag requested as separate field for each conv type
      if ( diag_list(i_diag)%pt % request % x_y_z_typ ) then
        ! Find lower-bounds of the output array
        lb(1:4) = lbound( diag_list(i_diag)%pt % field_4d )
        ub(1:4) = ubound( diag_list(i_diag)%pt % field_4d )
        ! Scatter data for each convection type
        do i_type = 1, n_conv_types
          if ( diags_super_mean1(i_type) % cmpr % n_points > 0 ) then
            call decompress( diags_super_mean1(i_type) % cmpr,                 &
                             diags_super_mean1(i_type) % super(:,i_super),     &
                             lb(1:2), ub(1:2),                                 &
                             diag_list(i_diag)%pt % field_4d(:,:,k,i_type) )
          end if
        end do
      end if

      ! If diag is requested as a mean over types
      if ( diag_list(i_diag)%pt % request % x_y_z ) then
        ! Find lower-bounds of the output array
        lb(1:3) = lbound( diag_list(i_diag)%pt % field_3d )
        ub(1:3) = ubound( diag_list(i_diag)%pt % field_3d )
        ! Scatter data
        call decompress( diags_super_mean2 % cmpr,                             &
                         diags_super_mean2 % super(:,i_super),                 &
                         lb(1:2), ub(1:2),                                     &
                         diag_list(i_diag)%pt % field_3d(:,:,k) )
      end if

    end do  ! i_diag = 1, n_diags

  end if  ! ( diags_super_mean2 % cmpr % n_points > 0 )


end do  ! k = k_bot_conv, k_top_conv


return
end subroutine diags_super_compute_means


!----------------------------------------------------------------
! Subroutine to expand compressed diagnostics within the existing
! arrays, to put them on a new, larger compression list
!----------------------------------------------------------------
subroutine diags_super_expand( n_diags_super, n_points_new, index_ic,          &
                               diags_super )

use comorph_constants_mod, only: zero

implicit none

! Number of diags in the super-array
integer, intent(in) :: n_diags_super

! Number of points in the new larger compression list
integer, intent(in) :: n_points_new

! Diagnostics super-array to be expanded
type(diags_super_type), intent(in out) :: diags_super

! Indices for transferring data from the old to the new compression list
integer, intent(in) :: index_ic( diags_super % cmpr % n_points )

! Number of points in the new compression list which will be vacant,
! and their indices
integer :: n_vacant
integer :: index_vacant(n_points_new)

! Loop counters
integer :: ic, ic2, ic_first, i_field

! See if any points need to be moved
ic_first = 0
over_n_points: do ic = 1, diags_super % cmpr % n_points
  ! Search until we find at least one point whose index changes
  ! then exit the loop.
  if ( .not. index_ic(ic) == ic ) then
    ! Save the index of the first point to move
    ic_first = ic
    exit over_n_points
  end if
end do over_n_points

! If any points need to be moved
if ( ic_first > 0 ) then

  ! For each super-array, loop through the points that need
  ! to be moved (only from ic_first onwards) and use index_ic to
  ! transfer the data to its new location.
  ! Note we need to loop backwards through the points to avoid
  ! overwriting data that hasn't been transfered yet.
  do i_field = 1, n_diags_super
    do ic = diags_super % cmpr % n_points, ic_first, -1
      diags_super % super(index_ic(ic),i_field)                                &
        = diags_super % super(ic,i_field)
    end do
  end do

end if  ! ( ic_first > 0 )

! Find the indices of points left vacant by the expansion.
n_vacant = 0
do ic2 = 1, n_points_new
  index_vacant(ic2) = 0
end do
do ic = 1, diags_super % cmpr % n_points
  index_vacant(index_ic(ic)) = 1
end do
do ic2 = 1, n_points_new
  if ( index_vacant(ic2) == 0 ) then
    n_vacant = n_vacant + 1
    index_vacant(n_vacant) = ic2
  end if
end do

! If any vacant points
if ( n_vacant > 0 ) then
  ! Zero all data in the vacant points, for safety.  They will presumably
  ! be populated with new data after the call to this routine.

  do i_field = 1, n_diags_super
    do ic = 1, n_vacant
      diags_super % super(index_vacant(ic),i_field) = zero
    end do
  end do

end if  ! ( n_vacant > 0 )

return
end subroutine diags_super_expand


!----------------------------------------------------------------
! Subroutine to add a compressed diagnostics super-array onto another one
!----------------------------------------------------------------
subroutine diags_super_combine( n_diags_super, l_weight, index_ic,             &
                                diags_super_a, diags_super_m )

use comorph_constants_mod, only: real_cvprec, one, min_float

implicit none

! Number of diags in the super-array
integer, intent(in) :: n_diags_super

! List of flags indicating which diagnostic(s) in the super-array
! correspond to masses which are summed and used as the weighting
! for computing means of subsequent diagnostics
logical, intent(in) :: l_weight(n_diags_super)

! Diags super-array to be added on
type(diags_super_type), intent(in) :: diags_super_a

! Accumulated mean diags super-array to be added onto
type(diags_super_type), intent(in out) :: diags_super_m

! Index list for referencing the diags_super_m compression list
! from diags_super_a
integer, intent(in) :: index_ic( diags_super_a % cmpr % n_points )

! Weight for computing means of diagnostics
real(kind=real_cvprec) :: weight_a( diags_super_a % cmpr % n_points )

! Loop counters
integer :: ic, ic2, i_ds

! Loop over diags in the super-array
do i_ds = 1, n_diags_super
  if ( l_weight(i_ds) ) then
    ! If this field is to be used as a mass-weight

    do ic = 1, diags_super_a % cmpr % n_points
      ic2 = index_ic(ic)
      ! Add contribution from diags_super_a to the combined total mass
      diags_super_m % super(ic2,i_ds) = diags_super_m % super(ic2,i_ds)        &
                                      + diags_super_a % super(ic,i_ds)
      ! Calculate mass-fraction weight to be applied to subsequent diagnostics
      weight_a(ic) = diags_super_a % super(ic,i_ds)                            &
              / max( diags_super_m % super(ic2,i_ds), min_float )
    end do
    do ic = 1, diags_super_a % cmpr % n_points
      ic2 = index_ic(ic)
      ! Force weight to be exactly 1.0 at points where the mass-fluxes are
      ! equal (needed for reproducibility on different decompositions).
      if ( diags_super_a%super(ic,i_ds) == diags_super_m%super(ic2,i_ds) ) then
        weight_a(ic) = one
      end if
    end do

  else  ! ( l_weight(i_ds) )
    ! Field is not a weight

    do ic = 1, diags_super_a % cmpr % n_points
      ic2 = index_ic(ic)
      ! Calculate weighted-mean
      diags_super_m % super(ic2,i_ds)                                          &
        = (one-weight_a(ic)) * diags_super_m % super(ic2,i_ds)                 &
             + weight_a(ic)  * diags_super_a % super(ic,i_ds)
    end do

  end if  ! ( l_weight(i_ds) )
end do  ! i_ds = 1, n_diags_super

return
end subroutine diags_super_combine


end module diags_super_type_mod
