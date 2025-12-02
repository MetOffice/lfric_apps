! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module par_gen_distinct_layers_mod

implicit none

contains


! Subroutine to identify distinct layers of convective
! mass sources (i.e. multiple adjacent layers in the same
! column with nonzero mass source),
! and rearrange the initial parcel data so that each distinct
! layer is in its own compression list
subroutine par_gen_distinct_layers( l_down, l_tracer,                          &
                                    ij_first, ij_last,                         &
                                    n_conv_types, n_conv_layers,               &
                                    par_gen_tmp, par_gen )

use comorph_constants_mod, only: nx_full, k_bot_conv, k_top_conv, k_top_init
use parcel_type_mod, only: parcel_type, parcel_alloc, parcel_copy

implicit none

! Flag for whether testing downwards
! (i.e. calculating mass sources for downdrafts)
! as opposed to testing upwards for updrafts.
logical, intent(in) :: l_down

! Flag for whether tracers are used
logical, intent(in) :: l_tracer

! i,j indices of the first and last point on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Number of convection types
integer, intent(in) :: n_conv_types

! Counter for the max number of convection mass-source
! layers found in a column
integer, intent(out) :: n_conv_layers

! Arrays used to initially calculate the initiating mass source
! properties, before divying up into distinct layers
type(parcel_type), intent(in) :: par_gen_tmp                                   &
                         ( n_conv_types, k_bot_conv:k_top_init )
! Array of structures of output initiation mass source properties
type(parcel_type), allocatable, intent(in out) :: par_gen(:,:,:)

! Array storing current convective layer count in each column
integer :: i_layr_ij(ij_first:ij_last)

! Ragged array storing the compressed layer indices at each
! model-level
type :: i_layr_type
  integer, allocatable :: i(:)
end type i_layr_type
type(i_layr_type) :: i_layr_cmpr( n_conv_types,                                &
                                  k_bot_conv:k_top_conv )

! Array storing flags for mass sources found at current level
logical :: l_mass_source(ij_first:ij_last)

! Loop bounds and stride to use in the level-loop
integer :: k_first, k_last, dk

! k-index of next model-level
integer :: k_next

! Flag for whether the data on one model-level needs to be split
! between multiple different distinct convecting layers
logical :: l_multi_layer

! Compression list info for moving data into separate
! layer arrays
integer :: n_points
integer, allocatable :: index_ic(:)

! Loop counters
integer :: i, j, ij, k, ic, i_layr, i_type


! Initialise flag for mass-sources found at each point
do ij = ij_first, ij_last
  l_mass_source(ij) = .false.
end do

! Set vertical level loop indices based on whether
! going up or down
if ( l_down ) then
  k_first = k_top_init
  k_last = k_bot_conv
  dk = -1
else
  k_first = k_bot_conv
  k_last = k_top_init
  dk = 1
end if

! Initialise number of convecting layers to zero
n_conv_layers = 0

! Loop over convection types
do i_type = 1, n_conv_types

  ! Initialise convecting layer counter
  do ij = ij_first, ij_last
    i_layr_ij(ij) = 1
  end do

  ! First loop over levels to find the layer index for each
  ! initiating point.
  do k = k_first, k_last, dk

    ! If any mass sources from this level
    if ( par_gen_tmp(i_type,k) % cmpr % n_points > 0 ) then

      ! Set k-index of next model-level
      k_next = k + dk

      ! Allocate layer indices for this level
      allocate( i_layr_cmpr(i_type,k) % i                                      &
                ( par_gen_tmp(i_type,k) % cmpr % n_points ) )

      ! Set layer indices for current level and set flag
      do ic = 1, par_gen_tmp(i_type,k) % cmpr % n_points
        i = par_gen_tmp(i_type,k) % cmpr % index_i(ic)
        j = par_gen_tmp(i_type,k) % cmpr % index_j(ic)
        ij = (j-1)*nx_full + i
        ! Set layer index for this point
        i_layr_cmpr(i_type,k)%i(ic) = i_layr_ij(ij)
        ! Set flag
        l_mass_source(ij) = .true.
      end do

      ! If any mass sources from the next level
      ! (can't have any if we're at the last initiating level!)
      if ( .not. k==k_last ) then
        if ( par_gen_tmp(i_type,k_next) % cmpr % n_points > 0 ) then
          ! Reset flag to false at columns where there is also
          ! a mass source at the next level
          do ic = 1, par_gen_tmp(i_type,k_next) % cmpr % n_points
            i = par_gen_tmp(i_type,k_next) % cmpr % index_i(ic)
            j = par_gen_tmp(i_type,k_next) % cmpr % index_j(ic)
            ij = (j-1)*nx_full + i
            l_mass_source(ij) = .false.
          end do
        end if
      end if

      ! Any points where l_mass_source is still set to true
      ! must have mass sources at the current level but not
      ! at the next level.
      ! At these points, increment the layer counter
      ! ready for the next layer to be found
      do ic = 1, par_gen_tmp(i_type,k) % cmpr % n_points
        i = par_gen_tmp(i_type,k) % cmpr % index_i(ic)
        j = par_gen_tmp(i_type,k) % cmpr % index_j(ic)
        ij = (j-1)*nx_full + i
        if ( l_mass_source(ij) ) then
          i_layr_ij(ij) = i_layr_ij(ij) + 1
          l_mass_source(ij) = .false.
        end if
      end do

    end if  ! ( par_gen_tmp(i_type,k) % cmpr % n_points > 0 )

  end do  ! k = k_first, k_last, dk

  ! See if any columns have the highest number of
  ! mass-source layers layers yet found
  n_conv_layers = max( n_conv_layers, maxval(i_layr_ij)-1 )
  ! Note: number of layers in each column is actually
  ! i_layr_ij - 1
  ! because it gets incremented at the top of each layer found.

end do  ! i_type = 1, n_conv_types


! Only anything else to do if at least one convective source
! layer has been found
if ( n_conv_layers > 0 ) then

  ! Now that we know how many convecting layers there are,
  ! allocate the main par_gen array
  allocate( par_gen( n_conv_types, n_conv_layers,                              &
                     k_bot_conv:k_top_conv ) )

  ! Initialise number of mass-source points to zero for all
  ! levels, layers and types
  do k = k_bot_conv, k_top_conv
    do i_layr = 1, n_conv_layers
      do i_type = 1, n_conv_types
        par_gen(i_type,i_layr,k) % cmpr % n_points = 0
      end do
    end do
  end do

  ! NOTE: par_gen is allocated up to level k_top_conv even though
  ! it is only used up to level k_top_init (the levels between
  ! k_top_init and k_top_conv will just be left with n_points=0
  ! so that they don't do anything).  This is to avoid subscripting
  ! par_gen out-of-bounds in various places.

  ! 2nd loop over levels to actually divvy-up the initiation
  ! mass source data between different layers where needed
  do k = k_first, k_last, dk

    ! Loop over convection types
    do i_type = 1, n_conv_types

      ! If any mass sources from this level
      if ( par_gen_tmp(i_type,k) % cmpr % n_points > 0 ) then

        ! See if there are any compression lists that need to be
        ! split between multiple different distinct layers
        l_multi_layer = .false.
        over_n_points: do ic = 1, par_gen_tmp(i_type,k) % cmpr % n_points
          if ( .not. i_layr_cmpr(i_type,k)%i(ic)                               &
                  == i_layr_cmpr(i_type,k)%i(1) ) then
            l_multi_layer = .true.
            exit over_n_points
          end if
        end do over_n_points

        ! If compression list needs splitting up...
        if ( l_multi_layer ) then

          allocate( index_ic( par_gen_tmp(i_type,k) % cmpr                     &
                              % n_points ) )

          ! Loop over layers
          do i_layr = minval(i_layr_cmpr(i_type,k)%i),                         &
                      maxval(i_layr_cmpr(i_type,k)%i)

            ! Find indices of points which belong on this layer
            n_points = 0
            do ic = 1, par_gen_tmp(i_type,k) % cmpr % n_points
              if ( i_layr_cmpr(i_type,k)%i(ic) == i_layr ) then
                n_points = n_points + 1
                index_ic(n_points) = ic
              end if
            end do

            if ( n_points > 0 ) then

              ! Allocate arrays in output multi-layer structure
              par_gen(i_type,i_layr,k) % cmpr%n_points = n_points
              call parcel_alloc( l_tracer, n_points,                           &
                                 par_gen(i_type,i_layr,k) )

              ! Copy data into the right address in the output
              ! multi-layer structure
              call parcel_copy( l_tracer,                                      &
                                par_gen_tmp(i_type,k),                         &
                                par_gen(i_type,i_layr,k),                      &
                                index_ic=index_ic )

            end if

          end do  ! i_layr = minval(i_layr_cmpr(i_type,k)%i),
                  !          maxval(i_layr_cmpr(i_type,k)%i)

          deallocate( index_ic )

        else  ! ( l_multi_layer )

          ! All points belong on the same layer...

          i_layr = i_layr_cmpr(i_type,k)%i(1)
          n_points = par_gen_tmp(i_type,k) % cmpr % n_points

          ! Allocate arrays in the output multi-layer structure
          par_gen(i_type,i_layr,k) % cmpr % n_points = n_points
          call parcel_alloc( l_tracer, n_points,                               &
                             par_gen(i_type,i_layr,k) )

          ! Copy data into the right address in the output
          ! multi-layer structure
          call parcel_copy( l_tracer,                                          &
                            par_gen_tmp(i_type,k),                             &
                            par_gen(i_type,i_layr,k) )

        end if  ! ( l_multi_layer )

        ! Deallocate compressed layer indices as finished
        deallocate( i_layr_cmpr(i_type,k)%i )

      end if  ! ( par_gen_tmp(i_type,k) % cmpr % n_points > 0 )

    end do  ! i_type = 1, n_conv_types

  end do  ! k = k_first, k_last, dk

end if  ! ( n_conv_layers > 0 )


return
end subroutine par_gen_distinct_layers

end module par_gen_distinct_layers_mod
