! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module homog_conv_bl_ctl_mod

implicit none

contains


! Subroutine to do higher-level data-management and call
! the routine homog_conv_bl to homogenize the convective
! resolved-scale source terms within the boundary-layer
subroutine homog_conv_bl_ctl( n_conv_types, n_conv_layers,                     &
                              max_points, ij_first, ij_last,                   &
                              n_fields_tot, l_down,                            &
                              grid, fields_np1, layer_mass,                    &
                              par_bl_top, turb, res_source )

use comorph_constants_mod, only: real_hmprec,                                  &
                                 nx_full, ny_full, k_bot_conv, k_top_conv
use cmpr_type_mod, only: cmpr_dealloc
use grid_type_mod, only: grid_type
use fields_type_mod, only: fields_type
use turb_type_mod, only: turb_type
use parcel_type_mod, only: parcel_type, i_massflux_d
use res_source_type_mod, only: res_source_type, res_source_compress
use homog_conv_bl_mod, only: homog_conv_bl
use set_l_within_bl_mod, only: set_l_within_bl

implicit none

! Number of convection types
integer, intent(in) :: n_conv_types

! Number of distinct convecting layers
integer, intent(in) :: n_conv_layers

! Total number of points on the current segment
integer, intent(in) :: max_points

! ij = nx_full*(j-1)+i, for the first and last point where
! convection is occurring
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last


! Number of transported fields (including tracers)
integer, intent(in) :: n_fields_tot

! Flag for downdraft versus updraft
logical, intent(in) :: l_down

! Structure containing model grid-fields
type(grid_type), intent(in) :: grid

! Structure containing full 3-D primary fields
type(fields_type), intent(in) :: fields_np1

! Full 3-D array of dry-mass per unit surface area
! contained in each grid-cell
real(kind=real_hmprec), intent(in) :: layer_mass                               &
       ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Array of structures containing parcel properties at the
! first level above the boundary-layer top, for each convection
! type and layer
type(parcel_type), intent(in out) :: par_bl_top                                &
          ( n_conv_types, n_conv_layers, k_bot_conv:k_top_conv )
! intent inout as homog_conv_bl modifies the contained fields to
! interpolate the parcel properties to the accurate BL-height.

! Structure containing turbulence fields, including BL-top height
type(turb_type), intent(in) :: turb

! Resolved-scale source terms to be modified
type(res_source_type), intent(in out) :: res_source                            &
          ( n_conv_types, n_conv_layers, k_bot_conv:k_top_conv )

! Highest model-level where convection crossed the BL-top
integer :: k_max

! Array to store compression list indices of the resolved-scale
! source term arrays on a grid, to enable addressing res_source
! arrays from the par_bl_top compression lists
integer, allocatable :: index_ic(:,:)

! Flag for points that are within the boundary-layer
logical :: l_within_bl(max_points)

! Grid of flags for points where at least some convection occured
! above the boundary-layer top
logical :: l_above_bl( ij_first:ij_last )

! Flag for whether any convection of each layer / type
! crossed the BL-top
logical :: l_any( n_conv_types, n_conv_layers )

! Number of points where convection not entirely below BL-top
integer :: n_above
! Indices of those points in the res_source arrays
integer :: index_ic_above(max_points)

! Array bounds for BL-height and grid heights
integer :: lbz(2), ubz(2)
integer :: lbh(3), ubh(3)

! Model-level up-to-which to homogenise source-terms
integer :: k_bl_top

! Loop counters
integer :: i, j, ij, k, ic, i_type, i_layr


! Find array bounds for the BL-top height array and half-level heights
lbz = lbound( turb % z_bl_top )
ubz = ubound( turb % z_bl_top )
lbh = lbound( grid % height_half )
ubh = ubound( grid % height_half )

! Find highest model-level where convection passed
! the boundary-layer top
k_max = k_bot_conv - 1
do i_layr = 1, n_conv_layers
  do i_type = 1, n_conv_types
    l_any(i_type,i_layr) = .false.
    do k = k_bot_conv, k_top_conv
      if ( par_bl_top(i_type,i_layr,k) % cmpr % n_points > 0 ) then
        ! Update highest model-level
        k_max = max( k_max, k )
        ! Set flag for convection of this type / layer
        l_any(i_type,i_layr) = .true.
      end if
    end do
  end do
end do


! If any convection crossed the BL-top
if ( k_max >= k_bot_conv ) then

  ! Allocate array to store indices of res_source points
  allocate( index_ic( ij_first:ij_last, k_bot_conv:k_max ) )
  do k = k_bot_conv, k_max
    do ij = ij_first, ij_last
      index_ic(ij,k) = 0
    end do
  end do

  ! Loop over layers and types
  do i_layr = 1, n_conv_layers
    do i_type = 1, n_conv_types
      ! If any convection of this type / layer crossed the BL-top
      if ( l_any(i_type,i_layr) ) then

        ! Scatter res_source compression indices for the current
        ! layer / type into the index array
        do k = k_bot_conv, k_max
          if ( res_source(i_type,i_layr,k) % cmpr % n_points > 0 ) then
            do ic = 1, res_source(i_type,i_layr,k) % cmpr % n_points
              i = res_source(i_type,i_layr,k) % cmpr %index_i(ic)
              j = res_source(i_type,i_layr,k) % cmpr %index_j(ic)
              index_ic( nx_full*(j-1)+i, k ) = ic
            end do
          end if
        end do

        ! Loop over possible levels where parcel crossed BL-top
        do k = k_bot_conv, k_max
          if ( par_bl_top(i_type,i_layr,k) % cmpr % n_points > 0 ) then

            ! TEMPORARY CODE TO PRESERVE KGO; TO BE REMOVED SOON
            ! (only homogenize up to and including k-1 for downdrafts)
            if ( l_down ) then
              k_bl_top = k - 1
            else
              k_bl_top = k
            end if

            ! Call routine to vertically homogenise the
            ! resolved-scale source terms in columns which
            ! hit the BL-top at the current model-level
            call homog_conv_bl( par_bl_top(i_type,i_layr,k) % cmpr % n_points, &
                                n_conv_types, n_conv_layers,                   &
                                n_fields_tot, l_down,                          &
                                i_type, i_layr, k_bl_top,                      &
                                ij_first, ij_last, index_ic,                   &
                                grid, fields_np1, layer_mass,                  &
                                par_bl_top(i_type,i_layr,k) % cmpr,            &
                                par_bl_top(i_type,i_layr,k) % par_super        &
                                                              (:,i_massflux_d),&
                                par_bl_top(i_type,i_layr,k) % mean_super,      &
                                res_source )

          end if
        end do

        ! Reset indices to zero ready for next type / layer
        do k = k_bot_conv, k_max
          if ( res_source(i_type,i_layr,k) % cmpr % n_points > 0 ) then
            do ic = 1, res_source(i_type,i_layr,k) % cmpr % n_points
              i = res_source(i_type,i_layr,k) % cmpr %index_i(ic)
              j = res_source(i_type,i_layr,k) % cmpr %index_j(ic)
              index_ic( nx_full*(j-1)+i, k ) = 0
            end do
          end if
        end do

      end if  ! ( l_any(i_type,i_layr) )
    end do  ! i_type = 1, n_conv_types
  end do  ! i_layr = 1, n_conv_layers

  ! Deallocate index array
  deallocate( index_ic )

end if  ! ( k_max >= k_bot_conv )


! Any convective columns that remained entirely inside the
! boundary-layer and never passed through the BL-top ought to
! be homogenized entirely, though they will not be included in
! the compression lists of points where convection crossed
! the BL-top.
! Total homogenisation of these columns amounts to zeroing the
! resolved-scale source terms.

! Set flag for any convection of each type / layer
do i_layr = 1, n_conv_layers
  do i_type = 1, n_conv_types
    l_any(i_type,i_layr) = .false.
    do k = k_bot_conv, k_top_conv
      if ( res_source(i_type,i_layr,k) % cmpr % n_points > 0 ) then
        l_any(i_type,i_layr) = .true.
      end if
    end do
  end do
end do

! Loop over active layers and types
do i_layr = 1, n_conv_layers
  do i_type = 1, n_conv_types
    if ( l_any(i_type,i_layr) ) then

      ! Initialise flags to false
      do ij = ij_first, ij_last
        l_above_bl(ij) = .false.
      end do

      ! Change flag to true at points where convection reached
      ! above the boundary-layer top
      do k = k_bot_conv, k_top_conv
        if ( res_source(i_type,i_layr,k) % cmpr % n_points > 0 ) then
          ! Call routine to flag points which are within the BL
          call set_l_within_bl( res_source(i_type,i_layr,k) % cmpr, k,         &
                                lbh, ubh, grid % height_half,                  &
                                lbz, ubz, turb % z_bl_top,                     &
                                l_within_bl )
          ! Set l_above_bl to true at points where flag is false
          do ic = 1, res_source(i_type,i_layr,k) % cmpr % n_points
            if ( .not. l_within_bl(ic) ) then
              i = res_source(i_type,i_layr,k) % cmpr % index_i(ic)
              j = res_source(i_type,i_layr,k) % cmpr % index_j(ic)
              ij = nx_full*(j-1)+i
              l_above_bl(ij) = .true.
            end if
          end do
        end if
      end do

      ! Loop over levels
      do k = k_bot_conv, k_top_conv
        if ( res_source(i_type,i_layr,k) % cmpr % n_points > 0 ) then

          ! Count points which are in columns where convection extended
          ! above the BL-top
          n_above = 0
          do ic = 1, res_source(i_type,i_layr,k) % cmpr % n_points
            i = res_source(i_type,i_layr,k) % cmpr % index_i(ic)
            j = res_source(i_type,i_layr,k) % cmpr % index_j(ic)
            ij = nx_full*(j-1)+i
            if ( l_above_bl(ij) ) then
              n_above = n_above + 1
              index_ic_above(n_above) = ic
            end if
          end do

          ! If at least one point not above, we have work to do
          if ( n_above < res_source(i_type,i_layr,k) % cmpr % n_points ) then

            ! Set new number of points in the res_source
            ! compression list
            res_source(i_type,i_layr,k) % cmpr % n_points = n_above

            ! If any above points found, copy them backwards to overwrite
            ! the below points
            if ( n_above > 0 ) then
              call res_source_compress( n_fields_tot, n_above, index_ic_above, &
                                        res_source(i_type,i_layr,k) )
            end if

          end if  ! ( n_above < res_source(i_type,i_layr,k) % cmpr % n_points )

        end if  ! ( res_source(i_type,i_layr,k) % cmpr % n_points > 0 )
      end do  ! k = k_bot_conv, k_top_conv

    end if  ! ( l_any(i_type,i_layr) )
  end do  ! i_type = 1, n_conv_types
end do  ! i_layr = 1, n_conv_layers


return
end subroutine homog_conv_bl_ctl


end module homog_conv_bl_ctl_mod
