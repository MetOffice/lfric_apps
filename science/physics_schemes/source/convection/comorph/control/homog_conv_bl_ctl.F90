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
                              ij_first, ij_last,                               &
                              n_fields_tot, l_down,                            &
                              grid, fields_np1, layer_mass,                    &
                              par_bl_top, turb, res_source )

use comorph_constants_mod, only: real_hmprec,                                  &
                     nx_full, ny_full, k_bot_conv, k_top_conv
use cmpr_type_mod, only: cmpr_dealloc
use grid_type_mod, only: grid_type
use fields_type_mod, only: fields_type
use turb_type_mod, only: turb_type
use parcel_type_mod, only: parcel_type
use res_source_type_mod, only: res_source_type, n_res
use cloudfracs_type_mod, only: n_convcloud
use homog_conv_bl_mod, only: homog_conv_bl

implicit none

! Number of convection types
integer, intent(in) :: n_conv_types

! Number of distinct convecting layers
integer, intent(in) :: n_conv_layers

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
! intent inout as we deallocate the contained arrays in here.

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

! Grid of flags for points where at least some convection occured
! above the boundary-layer top
logical :: l_above_bl( ij_first:ij_last )

! Work store for number of points
integer :: nc

! Flag for whether any convection of each layer / type
! crossed the BL-top
logical :: l_any( n_conv_types, n_conv_layers )

! Number of points where convection not entirely below BL-top
integer :: n_above
! Indices of those points in the res_source arrays
integer, allocatable :: index_ic_above(:)

! First index found where convection entirely below BL-top
integer :: ic_first

! Loop counters
integer :: i, j, ij, k, ic, ic2, i_type, i_layr, i_field


! Find highest model-level where convection passed
! the boundary-layer top
k_max = k_bot_conv - 1
do i_layr = 1, n_conv_layers
  do i_type = 1, n_conv_types
    l_any(i_type,i_layr) = .false.
    do k = k_bot_conv, k_top_conv
      nc = par_bl_top(i_type,i_layr,k) % cmpr % n_points
      if ( nc > 0 ) then
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
          nc = res_source(i_type,i_layr,k) % cmpr % n_points
          if ( nc > 0 ) then
            do ic = 1, nc
              i = res_source(i_type,i_layr,k) % cmpr %index_i(ic)
              j = res_source(i_type,i_layr,k) % cmpr %index_j(ic)
              index_ic( nx_full*(j-1)+i, k ) = ic
            end do
          end if
        end do

        ! Loop over possible levels where parcel crossed BL-top
        do k = k_bot_conv, k_max
          nc = par_bl_top(i_type,i_layr,k) % cmpr % n_points
          if ( nc > 0 ) then

            ! Call routine to vertically homogenise the
            ! resolved-scale source terms in columns which
            ! hit the BL-top at the current model-level
            call homog_conv_bl( n_conv_types, n_conv_layers,                   &
                                n_fields_tot, l_down,                          &
                                i_type, i_layr, k,                             &
                                ij_first, ij_last, index_ic,                   &
                                grid, fields_np1, layer_mass,                  &
                                par_bl_top(i_type,i_layr,k),                   &
                                res_source )

            ! Deallocate stored parcel properties at BL-top
            call cmpr_dealloc( par_bl_top(i_type,i_layr,k)%cmpr )
            deallocate( par_bl_top(i_type,i_layr,k) % par_super )
            deallocate( par_bl_top(i_type,i_layr,k) % mean_super)

          end if
        end do

        ! Reset indices to zero ready for next type / layer
        do k = k_bot_conv, k_max
          nc = res_source(i_type,i_layr,k) % cmpr % n_points
          if ( nc > 0 ) then
            do ic = 1, nc
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

! Find max number of points undergoing convection of a single
! type / layer / level
n_above = 0
do i_layr = 1, n_conv_layers
  do i_type = 1, n_conv_types
    l_any(i_type,i_layr) = .false.
    do k = k_bot_conv, k_top_conv
      nc = res_source(i_type,i_layr,k) % cmpr % n_points
      if ( nc > 0 ) then
        ! Save max number of points
        n_above = max( n_above, nc )
        ! Set flag for convection of this type / layer
        l_any(i_type,i_layr) = .true.
      end if
    end do
  end do
end do


! Allocate array to store points where convection is
! entirely below the BL-top
allocate( index_ic_above(n_above) )


! Loop over layrs and types
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
        nc = res_source(i_type,i_layr,k) % cmpr % n_points
        if ( nc > 0 ) then
          do ic = 1, nc
            i = res_source(i_type,i_layr,k) % cmpr % index_i(ic)
            j = res_source(i_type,i_layr,k) % cmpr % index_j(ic)
            ij = nx_full*(j-1)+i
            if ( .not. l_above_bl(ij) ) then
              ! Model-level is considered above the BL-top if its
              ! lower interface height_half(k) is above z_bl_top:
              if ( grid % height_half(i,j,k)                                   &
                 > turb % z_bl_top(i,j) ) l_above_bl(ij) = .true.
            end if
          end do
        end if
      end do

      ! Loop over levels
      do k = k_bot_conv, k_top_conv
        nc = res_source(i_type,i_layr,k) % cmpr % n_points
        if ( nc > 0 ) then

          ! See if any of the res_sources are in columns which
          ! are entirely below the BL-top
          ic_first = 0
          over_nc: do ic = 1, nc
            i = res_source(i_type,i_layr,k) % cmpr % index_i(ic)
            j = res_source(i_type,i_layr,k) % cmpr % index_j(ic)
            ij = nx_full*(j-1)+i
            if ( .not. l_above_bl(ij) ) then
              ! Store index and exit the loop when we find the
              ! first point that needs to be excluded.
              ic_first = ic
              exit over_nc
            end if
          end do over_nc

          ! If at least one point not above, we have work to do
          if ( ic_first > 0 ) then

            ! All points before ic_first do go above BL-top
            n_above = ic_first - 1

            if ( ic_first < nc ) then

              ! See if there are any more...
              do ic = ic_first+1, nc
                i = res_source(i_type,i_layr,k)%cmpr%index_i(ic)
                j = res_source(i_type,i_layr,k)%cmpr%index_j(ic)
                ij = nx_full*(j-1)+i
                if ( l_above_bl(ij) ) then
                  n_above = n_above + 1
                  index_ic_above(n_above) = ic
                end if
              end do

            end if

            ! Set new number of points in the res_source
            ! compression list
            res_source(i_type,i_layr,k) % cmpr % n_points                      &
              = n_above

            ! If any more above points found after below
            ! points, copy them backwards to overwrite the
            ! below points
            if ( n_above >= ic_first ) then

              ! Compress the arrays in-situ, leaving junk in the
              ! elements with ic > n_above...

              ! Note that ic_first (saved earlier) is the index
              ! of the first terminated point in the existing
              ! list, and therefore corresponds to the index of
              ! the first point that needs to be moved to form
              ! the recompressed list.

              ! Indices of points in the full 2-D grid
              do ic2 = ic_first, n_above
                ic = index_ic_above(ic2)
                res_source(i_type,i_layr,k) % cmpr%index_i(ic2)                &
                 = res_source(i_type,i_layr,k) % cmpr%index_i(ic)
                res_source(i_type,i_layr,k) % cmpr%index_j(ic2)                &
                 = res_source(i_type,i_layr,k) % cmpr%index_j(ic)
              end do

              ! Dry-mass increments from ent/det-rainment,
              ! and convective cloud fields
              do i_field = 1, n_res
                do ic2 = ic_first, n_above
                  ic = index_ic_above(ic2)
                  res_source(i_type,i_layr,k)                                  &
                        % res_super(ic2,i_field)                               &
                      = res_source(i_type,i_layr,k)                            &
                        % res_super(ic,i_field)
                end do
              end do

              ! Primary field source terms
              do i_field = 1, n_fields_tot
                do ic2 = ic_first, n_above
                  ic = index_ic_above(ic2)
                  res_source(i_type,i_layr,k)                                  &
                        % fields_super(ic2,i_field)                            &
                      = res_source(i_type,i_layr,k)                            &
                        % fields_super(ic,i_field)
                end do
              end do

              ! Convective cloud fields
              if ( n_convcloud > 0 ) then
                do i_field = 1, n_convcloud
                  do ic2 = ic_first, n_above
                    ic = index_ic_above(ic2)
                    res_source(i_type,i_layr,k)                                &
                        % convcloud_super(ic2,i_field)                         &
                      = res_source(i_type,i_layr,k)                            &
                        % convcloud_super(ic,i_field)
                  end do
                end do
              end if

            end if  ! ( n_above >= ic_first )

          end if  ! ( ic_first > 0 )

        end if  ! ( nc > 0 )
      end do  ! k = k_bot_conv, k_top_conv

    end if  ! ( l_any(i_type,i_layr) )
  end do  ! i_type = 1, n_conv_types
end do  ! i_layr = 1, n_conv_layers


! Deallocate work array
deallocate( index_ic_above )


return
end subroutine homog_conv_bl_ctl


end module homog_conv_bl_ctl_mod
