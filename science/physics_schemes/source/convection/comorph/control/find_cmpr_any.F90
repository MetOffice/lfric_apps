! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module find_cmpr_any_mod

implicit none

contains

! Subroutine to find the combined compression list of points
! where any type of convection is active, on each level
subroutine find_cmpr_any( n_updraft_layers, n_dndraft_layers,                  &
                          ij_first, ij_last,                                   &
                          updraft_res_source,                                  &
                          updraft_fallback_res_source,                         &
                          dndraft_res_source,                                  &
                          dndraft_fallback_res_source,                         &
                          cmpr_any )

use comorph_constants_mod, only: nx_full, k_bot_conv, k_top_conv,              &
                                 n_updraft_types, n_dndraft_types,             &
                                 l_updraft_fallback, l_dndraft_fallback
use cmpr_type_mod, only: cmpr_type, cmpr_alloc
use res_source_type_mod, only: res_source_type

implicit none

! Number of mass-source layers found for updrafts and downdrafts
integer, intent(in) :: n_updraft_layers
integer, intent(in) :: n_dndraft_layers

! Position indices ( nx*(j-1) + i ) of the first and last
! convecting point on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

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

! Structure storing compression list indices for the points where
! any kind of convection is occuring
type(cmpr_type), intent(in out) :: cmpr_any( k_bot_conv:k_top_conv )
! intent(inout) so that non-allocated status of the contained
! compression index arrays is preserved on input.

! Temporary store for ij indices of the combined
! compression list points
integer :: index_ij(ij_last-ij_first+1)

! Work array set to true at points where any kind of convection
! is occuring
logical :: l_mask(ij_first:ij_last)

! Flag for whether any convection on current level
logical :: l_any

! Loop counters
integer :: i, j, k, ic, nc, ij, i_type, i_layr


! Initialise mask to false
do ij = ij_first, ij_last
  l_mask(ij) = .false.
end do


! Loop over levels
do k = k_bot_conv, k_top_conv


  ! Initialise flag to false.  It will beset to true if any
  ! resolved-scale source terms on this level
  l_any = .false.


  ! Set mask to true at all points where any kind of
  ! convection is occuring

  if ( n_updraft_types > 0 .and. n_updraft_layers > 0 ) then

    do i_layr = 1, n_updraft_layers
      do i_type = 1, n_updraft_types
        nc = updraft_res_source(i_type,i_layr,k)                               &
             % cmpr % n_points
        if ( nc > 0 ) then
          l_any = .true.
          do ic = 1, nc
            i = updraft_res_source(i_type,i_layr,k)                            &
                % cmpr % index_i(ic)
            j = updraft_res_source(i_type,i_layr,k)                            &
                % cmpr % index_j(ic)
            l_mask( nx_full*(j-1)+i ) = .true.
          end do
        end if
      end do
    end do

    if ( l_updraft_fallback ) then
      do i_layr = 1, n_updraft_layers
        do i_type = 1, n_updraft_types
          nc = updraft_fallback_res_source(i_type,i_layr,k)                    &
               % cmpr % n_points
          if ( nc > 0 ) then
            l_any = .true.
            do ic = 1, nc
              i = updraft_fallback_res_source(i_type,i_layr,k)                 &
                  % cmpr % index_i(ic)
              j = updraft_fallback_res_source(i_type,i_layr,k)                 &
                  % cmpr % index_j(ic)
              l_mask( nx_full*(j-1)+i ) = .true.
            end do
          end if
        end do
      end do
    end if

  end if  ! ( n_updraft_types > 0 .and. n_updraft_layers > 0 )

  if ( n_dndraft_types > 0 .and. n_dndraft_layers > 0 ) then

    do i_layr = 1, n_dndraft_layers
      do i_type = 1, n_dndraft_types
        nc = dndraft_res_source(i_type,i_layr,k)                               &
             % cmpr % n_points
        if ( nc > 0 ) then
          l_any = .true.
          do ic = 1, nc
            i = dndraft_res_source(i_type,i_layr,k)                            &
                % cmpr % index_i(ic)
            j = dndraft_res_source(i_type,i_layr,k)                            &
                % cmpr % index_j(ic)
            l_mask( nx_full*(j-1)+i ) = .true.
          end do
        end if
      end do
    end do

    if ( l_dndraft_fallback ) then
      do i_layr = 1, n_dndraft_layers
        do i_type = 1, n_dndraft_types
          nc = dndraft_fallback_res_source(i_type,i_layr,k)                    &
               % cmpr % n_points
          if ( nc > 0 ) then
            l_any = .true.
            do ic = 1, nc
              i = dndraft_fallback_res_source(i_type,i_layr,k)                 &
                  % cmpr % index_i(ic)
              j = dndraft_fallback_res_source(i_type,i_layr,k)                 &
                  % cmpr % index_j(ic)
              l_mask( nx_full*(j-1)+i ) = .true.
            end do
          end if
        end do
      end do
    end if

  end if  ! ( n_dndraft_types > 0 .and. n_dndraft_layers > 0 )


  ! Initialise count for number of convecting points to zero
  cmpr_any(k) % n_points = 0

  ! If any convection occuring on this level
  if ( l_any ) then

    ! Loop over all points on current segment
    do ij = ij_first, ij_last
      ! If mask is true, add point to list
      if ( l_mask(ij) ) then
        cmpr_any(k) % n_points = cmpr_any(k) % n_points + 1
        index_ij(cmpr_any(k) % n_points) = ij
        ! Reset mask to false ready for next level
        l_mask(ij) = .false.
      end if
    end do

    ! Set the cmpr_any fields using the above
    call cmpr_alloc( cmpr_any(k), cmpr_any(k) % n_points )
    do ic = 1, cmpr_any(k) % n_points
      ! Reverse ij = nx*(j-1)+i to get back i, j
      j = 1 + (index_ij(ic)-1) / nx_full
      i = index_ij(ic) - nx_full * (j-1)
      cmpr_any(k) % index_i(ic) = i
      cmpr_any(k) % index_j(ic) = j
    end do

  end if  ! ( l_any )


end do  ! k = k_bot_conv, k_top_conv


return
end subroutine find_cmpr_any


end module find_cmpr_any_mod
