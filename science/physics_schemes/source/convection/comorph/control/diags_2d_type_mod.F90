! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module diags_2d_type_mod

use diag_type_mod, only: diag_type, diag_list_type

implicit none

!----------------------------------------------------------------
! Type definition for a structure to store 2D diagnostics for
! updrafts or downdrafts
!----------------------------------------------------------------
! See the structure of addresses in fields_2d_mod
type :: diags_2d_type

  ! Number of diagnostics in this structure requested for output
  integer :: n_diags = 0

  ! List of pointers to active diagnostic structures
  ! (points to subset of the main list stored in the
  !  comorph_diags_type structure)
  type(diag_list_type), pointer :: list(:) => null()

  ! CAPE
  type(diag_type) :: cape

  ! Mass-flux-weighted CAPE
  type(diag_type) :: mfw_cape

end type diags_2d_type


contains


!----------------------------------------------------------------
! Subroutine to count the number of active diagnostics in a
! diags_2d_type structure and store a list of pointers to them
!----------------------------------------------------------------
! Note: this routine gets called twice
!   1st call (l_count_diags = .true.):
!     Check whether each diag is requested and count how many
!   2nd call (l_count_diags = .false.):
!     Set other properties for requested diags, and assign
!     pointers from active diags list.
subroutine diags_2d_assign( parent_name, l_count_diags, diags_2d,              &
                            parent_list, parent_i_diag )

use comorph_constants_mod, only: name_length
use diag_type_mod, only: diag_list_type, diag_assign, dom_type
use fields_2d_mod, only: i_cape, i_mfw_cape

implicit none

! Character string to prepend onto all the diagnostic names
! in here; e.g. "updraft" or "downdraft"
character(len=name_length), intent(in) :: parent_name

! Flag for 1st call to this routine, in which we just check
! whether each diag is requested and count up number requested
logical, intent(in) :: l_count_diags

! Structure storing meta-data for and pointers to the 2D diagnostics
type(diags_2d_type), intent(in out) :: diags_2d

! List of active diagnostics in the parent structure
type(diag_list_type), target, intent(in out) :: parent_list(:)

! Counter for active diagnostics in the parent structure
integer, intent(in out) :: parent_i_diag

! Structure containing the allowed doms (2d,3d,4d) for the diags
type(dom_type) :: doms

! Name for each diagnostic
character(len=name_length) :: diag_name

! Counters
integer :: i_diag, i_dummy


! Set fields_diags list to point at a section of the parent list
if ( l_count_diags ) then
  ! Not used in the first call, but needs to be assigned
  ! to avoid error when passing it into diag_assign
  diags_2d % list => parent_list
else
  ! Assign to subset address in parent:
  diags_2d % list => parent_list                                               &
    ( parent_i_diag + 1 :                                                      &
      parent_i_diag + diags_2d % n_diags )
end if

! Set allowed domain profiles for the diagnostics;
! Horizontal 2D domain, optionally separate for each convection type/layer
doms % x_y = .true.
doms % x_y_typ = .true.
doms % x_y_lay_typ = .true.

! Initialise local requested diagnostics counter to zero
i_diag = 0
! Initialise dummy for unused super-array address
i_dummy = 0

! Check for / process diagnostic requests...

! CAPE diagnostic
diag_name = trim(adjustl(parent_name)) // "_cape"
call diag_assign( diag_name, l_count_diags, doms, diags_2d%cape,               &
                  diags_2d%list, i_diag, i_dummy, i_field=i_cape )

! Mass-flux-weighted CAPE diagnostic
diag_name = trim(adjustl(parent_name)) // "_mfw_cape"
call diag_assign( diag_name, l_count_diags, doms, diags_2d%mfw_cape,           &
                  diags_2d%list, i_diag, i_dummy, i_field=i_mfw_cape )

! Increment parent active diagnostics counter
parent_i_diag = parent_i_diag + i_diag

! If this is the initial call to this routine to count the
! number of requested diags
if ( l_count_diags ) then

  ! Set number of active diags in this fields_diags structure
  diags_2d % n_diags = i_diag

  ! Nullify the list pointer ready to be reassigned next call
  diags_2d % list => null()

end if

return
end subroutine diags_2d_assign


!----------------------------------------------------------------
! Subroutine to compute means of 2D diagnostics over
! layers / types, and copy them into the output arrays
!----------------------------------------------------------------
subroutine diags_2d_compute_means( n_conv_types, n_conv_layers,                &
                                   ij_first, ij_last,                          &
                                   diags_2d, fields_2d )

use comorph_constants_mod, only: real_cvprec, real_hmprec, zero, min_float,    &
                                 nx_full, n_conv_layers_diag
use fields_2d_mod, only: n_fields_2d, i_m_ref, i_mfw_cape, i_cape

implicit none

! Number of convection types
integer, intent(in) :: n_conv_types
! Number of vertically distinct convecting layers
integer, intent(in) :: n_conv_layers

! i,j indices of the first and last point on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last

! Structure storing meta-data for and pointers to the 2D diagnostics
type(diags_2d_type), intent(in out) :: diags_2d

! Super-array storing 2D work arrays
real(kind=real_cvprec), intent(in) :: fields_2d                                &
                                      ( ij_first:ij_last, n_fields_2d,         &
                                        n_conv_types, n_conv_layers )

! Mean of CAPE over convecting layers, and over layers and types
real(kind=real_cvprec) :: cape_mean1 ( ij_first:ij_last, n_conv_types )
real(kind=real_cvprec) :: cape_mean2 ( ij_first:ij_last )
! Sums of reference mass-fluxes (used for computing mass-flux-weighted means)
real(kind=real_cvprec) :: m_ref_sum1 ( ij_first:ij_last, n_conv_types )
real(kind=real_cvprec) :: m_ref_sum2 ( ij_first:ij_last )

! Loop counters
integer :: i, j, ij, i_type, i_layr, i_diag, i_field


! Loop over requested diagnostics...
do i_diag = 1, diags_2d % n_diags

  ! Select stored address of this diagnostic in the fields_2d array
  i_field = diags_2d % list(i_diag)%pt % i_field

  if ( diags_2d % list(i_diag)%pt % request % x_y_lay_typ ) then
    ! If diag is requested separately for each layer and type...

    ! Copy to output array
    do i_layr = 1, min( n_conv_layers, n_conv_layers_diag )
      do i_type = 1, n_conv_types
        do ij = ij_first, ij_last
          ! Reverse ij = nx*(j-1)+i to get back i, j
          j = 1 + (ij-1) / nx_full
          i = ij - nx_full * (j-1)
          diags_2d % list(i_diag)%pt % field_4d(i,j,i_layr,i_type)             &
            = real( fields_2d(ij,i_field,i_type,i_layr), real_hmprec )
        end do
      end do
    end do

  end if

  if ( diags_2d % list(i_diag)%pt % request % x_y .or.                         &
       diags_2d % list(i_diag)%pt % request % x_y_typ ) then
    ! If diag is requested combined over multiple types or layers...

    if ( i_field == i_cape ) then
      ! If this is the CAPE diagnostic...

      do i_type = 1, n_conv_types
        ! For each convection type...

        ! Set combined CAPE to max over all layers
        do ij = ij_first, ij_last
          cape_mean1(ij,i_type) = zero
        end do
        do i_layr = 1, n_conv_layers
          do ij = ij_first, ij_last
            cape_mean1(ij,i_type) = max( cape_mean1(ij,i_type),                &
                                         fields_2d(ij,i_cape,i_type,i_layr) )
          end do
        end do

      end do  ! i_type = 1, n_conv_types

    else if ( i_field == i_mfw_cape ) then
      ! If this is the mass-flux-weighted CAPE diagnostic...

      do i_type = 1, n_conv_types
        ! For each convection type...

        ! Set combined CAPE to mass-flux-weighted mean over all layers
        do ij = ij_first, ij_last
          cape_mean1(ij,i_type) = zero
          m_ref_sum1(ij,i_type) = zero
        end do
        do i_layr = 1, n_conv_layers
          do ij = ij_first, ij_last
            cape_mean1(ij,i_type) = cape_mean1(ij,i_type)                      &
                                  + fields_2d(ij,i_mfw_cape,i_type,i_layr)     &
                                  * fields_2d(ij,i_m_ref,i_type,i_layr)
            m_ref_sum1(ij,i_type) = m_ref_sum1(ij,i_type)                      &
                                  + fields_2d(ij,i_m_ref,i_type,i_layr)
          end do
        end do
        do ij = ij_first, ij_last
          cape_mean1(ij,i_type) = cape_mean1(ij,i_type)                        &
                                / max( m_ref_sum1(ij,i_type), min_float )
        end do

      end do  ! i_type = 1, n_conv_types

    end if  ! ( i_field )

    if ( diags_2d % list(i_diag)%pt % request % x_y_typ ) then
      ! If diag is requested separately for each type...

      ! Copy to output array
      do i_type = 1, n_conv_types
        do ij = ij_first, ij_last
          ! Reverse ij = nx*(j-1)+i to get back i, j
          j = 1 + (ij-1) / nx_full
          i = ij - nx_full * (j-1)
          diags_2d % list(i_diag)%pt % field_3d(i,j,i_type)                    &
            = real( cape_mean1(ij,i_type), real_hmprec )
        end do
      end do

    end if

    if ( diags_2d % list(i_diag)%pt % request % x_y ) then
      ! If diag is requested combined over all types...

      if ( i_field == i_cape ) then
        ! If this is the CAPE diagnostic...

        ! Set combined CAPE to max over all types
        do ij = ij_first, ij_last
          cape_mean2(ij) = zero
        end do
        do i_type = 1, n_conv_types
          do ij = ij_first, ij_last
            cape_mean2(ij) = max( cape_mean2(ij),                              &
                                  cape_mean1(ij,i_type) )
          end do
        end do

      else if ( i_field == i_mfw_cape ) then
        ! If this is the mass-flux-weighted CAPE diagnostic...

        ! Set combined CAPE to mass-flux-weighted mean over all types
        do ij = ij_first, ij_last
          cape_mean2(ij) = zero
          m_ref_sum2(ij) = zero
        end do
        do i_type = 1, n_conv_types
          do ij = ij_first, ij_last
            cape_mean2(ij) = cape_mean2(ij)                                    &
                           + cape_mean1(ij,i_type)                             &
                           * m_ref_sum1(ij,i_type)
            m_ref_sum2(ij) = m_ref_sum2(ij)                                    &
                           + m_ref_sum1(ij,i_type)
          end do
        end do
        do ij = ij_first, ij_last
          cape_mean2(ij) = cape_mean2(ij)                                      &
                         / max( m_ref_sum2(ij), min_float )
        end do

      end if  ! ( i_field )

      ! Copy to output array
      do ij = ij_first, ij_last
        ! Reverse ij = nx*(j-1)+i to get back i, j
        j = 1 + (ij-1) / nx_full
        i = ij - nx_full * (j-1)
        diags_2d % list(i_diag)%pt % field_2d(i,j)                             &
          = real( cape_mean2(ij), real_hmprec )
      end do

    end if  ! ( diags_2d % list(i_diag)%pt % request % x_y )

  end if  ! ( diags_2d % list(i_diag)%pt % request % x_y .or.
          !   diags_2d % list(i_diag)%pt % request % x_y_typ )

end do  ! i_diag = 1, diags_2d % n_diags

return
end subroutine diags_2d_compute_means


end module diags_2d_type_mod
