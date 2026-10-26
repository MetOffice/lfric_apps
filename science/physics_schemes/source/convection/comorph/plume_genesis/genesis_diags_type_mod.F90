! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module genesis_diags_type_mod

use diag_type_mod, only: diag_type, diag_list_type
use subregion_diags_type_mod, only: subregion_diags_type
use subregion_mod, only: n_regions

implicit none


! Type definition to store the diagnostics available from comorph's
! initiation mass-source calculations
type :: genesis_diags_type

  ! Number of diagnostics in this structure requested for output
  integer :: n_diags = 0

  ! List of pointers to active diagnostic structures
  ! (points to subset of the main list in the
  !  comorph_diags_type structure)
  type(diag_list_type), pointer :: list(:) => null()

  ! Properties of each of the grid-box sub-regions;
  ! dry (no cloud or precep)
  ! liq (liquid-only cloud)
  ! mph (mixed-phased cloud)
  ! icr (ice and/or rain but no liquid cloud)
  type(subregion_diags_type) :: subregion_diags(n_regions)

end type genesis_diags_type


contains


!----------------------------------------------------------------
! Subroutine to count the number of active diagnostics in the
! genesis_diags structure and store a list of pointers to them
!----------------------------------------------------------------
! Note: this routine gets called twice
!   1st call (l_count_diags = .TRUE.):
!     Check whether each diag is requested and count how many
!   2nd call (l_count_diags = .FALSE.):
!     Set other properties for requested diags, and assign
!     pointers from active diags list.
subroutine genesis_diags_assign( l_count_diags, genesis_diags,                 &
                                 parent_list, parent_i_diag )

use diag_type_mod, only: diag_list_type
use subregion_mod, only: n_regions, region_names
use subregion_diags_type_mod, only: subregion_diags_assign

implicit none

! Flag for 1st call to this routine, in which we just check
! whether each diag is requested and count up number requested
logical, intent(in) :: l_count_diags

! The main diagnostics structure, containing pointers
! for all diagnostics, including those not in use
type(genesis_diags_type), target, intent(in out) :: genesis_diags

! List of active diagnostics in the parent structure
type(diag_list_type), target, intent(in out) :: parent_list(:)

! Counter for active diagnostics in the parent structure
integer, intent(in out) :: parent_i_diag

! Counters
integer :: i_diag, i_super, i_region


! Set genesis_diags list to point at a section of parent list
if ( l_count_diags ) then
  ! Not used in the init call, but needs to be assigned
  ! to avoid error when passing it into diag_3d_assign
  genesis_diags % list => parent_list
else
  ! Assign to subset address in parent:
  genesis_diags % list => parent_list                                          &
    ( parent_i_diag + 1 :                                                      &
      parent_i_diag + genesis_diags % n_diags )
end if

! Initialise active 3-D diagnostic counter to zero
i_diag = 0

! Initialise super-array address for diagnostics
i_super = 0

! Fractions and fields for sub-regions of the grid-box
do i_region = 1, n_regions
  call subregion_diags_assign( region_names(i_region), l_count_diags,          &
         genesis_diags % subregion_diags(i_region),                            &
         genesis_diags % list, i_diag, i_super )
end do

! Increment parent active diagnostics counter
parent_i_diag = parent_i_diag + i_diag


! If this is the initial call to this routine to count the
! number of active diags
if ( l_count_diags ) then

  ! Set number of active diags in this subregion_diags structure
  genesis_diags % n_diags = i_diag

  ! Nullify the list pointer ready to be reassigned next call
  genesis_diags % list => null()

end if


return
end subroutine genesis_diags_assign


end module genesis_diags_type_mod
