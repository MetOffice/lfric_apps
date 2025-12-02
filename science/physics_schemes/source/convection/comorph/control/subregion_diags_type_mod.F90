! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module subregion_diags_type_mod

use diag_type_mod, only: diag_type, diag_list_type
use fields_diags_type_mod, only: fields_diags_type

implicit none


! Type definition to store the diags available for
! each sub-region
type :: subregion_diags_type

  ! Number of diagnostics in this structure requested for output
  integer :: n_diags = 0

  ! List of pointers to active diagnostic structures
  ! (points to subset of the main list in the
  !  comorph_diags_type structure)
  type(diag_list_type), pointer :: list(:) => null()

  ! Fractional area
  type(diag_type) :: frac

  ! Primary fields
  type(fields_diags_type) :: fields

end type subregion_diags_type


contains


!----------------------------------------------------------------
! Subroutine to count the number of active diagnostics in a
! subregion_diags structure and store a list of pointers to them
!----------------------------------------------------------------
! Note: this routine gets called twice
!   1st call (l_count_diags = .true.):
!     Check whether each diag is requested and count how many
!   2nd call (l_count_diags = .false.):
!     Set other properties for requested diags, and assign
!     pointers from active diags list.
subroutine subregion_diags_assign( region_name, l_count_diags,                 &
                                   subregion_diags,                            &
                                   parent_list, parent_i_diag )

use comorph_constants_mod, only: name_length
use diag_type_mod, only: diag_list_type, diag_assign, dom_type
use fields_diags_type_mod, only: fields_diags_assign

implicit none

! Name of current subgrid region
character(len=name_length), intent(in) :: region_name

! Flag for 1st call to this routine, in which we just check
! whether each diag is requested and count up number requested
logical, intent(in) :: l_count_diags

! The main diagnostics structure, containing pointers
! for all diagnostics, including those not in use
type(subregion_diags_type), target, intent(in out) ::                          &
                                               subregion_diags

! List of active diagnostics in the parent structure
type(diag_list_type), target, intent(in out) :: parent_list(:)

! Counter for active diagnostics in the parent structure
integer, intent(in out) :: parent_i_diag

! Name for each diagnostic
character(len=name_length) :: diag_name

! Structure containing the allowed domain profiles for the diags
type(dom_type) :: doms

! Counters
integer :: i_diag, i_super


! Set subregion_diags list to point at a section of parent list
if ( l_count_diags ) then
  ! Not used in the init call, but needs to be assigned
  ! to avoid error when passing it into diag_3d_assign
  subregion_diags % list => parent_list
else
  ! Assign to subset address in parent:
  subregion_diags % list => parent_list                                        &
    ( parent_i_diag + 1 :                                                      &
      parent_i_diag + subregion_diags % n_diags )
end if

! Initialise active 3-D diagnostic counter to zero
i_diag = 0

! Initialise super-array address for diagnostics
i_super = 0

! Allowed domain profile for the following diags is 3-D only.
doms % x_y_z = .true.


! Fractional area
diag_name = trim(adjustl(region_name)) // "_frac"
call diag_assign( diag_name, l_count_diags, doms,                              &
                  subregion_diags % frac,                                      &
                  subregion_diags % list, i_diag, i_super )

! Primary fields
call fields_diags_assign( region_name, l_count_diags, doms,                    &
                          subregion_diags % fields,                            &
                          subregion_diags % list,                              &
                          i_diag, i_super )


! Increment parent active diagnostics counter
parent_i_diag = parent_i_diag + i_diag


! If this is the initial call to this routine to count the
! number of active diags
if ( l_count_diags ) then

  ! Set number of active diags in this subregion_diags structure
  subregion_diags % n_diags = i_diag

  ! Nullify the list pointer ready to be reassigned next call
  subregion_diags % list => null()

end if


return
end subroutine subregion_diags_assign


end module subregion_diags_type_mod
