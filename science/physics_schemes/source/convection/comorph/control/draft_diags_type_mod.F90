! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module draft_diags_type_mod

use comorph_constants_mod, only: real_cvprec
use cmpr_type_mod, only: cmpr_type
use diag_type_mod, only: diag_type, diag_list_type
use parcel_diags_type_mod, only: parcel_diags_type
use plume_model_diags_type_mod, only: plume_model_diags_type
use diags_super_type_mod, only: diags_super_type

implicit none

!----------------------------------------------------------------
! Type definition for a structure to store all diagnostics for
! a given type of draft (updraft, downdraft or fall-back flow)
!----------------------------------------------------------------
type :: draft_diags_type

  ! Number of diagnostics in this structure requested for output
  integer :: n_diags = 0

  ! List of pointers to active diagnostic structures
  ! (points to subset of the main list in the
  !  comorph_diags_type structure)
  type(diag_list_type), pointer :: list(:) => null()


  ! Updraft / downdraft parcel properties
  type(parcel_diags_type) :: par

  ! Initiation mass-source parcel properties
  type(parcel_diags_type) :: gen

  ! Diagnostics from the plume-model
  ! (entrainment, detrainment, moist processes etc)
  type(plume_model_diags_type) :: plume_model

end type draft_diags_type


!----------------------------------------------------------------
! Type definition containing super-arrays to store the above
! diagnostics in compressed form on each model-level
!----------------------------------------------------------------
type :: draft_diags_super_type
  type(diags_super_type), allocatable :: par(:,:,:)
  type(diags_super_type), allocatable :: gen(:,:,:)
  type(diags_super_type), allocatable :: plume_model(:,:,:)
end type draft_diags_super_type
! One of these structures is output by conv_sweep_ctl,
! containing all requested diagnostics in ragged-array-form
! for each convection type / layer and on each level.


contains


!----------------------------------------------------------------
! Subroutine to count the number of active diagnostics in a
! draft_diags structure and store a list of pointers to them
!----------------------------------------------------------------
! Note: this routine gets called twice
!   1st call (l_count_diags = .true.):
!     Check whether each diag is requested and count how many
!   2nd call (l_count_diags = .false.):
!     Set other properties for requested diags, and assign
!     pointers from active diags list.
subroutine draft_diags_assign( parent_name, l_count_diags,                     &
                               draft_diags, parent_list,                       &
                               parent_i_diag )

use comorph_constants_mod, only: name_length
use diag_type_mod, only: diag_list_type, dom_type
use parcel_diags_type_mod, only: parcel_diags_assign
use plume_model_diags_type_mod, only: plume_model_diags_assign

implicit none

! Character string to prepend onto all the diagnostic names
! in here; e.g. "updraft" or "downdraft"
character(len=name_length), intent(in) :: parent_name

! Flag for 1st call to this routine, in which we just check
! whether each diag is requested and count up number requested
logical, intent(in) :: l_count_diags

! The main diagnostics structure, containing pointers
! for all diagnostics, including those not in use
type(draft_diags_type), intent(in out) :: draft_diags

! List of active diagnostics in the parent structure
type(diag_list_type), target, intent(in out) :: parent_list(:)

! Counter for active diagnostics in the parent structure
integer, intent(in out) :: parent_i_diag

! Structure containing the allowed doms (2d,3d,4d) for the diags
type(dom_type) :: doms

! Name for diagnostics belonging to each parcel
character(len=name_length) :: par_name

! Counters
integer :: i_diag


! Set draft_diags list to point at a section of the parent list
if ( l_count_diags ) then
  ! Not used in the init call, but needs to be assigned
  ! to avoid error when passing it into diag_3d_assign
  draft_diags % list => parent_list
else
  ! Assign to subset address in parent:
  draft_diags % list => parent_list                                            &
    ( parent_i_diag + 1 :                                                      &
      parent_i_diag + draft_diags % n_diags )
end if

! Initialise local requested diagnostics counter to zero
i_diag = 0

! Process each diagnostic within the draft_diags structure...

! Set allowed domain profiles for the diagnostics;
! parcel properties, initiating parcel properties,
! and plume-model diagnostics can all be 3-D, 4-D or 5-D
! (mean over all types and layers, separate field for each type
!  but mean over layers, or separate field for every layer and type).
doms % x_y_z = .true.
doms % x_y_z_typ = .true.
doms % x_y_z_lay_typ = .true.

if ( draft_diags % par % n_diags > 0 .or. l_count_diags ) then
  ! Call routine to count number of active diags in
  ! draft_diags % par
  par_name = trim(adjustl(parent_name)) // "_par"
  call parcel_diags_assign( par_name, l_count_diags, doms,                     &
                            draft_diags % par,                                 &
                            draft_diags % list, i_diag )
end if

if ( draft_diags % gen % n_diags > 0 .or. l_count_diags ) then
  ! Call routine to count number of active diags in
  ! draft_diags % gen
  par_name = trim(adjustl(parent_name)) // "_gen"
  call parcel_diags_assign( par_name, l_count_diags, doms,                     &
                            draft_diags % gen,                                 &
                            draft_diags % list, i_diag )
end if

if ( draft_diags % plume_model % n_diags > 0 .or. l_count_diags ) then
  ! Call routine to count number of active diags in
  ! draft_diags % plume_model
  call plume_model_diags_assign( parent_name, l_count_diags, doms,             &
                                 draft_diags % plume_model,                    &
                                 draft_diags % list, i_diag )
end if


! Increment parent active diagnostics counter
parent_i_diag = parent_i_diag + i_diag


! If this is the initial call to this routine to set the
! number of active diags
if ( l_count_diags ) then

  ! Set number of active diags in this parcel_diags structure
  draft_diags % n_diags = i_diag

  ! Nullify the list pointer ready to be reassigned next call
  draft_diags % list => null()

end if


return
end subroutine draft_diags_assign


!----------------------------------------------------------------
! Subroutine to apply closure scalings to draft diagnostics
!----------------------------------------------------------------
subroutine draft_diags_scaling( n_conv_types, n_conv_layers,                   &
                                ij_first, ij_last, draft_scaling,              &
                                draft_diags, draft_diags_super )

use comorph_constants_mod, only: real_cvprec
use diags_super_type_mod, only: diags_super_scaling

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

! Structure storing meta-data for the diagnostics, and
! pointers to the output arrays
type(draft_diags_type), intent(in) :: draft_diags

! Structure containing compressed diags super-arrays
type(draft_diags_super_type), intent(in out) :: draft_diags_super

! Process parcel property diagnostics
if ( draft_diags % par % n_diags > 0 ) then
  call diags_super_scaling( n_conv_types, n_conv_layers,                       &
                            ij_first, ij_last, draft_scaling,                  &
                            draft_diags % par % n_diags_super,                 &
                            draft_diags % par % l_weight,                      &
                            draft_diags_super % par )
end if

! Process initiating mass-source property diagnostics
if ( draft_diags % gen % n_diags > 0 ) then
  call diags_super_scaling( n_conv_types, n_conv_layers,                       &
                            ij_first, ij_last, draft_scaling,                  &
                            draft_diags % gen % n_diags_super,                 &
                            draft_diags % gen % l_weight,                      &
                            draft_diags_super % gen )
end if

! Process plume-model diagnostics
if ( draft_diags % plume_model % n_diags > 0 ) then
  call diags_super_scaling( n_conv_types, n_conv_layers,                       &
                            ij_first, ij_last, draft_scaling,                  &
                            draft_diags % plume_model % n_diags_super,         &
                            draft_diags % plume_model % l_weight,              &
                            draft_diags_super % plume_model )
end if

return
end subroutine draft_diags_scaling


!----------------------------------------------------------------
! Subroutine to compute means of draft diagnostics over
! layers / types, and copy them into the output arrays
!----------------------------------------------------------------
subroutine draft_diags_compute_means( n_conv_types, n_conv_layers,             &
                                      max_points, ij_first, ij_last,           &
                                      draft_diags, draft_diags_super )

use diags_super_type_mod, only: diags_super_compute_means

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

! Structure storing meta-data for the diagnostics, and
! pointers to the output arrays
type(draft_diags_type), intent(in out) :: draft_diags

! Structure containing compressed diags super-arrays
type(draft_diags_super_type), intent(in out) :: draft_diags_super

! Process parcel property diagnostics
if ( draft_diags % par % n_diags > 0 ) then
  call diags_super_compute_means( n_conv_types, n_conv_layers,                 &
                                  max_points, ij_first, ij_last,               &
                                  draft_diags % par % n_diags,                 &
                                  draft_diags % par % n_diags_super,           &
                                  draft_diags % par % l_weight,                &
                                  draft_diags % par % list,                    &
                                  draft_diags_super % par,                     &
                                  fields_diags1 = draft_diags % par % mean,    &
                                  fields_diags2 = draft_diags % par % core )
end if

! Process initiating mass-source property diagnostics
if ( draft_diags % gen % n_diags > 0 ) then
  call diags_super_compute_means( n_conv_types, n_conv_layers,                 &
                                  max_points, ij_first, ij_last,               &
                                  draft_diags % gen % n_diags,                 &
                                  draft_diags % gen % n_diags_super,           &
                                  draft_diags % gen % l_weight,                &
                                  draft_diags % gen % list,                    &
                                  draft_diags_super % gen,                     &
                                  fields_diags1 = draft_diags % gen % mean,    &
                                  fields_diags2 = draft_diags % gen % core )
end if

! Process plume-model diagnostics
if ( draft_diags % plume_model % n_diags > 0 ) then
  call diags_super_compute_means( n_conv_types, n_conv_layers,                 &
                                  max_points, ij_first, ij_last,               &
                                  draft_diags % plume_model % n_diags,         &
                                  draft_diags % plume_model % n_diags_super,   &
                                  draft_diags % plume_model % l_weight,        &
                                  draft_diags % plume_model % list,            &
                                  draft_diags_super % plume_model,             &
                                  fields_diags1 = draft_diags                  &
                                                  % plume_model % ent_fields,  &
                                  fields_diags2 = draft_diags                  &
                                                  % plume_model % det_fields )
end if

! Deallocate the arrays of structures for compressed diags,
! which we have now processed and no-longer need
deallocate( draft_diags_super % plume_model )
deallocate( draft_diags_super % gen )
deallocate( draft_diags_super % par )

return
end subroutine draft_diags_compute_means


end module draft_diags_type_mod
