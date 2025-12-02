! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module parcel_diags_type_mod

use diag_type_mod, only: diag_type, diag_list_type
use fields_diags_type_mod, only: fields_diags_type

implicit none

!----------------------------------------------------------------
! Type definition for a structure to store diagnostics from
! the updraft and downdraft parcels
!----------------------------------------------------------------
! This follows the structure of parcel_type in parcel_type_mod.

! Note: during the main ascent / descent, these diagnostics are
! valid at the next half-level, i.e. k+1/2 for updrafts,
! k-1/2 for downdrafts.
type :: parcel_diags_type

  ! Number of diagnostics in this structure requested for output
  integer :: n_diags = 0

  ! Number of diagnostics that need to be stored in the
  ! compressed super-array, accounting for any fields that
  ! aren't requested but are needed for computing other diags
  integer :: n_diags_super = 0

  ! Number of active diagnostics from only the below 4 fields
  ! (i.e. number that correspond to fields in the parcel_type
  !  super-array)
  integer :: n_diags_par = 0

  ! Flags for which fields in the super-array to use as mass-weights
  ! for computing means of other fields
  logical, allocatable :: l_weight(:)

  ! List of pointers to active diagnostic structures
  ! (points to subset of the main list in the
  !  comorph_diags_type structure)
  type(diag_list_type), pointer :: list(:) => null()


  ! Dry-mass flux / kg m-2 s-1
  ! (kg per s per unit surface area)
  type(diag_type) :: massflux_d
  ! Parcel radius / m
  type(diag_type) :: radius
  ! Edge virtual temperature / K
  type(diag_type) :: edge_virt_temp

  ! In-plume means of primary fields
  type(fields_diags_type) :: mean

  ! In-plume core properties
  type(fields_diags_type) :: core

end type parcel_diags_type


contains


!----------------------------------------------------------------
! Subroutine to count the number of active diagnostics in a
! parcel structure and store a list of pointers to them
!----------------------------------------------------------------
! Note: this routine gets called twice
!   1st call (l_count_diags = .true.):
!     Check whether each diag is requested and count how many
!   2nd call (l_count_diags = .false.):
!     Set other properties for requested diags, and assign
!     pointers from active diags list.
subroutine parcel_diags_assign( parent_name, l_count_diags, doms,              &
                                parcel_diags, parent_list,                     &
                                parent_i_diag )

use comorph_constants_mod, only: name_length, newline, l_par_core
use fields_diags_type_mod, only: fields_diags_assign
use diag_type_mod, only: diag_list_type, diag_assign, dom_type
use parcel_type_mod, only: i_massflux_d, i_radius, i_edge_virt_temp
use raise_error_mod, only: raise_fatal

implicit none

! Character string to prepend onto all the diagnostic names
! in here; e.g. "updraft" or "downdraft"
character(len=name_length), intent(in) :: parent_name

! Flag for 1st call to this routine, in which we just check
! whether each diag is requested and count up number requested
logical, intent(in) :: l_count_diags

! Structure containing the allowed domain profiles for the diags
type(dom_type), intent(in) :: doms
! Note: all diags in the parcel_diags structure are assumed
! to have the same allowed domain profiles.

! The main diagnostics structure, containing diag structures
! for all diagnostics, including those not in use
type(parcel_diags_type), intent(in out) :: parcel_diags

! List of active diagnostics in the parent structure
type(diag_list_type), target, intent(in out) :: parent_list(:)

! Counter for active diagnostics in the parent structure
integer, intent(in out) :: parent_i_diag

! Counter for addresses of fields in compressed diags super-array
integer :: i_super

! Flag to indicate that means of any requested primary fields
! need to be computed, averaging over multiple layers
! or types of convection (passed into fields_diags_assign)
logical, parameter :: l_mean_true = .true.

! Name for each diagnostic
character(len=name_length) :: diag_name

! Counter
integer :: i_diag

character(len=*), parameter :: routinename                                     &
                               = "PARCEL_DIAGS_ASSIGN"


! Set parcel_diags list to point at a section of the parent list
if ( l_count_diags ) then
  ! Not used in the init call, but needs to be assigned
  ! to avoid error when passing it into diag_3d_assign
  parcel_diags % list => parent_list
else
  ! Assign to subset address in parent:
  parcel_diags % list => parent_list                                           &
    ( parent_i_diag + 1 :                                                      &
      parent_i_diag + parcel_diags % n_diags )
  ! Initialise counter for super-array addresses
  i_super = 0
end if

! Initialise local requested diagnostics counter to zero
i_diag = 0

! Process each diagnostic within the parcel_diags structure...

if ( parcel_diags % n_diags > 0 .or. l_count_diags ) then

  diag_name = trim(adjustl(parent_name)) // "_massflux_d"
  call diag_assign( diag_name, l_count_diags, doms,                            &
                    parcel_diags % massflux_d,                                 &
                    parcel_diags % list, i_diag, i_super,                      &
                    i_field=i_massflux_d )
  diag_name = trim(adjustl(parent_name)) // "_radius"
  call diag_assign( diag_name, l_count_diags, doms,                            &
                    parcel_diags % radius,                                     &
                    parcel_diags % list, i_diag, i_super,                      &
                    i_field=i_radius )
  diag_name = trim(adjustl(parent_name)) // "_edge_virt_temp"
  call diag_assign( diag_name, l_count_diags, doms,                            &
                    parcel_diags % edge_virt_temp,                             &
                    parcel_diags % list, i_diag, i_super,                      &
                    i_field=i_edge_virt_temp )

end if

! Save count of number of the above parcel super-array fields
! which have diagnostics requested
if ( l_count_diags )  parcel_diags % n_diags_par = i_diag


! If any parcel mean primary field diagnostics requested
! (or this is the first counting call, in which case we
!  need to check whether any are requested)
if ( parcel_diags % mean % n_diags > 0 .or. l_count_diags ) then
  ! Call routine to count number of active diags in
  ! parcel_diags % mean
  diag_name = trim(adjustl(parent_name)) // "_mean"
  call fields_diags_assign( diag_name, l_count_diags, doms,                    &
                            parcel_diags%mean, parcel_diags%list,              &
                            i_diag, i_super, l_mean=l_mean_true )
end if

! Ditto for parcel core properties
if ( parcel_diags % core % n_diags > 0 .or. l_count_diags ) then
  ! Call routine to count number of active diags in
  ! parcel_diags % core
  diag_name = trim(adjustl(parent_name)) // "_core"
  call fields_diags_assign( diag_name, l_count_diags, doms,                    &
                            parcel_diags%core, parcel_diags%list,              &
                            i_diag, i_super, l_mean=l_mean_true )
end if

! Fatal error if parcel core diagnostics are requested,
! but parcel core properties are not in use.
if ( parcel_diags % core % n_diags > 0                                         &
     .and. (.not. l_par_core ) ) then
  call raise_fatal( routinename,                                               &
         "Diagnostics have been requested for parcel core "   //               &
         "properties, but l_par_core is set to "     //newline//               &
         "false, so the parcel core properties don't exist!" )
end if

! Increment parent active diagnostics counter
parent_i_diag = parent_i_diag + i_diag


! If this is the initial call to this routine to count the
! number of active diags
if ( l_count_diags ) then

  ! Set number of active diags in this parcel_diags structure
  parcel_diags % n_diags = i_diag

  ! Nullify the list pointer ready to be reassigned next call
  parcel_diags % list => null()

  ! Set flag true for diags needed for calculating other diags

  ! If any other parcel diags were requested, then the
  ! mass-flux is also needed for computing mass-flux-weighted
  ! means over multiple convecting layers.
  if ( parcel_diags % n_diags > 0 ) then
    parcel_diags % massflux_d % flag = .true.
  end if

  ! Note: the 2nd sweep with l_count_diags=.false. will now
  ! account for flag when setting super-array addresses
  ! in diag_assign.

else  ! ( l_count_diags )

  ! End of 2nd sweep; set final count for number of diags
  ! in the compressed super-array
  parcel_diags % n_diags_super = i_super

  ! Set flags for whether each field in the super-array is a mass-weight
  if ( allocated( parcel_diags % l_weight ) ) then
    deallocate( parcel_diags % l_weight )
  end if
  allocate( parcel_diags % l_weight( i_super ) )
  do i_super = 1, parcel_diags % n_diags_super
    ! Initialise to false everywhere
    parcel_diags % l_weight(i_super) = .false.
  end do
  ! Set to true for dry-mass-flux
  if ( parcel_diags % massflux_d % flag ) then
    i_super = parcel_diags % massflux_d % i_super
    parcel_diags % l_weight(i_super) = .true.
  end if
  ! At present, diags_super_compute_means (which performs the averaging)
  ! assumes that the weight comes before all the fields which
  ! need to be averaged in the super-array.  Therefore, if one of
  ! the above fields isn't first (or isn't there), throw a wobbly
  if ( .not. parcel_diags % l_weight(1) ) then
    call raise_fatal( routinename,                                             &
      "Error: trying to compute means of compressed "           //             &
      "parcel diagnostics, but the mass-flux "         //newline//             &
      "field needed for use as the weight in the averaging "    //             &
      "seems to have been misplaced." )
  end if

end if  ! ( l_count_diags )


return
end subroutine parcel_diags_assign


!----------------------------------------------------------------
! Subroutine to copy parcel properties into a super-array
! for diagnostic output
!----------------------------------------------------------------
subroutine parcel_diags_copy( n_fields_tot,                                    &
                              parcel_diags, parcel, pressure, diags_super )

use comorph_constants_mod, only: real_cvprec, l_par_core
use parcel_type_mod, only: parcel_type, i_massflux_d, i_edge_virt_temp
use fields_diags_type_mod, only: fields_diags_copy,                            &
                                 fields_diags_conserved_vars

implicit none

! Number of fields in parcel % fields_super
integer, intent(in) :: n_fields_tot

! Structure containing meta-data for parcel diags
type(parcel_diags_type), intent(in) :: parcel_diags

! Input structure storing parcel properties to be copied
type(parcel_type), intent(in) :: parcel

! Pressure, needed for computing RH diagnostic
real(kind=real_cvprec), intent(in) :: pressure( parcel%cmpr%n_points )

! Diagnostics super-array to copy fields into
real(kind=real_cvprec), intent(in out) :: diags_super                          &
            ( parcel%cmpr%n_points, parcel_diags%n_diags_super )

! Size of the parcel super-arrays (maybe bigger than needed,
! to save having to deallocate and reallocate)
integer :: n_points_super

! Flag input to fields_diags_copy
logical :: l_conserved_form

! Flag to pass into fields_diags_conserved_vars
logical, parameter :: l_reverse = .false.

! Loop counters
integer :: ic, i_req, i_diag, i_field

! Loop over requested parcel fields
! (excluding the primary fields)
do i_req = 1, parcel_diags % n_diags_par
  ! Extract diag-array index and fields-array index
  i_diag  = parcel_diags % list(i_req)%pt % i_super
  i_field = parcel_diags % list(i_req)%pt % i_field
  ! Copy field
  do ic = 1, parcel % cmpr % n_points
    diags_super(ic,i_diag) = parcel%par_super(ic,i_field)
  end do
end do

! Dry-mass flux might be required even if not requested
if ( parcel_diags % massflux_d % flag ) then
  if ( .not. parcel_diags % massflux_d % l_req ) then
    i_diag = parcel_diags % massflux_d % i_super
    do ic = 1, parcel % cmpr % n_points
      diags_super(ic,i_diag) = parcel%par_super(ic,i_massflux_d)
    end do
  end if
end if


n_points_super = size( parcel % mean_super, 1 )

! Copy parcel-mean primary fields
! (note: these are already in conserved variable form)
if ( parcel_diags % mean % n_diags > 0 ) then
  l_conserved_form = .true.
  call fields_diags_copy( parcel%cmpr%n_points, n_points_super,                &
                          parcel%cmpr%n_points, parcel_diags%n_diags_super,    &
                          n_fields_tot, l_conserved_form,                      &
                          parcel_diags%mean, parcel%mean_super,                &
                          parcel%par_super(:,i_edge_virt_temp), pressure,      &
                          diags_super )
end if

! Copy parcel core primary fields, if in use
if ( parcel_diags % core % n_diags > 0 .and. l_par_core ) then
  l_conserved_form = .false.
  call fields_diags_copy( parcel%cmpr%n_points, n_points_super,                &
                          parcel%cmpr%n_points, parcel_diags%n_diags_super,    &
                          n_fields_tot, l_conserved_form,                      &
                          parcel_diags%core, parcel%core_super,                &
                          parcel%par_super(:,i_edge_virt_temp), pressure,      &
                          diags_super )
  ! Convert core field properties to conserved variable form
  ! ready for calculating means over types and layers.
  call fields_diags_conserved_vars( parcel%cmpr%n_points,                      &
                                    parcel%cmpr%n_points,                      &
                                    parcel_diags%n_diags_super,                &
                                    l_reverse, parcel_diags%core,              &
                                    diags_super )
end if

return
end subroutine parcel_diags_copy


end module parcel_diags_type_mod
