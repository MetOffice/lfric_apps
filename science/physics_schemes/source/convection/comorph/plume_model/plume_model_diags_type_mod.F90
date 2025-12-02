! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module plume_model_diags_type_mod

use comorph_constants_mod, only: real_cvprec
use diag_type_mod, only: diag_type, diag_list_type
use moist_proc_diags_type_mod, only: moist_proc_diags_type
use fields_diags_type_mod, only: fields_diags_type

implicit none

!----------------------------------------------------------------
! Type definition for a structure to store diagnostics from
! the updraft and downdraft diagnostic ascent / descent model.
!----------------------------------------------------------------

! Note: these diagnostics are valid at the full-level k, i.e.
! over the updraft level-step from k-1/2 to k+1/2 (or the reverse
! for downdrafts).
type :: plume_model_diags_type

  ! Number of diagnostics in this structure requested for output
  integer :: n_diags = 0

  ! Number of diagnostics that need to be stored in the
  ! compressed super-array, accounting for any fields that
  ! aren't requested but are needed for computing other diags
  integer :: n_diags_super = 0

  ! Flags for which fields in the super-array to use as mass-weights
  ! for computing means of other fields
  logical, allocatable :: l_weight(:)

  ! List of pointers to active diagnostic structures
  ! (points to subset of the main list in the
  !  comorph_diags_type structure)
  type(diag_list_type), pointer :: list(:) => null()

  ! Dry-mass-flux in which the moist_proc diags occur
  ! (currently the same as the mass-flux at the next half-level)
  type(diag_type) :: massflux_d_k

  ! Time taken for parcel to traverse the current model-level
  type(diag_type) :: delta_t

  ! Moist process diagnostics
  type(moist_proc_diags_type) :: moist_proc

  ! Entrained dry-mass and air properties
  type(diag_type) :: ent_mass_d
  type(fields_diags_type) :: ent_fields

  ! Detrained dry-mass and air properties
  type(diag_type) :: det_mass_d
  type(fields_diags_type) :: det_fields

  ! Ratio between core and mean parcel properties, used for detrainment calc.
  type(diag_type) :: core_mean_ratio

  ! Note: conv_level_step does not directly populate the full 3-D
  ! arrays of these diagnostics; rather it calculates the
  ! diags in a compressed super-array which is passed out
  ! and decompressed into the full 3-D arrays at a higher level.
  ! That super-array is not stored in here; it can't be,
  ! because for each instance of this structure, there
  ! might be multiple parallel calls to conv_level_step within
  ! the OMP loop over convection segments.

end type plume_model_diags_type


contains


!----------------------------------------------------------------
! Subroutine to count the number of active diagnostics in a
! plume_model_diags structure / store a list of pointers to them
!----------------------------------------------------------------
! Note: this routine gets called twice
!   1st call (l_count_diags = .true.):
!     Check whether each diag is requested and count how many
!   2nd call (l_count_diags = .false.):
!     Set other properties for requested diags, and assign
!     pointers from active diags list.
subroutine plume_model_diags_assign( parent_name, l_count_diags, doms,         &
                                     plume_model_diags,                        &
                                     parent_list,                              &
                                     parent_i_diag )

use comorph_constants_mod, only: name_length, newline
use diag_type_mod, only: diag_list_type, diag_assign, dom_type
use moist_proc_diags_type_mod, only: moist_proc_diags_assign
use fields_diags_type_mod, only: fields_diags_assign
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
! Note: all diags in the plume_model_diags structure are assumed
! to have the same allowed domain profiles.

! The main diagnostics structure, containing pointers
! for all diagnostics, including those not in use
type(plume_model_diags_type), target, intent(in out) ::                        &
                                 plume_model_diags

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
                               = "PLUME_MODEL_DIAGS_ASSIGN"


! Set plume_model_diags list to point at a section of parent list
if ( l_count_diags ) then
  ! Not used in the init call, but needs to be assigned
  ! to avoid error when passing it into diag_3d_assign
  plume_model_diags % list => parent_list
else
  ! Assign to subset address in parent:
  plume_model_diags % list => parent_list                                      &
    ( parent_i_diag + 1 :                                                      &
      parent_i_diag + plume_model_diags % n_diags )
  ! Initialise counter for super-array addresses
  i_super = 0
end if

! Initialise local requested diagnostics counter to zero
i_diag = 0

! Process each diagnostic in the plume_model_diags structure...

! Flux of dry-mass at k in which moist processes occur
diag_name = trim(adjustl(parent_name)) // "_massflux_d_k"
call diag_assign( diag_name, l_count_diags, doms,                              &
                  plume_model_diags % massflux_d_k,                            &
                  plume_model_diags % list, i_diag, i_super )

! Time interval for level-step
diag_name = trim(adjustl(parent_name)) // "_delta_t"
call diag_assign( diag_name, l_count_diags, doms,                              &
                  plume_model_diags % delta_t,                                 &
                  plume_model_diags % list, i_diag, i_super )

! Moist process diagnostics
if ( plume_model_diags % moist_proc % n_diags > 0                              &
     .or. l_count_diags ) then
  call moist_proc_diags_assign( parent_name, l_count_diags, doms,              &
                                plume_model_diags % moist_proc,                &
                                plume_model_diags % list,                      &
                                i_diag, i_super )
end if

! Entrained and detrained mass and air properties
diag_name = trim(adjustl(parent_name)) // "_ent_mass_d"
call diag_assign( diag_name, l_count_diags, doms,                              &
                  plume_model_diags % ent_mass_d,                              &
                  plume_model_diags % list, i_diag, i_super )
diag_name = trim(adjustl(parent_name)) // "_ent"
call fields_diags_assign( diag_name, l_count_diags, doms,                      &
                          plume_model_diags % ent_fields,                      &
                          plume_model_diags % list,                            &
                          i_diag, i_super, l_mean_true )
diag_name = trim(adjustl(parent_name)) // "_det_mass_d"
call diag_assign( diag_name, l_count_diags, doms,                              &
                  plume_model_diags % det_mass_d,                              &
                  plume_model_diags % list, i_diag, i_super )
diag_name = trim(adjustl(parent_name)) // "_det"
call fields_diags_assign( diag_name, l_count_diags, doms,                      &
                          plume_model_diags % det_fields,                      &
                          plume_model_diags % list,                            &
                          i_diag, i_super, l_mean_true )

! Ratio between core and mean parcel properties
diag_name = trim(adjustl(parent_name)) // "_core_mean_ratio"
call diag_assign( diag_name, l_count_diags, doms,                              &
                  plume_model_diags % core_mean_ratio,                         &
                  plume_model_diags % list, i_diag, i_super )


! Increment parent active diagnostics counter
parent_i_diag = parent_i_diag + i_diag


! If this is the initial call to this routine to set the
! number of active diags
if ( l_count_diags ) then

  ! Set number of active diags in this moist_proc_diags structure
  plume_model_diags % n_diags = i_diag

  ! Nullify the list pointer ready to be reassigned next call
  plume_model_diags % list => null()

  ! The mass-flux at k after entrainment is needed for computing
  ! mass-weighted-means of the moist process diagnostics,
  ! and so needs to be assigned an address in the compressed
  ! super-array even if not requested for output...
  if ( plume_model_diags % moist_proc % n_diags > 0 .or.                       &
       plume_model_diags % delta_t % flag .or.                                 &
       plume_model_diags % core_mean_ratio % flag ) then
    plume_model_diags % massflux_d_k % flag = .true.
  end if

  ! The entrained mass is needed for computing mean entrained
  ! air properties
  if ( plume_model_diags % ent_fields % n_diags > 0 ) then
    plume_model_diags % ent_mass_d % flag = .true.
  end if

  ! Ditto for the detrained mass
  if ( plume_model_diags % det_fields % n_diags > 0 ) then
    plume_model_diags % det_mass_d % flag = .true.
  end if

  ! Note: the 2nd sweep with l_count_diags=.false. will now
  ! account for flag when setting super-array addresses
  ! in diag_assign.

else  ! ( l_count_diags )

  ! End of 2nd sweep; set final count for number of diags
  ! in the compressed super-array
  plume_model_diags % n_diags_super = i_super

  ! Set flags for whether each field in the super-array is a mass-weight
  if ( allocated( plume_model_diags % l_weight ) ) then
    deallocate( plume_model_diags % l_weight )
  end if
  allocate( plume_model_diags % l_weight( i_super ) )
  do i_super = 1, plume_model_diags % n_diags_super
    ! Initialise to false everywhere
    plume_model_diags % l_weight(i_super) = .false.
  end do
  ! Set to true for dry-mass-flux
  if ( plume_model_diags % massflux_d_k % flag ) then
    i_super = plume_model_diags % massflux_d_k % i_super
    plume_model_diags % l_weight(i_super) = .true.
  end if
  ! Set to true for entrained dry-mass
  if ( plume_model_diags % ent_mass_d % flag ) then
    i_super = plume_model_diags % ent_mass_d % i_super
    plume_model_diags % l_weight(i_super) = .true.
  end if
  ! Set true for detrained dry-mass
  if ( plume_model_diags % det_mass_d % flag ) then
    i_super = plume_model_diags % det_mass_d % i_super
    plume_model_diags % l_weight(i_super) = .true.
  end if
  ! At present, diags_super_compute_means (which performs the averaging)
  ! assumes that the weight comes before all the fields which
  ! need to be averaged in the super-array.  Therefore, if one of
  ! the above fields isn't first (or isn't there), throw a wobbly
  if ( .not. plume_model_diags % l_weight(1) ) then
    call raise_fatal( routinename,                                             &
      "Error: trying to compute means of compressed "           //             &
      "plume-model diagnostics, but the mass-flux "    //newline//             &
      "field needed for use as the weight in the averaging "    //             &
      "seems to have been misplaced." )
  end if

end if  ! ( l_count_diags )


return
end subroutine plume_model_diags_assign


end module plume_model_diags_type_mod
