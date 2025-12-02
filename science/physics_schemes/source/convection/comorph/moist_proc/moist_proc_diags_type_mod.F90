! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module moist_proc_diags_type_mod

use diag_type_mod, only: diag_type, diag_list_type

implicit none


! Type definition to store the diagnostics available for
! each condensed water species.
type :: diags_cond_type
  ! Increment due to sedimentation
  type(diag_type) :: dq_fall
  ! Increment due to condensation / evaporation
  type(diag_type) :: dq_cond
  ! Increment due to freezing and melting
  type(diag_type) :: dq_frzmlt
  ! Increment due to collision processes (accretion / riming)
  type(diag_type) :: dq_col
  ! Increment due to autoconversion
  type(diag_type) :: dq_aut
  ! Fall-speed
  type(diag_type) :: fall_speed
  ! Diagnosed hydrometeor temperature
  type(diag_type) :: cond_temp
end type diags_cond_type


! Type for defining a list of pointers to the above, assigned
! only for condensed water species that are in use.
type :: diags_cond_list_type
  type(diags_cond_type), pointer :: pt => null()
end type diags_cond_list_type


! Type for a structure to store diagnostics from the moist
! processes subroutine
type :: moist_proc_diags_type

  ! Number of diagnostics in this structure requested for output
  integer :: n_diags = 0

  ! List of pointers to active diagnostic structures
  ! (points to subset of the main list in the
  !  comorph_diags_type structure)
  type(diag_list_type), pointer :: list(:) => null()

  ! Diagnostics for each condensed water species
  type(diags_cond_type) :: diags_cl
  type(diags_cond_type) :: diags_rain
  type(diags_cond_type) :: diags_cf
  type(diags_cond_type) :: diags_snow
  type(diags_cond_type) :: diags_graup

  ! List of pointers to the above structures for only those
  ! condensed water species that are in use
  type(diags_cond_list_type), allocatable :: diags_cond(:)

  ! Note: moist_proc does not directly populate the full 3-D
  ! arrays of these diagnostics; rather it calculates the
  ! diags in a compressed super-array which is passed out
  ! and decompressed into the full 3-D arrays at a higher level.
  ! That super-array is not stored in here; it can't be,
  ! because for each instance of this structure, there
  ! might be multiple parallel calls to moist_proc within
  ! the OMP loop over convection segments.

end type moist_proc_diags_type


contains


!----------------------------------------------------------------
! Subroutine to count the number of active diagnostics in a
! moist_proc_diags structure and store a list of pointers to them
!----------------------------------------------------------------
! Note: this routine gets called twice
!   1st call (l_count_diags = .true.):
!     Check whether each diag is requested and count how many
!   2nd call (l_count_diags = .false.):
!     Set other properties for requested diags, and assign
!     pointers from active diags list.
subroutine moist_proc_diags_assign( parent_name, l_count_diags, doms,          &
                                    moist_proc_diags,                          &
                                    parent_list,                               &
                                    parent_i_diag, i_super )

use comorph_constants_mod, only: n_cond_species, cond_params,                  &
                     i_cond_cl, i_cond_rain,                                   &
                     i_cond_cf, i_cond_snow, i_cond_graup,                     &
                     l_cv_rain, l_cv_cf, l_cv_snow, l_cv_graup,                &
                     name_length
use diag_type_mod, only: diag_list_type, diag_assign, dom_type

implicit none

! Character string to prepend onto all the diagnostic names
! in here; e.g. "updraft" or "downdraft"
character(len=name_length), intent(in) :: parent_name

! Flag for 1st call to this routine, in which we just check
! whether each diag is requested and count up number requested
logical, intent(in) :: l_count_diags

! Structure containing the allowed domain profiles for the diags
type(dom_type), intent(in) :: doms
! Note: all diags in the moist_proc_diags structure are assumed
! to have the same allowed domain profiles.

! The main diagnostics structure, containing pointers
! for all diagnostics, including those not in use
type(moist_proc_diags_type), target, intent(in out) ::                         &
                                 moist_proc_diags

! List of active diagnostics in the parent structure
type(diag_list_type), target, intent(in out) :: parent_list(:)

! Counter for active diagnostics in the parent structure
integer, intent(in out) :: parent_i_diag

! Counter for addresses of fields in compressed diags super-array
integer, intent(in out) :: i_super

! Name for each diagnostic
character(len=name_length) :: name_tmp
character(len=name_length) :: diag_name

! Counters
integer :: i_diag, i_cond


! Setup list of cond_diag structures to point at those
! condensed water species that are actually in use
if ( l_count_diags ) then
  allocate( moist_proc_diags % diags_cond(n_cond_species) )
  moist_proc_diags % diags_cond(i_cond_cl)%pt                                  &
                     => moist_proc_diags % diags_cl
  if ( l_cv_rain )   moist_proc_diags                                          &
                   % diags_cond(i_cond_rain)%pt                                &
                     => moist_proc_diags % diags_rain
  if ( l_cv_cf )     moist_proc_diags                                          &
                   % diags_cond(i_cond_cf)%pt                                  &
                     => moist_proc_diags % diags_cf
  if ( l_cv_snow )   moist_proc_diags                                          &
                   % diags_cond(i_cond_snow)%pt                                &
                     => moist_proc_diags % diags_snow
  if ( l_cv_graup )  moist_proc_diags                                          &
                   % diags_cond(i_cond_graup)%pt                               &
                     => moist_proc_diags % diags_graup
end if


! Set moist_proc_diags list to point at a section of parent list
if ( l_count_diags ) then
  ! Not used in the init call, but needs to be assigned
  ! to avoid error when passing it into diag_3d_assign
  moist_proc_diags % list => parent_list
else
  ! Assign to subset address in parent:
  moist_proc_diags % list => parent_list                                       &
    ( parent_i_diag + 1 :                                                      &
      parent_i_diag + moist_proc_diags % n_diags )
end if

! Initialise active 3-D diagnostic counter to zero
i_diag = 0

! Loop over condensed water species
do i_cond = 1, n_cond_species
  ! Count number of diags requested for this species...

  ! Construct name of q increment
  name_tmp = adjustl( trim(adjustl(parent_name)) // "_dq_" //                  &
                  trim(adjustl(cond_params(i_cond)%pt % cond_name)) )

  ! Sedimentation increment
  diag_name = trim(name_tmp) // "_fall"
  call diag_assign( diag_name, l_count_diags, doms, moist_proc_diags           &
                    % diags_cond(i_cond)%pt % dq_fall,                         &
                    moist_proc_diags%list, i_diag, i_super )

  ! Condensation / evaporation increment
  diag_name = trim(name_tmp) // "_cond"
  call diag_assign( diag_name, l_count_diags, doms, moist_proc_diags           &
                    % diags_cond(i_cond)%pt % dq_cond,                         &
                    moist_proc_diags%list, i_diag, i_super )

  ! Freezing / melting increment
  diag_name = trim(name_tmp) // "_frzmlt"
  call diag_assign( diag_name, l_count_diags, doms, moist_proc_diags           &
                    % diags_cond(i_cond)%pt % dq_frzmlt,                       &
                    moist_proc_diags%list, i_diag, i_super )

  ! Collision process increment
  diag_name = trim(name_tmp) // "_col"
  call diag_assign( diag_name, l_count_diags, doms, moist_proc_diags           &
                    % diags_cond(i_cond)%pt % dq_col,                          &
                    moist_proc_diags%list, i_diag, i_super )

  ! Autoconversion increment
  diag_name = trim(name_tmp) // "_aut"
  call diag_assign( diag_name, l_count_diags, doms, moist_proc_diags           &
                    % diags_cond(i_cond)%pt % dq_aut,                          &
                    moist_proc_diags%list, i_diag, i_super )

  ! Construct name to prepend onto other field
  name_tmp = adjustl( trim(adjustl(parent_name)) // "_" //                     &
                  trim(adjustl(cond_params(i_cond)%pt % cond_name)) )

  ! Fall-speed
  diag_name = trim(name_tmp) // "_fall_speed"
  call diag_assign( diag_name, l_count_diags, doms, moist_proc_diags           &
                    % diags_cond(i_cond)%pt % fall_speed,                      &
                    moist_proc_diags%list, i_diag, i_super )

  ! Hydrometeor temperature
  diag_name = trim(name_tmp) // "_temp"
  call diag_assign( diag_name, l_count_diags, doms, moist_proc_diags           &
                    % diags_cond(i_cond)%pt % cond_temp,                       &
                    moist_proc_diags%list, i_diag, i_super )

end do  ! i_cond = 1, n_cond_species


! Increment parent active diagnostics counter
parent_i_diag = parent_i_diag + i_diag


! If this is the initial call to this routine to count the
! number of active diags
if ( l_count_diags ) then

  ! Set number of active diags in this moist_proc_diags structure
  moist_proc_diags % n_diags = i_diag

  ! Nullify the list pointer ready to be reassigned next call
  moist_proc_diags % list => null()

end if


return
end subroutine moist_proc_diags_assign


end module moist_proc_diags_type_mod
