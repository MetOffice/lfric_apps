! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module comorph_diags_type_mod

use diag_type_mod, only: diag_type, diag_list_type
use fields_diags_type_mod, only: fields_diags_type
use draft_diags_type_mod, only: draft_diags_type
use diags_2d_type_mod, only: diags_2d_type
use subregion_diags_type_mod, only: subregion_diags_type
use subregion_mod, only: n_regions

implicit none

!---------------------------------------------------------------
! Type definition for a structure to store diagnostics.
!---------------------------------------------------------------
! This structure should be declared in the routine that calls
! comorph, and is in the comorp_ctl argument list for outputting
! the diagnostics.  Crucially, the contained pointers should be
! assigned to the desired memory addresses in the calling
! routine.
type :: comorph_diags_type

  ! Number of diagnostics in this structure requested for output
  integer :: n_diags = 0

  ! List of pointers to active diagnostic structures
  type(diag_list_type), allocatable :: list(:)


  ! Sub-structures containing meta-data for each diagnostic,
  ! and pointers to the output fields...


  ! Increments to primary fields:

  type(fields_diags_type) :: fields_incr


  ! Miscellanious diagnostics:

  ! Layer-masses used for conservation
  type(diag_type) :: layer_mass

  ! Lagrangian pressure increment (following environment parcels)
  ! due to compensating vertical mass rearrangement within
  ! the column.
  ! This is needed by the prognostic cloud scheme to calculate
  ! the homogeneous forcing of liquid cloud by the compensating
  ! subsidence term.
  type(diag_type) :: pressure_incr_env

  ! Total net upwards mass-flux (dry-mass),
  ! with negative contributions from downdrafts,
  ! on model-layer interfaces
  type(diag_type) :: massflux_d_net

  ! Total net upwards mass-flux (wet-mass)...
  type(diag_type) :: massflux_w_net

  ! Total convective cloud fractions
  type(diag_type) :: conv_cf_liq
  type(diag_type) :: conv_cf_ice
  type(diag_type) :: conv_cf_bulk

  ! Convective cloud liquid and ice mixing ratios
  ! (grid-mean values)
  type(diag_type) :: conv_q_cl
  type(diag_type) :: conv_q_cf

  ! Turbulence-based parcel perturbations
  type(fields_diags_type) :: turb_fields_pert
  ! Turbulence-based parcel radius
  type(diag_type) :: turb_radius

  ! Properties of each of the grid sub-regions;
  ! dry (no cloud or precep)
  ! liq (liquid-only cloud)
  ! mph (mixed-phased cloud)
  ! icr (ice and/or rain but no liquid cloud)
  type(subregion_diags_type) :: subregion_diags(n_regions)
  ! Number of requested diagnostics over all subregions
  integer :: n_diags_subregions = 0


  ! Diagnostics from the convective updrafts and downdrafts
  !
  ! There are separate diagnostic structures for the parcel
  ! ascent / descent for updrafts, downdraft,
  ! and their corresponding fall-back flows.
  !
  ! Where distinct convecting layers of the same type overlap
  ! in the vertical, the diagnostics are mass-flux-weighted means
  ! over the overlapping layers
  ! (apart from the contained mass-flux and area-fraction,
  !  which are summed over the layers).
  !
  ! Each contained diagnostic maybe output in either 4-D
  ! (separate field for each convection type) and/or 3-D
  ! (single field containing mass-flux-weighted mean over
  ! all convection types).

  ! Primary updraft
  type(draft_diags_type) :: updraft

  ! Primary downdraft
  type(draft_diags_type) :: dndraft

  ! Updraft fall-back flows
  type(draft_diags_type) :: updraft_fallback

  ! Downdraft fall-back flows
  type(draft_diags_type) :: dndraft_fallback


  ! 2D diagnostics (e.g. CAPE etc)
  !
  ! There are separate diagnostic structures for updrafts and downdrafts
  ! (each includes combined contributions from primary drafts and fall-backs)
  type(diags_2d_type) :: updraft_diags_2d
  type(diags_2d_type) :: dndraft_diags_2d

end type comorph_diags_type


contains


!----------------------------------------------------------------
! Subroutine to count the number of active diagnostics and
! store a list of pointers to them
!----------------------------------------------------------------
! Note: this routine gets called twice
!   1st call (l_count_diags = .true.):
!     Check whether each diag is requested and count how many
!   2nd call (l_count_diags = .false.):
!     Set other properties for requested diags, and assign
!     pointers from active diags list.
subroutine comorph_diags_assign( l_count_diags, comorph_diags )

use comorph_constants_mod, only: n_updraft_types, n_dndraft_types,             &
                     l_updraft_fallback, l_dndraft_fallback,                   &
                     name_length
use fields_diags_type_mod, only: fields_diags_assign
use draft_diags_type_mod, only: draft_diags_assign
use diags_2d_type_mod, only: diags_2d_assign
use subregion_mod, only: n_regions, region_names
use subregion_diags_type_mod, only: subregion_diags_assign
use diag_type_mod, only: diag_assign, dom_type
use init_diag_array_mod, only: init_diag_array

implicit none

! Flag for 1st call to this routine, in which we just check
! whether each diag is requested and count up number requested
logical, intent(in) :: l_count_diags

! The main diagnostics structure, containing pointers
! for all diagnostics, including those not in use
type(comorph_diags_type), intent(in out) :: comorph_diags

! Structure containing the allowed doms (2d,3d,4d) for the diags
type(dom_type) :: doms

! Name for each diagnostic
character(len=name_length) :: diag_name

! Counters
integer :: i_diag, i_super, i_region


if ( l_count_diags ) then
  ! Allocate the diags list to size 1 just to avoid illegal
  ! passing of an unallocated array into a subroutine
  allocate( comorph_diags % list(1) )
else
  ! Allocate the diags list to full size
  allocate( comorph_diags % list(comorph_diags % n_diags) )
end if

! Initialise requested diagnostics counter to zero
i_diag = 0


! Argument to diag_assign which is redundant here
i_super = 0

! Allowed domain profile for the following diags is 3-D only.
doms % x_y_z = .true.


! If any primary field increment diagnostics requested
! (or this is the initialisation call, in which case we
!  need to check whether any are requested)
if ( comorph_diags%fields_incr%n_diags > 0 .or. l_count_diags ) then
  ! Call routine to count number of active diags in fields_incr
  diag_name = "conv_incr"
  call fields_diags_assign( diag_name, l_count_diags, doms,                    &
                            comorph_diags % fields_incr,                       &
                            comorph_diags % list,                              &
                            i_diag, i_super )
end if


! Check whether various miscellaneous diagnostics are requested
diag_name = "layer_mass"
call diag_assign( diag_name, l_count_diags, doms,                              &
                  comorph_diags % layer_mass,                                  &
                  comorph_diags % list, i_diag, i_super )
diag_name = "pressure_incr_env"
call diag_assign( diag_name, l_count_diags, doms,                              &
                  comorph_diags % pressure_incr_env,                           &
                  comorph_diags % list, i_diag, i_super )
diag_name = "massflux_d_net"
call diag_assign( diag_name, l_count_diags, doms,                              &
                  comorph_diags % massflux_d_net,                              &
                  comorph_diags % list, i_diag, i_super )
diag_name = "massflux_w_net"
call diag_assign( diag_name, l_count_diags, doms,                              &
                  comorph_diags % massflux_w_net,                              &
                  comorph_diags % list, i_diag, i_super )


! Convective cloud diagnostics
diag_name = "conv_cf_liq"
call diag_assign( diag_name, l_count_diags, doms,                              &
                  comorph_diags % conv_cf_liq,                                 &
                  comorph_diags % list, i_diag, i_super )
diag_name = "conv_cf_ice"
call diag_assign( diag_name, l_count_diags, doms,                              &
                  comorph_diags % conv_cf_ice,                                 &
                  comorph_diags % list, i_diag, i_super )
diag_name = "conv_cf_bulk"
call diag_assign( diag_name, l_count_diags, doms,                              &
                  comorph_diags % conv_cf_bulk,                                &
                  comorph_diags % list, i_diag, i_super )
diag_name = "conv_q_cl"
call diag_assign( diag_name, l_count_diags, doms,                              &
                  comorph_diags % conv_q_cl,                                   &
                  comorph_diags % list, i_diag, i_super )
diag_name = "conv_q_cf"
call diag_assign( diag_name, l_count_diags, doms,                              &
                  comorph_diags % conv_q_cf,                                   &
                  comorph_diags % list, i_diag, i_super )


! If any turbulent perturbation diagnostics requested
! (or this is the initialisation call, in which case we
!  need to check whether any are requested)
if ( comorph_diags%turb_fields_pert%n_diags > 0                                &
     .or. l_count_diags ) then
  ! Call routine to count number of active diags in fields_incr
  diag_name = "turb_pert"
  call fields_diags_assign( diag_name, l_count_diags, doms,                    &
                            comorph_diags % turb_fields_pert,                  &
                            comorph_diags % list,                              &
                            i_diag, i_super )
end if

! Turbulent radius
diag_name = "turb_radius"
call diag_assign( diag_name, l_count_diags, doms,                              &
                  comorph_diags % turb_radius,                                 &
                  comorph_diags % list, i_diag, i_super )


! Fractions and fields for sub-regions of the grid-box
do i_region = 1, n_regions
  call subregion_diags_assign( region_names(i_region), l_count_diags,          &
         comorph_diags % subregion_diags(i_region),                            &
         comorph_diags % list,  i_diag )
  if ( l_count_diags ) then
    ! Sum requested diags over all subregions
    comorph_diags % n_diags_subregions                                         &
       = comorph_diags % n_diags_subregions                                    &
       + comorph_diags % subregion_diags(i_region) % n_diags
  end if
end do


! If any updraft or downdraft diagnostics are requested,
! (or this is the initialisation call, in which case we
!  need to check whether any are requested)...

! If any updraft types are in use
if ( n_updraft_types > 0 ) then

  if ( comorph_diags % updraft % n_diags > 0 .or. l_count_diags ) then
    ! Call to count number of active diags in updraft
    diag_name = "updraft"
    call draft_diags_assign( diag_name, l_count_diags,                         &
                             comorph_diags % updraft,                          &
                             comorph_diags % list, i_diag )
  end if

  ! If updraft fallbacks are on
  if ( l_updraft_fallback ) then
    if ( comorph_diags % updraft_fallback % n_diags > 0                        &
         .or. l_count_diags ) then
      ! Call to count number of active diags in updraft_fallback
      diag_name = "updraft_fallback"
      call draft_diags_assign( diag_name, l_count_diags,                       &
                               comorph_diags % updraft_fallback,               &
                               comorph_diags % list, i_diag )
    end if
  end if

  if ( comorph_diags % updraft_diags_2d % n_diags > 0 .or. l_count_diags ) then
    ! Call to count number of active diags in updraft_diags_2d
    diag_name = "updraft"
    call diags_2d_assign( diag_name, l_count_diags,                            &
                          comorph_diags % updraft_diags_2d,                    &
                          comorph_diags % list, i_diag )
  end if

end if  ! ( n_updraft_types > 0 )

! If any downdraft types are in use
if ( n_dndraft_types > 0 ) then

  if ( comorph_diags % dndraft % n_diags > 0 .or. l_count_diags ) then
    ! Call to count number of active diags in dndraft
    diag_name = "dndraft"
    call draft_diags_assign( diag_name, l_count_diags,                         &
                             comorph_diags % dndraft,                          &
                             comorph_diags % list, i_diag )
  end if

  ! If dndraft fallbacks are on
  if ( l_dndraft_fallback ) then
    if ( comorph_diags % dndraft_fallback % n_diags > 0                        &
         .or. l_count_diags ) then
      ! Call to count number of active diags in dndraft_fallback
      diag_name = "dndraft_fallback"
      call draft_diags_assign( diag_name, l_count_diags,                       &
                               comorph_diags % dndraft_fallback,               &
                               comorph_diags % list, i_diag )
    end if
  end if

  if ( comorph_diags % dndraft_diags_2d % n_diags > 0 .or. l_count_diags ) then
    ! Call to count number of active diags in dndraft_diags_2d
    diag_name = "dndraft"
    call diags_2d_assign( diag_name, l_count_diags,                            &
                          comorph_diags % dndraft_diags_2d,                    &
                          comorph_diags % list, i_diag )
  end if

end if  ! ( n_dndraft_types > 0 )


! If this is the initial call to this routine to set the
! number of active diags
if ( l_count_diags ) then

  ! Set it to the number of diags counted
  comorph_diags % n_diags = i_diag

  ! Deallocate the list ready to be allocated to full size
  ! at the next call
  deallocate( comorph_diags % list )

  ! If this is the 2nd call which has fully setup the list
else  ! ( l_count_diags )

  ! Now that the list is all setup, loop over all the requested
  ! diagnostics...
  do i_diag = 1, comorph_diags % n_diags

    ! Check that each diagnostic requested
    ! has a pointer pointing at some memory tested
    ! Its shape and extent must match, and if so
    ! The spell is cast; initialise to zero!
    call init_diag_array( comorph_diags % list(i_diag)%pt )

  end do  ! i_diag = 1, comorph_diags % n_diags

end if  ! ( l_count_diags )


return
end subroutine comorph_diags_assign


!---------------------------------------------------------------
! Routine to deallocate diags stuff at the end
!---------------------------------------------------------------
subroutine comorph_diags_dealloc( comorph_diags )

use comorph_constants_mod, only: n_updraft_types, n_dndraft_types,             &
                     l_updraft_fallback, l_dndraft_fallback

implicit none

! Diagnostic structure
type(comorph_diags_type), intent(in out) :: comorph_diags

! Deallocate condensed water species diagnostics pointer lists
if ( n_dndraft_types > 0 ) then
  if ( l_dndraft_fallback ) then
    deallocate( comorph_diags % dndraft_fallback % plume_model                 &
                % moist_proc % diags_cond )
  end if
  deallocate( comorph_diags % dndraft % plume_model                            &
              % moist_proc % diags_cond )
end if
if ( n_updraft_types > 0 ) then
  if ( l_updraft_fallback ) then
    deallocate( comorph_diags % updraft_fallback % plume_model                 &
                % moist_proc % diags_cond )
  end if
  deallocate( comorph_diags % updraft % plume_model                            &
              % moist_proc % diags_cond )
end if

! Deallocate the pointer list of active diagnostics
if ( comorph_diags % n_diags > 0 )                                             &
  deallocate( comorph_diags % list )
! (was only allocated if at least one diagnostic was requested).

return
end subroutine comorph_diags_dealloc


end module comorph_diags_type_mod
