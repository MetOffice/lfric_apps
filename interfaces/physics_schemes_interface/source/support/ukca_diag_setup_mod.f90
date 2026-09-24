! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! -----------------------------------------------------------------------
!> @brief Module containing variables and routines for setting up UKCA diagnostics
!
! ----------------------------------------------------------------------

module ukca_diag_setup_mod

use ukca_api_mod,         only : ukca_maxlen_diagname, ukca_maxlen_message,    &
                                 ukca_diag_status_requested,                   &
                                 ukca_diag_status_inactive,                    &
                                 ukca_set_diagnostic_requests,                 &
                                 ukca_diagname_rxnflux_oh_ch4_trop,            &
                                 ukca_diagname_o3_column_du

use key_value_collection_mod,   only: key_value_collection_type
use lfric_xios_diag_mod,  only: field_is_active
use constants_mod,        only: imdi, i_def, l_def, str_def, i_um
use log_mod,              only: log_scratch_space, log_event, LOG_LEVEL_ERROR, & 
                                LOG_LEVEL_DEBUG
  
implicit none

private

integer(i_def), parameter :: n_diag_group = 2_i_def  ! Number of UKCA
                                                     ! diagnostic groups used

integer(i_def), parameter :: i_dgroup_2d = 1_i_def   ! Index used for
                                                   ! 2D group (not active yet)
integer(i_def), parameter :: i_dgroup_3d = 2_i_def   ! Index used for 
                                                     ! 3D group requests

integer(i_def), parameter :: max_ukca_diags = 2_i_def  ! Maximum number
                                                       ! of UKCA diagnostics
                                                       ! currently supported

integer(i_def), parameter :: n_req_max(n_diag_group) = [2_i_def, max_ukca_diags]
                             ! Max. no. of fields currently supported, by group

! Diagnostic short-names/ xios ids
character(len=ukca_maxlen_diagname), parameter, public ::   &
  nm_rxnflux_oh_ch4_trop = 'rxnflux_oh_ch4_trop'
character(len=ukca_maxlen_diagname), parameter, public ::   &
  nm_o3_column_du = 'o3_column_du'
character(len=ukca_maxlen_diagname), parameter, public ::   &
  key_diagname_base = 'ukca_req_diagnames_fullht_'

public :: ukca_diag_setup

contains

! ----------------------------------------------------------------------
subroutine ukca_diag_setup( values )
! ----------------------------------------------------------------------
! Description:
!   Set up the diagnostic request information be passed to UKCA
!   Note, currently this handles all defined diagnostics, but in future
!   it should only consider diagnostics requested via XIOS (active or inactive
!   on a given timesteps )

! ----------------------------------------------------------------------
implicit none

type(key_value_collection_type), intent(inout) :: values

! Local variables

! Dictionary mapping the short/ XIOS id of diagnostic to full name in UKCA
character(len=ukca_maxlen_diagname) :: ukca_diagnames_map(max_ukca_diags, 2)
  data ukca_diagnames_map(1,:) /nm_rxnflux_oh_ch4_trop,   &
                                ukca_diagname_rxnflux_oh_ch4_trop/
  data ukca_diagnames_map(2,:) /nm_o3_column_du, ukca_diagname_o3_column_du/

integer(i_def) :: n_req_ukca_diags_3d  ! Requested UKCA diagnostics (3D)

! Arrays to hold names and status flags of requested diagnostics, currently 3-D
character(len=ukca_maxlen_diagname), allocatable :: diagnames_fullht_real(:)
                                  ! In UKCA format (CF/ long-names), for API
character(len=ukca_maxlen_diagname), allocatable :: tmp_diagnames_fullht_real(:)
character(len=ukca_maxlen_diagname), allocatable ::  &   ! In XIOS id form
                                    req_diagnames_fullht(:)
integer(i_um), allocatable :: tmp_diag_status_3d(:)
integer(i_um), allocatable :: idiag_status_3d(:)

integer(i_def) :: i
logical :: l_diag_requested
character(len=ukca_maxlen_diagname):: diagname_key

! Error handling variables
integer(i_um) :: errcode
character(len=ukca_maxlen_message) :: ukca_errmsg

! End of header

errcode = 0_i_um
allocate(tmp_diagnames_fullht_real(n_req_max(i_dgroup_3d)))
allocate(req_diagnames_fullht(n_req_max(i_dgroup_3d)))
allocate(tmp_diag_status_3d(n_req_max(i_dgroup_3d)))
tmp_diagnames_fullht_real(:) = ''
req_diagnames_fullht(:) = ''
tmp_diag_status_3d(:) = ukca_diag_status_inactive

! Counter for active requests
n_req_ukca_diags_3d = 0_i_def

! Check if field is requested via XIOS configuration, irrespective of whether
! it is active on this timestep
do i = 1, n_req_max(i_dgroup_3d)
  if ( field_is_active('chemistry__'//ukca_diagnames_map(i, 1),               &
                        at_current_timestep=.false.) ) then
    n_req_ukca_diags_3d = n_req_ukca_diags_3d + 1
    tmp_diagnames_fullht_real(n_req_ukca_diags_3d) = ukca_diagnames_map(i, 2)
    req_diagnames_fullht(n_req_ukca_diags_3d) = ukca_diagnames_map(i, 1)
    tmp_diag_status_3d(n_req_ukca_diags_3d) = ukca_diag_status_requested
  end if  
end do

! Populate the allocatable arrays with the requested diagnostics
allocate(diagnames_fullht_real(n_req_ukca_diags_3d))
allocate(idiag_status_3d(n_req_ukca_diags_3d))
diagnames_fullht_real(:) = tmp_diagnames_fullht_real(1:n_req_ukca_diags_3d)
idiag_status_3d(:) = tmp_diag_status_3d(1:n_req_ukca_diags_3d)

deallocate(tmp_diag_status_3d)
deallocate(tmp_diagnames_fullht_real)

! Pass on active diagnostic information to UKCA via API - if any requested
if ( n_req_ukca_diags_3d > 0_i_def ) then
  CALL ukca_set_diagnostic_requests(                                           &
         errcode,                                                              &
         names_fullht_real=diagnames_fullht_real,                              &
         dreq_status_fullht_real=idiag_status_3d,                              &
         error_message=ukca_errmsg )
  if (errcode > 0_i_um)  then
    write(log_scratch_space, '(A,I0,A,A)')'Error in UKCA_SET_DIAG_REQUESTS: ', &
      errcode, ': ', TRIM(ukca_errmsg)
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if
end if

! Upload the list of requested diagnostics (xios_id) to the values dictionary
! NOTE: The key_values currently do not accept arrays, so store each value separately.
call values%add_key_value('n_req_ukca_diags_3d', n_req_ukca_diags_3d)
do i = 1, n_req_ukca_diags_3d
  write(diagname_key, '(A,I0)') trim(key_diagname_base), i
  call values%add_key_value(diagname_key, trim(req_diagnames_fullht(i)))
end do

deallocate(idiag_status_3d)
deallocate(diagnames_fullht_real)
deallocate(req_diagnames_fullht)

end subroutine ukca_diag_setup
! ----------------------------------------------------------------------

end module ukca_diag_setup_mod
