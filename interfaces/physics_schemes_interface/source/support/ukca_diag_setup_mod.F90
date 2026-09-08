! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
! -----------------------------------------------------------------------
!> @brief Module containing variables and routines for setting up UKCA diagnostics
!
! ----------------------------------------------------------------------

module ukca_diag_setup_mod

use ukca_api_mod, only :                                                       &    
    ukca_maxlen_diagname, ukca_maxlen_message,                                 &
    ukca_diag_status_requested, ukca_set_diagnostic_requests,                  &
    ukca_diagname_rxnflux_oh_ch4_trop,                                         &
    ukca_diagname_o3_column_du    

use constants_mod, only : imdi, i_def, l_def, str_def, i_um
use log_mod,       only : log_scratch_space, log_event, LOG_LEVEL_ERROR,    & 
                          LOG_LEVEL_DEBUG
  
implicit none

private

public :: ukca_diag_setup

integer(i_def), parameter, public :: n_diag_group = 2_i_def  ! Number of UKCA
                                                     ! diagnostic groups used

integer(i_def), parameter, public :: i_dgroup_2d = 1_i_def   ! Index used for
                                                   ! 2D group (not active yet)
integer(i_def), parameter, public :: i_dgroup_3d = 2_i_def   ! Index used for 
                                                     ! 3D group requests

integer(i_def), parameter :: max_ukca_diags = 2_i_def  ! Maximum number of UKCA 
                                                       ! diagnostics currently
                                                       ! supported for output

integer(i_def), parameter :: n_req_max(n_diag_group) = [2_i_def, max_ukca_diags]
                             ! Max. no. of fields currently supported, by group

! Diagnostic short-names/ xios ids
character(len=ukca_maxlen_diagname), parameter, public ::   &
  nm_rxnflux_oh_ch4_trop = 'rxnflux_oh_ch4_trop'
character(len=ukca_maxlen_diagname), parameter, public ::   &
  nm_o3_column_du = 'o3_column_du'

! Dictionary mapping the short/ XIOS id of diagnostic to full name in UKCA
character(len=ukca_maxlen_diagname), public :: diagnames_map(max_ukca_diags, 2)
  data diagnames_map(1,:) /nm_rxnflux_oh_ch4_trop,   &
                                ukca_diagname_rxnflux_oh_ch4_trop/
  data diagnames_map(2,:) /nm_o3_column_du, ukca_diagname_o3_column_du/

integer(i_def), public :: n_ukca_diags_3d  ! Counter for active requests

! Names of requested diagnostics -currently only 3D
character(len=ukca_maxlen_diagname), allocatable, public :: diagnames_fullht_real(:)

! Integer status flag for diagnostics to match API arguments
integer(i_um), allocatable, public :: idiag_status_3d(:)

contains

! ----------------------------------------------------------------------
subroutine ukca_diag_setup( )
! ----------------------------------------------------------------------
! Description:
!   Set up the diagnostic request information be passed to UKCA
!   Note, currently this handles all defined diagnostics, but in future
!   it should only consider diagnostics requested via XIOS (active or inactive
!   on a given timesteps )

! ----------------------------------------------------------------------
implicit none

! Local variables
integer(i_def) :: n(n_diag_group)
integer(i_def) :: i_dgroup  
integer(i_def) :: i
logical :: l_diag_requested
character(len=ukca_maxlen_diagname) :: diagname

! Array to temporarily hold names and status flags of requested diagnostics
character(len=ukca_maxlen_diagname) :: tmp_diagnames_fullht_real(n_req_max(i_dgroup_3d))
integer(i_um) :: tmp_diag_status_3d(n_req_max(i_dgroup_3d))

! Error handling variables
integer(i_um) :: errcode
character(len=ukca_maxlen_message) :: ukca_errmsg

! End of header

errcode = 0_i_um
tmp_diagnames_fullht_real(:) = ''
tmp_diag_status_3d(:) = 0_i_um

! Counter for active requests
n_ukca_diags_3d = 0_i_def

do i = 1, n_req_max(i_dgroup_3d)
  n_ukca_diags_3d = n_ukca_diags_3d + 1
  tmp_diagnames_fullht_real(i) = diagnames_map(i, 2)
  tmp_diag_status_3d(i) = ukca_diag_status_requested
end do

! Populate the allocatable arrays with the requested diagnostics
allocate(diagnames_fullht_real(n_ukca_diags_3d))
allocate(idiag_status_3d(n_ukca_diags_3d))
diagnames_fullht_real = tmp_diagnames_fullht_real(1:n_ukca_diags_3d)
idiag_status_3d = tmp_diag_status_3d(1:n_ukca_diags_3d)

! Pass on active diagnostic information to UKCA via API - if any requested
if ( n_ukca_diags_3d > 0_i_def ) then
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

end subroutine ukca_diag_setup
! ----------------------------------------------------------------------

end module ukca_diag_setup_mod
