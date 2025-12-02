! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module init_diag_array_mod

implicit none

contains

!----------------------------------------------------------------
! Subroutine to check that a requested diagnostic has the
! appropriate pointer assigned to some memory of the appropriate
! size and shape, and initialise it to zero if so
!----------------------------------------------------------------
subroutine init_diag_array( diag )

use comorph_constants_mod, only: n_updraft_types, n_dndraft_types,             &
                                 n_conv_layers_diag,                           &
                                 nx_full, ny_full, k_bot_conv, k_top_conv,     &
                                 newline
use diag_type_mod, only: diag_type

use init_zero_mod, only: init_zero_2d, init_zero_3d
use raise_error_mod, only: raise_fatal

implicit none

! Structure containing diagnostic meta-data and pointer to
! output array
type(diag_type), intent(in out) :: diag

! Number of convection types (for updrafts or downdrafts)
integer :: n_conv_types

! Lower and upper bounds of array pointers
integer :: lb(5), ub(5)

! Flags for possible errors
logical :: error_1, error_2, error_3

! Loop counters
integer :: i_type, i_layr

character(len=*), parameter :: routinename = "INIT_DIAG_ARRAY"


! Initialise error flags to false
error_1 = .false.
error_2 = .false.
error_3 = .false.

! Diagnostic requested in x,y domain
if ( diag % request % x_y ) then
  ! Check the appropriate pointer is assigned
  if ( .not. associated( diag % field_2d ) ) then
    error_1 = .true.
  else
    ! Check that the array has sufficient extent
    lb(1:2) = lbound( diag % field_2d )
    ub(1:2) = ubound( diag % field_2d )
    if ( .not. ( lb(1) <= 1 .and. ub(1) >= nx_full                             &
           .and. lb(2) <= 1 .and. ub(2) >= ny_full ) ) then
      error_2 = .true.
    else
      ! If we passed the above 2 checks, initialise the
      ! diagnostic to zero.
      call init_zero_2d( lb(1:2), ub(1:2), diag % field_2d )
    end if
  end if
end if

! Diagnostic requested in x,y domain for all
! updraft or downdraft types
if ( diag % request % x_y_typ ) then
  ! Check the appropriate pointer is assigned
  if ( .not. associated( diag % field_3d ) ) then
    error_1 = .true.
  else
    ! Number of convection types depends on whether updraft
    ! or downdraft
    if ( diag % diag_name(1:7) == "updraft" ) then
      n_conv_types = n_updraft_types
    else if ( diag % diag_name(1:7) == "dndraft" ) then
      n_conv_types = n_dndraft_types
    else
      ! Problem if not an updraft or downdraft diagnostic
      n_conv_types = 0
      error_3 = .true.
    end if
    ! Check that the array has sufficient extent
    lb(1:3) = lbound( diag % field_3d )
    ub(1:3) = ubound( diag % field_3d )
    if ( .not. ( lb(1) <= 1 .and. ub(1) >= nx_full                             &
           .and. lb(2) <= 1 .and. ub(2) >= ny_full                             &
           .and. lb(3) <= 1 .and. ub(3) >= n_conv_types ) ) then
      error_2 = .true.
    else
      ! If we passed the above 2 checks, initialise the
      ! diagnostic to zero.
      do i_type = 1, n_conv_types
        call init_zero_2d( lb(1:2), ub(1:2), diag % field_3d(:,:,i_type) )
      end do
    end if
  end if
end if

! Diagnostic requested in x,y domain for all
! updraft or downdraft types and layers
if ( diag % request % x_y_lay_typ ) then
  ! Check the appropriate pointer is assigned
  if ( .not. associated( diag % field_4d ) ) then
    error_1 = .true.
  else
    ! Number of convection types depends on whether updraft
    ! or downdraft
    if ( diag % diag_name(1:7) == "updraft" ) then
      n_conv_types = n_updraft_types
    else if ( diag % diag_name(1:7) == "dndraft" ) then
      n_conv_types = n_dndraft_types
    else
      ! Problem if not an updraft or downdraft diagnostic
      n_conv_types = 0
      error_3 = .true.
    end if
    ! Check that the array has sufficient extent
    lb(1:4) = lbound( diag % field_4d )
    ub(1:4) = ubound( diag % field_4d )
    if ( .not. ( lb(1) <= 1 .and. ub(1) >= nx_full                             &
           .and. lb(2) <= 1 .and. ub(2) >= ny_full                             &
           .and. lb(3) <= 1 .and. ub(3) >= n_conv_layers_diag                  &
           .and. lb(4) <= 1 .and. ub(4) >= n_conv_types ) ) then
      error_2 = .true.
    else
      ! If we passed the above 2 checks, initialise the
      ! diagnostic to zero.
      do i_type = 1, n_conv_types
        do i_layr = 1, n_conv_layers_diag
          call init_zero_2d( lb(1:2), ub(1:2),                                 &
                             diag % field_4d(:,:,i_layr,i_type) )
        end do
      end do
    end if
  end if
end if

! Diagnostic requested in x,y,z domain
if ( diag % request % x_y_z ) then
  ! Check the appropriate pointer is assigned
  if ( .not. associated( diag % field_3d ) ) then
    error_1 = .true.
  else
    ! Check that the array has sufficient extent
    lb(1:3) = lbound( diag % field_3d )
    ub(1:3) = ubound( diag % field_3d )
    if ( .not. ( lb(1) <= 1 .and. ub(1) >= nx_full                             &
           .and. lb(2) <= 1 .and. ub(2) >= ny_full                             &
           .and. lb(3) <= k_bot_conv .and. ub(3) >= k_top_conv ) ) then
      error_2 = .true.
    else
      ! If we passed the above 2 checks, initialise the
      ! diagnostic to zero.
      call init_zero_3d( lb(1:3), ub(1:3), diag % field_3d )
    end if
  end if
end if

! Diagnostic requested in x,y,z domain for all
! updraft or downdraft types
if ( diag % request % x_y_z_typ ) then
  ! Check the appropriate pointer is assigned
  if ( .not. associated( diag % field_4d ) ) then
    error_1 = .true.
  else
    ! Number of convection types depends on whether updraft
    ! or downdraft
    if ( diag % diag_name(1:7) == "updraft" ) then
      n_conv_types = n_updraft_types
    else if ( diag % diag_name(1:7) == "dndraft" ) then
      n_conv_types = n_dndraft_types
    else
      ! Problem if not an updraft or downdraft diagnostic
      n_conv_types = 0
      error_3 = .true.
    end if
    ! Check that the array has sufficient extent
    lb(1:4) = lbound( diag % field_4d )
    ub(1:4) = ubound( diag % field_4d )
    if ( .not. ( lb(1) <= 1 .and. ub(1) >= nx_full                             &
           .and. lb(2) <= 1 .and. ub(2) >= ny_full                             &
           .and. lb(3) <= k_bot_conv .and. ub(3) >= k_top_conv                 &
           .and. lb(4) <= 1 .and. ub(4) >= n_conv_types ) ) then
      error_2 = .true.
    else
      ! If we passed the above 2 checks, initialise the
      ! diagnostic to zero.
      do i_type = 1, n_conv_types
        call init_zero_3d( lb(1:3), ub(1:3), diag%field_4d(:,:,:,i_type) )
      end do
    end if
  end if
end if

! Diagnostic requested in x,y,z domain for all
! updraft or downdraft types and layers
if ( diag % request % x_y_z_lay_typ ) then
  ! Check the appropriate pointer is assigned
  if ( .not. associated( diag % field_5d ) ) then
    error_1 = .true.
  else
    ! Number of convection types depends on whether updraft
    ! or downdraft
    if ( diag % diag_name(1:7) == "updraft" ) then
      n_conv_types = n_updraft_types
    else if ( diag % diag_name(1:7) == "dndraft" ) then
      n_conv_types = n_dndraft_types
    else
      ! Problem if not an updraft or downdraft diagnostic
      n_conv_types = 0
      error_3 = .true.
    end if
    ! Check that the array has sufficient extent
    lb(1:5) = lbound( diag % field_5d )
    ub(1:5) = ubound( diag % field_5d )
    if ( .not. ( lb(1) <= 1 .and. ub(1) >= nx_full                             &
           .and. lb(2) <= 1 .and. ub(2) >= ny_full                             &
           .and. lb(3) <= k_bot_conv .and. ub(3) >= k_top_conv                 &
           .and. lb(4) <= 1 .and. ub(4) >= n_conv_layers_diag                  &
           .and. lb(5) <= 1 .and. ub(5) >= n_conv_types ) ) then
      error_2 = .true.
    else
      ! If we passed the above 2 checks, initialise the
      ! diagnostic to zero.
      do i_type = 1, n_conv_types
        do i_layr = 1, n_conv_layers_diag
          call init_zero_3d( lb(1:3), ub(1:3),                                 &
                             diag%field_5d(:,:,:,i_layr,i_type) )
        end do
      end do
    end if
  end if
end if

! Print error message if problem detected
if ( error_1 ) then
  call raise_fatal( routinename,                                               &
         "Diagnostic " // trim(adjustl( diag % diag_name )) // " " //          &
         "has been requested but not assigned to "   //newline//               &
         "any memory." )
end if
if ( error_2 ) then
  call raise_fatal( routinename,                                               &
         "Diagnostic " // trim(adjustl( diag % diag_name )) // " " //          &
         "has been assigned to an array which has "  //newline//               &
         "insufficient extent / the wrong shape for the "     //               &
         "field to be output in it." )
end if
if ( error_3 ) then
  call raise_fatal( routinename,                                               &
         "Diagnostic " // trim(adjustl( diag % diag_name )) // " " //          &
         "has been requested with a separate field " //newline//               &
         "for each convection type, but is not available in " //               &
         "that form." )
end if

return
end subroutine init_diag_array

end module init_diag_array_mod
