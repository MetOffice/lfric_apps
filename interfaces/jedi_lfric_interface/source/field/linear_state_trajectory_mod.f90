!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing the linear state trajectory object.
!>
!> @details This module provides the linear state trajectory object for use
!>          with the linear model. Its purpose is to store an array of field
!>          collections containing the linear state fields. The fields are
!>          created and populated by the non-linear model. The linear model
!>          uses the data stored in the trajectory as part of the linear model
!>          time step. Because the time-resolution of the linear states can be
!>          less than the linear time-step. Linear time interpolation is
!>          applied when retrieving the linear state for a given time-step.
!>
module linear_state_trajectory_mod

  use jedi_lfric_datetime_mod,         only : jedi_datetime_type
  use jedi_lfric_duration_mod,         only : jedi_duration_type
  use log_mod,                         only : log_event,          &
                                              log_scratch_space,  &
                                              LOG_LEVEL_ERROR
  use constants_mod,                   only : i_def, str_def, r_def, &
                                              imdi, rmdi
  use field_collection_mod,            only : field_collection_type
  use trajectory_field_utils_mod,      only : copy_fields,        &
                                              copy_create_fields, &
                                              interpolate_fields

  implicit none

  private

type, public :: linear_state_trajectory_type
  private

  !> An array of linear state field collections that store the 4D trajectory
  type( field_collection_type ), allocatable :: linear_state_fields(:)

  !> An array of datetime objects the correspond to the linear state array
  type( jedi_datetime_type ), allocatable    :: linear_state_time(:)

  !> The index of the next item to be stored in the trajectory
  integer( kind=i_def)                       :: next_item_index

  !> The number of items stored in the trajectory
  integer( kind=i_def)                       :: number_of_items

  !> The length of the time-step
  type( jedi_duration_type )                 :: time_step

contains

  !> Initialiser for the linear_state_trajectory_trajectory.
  procedure          :: initialise

  !> Add the next linear state and associate with a datetime
  procedure          :: add_linear_state

  !> Get the linear state for a given datetime
  procedure          :: get_linear_state

  !> Finalizer
  final              :: linear_state_trajectory_destructor

end type linear_state_trajectory_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> @brief    Initialiser for linear_state_trajectory_type
!>
!> @param [in] forecast_length  The duration of the forecast to store
!> @param [in] time_step        The duration of a time step
subroutine initialise( self, forecast_length, time_step )

  implicit none

  class( linear_state_trajectory_type ), intent(inout) :: self
  type( jedi_duration_type ),               intent(in) :: forecast_length
  type( jedi_duration_type ),               intent(in) :: time_step

  ! Local
  integer :: forecast_length_secs
  integer :: time_step_secs

  ! Get the number of slots
  call forecast_length%get_duration(forecast_length_secs)
  call time_step%get_duration(time_step_secs)
  self%number_of_items = forecast_length_secs/time_step_secs + 1

  ! Setup the trajectory
  allocate( self%linear_state_fields(self%number_of_items) )
  allocate( self%linear_state_time(self%number_of_items) )

  self%next_item_index = 1_i_def
  self%time_step = time_step

end subroutine initialise

!> @brief    Add the next next linear state to the trajectory
!>
!> @param [in]  next_datetime      The datetime to associated  with the next
!>                                 linear state field collection
!> @param [out] next_linear_state  The next linear state field collection
subroutine add_linear_state( self, next_datetime, next_linear_state )

  implicit none

  class( linear_state_trajectory_type ), intent(inout) :: self
  type( jedi_datetime_type ),               intent(in) :: next_datetime
  type( field_collection_type ),            intent(in) :: next_linear_state

  ! Local
  type(jedi_duration_type) :: time_difference
  integer                  :: i_item
  integer                  :: next_item_index

  ! Create and add new item
  next_item_index = self%next_item_index

  ! If no slots available
  if (next_item_index > self%number_of_items) then
    write(log_scratch_space, '(2A)' ) &
        "There are no slots available, ", &
        "the object was not initialised with enough model_data slots."
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! Check duplicates are not added
  do i_item=1,next_item_index-1

    time_difference = next_datetime - self%linear_state_time(i_item)

    ! The next_datetime should be current_datetime + time_step: abort if it is not
    if (i_item==next_item_index-1) then
      if ( time_difference/=self%time_step ) then
        log_scratch_space = 'The input datetime is not one time step greater than the previous datetime.'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      endif
    endif

    ! If slot has been filled then check the date is not already stored
    if (time_difference==0) then
      log_scratch_space = 'The input datetime has already been used.'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    endif

  enddo

  ! Copy create the linear_state_fields stored internally from the collection
  ! passed in.
  call copy_create_fields( self%linear_state_fields(next_item_index), &
                                      next_linear_state, &
                                      prefix_name="ls_")
  self%linear_state_time(next_item_index) = next_datetime
  self%next_item_index = self%next_item_index + 1_i_def

end subroutine add_linear_state

!> @brief    Get the linear state field collection for a given datetime
!>
!> @param [in]  req_time          The datetime to associated with the linear
!>                                state field collection
!> @param [out] req_linear_state  The linear state field collection to be
!>                                updated at the requested datetime
subroutine get_linear_state( self, req_time, req_linear_state )

  implicit none

  class( linear_state_trajectory_type ), intent(inout) :: self
  type( jedi_datetime_type ),               intent(in) :: req_time
  type( field_collection_type ),         intent(inout) :: req_linear_state

  ! Local
  type(jedi_duration_type) :: time_difference
  integer( kind=i_def )    :: time_difference_secs
  integer( kind=i_def )    :: time_step_secs
  integer( kind=i_def )    :: i_item
  integer( kind=i_def )    :: time_index(2)
  real( kind=r_def )       :: time_weight(2)
  character( len=str_def ) :: req_time_char
  character( len=str_def ) :: traj_time_char

  ! Check the requested time is not before the start of the trajectory
  time_difference = req_time - self%linear_state_time(1)
  if (time_difference < 0) then
    call self%linear_state_time(1)%to_string(traj_time_char)
    call req_time%to_string(req_time_char)

    write(log_scratch_space, '(4A)' )                      &
          'The requested datetime: ', trim(req_time_char), &
          ', is before the first trajectory datetime: ',   &
          trim(traj_time_char)
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  endif

  time_weight = rmdi
  time_index = imdi
  do i_item=1,self%number_of_items
    time_difference = req_time - self%linear_state_time(i_item)

    ! Get the time weight and index
    if (time_difference == 0) then
      time_weight = [1.0_r_def , 0.0_r_def]
      time_index = [i_item, i_item]
      exit
    elseif (time_difference < 0) then
      call time_difference%get_duration(time_difference_secs)
      call self%time_step%get_duration(time_step_secs)
      time_weight(1)= real(-1*time_difference_secs)/real(time_step_secs)
      time_weight(2)= 1.0_r_def - time_weight(1)

      time_index = [i_item-1, i_item]
      exit
    endif
  enddo

  ! Trajectory index equals imdi if not found
  if (time_index(1)==imdi) then
    call self%linear_state_time(self%number_of_items)%to_string(traj_time_char)
    call req_time%to_string(req_time_char)

    write(log_scratch_space, '(4A)' )                      &
          'The requested datetime: ', trim(req_time_char), &
          ', is after the last trajectory datetime: ',     &
          trim(traj_time_char)
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  endif

  ! Interpolate
  if (time_index(1)==time_index(2)) then
    ! No interpolation required, just copy the correct collection
    call copy_fields( req_linear_state, &
                      self%linear_state_fields(time_index(1)) )

  else
    ! Do the interpolation using weights and indexes and put result into the output collection
    call interpolate_fields( req_linear_state,                        &
                             self%linear_state_fields(time_index(1)), &
                             self%linear_state_fields(time_index(2)), &
                             time_weight(1),                          &
                             time_weight(2) )
  endif

end subroutine get_linear_state

!> @brief    The linear_state_trajectory_type finalizer
!>
subroutine linear_state_trajectory_destructor( self )

  implicit none

  type( linear_state_trajectory_type ), intent(inout) :: self

  if ( allocated(self%linear_state_fields) ) &
                                  deallocate(self%linear_state_fields )
  if ( allocated(self%linear_state_time) ) deallocate(self%linear_state_time )

end subroutine linear_state_trajectory_destructor

end module linear_state_trajectory_mod
