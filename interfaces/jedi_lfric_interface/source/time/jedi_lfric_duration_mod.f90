!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief   A module providing JEDI/OOPS like duration functionality
!>
!> @details This module emulates the duration object used in the
!!          LFRic-JEDI interface, see:
!!          https://github.com/JCSDA-internal/oops/blob/develop/src/oops/util/Duration.h
module jedi_lfric_duration_mod

  use constants_mod,               only : i_timestep, str_def, l_def
  use log_mod,                     only : log_event,         &
                                          log_scratch_space, &
                                          LOG_LEVEL_DEBUG,   &
                                          LOG_LEVEL_ERROR,   &
                                          LOG_LEVEL_WARNING

  implicit none
  private

  integer(i_timestep), parameter :: seconds_in_day    = 86400
  integer(i_timestep), parameter :: seconds_in_hour   = 3600
  integer(i_timestep), parameter :: seconds_in_minute = 60

  !> @brief Stores the duration in seconds between two jedi_datetime objects
  type, public :: jedi_duration_type
    private

    integer(i_timestep) :: seconds  !< duration between two times in seconds

  contains

    procedure, private :: init_empty
    procedure, private :: init_iso_string
    procedure, private :: init_iso_string_err
    procedure, private :: init_seconds
    generic,   public  :: init => init_empty,      &
                                  init_iso_string, &
                                  init_seconds

    procedure, public  :: get_duration

    procedure, public  :: to_string

    procedure, private, pass(self) :: add
    procedure, private, pass(self) :: subtract
    procedure, private, pass(self) :: multiply
    procedure, private, pass(self) :: divide

    procedure, private, pass(self) :: add_int
    procedure, private, pass(self) :: subtract_int
    procedure, private, pass(self) :: multiply_int
    procedure, private, pass(self) :: divide_int

    generic,   public              :: operator(+) => add, add_int
    generic,   public              :: operator(-) => subtract, subtract_int
    generic,   public              :: operator(*) => multiply, multiply_int
    generic,   public              :: operator(/) => divide, divide_int

    procedure, private, pass(self) :: eq
    procedure, private, pass(self) :: ne
    procedure, private, pass(self) :: gt
    procedure, private, pass(self) :: lt
    procedure, private, pass(self) :: ge
    procedure, private, pass(self) :: le

    procedure, private, pass(self) :: eq_int
    procedure, private, pass(self) :: ne_int
    procedure, private, pass(self) :: gt_int
    procedure, private, pass(self) :: lt_int
    procedure, private, pass(self) :: ge_int
    procedure, private, pass(self) :: le_int

    generic,   public              :: operator(==) => eq, eq_int
    generic,   public              :: operator(/=) => ne, ne_int
    generic,   public              :: operator(>)  => gt, gt_int
    generic,   public              :: operator(<)  => lt, lt_int
    generic,   public              :: operator(>=) => ge, ge_int
    generic,   public              :: operator(<=) => le, le_int

  end type jedi_duration_type

contains

  !> @brief Initialise a duration instance with 0 seconds
  subroutine init_empty( self )

    implicit none

    class( jedi_duration_type ), intent(inout) :: self

    call log_event( 'Initialsing JEDI duration = 0 seconds', LOG_LEVEL_DEBUG )

    self%seconds = 0

  end subroutine init_empty

  !> @brief Initialise a duration instance with an ISO string
  !!
  !> @param [in] iso_duration An ISO 8601 duration string of
  !!                          the form [-]P[nD][T][nH][nM][nS],
  !!                          one or both of D and T must be present
  subroutine init_iso_string( self, iso_duration )

    implicit none

    class( jedi_duration_type ), intent(inout) :: self
    character(*),                intent(in)    :: iso_duration

    integer(i_timestep)     :: err
    character(str_def) :: message

    integer(i_timestep)     :: is_negative
    character(str_def) :: absolute_duration

    integer(i_timestep)     :: last_index

    integer(i_timestep)     :: D_index
    integer(i_timestep)     :: days

    integer(i_timestep)     :: T_index

    integer(i_timestep)     :: H_index
    integer(i_timestep)     :: hours
    integer(i_timestep)     :: M_index
    integer(i_timestep)     :: minutes
    integer(i_timestep)     :: S_index
    integer(i_timestep)     :: seconds

    write ( log_scratch_space, '(A)' )                     &
            'Initialising JEDI duration from iso_string: ' &
            // iso_duration
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG)

    days    = 0
    hours   = 0
    minutes = 0
    seconds = 0

    self%seconds = 0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Check if first character is '-' for negative durations
    is_negative = scan( iso_duration, '-', kind=i_timestep )
    if ( is_negative > 1 ) then       ! found '-' sign not at string start
      message = 'Failed negative check'
      call self%init_iso_string_err( message )
    else if ( is_negative == 1 ) then
      absolute_duration = iso_duration(2:)  ! remove first char ( '-' sign )
    else
      absolute_duration = iso_duration      ! number is positive
    end if

    last_index = 1

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Scan for the D, T, H, M, S chars

    D_index = scan( absolute_duration, 'D', kind=i_timestep )
    T_index = scan( absolute_duration, 'T', kind=i_timestep )

    ! raise error if T and D aren't present
    if ( (T_index == 0) .and. (D_index == 0) ) then
      message = 'No D or T present'
      call self%init_iso_string_err( message )
    end if

    H_index = scan( absolute_duration, 'H', kind=i_timestep )
    M_index = scan( absolute_duration, 'M', kind=i_timestep )
    S_index = scan( absolute_duration, 'S', kind=i_timestep )

    ! raise error if no D, H, M, or S is present
    if ( (D_index + H_index + M_index + S_index) == 0 ) then
      message = 'No D, H, M, or S present'
      call self%init_iso_string_err( message )
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! read days if present
    if ( D_index /= 0 ) then
      read ( absolute_duration( last_index+1:D_index-1 ), *, iostat=err ) days
      if ( err /= 0 ) then
        message = 'Failed to read D'
        call self%init_iso_string_err( message )
      end if
      last_index = D_index
    end if

    ! T is present after D or P so increase last_index
    last_index = last_index + 1

    ! read hours if present
    if ( H_index /= 0 ) then
      read ( absolute_duration( last_index+1:H_index-1 ), *, iostat=err ) hours
      if ( err /= 0 ) then
        message = 'Failed to read H'
        call self%init_iso_string_err( message )
      end if
      last_index = H_index
    end if

    ! read minutes if present
    if ( M_index /= 0 ) then
      read ( absolute_duration( last_index+1:M_index-1 ), *, iostat=err ) minutes
      if ( err /= 0 ) then
        message = 'Failed to read M'
        call self%init_iso_string_err( message )
      end if
      last_index = M_index
    end if

    ! read seconds if present
    if ( S_index /= 0 ) then
      read ( absolute_duration( last_index+1:S_index-1 ), *, iostat=err ) seconds
      if ( err /= 0 ) then
        message = 'Failed to read S'
        call self%init_iso_string_err( message )
      end if
    end if

    self%seconds = seconds                        &
                   + minutes * seconds_in_minute  &
                   + hours * seconds_in_hour      &
                   + days * seconds_in_day

    if ( is_negative == 1 ) then
      self%seconds = - self%seconds
    end if

  end subroutine init_iso_string

  !> @brief Logs an error if initialising a duration instance fails
  !!
  !> @param [in] extra_message An extra message string to append
  !!                           to the main error message
  subroutine init_iso_string_err( self, extra_message )

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    character(*),                intent(in) :: extra_message

    write ( log_scratch_space, '(A)' ) &
              'Creating JEDI duration FAIL: Expecting string &
              & of form [-]P[nD][T][nH][nM][nS] - '          &
              // extra_message
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )

  end subroutine init_iso_string_err

  !> @brief Initialise a duration from a time in seconds
  !!
  !> @param [in] seconds Time duration in seconds
  subroutine init_seconds( self, seconds )

    implicit none

    class( jedi_duration_type ), intent(inout) :: self
    integer(i_timestep),         intent(in)    :: seconds

    write ( log_scratch_space, * ) seconds
    call log_event( 'Initialising JEDI duration with: ' &
                    //trim(log_scratch_space)//' seconds', LOG_LEVEL_DEBUG )

    self%seconds = seconds

  end subroutine init_seconds

  !> @brief Returns the duration in seconds
  !!
  !> @param [out] duration Seconds
  subroutine get_duration( self, duration )

    implicit none

    class( jedi_duration_type ), intent(in)  :: self
    integer(i_timestep),         intent(out) :: duration

    duration = self%seconds

  end subroutine get_duration

  !> @brief   Returns the duration as an iso string
  !!
  !> @details Returns in the form [-]P[nD][T][nH][nM][nS],
  !!          output for 0 seconds is PT0S
  !!
  !> @param [inout] iso_duration The string to return
  subroutine to_string( self, iso_duration )

    implicit none

    class( jedi_duration_type ), intent(in)    :: self
    character(str_def),          intent(inout) :: iso_duration

    character(str_def) :: temp_char

    integer(i_timestep) :: remainder

    integer(i_timestep) :: days
    integer(i_timestep) :: hours
    integer(i_timestep) :: minutes

    remainder = self%seconds

    if ( remainder < 0 ) then
      iso_duration = '-P'
      remainder = abs(remainder)
    else
      iso_duration = 'P'
    end if

    days = remainder / seconds_in_day
    if ( days /= 0 ) then
      write ( temp_char, * ) days
      iso_duration = trim(iso_duration) // trim(adjustl(temp_char)) // 'D'
      remainder = remainder - days * seconds_in_day
    end if

    if ( (remainder /= 0) .or. (days == 0) ) then
      iso_duration = trim(iso_duration) // 'T'
    end if

    hours = remainder / seconds_in_hour
    if ( hours /= 0 ) then
      write ( temp_char, * ) hours
      iso_duration = trim(iso_duration) // trim(adjustl(temp_char)) // 'H'
      remainder = remainder - hours * seconds_in_hour
    end if

    minutes = remainder / seconds_in_minute
    if ( minutes /= 0 ) then
      write ( temp_char, * ) minutes
      iso_duration = trim(iso_duration) // trim(adjustl(temp_char)) // 'M'
      remainder = remainder - minutes * seconds_in_minute
    end if

    if ( (remainder /= 0) .or. (self%seconds == 0) ) then
      write ( temp_char, * ) remainder
      iso_duration = trim(iso_duration) // trim(adjustl(temp_char)) // 'S'
    end if

  end subroutine to_string

  !> @brief Overload + operator for duration + duration
  !!
  !> @param [in] duration     duration instance to add
  !> @result     new_duration New duration instance
  function add( self, duration ) result(new_duration)

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    type( jedi_duration_type ),  intent(in) :: duration

    type( jedi_duration_type ) :: new_duration

    integer(i_timestep) :: seconds

    seconds = self%seconds + duration%seconds

    call new_duration%init( seconds )

  end function add

  !> @brief Overload + operator for duration + integer
  !!
  !> @param [in] duration     duration in seconds to add
  !> @result     new_duration New duration instance
  function add_int( self, duration ) result(new_duration)

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    integer(i_timestep),         intent(in) :: duration

    type( jedi_duration_type ) :: new_duration

    integer(i_timestep) :: seconds

    seconds = self%seconds + duration

    call new_duration%init( seconds )

  end function add_int

  !> @brief Overload - operator for duration - duration
  !!
  !> @param [in] duration     duration instance to subtract
  !> @result     new_duration New duration instance
  function subtract( self, duration ) result(new_duration)

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    type( jedi_duration_type ),  intent(in) :: duration

    type( jedi_duration_type ) :: new_duration

    integer(i_timestep) :: seconds

    seconds = self%seconds - duration%seconds

    call new_duration%init( seconds )

  end function subtract

  !> @brief Overload - operator for duration - integer
  !!
  !> @param [in] duration     duration in seconds to subtract
  !> @result     new_duration New duration instance
  function subtract_int( self, duration ) result(new_duration)

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    integer(i_timestep),         intent(in) :: duration

    type( jedi_duration_type ) :: new_duration

    integer(i_timestep) :: seconds

    seconds = self%seconds - duration

    call new_duration%init( seconds )

  end function subtract_int

  !> @brief Overload * operator for duration * duration
  !!
  !> @param [in] duration     duration instance to multiply by
  !> @result     new_duration New duration instance
  function multiply( self, duration ) result(new_duration)

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    type( jedi_duration_type ),  intent(in) :: duration

    type( jedi_duration_type ) :: new_duration

    integer(i_timestep) :: seconds

    seconds = self%seconds * duration%seconds

    call new_duration%init( seconds )

  end function multiply

  !> @brief Overload * operator for integer * duration
  !!
  !> @param [in] duration     duration in seconds to multiply by
  !> @result     new_duration New duration instance
  function multiply_int( self, duration ) result(new_duration)

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    integer(i_timestep),         intent(in) :: duration

    type( jedi_duration_type ) :: new_duration

    integer(i_timestep) :: seconds

    seconds = duration * self%seconds

    call new_duration%init( seconds )

  end function multiply_int

  !> @brief Overload / operator for duration / duration
  !!
  !> @param [in] duration     duration instance to divide by
  !> @result     new_duration New duration instance
  function divide( self, duration ) result(new_duration)

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    type( jedi_duration_type ),  intent(in) :: duration

    type( jedi_duration_type ) :: new_duration

    integer(i_timestep) :: seconds

    if ( duration%seconds == 0 ) then
      write ( log_scratch_space, '(A)' ) &
        'Cannot divide duration by 0'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    else if ( mod(self%seconds, duration%seconds) /= 0 ) then
      write ( log_scratch_space, '(A)' ) &
        'Durations not evenly divisible'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    if ( (self%seconds > 0) .and. (duration%seconds < 0) ) then
      write ( log_scratch_space, '(A)' ) &
        'Dividing + duration by - duration'
      call log_event( log_scratch_space, LOG_LEVEL_WARNING )
    end if

    seconds = self%seconds / duration%seconds

    call new_duration%init( seconds )

  end function divide

  !> @brief Overload / operator for duration / integer
  !!
  !> @param [in] duration     duration in seconds to divide by
  !> @result     new_duration New duration instance
  function divide_int( self, duration ) result(new_duration)

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    integer(i_timestep),         intent(in) :: duration

    type( jedi_duration_type ) :: new_duration

    integer(i_timestep) :: seconds

    if ( duration == 0 ) then
      write ( log_scratch_space, '(A)' ) &
        'Cannot divide duration by 0'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    else if ( mod(self%seconds, duration) /= 0 ) then
      write ( log_scratch_space, '(A)' ) &
        'Durations not evenly divisible'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    if ( (self%seconds > 0) .and. (duration < 0) ) then
      write ( log_scratch_space, '(A)' ) &
        'Dividing + duration by - duration'
      call log_event( log_scratch_space, LOG_LEVEL_WARNING )
    end if

    seconds = self%seconds / duration

    call new_duration%init( seconds )

  end function divide_int

  !> @brief Overload == comparator
  !!
  !> @param [in] duration  duration to compare to
  function eq( self, duration )

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    type( jedi_duration_type ),  intent(in) :: duration

    logical(l_def) :: eq

    eq = ( self%seconds == duration%seconds )

  end function eq

  !> @brief Overload == comparator
  !!
  !> @param [in] seconds  seconds to compare to
  function eq_int( self, seconds )

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    integer(i_timestep),         intent(in) :: seconds

    logical(l_def) :: eq_int

    eq_int = ( self%seconds == seconds )

  end function eq_int

  !> @brief Overload /= comparator
  !!
  !> @param [in] duration  duration to compare to
  function ne( self, duration )

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    type( jedi_duration_type ),  intent(in) :: duration

    logical(l_def) :: ne

    ne = ( self%seconds /= duration%seconds )

  end function ne

  !> @brief Overload /= comparator
  !!
  !> @param [in] seconds  seconds to compare to
  function ne_int( self, seconds )

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    integer(i_timestep),         intent(in) :: seconds

    logical(l_def) :: ne_int

    ne_int = ( self%seconds /= seconds )

  end function ne_int

  !> @brief Overload > comparator
  !!
  !> @param [in] duration  duration to compare to
  function gt( self, duration )

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    type( jedi_duration_type ),  intent(in) :: duration

    logical(l_def) :: gt

    gt = ( self%seconds > duration%seconds )

  end function gt

  !> @brief Overload > comparator
  !!
  !> @param [in] seconds  seconds to compare to
  function gt_int( self, seconds )

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    integer(i_timestep),         intent(in) :: seconds

    logical(l_def) :: gt_int

    gt_int = ( self%seconds > seconds )

  end function gt_int

  !> @brief Overload < comparator
  !!
  !> @param [in] duration  duration to compare to
  function lt( self, duration )

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    type( jedi_duration_type ),  intent(in) :: duration

    logical(l_def) :: lt

    lt = ( self%seconds < duration%seconds )

  end function lt

  !> @brief Overload < comparator
  !!
  !> @param [in] seconds  seconds to compare to
  function lt_int( self, seconds )

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    integer(i_timestep),         intent(in) :: seconds

    logical(l_def) :: lt_int

    lt_int = ( self%seconds < seconds )

  end function lt_int

  !> @brief Overload >= comparator
  !!
  !> @param [in] duration  duration to compare to
  function ge( self, duration )

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    type( jedi_duration_type ),  intent(in) :: duration

    logical(l_def) :: ge

    ge = ( self%seconds >= duration%seconds )

  end function ge

  !> @brief Overload >= comparator
  !!
  !> @param [in] seconds  seconds to compare to
  function ge_int( self, seconds )

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    integer(i_timestep),         intent(in) :: seconds

    logical(l_def) :: ge_int

    ge_int = ( self%seconds >= seconds )

  end function ge_int

  !> @brief Overload <= comparator
  !!
  !> @param [in] duration  duration to compare to
  function le( self, duration )

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    type( jedi_duration_type ),  intent(in) :: duration

    logical(l_def) :: le

    le = ( self%seconds <= duration%seconds )

  end function le

  !> @brief Overload <= comparator
  !!
  !> @param [in] seconds  seconds to compare to
  function le_int( self, seconds )

    implicit none

    class( jedi_duration_type ), intent(in) :: self
    integer(i_timestep),         intent(in) :: seconds

    logical(l_def) :: le_int

    le_int = ( self%seconds <= seconds )

  end function le_int

end module jedi_lfric_duration_mod
