!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief   A module providing JEDI/OOPS like date time functionality
!>
!> @details This module emulates the date time object used in the
!!          LFRic-JEDI interface
module jedi_lfric_datetime_mod

  use constants_mod,                     only : i_timestep, r_def, str_def, l_def
  use jedi_lfric_datetime_functions_mod, only : YYYYMMDD_to_JDN,   &
                                                JDN_to_YYYYMMDD,   &
                                                hhmmss_to_seconds, &
                                                seconds_to_hhmmss, &
                                                is_valid_datetime
  use jedi_lfric_duration_mod,           only : jedi_duration_type
  use log_mod,                           only : log_event,         &
                                                log_scratch_space, &
                                                LOG_LEVEL_INFO,    &
                                                LOG_LEVEL_ERROR

  implicit none
  private

  integer(i_timestep), parameter  :: seconds_in_day = 86400

  !> @brief This type stores the JEDI date and time
  type, public :: jedi_datetime_type
    private

    integer(i_timestep) :: date    !< Julian day number
    integer(i_timestep) :: time    !< seconds since start of day
    logical(l_def)      :: initialised = .false.

  contains

    procedure, private :: init_iso_string
    procedure, private :: init_string
    procedure, private :: init_YYMMDD_hhmmss
    procedure, private :: init_YYYYMMDDhhmmss
    generic,   public  :: init => init_iso_string,     &
                                  init_YYMMDD_hhmmss,  &
                                  init_YYYYMMDDhhmmss

    procedure, public  :: get_date
    procedure, public  :: get_time

    procedure, public  :: add_seconds

    procedure, public  :: to_string
    procedure, public  :: print

    ! new_datetime = datetime
    procedure, private, pass(self) :: copy
    generic,   public              :: assignment(=) => copy

    ! new_datetime = datetime + duration
    procedure, private, pass(self) :: add_duration
    generic,   public              :: operator(+) => add_duration

    ! duration = datetime1 - datetime2
    procedure, private, pass(self) :: seconds_between
    generic,   public              :: operator(-) => seconds_between

    ! logical = datetime1 > datetime2
    procedure, private, pass(self) :: is_greater_than
    generic,   public              :: operator(>) => is_greater_than

    ! logical = datetime1 < datetime2
    procedure, private, pass(self) :: is_less_than
    generic,   public              :: operator(<) => is_less_than


  end type jedi_datetime_type

contains

  !> @brief   Initialise a datetime using an iso string (UTC)
  !> @details Initialise a jedi datetime using an iso string
  !!          of the form 2023-04-05T11:41:38Z. The T may be
  !!          a space and the timezone, Z is ignored
  !!
  !> @param [in] iso_datetime ISO datetime string
  subroutine init_iso_string( self, iso_datetime )

    implicit none

    class( jedi_datetime_type ), intent(inout) :: self
    character(*),                intent(in)    :: iso_datetime

    write ( log_scratch_space, '(2(A))' ) &
      'Creating JEDI datetime using iso datetime string: ', iso_datetime
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    call self%init_string( iso_datetime )

  end subroutine init_iso_string

  !> @brief   Initialise a datetime using an iso string (UTC)
  !!
  !> @param [in] iso_datetime ISO datetime string
  subroutine init_string( self, iso_datetime )

    implicit none

    class( jedi_datetime_type ), intent(inout) :: self
    character(*),                intent(in)    :: iso_datetime

    integer(i_timestep) :: year
    integer(i_timestep) :: month
    integer(i_timestep) :: day

    integer(i_timestep) :: hour
    integer(i_timestep) :: minute
    integer(i_timestep) :: second

    integer(i_timestep) :: err

    read ( iso_datetime(1:4), '(I4)', iostat=err ) year
    read ( iso_datetime(6:7), '(I2)', iostat=err ) month
    read ( iso_datetime(9:10), '(I2)', iostat=err  ) day

    read ( iso_datetime(12:13), '(I2)', iostat=err  ) hour
    read ( iso_datetime(15:16), '(I2)', iostat=err  ) minute
    read ( iso_datetime(18:19), '(I2)', iostat=err  ) second

    if ( err /= 0 ) then
      write ( log_scratch_space, '(A)' ) &
              'Creating JEDI datetime FAIL: Failed to read iso_datetime string'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    call self%init_YYYYMMDDhhmmss( year, month, day, hour, minute, second )

  end subroutine init_string

  !> @brief   Initialise a datetime with two integers,
  !!          one for the date and one for time
  !> @details Initialise a datetime with a YYYYMMDD integer
  !!          for the year, month, and day eg. 20230405,
  !!          and a hhmmss integer for the hours, minutes,
  !!          and seconds eg. 130210.
  !!
  !> @param [in] YYYYMMDD Integer representing the year, month, day
  !> @param [in] hhmmss   Integer representing the hour, minute, second
  subroutine init_YYMMDD_hhmmss( self, YYYYMMDD, hhmmss )

    implicit none

    class( jedi_datetime_type ), intent(inout) :: self
    integer(i_timestep), intent(in) :: YYYYMMDD !< year, month, and day eg 20230405
    integer(i_timestep), intent(in) :: hhmmss   !< hour, minute, second eg 130210

    integer(i_timestep) :: year
    integer(i_timestep) :: month
    integer(i_timestep) :: day

    integer(i_timestep) :: hour
    integer(i_timestep) :: minute
    integer(i_timestep) :: second

    integer(i_timestep) :: temp_int

    year     = YYYYMMDD / 10000        ! remove lower 4 digits, yields YYYY
    temp_int = mod( YYYYMMDD, 10000 )  ! keep lower 4 digits, yields MMDD
    month    = temp_int / 100          ! remove lower 2 digits, yields MM
    day      = mod( temp_int, 100 )    ! keep lower 2 digits, yields DD

    hour     = int(hhmmss, i_timestep) / 10000        ! remove lower 4 digits, yields hh
    temp_int = mod( int(hhmmss, i_timestep), 10000 )  ! keep lower 4 digits, yields mmss
    minute   = temp_int / 100                    ! remove lower 2 digits, yields mm
    second   = mod( temp_int, 100 )              ! keep lower 2 digits, yeilds ss

    call self%init_YYYYMMDDhhmmss( year, month, day, hour, minute, second )

  end subroutine init_YYMMDD_hhmmss

  !> @brief      Initialise a datetime with six integers for the year, month,
  !!             day, hour, minute, and second
  !!
  !> @param [in] year   Year   - YYYY
  !> @param [in] month  Month  - MM
  !> @param [in] day    Day    - DD
  !> @param [in] hour   Hour   - hh
  !> @param [in] minute Minute - mm
  !> @param [in] second Second - ss
  subroutine init_YYYYMMDDhhmmss( self, year, month, day, hour, minute, second )

    implicit none

    class( jedi_datetime_type ), intent(inout) :: self
    integer(i_timestep), intent(in) :: year    !< in YYYY format, eg 2020
    integer(i_timestep), intent(in) :: month   !< in MM format, 01 (Jan) through 12 (Dec)
    integer(i_timestep), intent(in) :: day     !< in DD format, 01 through 28, 29, 30, or 31

    integer(i_timestep), intent(in) :: hour    !< in hh format, 00 through 23
    integer(i_timestep), intent(in) :: minute  !< in mm format, 00 through 59
    integer(i_timestep), intent(in) :: second  !< in ss format, 00 through 59

    call log_event( 'Initialising JEDI datetime', LOG_LEVEL_INFO )

    ! Also checks the date and time are valid
    call YYYYMMDD_to_JDN( year, month, day, self%date )
    call hhmmss_to_seconds( hour, minute, second, self%time )

    self%initialised = .true.

    write ( log_scratch_space, '(A, I8, A, I5)' )                   &
           'Initialising JEDI datetime SUCCESS: JDN = ', self%date, &
           ' Time (s) = ', self%time
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine init_YYYYMMDDhhmmss

  !> @brief Gets the date from the datetime instance
  !!
  !> @param [out] date The Julian Day Number
  subroutine get_date( self, date )

    implicit none

    class( jedi_datetime_type ), intent(in)  :: self
    integer(i_timestep),         intent(out) :: date

    date = self%date

  end subroutine get_date

  !> @brief Gets the time from the datetime instance
  !!
  !> @param [out] time Seconds since the start of the day
  subroutine get_time( self, time )

    implicit none

    class( jedi_datetime_type ), intent(in)  :: self
    integer(i_timestep),         intent(out) :: time

    time = self%time

  end subroutine get_time

  !> @brief Adds a time in seconds to the datetime instance
  !!
  !> @param [in] seconds Time in seconds to add, can be negative
  subroutine add_seconds( self, seconds )

    implicit none

    class( jedi_datetime_type ), intent(inout) :: self
    integer(i_timestep),         intent(in)    :: seconds

    integer(i_timestep) :: days
    integer(i_timestep) :: new_time

    new_time = seconds
    days = new_time / seconds_in_day

    self%date = self%date + days
    new_time = mod( new_time, seconds_in_day )

    self%time = self%time + new_time

    if ( self%time < 0 ) then
      self%date = self%date - 1
      self%time = self%time + seconds_in_day
    else if ( self%time >= seconds_in_day ) then
      self%date = self%date + 1
      self%time = self%time - seconds_in_day
    end if

  end subroutine

  !> @brief Returns true if the datetime is greater than the passed datetime
  !!
  !> @param [in] datetime The jedi_datetime instance to compare with
  function is_greater_than( self, datetime )

    implicit none

    class(jedi_datetime_type), intent(in) :: self
    type(jedi_datetime_type),  intent(in) :: datetime

    logical(l_def) :: is_greater_than

    type( jedi_duration_type ) :: duration

    is_greater_than = .false.

    duration = datetime - self

    if ( duration < 0 ) is_greater_than = .true.

  end function is_greater_than

  !> @brief Returns true if the datetime is less than the passed datetime
  !!
  !> @param [in] datetime The jedi_datetime instance to compare with
  function is_less_than( self, datetime )

    implicit none

    class(jedi_datetime_type), intent(in) :: self
    type(jedi_datetime_type),  intent(in) :: datetime

    logical(l_def) :: is_less_than

    type( jedi_duration_type ) :: duration

    is_less_than = .false.

    duration = datetime - self

    if ( duration > 0 ) is_less_than = .true.

  end function is_less_than

  !> @brief Returns the datetime as an iso string (UTC)
  !!
  !> @param [inout] iso_datetime The string to return
  subroutine to_string( self, iso_datetime )

    implicit none

    class( jedi_datetime_type ),    intent(in) :: self
    character(str_def),          intent(inout) :: iso_datetime

    integer(i_timestep) :: year
    integer(i_timestep) :: month
    integer(i_timestep) :: day

    integer(i_timestep) :: hour
    integer(i_timestep) :: minute
    integer(i_timestep) :: second

    character(len=4) :: temp_str_4
    character(len=2) :: temp_str_2
    character(len=1) :: dash, space, colon

    call JDN_to_YYYYMMDD( self%date, year, month, day )
    call seconds_to_hhmmss( self%time, hour, minute, second )

    dash  = '-'
    space = ' '
    colon = ':'

    write ( temp_str_4, '(I4)' ) year
    iso_datetime = temp_str_4 // dash
    write ( temp_str_2, '(I2.2)' ) month
    iso_datetime = trim(iso_datetime) // temp_str_2 // dash
    write ( temp_str_2, '(I2.2)' ) day
    iso_datetime = trim(iso_datetime) // temp_str_2

    write ( temp_str_2, '(I2.2)' ) hour
    iso_datetime = trim(iso_datetime) // space // temp_str_2 // colon
    write ( temp_str_2, '(I2.2)' ) minute
    iso_datetime = trim(iso_datetime) // temp_str_2 // colon
    write ( temp_str_2, '(I2.2)' ) second
    iso_datetime = trim(iso_datetime) // temp_str_2

  end subroutine to_string

  !> @brief Writes the curent dateime via log_event
  subroutine print( self )

    implicit none

    class( jedi_datetime_type ), intent(in) :: self

    ! Local
    character(str_def) :: iso_datetime

    call self%to_string(iso_datetime)

    write ( log_scratch_space, '(2(A))' ) &
      'JEDI iso datetime: ', iso_datetime
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine print

  !> @brief Overload assignment
  !!
  !> @param [in] from datetime instance to copy from
  subroutine copy( self, from )

    implicit none

    class( jedi_datetime_type ), intent(inout) :: self
    class( jedi_datetime_type ), intent(in)    :: from

    if ( .not. from%initialised ) then
      write ( log_scratch_space, '(A)' ) &
        'Cannot initialise this datetime from an uninitialised datetime'
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end if

    self%date = from%date
    self%time = from%time
    self%initialised = .true.

  end subroutine copy

  !> @brief Overload + operator for datetime + duration
  !!
  !> @param [in] duration     duration instance to add to datetime
  !> @result     new_datetime New datetime instance
  function add_duration( self, duration ) result(new_datetime)

    implicit none

    class( jedi_datetime_type ), intent(in) :: self
    type( jedi_duration_type ),  intent(in) :: duration

    type( jedi_datetime_type ) :: new_datetime

    integer(i_timestep) :: seconds

    new_datetime%date = self%date
    new_datetime%time = self%time

    call duration%get_duration( seconds )
    call new_datetime%add_seconds( seconds )

    new_datetime%initialised = .true.

  end function add_duration

  !> @brief Calculates the number of seconds between two datetimes
  !!
  !> @param [in] datetime  The jedi_datetime instance to subtract
  !> @result     duration  Seconds between the two datetimes
  function seconds_between( self, datetime ) result(duration)

    implicit none

    class( jedi_datetime_type ), intent(in)    :: self
    type( jedi_datetime_type ),  intent(in)    :: datetime

    type( jedi_duration_type ) :: duration

    integer(i_timestep)        :: diff_date
    integer(i_timestep)        :: diff_time
    integer(i_timestep)        :: diff_seconds

    diff_date = self%date - datetime%date
    diff_time = self%time - datetime%time

    diff_seconds = diff_date * seconds_in_day + diff_time
    call duration%init( diff_seconds )

  end function seconds_between

end module jedi_lfric_datetime_mod
