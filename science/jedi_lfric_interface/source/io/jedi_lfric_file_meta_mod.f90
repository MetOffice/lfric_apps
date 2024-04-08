!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing a IO meta-data type for JEDI-LFRIC.
!>
!> @details Defines an object to store the meta-data associated with
!>          setting up IO for JEDI-LFRIC required to setup the file_list
!>          used by the lfric-xios context initialiser.
!>
module jedi_lfric_file_meta_mod

  use constants_mod,  only : str_def, l_def, i_def
  use file_mod,       only : FILE_MODE_READ, FILE_MODE_WRITE
  use log_mod,        only : log_event, log_scratch_space, &
                             LOG_LEVEL_ERROR

  implicit none

  private

type, public :: jedi_lfric_file_meta_type

    !> The name of the file including full path if it is not located in the current directory
    character(str_def)  :: file_name

    !> A unique string that will be used to access the file within a given context
    character(str_def)  :: xios_id

    !> The mode to perform IO with defined as an integer value
    integer(i_def)      :: io_mode

    !> The frequency in timesteps that the file is read-from/written-to
    integer(i_def)      :: freq

    !> The XIOS ID for the group that defines the fields to be used in the IO
    character(str_def)  :: field_group_id

contains

  !> Initialiser
  procedure, public :: initialise

  !> jedi_lfric_file_meta_type finalizer
  final :: jedi_lfric_file_meta_destructor

end type jedi_lfric_file_meta_type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @brief    Initialiser for jedi_lfric_file_meta_type
!>
!> @param[in] file_name      The name of the file including full path if
!>                           not located in pwd
!> @param[in] xios_id        A unique string that will be used to access the
!>                           file within a given context
!> @param[in] io_mode_str    A string either "read" or "write" that defines the
!>                           mode to perform IO
!> @param[in] freq           The frequency in timesteps that the file is
!>                           read-from/written-to
!> @param[in] field_group_id The XIOS ID for the group that defines the fields
!>                           to be used in the IO
subroutine initialise( self, file_name, xios_id, io_mode_str, &
                       freq, field_group_id )

  implicit none

  class( jedi_lfric_file_meta_type ), intent(inout) :: self

  character(str_def), intent(in) :: file_name
  character(str_def), intent(in) :: xios_id
  character(str_def), intent(in) :: io_mode_str
  integer(i_def),     intent(in) :: freq
  character(str_def), intent(in) :: field_group_id

  self%file_name = file_name
  self%xios_id = xios_id
  ! Set the io_mode based on the supplied value
  select case (convert_to_lower(io_mode_str))
    case ("read")
      self%io_mode = FILE_MODE_READ
    case ("write")
      self%io_mode = FILE_MODE_WRITE
    case default
       write( log_scratch_space, '(3A)' ) &
             'Received: io_mode_str = ', trim( io_mode_str ), &
             '. Expected "read" or "write", check the input configuration.'
       call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end select
  self%freq = freq
  self%field_group_id = field_group_id

end subroutine initialise

!> @brief    jedi_lfric_file_meta_type finalizer
!>
subroutine jedi_lfric_file_meta_destructor(self)

  implicit none

  type( jedi_lfric_file_meta_type ), intent(inout) :: self

end subroutine jedi_lfric_file_meta_destructor

!=============================================================================!
! Local methods
!=============================================================================!

!> @brief Changes a string to lower case
!> @param[in] str Input string to convert
!> @result string Lower case string
pure function convert_to_lower (str) Result (string)

  implicit none

  character(*), intent(in) :: str
  character(len(str))      :: string

  integer(i_def) :: ic
  integer(i_def) :: i

  character(26), parameter :: cap = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
  character(26), parameter :: low = 'abcdefghijklmnopqrstuvwxyz'

  string = str
  do i = 1, len_trim(str)
      ic = index(cap, str(i:i))
      if (ic > 0) string(i:i) = low(ic:ic)
  end do

end function convert_to_lower

end module jedi_lfric_file_meta_mod
