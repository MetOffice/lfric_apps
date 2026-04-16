!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Module containing the lfric_apps version string components and
!> character return function (for printing & labelling).

module lfric_apps_version_mod

  use constants_mod,                 only: i_def, l_def, str_def

  implicit none

  private
  public lfric_apps_version_char

  integer(i_def), public  :: lfric_apps_major_version = 3
  integer(i_def), public  :: lfric_apps_minor_version = 1
  integer(i_def), public  :: lfric_apps_patch_version = 0
  logical(l_def), public  :: lfric_apps_release_version = .false.
  character(2), parameter :: prefix = 'vn'
  character(4), parameter :: dev_suffix = '_dev'

contains

  !> Return a character representation of the current version
  !> of lfric_apps.
  !> This will only include the patch version if it is greater than zero
  !> and will only skip the '_dev' suffix if the release logical is .true.
  !> e.g. 3.1, 3.1_dev, 3.1.1_dev, 3.2, 3.2_dev
  function lfric_apps_version_char() result(apps_version_char)
    character(str_def) :: apps_version_char
    if (lfric_apps_patch_version > 0) then
      if (lfric_apps_release_version) then
        write (unit=apps_version_char,fmt='(A,I0,A,I0,A,I0)') prefix, &
               lfric_apps_major_version, '.', lfric_apps_minor_version, '.', &
               lfric_apps_patch_version
      else
        write (unit=apps_version_char,fmt='(A,I0,A,I0,A,I0,A)') prefix, &
               lfric_apps_major_version, '.', lfric_apps_minor_version, '.', &
               lfric_apps_patch_version, dev_suffix
      end if
    else
      if (lfric_apps_release_version) then
        write (unit=apps_version_char,fmt='(A,I0,A,I0)') prefix, &
               lfric_apps_major_version, '.', lfric_apps_minor_version
      else
        write (unit=apps_version_char,fmt='(A,I0,A,I0,A)') prefix, &
               lfric_apps_major_version, '.', lfric_apps_minor_version, &
               dev_suffix
      end if
    end if

  end function lfric_apps_version_char

end module lfric_apps_version_mod
