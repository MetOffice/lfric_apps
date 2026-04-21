!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Module containing the lfric_apps version string components and
!> character return function (for printing & labelling).

module lfric_apps_version_mod

  use constants_mod,          only: i_def, l_def, str_def
  use lfric_core_version_mod, only: lfric_core_version_char

  implicit none

  private
  public lfric_apps_version_char

  integer(i_def), public, parameter   :: lfric_apps_major_version = 3
  integer(i_def), public, parameter   :: lfric_apps_minor_version = 1
  integer(i_def), public, parameter   :: lfric_apps_patch_version = 0
  logical(l_def), public, parameter   :: lfric_apps_release_version = .false.
  character(2), parameter             :: prefix = 'vn'
  character(4), parameter             :: dev_suffix = '_dev'

contains

  !> Return a character representation of the current version
  !> of lfric_apps.
  !> This will only include the patch version if it is greater than zero
  !> and will only skip the '_dev' suffix if the release logical is .true.
  !> e.g. 3.1, 3.1_dev, 3.1.1_dev, 3.2, 3.2_dev
  function lfric_apps_version_char(input_major_version, input_minor_version, &
                                   input_patch_version, input_release_version) &
                                   result(apps_version_char)

    character(str_def)                   :: apps_version_char
    integer(i_def), intent(in), optional :: input_major_version, &
                                            input_minor_version, &
                                            input_patch_version
    logical(l_def), intent(in), optional :: input_release_version
    integer(i_def)                       :: lfric_major_version, &
                                            lfric_minor_version, &
                                            lfric_patch_version
    logical(l_def)                       :: lfric_release_version

    if (present(input_major_version)) then
      lfric_major_version = input_major_version
    else
      lfric_major_version = lfric_apps_major_version
    end if
    if (present(input_minor_version)) then
      lfric_minor_version = input_minor_version
    else
      lfric_minor_version = lfric_apps_minor_version
    end if
    if (present(input_patch_version)) then
      lfric_patch_version = input_patch_version
    else
      lfric_patch_version = lfric_apps_patch_version
    end if
    if (present(input_release_version)) then
      lfric_release_version = input_release_version
    else
      lfric_release_version = lfric_apps_release_version
    end if

    apps_version_char = lfric_core_version_char(lfric_major_version, &
                                                lfric_minor_version, &
                                                lfric_patch_version, &
                                                lfric_release_version)

  end function lfric_apps_version_char

end module lfric_apps_version_mod
