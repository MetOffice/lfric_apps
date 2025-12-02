!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Module for printing mean and root mean square of fields.

module print_meanrms_field_mod
  use constants_mod,                   only : i_def, str_def
  use field_mod,                       only : field_type
  use log_mod,                         only : log_event,         &
                                              log_level,         &
                                              log_scratch_space
  use norm_alg_mod,                    only : mean_alg,         &
                                              root_mean_square_alg

  implicit none

  private
  public :: print_meanrms_field

  contains

  !> @brief Prints the mean and root mean square of a field.
  !> @param[in] field   Field to be printed
  !> @param[in] level   Level of logging. If the configured
  !>                    log_level is less than or equal to
  !>                    level, output will be shown.
  !> @param[in] fname   Optional name of field.
  subroutine print_meanrms_field( field, level, fname )

    implicit none

    type( field_type ),    intent( in )           :: field
    integer( kind=i_def ), intent( in )           :: level
    character (*),         intent( in ), optional :: fname

    character ( str_def ) :: field_name

    if ( log_level() <= level ) then

      if ( present( fname ) ) then
        field_name = fname
      else
        field_name = field % get_name()
      end if

      write( log_scratch_space, '( A, A, A, 2E16.8 )' ) &
      "mean/rms ", trim( field_name ), " = ",           &
      mean_alg( field ), root_mean_square_alg( field )
      call log_event( log_scratch_space, level )

    end if

  end subroutine print_meanrms_field

end module print_meanrms_field_mod
