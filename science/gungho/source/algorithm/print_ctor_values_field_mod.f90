!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Module for printing values used in constructor of a field.
!> @details Key issue is not knowing the proper constructor values of fields
!>          in live environments, so this tool can be used to determine the
!>          constructor values so e.g. adjoint tests can be developed more easily.
module print_ctor_values_field_mod

  use constants_mod,                   only : i_def, str_def
  use field_mod,                       only : field_type
  use function_space_mod,              only : function_space_type
  use log_mod,                         only : log_event,         &
                                              log_level,         &
                                              log_scratch_space
  use fs_continuity_mod,               only : name_from_functionspace

  implicit none

  private
  public :: print_ctor_values_field

  contains

  !=============================================================================
  !> @brief Prints the values that were used to construct a field.
  !> @param[in] field   Field to be printed
  !> @param[in] level   Level of logging. If the configured
  !>                    log_level is less than or equal to
  !>                    level, output will be shown.
  !> @param[in] fname   Optional name of field.
  subroutine print_ctor_values_field( field, level, fname )

    implicit none

    ! Arguments
    type(field_type),         intent(in) :: field
    integer(kind=i_def),      intent(in) :: level
    character(*), intent(in), optional   :: fname

    ! Internal variables
    type(function_space_type), pointer   :: fs_ptr
    character(str_def)                   :: field_name

    nullify(fs_ptr)

    if( log_level() <= level ) then
      if ( present( fname ) ) then
        field_name = fname
      else
        field_name = field%get_name()
      end if

      write(log_scratch_space, "(a, a)") "Printing constructor values for: ", trim( field_name )
      call log_event( log_scratch_space, level )

      fs_ptr => field%get_function_space()

      ! Function space
      write(log_scratch_space,"(1x, a, a)") "Function space name = ", trim( name_from_functionspace( fs_ptr%which() ) )
      call log_event( log_scratch_space, level )

      ! Element order
      write(log_scratch_space, "(1x, a, i3)") "Horizontal element order = ", fs_ptr%get_element_order_h()
      call log_event( log_scratch_space, level )
      write(log_scratch_space, "(1x, a, i3)") "Vertical element order = ", fs_ptr%get_element_order_v()
      call log_event( log_scratch_space, level )

      ! Multidata
      write(log_scratch_space, "(1x, a, i3)") "ndata = ", fs_ptr%get_ndata()
      call log_event( log_scratch_space, level )
      write(log_scratch_space, "(1x, a, l1)") "ndata_first = ", fs_ptr%is_ndata_first()
      call log_event( log_scratch_space, level )

      ! Halo depth
      write(log_scratch_space, "(1x, a, i3)") "halo_depth = ", field%get_field_halo_depth()
      call log_event( log_scratch_space, level )

    end if

    end subroutine print_ctor_values_field

end module print_ctor_values_field_mod

