!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Adjoint of mixed precision field copy routines for real32 and real64 fields.
module adj_copy_field_mod

  use, intrinsic :: iso_fortran_env, only: real32, real64

  ! Missing invoke_adj_copy_field_64_32 since this is non
  ! linear. Source must not have greater precision than destination.
  use adj_sci_psykal_builtin_light_mod, only : &
                          invoke_adj_copy_field_32_64, &
                          invoke_adj_copy_field_32_32, &
                          invoke_adj_copy_field_64_64

  use field_real32_mod,      only : field_real32_type
  use field_real64_mod,      only : field_real64_type
  use log_mod,               only : log_event, LOG_LEVEL_ERROR

  implicit none

  ! We keep adj_copy_field_64_32 in the interface to throw
  ! an error in case someone tries to use it.
  interface adj_copy_field
     module procedure &
        adj_copy_field_32_32, adj_copy_field_64_64, &
        adj_copy_field_32_64, adj_copy_field_64_32
  end interface adj_copy_field

  contains
  !--------------------------------------------------------------
  !> Adjoint of field copy routines for all exisiting combinations
  !> of field real types.
  subroutine adj_copy_field_32_64(fsrce_32, fdest_64)
    implicit none
    type(field_real32_type), intent(inout)  :: fsrce_32
    type(field_real64_type), intent(inout)  :: fdest_64
    call invoke_adj_copy_field_32_64(fsrce_32, fdest_64)
  end subroutine adj_copy_field_32_64

  subroutine adj_copy_field_64_32(fsrce_64, fdest_32)
    implicit none
    type(field_real64_type), intent(in)     :: fsrce_64
    type(field_real32_type), intent(in)     :: fdest_32
    call log_event( "Cannot perform non-linear adjoint copy of field, source is more precise than destination.", &
                    LOG_LEVEL_ERROR )
  end subroutine adj_copy_field_64_32

  subroutine adj_copy_field_32_32(fsrce_32, fdest_32)
    implicit none
    type(field_real32_type), intent(inout)  :: fsrce_32
    type(field_real32_type), intent(inout)  :: fdest_32
    call invoke_adj_copy_field_32_32(fsrce_32, fdest_32)
  end subroutine adj_copy_field_32_32

  subroutine adj_copy_field_64_64(fsrce_64, fdest_64)
    implicit none
    type(field_real64_type), intent(inout)  :: fsrce_64
    type(field_real64_type), intent(inout)  :: fdest_64
    call invoke_adj_copy_field_64_64(fsrce_64, fdest_64)
  end subroutine adj_copy_field_64_64

end module adj_copy_field_mod
