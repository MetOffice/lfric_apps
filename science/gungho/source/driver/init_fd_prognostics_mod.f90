!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>@brief Initialisation of finite difference prognostic fields
!>@details At present this means reading from a (UM2LFRic) start dump
module init_fd_prognostics_mod

  use log_mod,                         only: log_event,         &
                                             log_scratch_space, &
                                             LOG_LEVEL_INFO,    &
                                             LOG_LEVEL_TRACE
  use lfric_xios_read_mod,             only: read_state

  ! Derived Types
  use field_mod,                       only: field_type
  use field_collection_mod,            only: field_collection_type

  implicit none

  private
  public :: init_fd_prognostics_dump

contains

  !> @details Initialise FD prognostic fields
  !>          from UM2LFric dump (checkpoint format)
  !> @param[in,out] fd_field_collection The collection of FD fields
  subroutine init_fd_prognostics_dump( fd_field_collection )

    implicit none

    ! FD Prognostic fields collection

    type(field_collection_type), intent(inout) :: fd_field_collection

    call read_state(fd_field_collection,prefix='read_')


    call log_event( "Physics: Initialised FD prognostic fields from UM2LFRic dump", LOG_LEVEL_INFO )


  end subroutine init_fd_prognostics_dump


end module init_fd_prognostics_mod
