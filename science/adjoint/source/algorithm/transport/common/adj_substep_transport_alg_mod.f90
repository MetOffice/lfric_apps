!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Adjoint routine of the setup of the substepping of the transport

module adj_substep_transport_alg_mod

  use constants_mod,        only : i_def, r_tran
  use log_mod,              only : log_event,         &
                                   log_scratch_space, &
                                   LOG_LEVEL_ERROR,   &
                                   LOG_LEVEL_INFO
  use r_tran_field_mod,     only : r_tran_field_type
  use transport_config_mod, only : substep_transport,          &
                                   substep_transport_adaptive, &
                                   substep_transport_off,      &
                                   substep_transport_two,      &
                                   substep_transport_four

  implicit none

  private
  public :: adj_substep_transport_alg

  contains

  !=============================================================================
  !> @brief   Set up the number of substeps when substepping transport
  !> @details This code currently does not have a test since the active fields
  !!          are not modified. It solely exists as a way of preventing developers
  !!          from using the incorrect configuration for the adjoint.
  !> @param[in,out] substeps   Number of substeps
  !> @param[in]     wind       The advecting wind
  !> @param[in]     dt         Transport time step
  subroutine adj_substep_transport_alg(substeps, &
                                       wind,     &
                                       dt)

    implicit none

    integer(kind=i_def),     intent(inout) :: substeps
    type(r_tran_field_type), intent(in)    :: wind
    real(kind=r_tran),       intent(in)    :: dt

    select case ( substep_transport )

    case ( substep_transport_adaptive )
      call log_event('substep_transport_adaptive has non-linear logic in forward routine, cannot be adjointed', &
                      LOG_LEVEL_ERROR)

    case ( substep_transport_off )
      ! Set number of substeps to 1
      substeps = 1_i_def

    case ( substep_transport_two )
      ! Set number of substeps to 2 for debugging
      substeps = 2_i_def

    case ( substep_transport_four )
      ! Set number of substeps to 4 for debugging
      substeps = 4_i_def

    case default
      call log_event('Unrecognised substep_transport option', &
                      LOG_LEVEL_ERROR)

    end select

    if (substeps > 1) then
      write( log_scratch_space, '(A,I4)' ) &
      'Transport: number of gungho transport substeps = ', substeps
      call log_event(log_scratch_space, LOG_LEVEL_INFO)
    end if

  end subroutine adj_substep_transport_alg

end module adj_substep_transport_alg_mod
