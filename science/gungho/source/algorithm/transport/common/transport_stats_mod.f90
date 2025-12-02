!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Contains a routine for writing out statistics for transport tests
!> @details This file contains a routine for calculating and then logging
!!          statistics used in transport tests for evaluating the quality of
!!          transport schemes, such as field minima and maxima and L2 errors.
module transport_stats_mod

  use constants_mod,                     only: r_def, str_def, EPS
  use field_mod,                         only: field_type
  use log_mod,                           only: log_event,           &
                                               log_scratch_space,   &
                                               LOG_LEVEL_INFO
  use sci_field_minmax_alg_mod,          only: get_field_minmax

  implicit none

  private

  public :: write_transport_stats

contains


  !> @brief Routine to calculate and write out various transport statistics
  !> @param[in] field       The transported field
  !> @param[in] true_field  The true field to compare against
  !> @param[in] field_name  String to print out next to statistics
  subroutine write_transport_stats( field, true_field, field_name )

    use norm_alg_mod,     only: l2_norm_alg,              &
                                rel_l2_error_alg,         &
                                volume_normalisation_alg, &
                                dispersion_error_alg,     &
                                dissipation_error_alg

    implicit none

    type(field_type),       intent(in) :: field, true_field
    character(len=str_def), intent(in) :: field_name
    real(kind=r_def)                   :: min_field, max_field, volume
    real(kind=r_def)                   :: min_field0, max_field0
    real(kind=r_def)                   :: l2_field, l2_field0
    real(kind=r_def)                   :: l2_error, diss, disp

    ! Compute statistics
    ! field min and max
    call get_field_minmax( field, min_field, max_field )

    ! true_field min and max
    call get_field_minmax( true_field, min_field0, max_field0 )

    volume = volume_normalisation_alg(field)
    l2_field = l2_norm_alg(field)
    l2_field0 = l2_norm_alg(true_field)
    l2_error = rel_l2_error_alg(field, true_field)
    diss = dissipation_error_alg(field, true_field)
    disp = dispersion_error_alg(field, true_field)

    ! Normalise disp and diss errors to compare them with relative l2 error
    if ((l2_field0 / sqrt(volume)) > EPS) then
      diss = diss / (l2_field0 / sqrt(volume))
      disp = disp / (l2_field0 / sqrt(volume))
    end if

    ! Write out statistics to log
    write( log_scratch_space, '(A, E15.8)' ) &
      "Transport stats: Min-initial " // trim(field_name) // " =", min_field0
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(A, E15.8)' ) &
      "Transport stats: Max-initial " // trim(field_name) // " =", max_field0
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(A, E15.8)' ) &
      "Transport stats: Min-final " // trim(field_name) // " =", min_field
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(A, E15.8)' ) &
      "Transport stats: Max-final " // trim(field_name) // " =", max_field
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(A, E15.8)' ) &
      "Transport stats: L2-initial " // trim(field_name) // " =", l2_field0
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(A, E15.8)' ) &
      "Transport stats: L2-final " // trim(field_name) // " =", l2_field
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(A, E15.8)' ) &
      "Transport stats: Rel-L2-error " // trim(field_name) // " =", l2_error
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(A, E15.8)' ) &
      "Transport stats: Dissipation " // trim(field_name) // " =", diss
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(A, E15.8)' ) &
      "Transport stats: Dispersion " // trim(field_name) // " =", disp
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine write_transport_stats

end module transport_stats_mod
