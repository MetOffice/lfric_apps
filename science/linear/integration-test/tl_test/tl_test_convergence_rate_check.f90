!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>@brief   Test the convergence rate (Taylor remainder convergence).
module tl_test_convergence_rate_check

  use constants_mod,                  only: r_def, str_def
  use log_mod,                        only: log_event,         &
                                            log_scratch_space, &
                                            LOG_LEVEL_INFO

  implicit none

  private
  public convergence_rate_check

  contains

  !> @brief   Calculate and test the convergence rate.
  !> @details Calculate the convergence rate based on the norms
  !>          at two different iterations, and compare with the
  !>          expected value. Print out either PASS or FAIL, which
  !>          will then be used by the top level integration test.
  !> @param[in] norm_diff       Norm at first iteration
  !> @param[in] norm_diff_prev  Norm at second iteration
  !> @param[in] label           Name of the code being tested
  !> @param[in] tol             Tolerance value
  subroutine convergence_rate_check( norm_diff, norm_diff_prev, label, tol )

    implicit none

    real(r_def),           intent(in) :: norm_diff, norm_diff_prev
    character(str_def),    intent(in) :: label
    real(r_def), optional, intent(in) :: tol
    real(r_def) :: tolerance
    real(r_def) :: conv_rate
    character(len=4) :: pass_str

    if ( present(tol) ) then
      tolerance = tol
    else
      tolerance = 1.0E-8_r_def
    end if

    conv_rate =  norm_diff_prev/ norm_diff



    write( log_scratch_space, '(A)' ) &
           "TL Test: " // trim(label)
    call log_event( log_scratch_space, LOG_LEVEL_INFO )
    write( log_scratch_space, '(A, E16.8)') &
           "Convergence rate: ", conv_rate
    call log_event( log_scratch_space, LOG_LEVEL_INFO )


    if ( abs(conv_rate - 4.0_r_def ) < tolerance ) then
      pass_str = "PASS"
    else
      pass_str = "FAIL"
    end if

    write(log_scratch_space,'("   test",A32," : ",A4)') trim(label), pass_str
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine convergence_rate_check

end module tl_test_convergence_rate_check
