!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>@brief   Test the convergence rate (Taylor remainder convergence).
module tl_test_convergence_rate_check

  use constants_mod,                  only: r_def, str_def, i_def
  use log_mod,                        only: log_event,         &
                                            log_scratch_space, &
                                            LOG_LEVEL_INFO

  implicit none

  private
  public convergence_pass_string,      &
         array_convergence_rate_check, &
         convergence_rate_check

contains

  !> @brief   Test the convergence rate, with application of square root.
  !> @details If the convergence rate is not close to 4, within the
  !!          specified tolerance, set the pass string to FAIL.
  !> @param[in] name         Variable or test name
  !> @param[in] norm         Current norm
  !> @param[in] norm_prev    Previous norm
  !> @param[in,out] pass_str Pass string (either PASS or FAIL)
  !> @param[in] tol          Tolerance
  subroutine convergence_pass_string( name, norm, norm_prev, pass_str, tol )
    implicit none

    real(r_def),           intent(in)    :: norm, norm_prev
    character(len=4),      intent(inout) :: pass_str
    real(r_def),           intent(in)    :: tol
    character(str_def),    intent(in)    :: name

    real(r_def) :: conv_rate

    ! Let the error between the nonlinear (N) difference and the linear (L) be
    ! E(g) = || N(x +  g dx) - N(x) - g Ldx || = O(g^2)
    ! where g is a scalar, x and dx are vectors, O is the order
    ! and || . || = (x^T x)^1/2 is the L2 norm.
    !
    ! Then the ratio
    ! E(2g) / E(g) = O(4 g^2) / O(g^2) = 4
    ! i.e. we need to check whether the convergence rate is close to 4.

    ! The norms have not had a square root applied yet.
    conv_rate =  sqrt(norm_prev / norm)

    if ( abs( conv_rate - 4.0_r_def ) >= tol ) then
      pass_str = "FAIL"
    else
      pass_str = "PASS"
    end if

    write( log_scratch_space, '(A, A, A, E16.8, A, E16.8)') &
     name, pass_str, " Convergence rate: ", conv_rate, "  Tolerance: ", tol
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine convergence_pass_string

  !> @brief   Test the convergence rate for individual variables and the sum.
  !> @details Calculate the convergence rate based on the norms
  !>          at two different iterations, and compare with the
  !>          expected value. Print out either PASS or FAIL, which
  !>          will then be used by the top level integration test.
  !> @param[in] array_norm       Norm at second iteration
  !> @param[in] array_norm_prev  Norm at first iteration
  !> @param[in] array_names      Name of the variable being tested
  !> @param[in] n_variables      Array lengths (number of variables)
  !> @param[in] label            Test name
  !> @param[in] tol              Tolerance value
  subroutine array_convergence_rate_check( array_norm, array_norm_prev, array_names, n_variables, label, tol, indiv_tol)
    integer(i_def), intent(in) :: n_variables
    real(r_def),    intent(in) :: array_norm(n_variables)
    real(r_def),    intent(in) :: array_norm_prev(n_variables)
    character(str_def),    intent(in) :: array_names(n_variables)
    character(str_def),    intent(in) :: label
    real(r_def), optional, intent(in) :: tol
    real(r_def), optional, intent(in) :: indiv_tol
    real(r_def) :: tolerance, individual_tolerance
    character(len=4) :: pass_str_arr(n_variables)
    character(len=4) :: sum_pass_str, pass_str
    character(str_def), parameter :: sum_name = "sum"
    integer(i_def) :: i
    real(r_def) :: sum, sum_prev

    call log_event( "Checking convergence rate", LOG_LEVEL_INFO )

    if ( present(tol) ) then
      tolerance = tol
    else
      tolerance = 1.0E-8_r_def
    end if

    if ( present(indiv_tol) ) then
      individual_tolerance = indiv_tol
    else
      individual_tolerance = 1.0E-8_r_def
    end if

    sum = 0.0_r_def
    sum_prev = 0.0_r_def
    do i= 1, n_variables
      ! Check individual convergence rates
       call convergence_pass_string( &
            array_names(i),  array_norm(i),  &
            array_norm_prev(i), pass_str_arr(i), &
            individual_tolerance )

      ! Weighted sum
      sum = sum + array_norm(i) / array_norm_prev(i)
      sum_prev = sum_prev + array_norm_prev(i) / array_norm_prev(i)
    end do

    call convergence_pass_string( &
         sum_name, sum, sum_prev, sum_pass_str, tolerance)

    pass_str = "PASS"
    do i= 1, n_variables
      if ( pass_str_arr(i) == "FAIL" ) pass_str = "FAIL"
    end do
    if ( sum_pass_str == "FAIL" ) pass_str = "FAIL"

    write(log_scratch_space,'("   test",A32," : ",A4)') trim(label), pass_str
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine array_convergence_rate_check

  !> @brief   Calculate and test the convergence rate.
  !> @details Calculate the convergence rate based on the norms
  !>          at two different iterations, and compare with the
  !>          expected value. Print out either PASS or FAIL, which
  !>          will then be used by the top level integration test.
  !> @param[in] norm_diff       Norm at second iteration
  !> @param[in] norm_diff_prev  Norm at first iteration
  !> @param[in] label           Name of the code being tested
  !> @param[in] tol             Tolerance value
  subroutine convergence_rate_check( norm_diff, norm_diff_prev, label, tol )

    implicit none

    real(r_def),           intent(in) :: norm_diff, norm_diff_prev
    character(str_def),    intent(in) :: label
    real(r_def), optional, intent(in) :: tol
    real(r_def) :: tolerance
    character(len=4) :: pass_str

    if ( present(tol) ) then
      tolerance = tol
    else
      tolerance = 1.0E-8_r_def
    end if

    call convergence_pass_string( &
         label, norm_diff, norm_diff_prev, pass_str, tolerance)

    write(log_scratch_space,'("   test",A32," : ",A4)') trim(label), pass_str
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

  end subroutine convergence_rate_check

end module tl_test_convergence_rate_check
