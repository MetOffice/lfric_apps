module create_pressure_solver_mod

  !=============================================================================!
  !> @file create_pressure_solver_mod.F90
  !> @brief Create iterative solver for (Helmholtz) pressure problem
  !> @details This module provides a subroutine to create an iterative solver for the
  !!          (Helmholtz) pressure problem, which is used in the semi-implicit solver.
  !=============================================================================!

  use sci_iterative_solver_mod, only: abstract_iterative_solver_type, &
                                      conjugate_gradient_type, bicgstab_type, &
                                      gmres_type, fgmres_type, gcr_type, &
                                      precondition_only_type, jacobi_type

  use sci_preconditioner_mod, only: abstract_preconditioner_type
  use pressure_operator_alg_mod, only: pressure_operator_type
  use log_mod, only: log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR

  implicit none

  private
  public  :: create_pressure_solver

  contains

!=============================================================================!
  !> @brief Create iterative solver for (Helmholtz) pressure problem
  !> @details Called by init method of this module, but also by
  !!          adjt_mixed_schur_preconditioner_alg_mod and adjt_mixed_solver_alg_mod
  !> @param[in]  pressure_operator_in       Input (Helmholtz) pressure operator
  !> @param[in]  pressure_preconditioner_in Input (Helmholtz) pressure preconditioner
  !> @param[out] pressure_solver_out        Output (Helmholtz) pressure solver
  subroutine create_pressure_solver( pressure_operator_in, pressure_preconditioner_in, pressure_solver_out )

    use helmholtz_solver_config_mod,   only: si_pressure_maximum_iterations,         &
                                             helmholtz_gcrk => gcrk,                 &
                                             si_pressure_tolerance,                  &
                                             si_pressure_a_tol,                      &
                                             helmholtz_method => method,             &
                                             method_cg,                              &
                                             method_bicgstab,                        &
                                             method_gmres,                           &
                                             method_fgmres,                          &
                                             method_gcr,                             &
                                             method_prec_only,                       &
                                             method_jacobi,                          &
                                             si_pressure_monitor_convergence =>      &
                                                              monitor_convergence,   &
                                             si_pressure_fail_on_non_converged =>    &
                                                              fail_on_non_converged, &
                                             si_pressure_jacobi_relaxation =>        &
                                                             jacobi_relaxation

    implicit none

    ! Input operator and preconditioner for (Helmholtz) pressure problem
    type(pressure_operator_type),                     intent(in) :: pressure_operator_in
    class(abstract_preconditioner_type), allocatable, intent(in) :: pressure_preconditioner_in

    ! Output iterative solver for (Helmholtz) pressure problem
    class(abstract_iterative_solver_type), allocatable, intent(out) :: pressure_solver_out

    call log_event( "create_pressure_solver: starting", LOG_LEVEL_INFO )

    ! Allocate pressure solver of correct type
    select case( helmholtz_method )
    case (METHOD_BICGSTAB)
      allocate( pressure_solver_out,                                     &
                source = bicgstab_type( pressure_operator_in,            &
                                        pressure_preconditioner_in,      &
                                        si_pressure_tolerance,           &
                                        si_pressure_a_tol,               &
                                        si_pressure_maximum_iterations,  &
                                        si_pressure_monitor_convergence, &
                                        si_pressure_fail_on_non_converged ) )
    case(METHOD_CG)
      allocate( pressure_solver_out,                                               &
                source = conjugate_gradient_type( pressure_operator_in,            &
                                                  pressure_preconditioner_in,      &
                                                  si_pressure_tolerance,           &
                                                  si_pressure_a_tol,               &
                                                  si_pressure_maximum_iterations,  &
                                                  si_pressure_monitor_convergence, &
                                                  si_pressure_fail_on_non_converged ) )
    case(METHOD_GMRES)
      allocate( pressure_solver_out,                                  &
                source = gmres_type( pressure_operator_in,            &
                                     pressure_preconditioner_in,      &
                                     helmholtz_gcrk,                  &
                                     si_pressure_tolerance,           &
                                     si_pressure_a_tol,               &
                                     si_pressure_maximum_iterations,  &
                                     si_pressure_monitor_convergence, &
                                     si_pressure_fail_on_non_converged ) )
    case(METHOD_FGMRES)
      allocate( pressure_solver_out,                                   &
                source = fgmres_type( pressure_operator_in,            &
                                      pressure_preconditioner_in,      &
                                      helmholtz_gcrk,                  &
                                      si_pressure_tolerance,           &
                                      si_pressure_a_tol,               &
                                      si_pressure_maximum_iterations,  &
                                      si_pressure_monitor_convergence, &
                                      si_pressure_fail_on_non_converged ) )
    case(METHOD_GCR)
      allocate( pressure_solver_out,                                &
                source = gcr_type( pressure_operator_in,            &
                                   pressure_preconditioner_in,      &
                                   helmholtz_gcrk,                  &
                                   si_pressure_tolerance,           &
                                   si_pressure_a_tol,               &
                                   si_pressure_maximum_iterations,  &
                                   si_pressure_monitor_convergence, &
                                   si_pressure_fail_on_non_converged ) )
    case(METHOD_PREC_ONLY)
      allocate( pressure_solver_out,                                         &
                source = precondition_only_type( pressure_operator_in,       &
                                                 pressure_preconditioner_in, &
                                                 si_pressure_monitor_convergence) )
    case(METHOD_JACOBI)
      allocate( pressure_solver_out,                                     &
                source = jacobi_type( pressure_operator_in,              &
                                      pressure_preconditioner_in,        &
                                      si_pressure_tolerance,             &
                                      si_pressure_a_tol,                 &
                                      si_pressure_maximum_iterations,    &
                                      si_pressure_monitor_convergence,   &
                                      si_pressure_fail_on_non_converged, &
                                      si_pressure_jacobi_relaxation ) )
    case default
      call log_event("Unknown pressure solver specified",LOG_LEVEL_ERROR)
    end select

    call log_event( "create_pressure_solver: done", LOG_LEVEL_INFO )

  end subroutine create_pressure_solver

end module create_pressure_solver_mod
