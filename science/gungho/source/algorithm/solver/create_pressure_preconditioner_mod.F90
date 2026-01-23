module create_pressure_preconditioner_mod

  !=============================================================================!
  !> @file create_pressure_preconditioner_mod.F90
  !> @brief Create iterative preconditioner for (Helmholtz) pressure problem
  !> @details This module provides a subroutine to create an iterative preconditioner for the
  !!          (Helmholtz) pressure problem, which is used in the semi-implicit preconditioner.
  !=============================================================================!

  use log_mod,                          only: log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR
  use constants_mod,                    only: i_def
  use pressure_operator_alg_mod,        only: pressure_operator_type
  use sci_preconditioner_mod,           only: abstract_preconditioner_type
  use pressure_precon_alg_mod,          only: pressure_preconditioner_type
  use pressure_diag_precon_alg_mod,     only: pressure_diag_preconditioner_type
  use multigrid_preconditioner_alg_mod, only: multigrid_preconditioner_type
  use sci_null_preconditioner_alg_mod,  only: null_preconditioner_type

  implicit none

  private
  public  :: create_pressure_preconditioner

  contains

!=============================================================================!
  !> @brief Create operator and preconditioner for (Helmholtz) pressure problem
  !> @details Called by init method of this module, but also by adj_semi_implicit_solver_alg_mod,
  !!          adjt_mixed_schur_preconditioner_alg_mod and adjt_mixed_solver_alg_mod
  !> @param[in]  fs_id                       ID of the pressure function space
  !> @param[out] pressure_operator_out       Output (Helmholtz) pressure operator
  !> @param[out] pressure_preconditioner_out Output (Helmholtz) pressure preconditioner
  !> @param[in]  level                       Multigrid level to create preconditioner on
  subroutine create_pressure_preconditioner( fs_id, pressure_operator_out, pressure_preconditioner_out, level )

    use helmholtz_solver_config_mod,   only: helmholtz_preconditioner => preconditioner, &
                                             preconditioner_none,                        &
                                             preconditioner_diagonal,                    &
                                             preconditioner_tridiagonal,                 &
                                             preconditioner_multigrid
    use function_space_mod,            only: function_space_type

    implicit none

    integer(kind=i_def), intent(in) :: fs_id
    integer(kind=i_def), intent(in) :: level

    ! Output operator and preconditioner for (Helmholtz) pressure problem
    type(pressure_operator_type),                     intent(out) :: pressure_operator_out
    class(abstract_preconditioner_type), allocatable, intent(out) :: pressure_preconditioner_out

    ! Vertical pressure preconditioner
    type(pressure_preconditioner_type) :: Hz_preconditioner

    integer(kind=i_def) :: n_fields

    pressure_operator_out = pressure_operator_type(level)

    call log_event( "create_pressure_preconditioner: starting", LOG_LEVEL_INFO )

    ! Allocate pressure preconditioner of correct type
    select case(helmholtz_preconditioner)
    case(PRECONDITIONER_NONE)
      call log_event( "No pressure preconditioner specified", LOG_LEVEL_INFO )
      allocate( pressure_preconditioner_out, &
                source = null_preconditioner_type() )
    case(PRECONDITIONER_DIAGONAL)
      call log_event( "Creating diagonal pressure preconditioner", LOG_LEVEL_INFO )
      allocate( pressure_preconditioner_out, &
                source = pressure_diag_preconditioner_type() )
    case(PRECONDITIONER_TRIDIAGONAL)
      call log_event( "Creating tridiagonal pressure preconditioner", LOG_LEVEL_INFO )
      allocate( pressure_preconditioner_out, &
                source = pressure_preconditioner_type(level) )
    case(PRECONDITIONER_MULTIGRID)
      call log_event( "Creating multigrid pressure preconditioner", LOG_LEVEL_INFO )
      if (level /= 1_i_def) &
        call log_event( "Pressure Multigrid has to be constructed on level 1", LOG_LEVEL_ERROR)

      Hz_preconditioner = pressure_preconditioner_type(level=1_i_def)
      n_fields = 1_i_def
      allocate( pressure_preconditioner_out, &
                source = multigrid_preconditioner_type( (/fs_id/),             &
                                                        n_fields,              &
                                                        pressure_operator_out, &
                                                        Hz_preconditioner ) )
    case default
      call log_event( "Unknown pressure preconditioner specified", LOG_LEVEL_ERROR)
    end select

    call log_event( "create_pressure_preconditioner: done", LOG_LEVEL_INFO )

  end subroutine create_pressure_preconditioner

  end module create_pressure_preconditioner_mod
