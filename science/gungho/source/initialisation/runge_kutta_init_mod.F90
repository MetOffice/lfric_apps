!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------

!>@brief Intialisation module for the Runge-Kutta method, sets the
!!       number of stages and coefficients
module runge_kutta_init_mod

  use constants_mod, only: r_def, i_def, r_tran, r_double

implicit none

  private
  real(kind=r_double), dimension(1,1) :: euler_ak
  real(kind=r_double), dimension(2,2) :: ssp2_ak
  real(kind=r_double), dimension(3,3) :: ssp3_ak
  real(kind=r_double), dimension(4,4) :: ssp4_ak
  real(kind=r_double), dimension(5,5) :: ssp5_ak

  public :: runge_kutta_init
  public :: get_rk_timestepping_weights
  public :: get_rk_transport_weights

contains

  !> @brief Compute the weights for the set of allowed Runge-Kutta schemes

  subroutine runge_kutta_init()

    implicit none

    ! Weights for each defined scheme
    euler_ak(1,1) = 1.0_r_double

    ssp2_ak(1,:) = (/ 1.0_r_double, 0.0_r_double /)
    ssp2_ak(2,:) = (/ 0.5_r_double, 0.5_r_double /)

    ssp3_ak(1,:) = (/ 1.0_r_double,  0.0_r_double,  0.0_r_double /)
    ssp3_ak(2,:) = (/ 0.25_r_double, 0.25_r_double, 0.0_r_double /)
    ssp3_ak(3,:) = (/ 1.0_r_double,  1.0_r_double,  4.0_r_double /)/6.0_r_double

    ssp4_ak(1,:) = (/ 0.5_r_double, 0.0_r_double, 0.0_r_double, 0.0_r_double /)
    ssp4_ak(2,:) = (/ 0.5_r_double, 0.5_r_double, 0.0_r_double, 0.0_r_double /)
    ssp4_ak(3,:) = (/ 1.0_r_double, 1.0_r_double, 1.0_r_double, 0.0_r_double /)/6.0_r_double
    ssp4_ak(4,:) = (/ 1.0_r_double, 1.0_r_double, 1.0_r_double, 3.0_r_double /)/6.0_r_double

    ssp5_ak(1,:) = (/ 0.377268915331368_r_double, 0.0_r_double,               &
                      0.0_r_double,               0.0_r_double,              0.0_r_double /)
    ssp5_ak(2,:) = (/ 0.377268915331368_r_double, 0.377268915331368_r_double, &
                      0.0_r_double,               0.0_r_double,              0.0_r_double /)
    ssp5_ak(3,:) = (/ 0.242995220537395_r_double, 0.242995220537395_r_double, &
                      0.242995220537395_r_double, 0.0_r_double,              0.0_r_double /)
    ssp5_ak(4,:) = (/ 0.153589067695126_r_double, 0.153589067695126_r_double, &
                      0.153589067695126_r_double, 0.23845893284629_r_double, 0.0_r_double /)
    ssp5_ak(5,:) = (/ 0.206734020864804_r_double, 0.206734020864804_r_double, &
                      0.117097251841844_r_double, 0.181802560120140_r_double, &
                      0.287632146308408_r_double /)

  end subroutine runge_kutta_init

  !> @brief Get the Butcher tableau and number of stages for a Runge-Kutta (RK)
  !!        timestepping scheme
  !> @param[out] num_rk_stage       Number of stages of RK scheme
  !> @param[out] ak                 Butcher tableau of RK weights
  !> @param[in]  runge_kutta_method Which RK method to choose
  subroutine get_rk_timestepping_weights(num_rk_stage, ak, runge_kutta_method)

    use log_mod,                 only: log_event,         &
                                       log_scratch_space, &
                                       LOG_LEVEL_ERROR
    use timestepping_config_mod, only: runge_kutta_method_forward_euler, &
                                       runge_kutta_method_ssp2,          &
                                       runge_kutta_method_ssp3,          &
                                       runge_kutta_method_ssp4,          &
                                       runge_kutta_method_ssp5
    implicit none

    integer(kind=i_def),           intent(out) :: num_rk_stage
    real(kind=r_def), allocatable, intent(out) :: ak(:,:)
    integer(kind=i_def),           intent(in)  :: runge_kutta_method

    select case(runge_kutta_method)

      case(runge_kutta_method_forward_euler)
        num_rk_stage = 1
        allocate ( ak (num_rk_stage,num_rk_stage) )
        ak = real( euler_ak, r_def )

      case(runge_kutta_method_ssp2)
        num_rk_stage = 2
        allocate ( ak (num_rk_stage,num_rk_stage) )
        ak = real( ssp2_ak, r_def )

      case(runge_kutta_method_ssp3)
        num_rk_stage = 3
        allocate ( ak (num_rk_stage,num_rk_stage) )
        ak = real( ssp3_ak, r_def )

      case(runge_kutta_method_ssp4)
        num_rk_stage = 4
        allocate ( ak (num_rk_stage,num_rk_stage) )
        ak = real( ssp4_ak, r_def )

      case(runge_kutta_method_ssp5)
        num_rk_stage = 5
        allocate ( ak (num_rk_stage,num_rk_stage) )
        ak = real( ssp5_ak, r_def )

      case default
        write( log_scratch_space, '(A,I3)' )  &
          'Invalid Runge Kutta timestepping method, stopping', runge_kutta_method
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select

  end subroutine get_rk_timestepping_weights

  !> @brief Get the Butcher tableau and number of stages for a Runge-Kutta (RK)
  !!        transport scheme
  !> @param[out] num_rk_stage       Number of stages of RK scheme
  !> @param[out] ak                 Butcher tableau of RK weights
  !> @param[in]  runge_kutta_method Which RK method to choose
  subroutine get_rk_transport_weights(num_rk_stage, ak, runge_kutta_method)

    use log_mod,              only: log_event,         &
                                    log_scratch_space, &
                                    LOG_LEVEL_ERROR
    use transport_config_mod, only: runge_kutta_method_forward_euler, &
                                    runge_kutta_method_ssp2,          &
                                    runge_kutta_method_ssp3,          &
                                    runge_kutta_method_ssp4,          &
                                    runge_kutta_method_ssp5
    implicit none

    integer(kind=i_def),           intent(out) :: num_rk_stage
    real(kind=r_tran), allocatable, intent(out) :: ak(:,:)
    integer(kind=i_def),           intent(in)  :: runge_kutta_method

    select case(runge_kutta_method)

      case(runge_kutta_method_forward_euler)
        num_rk_stage = 1
        allocate ( ak (num_rk_stage,num_rk_stage) )
        ak = real( euler_ak, r_tran)

      case(runge_kutta_method_ssp2)
        num_rk_stage = 2
        allocate ( ak (num_rk_stage,num_rk_stage) )
        ak = real( ssp2_ak, r_tran)

      case(runge_kutta_method_ssp3)
        num_rk_stage = 3
        allocate ( ak (num_rk_stage,num_rk_stage) )
        ak = real( ssp3_ak, r_tran)

      case(runge_kutta_method_ssp4)
        num_rk_stage = 4
        allocate ( ak (num_rk_stage,num_rk_stage) )
        ak = real( ssp4_ak, r_tran)

      case(runge_kutta_method_ssp5)
        num_rk_stage = 5
        allocate ( ak (num_rk_stage,num_rk_stage) )
        ak = real( ssp5_ak, r_tran)

      case default
        write( log_scratch_space, '(A,I3)' )  &
          'Invalid Runge Kutta transport method, stopping', runge_kutta_method
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select

  end subroutine get_rk_transport_weights
end module runge_kutta_init_mod
