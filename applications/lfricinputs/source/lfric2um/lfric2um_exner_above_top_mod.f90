! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfric2um_exner_above_top_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env,     only: int64, real64
use constants_mod,                     only: r_def

! lfric2um modules
use lfricinp_lfric_driver_mod,         only: local_rank
use lfric2um_initialise_um_mod,        only: um_output_file
use lfricinp_um_grid_mod,              only: um_grid

! lfricinp modules
use lfricinp_um_level_codes_mod,       only: lfricinp_get_num_levels
use lfricinp_check_shumlib_status_mod, only: shumlib
use lfricinp_stashmaster_mod,          only: stashcode_theta

! lfric modules
use log_mod, only: log_event, log_scratch_space, LOG_LEVEL_INFO, &
                   LOG_LEVEL_ERROR

implicit none

private
public :: lfric2um_exner_above_top

contains

subroutine exner_above_top_level_cell(theta, pi_km1, q, gravity, cp, repsilon, dz, tol, pi)

implicit none

real(kind=real64), intent(in)  :: theta    ! Potential temperature at top model level
real(kind=real64), intent(in)  :: pi_km1   ! Exner pressure in top model half-level
real(kind=real64), intent(in)  :: q        ! Specific humidity at top model level
real(kind=real64), intent(in)  :: gravity
real(kind=real64), intent(in)  :: cp       ! Specific heat at constant pressure
real(kind=real64), intent(in)  :: repsilon !
real(kind=real64), intent(in)  :: dz       ! 2 x (top model level - top model half-level)
real(kind=real64), intent(in)  :: tol      ! Solver tolerance
real(kind=real64), intent(out) :: pi       ! Exner pressure at half-level above model top

integer(kind=int64)               :: iteration, max_its
real(kind=real64)                 :: R_T, R_pi, d_T, d_pi, T, tmp1, tmp2, detInv
real(kind=real64), dimension(2,2) :: J, Jinv

pi = pi_km1
T  = theta*pi_km1

tmp1 = cp*(1.0_r_def + (1.0_r_def/repsilon - 1.0_r_def)*q)
tmp2 = 0.5_r_def*gravity*dz

max_its = 100
do iteration = 1, max_its

  R_pi = (tmp1*T + tmp2)*pi - (tmp1*T - tmp2)*pi_km1
  R_T  = T - 0.5_r_def*theta*(pi + pi_km1)

  J(1,1) = tmp1*T + tmp2       ! d_R_pi / d_pi
  J(1,2) = tmp1*(pi - pi_km1)  ! d_R_pi / d_T
  J(2,1) = -0.5_r_def*theta          ! d_R_T  / d_pi
  J(2,2) = 1.0_r_def                 ! d_R_T  / d_T

  detInv = 1.0_r_def / (J(1,1)*J(2,2) - J(1,2)*J(2,1))

  Jinv(1,1) =  detInv*J(2,2)
  Jinv(1,2) = -detInv*J(1,2)
  Jinv(2,1) = -detInv*J(2,1)
  Jinv(2,2) =  detInv*J(1,1)

  d_pi = Jinv(1,1)*R_pi + Jinv(1,2)*R_T
  d_T  = Jinv(2,1)*R_pi + Jinv(2,2)*R_T

  pi = pi - d_pi
  T  = T  - d_T

  if (abs(d_pi/pi) < tol) then
    exit
  end if

  if (iteration == max_its) then
    write(log_scratch_space, '(A)') "Error! Convergence failure in calculation of Exner above model top."
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

end do

end subroutine exner_above_top_level_cell

subroutine lfric2um_exner_above_top(theta, pi_km1, q, pi)
! Description:
!  Compute the Exner pressure at a half level above the model top, so that it
!  can get passed into the UM fieldsfile as required by VAR.
!  This is done by solving a 2x2 Newton problem for the Exner a half level
!  above the model top, and the temperature at the model top, from the potential
!  temperature and the humidity at the model top and the Exner pressure a half
!  level below the model top.
!  Assumption: The top model level is of constant thickness everywhere.
use planet_config_mod,          only: cp, gravity, epsilon
use lfricinp_um_parameters_mod, only: ldc_zsea_theta, ldc_zsea_rho

implicit none

real(kind=real64), dimension(:,:), intent(in)  :: theta, pi_km1, q
real(kind=real64), dimension(:,:), intent(out) :: pi

real(kind=real64)   :: tol, dz
integer(kind=int64) :: ix, iy, nx, ny, num_levels_theta
real(kind=real64), allocatable :: level_dep_c(:,:)
character(len=*), parameter :: routinename='lfric2um_exner_above_top'

! Compute the thickness between the upper layer half level and the half level
! immediately above the upper level. We assume that this upper level thickness
! is constant for all horizontal cells.
call shumlib(routinename// &
     '::get_level_dependent_constants', &
     um_output_file % get_level_dependent_constants(level_dep_c))

num_levels_theta = lfricinp_get_num_levels(um_output_file, stashcode_theta)
dz = 2.0_r_def*( level_dep_c(num_levels_theta, ldc_zsea_theta) - &
        level_dep_c(num_levels_theta-1, ldc_zsea_rho) )

nx = um_grid%num_p_points_x
ny = um_grid%num_p_points_y

tol = 1.0e-12_r_def

! Regridding is all done on the root process, so compute Exner above top level
! on this processor only as well
if (local_rank == 0) then
  do ix = 1, nx
    do iy = 1, ny
      call exner_above_top_level_cell(theta(ix,iy), pi_km1(ix,iy), q(ix,iy),   &
                                      gravity, cp, epsilon, dz, tol, pi(ix,iy))
    end do
  end do
end if

end subroutine lfric2um_exner_above_top

end module lfric2um_exner_above_top_mod
