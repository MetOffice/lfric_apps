!----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief   Parameters used in the adjoint tests.
!> @details Contains ranges for the linearisation state to be randomised in.
!----------------------------------------------------------------------------

module adjoint_test_parameters_mod

  use constants_mod, only : r_def

  implicit none

  private
  ! ls ranges
  public :: ls_u_range
  public :: ls_theta_range
  public :: ls_rho_range
  public :: ls_exner_range
  public :: ls_md1_range
  public :: ls_md2_range
  public :: ls_md3_range

  ! Ranges for ls.
  ! Some kernels produce wildly divergent values
  ! if the ls is not realistic. ls can still be
  ! randomly assigned, but in a sensible range
  ! to prevent these issues.
  real(r_def), dimension(2), parameter :: ls_u_range = (/ 0.0_r_def, 10.0_r_def /)
  real(r_def), dimension(2), parameter :: ls_theta_range = (/ 280.0_r_def, 340.0_r_def /)
  real(r_def), dimension(2), parameter :: ls_rho_range = (/ 0.0_r_def, 1.0_r_def /)
  real(r_def), dimension(2), parameter :: ls_exner_range = (/ 0.0_r_def, 1.0_r_def /)
  real(r_def), dimension(2), parameter :: ls_md1_range = (/ 0.0_r_def, 1.0_r_def /)
  real(r_def), dimension(2), parameter :: ls_md2_range = (/ 0.0_r_def, 1.0_r_def /)
  real(r_def), dimension(2), parameter :: ls_md3_range = (/ 0.0_r_def, 1.0_r_def /)

end module adjoint_test_parameters_mod
