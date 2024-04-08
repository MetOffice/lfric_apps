!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @brief Updates time-correlated random perturbation

module stph_time_corr_pert_alg_mod

  use constants_mod, only: r_def, l_def

  implicit none

  private

  public stph_time_corr_pert_alg

contains
  !>@brief     Updates a random perturbation according to a stochastic process
  !>
  !>@details   This routine updates a given perturbation using a simple
  !>           AR1 process.  This results in an updated perturbation
  !>           correlated in time with the original perturbation.
  !>
  !>           The perturbation is assumed to have zero mean with max and min
  !>           values specified at input.  It's assumed that the perturbed value
  !>           has already been set randomly between the min and max values.
  !>           In subsequent calls, the perturbation is updated and any perturbation
  !>           calculated outside of the max or min is reflected back into the specified range.
  !>
  !>           The algorithm to update the perturbation includes a random shock
  !>           term with a specified amplitude.  For the case when a time-correlated
  !>           random number sequence is required, the shock amplitude should be set
  !>           to 1.
  !>
  !>           The auto-correlation coefficient, mu, should be calculated outside of
  !>           this module and typically takes the form EXP(-delta_t/delta_tau)
  !>           where delta_t is the timestep between perturbation updates, and
  !>           delta_tau is the decorrelation timescale.
  !>
  !>           This routine is called by stochastic physics schemes.
  !>
  !>@param[in]      mu          Auto-correlation coefficient
  !>@param[in]      rand_numb   Random number
  !>@param[in]      shock_amp   Amplitude of shock term
  !>@param[in]      pert_min    Min value of parameter
  !>@param[in]      pert_max    Max value of parameter
  !>@param[in, out] pert        Parameter to update

  subroutine stph_time_corr_pert_alg(mu, rand_numb, shock_amp, &
                                     pert_min, pert_max, pert)

  implicit none

  ! Arguments
  real(kind=r_def), intent(in) :: mu
  real(kind=r_def), intent(in) :: rand_numb
  real(kind=r_def), intent(in) :: shock_amp
  real(kind=r_def), intent(in) :: pert_min
  real(kind=r_def), intent(in) :: pert_max
  real(kind=r_def), intent(inout) :: pert


  real(kind=r_def) :: rand_mult ! multiplier based on random number

  ! Adjust random number to vary between -1 and 1
  rand_mult = 2.0_r_def*rand_numb - 1.0_r_def

  pert = mu*pert + SQRT(1.0_r_def-mu*mu)*rand_mult*shock_amp
  if (pert > pert_max) pert = 2.0_r_def*pert_max - pert
  if (pert < pert_min) pert = 2.0_r_def*pert_min - pert

  end subroutine stph_time_corr_pert_alg
end module stph_time_corr_pert_alg_mod
