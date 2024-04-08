!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!> @brief Applies the Random Parameter (RP) scheme

module stph_rp_main_alg_mod

  use constants_mod,        only: r_def, i_def, l_def
  use gungho_modeldb_mod,   only: modeldb_type

  implicit none

  private

  public stph_rp_init_alg
  public stph_rp_main_alg
  private rp_update_alg

contains
  !>@brief     Initialise the Random Parameter (RP) scheme
  !>@details   Set the initial values of the random numbers used
  !>           used to perturb the parameters in the RP scheme
  !>           See UMDP-081 for full scheme details
  !>
  !>@param[inout] modeldb  The structure that holds the model state
  subroutine stph_rp_init_alg(modeldb)

  use stochastic_physics_config_mod, only: rp_ran_max

  implicit none

  ! Arguments
  type(modeldb_type), intent(inout), target  :: modeldb

  ! Time-correlated random numbers used to perturb
  ! boundary layer parameter values.
  real(kind=r_def) :: rp_bl_a_ent_1_rand
  real(kind=r_def) :: rp_bl_cld_top_diffusion_rand
  real(kind=r_def) :: rp_bl_mix_rand
  real(kind=r_def) :: rp_bl_stab_rand

  ! Land-surface
  real(kind=r_def) :: rp_lsfc_lai_mult_rand
  real(kind=r_def) :: rp_lsfc_z0hm_pft_rand

  ! Micro-physics
  real(kind=r_def) :: rp_mp_ndrop_surf_rand

  ! Random numbers for each parameter
  real(kind=r_def) :: rp_rand(rp_ran_max)
  integer(kind=i_def) :: rp_count

  ! Initialise rp_count
  rp_count = 1

  ! Get set of random numbers and store in rp_rand
  call random_number(rp_rand)

  ! For each random number, initialise and adjust to
  ! vary between -1 and 1
  rp_bl_a_ent_1_rand = 2.0_r_def*rp_rand(rp_count) - 1.0_r_def
  rp_count = rp_count + 1

  rp_bl_cld_top_diffusion_rand = 2.0_r_def*rp_rand(rp_count) - 1.0_r_def
  rp_count = rp_count + 1

  rp_bl_mix_rand = 2.0_r_def*rp_rand(rp_count) - 1.0_r_def
  rp_count = rp_count + 1

  rp_bl_stab_rand = 2.0_r_def*rp_rand(rp_count) - 1.0_r_def
  rp_count = rp_count + 1

  rp_lsfc_lai_mult_rand = 2.0_r_def*rp_rand(rp_count) - 1.0_r_def
  rp_count = rp_count + 1

  rp_lsfc_z0hm_pft_rand = 2.0_r_def*rp_rand(rp_count) - 1.0_r_def
  rp_count = rp_count + 1

  rp_mp_ndrop_surf_rand = 2.0_r_def*rp_rand(rp_count) - 1.0_r_def

  ! Save to the modeldb for future use
  call modeldb%values%add_key_value('rp_bl_a_ent_1_rand', rp_bl_a_ent_1_rand)
  call modeldb%values%add_key_value('rp_bl_cld_top_diffusion_rand', &
                                     rp_bl_cld_top_diffusion_rand)
  call modeldb%values%add_key_value('rp_bl_mix_rand', rp_bl_mix_rand)
  call modeldb%values%add_key_value('rp_bl_stab_rand', rp_bl_stab_rand)
  call modeldb%values%add_key_value('rp_lsfc_lai_mult_rand', &
                                     rp_lsfc_lai_mult_rand)
  call modeldb%values%add_key_value('rp_lsfc_z0hm_pft_rand', &
                                     rp_lsfc_z0hm_pft_rand)
  call modeldb%values%add_key_value('rp_mp_ndrop_surf_rand', &
                                     rp_mp_ndrop_surf_rand)

  end subroutine stph_rp_init_alg

  !>@brief     Run the Random Parameter scheme
  !>@details   Update each parameter according to a stochastic process
  !>           See UMDP-081 for full scheme details
  !>
  !>@param[inout] modeldb  The structure that holds the model state

  subroutine stph_rp_main_alg(modeldb)



  ! Stochastic Physics namelist
  use stochastic_physics_config_mod, only: rp_bl_a_ent_1,                      &
                                           rp_bl_cbl_mix_fac,                  &
                                           rp_bl_cld_top_diffusion,            &
                                           rp_bl_stable_ri_coef,               &
                                           rp_bl_min_mix_length,               &
                                           rp_bl_neutral_mix_length,           &
                                           rp_bl_ricrit,                       &
                                           rp_bl_smag_coef,                    &
                                           rp_callfreq,                        &
                                           rp_cycle_in,                        &
                                           rp_cycle_out,                       &
                                           rp_cycle_tm,                        &
                                           rp_decorr_ts,                       &
                                           rp_lsfc_lai_mult_max,               &
                                           rp_lsfc_lai_mult_min,               &
                                           rp_lsfc_lai_mult,                   &
                                           rp_lsfc_z0hm_pft_max,               &
                                           rp_lsfc_z0hm_pft_min,               &
                                           rp_lsfc_z0hm_pft,                   &
                                           rp_mp_ndrop_surf,                   &
                                           rp_ran_max

  use microphysics_config_mod,   only : droplet_tpr

  use jules_surface_types_mod, only: npft
  use stph_time_corr_pert_alg_mod, only: stph_time_corr_pert_alg
  use extrusion_config_mod, only : number_of_layers

  use log_mod, only: log_event, log_scratch_space, LOG_LEVEL_INFO, LOG_LEVEL_DEBUG

  ! UM Namelists to pass through perturbed values
  use stochastic_physics_run_mod, only: a_ent_1_rp, cbl_mix_fac_rp, cs_rp,     &
       g0_rp, g1_rp, lai_mult_rp, lambda_min_rp, ndrop_surf_rp,                &
       par_mezcla_rp, ricrit_rp, z0hm_pft_rp

  implicit none

  ! Arguments
  type(modeldb_type), intent(inout), target  :: modeldb

  ! Model timestep
  integer(kind=i_def) :: timestep_number
  real(kind=r_def) :: dt

  ! Local variables
  integer(kind=i_def) :: lead_time
  real(kind=r_def) :: mu
  integer(kind=i_def) :: n, rp_count
  real(kind=r_def), parameter :: shock_amp = 2.0
  real(kind=r_def), parameter :: rand_max = 1.0
  real(kind=r_def), parameter :: rand_min = -1.0

  ! The random parameter is an array of length three
  ! The first value is the minimum, the second the parameter
  ! value that will be used in the physics schemes and the
  ! third is the maximum value
  integer(kind=i_def), parameter :: min_idx = 1  ! Index of min_rp
  integer(kind=i_def), parameter :: idx     = 2  ! Index of rp
  integer(kind=i_def), parameter :: max_idx = 3  ! Index of max_rp

  ! Random numbers for each parameter
  real(kind=r_def) :: rp_rand(rp_ran_max)

  ! Time-correlated random numbers used to perturb
  ! boundary layer parameter values.  Note that some parameters
  ! are perturbed together with the same random number
  real(kind=r_def), pointer :: rp_bl_a_ent_1_rand
  real(kind=r_def), pointer :: rp_bl_cld_top_diffusion_rand
  real(kind=r_def), pointer :: rp_bl_mix_rand
  real(kind=r_def), pointer :: rp_bl_stab_rand

  ! Perturbed BL parameter values to pass through to the physics schemes
  real(kind=r_def) :: rp_bl_a_ent_1_pert
  real(kind=r_def) :: rp_bl_cbl_mix_fac_pert
  real(kind=r_def) :: rp_bl_cld_top_diffusion_pert
  real(kind=r_def) :: rp_bl_stable_ri_coef_pert
  real(kind=r_def) :: rp_bl_min_mix_length_pert
  real(kind=r_def) :: rp_bl_neutral_mix_length_pert
  real(kind=r_def) :: rp_bl_ricrit_pert
  real(kind=r_def) :: rp_bl_smag_coef_pert

  ! Land-surface variables
  real(kind=r_def), pointer :: rp_lsfc_lai_mult_rand
  real(kind=r_def), pointer :: rp_lsfc_z0hm_pft_rand

  real(kind=r_def), allocatable :: rp_lsfc_lai_mult_pert(:)
  real(kind=r_def), allocatable :: rp_lsfc_z0hm_pft_pert(:)

  ! Micro-physics variables
  real(kind=r_def), pointer :: rp_mp_ndrop_surf_rand

  real(kind=r_def) :: rp_mp_ndrop_surf_pert

  ! Call scheme depending on rp_callfreq
  timestep_number = modeldb%clock%get_step()
  dt = modeldb%clock%get_seconds_per_step()
  lead_time = int(dt * (timestep_number - 1))

  if (mod(lead_time, rp_callfreq) == 0) then
    call log_event( 'gungho_driver: Running Random Parameter Scheme',  LOG_LEVEL_INFO )
  else
    ! Will need code in here to output parameter values to a dump or
    ! other output file, even if the time doesn't coincide with the
    ! RP scheme running
    return
  end if

  ! Set up allocatable arrays and initialise values where necessary
  if (.not. allocated(rp_lsfc_lai_mult_pert)) then
    allocate(rp_lsfc_lai_mult_pert(1:npft))
  end if
  if (.not. allocated(rp_lsfc_z0hm_pft_pert)) then
    allocate(rp_lsfc_z0hm_pft_pert(1:npft))
  end if

  ! Calculate auto-correlation coefficient
  mu = 1.0_r_def - (rp_callfreq / rp_decorr_ts)

  ! Initialise rp_count
  rp_count = 1

  ! Get set of random numbers and store in rp_rand
  call random_number(rp_rand)

  ! Extract the time-correlated random numbers from the modeldb
  call modeldb%values%get_value('rp_bl_a_ent_1_rand', rp_bl_a_ent_1_rand)
  call modeldb%values%get_value('rp_bl_cld_top_diffusion_rand', &
                                 rp_bl_cld_top_diffusion_rand)
  call modeldb%values%get_value('rp_bl_mix_rand', rp_bl_mix_rand)
  call modeldb%values%get_value('rp_bl_stab_rand', rp_bl_stab_rand)
  call modeldb%values%get_value('rp_lsfc_lai_mult_rand', &
                                 rp_lsfc_lai_mult_rand)
  call modeldb%values%get_value('rp_lsfc_z0hm_pft_rand', &
                                 rp_lsfc_z0hm_pft_rand)
  call modeldb%values%get_value('rp_mp_ndrop_surf_rand', &
                                 rp_mp_ndrop_surf_rand)

  ! TODO in future ticket:
  ! - running from CRUNS
  ! - reading in values from a file if requested
  ! - writing out values to a file if requested
  ! - hook jobs

  ! Update parameters

  !
  ! Boundary layer parameters
  !
  call stph_time_corr_pert_alg(mu, rp_rand(rp_count), shock_amp,               &
                               rand_min, rand_max, rp_bl_a_ent_1_rand)
  rp_count = rp_count + 1
  call rp_update_alg(rp_bl_a_ent_1_rand, rp_bl_a_ent_1(min_idx),               &
                     rp_bl_a_ent_1(max_idx), rp_bl_a_ent_1(idx),               &
                     rp_bl_a_ent_1_pert)
  write( log_scratch_space, '("stph_rp_main: rp_bl_a_ent_1: ", 2ES16.8 )' )    &
         rp_bl_a_ent_1_pert, rp_bl_a_ent_1(idx)
  call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

  call stph_time_corr_pert_alg(mu, rp_rand(rp_count), shock_amp,               &
                               rand_min, rand_max,                             &
                               rp_bl_cld_top_diffusion_rand)
  rp_count = rp_count + 1
  call rp_update_alg(rp_bl_cld_top_diffusion_rand,                             &
                     rp_bl_cld_top_diffusion(min_idx),                         &
                     rp_bl_cld_top_diffusion(max_idx),                         &
                     rp_bl_cld_top_diffusion(idx),                             &
                     rp_bl_cld_top_diffusion_pert)
  write( log_scratch_space,                                                    &
         '("stph_rp_main: rp_bl_cld_top_diffusion: ", 2ES16.8 )' )             &
         rp_bl_cld_top_diffusion_pert, rp_bl_cld_top_diffusion(idx)
  call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

  ! Perturb all mixing lengths together
  call stph_time_corr_pert_alg(mu, rp_rand(rp_count), shock_amp,               &
                               rand_min, rand_max, rp_bl_mix_rand)
  rp_count = rp_count + 1
  call rp_update_alg(rp_bl_mix_rand, rp_bl_cbl_mix_fac(min_idx),               &
                     rp_bl_cbl_mix_fac(max_idx), rp_bl_cbl_mix_fac(idx),       &
                     rp_bl_cbl_mix_fac_pert)
  write( log_scratch_space, '("stph_rp_main: rp_bl_cbl_mix_fac: ", 2ES16.8 )' )&
         rp_bl_cbl_mix_fac_pert, rp_bl_cbl_mix_fac(idx)
  call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

  call rp_update_alg(rp_bl_mix_rand, rp_bl_min_mix_length(min_idx),            &
                     rp_bl_min_mix_length(max_idx), rp_bl_min_mix_length(idx), &
                     rp_bl_min_mix_length_pert)
  write( log_scratch_space,                                                    &
         '("stph_rp_main: rp_bl_min_mix_length: ", 2ES16.8 )' )                &
         rp_bl_min_mix_length_pert, rp_bl_min_mix_length(idx)
  call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

  call rp_update_alg(rp_bl_mix_rand, rp_bl_neutral_mix_length(min_idx),        &
                     rp_bl_neutral_mix_length(max_idx),                        &
                     rp_bl_neutral_mix_length(idx),                            &
                     rp_bl_neutral_mix_length_pert)
  write( log_scratch_space,                                                    &
         '("stph_rp_main: rp_bl_neutral_mix_length: ", 2ES16.8 )' )            &
         rp_bl_neutral_mix_length_pert, rp_bl_neutral_mix_length(idx)
  call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

  call rp_update_alg(rp_bl_mix_rand, rp_bl_smag_coef(min_idx),                 &
                     rp_bl_smag_coef(max_idx), rp_bl_smag_coef(idx),           &
                     rp_bl_smag_coef_pert)
  write( log_scratch_space, '("stph_rp_main: rp_bl_smag_coef: ", 2ES16.8 )' )  &
         rp_bl_smag_coef_pert, rp_bl_smag_coef(idx)
  call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

  ! Perturb stability functions together
  call stph_time_corr_pert_alg(mu, rp_rand(rp_count), shock_amp,               &
                               rand_min, rand_max, rp_bl_stab_rand)
  rp_count = rp_count + 1
  call rp_update_alg(rp_bl_stab_rand, rp_bl_stable_ri_coef(min_idx),           &
                     rp_bl_stable_ri_coef(max_idx), rp_bl_stable_ri_coef(idx), &
                     rp_bl_stable_ri_coef_pert)
  write( log_scratch_space,                                                    &
         '("stph_rp_main: rp_bl_stable_ri_coef: ", 2ES16.8 )' )                &
         rp_bl_stable_ri_coef_pert, rp_bl_stable_ri_coef(idx)
  call log_event( log_scratch_space, LOG_LEVEL_DEBUG )

  call rp_update_alg(rp_bl_stab_rand, rp_bl_ricrit(min_idx),                   &
                     rp_bl_ricrit(max_idx), rp_bl_ricrit(idx),                 &
                     rp_bl_ricrit_pert)
  write( log_scratch_space, '("stph_rp_main: rp_bl_ricrit: ", 2ES16.8 )' )     &
         rp_bl_ricrit_pert, rp_bl_ricrit(idx)
  call log_event( log_scratch_space, LOG_LEVEL_DEBUG )


  !
  ! Micro-phsyics parameters
  !
  if ( droplet_tpr ) then
    call stph_time_corr_pert_alg(mu, rp_rand(rp_count), shock_amp,             &
                                 rand_min, rand_max, rp_mp_ndrop_surf_rand)
    rp_count = rp_count + 1
    call rp_update_alg(rp_mp_ndrop_surf_rand, rp_mp_ndrop_surf(min_idx),       &
                       rp_mp_ndrop_surf(max_idx), rp_mp_ndrop_surf(idx),       &
                       rp_mp_ndrop_surf_pert)
    write( log_scratch_space,                                                  &
           '("stph_rp_main: rp_mp_ndrop_surf: ", 2ES16.8 )' )                  &
           rp_mp_ndrop_surf_pert, rp_mp_ndrop_surf(idx)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end if

  !
  ! Land-surface parameters
  !
  ! Land-surface parameters are input as an array for each plant function type
  ! (1:npft).  To ensure the value for each plant function
  ! type (PFT) is changing in the same way (for example, making the
  ! forecast smoother/rougher overall) we use the same random number
  ! to perturb each PFT for a given parameter.

  call stph_time_corr_pert_alg(mu, rp_rand(rp_count), shock_amp,               &
                               rand_min, rand_max, rp_lsfc_lai_mult_rand)
  rp_count = rp_count + 1
  do n = 1, npft
    call rp_update_alg(rp_lsfc_lai_mult_rand, rp_lsfc_lai_mult_min(n),         &
                       rp_lsfc_lai_mult_max(n), rp_lsfc_lai_mult(n),           &
                       rp_lsfc_lai_mult_pert(n))
    write( log_scratch_space,                                                  &
           '("stph_rp_main: rp_lsfc_lai_mult: ", I0, 2ES16.8 )' )              &
           n, rp_lsfc_lai_mult_pert(n), rp_lsfc_lai_mult(n)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end do

  call stph_time_corr_pert_alg(mu, rp_rand(rp_count), shock_amp,               &
                               rand_min, rand_max, rp_lsfc_z0hm_pft_rand)
  rp_count = rp_count + 1
  do n = 1, npft
    call rp_update_alg(rp_lsfc_z0hm_pft_rand, rp_lsfc_z0hm_pft_min(n),         &
                       rp_lsfc_z0hm_pft_max(n), rp_lsfc_z0hm_pft(n),           &
                       rp_lsfc_z0hm_pft_pert(n))
    write( log_scratch_space,                                                  &
           '("stph_rp_main: rp_lsfc_z0hm_pft: ", I0, 2ES16.8 )' )              &
           n, rp_lsfc_z0hm_pft_pert(n), rp_lsfc_z0hm_pft(n)
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end do

  ! Pass the values through to the UM schemes to use
  a_ent_1_rp = rp_bl_a_ent_1_pert
  cbl_mix_fac_rp = rp_bl_cbl_mix_fac_pert
  cs_rp = rp_bl_smag_coef_pert
  g0_rp = rp_bl_stable_ri_coef_pert
  g1_rp = rp_bl_cld_top_diffusion_pert
  lambda_min_rp = rp_bl_min_mix_length_pert
  ndrop_surf_rp = rp_mp_ndrop_surf_pert
  par_mezcla_rp = rp_bl_neutral_mix_length_pert
  ricrit_rp = rp_bl_ricrit_pert
  do n = 1, npft
    lai_mult_rp(n) = rp_lsfc_lai_mult_pert(n)
    z0hm_pft_rp(n) = rp_lsfc_z0hm_pft_pert(n)
  end do

  end subroutine stph_rp_main_alg


  !>@brief     Updates a given random parameter
  !>
  !>@details   This routine updates a given random parameter by scaling it
  !>           by a time correlated random perturbation
  !>
  !>@param[in]      rand_no     Time correlated random number
  !>@param[in]      param_min   Min value of parameter
  !>@param[in]      param_max   Max value of parameter
  !>@param[in]      param_0     Default value of parameter
  !>@param[in, out] param_pert  Perturbed parameter to update

  subroutine rp_update_alg(rand_no, param_min, param_max, param_0, param_pert)

  implicit none

  ! Arguments
  real(kind=r_def), intent(in) :: rand_no
  real(kind=r_def), intent(in) :: param_min
  real(kind=r_def), intent(in) :: param_max
  real(kind=r_def), intent(in) :: param_0
  real(kind=r_def), intent(out) :: param_pert

  if ( rand_no < 0.0 ) then
    param_pert = param_0 + (param_0 - param_min) * rand_no
  else
    param_pert = param_0 + (param_max - param_0) * rand_no
  end if

  end subroutine rp_update_alg

end module stph_rp_main_alg_mod
