!-------------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Radiative properties are calculated
!>        from GLOMAP aerosol derived fields.
!>        These are required by the Socrates radiation scheme.

module radaer_kernel_mod

use argument_mod,      only: arg_type,                                         &
                             GH_FIELD, GH_REAL, GH_READ, GH_WRITE,             &
                             CELL_COLUMN, GH_INTEGER,                          &
                             ANY_DISCONTINUOUS_SPACE_1,                        &
                             ANY_DISCONTINUOUS_SPACE_2,                        &
                             ANY_DISCONTINUOUS_SPACE_3,                        &
                             ANY_DISCONTINUOUS_SPACE_4,                        &
                             ANY_DISCONTINUOUS_SPACE_5

use empty_data_mod,    only: empty_real_data

use fs_continuity_mod, only: WTHETA, W3

use kernel_mod,        only: kernel_type

use log_mod,           only: log_event, log_scratch_space, LOG_LEVEL_ERROR

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: radaer_kernel_type
  private
  type(arg_type) :: meta_args(80) = (/                &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! theta_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! exner_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),     & ! exner_in_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! rho_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! dz_in_wth
       ! trop_level
       arg_type(GH_FIELD, GH_INTEGER, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
       ! lit_fraction
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! n_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! ait_sol_su
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! ait_sol_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! ait_sol_om
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! n_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! acc_sol_su
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! acc_sol_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! acc_sol_om
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! acc_sol_ss
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! n_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cor_sol_su
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cor_sol_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cor_sol_om
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cor_sol_ss
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! n_ait_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! ait_ins_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! ait_ins_om
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! n_acc_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! acc_ins_du
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! n_cor_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! cor_ins_du
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! drydp_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! drydp_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! drydp_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! drydp_ait_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! drydp_acc_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! drydp_cor_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! wetdp_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! wetdp_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! wetdp_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! rhopar_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! rhopar_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! rhopar_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! rhopar_ait_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! rhopar_acc_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! rhopar_cor_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_wat_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_wat_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_wat_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_su_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_bc_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_om_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_su_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_bc_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_om_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_ss_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_su_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_bc_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_om_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_ss_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_bc_ait_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_om_ait_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_du_acc_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA), & ! pvol_du_cor_ins
       ! aer_mix_ratio
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
       ! aer_sw_absorption
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3), &
       ! aer_sw_scattering
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3), &
       ! aer_sw_asymmetry
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_3), &
       ! aer_lw_absorption
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_4), &
       ! aer_lw_scattering
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_4), &
       ! aer_lw_asymmetry
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_4), &
       ! aod_ukca_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), &
       ! aaod_ukca_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), &
       ! aod_ukca_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), &
       ! aaod_ukca_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), &
       ! aod_ukca_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), &
       ! aaod_ukca_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), &
       ! aod_ukca_ait_ins
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), &
       ! aaod_ukca_ait_ins
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), &
       ! aod_ukca_acc_ins
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), &
       ! aaod_ukca_acc_ins
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), &
       ! aod_ukca_cor_ins
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), &
       ! aaod_ukca_cor_ins
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5) &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: radaer_code
end type radaer_kernel_type

public :: radaer_code

contains

!> @brief Interface to glomap aerosol climatology scheme.
!> @param[in]     nlayers            The number of layers
!> @param[in]     theta_in_wth       Potential temperature field
!> @param[in]     exner_in_wth       Exner pressure
!>                                    in potential temperature space
!> @param[in]     exner_in_w3        Exner pressure
!>                                    in density space
!> @param[in]     rho_in_wth         Density field
!>                                    in potential temperature space
!> @param[in]     dz_in_wth          Depth of temperature space levels
!> @param[in]     trop_level         Level of tropopause
!> @param[in]     lit_fraction       Fraction of radiation timestep lit by SW
!> @param[in]     n_ait_sol          Climatology aerosol field
!> @param[in]     ait_sol_su         Climatology aerosol field
!> @param[in]     ait_sol_bc         Climatology aerosol field
!> @param[in]     ait_sol_om         Climatology aerosol field
!> @param[in]     n_acc_sol          Climatology aerosol field
!> @param[in]     acc_sol_su         Climatology aerosol field
!> @param[in]     acc_sol_bc         Climatology aerosol field
!> @param[in]     acc_sol_om         Climatology aerosol field
!> @param[in]     acc_sol_ss         Climatology aerosol field
!> @param[in]     n_cor_sol          Climatology aerosol field
!> @param[in]     cor_sol_su         Climatology aerosol field
!> @param[in]     cor_sol_bc         Climatology aerosol field
!> @param[in]     cor_sol_om         Climatology aerosol field
!> @param[in]     cor_sol_ss         Climatology aerosol field
!> @param[in]     n_ait_ins          Climatology aerosol field
!> @param[in]     ait_ins_bc         Climatology aerosol field
!> @param[in]     ait_ins_om         Climatology aerosol field
!> @param[in]     n_acc_ins          Climatology aerosol field
!> @param[in]     acc_ins_du         Climatology aerosol field
!> @param[in]     n_cor_ins          Climatology aerosol field
!> @param[in]     cor_ins_du         Climatology aerosol field
!> @param[in]     drydp_ait_sol      Median particle dry diameter (Ait_Sol)
!> @param[in]     drydp_acc_sol      Median particle dry diameter (Acc_Sol)
!> @param[in]     drydp_cor_sol      Median particle dry diameter (Cor_Sol)
!> @param[in]     drydp_ait_ins      Median particle dry diameter (Ait_Ins)
!> @param[in]     drydp_acc_ins      Median particle dry diameter (Acc_Ins)
!> @param[in]     drydp_cor_ins      Median particle dry diameter (Cor_Ins)
!> @param[in]     wetdp_ait_sol      Avg wet diameter (Ait_Sol)
!> @param[in]     wetdp_acc_sol      Avg wet diameter (Acc_Sol)
!> @param[in]     wetdp_cor_sol      Avg wet diameter (Cor_Sol)
!> @param[in]     rhopar_ait_sol     Particle density (Ait_Sol)
!> @param[in]     rhopar_acc_sol     Particle density (Acc_Sol)
!> @param[in]     rhopar_cor_sol     Particle density (Cor_Sol)
!> @param[in]     rhopar_ait_ins     Particle density (Ait_Ins)
!> @param[in]     rhopar_acc_ins     Particle density (Acc_Ins)
!> @param[in]     rhopar_cor_ins     Particle density (Cor_Ins)
!> @param[in]     pvol_wat_ait_sol   Partial volume of water (Ait_Sol)
!> @param[in]     pvol_wat_acc_sol   Partial volume of water (Acc_Sol)
!> @param[in]     pvol_wat_cor_sol   Partial volume of water (Cor_Sol)
!> @param[in]     pvol_su_ait_sol    Partial volume (Ait_Sol h2so4)
!> @param[in]     pvol_bc_ait_sol    Partial volume (Ait_Sol black carbon)
!> @param[in]     pvol_om_ait_sol    Partial volume (Ait_Sol organic matter)
!> @param[in]     pvol_su_acc_sol    Partial volume (Acc_Sol h2so4)
!> @param[in]     pvol_bc_acc_sol    Partial volume (Acc_Sol black carbon)
!> @param[in]     pvol_om_acc_sol    Partial volume (Acc_Sol organic matter)
!> @param[in]     pvol_ss_acc_sol    Partial volume (Acc_Sol sea salt)
!> @param[in]     pvol_su_cor_sol    Partial volume (Cor_Sol h2so4)
!> @param[in]     pvol_bc_cor_sol    Partial volume (Cor_Sol black carbon)
!> @param[in]     pvol_om_cor_sol    Partial volume (Cor_Sol organic matter)
!> @param[in]     pvol_ss_cor_sol    Partial volume (Cor_Sol sea salt)
!> @param[in]     pvol_bc_ait_ins    Partial volume (Ait_Ins black carbon)
!> @param[in]     pvol_om_ait_ins    Partial volume (Ait_Ins organic matter)
!> @param[in]     pvol_du_acc_ins    Partial volume (Acc_Ins dust)
!> @param[in]     pvol_du_cor_ins    Partial volume (Cor_Ins dust)
!> @param[in,out] aer_mix_ratio      MODE aerosol mixing ratios
!> @param[in,out] aer_sw_absorption  MODE aerosol SW absorption
!> @param[in,out] aer_sw_scattering  MODE aerosol SW scattering
!> @param[in,out] aer_sw_asymmetry   MODE aerosol SW asymmetry
!> @param[in,out] aer_lw_absorption  MODE aerosol LW absorption
!> @param[in,out] aer_lw_scattering  MODE aerosol LW scattering
!> @param[in,out] aer_lw_asymmetry   MODE aerosol LW asymmetry
!> @param[in,out] aod_ukca_ait_sol   Modal extinction aerosol opt depth
!> @param[in,out] aaod_ukca_ait_sol  Modal absorption aerosol opt depth
!> @param[in,out] aod_ukca_acc_sol   Modal extinction aerosol opt depth
!> @param[in,out] aaod_ukca_acc_sol  Modal absorption aerosol opt depth
!> @param[in,out] aod_ukca_cor_sol   Modal extinction aerosol opt depth
!> @param[in,out] aaod_ukca_cor_sol  Modal absorption aerosol opt depth
!> @param[in,out] aod_ukca_ait_ins   Modal extinction aerosol opt depth
!> @param[in,out] aaod_ukca_ait_ins  Modal absorption aerosol opt depth
!> @param[in,out] aod_ukca_acc_ins   Modal extinction aerosol opt depth
!> @param[in,out] aaod_ukca_acc_ins  Modal absorption aerosol opt depth
!> @param[in,out] aod_ukca_cor_ins   Modal extinction aerosol opt depth
!> @param[in,out] aaod_ukca_cor_ins  Modal absorption aerosol opt depth
!> @param[in]     ndf_wth            Number of degrees of freedom per cell for
!>                                    potential temperature space
!> @param[in]     undf_wth           Unique number of degrees of freedom for
!>                                    potential temperature space
!> @param[in]     map_wth            Dofmap for the cell at the base of the
!>                                    column for potential temperature space
!> @param[in]     ndf_w3             Number of degrees of freedom per cell for
!>                                    density space
!> @param[in]     undf_w3            Unique number of degrees of freedom for
!>                                    density space
!> @param[in]     map_w3             Dofmap for the cell at the base of the
!>                                    column for density space
!> @param[in]     ndf_2d             No. DOFs per cell for 2D space
!> @param[in]     undf_2d            No. unique DOFs for 2D space
!> @param[in]     map_2d             Dofmap for 2D space column base cell
!> @param[in]     ndf_mode           No. of DOFs per cell for mode space
!> @param[in]     undf_mode          No. unique of DOFs for mode space
!> @param[in]     map_mode           Dofmap for mode space column base cell
!> @param[in]     ndf_rmode_sw       No. of DOFs per cell for rmode_sw space
!> @param[in]     undf_rmode_sw      No. unique of DOFs for rmode_sw space
!> @param[in]     map_rmode_sw       Dofmap for rmode_sw space column base cell
!> @param[in]     ndf_rmode_lw       No. of DOFs per cell for rmode_lw space
!> @param[in]     undf_rmode_lw      No. unique of DOFs for rmode_lw space
!> @param[in]     map_rmode_lw       Dofmap for rmode_lw space column base cell
!> @param[in]     ndf_aod_wavel      No. DOFs per cell for aod_wavel
!> @param[in]     undf_aod_wavel     No. unique DOFs for aod_wavel
!> @param[in]     map_aod_wavel      Dofmap for the cell at the base of the
!>                                    column for aod_wavel

subroutine radaer_code( nlayers,                                               &
                        theta_in_wth,                                          &
                        exner_in_wth,                                          &
                        exner_in_w3,                                           &
                        rho_in_wth,                                            &
                        dz_in_wth,                                             &
                        trop_level,                                            &
                        lit_fraction,                                          &
                        n_ait_sol,                                             &
                        ait_sol_su,                                            &
                        ait_sol_bc,                                            &
                        ait_sol_om,                                            &
                        n_acc_sol,                                             &
                        acc_sol_su,                                            &
                        acc_sol_bc,                                            &
                        acc_sol_om,                                            &
                        acc_sol_ss,                                            &
                        n_cor_sol,                                             &
                        cor_sol_su,                                            &
                        cor_sol_bc,                                            &
                        cor_sol_om,                                            &
                        cor_sol_ss,                                            &
                        n_ait_ins,                                             &
                        ait_ins_bc,                                            &
                        ait_ins_om,                                            &
                        n_acc_ins,                                             &
                        acc_ins_du,                                            &
                        n_cor_ins,                                             &
                        cor_ins_du,                                            &
                        drydp_ait_sol,                                         &
                        drydp_acc_sol,                                         &
                        drydp_cor_sol,                                         &
                        drydp_ait_ins,                                         &
                        drydp_acc_ins,                                         &
                        drydp_cor_ins,                                         &
                        wetdp_ait_sol,                                         &
                        wetdp_acc_sol,                                         &
                        wetdp_cor_sol,                                         &
                        rhopar_ait_sol,                                        &
                        rhopar_acc_sol,                                        &
                        rhopar_cor_sol,                                        &
                        rhopar_ait_ins,                                        &
                        rhopar_acc_ins,                                        &
                        rhopar_cor_ins,                                        &
                        pvol_wat_ait_sol,                                      &
                        pvol_wat_acc_sol,                                      &
                        pvol_wat_cor_sol,                                      &
                        pvol_su_ait_sol,                                       &
                        pvol_bc_ait_sol,                                       &
                        pvol_om_ait_sol,                                       &
                        pvol_su_acc_sol,                                       &
                        pvol_bc_acc_sol,                                       &
                        pvol_om_acc_sol,                                       &
                        pvol_ss_acc_sol,                                       &
                        pvol_su_cor_sol,                                       &
                        pvol_bc_cor_sol,                                       &
                        pvol_om_cor_sol,                                       &
                        pvol_ss_cor_sol,                                       &
                        pvol_bc_ait_ins,                                       &
                        pvol_om_ait_ins,                                       &
                        pvol_du_acc_ins,                                       &
                        pvol_du_cor_ins,                                       &
                        aer_mix_ratio,                                         &
                        aer_sw_absorption,                                     &
                        aer_sw_scattering,                                     &
                        aer_sw_asymmetry,                                      &
                        aer_lw_absorption,                                     &
                        aer_lw_scattering,                                     &
                        aer_lw_asymmetry,                                      &
                        aod_ukca_ait_sol,                                      &
                        aaod_ukca_ait_sol,                                     &
                        aod_ukca_acc_sol,                                      &
                        aaod_ukca_acc_sol,                                     &
                        aod_ukca_cor_sol,                                      &
                        aaod_ukca_cor_sol,                                     &
                        aod_ukca_ait_ins,                                      &
                        aaod_ukca_ait_ins,                                     &
                        aod_ukca_acc_ins,                                      &
                        aaod_ukca_acc_ins,                                     &
                        aod_ukca_cor_ins,                                      &
                        aaod_ukca_cor_ins,                                     &
                        ndf_wth, undf_wth, map_wth,                            &
                        ndf_w3, undf_w3, map_w3,                               &
                        ndf_2d, undf_2d, map_2d,                               &
                        ndf_mode, undf_mode, map_mode,                         &
                        ndf_rmode_sw, undf_rmode_sw, map_rmode_sw,             &
                        ndf_rmode_lw, undf_rmode_lw, map_rmode_lw,             &
                        ndf_aod_wavel, undf_aod_wavel, map_aod_wavel )


  use constants_mod,                     only: r_def, i_def, r_um, i_um
  use aerosol_config_mod,                only: n_radaer_step
  use socrates_init_mod,                 only: n_sw_band,                      &
                                               sw_n_band_exclude,              &
                                               sw_index_exclude,               &
                                               n_lw_band,                      &
                                               lw_n_band_exclude,              &
                                               lw_index_exclude

  use um_physics_init_mod,               only: n_radaer_mode,                  &
                                               n_aer_mode_sw, n_aer_mode_lw

  use nlsizes_namelist_mod,              only: row_length, rows

  use ukca_mode_setup,                   only: nmodes, ncp_max,                &
                                               mode_ait_sol, mode_acc_sol,     &
                                               mode_cor_sol, mode_ait_insol,   &
                                               mode_acc_insol, mode_cor_insol, &
                                               cp_su,  cp_bc, cp_oc,           &
                                               cp_cl,  cp_du, cp_so,           &
                                               cp_no3, cp_nn, cp_nh4,          &
                                               ip_ukca_mode_aitken,            &
                                               ip_ukca_mode_accum,             &
                                               ip_ukca_mode_coarse

  use ukca_radaer_band_average_mod,      only: ukca_radaer_band_average

  use ukca_radaer_prepare_mod,           only: ukca_radaer_prepare

  use ukca_radaer_compute_aod_mod,       only: ukca_radaer_compute_aod

  use planet_config_mod,                 only: p_zero, kappa, gravity

  use ukca_radaer_precalc,               only: npd_ukca_aod_wavel

  use ukca_option_mod,                   only: do_not_prescribe

  use ukca_radaer_lfric_api_mod,         only: ukca_radaer_lfric_interface

  implicit none

  ! Arguments

  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_wth
  integer(kind=i_def), intent(in) :: undf_wth
  integer(kind=i_def), dimension(ndf_wth), intent(in) :: map_wth

  integer(kind=i_def), intent(in) :: ndf_w3
  integer(kind=i_def), intent(in) :: undf_w3
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3

  integer(kind=i_def), intent(in) :: ndf_2d
  integer(kind=i_def), intent(in) :: undf_2d
  integer(kind=i_def), dimension(ndf_2d), intent(in) :: map_2d

  integer(kind=i_def), intent(in) :: ndf_mode
  integer(kind=i_def), intent(in) :: undf_mode
  integer(kind=i_def), dimension(ndf_mode), intent(in) :: map_mode

  integer(kind=i_def), intent(in) :: ndf_rmode_sw
  integer(kind=i_def), intent(in) :: undf_rmode_sw
  integer(kind=i_def), dimension(ndf_rmode_sw), intent(in) :: map_rmode_sw

  integer(kind=i_def), intent(in) :: ndf_rmode_lw
  integer(kind=i_def), intent(in) :: undf_rmode_lw
  integer(kind=i_def), dimension(ndf_rmode_lw), intent(in) :: map_rmode_lw

  integer(kind=i_def), intent(in) :: ndf_aod_wavel
  integer(kind=i_def), intent(in) :: undf_aod_wavel
  integer(kind=i_def), dimension(ndf_aod_wavel), intent(in) :: map_aod_wavel

  real(kind=r_def), intent(in),    dimension(undf_wth)   :: theta_in_wth
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: exner_in_wth

  real(kind=r_def), intent(in),    dimension(undf_w3)    :: exner_in_w3
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: rho_in_wth
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: dz_in_wth

  integer(kind=i_def), intent(in), dimension(undf_2d)    :: trop_level
  real(kind=r_def), intent(in),    dimension(undf_2d)    :: lit_fraction
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: n_ait_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: ait_sol_su
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: ait_sol_bc
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: ait_sol_om
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: n_acc_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: acc_sol_su
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: acc_sol_bc
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: acc_sol_om
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: acc_sol_ss
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: n_cor_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: cor_sol_su
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: cor_sol_bc
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: cor_sol_om
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: cor_sol_ss
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: n_ait_ins
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: ait_ins_bc
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: ait_ins_om
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: n_acc_ins
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: acc_ins_du
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: n_cor_ins
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: cor_ins_du
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: drydp_ait_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: drydp_acc_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: drydp_cor_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: drydp_ait_ins
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: drydp_acc_ins
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: drydp_cor_ins
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: wetdp_ait_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: wetdp_acc_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: wetdp_cor_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: rhopar_ait_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: rhopar_acc_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: rhopar_cor_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: rhopar_ait_ins
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: rhopar_acc_ins
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: rhopar_cor_ins
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_wat_ait_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_wat_acc_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_wat_cor_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_su_ait_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_bc_ait_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_om_ait_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_su_acc_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_bc_acc_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_om_acc_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_ss_acc_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_su_cor_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_bc_cor_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_om_cor_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_ss_cor_sol
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_bc_ait_ins
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_om_ait_ins
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_du_acc_ins
  real(kind=r_def), intent(in),    dimension(undf_wth)   :: pvol_du_cor_ins
  real(kind=r_def), intent(inout), dimension(undf_mode)  :: aer_mix_ratio
  real(kind=r_def), intent(inout), dimension(undf_rmode_sw) :: aer_sw_absorption
  real(kind=r_def), intent(inout), dimension(undf_rmode_sw) :: aer_sw_scattering
  real(kind=r_def), intent(inout), dimension(undf_rmode_sw) :: aer_sw_asymmetry
  real(kind=r_def), intent(inout), dimension(undf_rmode_lw) :: aer_lw_absorption
  real(kind=r_def), intent(inout), dimension(undf_rmode_lw) :: aer_lw_scattering
  real(kind=r_def), intent(inout), dimension(undf_rmode_lw) :: aer_lw_asymmetry

  ! Diagnostic arguments
  real(kind=r_def), pointer, intent(inout) :: aod_ukca_ait_sol(:)
  real(kind=r_def), pointer, intent(inout) :: aaod_ukca_ait_sol(:)
  real(kind=r_def), pointer, intent(inout) :: aod_ukca_acc_sol(:)
  real(kind=r_def), pointer, intent(inout) :: aaod_ukca_acc_sol(:)
  real(kind=r_def), pointer, intent(inout) :: aod_ukca_cor_sol(:)
  real(kind=r_def), pointer, intent(inout) :: aaod_ukca_cor_sol(:)
  real(kind=r_def), pointer, intent(inout) :: aod_ukca_ait_ins(:)
  real(kind=r_def), pointer, intent(inout) :: aaod_ukca_ait_ins(:)
  real(kind=r_def), pointer, intent(inout) :: aod_ukca_acc_ins(:)
  real(kind=r_def), pointer, intent(inout) :: aaod_ukca_acc_ins(:)
  real(kind=r_def), pointer, intent(inout) :: aod_ukca_cor_ins(:)
  real(kind=r_def), pointer, intent(inout) :: aaod_ukca_cor_ins(:)

  ! Local variables for the kernel

  ! Note - n_ukca_mode excludes the GLOMAP nucleation mode
  ! Since nucleation is the first mode in GLOMAP, the subsequent modes 2-7
  ! have been reordered in RADAER as modes 1-6
  integer(i_um), parameter :: n_ukca_mode = 6
  integer(i_um), parameter :: n_ukca_cpnt = 17

  integer(i_um) :: npd_exclude_lw
  integer(i_um) :: npd_exclude_sw
  logical, parameter       :: l_exclude_sw = .true.
  logical, parameter       :: l_exclude_lw = .true.

  ! Prescribed single-scattering albedo dummy variables
  ! Make these namelist options later
  integer, parameter       :: i_ukca_radaer_prescribe_ssa = do_not_prescribe
  integer(i_um), parameter :: nd_prof_ssa = 1
  integer(i_um), parameter :: nd_layr_ssa = 1
  integer(i_um), parameter :: nd_band_ssa = 1
  real(r_um),dimension( nd_prof_ssa, nd_layr_ssa, nd_band_ssa ) ::             &
                                                          ukca_radaer_presc_ssa

  ! Loop counters
  integer(i_um) :: k, i, i_band, i_mode, i_rmode
  integer(i_um) :: mm, m, n_fields, nn_fields

  ! Need to pass through argument list
  integer(i_um) :: ncp_max_x_nmodes

  ! pressure on theta levels
  real(r_um),dimension( row_length, rows, nlayers ) :: p_theta_levels_um

  ! temperature on theta levels
  real(r_um),dimension( row_length, rows, nlayers ) :: t_theta_levels_um

  ! d_mass on theta levels
  real(r_um),dimension( row_length, rows, nlayers ) :: d_mass_theta_levels_um

  ! Allocation of arrays

  real(r_um), allocatable :: ukca_comp_vol_um(:,:,:)
  real(r_um), allocatable :: ukca_mix_ratio_um(:,:,:)
  real(r_um), allocatable :: ukca_dry_diam_um(:,:,:)
  real(r_um), allocatable :: ukca_wet_diam_um(:,:,:)
  real(r_um), allocatable :: ukca_modal_nbr_um(:,:,:)
  real(r_um), allocatable :: ukca_modal_rho_um(:,:,:)
  real(r_um), allocatable :: ukca_modal_vol_um(:,:,:)
  real(r_um), allocatable :: ukca_modal_wtv_um(:,:,:)

  real(r_um), allocatable :: ukca_mode_mix_ratio_um(:,:,:)

  real(r_um), allocatable :: aer_lw_absorption_um(:,:,:,:)
  real(r_um), allocatable :: aer_lw_scattering_um(:,:,:,:)
  real(r_um), allocatable :: aer_lw_asymmetry_um(:,:,:,:)
  real(r_um), allocatable :: aer_sw_absorption_um(:,:,:,:)
  real(r_um), allocatable :: aer_sw_scattering_um(:,:,:,:)
  real(r_um), allocatable :: aer_sw_asymmetry_um(:,:,:,:)
  real(r_um), allocatable :: aod_ukca_all_modes_um(:,:,:)
  real(r_um), allocatable :: sod_ukca_all_modes_um(:,:,:)
  real(r_um), allocatable :: aaod_ukca_all_modes_um(:,:,:)

  ! By convention, arrays are inverted in UM radiation code
  ! Since we are calling from LFRic, arrays will not be inverted
  ! This matters for determining whether a level is above the tropopause
  logical, parameter :: l_inverted = .false.
  integer(i_um) :: trindxrad_um( row_length * rows )

  ! Variables close to but not exactly 1 or -1 for bounding asymmetry
  real(r_def), parameter :: one_minus_eps = 1.0_r_def - epsilon(1.0_r_def)
  real(r_def), parameter :: minus1_plus_eps = -1.0_r_def + epsilon(1.0_r_def)

  !-----------------------------------------------------------------------

  ! UKCA modal optical depth diagnostics: full column
  real(r_um) :: aod_ukca_this_mode_um(  row_length*rows, npd_ukca_aod_wavel )
  ! Not yet included as diagnostic
  ! UKCA modal optical depth diagnostics: stratosphere
  real(r_um) :: sod_ukca_this_mode_um(  row_length*rows, npd_ukca_aod_wavel )
  ! UKCA modal absorption optical depth diagnostics: full column
  real(r_um) :: aaod_ukca_this_mode_um( row_length*rows, npd_ukca_aod_wavel )

  !-----------------------------------------------------------------------

  logical :: l_aod_ukca_ait_sol
  logical :: l_aaod_ukca_ait_sol
  logical :: l_aod_ukca_acc_sol
  logical :: l_aaod_ukca_acc_sol
  logical :: l_aod_ukca_cor_sol
  logical :: l_aaod_ukca_cor_sol
  logical :: l_aod_ukca_ait_ins
  logical :: l_aaod_ukca_ait_ins
  logical :: l_aod_ukca_acc_ins
  logical :: l_aaod_ukca_acc_ins
  logical :: l_aod_ukca_cor_ins
  logical :: l_aaod_ukca_cor_ins
  logical :: l_any_lit_points_um

  !-----------------------------------------------------------------------

  integer, parameter :: max_fldname_len = 40

  character(len=max_fldname_len), parameter, dimension(17) :: &
       pvol_comp_names = [ 'fldname_pvol_su_ait_sol' , &
                           'fldname_pvol_bc_ait_sol' , &
                           'fldname_pvol_om_ait_sol' , &
                           'fldname_pvol_su_acc_sol' , &
                           'fldname_pvol_bc_acc_sol' , &
                           'fldname_pvol_om_acc_sol' , &
                           'fldname_pvol_ss_acc_sol' , &
                           'fldname_pvol_du_acc_sol' , &
                           'fldname_pvol_su_cor_sol' , &
                           'fldname_pvol_bc_cor_sol' , &
                           'fldname_pvol_om_cor_sol' , &
                           'fldname_pvol_ss_cor_sol' , &
                           'fldname_pvol_du_cor_sol' , &
                           'fldname_pvol_bc_ait_ins' , &
                           'fldname_pvol_om_ait_ins' , &
                           'fldname_pvol_du_acc_ins' , &
                           'fldname_pvol_du_cor_ins' ]

  character(len=max_fldname_len), parameter, dimension(17) :: &
       comp_names = [ 'fldname_ait_sol_su' , &
                      'fldname_ait_sol_bc' , &
                      'fldname_ait_sol_om' , &
                      'fldname_acc_sol_su' , &
                      'fldname_acc_sol_bc' , &
                      'fldname_acc_sol_om' , &
                      'fldname_acc_sol_ss' , &
                      'fldname_acc_sol_du' , &
                      'fldname_cor_sol_su' , &
                      'fldname_cor_sol_bc' , &
                      'fldname_cor_sol_om' , &
                      'fldname_cor_sol_ss' , &
                      'fldname_cor_sol_du' , &
                      'fldname_ait_ins_bc' , &
                      'fldname_ait_ins_om' , &
                      'fldname_acc_ins_du' , &
                      'fldname_cor_ins_du' ]

  character(len=max_fldname_len), parameter, dimension(6) :: &
       mode_names = [ 'fldname_n_ait_sol' , &
                      'fldname_n_acc_sol' , &
                      'fldname_n_cor_sol' , &
                      'fldname_n_ait_ins' , &
                      'fldname_n_acc_ins' , &
                      'fldname_n_cor_ins' ]

  character(len=max_fldname_len), parameter, dimension(6) :: &
       rhopar_mode_names = [ 'fldname_rhopar_ait_sol' , &
                             'fldname_rhopar_acc_sol' , &
                             'fldname_rhopar_cor_sol' , &
                             'fldname_rhopar_ait_ins' , &
                             'fldname_rhopar_acc_ins' , &
                             'fldname_rhopar_cor_ins' ]

  character(len=max_fldname_len), parameter, dimension(6) :: &
       dry_diam_mode_names = [ 'fldname_drydp_ait_sol' , &
                               'fldname_drydp_acc_sol' , &
                               'fldname_drydp_cor_sol' , &
                               'fldname_drydp_ait_ins' , &
                               'fldname_drydp_acc_ins' , &
                               'fldname_drydp_cor_ins' ]

  character(len=max_fldname_len), parameter, dimension(6) :: &
       modal_volume_names = [ 'fldname_mod_vol_ait_sol' , &
                              'fldname_mod_vol_acc_sol' , &
                              'fldname_mod_vol_cor_sol' , &
                              'fldname_mod_vol_ait_ins' , &
                              'fldname_mod_vol_acc_ins' , &
                              'fldname_mod_vol_cor_ins' ]

  character(len=max_fldname_len), parameter, dimension(4) :: &
       ait_sol_volume_names = [ 'pvol_wat_ait_sol' , &
                                'pvol_su_ait_sol ' , &
                                'pvol_bc_ait_sol ' , &
                                'pvol_om_ait_sol ' ]

  character(len=max_fldname_len), parameter,dimension(5) :: &
       acc_sol_volume_names = [ 'pvol_wat_acc_sol' , &
                                'pvol_su_acc_sol ' , &
                                'pvol_bc_acc_sol ' , &
                                'pvol_om_acc_sol ' , &
                                'pvol_ss_acc_sol ' ]

  character(len=max_fldname_len), parameter,dimension(5) :: &
       cor_sol_volume_names = [ 'pvol_wat_cor_sol' , &
                                'pvol_su_cor_sol ' , &
                                'pvol_bc_cor_sol ' , &
                                'pvol_om_cor_sol ' , &
                                'pvol_ss_cor_sol ' ]

  character(len=max_fldname_len), parameter,dimension(2) :: &
       ait_ins_volume_names = [ 'pvol_bc_ait_ins' , &
                                'pvol_om_ait_ins' ]

  character(len=max_fldname_len), parameter, dimension(1) :: &
       acc_ins_volume_names = [ 'pvol_du_acc_ins' ]

  character(len=max_fldname_len), parameter, dimension(1) :: &
       cor_ins_volume_names = [ 'pvol_du_cor_ins' ]

  character(len=max_fldname_len), parameter, dimension(6) :: &
       modal_wtv_names = [ 'fldname_pvol_wat_ait_sol' , &
                           'fldname_pvol_wat_acc_sol' , &
                           'fldname_pvol_wat_cor_sol' , &
                           'fldname_pvol_wat_ait_ins' , &
                           'fldname_pvol_wat_acc_ins' , &
                           'fldname_pvol_wat_cor_ins' ]

  character(len=max_fldname_len), parameter, dimension(6) :: &
       wet_diam_mode_names = [ 'fldname_wetdp_ait_sol' , &
                               'fldname_wetdp_acc_sol' , &
                               'fldname_wetdp_cor_sol' , &
                               'fldname_wetdp_ait_ins' , &
                               'fldname_wetdp_acc_ins' , &
                               'fldname_wetdp_cor_ins' ]

  !-----------------------------------------------------------------------

  ncp_max_x_nmodes = ncp_max * nmodes

  npd_exclude_lw = SIZE( lw_index_exclude, 1 )
  npd_exclude_sw = SIZE( sw_index_exclude, 1 )

  ! Note that this is inverted compared to the UM
  ! This will be dealt with in ukca_radaer_band_average
  trindxrad_um(1) = trop_level( map_2d(1) )

  !-----------------------------------------------------------------------
  ! Populate ukca_radaer element arrays
  ! Note that nucleation mode gets ignored in some of these

  l_aod_ukca_ait_sol = ( .not. associated( aod_ukca_ait_sol, empty_real_data ))
  l_aaod_ukca_ait_sol= ( .not. associated(aaod_ukca_ait_sol, empty_real_data ))
  l_aod_ukca_acc_sol = ( .not. associated( aod_ukca_acc_sol, empty_real_data ))
  l_aaod_ukca_acc_sol= ( .not. associated(aaod_ukca_acc_sol, empty_real_data ))
  l_aod_ukca_cor_sol = ( .not. associated( aod_ukca_cor_sol, empty_real_data ))
  l_aaod_ukca_cor_sol= ( .not. associated(aaod_ukca_cor_sol, empty_real_data ))
  l_aod_ukca_ait_ins = ( .not. associated( aod_ukca_ait_ins, empty_real_data ))
  l_aaod_ukca_ait_ins= ( .not. associated(aaod_ukca_ait_ins, empty_real_data ))
  l_aod_ukca_acc_ins = ( .not. associated( aod_ukca_acc_ins, empty_real_data ))
  l_aaod_ukca_acc_ins= ( .not. associated(aaod_ukca_acc_ins, empty_real_data ))
  l_aod_ukca_cor_ins = ( .not. associated( aod_ukca_cor_ins, empty_real_data ))
  l_aaod_ukca_cor_ins= ( .not. associated(aaod_ukca_cor_ins, empty_real_data ))

  !-----------------------------------------------------------------------
  ! Allocation of arrays

  allocate( ukca_mode_mix_ratio_um( 1, nlayers, n_radaer_mode ) )
  allocate( aer_lw_absorption_um( 1, nlayers, n_radaer_mode, n_lw_band ) )
  allocate( aer_lw_scattering_um( 1, nlayers, n_radaer_mode, n_lw_band ) )
  allocate( aer_lw_asymmetry_um( 1, nlayers, n_radaer_mode, n_lw_band ) )
  allocate( aer_sw_absorption_um( 1, nlayers, n_radaer_mode, n_sw_band ) )
  allocate( aer_sw_scattering_um( 1, nlayers, n_radaer_mode, n_sw_band ) )
  allocate( aer_sw_asymmetry_um( 1, nlayers, n_radaer_mode, n_sw_band ) )
  allocate( aod_ukca_all_modes_um( 1, npd_ukca_aod_wavel, n_ukca_mode ) )
  allocate( sod_ukca_all_modes_um( 1, npd_ukca_aod_wavel, n_ukca_mode ) )
  allocate( aaod_ukca_all_modes_um( 1, npd_ukca_aod_wavel, n_ukca_mode ) )

  !-----------------------------------------------------------------------
  ! Segmentation and openmp would start here
  !-----------------------------------------------------------------------

  ! Whether we need to run shortwave band_average because lit or not
  l_any_lit_points_um = .false.
  if ( n_radaer_step > 1 ) then
    l_any_lit_points_um = .true.
  else
    if ( lit_fraction( map_2d(1) ) > 0.0_r_def ) then
      l_any_lit_points_um = .true.
    end if
  end if

  ! Note that this is inverted compared to the UM
  ! This will be dealt with in ukca_radaer_band_average
  trindxrad_um(1) = trop_level( map_2d(1) )

  !-----------------------------------------------------------------------
  ! Initialisation of prognostic variables and arrays
  !-----------------------------------------------------------------------

  do k = 1, nlayers
    p_theta_levels_um(1,1,k) = p_zero *                                        &
                             ( exner_in_wth(map_wth(1) + k) )**(1.0_r_um/kappa)
  end do

  do k = 1, nlayers
    t_theta_levels_um(1,1,k) = exner_in_wth(map_wth(1) + k) *                  &
                               theta_in_wth(map_wth(1) + k)
  end do

  !-----------------------------------------------------------------------

  ! -- ukca_comp_vol_um --
  n_fields = size(pvol_comp_names)
  allocate(ukca_comp_vol_um(n_fields, 1, nlayers))
  ukca_comp_vol_um = 0.0_r_um

  do m = 1, n_fields
    select case(trim(pvol_comp_names(m)))
    case('fldname_pvol_su_ait_sol')
      do k = 1, nlayers
        ukca_comp_vol_um(m, 1, k)= real( pvol_su_ait_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_pvol_bc_ait_sol')
      do k = 1, nlayers
        ukca_comp_vol_um(m, 1, k)= real( pvol_bc_ait_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_pvol_om_ait_sol')
      do k = 1, nlayers
        ukca_comp_vol_um(m, 1, k)= real( pvol_om_ait_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_pvol_su_acc_sol')
      do k = 1, nlayers
        ukca_comp_vol_um(m, 1, k)= real( pvol_su_acc_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_pvol_bc_acc_sol')
      do k = 1, nlayers
        ukca_comp_vol_um(m, 1, k)= real( pvol_bc_acc_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_pvol_om_acc_sol')
      do k = 1, nlayers
        ukca_comp_vol_um(m, 1, k)= real( pvol_om_acc_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_pvol_ss_acc_sol')
      do k = 1, nlayers
        ukca_comp_vol_um(m, 1, k)= real( pvol_ss_acc_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_pvol_du_acc_sol')
      do k = 1, nlayers
        ukca_comp_vol_um(m, 1, k)= 0.0_r_um
      end do
    case('fldname_pvol_su_cor_sol')
      do k = 1, nlayers
        ukca_comp_vol_um(m, 1, k)= real( pvol_su_cor_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_pvol_bc_cor_sol')
      do k = 1, nlayers
        ukca_comp_vol_um(m, 1, k)= real( pvol_bc_cor_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_pvol_om_cor_sol')
      do k = 1, nlayers
        ukca_comp_vol_um(m, 1, k)= real( pvol_om_cor_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_pvol_ss_cor_sol')
      do k = 1, nlayers
        ukca_comp_vol_um(m, 1, k)= real( pvol_ss_cor_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_pvol_du_cor_sol')
      do k = 1, nlayers
        ukca_comp_vol_um(m, 1, k)= 0.0_r_um
      end do
    case('fldname_pvol_bc_ait_ins')
      do k = 1, nlayers
        ukca_comp_vol_um(m, 1, k)= real( pvol_bc_ait_ins(map_wth(1) + k), r_um )
      end do
    case('fldname_pvol_om_ait_ins')
      do k = 1, nlayers
        ukca_comp_vol_um(m, 1, k)= real( pvol_om_ait_ins(map_wth(1) + k), r_um )
      end do
    case('fldname_pvol_du_acc_ins')
      do k = 1, nlayers
        ukca_comp_vol_um(m, 1, k)= real( pvol_du_acc_ins(map_wth(1) + k), r_um )
      end do
    case('fldname_pvol_du_cor_ins')
      do k = 1, nlayers
        ukca_comp_vol_um(m, 1, k)= real( pvol_du_cor_ins(map_wth(1) + k), r_um )
      end do
    case default
      write( log_scratch_space, '(A,A)' )                                      &
           'Missing required ukca_comp_vol_um : ', pvol_comp_names(m)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  ! -- ukca_mix_ratio_um --
  n_fields = size(comp_names)
  allocate(ukca_mix_ratio_um(n_fields, 1, nlayers))
  ukca_mix_ratio_um = 0.0_r_um

  do m = 1, n_fields
    select case(trim(comp_names(m)))
    case('fldname_ait_sol_su')
      do k = 1, nlayers
        ukca_mix_ratio_um(m, 1, k) = real( ait_sol_su(map_wth(1) + k), r_um )
      end do
    case('fldname_ait_sol_bc')
      do k = 1, nlayers
        ukca_mix_ratio_um(m, 1, k) = real( ait_sol_bc(map_wth(1) + k), r_um )
      end do
    case('fldname_ait_sol_om')
      do k = 1, nlayers
        ukca_mix_ratio_um(m, 1, k) = real( ait_sol_om(map_wth(1) + k), r_um )
      end do
    case('fldname_acc_sol_su')
      do k = 1, nlayers
        ukca_mix_ratio_um(m, 1, k) = real( acc_sol_su(map_wth(1) + k), r_um )
      end do
    case('fldname_acc_sol_bc')
      do k = 1, nlayers
        ukca_mix_ratio_um(m, 1, k) = real( acc_sol_bc(map_wth(1) + k), r_um )
      end do
    case('fldname_acc_sol_om')
      do k = 1, nlayers
        ukca_mix_ratio_um(m, 1, k) = real( acc_sol_om(map_wth(1) + k), r_um )
      end do
    case('fldname_acc_sol_ss')
      do k = 1, nlayers
        ukca_mix_ratio_um(m, 1, k) = real( acc_sol_ss(map_wth(1) + k), r_um )
      end do
    case('fldname_acc_sol_du')
      do k = 1, nlayers
        ukca_mix_ratio_um(m, 1, k) = 0.0_r_um
      end do
    case('fldname_cor_sol_su')
      do k = 1, nlayers
        ukca_mix_ratio_um(m, 1, k) = real( cor_sol_su(map_wth(1) + k), r_um )
      end do
    case('fldname_cor_sol_bc')
      do k = 1, nlayers
        ukca_mix_ratio_um(m, 1, k) = real( cor_sol_bc(map_wth(1) + k), r_um )
      end do
    case('fldname_cor_sol_om')
      do k = 1, nlayers
        ukca_mix_ratio_um(m, 1, k) = real( cor_sol_om(map_wth(1) + k), r_um )
      end do
    case('fldname_cor_sol_ss')
      do k = 1, nlayers
        ukca_mix_ratio_um(m, 1, k) = real( cor_sol_ss(map_wth(1) + k), r_um )
      end do
    case('fldname_cor_sol_du')
      do k = 1, nlayers
        ukca_mix_ratio_um(m, 1, k) = 0.0_r_um
      end do
    case('fldname_ait_ins_bc')
      do k = 1, nlayers
        ukca_mix_ratio_um(m, 1, k) = real( ait_ins_bc(map_wth(1) + k), r_um )
      end do
    case('fldname_ait_ins_om')
      do k = 1, nlayers
        ukca_mix_ratio_um(m, 1, k) = real( ait_ins_om(map_wth(1) + k), r_um )
      end do
    case('fldname_acc_ins_du')
      do k = 1, nlayers
        ukca_mix_ratio_um(m, 1, k) = real( acc_ins_du(map_wth(1) + k), r_um )
      end do
    case('fldname_cor_ins_du')
      do k = 1, nlayers
        ukca_mix_ratio_um(m, 1, k) = real( cor_ins_du(map_wth(1) + k), r_um )
      end do
    case default
      write( log_scratch_space, '(A,A)' )                                      &
           'Missing required ukca_mix_ratio_um : ', comp_names(m)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  ! -- ukca_modal_nbr_um --
  n_fields = size(mode_names)
  allocate(ukca_modal_nbr_um(1, nlayers, n_fields))
  ukca_modal_nbr_um = 0.0_r_um

  do m = 1, n_fields
    select case(trim(mode_names(m)))
    case('fldname_n_nuc_sol')
      write( log_scratch_space, '(A,A)' )                                      &
           'Radaer should not request soluble nucleation mode: ', mode_names(m)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    case('fldname_n_ait_sol')
      do k = 1, nlayers
        ukca_modal_nbr_um(1, k, m) = real( n_ait_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_n_acc_sol')
      do k = 1, nlayers
        ukca_modal_nbr_um(1, k, m) = real( n_acc_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_n_cor_sol')
      do k = 1, nlayers
        ukca_modal_nbr_um(1, k, m) = real( n_cor_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_n_ait_ins')
      do k = 1, nlayers
        ukca_modal_nbr_um(1, k, m) = real( n_ait_ins(map_wth(1) + k), r_um )
      end do
    case('fldname_n_acc_ins')
      do k = 1, nlayers
        ukca_modal_nbr_um(1, k, m) = real( n_acc_ins(map_wth(1) + k), r_um )
      end do
    case('fldname_n_cor_ins')
      do k = 1, nlayers
        ukca_modal_nbr_um(1, k, m) = real( n_cor_ins(map_wth(1) + k), r_um )
      end do
    case default
      write( log_scratch_space, '(A,A)' )                                      &
           'Missing required ukca_modal_nbr_um : ', mode_names(m)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  ! -- ukca_modal_rho_um --
  n_fields = size(rhopar_mode_names)
  allocate(ukca_modal_rho_um(1, nlayers, n_fields))
  ukca_modal_rho_um = 0.0_r_um

  do m = 1, n_fields
    select case(trim(rhopar_mode_names(m)))
    case('fldname_rhopar_ait_sol')
      do k = 1, nlayers
        ukca_modal_rho_um(1, k, m)= real( rhopar_ait_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_rhopar_acc_sol')
      do k = 1, nlayers
        ukca_modal_rho_um(1, k, m)= real( rhopar_acc_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_rhopar_cor_sol')
      do k = 1, nlayers
        ukca_modal_rho_um(1, k, m)= real( rhopar_cor_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_rhopar_ait_ins')
      do k = 1, nlayers
        ukca_modal_rho_um(1, k, m)= real( rhopar_ait_ins(map_wth(1) + k), r_um )
      end do
    case('fldname_rhopar_acc_ins')
      do k = 1, nlayers
        ukca_modal_rho_um(1, k, m)= real( rhopar_acc_ins(map_wth(1) + k), r_um )
      end do
    case('fldname_rhopar_cor_ins')
      do k = 1, nlayers
        ukca_modal_rho_um(1, k, m)= real( rhopar_cor_ins(map_wth(1) + k), r_um )
      end do
    case default
      write( log_scratch_space, '(A,A)' )                                      &
           'Missing required ukca_modal_nbr_um : ', rhopar_mode_names(m)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  ! -- ukca_dry_diam_um --
  n_fields = size(dry_diam_mode_names)
  allocate(ukca_dry_diam_um(1, nlayers, n_fields))
  ukca_dry_diam_um = 0.0_r_um

  do m = 1, n_fields
    select case(trim(dry_diam_mode_names(m)))
    case('fldname_drydp_ait_sol')
      do k = 1, nlayers
        ukca_dry_diam_um(1, k, m) = real( drydp_ait_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_drydp_acc_sol')
      do k = 1, nlayers
        ukca_dry_diam_um(1, k, m) = real( drydp_acc_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_drydp_cor_sol')
      do k = 1, nlayers
        ukca_dry_diam_um(1, k, m) = real( drydp_cor_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_drydp_ait_ins')
      do k = 1, nlayers
        ukca_dry_diam_um(1, k, m) = real( drydp_ait_ins(map_wth(1) + k), r_um )
      end do
    case('fldname_drydp_acc_ins')
      do k = 1, nlayers
        ukca_dry_diam_um(1, k, m) = real( drydp_acc_ins(map_wth(1) + k), r_um )
      end do
    case('fldname_drydp_cor_ins')
      do k = 1, nlayers
        ukca_dry_diam_um(1, k, m) = real( drydp_cor_ins(map_wth(1) + k), r_um )
      end do
    case default
      write( log_scratch_space, '(A,A)' )                                      &
           'Missing required ukca_dry_diam_um : ', dry_diam_mode_names(m)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  ! -- ukca_modal_vol_um --
  n_fields = size(modal_volume_names)
  allocate(ukca_modal_vol_um(1, nlayers, n_fields))
  ukca_modal_vol_um = 0.0_r_um

  do m = 1, n_fields

    select case( trim( modal_volume_names(m) ) )

    case('fldname_mod_vol_ait_sol')

      nn_fields = size(ait_sol_volume_names)

      do mm = 1, nn_fields

        select case( trim( ait_sol_volume_names(mm) ) )

        case( 'pvol_wat_ait_sol' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_wat_ait_sol( map_wth(1) + k )
          end do

        case( 'pvol_su_ait_sol' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_su_ait_sol( map_wth(1) + k )
          end do

        case( 'pvol_bc_ait_sol' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_bc_ait_sol( map_wth(1) + k )
          end do

        case( 'pvol_om_ait_sol' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_om_ait_sol( map_wth(1) + k )
          end do

        case( 'null' )
          write( log_scratch_space, '(A,A)' )                                  &
               'This mode should not require a contribution to pvol : ',       &
               modal_volume_names(m)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )

        case default
          write( log_scratch_space, '(A,A)' )                                  &
           'Missing required ukca_modal_vol_um : ', ait_sol_volume_names(mm)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end select

      end do

    case('fldname_mod_vol_acc_sol')

      nn_fields = size(acc_sol_volume_names)

      do mm = 1, nn_fields
        select case( trim( acc_sol_volume_names(mm) ) )
        case( 'pvol_wat_acc_sol' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_wat_acc_sol( map_wth(1) + k )
          end do

        case( 'pvol_su_acc_sol' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_su_acc_sol( map_wth(1) + k )
          end do

        case( 'pvol_bc_acc_sol' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_bc_acc_sol( map_wth(1) + k )
          end do

        case( 'pvol_om_acc_sol' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_om_acc_sol( map_wth(1) + k )
          end do

        case( 'pvol_ss_acc_sol' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_ss_acc_sol( map_wth(1) + k )
          end do

        case( 'null' )
          write( log_scratch_space, '(A,A)' )                                  &
               'This mode should not require a contribution to pvol : ',       &
               modal_volume_names(m)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )

        case default
          write( log_scratch_space, '(A,A)' )                                  &
           'Missing required ukca_modal_vol_um : ', acc_sol_volume_names(mm)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end select

      end do

    case('fldname_mod_vol_cor_sol')

      nn_fields = size(cor_sol_volume_names)

      do mm = 1, nn_fields
        select case( trim( cor_sol_volume_names(mm) ) )
        case( 'pvol_wat_cor_sol' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_wat_cor_sol( map_wth(1) + k )
          end do

        case( 'pvol_su_cor_sol' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_su_cor_sol( map_wth(1) + k )
          end do

        case( 'pvol_bc_cor_sol' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_bc_cor_sol( map_wth(1) + k )
          end do

        case( 'pvol_om_cor_sol' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_om_cor_sol( map_wth(1) + k )
          end do

        case( 'pvol_ss_cor_sol' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_ss_cor_sol( map_wth(1) + k )
          end do

        case( 'null' )
          write( log_scratch_space, '(A,A)' )                                  &
               'This mode should not require a contribution to pvol : ',       &
               modal_volume_names(m)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )

        case default
          write( log_scratch_space, '(A,A)' )                                  &
           'Missing required ukca_modal_vol_um : ', cor_sol_volume_names(mm)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end select

      end do

    case('fldname_mod_vol_ait_ins')

      nn_fields = size(ait_ins_volume_names)

      do mm = 1, nn_fields

        select case( trim( ait_ins_volume_names(mm) ) )

        case( 'pvol_bc_ait_ins' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_bc_ait_ins( map_wth(1) + k )
          end do

        case( 'pvol_om_ait_ins' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_om_ait_ins( map_wth(1) + k )
          end do

        case( 'null' )
          write( log_scratch_space, '(A,A)' )                                  &
               'This mode should not require a contribution to pvol : ',       &
               modal_volume_names(m)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )

        case default
          write( log_scratch_space, '(A,A)' )                                  &
           'Missing required ukca_modal_vol_um : ', ait_ins_volume_names(mm)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end select

      end do

    case('fldname_mod_vol_acc_ins')

      nn_fields = size(acc_ins_volume_names)

      do mm = 1, nn_fields

        select case( trim( acc_ins_volume_names(mm) ) )

        case( 'pvol_du_acc_ins' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_du_acc_ins( map_wth(1) + k )
          end do

        case( 'null' )
          write( log_scratch_space, '(A,A)' )                                  &
               'This mode should not require a contribution to pvol : ',       &
               modal_volume_names(m)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )

        case default
          write( log_scratch_space, '(A,A)' )                                  &
           'Missing required ukca_modal_vol_um : ', acc_ins_volume_names(mm)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end select

      end do

    case('fldname_mod_vol_cor_ins')

      nn_fields = size(cor_ins_volume_names)

      do mm = 1, nn_fields

        select case( trim( cor_ins_volume_names(mm) ) )

        case( 'pvol_du_cor_ins' )
          do k = 1, nlayers
            ukca_modal_vol_um(1, k, m) = ukca_modal_vol_um(1, k, m) +          &
                                         pvol_du_cor_ins( map_wth(1) + k )
          end do

        case( 'null' )
          write( log_scratch_space, '(A,A)' )                                  &
               'This mode should not require a contribution to pvol : ',       &
               modal_volume_names(m)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )

        case default
          write( log_scratch_space, '(A,A)' )                                  &
           'Missing required ukca_modal_vol_um : ', cor_ins_volume_names(mm)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end select

      end do

    case default
      write( log_scratch_space, '(A,A)' )                                      &
           'Missing required ukca_modal_vol_um : ', modal_volume_names(m)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select

  end do

  ! -- ukca_modal_wtv_um --
  n_fields = size(modal_wtv_names)
  allocate(ukca_modal_wtv_um(1, nlayers, n_fields))
  ukca_modal_wtv_um = 0.0_r_um

  do m = 1, n_fields
    select case(trim(modal_wtv_names(m)))
    case('fldname_pvol_wat_ait_sol')
      do k = 1, nlayers
        ukca_modal_wtv_um(1,k,m) =real( pvol_wat_ait_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_pvol_wat_acc_sol')
      do k = 1, nlayers
        ukca_modal_wtv_um(1,k,m) =real( pvol_wat_acc_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_pvol_wat_cor_sol')
      do k = 1, nlayers
        ukca_modal_wtv_um(1,k,m) =real( pvol_wat_cor_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_pvol_wat_ait_ins')
      do k = 1, nlayers
        ukca_modal_wtv_um(1,k,m) = 0.0_r_um
      end do
    case('fldname_pvol_wat_acc_ins')
      do k = 1, nlayers
        ukca_modal_wtv_um(1,k,m) = 0.0_r_um
      end do
    case('fldname_pvol_wat_cor_ins')
      do k = 1, nlayers
        ukca_modal_wtv_um(1,k,m) = 0.0_r_um
      end do
    case default
      write( log_scratch_space, '(A,A)' )                                      &
           'Missing required ukca_modal_wtv_um : ', modal_wtv_names(m)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  ! -- ukca_wet_diam_um --
  n_fields = size(wet_diam_mode_names)
  allocate(ukca_wet_diam_um(1, nlayers, n_fields))
  ukca_wet_diam_um = 0.0_r_um

  do m = 1, n_fields
    select case(trim(wet_diam_mode_names(m)))
    case('fldname_wetdp_ait_sol')
      do k = 1, nlayers
        ukca_wet_diam_um(1, k, m) = real( wetdp_ait_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_wetdp_acc_sol')
      do k = 1, nlayers
        ukca_wet_diam_um(1, k, m) = real( wetdp_acc_sol(map_wth(1) + k), r_um )
      end do
    case('fldname_wetdp_cor_sol')
      do k = 1, nlayers
        ukca_wet_diam_um(1, k, m) = real( wetdp_cor_sol(map_wth(1) + k), r_um )
      end do
    !!!!! Note that wet and dry diameter are the same for insoluble modes
    case('fldname_wetdp_ait_ins')
      do k = 1, nlayers
        ukca_wet_diam_um(1, k, m) = real( drydp_ait_ins(map_wth(1) + k), r_um )
      end do
    !!!!! Note that wet and dry diameter are the same for insoluble modes
    case('fldname_wetdp_acc_ins')
      do k = 1, nlayers
        ukca_wet_diam_um(1, k, m) = real( drydp_acc_ins(map_wth(1) + k), r_um )
      end do
    !!!!! Note that wet and dry diameter are the same for insoluble modes
    case('fldname_wetdp_cor_ins')
      do k = 1, nlayers
        ukca_wet_diam_um(1, k, m) = real( drydp_cor_ins(map_wth(1) + k), r_um )
      end do
    case default
      write( log_scratch_space, '(A,A)' )                                      &
           'Missing required ukca_wet_diam_um : ', wet_diam_mode_names(m)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  !------------------------------------------------
  ! Calculate mass thickness of vertical levels
  ! This duplicates calculation of d_mass from set_thermodynamic_kernel_mod
  if ( ( .not. associated( aod_ukca_ait_sol, empty_real_data ) ) .or.          &
       ( .not. associated( aaod_ukca_ait_sol, empty_real_data ) ) .or.         &
       ( .not. associated( aod_ukca_acc_sol, empty_real_data ) ) .or.          &
       ( .not. associated( aaod_ukca_acc_sol, empty_real_data ) ) .or.         &
       ( .not. associated( aod_ukca_cor_sol, empty_real_data ) ) .or.          &
       ( .not. associated( aaod_ukca_cor_sol, empty_real_data ) ) .or.         &
       ( .not. associated( aod_ukca_ait_ins, empty_real_data ) ) .or.          &
       ( .not. associated( aaod_ukca_ait_ins, empty_real_data ) ) .or.         &
       ( .not. associated( aod_ukca_acc_ins, empty_real_data ) ) .or.          &
       ( .not. associated( aaod_ukca_acc_ins, empty_real_data ) ) .or.         &
       ( .not. associated( aod_ukca_cor_ins, empty_real_data ) ) .or.          &
       ( .not. associated( aaod_ukca_cor_ins, empty_real_data ) ) ) then

    d_mass_theta_levels_um(1,1,1) = rho_in_wth(  map_wth(2) ) *                &
                                    ( dz_in_wth( map_wth(2) ) +                &
                                    dz_in_wth( map_wth(1) ) )

    do k = 2, nlayers - 1
      d_mass_theta_levels_um(1,1,k) = rho_in_wth( map_wth(1) + k ) *           &
                                       dz_in_wth( map_wth(1) + k )
    end do

    d_mass_theta_levels_um(1,1,nlayers) = p_zero *                             &
                                        exner_in_w3( map_w3(1) + nlayers-1 )** &
                                        ( 1.0_r_def / kappa ) / gravity
  end if


    !-----------------------------------------------------------------------
    ! Call UKCA modules
    !-----------------------------------------------------------------------

    CALL ukca_radaer_lfric_interface(                                          &
      ! Fixed array dimensions (input)
      (row_length*rows),                                                       &
      nlayers,                                                                 &
      npd_exclude_lw,                                                          &
      npd_exclude_sw,                                                          &
      npd_ukca_aod_wavel,                                                      &
      ncp_max_x_nmodes,                                                        &
      ! Actual array dimensions (input)
      n_ukca_mode,                                                             &
      n_ukca_cpnt,                                                             &
      ! Prescribed SSA dimensions
      nd_prof_ssa,                                                             &
      nd_layr_ssa,                                                             &
      nd_band_ssa,                                                             &
      ! Variables related to waveband exclusion
      l_exclude_lw,                                                            &
      l_exclude_sw,                                                            &
      ! Modal diameters from UKCA module (input)
      ukca_dry_diam_um,                                                        &
      ukca_wet_diam_um,                                                        &
      ! Other inputs from UKCA module (input)
      ukca_comp_vol_um,                                                        &
      ukca_modal_vol_um,                                                       &
      ukca_modal_rho_um,                                                       &
      ukca_modal_wtv_um,                                                       &
      ! Logical to describe orientation
      l_inverted,                                                              &
      ! Control option for prescribed single scattering albedo array
      i_ukca_radaer_prescribe_ssa,                                             &
      ! Model level of the tropopause (input)
      trindxrad_um,                                                            &
      ! Whether we need to run shortwave band_average because lit or not
      l_any_lit_points_um,                                                     &
      ! Prescription of single-scattering albedo
      ukca_radaer_presc_ssa,                                                   &
      ! Input Component mass-mixing ratios
      ukca_mix_ratio_um,                                                       &
      ! Input modal number concentrations
      ukca_modal_nbr_um,                                                       &
      ! Input Pressure and temperature
      p_theta_levels_um, t_theta_levels_um,                                    &
      ! Which aerosol optical depth diagnostics to calculate
      l_aod_ukca_ait_sol, l_aaod_ukca_ait_sol,                                 &
      l_aod_ukca_acc_sol, l_aaod_ukca_acc_sol,                                 &
      l_aod_ukca_cor_sol, l_aaod_ukca_cor_sol,                                 &
      l_aod_ukca_ait_ins, l_aaod_ukca_ait_ins,                                 &
      l_aod_ukca_acc_ins, l_aaod_ukca_acc_ins,                                 &
      l_aod_ukca_cor_ins, l_aaod_ukca_cor_ins,                                 &
      ! Mass thickness of layers
      d_mass_theta_levels_um,                                                  &
      ! Modal mass-mixing ratios (input output)
      ukca_mode_mix_ratio_um,                                                  &
      ! Band-averaged optical properties (output)
      aer_lw_absorption_um,                                                    &
      aer_sw_absorption_um,                                                    &
      aer_lw_scattering_um,                                                    &
      aer_sw_scattering_um,                                                    &
      aer_lw_asymmetry_um,                                                     &
      aer_sw_asymmetry_um,                                                     &
      aod_ukca_all_modes_um,                                                   &
      aaod_ukca_all_modes_um )

    !-----------------------------------------------------------------------
    ! Convert back to LFRic arrays
    !-----------------------------------------------------------------------

    ! MODE aerosol mixing ratios
    do i_mode = 1, n_radaer_mode
      do k = 1, nlayers

          aer_mix_ratio( map_mode(1) + ( (i_mode-1)*(nlayers+1) ) + k ) =   &
                                         ukca_mode_mix_ratio_um( 1, k, i_mode )
      end do
    end do

    ! Socrates arrays filled with MODE aerosol optical properties in bands
    i_rmode = 0
    do i_band = 1, n_lw_band
      ! Fill the radaer modes within this band
      do i_mode = 1, n_radaer_mode
        i_rmode = i_rmode + 1
        do k = 1, nlayers

            aer_lw_absorption( map_rmode_lw(1) +                               &
                               ( (i_rmode-1)*(nlayers+1) ) + k ) =             &
                               aer_lw_absorption_um( 1, k, i_mode, i_band )

            aer_lw_scattering( map_rmode_lw(1) +                               &
                               ( (i_rmode-1)*(nlayers+1) ) + k ) =             &
                               aer_lw_scattering_um( 1, k, i_mode, i_band )

            aer_lw_asymmetry( map_rmode_lw(1)  +                               &
                              ( (i_rmode-1)*(nlayers+1) ) + k ) =              &
                              max(minus1_plus_eps, min(one_minus_eps,          &
                              aer_lw_asymmetry_um(   1, k, i_mode, i_band ) ) )

        end do ! n_layers
      end do ! n_radaer_mode

      ! If there are additional aerosol modes not associated with radaer
      ! (e.g. from easyaerosol) then i_rmode needs advancing past them
      ! before starting on the next radiation band.
      if (n_aer_mode_lw > n_radaer_mode) then
        i_rmode = i_rmode + n_aer_mode_lw - n_radaer_mode
      end if

    end do ! n_lw_band

    ! Only calculate SW on lit points
    ! If superstepping (n_radaer_step>1) then need to calculate on all points
    ! for use when the sun moves later

    ! Socrates arrays filled with MODE aerosol optical properties in bands
    i_rmode = 0
    do i_band = 1, n_sw_band

      ! Fill the radaer modes within this band
      do i_mode = 1, n_radaer_mode
        i_rmode = i_rmode + 1

        do k = 1, nlayers
            if ( ( lit_fraction(map_2d(1)) > 0.0_r_def  ) .or.             &
                 n_radaer_step > 1 ) then

              aer_sw_absorption( map_rmode_sw(1)     +                         &
                                 ( (i_rmode-1)*(nlayers+1) ) + k ) =           &
                                 aer_sw_absorption_um( 1, k, i_mode, i_band )

              aer_sw_scattering( map_rmode_sw(1)     +                         &
                                 ( (i_rmode-1)*(nlayers+1) ) + k ) =           &
                                 aer_sw_scattering_um( 1, k, i_mode, i_band )

              aer_sw_asymmetry( map_rmode_sw(1)      +                         &
                                ( (i_rmode-1)*(nlayers+1) ) + k ) =            &
                                max(minus1_plus_eps, min(one_minus_eps,        &
                                aer_sw_asymmetry_um(  1, k, i_mode, i_band ) ) )

            ! unlit points
            else

              aer_sw_absorption( map_rmode_sw(1) +                             &
                   ( (i_rmode-1)*(nlayers+1) ) + k ) =  1.0_r_def

              aer_sw_scattering( map_rmode_sw(1) +                             &
                   ( (i_rmode-1)*(nlayers+1) ) + k ) =  1.0_r_def

              aer_sw_asymmetry( map_rmode_sw(1)  +                             &
                   ( (i_rmode-1)*(nlayers+1) ) + k ) = one_minus_eps

            end if
        end do ! nlayers
      end do ! n_radaer_mode

      ! If there are additional aerosol modes not associated with radaer
      ! (e.g. from easyaerosol) then i_rmode needs advancing past them
      ! before starting on the next radiation band.
      if (n_aer_mode_sw > n_radaer_mode) then
        i_rmode = i_rmode + n_aer_mode_sw - n_radaer_mode
      end if

    end do ! n_sw_bands

    !------------------------------------------------
    ! Now calculate aod and aaod for Aitken Soluble mode

    if ( .not. associated( aod_ukca_ait_sol, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, row_length
          aod_ukca_ait_sol( map_aod_wavel(i) + k - 1 ) =                       &
                                      aod_ukca_all_modes_um(i,k,mode_ait_sol-1)
        end do
      end do
    end if

    if ( .not. associated( aaod_ukca_ait_sol, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, row_length
          aaod_ukca_ait_sol( map_aod_wavel(i) + k - 1 ) =                      &
                                     aaod_ukca_all_modes_um(i,k,mode_ait_sol-1)
        end do
      end do
    end if

    !------------------------------------------------
    ! Now calculate aod and aaod for Accumulation Soluble mode

    if ( .not. associated( aod_ukca_acc_sol, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, row_length
          aod_ukca_acc_sol( map_aod_wavel(i) + k - 1 ) =                       &
                                      aod_ukca_all_modes_um(i,k,mode_acc_sol-1)
        end do
      end do
    end if

    if ( .not. associated( aaod_ukca_acc_sol, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, row_length
          aaod_ukca_acc_sol( map_aod_wavel(i) + k - 1 ) =                      &
                                     aaod_ukca_all_modes_um(i,k,mode_acc_sol-1)
        end do
      end do
    end if

    !------------------------------------------------
    ! Now calculate aod and aaod for Coarse Soluble mode

    if ( .not. associated( aod_ukca_cor_sol, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, row_length
          aod_ukca_cor_sol( map_aod_wavel(i) + k - 1 ) =                       &
                                      aod_ukca_all_modes_um(i,k,mode_cor_sol-1)
        end do
      end do
    end if

    if ( .not. associated( aaod_ukca_cor_sol, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, row_length
          aaod_ukca_cor_sol( map_aod_wavel(i) + k - 1 ) =                      &
                                     aaod_ukca_all_modes_um(i,k,mode_cor_sol-1)
        end do
      end do
    end if

    !------------------------------------------------
    ! Now calculate aod and aaod for Aitken Insoluble mode

    if ( .not. associated( aod_ukca_ait_ins, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, row_length
          aod_ukca_ait_ins( map_aod_wavel(i) + k - 1 ) =                       &
                                    aod_ukca_all_modes_um(i,k,mode_ait_insol-1)
        end do
      end do
    end if

    if ( .not. associated( aaod_ukca_ait_ins, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, row_length
          aaod_ukca_ait_ins( map_aod_wavel(i) + k - 1 ) =                      &
                                   aaod_ukca_all_modes_um(i,k,mode_ait_insol-1)
        end do
      end do
    end if

    !------------------------------------------------
    ! Now calculate aod and aaod for Accumulation Insoluble mode

    if ( .not. associated( aod_ukca_acc_ins, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, row_length
          aod_ukca_acc_ins( map_aod_wavel(i) + k - 1 ) =                       &
                                    aod_ukca_all_modes_um(i,k,mode_acc_insol-1)
        end do
      end do
    end if

    if ( .not. associated( aaod_ukca_acc_ins, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, row_length
          aaod_ukca_acc_ins( map_aod_wavel(i) + k - 1 ) =                      &
                                   aaod_ukca_all_modes_um(i,k,mode_acc_insol-1)
        end do
      end do
    end if

    !------------------------------------------------
    ! Now calculate aod and aaod for Coarse Insoluble mode

    if ( .not. associated( aod_ukca_cor_ins, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, row_length
          aod_ukca_cor_ins( map_aod_wavel(i) + k - 1 ) =                       &
                                    aod_ukca_all_modes_um(i,k,mode_cor_insol-1)
        end do
      end do
    end if

    if ( .not. associated( aaod_ukca_cor_ins, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, row_length
           aaod_ukca_cor_ins( map_aod_wavel(i) + k - 1 ) =                     &
                                   aaod_ukca_all_modes_um(i,k,mode_cor_insol-1)
        end do
      end do
    end if

  !------------------------------------------------
  ! This is where we would close segmentation and openmp
  !------------------------------------------------

  deallocate( ukca_modal_wtv_um )
  deallocate( ukca_modal_vol_um )
  deallocate( ukca_modal_rho_um )
  deallocate( ukca_modal_nbr_um )
  deallocate( ukca_wet_diam_um )
  deallocate( ukca_dry_diam_um )

  deallocate( ukca_mix_ratio_um )
  deallocate( ukca_comp_vol_um )

  deallocate( aaod_ukca_all_modes_um )
  deallocate(  sod_ukca_all_modes_um )
  deallocate(  aod_ukca_all_modes_um )

  deallocate( aer_sw_asymmetry_um )
  deallocate( aer_sw_scattering_um )
  deallocate( aer_sw_absorption_um )
  deallocate( aer_lw_asymmetry_um )
  deallocate( aer_lw_scattering_um )
  deallocate( aer_lw_absorption_um )

  deallocate( ukca_mode_mix_ratio_um )

end subroutine radaer_code

end module radaer_kernel_mod
