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
                             DOMAIN, GH_INTEGER,                               &
                             ANY_DISCONTINUOUS_SPACE_1,                        &
                             ANY_DISCONTINUOUS_SPACE_2,                        &
                             ANY_DISCONTINUOUS_SPACE_3,                        &
                             ANY_DISCONTINUOUS_SPACE_4,                        &
                             ANY_DISCONTINUOUS_SPACE_5

use empty_data_mod,    only: empty_real_data

use fs_continuity_mod, only: WTHETA, W3

use kernel_mod,        only: kernel_type

implicit none

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
  integer :: operates_on = DOMAIN
contains
  procedure, nopass :: radaer_code
end type

public :: radaer_code

contains

!> @brief Interface to glomap aerosol climatology scheme.
!> @param[in]     nlayers            The number of layers
!> @param[in]     ncells             Number of horizontal cells in segment  
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
                        ncells,                                                &
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
                                               sw_index_exclude,               &
                                               n_lw_band,                      &
                                               lw_index_exclude

  use um_physics_init_mod,               only: n_radaer_mode,                  &
                                               n_aer_mode_sw, n_aer_mode_lw

  use ukca_mode_setup,                   only: nmodes, ncp_max,                &
                                               mode_nuc_sol,                   &
                                               mode_ait_sol, mode_acc_sol,     &
                                               mode_cor_sol, mode_ait_insol,   &
                                               mode_acc_insol, mode_cor_insol, &
                                               cp_su,  cp_bc, cp_oc,           &
                                               cp_cl,  cp_du, cp_so,           &
                                               cp_no3, cp_nn, cp_nh4,          &
                                               i_ukca_bc_tuned

  use planet_config_mod,                 only: p_zero, kappa, gravity

  use ukca_radaer_precalc,               only: npd_ukca_aod_wavel

  use ukca_option_mod,                   only: do_not_prescribe

  use ukca_radaer_lfric_mod,             only: ukca_radaer_lfric_interface

  implicit none

  ! Arguments

  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ncells

  integer(kind=i_def), intent(in) :: ndf_wth
  integer(kind=i_def), intent(in) :: undf_wth
  integer(kind=i_def), dimension(ndf_wth, ncells), intent(in) :: map_wth

  integer(kind=i_def), intent(in) :: ndf_w3
  integer(kind=i_def), intent(in) :: undf_w3
  integer(kind=i_def), dimension(ndf_w3, ncells), intent(in) :: map_w3

  integer(kind=i_def), intent(in) :: ndf_2d
  integer(kind=i_def), intent(in) :: undf_2d
  integer(kind=i_def), dimension(ndf_2d, ncells), intent(in) :: map_2d

  integer(kind=i_def), intent(in) :: ndf_mode
  integer(kind=i_def), intent(in) :: undf_mode
  integer(kind=i_def), dimension(ndf_mode, ncells), intent(in) :: map_mode

  integer(kind=i_def), intent(in) :: ndf_rmode_sw
  integer(kind=i_def), intent(in) :: undf_rmode_sw
  integer(kind=i_def), dimension(ndf_rmode_sw, ncells), intent(in) ::          &
                                                                   map_rmode_sw

  integer(kind=i_def), intent(in) :: ndf_rmode_lw
  integer(kind=i_def), intent(in) :: undf_rmode_lw
  integer(kind=i_def), dimension(ndf_rmode_lw, ncells), intent(in) ::          &
                                                                   map_rmode_lw

  integer(kind=i_def), intent(in) :: ndf_aod_wavel
  integer(kind=i_def), intent(in) :: undf_aod_wavel
  integer(kind=i_def), dimension(ndf_aod_wavel, ncells), intent(in) ::         &
                                                                  map_aod_wavel

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
  integer(i_um), parameter :: ip_solar = 1
  integer(i_um), parameter :: ip_infra_red = 2

  ! Prescribed single-scattering albedo dummy variables
  ! Make these namelist options later
  integer, parameter       :: i_ukca_radaer_prescribe_ssa = do_not_prescribe
  integer(i_um), parameter :: nd_prof_ssa = 1
  integer(i_um), parameter :: nd_layr_ssa = 1
  integer(i_um), parameter :: nd_band_ssa = 1
  real(r_um),dimension( nd_prof_ssa, nd_layr_ssa, nd_band_ssa ) ::             &
                                                          ukca_radaer_presc_ssa

  ! Loop counters
  integer(i_um) :: k, i, ii, i_band, i_mode, i_rmode , j, jj

  integer(i_um), dimension(:), allocatable :: trindxrad_um

  ! pressure on theta levels
  real(r_um), dimension(:,:), allocatable :: p_theta_levels_um

  ! temperature on theta levels
  real(r_um), dimension(:,:), allocatable :: t_theta_levels_um

  ! d_mass on theta levels
  real(r_um), dimension(:,:), allocatable :: d_mass_theta_levels_um

  real(r_um), dimension(:,:,:), allocatable :: ukca_comp_vol_um
  real(r_um), dimension(:,:,:), allocatable :: ukca_mix_ratio_um
  real(r_um), dimension(:,:,:), allocatable :: ukca_dry_diam_um
  real(r_um), dimension(:,:,:), allocatable :: ukca_wet_diam_um
  real(r_um), dimension(:,:,:), allocatable :: ukca_modal_nbr_um
  real(r_um), dimension(:,:,:), allocatable :: ukca_modal_rho_um
  real(r_um), dimension(:,:,:), allocatable :: ukca_modal_vol_um
  real(r_um), dimension(:,:,:), allocatable :: ukca_modal_wtv_um
  real(r_um), dimension(:,:,:), allocatable :: ukca_mode_mix_ratio_um

  real(r_um), dimension(:,:,:,:), allocatable :: aer_lw_absorption_um
  real(r_um), dimension(:,:,:,:), allocatable :: aer_lw_scattering_um
  real(r_um), dimension(:,:,:,:), allocatable :: aer_lw_asymmetry_um
  real(r_um), dimension(:,:,:,:), allocatable :: aer_sw_absorption_um
  real(r_um), dimension(:,:,:,:), allocatable :: aer_sw_scattering_um
  real(r_um), dimension(:,:,:,:), allocatable :: aer_sw_asymmetry_um

  ! UKCA modal optical depth diagnostics: full domain
  real(r_um), dimension(:,:,:), allocatable :: aod_ukca_all_modes_um

  ! Not yet included as diagnostic
  ! UKCA modal optical depth diagnostics: stratosphere
  real(r_um), dimension(:,:,:), allocatable :: sod_ukca_all_modes_um

  ! UKCA modal absorption optical depth diagnostics: full column
  real(r_um), dimension(:,:,:), allocatable :: aaod_ukca_all_modes_um

  integer, parameter :: i_ukca_tune_bc = i_ukca_bc_tuned
  integer, parameter :: i_glomap_clim_tune_bc = i_ukca_bc_tuned
  logical, parameter :: l_nitrate = .false. ! Make this a namelist option later
  logical, parameter :: l_sustrat = .true.  ! Make this a namelist option later
                                            ! l_sustrat=.true. for ga9
  logical, parameter :: l_cornarrow_ins = .false.
                                            ! Make this a namelist option later

  integer(i_um) :: ncp_max_x_nmodes
  integer(i_um) :: i_cpnt_index( ncp_max, nmodes )
  integer(i_um) :: i_cpnt_type( ncp_max * nmodes )
  integer(i_um) :: i_mode_type( nmodes )
  integer(i_um) :: n_cpnt_in_mode( nmodes )
  logical       :: l_soluble( nmodes )

  ! By convention, arrays are inverted in UM radiation code
  ! Since we are calling from LFRic, arrays will not be inverted
  ! This matters for determining whether a level is above the tropopause
  logical, parameter :: l_inverted = .false.

  ! Whether we need to run shortwave band_average because lit or not
  logical :: l_any_lit_points_um

  ! Variables close to but not exactly 1 or -1 for bounding asymmetry
  real(r_def), parameter :: one_minus_eps = 1.0_r_def - epsilon(1.0_r_def)
  real(r_def), parameter :: minus1_plus_eps = -1.0_r_def + epsilon(1.0_r_def)

  !-----------------------------------------------------------------------

  logical, parameter :: soluble_wanted   = .true.
  logical, parameter :: soluble_unwanted = .false.

  !-----------------------------------------------------------------------
  ! Segmetation

  integer(i_um) :: nml_segment_size , num_seg , seg_size

  !
  !
  ! Hard coded here, but need to make this a namelist item later
  !
  !
  integer, parameter :: ukca_radaer_segment_size = 4

  !-----------------------------------------------------------------------

  logical :: l_aod_ukca_ait_sol, l_aaod_ukca_ait_sol,                          &
             l_aod_ukca_acc_sol, l_aaod_ukca_acc_sol,                          &
             l_aod_ukca_cor_sol, l_aaod_ukca_cor_sol,                          &
             l_aod_ukca_ait_ins, l_aaod_ukca_ait_ins,                          &
             l_aod_ukca_acc_ins, l_aaod_ukca_acc_ins,                          &
             l_aod_ukca_cor_ins, l_aaod_ukca_cor_ins

  !-----------------------------------------------------------------------

  ncp_max_x_nmodes = ncp_max * nmodes

  !-----------------------------------------------------------------------
  ! Populate ukca_radaer element arrays
  ! Note that nucleation mode gets ignored in some of these
  !-----------------------------------------------------------------------

  ! No nucleation mode
  l_soluble(1:nmodes) =  (/.true., .true., .true., .false.,                    &
                           .false.,.false.,.false.,.false./)

  ! No nucleation mode
  n_cpnt_in_mode(1:nmodes) = (/ 3, 5, 5, 2, 1, 1, -1, -1 /)

  ! No nucleation mode
  i_mode_type(1:nmodes)    = (/ 1, 2, 3, 1, 2, 3, -1, -1 /)

  ! No nucleation mode
  i_cpnt_index(cp_su, 1:nmodes)=(/  1,  4,  9, 14, 16, 17, -1, -1 /)
  i_cpnt_index(cp_bc, 1:nmodes)=(/  2,  5, 10, 15, -1, -1, -1, -1 /)
  i_cpnt_index(cp_oc, 1:nmodes)=(/  3,  6, 11, -1, -1, -1, -1, -1 /)
  i_cpnt_index(cp_cl, 1:nmodes)=(/ -1,  7, 12, -1, -1, -1, -1, -1 /)
  i_cpnt_index(cp_du, 1:nmodes)=(/ -1,  8, 13, -1, -1, -1, -1, -1 /)
  i_cpnt_index(cp_so, 1:nmodes)=(/ -1, -1, -1, -1, -1, -1, -1, -1 /)
  i_cpnt_index(cp_no3,1:nmodes)=(/ -1, -1, -1, -1, -1, -1, -1, -1 /)
  i_cpnt_index(cp_nn, 1:nmodes)=(/ -1, -1, -1, -1, -1, -1, -1, -1 /)
  i_cpnt_index(cp_nh4,1:nmodes)=(/ -1, -1, -1, -1, -1, -1, -1, -1 /)

  i_cpnt_type(1:ncp_max_x_nmodes) = (/ 1,  2,  3,  1,  2,  3,  4,  5,  1,      &
                                       2,  3,  4,  5,  2,  3,  5,  5, -1,      &
                                      -1, -1, -1, -1, -1, -1, -1, -1, -1,      &
                                      -1, -1, -1, -1, -1, -1, -1, -1, -1,      &
                                      -1, -1, -1, -1, -1, -1, -1, -1, -1,      &
                                      -1, -1, -1, -1, -1, -1, -1, -1, -1,      &
                                      -1, -1, -1, -1, -1, -1, -1, -1, -1,      &
                                      -1, -1, -1, -1, -1, -1, -1, -1, -1 /)

  !-----------------------------------------------------------------------

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
  ! Start of segmetation

  nml_segment_size = min( ukca_radaer_segment_size, ncells )

  num_seg = ceiling( real(ncells,r_um) / real(nml_segment_size,r_um) )

  allocate( trindxrad_um( nml_segment_size ) )

  allocate( p_theta_levels_um( nml_segment_size, nlayers ) )
  allocate( t_theta_levels_um( nml_segment_size, nlayers ) )
  allocate( d_mass_theta_levels_um( nml_segment_size, nlayers ) )

  allocate( ukca_comp_vol_um( n_ukca_cpnt, nml_segment_size, nlayers ) )
  allocate( ukca_mix_ratio_um( n_ukca_cpnt, nml_segment_size, nlayers ) )

  allocate( ukca_dry_diam_um( nml_segment_size, nlayers, n_ukca_mode ) )
  allocate( ukca_wet_diam_um( nml_segment_size, nlayers, n_ukca_mode ) )
  allocate( ukca_modal_nbr_um( nml_segment_size, nlayers, n_ukca_mode ) )
  allocate( ukca_modal_rho_um( nml_segment_size, nlayers, n_ukca_mode ) )
  allocate( ukca_modal_vol_um( nml_segment_size, nlayers, n_ukca_mode ) )
  allocate( ukca_modal_wtv_um( nml_segment_size, nlayers, n_ukca_mode ) )

  allocate( ukca_mode_mix_ratio_um( nml_segment_size, nlayers, n_radaer_mode ) )

  allocate( aer_lw_absorption_um( nml_segment_size, nlayers, n_radaer_mode,    &
                                  n_lw_band ) )

  allocate( aer_lw_scattering_um( nml_segment_size, nlayers, n_radaer_mode,    &
                                  n_lw_band ) )

  allocate( aer_lw_asymmetry_um( nml_segment_size, nlayers, n_radaer_mode,     &
                                 n_lw_band ) )

  allocate( aer_sw_absorption_um( nml_segment_size, nlayers, n_radaer_mode,    &
                                  n_sw_band ) )

  allocate( aer_sw_scattering_um( nml_segment_size, nlayers, n_radaer_mode,    &
                                  n_sw_band ) )

  allocate( aer_sw_asymmetry_um( nml_segment_size, nlayers, n_radaer_mode,     &
                                 n_sw_band ) )

  allocate( aod_ukca_all_modes_um( nml_segment_size, npd_ukca_aod_wavel,       &
                                   n_ukca_mode ) )

  allocate( sod_ukca_all_modes_um( nml_segment_size, npd_ukca_aod_wavel,       &
                                   n_ukca_mode ) )

  allocate( aaod_ukca_all_modes_um(nml_segment_size, npd_ukca_aod_wavel,       &
                                   n_ukca_mode ) )

  !-----------------------------------------------------------------------

  seg_size = nml_segment_size

  npd_exclude_lw = SIZE( lw_index_exclude, 1 )
  npd_exclude_sw = SIZE( sw_index_exclude, 1 )

  !-----------------------------------------------------------------------

  ! The following appear to be parameters so dont need to be shared
  !  npd_ukca_aod_wavel, n_ukca_mode, n_ukca_cpnt,
  !  nd_prof_ssa, nd_layr_ssa, nd_band_ssa, nmodes, ncp_max,
  !  i_ukca_tune_bc, i_ukca_radaer_prescribe_ssa, i_glomap_clim_tune_bc,
  !  ip_infra_red, ip_solar,
  !  soluble_wanted, soluble_unwanted,
  !  l_exclude_lw, l_exclude_sw, l_nitrate, l_sustrat,
  !  l_cornarrow_ins, l_inverted,
  !  mode_ait_sol, mode_acc_sol, mode_cor_sol,
  !  mode_ait_insol, mode_acc_insol, mode_cor_insol,
  !  minus1_plus_eps, one_minus_eps,

!$OMP PARALLEL DEFAULT(none)                                                   &
!$OMP private( i, ii, j, jj, k, i_mode, i_rmode, i_band,                       &
!$OMP          l_any_lit_points_um, trindxrad_um,                              &
!$OMP          p_theta_levels_um, t_theta_levels_um, ukca_comp_vol_um,         &
!$OMP          ukca_mix_ratio_um, ukca_dry_diam_um, ukca_modal_nbr_um,         &
!$OMP          ukca_modal_rho_um, ukca_modal_vol_um, ukca_modal_wtv_um,        &
!$OMP          ukca_wet_diam_um, d_mass_theta_levels_um,                       &
!$OMP          ukca_mode_mix_ratio_um,                                         &
!$OMP          aer_lw_absorption_um, aer_lw_scattering_um, aer_lw_asymmetry_um,&
!$OMP          aer_sw_absorption_um, aer_sw_scattering_um, aer_sw_asymmetry_um,&
!$OMP          aod_ukca_all_modes_um, aaod_ukca_all_modes_um )                 &
!$OMP shared( num_seg, seg_size, n_radaer_step, nlayers, npd_exclude_lw,       &
!$OMP         npd_exclude_sw,                                                  &
!$OMP         ncp_max_x_nmodes, n_cpnt_in_mode, n_radaer_mode, n_lw_band,      &
!$OMP         n_sw_band, n_aer_mode_lw, n_aer_mode_sw,                         &
!$OMP         i_cpnt_index, i_cpnt_type, i_mode_type,                          &
!$OMP         trop_level, p_zero, kappa, lit_fraction, gravity,                &
!$OMP         theta_in_wth, exner_in_wth, ukca_radaer_presc_ssa,               &
!$OMP         map_2d, map_wth, map_w3, map_mode, map_rmode_lw, map_rmode_sw,   &
!$OMP         map_aod_wavel,                                                   &
!$OMP         l_soluble,                                                       &
!$OMP         l_aod_ukca_ait_sol, l_aod_ukca_acc_sol, l_aod_ukca_cor_sol,      &
!$OMP         l_aaod_ukca_ait_sol, l_aaod_ukca_acc_sol, l_aaod_ukca_cor_sol,   &
!$OMP         l_aod_ukca_ait_ins, l_aod_ukca_acc_ins, l_aod_ukca_cor_ins,      &
!$OMP         l_aaod_ukca_ait_ins, l_aaod_ukca_acc_ins, l_aaod_ukca_cor_ins,   &
!$OMP         pvol_wat_ait_sol, pvol_wat_acc_sol, pvol_wat_cor_sol,            &
!$OMP         pvol_su_ait_sol, pvol_bc_ait_sol, pvol_om_ait_sol,               &
!$OMP         pvol_su_acc_sol, pvol_bc_acc_sol, pvol_om_acc_sol,               &
!$OMP         pvol_ss_acc_sol,                                                 &
!$OMP         pvol_su_cor_sol, pvol_bc_cor_sol, pvol_om_cor_sol,               &
!$OMP         pvol_ss_cor_sol,                                                 &
!$OMP         pvol_bc_ait_ins, pvol_om_ait_ins,                                &
!$OMP         pvol_du_acc_ins, pvol_du_cor_ins,                                &
!$OMP         ait_sol_su, ait_sol_bc, ait_sol_om,                              &
!$OMP         acc_sol_su, acc_sol_bc, acc_sol_om, acc_sol_ss,                  &
!$OMP         cor_sol_su, cor_sol_bc, cor_sol_om, cor_sol_ss,                  &
!$OMP         ait_ins_bc, ait_ins_om, acc_ins_du, cor_ins_du,                  &
!$OMP         n_ait_sol, n_acc_sol, n_cor_sol, n_ait_ins, n_acc_ins, n_cor_ins,&
!$OMP         rhopar_ait_sol, rhopar_acc_sol, rhopar_cor_sol,                  &
!$OMP         rhopar_ait_ins, rhopar_acc_ins, rhopar_cor_ins,                  &
!$OMP         drydp_ait_sol, drydp_acc_sol, drydp_cor_sol,                     &
!$OMP         drydp_ait_ins, drydp_acc_ins, drydp_cor_ins,                     &
!$OMP         wetdp_ait_sol, wetdp_acc_sol, wetdp_cor_sol,                     &
!$OMP         aod_ukca_ait_sol, aod_ukca_acc_sol, aod_ukca_cor_sol,            &
!$OMP         aaod_ukca_ait_sol, aaod_ukca_acc_sol, aaod_ukca_cor_sol,         &
!$OMP         aod_ukca_ait_ins, aod_ukca_acc_ins, aod_ukca_cor_ins,            &
!$OMP         aaod_ukca_ait_ins, aaod_ukca_acc_ins, aaod_ukca_cor_ins,         &
!$OMP         rho_in_wth, dz_in_wth, exner_in_w3, aer_mix_ratio,               &
!$OMP         aer_lw_absorption, aer_lw_scattering, aer_lw_asymmetry,          &
!$OMP         aer_sw_absorption, aer_sw_scattering, aer_sw_asymmetry,          &
!$OMP         empty_real_data )

!$OMP do SCHEDULE(DYNAMIC)
  do j = 1, num_seg

    jj = ( (j-1)*seg_size ) + 1



    ! Whether we need to run shortwave band_average because lit or not
    l_any_lit_points_um = .false.
    if ( n_radaer_step > 1 ) then
      l_any_lit_points_um = .true.
    else
      do i = 1, seg_size
        ii = jj + i - 1
        if ( lit_fraction( map_2d(1, ii) ) > 0.0_r_def ) then
          l_any_lit_points_um = .true.
        end if   
      end do
    end if

    ! Note that this is inverted compared to the UM
    ! This will be dealt with in ukca_radaer_band_average
    do i = 1, seg_size
      ii = jj + i - 1
      trindxrad_um(i) = trop_level( map_2d(1, ii) )
    end do

    !-----------------------------------------------------------------------
    ! Initialisation of prognostic variables and arrays
    !-----------------------------------------------------------------------

    do i = 1, seg_size
      do k = 1, nlayers
        ii = jj + i - 1
        p_theta_levels_um(i,k) = p_zero *                                      &
                         ( exner_in_wth(map_wth(1, ii) + k) )**(1.0_r_um/kappa)
      end do
    end do

    do i = 1, seg_size
      do k = 1, nlayers
        ii = jj + i - 1
        t_theta_levels_um(i,k) = exner_in_wth(map_wth(1, ii) + k) *            &
                                 theta_in_wth(map_wth(1, ii) + k)
      end do
    end do

    ! note - zeroth level is redundant for these fields in UM
    do i = 1, seg_size
      do k = 1, nlayers

        ii = jj + i - 1

        ukca_comp_vol_um(1, i,k) = pvol_su_ait_sol(map_wth(1, ii) + k)
        ukca_comp_vol_um(2, i,k) = pvol_bc_ait_sol(map_wth(1, ii) + k)
        ukca_comp_vol_um(3, i,k) = pvol_om_ait_sol(map_wth(1, ii) + k)
        ukca_comp_vol_um(4, i,k) = pvol_su_acc_sol(map_wth(1, ii) + k)
        ukca_comp_vol_um(5, i,k) = pvol_bc_acc_sol(map_wth(1, ii) + k)
        ukca_comp_vol_um(6, i,k) = pvol_om_acc_sol(map_wth(1, ii) + k)
        ukca_comp_vol_um(7, i,k) = pvol_ss_acc_sol(map_wth(1, ii) + k)
        ukca_comp_vol_um(8, i,k) = 0.0_r_um ! no pvol_du_acc_sol prognostic
        ukca_comp_vol_um(9, i,k) = pvol_su_cor_sol(map_wth(1, ii) + k)
        ukca_comp_vol_um(10,i,k) = pvol_bc_cor_sol(map_wth(1, ii) + k)
        ukca_comp_vol_um(11,i,k) = pvol_om_cor_sol(map_wth(1, ii) + k)
        ukca_comp_vol_um(12,i,k) = pvol_ss_cor_sol(map_wth(1, ii) + k)
        ukca_comp_vol_um(13,i,k) = 0.0_r_um ! no pvol_du_cor_sol prognostic
        ukca_comp_vol_um(14,i,k) = pvol_bc_ait_ins(map_wth(1, ii) + k)
        ukca_comp_vol_um(15,i,k) = pvol_om_ait_ins(map_wth(1, ii) + k)
        ukca_comp_vol_um(16,i,k) = pvol_du_acc_ins(map_wth(1, ii) + k)
        ukca_comp_vol_um(17,i,k) = pvol_du_cor_ins(map_wth(1, ii) + k)

        ukca_mix_ratio_um(1, i,k) = ait_sol_su(map_wth(1, ii) + k)
        ukca_mix_ratio_um(2, i,k) = ait_sol_bc(map_wth(1, ii) + k)
        ukca_mix_ratio_um(3, i,k) = ait_sol_om(map_wth(1, ii) + k)
        ukca_mix_ratio_um(4, i,k) = acc_sol_su(map_wth(1, ii) + k)
        ukca_mix_ratio_um(5, i,k) = acc_sol_bc(map_wth(1, ii) + k)
        ukca_mix_ratio_um(6, i,k) = acc_sol_om(map_wth(1, ii) + k)
        ukca_mix_ratio_um(7, i,k) = acc_sol_ss(map_wth(1, ii) + k)
        ukca_mix_ratio_um(8, i,k) = 0.0_r_um ! no acc_sol_du prognostic
        ukca_mix_ratio_um(9, i,k) = cor_sol_su(map_wth(1, ii) + k)
        ukca_mix_ratio_um(10,i,k) = cor_sol_bc(map_wth(1, ii) + k)
        ukca_mix_ratio_um(11,i,k) = cor_sol_om(map_wth(1, ii) + k)
        ukca_mix_ratio_um(12,i,k) = cor_sol_ss(map_wth(1, ii) + k)
        ukca_mix_ratio_um(13,i,k) = 0.0_r_um ! no cor_sol_du prognostic
        ukca_mix_ratio_um(14,i,k) = ait_ins_bc(map_wth(1, ii) + k)
        ukca_mix_ratio_um(15,i,k) = ait_ins_om(map_wth(1, ii) + k)
        ukca_mix_ratio_um(16,i,k) = acc_ins_du(map_wth(1, ii) + k)
        ukca_mix_ratio_um(17,i,k) = cor_ins_du(map_wth(1, ii) + k)

        ukca_dry_diam_um(i,k,(mode_ait_sol-1))  =drydp_ait_sol(map_wth(1,ii)+ k)
        ukca_dry_diam_um(i,k,(mode_acc_sol-1))  =drydp_acc_sol(map_wth(1,ii)+ k)
        ukca_dry_diam_um(i,k,(mode_cor_sol-1))  =drydp_cor_sol(map_wth(1,ii)+ k)
        ukca_dry_diam_um(i,k,(mode_ait_insol-1))=drydp_ait_ins(map_wth(1,ii)+ k)
        ukca_dry_diam_um(i,k,(mode_acc_insol-1))=drydp_acc_ins(map_wth(1,ii)+ k)
        ukca_dry_diam_um(i,k,(mode_cor_insol-1))=drydp_cor_ins(map_wth(1,ii)+ k)

        ukca_modal_nbr_um(i,k,(mode_ait_sol-1))  = n_ait_sol(map_wth(1, ii) + k)
        ukca_modal_nbr_um(i,k,(mode_acc_sol-1))  = n_acc_sol(map_wth(1, ii) + k)
        ukca_modal_nbr_um(i,k,(mode_cor_sol-1))  = n_cor_sol(map_wth(1, ii) + k)
        ukca_modal_nbr_um(i,k,(mode_ait_insol-1))= n_ait_ins(map_wth(1, ii) + k)
        ukca_modal_nbr_um(i,k,(mode_acc_insol-1))= n_acc_ins(map_wth(1, ii) + k)
        ukca_modal_nbr_um(i,k,(mode_cor_insol-1))= n_cor_ins(map_wth(1, ii) + k)

        ukca_modal_rho_um(i,k,(mode_ait_sol-1))   =                            &
                                        rhopar_ait_sol(map_wth(1, ii) + k)

        ukca_modal_rho_um(i,k,(mode_acc_sol-1))   =                            &
                                        rhopar_acc_sol(map_wth(1, ii) + k)

        ukca_modal_rho_um(i,k,(mode_cor_sol-1))   =                            &
                                        rhopar_cor_sol(map_wth(1, ii) + k)

        ukca_modal_rho_um(i,k,(mode_ait_insol-1)) =                            &
                                        rhopar_ait_ins(map_wth(1, ii) + k)

        ukca_modal_rho_um(i,k,(mode_acc_insol-1)) =                            &
                                        rhopar_acc_ins(map_wth(1, ii) + k)

        ukca_modal_rho_um(i,k,(mode_cor_insol-1)) =                            &
                                        rhopar_cor_ins(map_wth(1, ii) + k)

        ukca_modal_vol_um(i,k,(mode_ait_sol-1)) =                              &
                                        pvol_wat_ait_sol(map_wth(1, ii) + k) + &
                                        pvol_su_ait_sol( map_wth(1, ii) + k) + &
                                        pvol_bc_ait_sol( map_wth(1, ii) + k) + &
                                        pvol_om_ait_sol( map_wth(1, ii) + k)

        ukca_modal_vol_um(i,k,(mode_acc_sol-1)) =                              &
                                        pvol_wat_acc_sol(map_wth(1, ii) + k) + &
                                        pvol_su_acc_sol( map_wth(1, ii) + k) + &
                                        pvol_bc_acc_sol( map_wth(1, ii) + k) + &
                                        pvol_om_acc_sol( map_wth(1, ii) + k) + &
                                        pvol_ss_acc_sol( map_wth(1, ii) + k)

        ukca_modal_vol_um(i,k,(mode_cor_sol-1)) =                              &
                                        pvol_wat_cor_sol(map_wth(1, ii) + k) + &
                                        pvol_su_cor_sol( map_wth(1, ii) + k) + &
                                        pvol_bc_cor_sol( map_wth(1, ii) + k) + &
                                        pvol_om_cor_sol( map_wth(1, ii) + k) + &
                                        pvol_ss_cor_sol( map_wth(1, ii) + k)

        ukca_modal_vol_um(i,k,(mode_ait_insol-1)) =                            &
                                        pvol_bc_ait_ins( map_wth(1, ii) + k) + &
                                        pvol_om_ait_ins( map_wth(1, ii) + k)

        ukca_modal_vol_um(i,k,(mode_acc_insol-1)) =                            &
                                        pvol_du_acc_ins( map_wth(1, ii) + k)

        ukca_modal_vol_um(i,k,(mode_cor_insol-1)) =                            &
                                        pvol_du_cor_ins( map_wth(1, ii) + k)

        ukca_modal_wtv_um(i,k,(mode_ait_sol-1)) =                              &
                                        pvol_wat_ait_sol(map_wth(1, ii) + k)

        ukca_modal_wtv_um(i,k,(mode_acc_sol-1)) =                              &
                                        pvol_wat_acc_sol(map_wth(1, ii) + k)

        ukca_modal_wtv_um(i,k,(mode_cor_sol-1)) =                              &
                                        pvol_wat_cor_sol(map_wth(1, ii) + k)

        ukca_modal_wtv_um(i,k,(mode_ait_insol-1)) = 0.0_r_um
        ukca_modal_wtv_um(i,k,(mode_acc_insol-1)) = 0.0_r_um
        ukca_modal_wtv_um(i,k,(mode_cor_insol-1)) = 0.0_r_um

        ukca_wet_diam_um(i,k,(mode_ait_sol-1))  =wetdp_ait_sol(map_wth(1,ii)+ k)
        ukca_wet_diam_um(i,k,(mode_acc_sol-1))  =wetdp_acc_sol(map_wth(1,ii)+ k)
        ukca_wet_diam_um(i,k,(mode_cor_sol-1))  =wetdp_cor_sol(map_wth(1,ii)+ k)
        ! Note that wet and dry diameter are the same for insoluble modes
        ukca_wet_diam_um(i,k,(mode_ait_insol-1))=drydp_ait_ins(map_wth(1,ii)+ k)
        ukca_wet_diam_um(i,k,(mode_acc_insol-1))=drydp_acc_ins(map_wth(1,ii)+ k)
        ukca_wet_diam_um(i,k,(mode_cor_insol-1))=drydp_cor_ins(map_wth(1,ii)+ k)
      end do
    end do

    !------------------------------------------------
    ! Calculate mass thickness of vertical levels
    ! This duplicates calculation of d_mass from set_thermodynamic_kernel_mod

    if ( ( .not. associated( aod_ukca_ait_sol,  empty_real_data ) ) .or.       &
         ( .not. associated( aaod_ukca_ait_sol, empty_real_data ) ) .or.       &
         ( .not. associated( aod_ukca_acc_sol,  empty_real_data ) ) .or.       &
         ( .not. associated( aaod_ukca_acc_sol, empty_real_data ) ) .or.       &
         ( .not. associated( aod_ukca_cor_sol,  empty_real_data ) ) .or.       &
         ( .not. associated( aaod_ukca_cor_sol, empty_real_data ) ) .or.       &
         ( .not. associated( aod_ukca_ait_ins,  empty_real_data ) ) .or.       &
         ( .not. associated( aaod_ukca_ait_ins, empty_real_data ) ) .or.       &
         ( .not. associated( aod_ukca_acc_ins,  empty_real_data ) ) .or.       &
         ( .not. associated( aaod_ukca_acc_ins, empty_real_data ) ) .or.       &
         ( .not. associated( aod_ukca_cor_ins,  empty_real_data ) ) .or.       &
         ( .not. associated( aaod_ukca_cor_ins, empty_real_data ) ) ) then

      do i = 1, seg_size
        ii = jj + i - 1
        d_mass_theta_levels_um(i,1) = rho_in_wth(  map_wth(2, ii) ) *          &
                                      ( dz_in_wth( map_wth(2, ii) ) +          &
                                        dz_in_wth( map_wth(1, ii) ) )
      end do

      do i = 1, seg_size
        do k = 2, nlayers - 1
          ii = jj + i - 1
          d_mass_theta_levels_um(i,k) = rho_in_wth( map_wth(1, ii) + k ) *     &
                                         dz_in_wth( map_wth(1, ii) + k )
        end do
      end do

      do i = 1, seg_size
        ii = jj + i - 1
        d_mass_theta_levels_um(i,nlayers) = p_zero *                           &
             exner_in_w3( map_w3(1, ii) + nlayers-1 )**( 1.0_r_def / kappa ) / &
                                            gravity
      end do

    end if

    !-----------------------------------------------------------------------
    ! Call UKCA modules
    !-----------------------------------------------------------------------

    CALL ukca_radaer_lfric_interface(                                          &
      ! Fixed array dimensions (input)
      seg_size,                                                                &
      nlayers,                                                                 &
      npd_exclude_lw,                                                          &
      npd_exclude_sw,                                                          &
      npd_ukca_aod_wavel,                                                      &
      ! Spectral information (input)
      ip_infra_red,                                                            &
      ip_solar,                                                                &
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
      ! UKCA_RADAER structure (input)
      nmodes,                                                                  &
      ncp_max,                                                                 &
      ncp_max_x_nmodes,                                                        &
      i_cpnt_index,                                                            &
      i_cpnt_type,                                                             &
      i_mode_type,                                                             &
      l_nitrate,                                                               &
      l_soluble,                                                               &
      l_sustrat,                                                               &
      l_cornarrow_ins,                                                         &
      n_cpnt_in_mode,                                                          &
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
      ! Maxwell-Garnett mixing approach logical control switches
      i_ukca_tune_bc, i_glomap_clim_tune_bc,                                   &
      ! Type selection
      soluble_wanted,                                                          &
      soluble_unwanted,                                                        &
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
        do i = 1, seg_size
          ii = jj + i - 1
          aer_mix_ratio( map_mode(1,ii) + ( (i_mode-1)*(nlayers+1) ) + k ) =   &
                                         ukca_mode_mix_ratio_um( i, k, i_mode )
        end do
      end do
    end do

    ! Socrates arrays filled with MODE aerosol optical properties in bands
    i_rmode = 0
    do i_band = 1, n_lw_band
      ! Fill the radaer modes within this band
      do i_mode = 1, n_radaer_mode
        i_rmode = i_rmode + 1
        do k = 1, nlayers
          do i = 1, seg_size
            ii = jj + i - 1
            aer_lw_absorption( map_rmode_lw(1, ii) +                           &
                               ( (i_rmode-1)*(nlayers+1) ) + k ) =             &
                               aer_lw_absorption_um( i, k, i_mode, i_band )

            aer_lw_scattering( map_rmode_lw(1, ii) +                           &
                               ( (i_rmode-1)*(nlayers+1) ) + k ) =             &
                               aer_lw_scattering_um( i, k, i_mode, i_band )

            aer_lw_asymmetry( map_rmode_lw(1, ii) +                            &
                              ( (i_rmode-1)*(nlayers+1) ) + k ) =              &
                              max(minus1_plus_eps, min(one_minus_eps,          &
                              aer_lw_asymmetry_um(  i, k, i_mode, i_band ) ) )

          end do ! seg_size
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
          do i = 1, seg_size
            ii = jj + i - 1
            ! lit points
            if ( ( lit_fraction(map_2d(1, ii)) > 0.0_r_def  ) .or.             &
                 n_radaer_step > 1 ) then

              aer_sw_absorption( map_rmode_sw(1, ii) +                         &
                                 ( (i_rmode-1)*(nlayers+1) ) + k ) =           &
                                 aer_sw_absorption_um( i, k, i_mode, i_band )

              aer_sw_scattering( map_rmode_sw(1, ii) +                         &
                                 ( (i_rmode-1)*(nlayers+1) ) + k ) =           &
                                 aer_sw_scattering_um( i, k, i_mode, i_band )

              aer_sw_asymmetry( map_rmode_sw(1, ii) +                          &
                                ( (i_rmode-1)*(nlayers+1) ) + k ) =            &
                                max(minus1_plus_eps, min(one_minus_eps,        &
                                aer_sw_asymmetry_um( i, k, i_mode, i_band ) ) )

            ! unlit points
            else

              aer_sw_absorption( map_rmode_sw(1, ii) +                         &
                   ( (i_rmode-1)*(nlayers+1) ) + k ) =  1.0_r_def

              aer_sw_scattering( map_rmode_sw(1, ii) +                         &
                   ( (i_rmode-1)*(nlayers+1) ) + k ) =  1.0_r_def

              aer_sw_asymmetry( map_rmode_sw(1, ii) +                          &
                   ( (i_rmode-1)*(nlayers+1) ) + k ) = one_minus_eps

            end if
          end do ! seg_size
        end do ! nlayers
      end do ! n_radaer_mode

      ! If there are additional aerosol modes not associated with radaer
      ! (e.g. from easyaerosol) then i_rmode needs advancing past them
      ! before starting on the next radiation band.
      if (n_aer_mode_sw > n_radaer_mode) then
        i_rmode = i_rmode + n_aer_mode_sw - n_radaer_mode
      end if

    end do ! n_sw_band

    !------------------------------------------------
    ! Now calculate aod and aaod for Aitken Soluble mode

    if ( .not. associated( aod_ukca_ait_sol, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, seg_size
          ii = jj + i - 1
          aod_ukca_ait_sol( map_aod_wavel(1, ii) + k - 1 ) =                   &
                                      aod_ukca_all_modes_um(i,k,mode_ait_sol-1)
        end do
      end do
    end if

    if ( .not. associated( aaod_ukca_ait_sol, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, seg_size
          ii = jj + i - 1
          aaod_ukca_ait_sol( map_aod_wavel(1, ii) + k - 1 ) =                  &
                                     aaod_ukca_all_modes_um(i,k,mode_ait_sol-1)
        end do
      end do
    end if

    !------------------------------------------------
    ! Now calculate aod and aaod for Accumulation Soluble mode

    if ( .not. associated( aod_ukca_acc_sol, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, seg_size
          ii = jj + i - 1
          aod_ukca_acc_sol( map_aod_wavel(1, ii) + k - 1 ) =                   &
                                      aod_ukca_all_modes_um(i,k,mode_acc_sol-1)
        end do
      end do
    end if

    if ( .not. associated( aaod_ukca_acc_sol, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, seg_size
          ii = jj + i - 1
          aaod_ukca_acc_sol( map_aod_wavel(1, ii) + k - 1 ) =                  &
                                     aaod_ukca_all_modes_um(i,k,mode_acc_sol-1)
        end do
      end do
    end if

    !------------------------------------------------
    ! Now calculate aod and aaod for Coarse Soluble mode

    if ( .not. associated( aod_ukca_cor_sol, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, seg_size
          ii = jj + i - 1
          aod_ukca_cor_sol( map_aod_wavel(1, ii) + k - 1 ) =                   &
                                      aod_ukca_all_modes_um(i,k,mode_cor_sol-1)
        end do
      end do
    end if

    if ( .not. associated( aaod_ukca_cor_sol, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, seg_size
          ii = jj + i - 1
          aaod_ukca_cor_sol( map_aod_wavel(1, ii) + k - 1 ) =                  &
                                     aaod_ukca_all_modes_um(i,k,mode_cor_sol-1)
        end do
      end do
    end if

    !------------------------------------------------
    ! Now calculate aod and aaod for Aitken Insoluble mode

    if ( .not. associated( aod_ukca_ait_ins, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, seg_size
          ii = jj + i - 1
          aod_ukca_ait_ins( map_aod_wavel(1, ii) + k - 1 ) =                   &
                                    aod_ukca_all_modes_um(i,k,mode_ait_insol-1)
        end do
      end do
    end if

    if ( .not. associated( aaod_ukca_ait_ins, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, seg_size
          ii = jj + i - 1
          aaod_ukca_ait_ins( map_aod_wavel(1, ii) + k - 1 ) =                  &
                                   aaod_ukca_all_modes_um(i,k,mode_ait_insol-1)
        end do
      end do
    end if

    !------------------------------------------------
    ! Now calculate aod and aaod for Accumulation Insoluble mode

    if ( .not. associated( aod_ukca_acc_ins, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, seg_size
          ii = jj + i - 1
          aod_ukca_acc_ins( map_aod_wavel(1, ii) + k - 1 ) =                   &
                                    aod_ukca_all_modes_um(i,k,mode_acc_insol-1)
        end do
      end do
    end if

    if ( .not. associated( aaod_ukca_acc_ins, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, seg_size
          ii = jj + i - 1
          aaod_ukca_acc_ins( map_aod_wavel(1, ii) + k - 1 ) =                  &
                                   aaod_ukca_all_modes_um(i,k,mode_acc_insol-1)
        end do
      end do
    end if

    !------------------------------------------------
    ! Now calculate aod and aaod for Coarse Insoluble mode

    if ( .not. associated( aod_ukca_cor_ins, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, seg_size
          ii = jj + i - 1
          aod_ukca_cor_ins( map_aod_wavel(1, ii) + k - 1 ) =                   &
                                    aod_ukca_all_modes_um(i,k,mode_cor_insol-1)
        end do
      end do
    end if

    if ( .not. associated( aaod_ukca_cor_ins, empty_real_data ) ) then
      do k = 1, npd_ukca_aod_wavel
        do i = 1, seg_size
           ii = jj + i - 1
           aaod_ukca_cor_ins( map_aod_wavel(1, ii) + k - 1 ) =                 &
                                   aaod_ukca_all_modes_um(i,k,mode_cor_insol-1)
        end do
      end do
    end if

    !------------------------------------------------

  end do ! Segmentation loop
!$OMP end do

!$OMP end PARALLEL  

  !------------------------------------------------

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

  deallocate( ukca_modal_wtv_um )
  deallocate( ukca_modal_vol_um )
  deallocate( ukca_modal_rho_um )
  deallocate( ukca_modal_nbr_um )
  deallocate( ukca_wet_diam_um )
  deallocate( ukca_dry_diam_um )

  deallocate( ukca_mix_ratio_um )
  deallocate( ukca_comp_vol_um )

  deallocate( d_mass_theta_levels_um )
  deallocate( t_theta_levels_um )
  deallocate( p_theta_levels_um )

  deallocate( trindxrad_um )

end subroutine radaer_code

end module radaer_kernel_mod
