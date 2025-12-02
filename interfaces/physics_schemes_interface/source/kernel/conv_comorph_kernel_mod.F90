!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the Comorph Convection scheme.
!>
module conv_comorph_kernel_mod

  use argument_mod,            only : arg_type,                  &
                                      GH_FIELD, GH_SCALAR,       &
                                      GH_INTEGER, GH_REAL,       &
                                      GH_READ, GH_WRITE,         &
                                      GH_READWRITE, DOMAIN,      &
                                      ANY_DISCONTINUOUS_SPACE_1, &
                                      ANY_DISCONTINUOUS_SPACE_2, &
                                      ANY_DISCONTINUOUS_SPACE_3
  use constants_mod,           only : i_def, i_um, r_def, r_um
  use empty_data_mod,          only : empty_real_data
  use fs_continuity_mod,       only : W3, Wtheta
  use kernel_mod,              only : kernel_type
  use timestepping_config_mod, only : outer_iterations
  use microphysics_config_mod, only : prog_tnuc, microphysics_casim

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: conv_comorph_kernel_type
    private
    type(arg_type) :: meta_args(195) = (/                                         &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                &! outer
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! rho_in_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! rho_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! wetrho_in_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! wetrho_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! exner_in_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! exner_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! delta
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! theta_n
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! u_in_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! v_in_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! w_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! theta_star
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! u_in_w3_star
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! v_in_w3_star
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! height_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! height_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dt_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dmv_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dmcl_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dms_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! du_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! dv_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! conv_prog_dtheta
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! conv_prog_dmv
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! mv_n
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! mcl_n
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! mci_n
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! mr_n
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! mg_n
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! ms_n
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_v
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_cl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_ci
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! m_r
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! m_g
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! m_s
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! cf_ice
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! cf_liq
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! cf_bulk
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cca
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! ccw
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! precfrac
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! cf_liq_n
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! cf_fro_n
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! cf_bulk_n
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! cape_diluted
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! cca_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! dd_mf_cb
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! massflux_up
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! massflux_down
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! tke_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     WTHETA),                   &! pressure_inc_env
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! zh_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! zh_nonloc
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! inv_depth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! zhsc_2d
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! bl_type
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! wvar
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! rhokm_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! moist_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! heat_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! taux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! tauy
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_3),&! surf_interp
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! ustar
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! ls_rain_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! ls_snow_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dcfl_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dcff_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dbcf_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! o3p
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! o1d
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! o3
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! nit
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! no
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! no3
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! lumped_n
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! n2o5
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! ho2no2
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! hono2
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! h2o2
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! ch4
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! co
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! hcho
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! meoo
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! meooh
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! h
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! oh
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! ho2
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! cl
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! cl2o2
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! clo
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! oclo
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! br
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! lumped_br
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! brcl
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! brono2
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! n2o
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! lumped_cl
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! hocl
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! hbr
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! hobr
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! clono2
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! cfcl3
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! cf2cl2
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! mebr
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! hono
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! c2h6
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! etoo
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! etooh
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! mecho
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! meco3
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! pan
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! c3h8
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! n_proo
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! i_proo
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! n_prooh
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! i_prooh
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! etcho
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! etco3
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! me2co
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! mecoch2oo
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! mecoch2ooh
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! ppan
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! meono2
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! c5h8
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! iso2
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! isooh
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! ison
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! macr
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! macro2
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! macrooh
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! mpan
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! hacet
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! mgly
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! nald
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! hcooh
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! meco3h
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! meco2h
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! h2
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! meoh
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! msa
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! nh3
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! cs2
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! csul
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! h2s
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! so3
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! passive_o3
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA ),                  & ! age_of_air
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dms
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! so2
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! h2so4
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dmso
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! monoterpene
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! secondary_organic
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! n_nuc_sol
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! nuc_sol_su
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! nuc_sol_om
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! n_ait_sol
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! ait_sol_su
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! ait_sol_bc
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! ait_sol_om
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! n_acc_sol
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! acc_sol_su
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! acc_sol_bc
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! acc_sol_om
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! acc_sol_ss
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! n_cor_sol
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cor_sol_su
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cor_sol_bc
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cor_sol_om
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cor_sol_ss
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! n_ait_ins
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! ait_ins_bc
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! ait_ins_om
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! n_acc_ins
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! acc_ins_du
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! n_cor_ins
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cor_ins_du
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! tnuc
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! lowest_cv_base
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! lowest_cv_top
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! cv_base
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! cv_top
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! pres_cv_base
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! pres_cv_top
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! pres_lowest_cv_base
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! pres_lowest_cv_top
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! lowest_cca_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! entrain_up
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! entrain_down
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! detrain_up
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! detrain_down
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3)                        &! massflux_up_half
        /)
    integer :: operates_on = DOMAIN
  contains
    procedure, nopass :: conv_comorph_code
  end type

  public :: conv_comorph_code

contains

  !> @brief Interface to the Comorph convection scheme
  !> @details The Comorph convection scheme does:
  !>             Vertical non-local transport of heat, momentum, moisture
  !>              and tracers
  !>             Conversion of moisture to precipitation and associated latent
  !>              heat release
  !> @param[in]     nlayers              Number of layers
  !> @param[in]     outer                Outer loop counter
  !> @param[in]     rho_in_w3            Density field in density space
  !> @param[in]     rho_in_wth           Density field in wth space
  !> @param[in]     wetrho_in_w3         Wet density field in density space
  !> @param[in]     wetrho_in_wth        Wet density field in wth space
  !> @param[in]     exner_in_w3          Exner pressure field in density space
  !> @param[in]     exner_in_wth         Exner pressure field in wth space
  !> @param[in]     delta                Edge length on wtheta points
  !> @param[in]     theta_n              Potential temperature at time n
  !> @param[in]     u_in_w3              'Zonal' wind at time n
  !> @param[in]     v_in_w3              'Meridional' wind at time n
  !> @param[in]     w_in_wth             'Vertical' wind in theta space
  !> @param[in]     theta_star           Potential temperature after advection
  !> @param[in]     u_in_w3_star         'Zonal' wind after advection
  !> @param[in]     v_in_w3_star         'Meridional' wind after advection
  !> @param[in]     height_w3            Height of density space above surface
  !> @param[in]     height_wth           Height of theta space above surface
  !> @param[in,out] dt_conv              Convection temperature increment
  !> @param[in,out] dmv_conv             Convection vapour increment
  !> @param[in,out] dmcl_conv            Convection cloud liquid increment
  !> @param[in,out] dms_conv             Convection cloud frozen increment
  !> @param[in,out] du_conv              Convection 'zonal' wind increment
  !> @param[in,out] dv_conv              Convection 'meridional' wind increment
  !> @param[in,out] conv_prog_dtheta     Time smoothed convective theta inc
  !> @param[in,out] conv_prog_dmv        Time smoothed convective humidity inc
  !> @param[in]     m_v                  Vapour mixing ratio after advection
  !> @param[in]     m_cl                 Cloud liq mixing ratio after advection
  !> @param[in]     m_ci                 Cloud ice mixing ratio after advection
  !> @param[in,out] m_r                  Rain mixing ratio after advection
  !> @param[in,out] m_g                  Graupel mixing ratio after advection
  !> @param[in,out] m_s                  Snow mixing ratio after advection
  !> @param[in]     cf_ice               Ice cloud fraction
  !> @param[in]     cf_liq               Liquid cloud fraction
  !> @param[in]     cf_bulk              Bulk cloud fraction
  !> @param[in,out] cca                  Convective cloud amount (fraction)
  !> @param[in,out] ccw                  Convective cloud water (kg/kg) (can be ice or liquid)
  !> @param[in,out] precfrac             3D precipitation fraction
  !> @param[in]     cf_liq_n             Liquid cloud fraction at start of timestep
  !> @param[in]     cf_fro_n             Ice cloud fraction at start of timestep
  !> @param[in]     cf_bulk_n            Bulk cloud fraction at start of timestep
  !> @param[in,out] cape_diluted         CAPE value
  !> @param[in,out] cca_2d               Convective cloud amout (2d) with no anvil
  !> @param[in,out] dd_mf_cb             Downdraft massflux at cloud base (Pa/s)
  !> @param[in,out] massflux_up          Convective upwards mass flux (Pa/s)
  !> @param[in,out] massflux_down        Convective downwards mass flux (Pa/s)
  !> @param[in,out] tke_bl               Turbulent kinetic energy (m2/s2)
  !> @param[out]    pressure_inv_env     Pressure increment to environment from convection
  !> @param[in]     zh_2d                Boundary layer depth
  !> @param[in]     zh_nonloc            Depth of non-local BL scheme
  !> @param[in]     inv_depth            Depth of BL top inversion layer
  !> @param[in]     zhsc_2d              Height of decoupled layer top
  !> @param[in]     bl_type_ind          Diagnosed BL types
  !> @param[in]     wvar                 Vertical velocity variance in wth
  !> @param[in]     rhokm_bl             Momentum eddy diffusivity on BL levels
  !> @param[in]     moist_flux           Vertical moisture flux on BL levels
  !> @param[in]     heat_flux            Vertical heat flux on BL levels
  !> @param[in]     taux                 Explicit u momentum flux at cell centres on BL levels
  !> @param[in]     tauy                 Explicit v momentum flux at cell centres on BL levels
  !> @param[in]     surf_interp          Surface variables for regridding
  !> @param[in]     ustar                Friction velocity
  !> @param[in]     ls_rain_2d           Large scale rain from twod_fields
  !> @param[in]     ls_snow_2d           Large scale snow from twod_fields
  !> @param[in,out] dcfl_conv            Increment to liquid cloud fraction from convection
  !> @param[in,out] dcff_conv            Increment to ice cloud fraction from convection
  !> @param[in,out] dbcf_conv            Increment to bulk cloud fraction from convection
  !> @param[in,out] o3p                  oxygen_ground_state m.m.r
  !> @param[in,out] o1d                  oxygen_excited_state m.m.r
  !> @param[in,out] o3                   ozone m.m.r
  !> @param[in,out] nit                  nitrogen_radical m.m.r
  !> @param[in,out] no                   nitric_oxide m.m.r
  !> @param[in,out] no3                  nitrate_radical m.m.r
  !> @param[in,out] lumped_n             lumped_n_as_nitrogen_dioxide m.m.r
  !> @param[in,out] n2o5                 dinitrogen_pentoxide m.m.r
  !> @param[in,out] ho2no2               peroxynitric_acid m.m.r
  !> @param[in,out] hono2                nitric_acid m.m.r
  !> @param[in,out] h2o2                 Hydrogen peroxide m.m.r.
  !> @param[in,out] ch4                  methane m.m.r
  !> @param[in,out] co                   carbon_monoxide m.m.r
  !> @param[in,out] hcho                 formaldehyde m.m.r
  !> @param[in,out] meoo                 methyl_peroxy_radical m.m.r
  !> @param[in,out] meooh                methyl_hydroperoxide m.m.r
  !> @param[in,out] h                    hydrogen_radical m.m.r
  !> @param[in,out] oh                   hydroxyl_radical m.m.r
  !> @param[in,out] ho2                  hydroxy_peroxyl_radical m.m.r
  !> @param[in,out] cl                   chloride_radical m.m.r
  !> @param[in,out] cl2o2                chlorine_monoxide_dimer m.m.r
  !> @param[in,out] clo                  chlorine_monoxide m.m.r
  !> @param[in,out] oclo                 chlorine_dioxide m.m.r
  !> @param[in,out] br                   bromine m.m.r
  !> @param[in,out] lumped_br            lumped_br_as_bromine_monoxide m.m.r
  !> @param[in,out] brcl                 bromine_chloride m.m.r
  !> @param[in,out] brono2               bromine_nitrate m.m.r
  !> @param[in,out] n2o                  nitrous_oxide m.m.r
  !> @param[in,out] lumped_cl            lumped_cl_as_hydrogen_chloride m.m.r
  !> @param[in,out] hocl                 hydrochlorous_acid m.m.r
  !> @param[in,out] hbr                  hydrogen_bromide m.m.r
  !> @param[in,out] hobr                 hydrobromous_acid m.m.r
  !> @param[in,out] clono2               chlorine_nitrate m.m.r
  !> @param[in,out] cfcl3                cfc11 m.m.r
  !> @param[in,out] cf2cl2               cfc12 m.m.r
  !> @param[in,out] mebr                 methyl_bromide m.m.r
  !> @param[in,out] hono                 nitrous_acid m.m.r
  !> @param[in,out] c2h6                 ethane m.m.r
  !> @param[in,out] etoo                 ethyl_peroxy_radical m.m.r
  !> @param[in,out] etooh                ethyl_hydroperoxide m.m.r
  !> @param[in,out] mecho                acetaldehyde m.m.r
  !> @param[in,out] meco3                peroxy_acetic_acide m.m.r
  !> @param[in,out] pan                  peroxy_acetyl_nitrate m.m.r
  !> @param[in,out] c3h8                 propane m.m.r
  !> @param[in,out] n_proo               n-propyl_peroxy_radical m.m.r
  !> @param[in,out] i_proo               i-propyl_peroxy_radical m.m.r
  !> @param[in,out] n_prooh              n-propyl_alcohol m.m.r
  !> @param[in,out] i_prooh              i-propyl_alcohol m.m.r
  !> @param[in,out] etcho                propionaldehyde m.m.r
  !> @param[in,out] etco3                ethyl_carbonate m.m.r
  !> @param[in,out] me2co                acetone m.m.r
  !> @param[in,out] mecoch2oo            acetone_peroxy_radical m.m.r
  !> @param[in,out] mecoch2ooh           acetone_hydroperoxide m.m.r
  !> @param[in,out] ppan                 peroxypropionyl_nitrate m.m.r
  !> @param[in,out] meono2               methyl_nitrate m.m.r
  !> @param[in,out] c5h8                 isoprene m.m.r
  !> @param[in,out] iso2                 isoprene_peroxy_radical m.m.r
  !> @param[in,out] isooh                isoprene_hydroperoxide m.m.r
  !> @param[in,out] ison                 isoprene_peroxy_acetyl_nitrate m.m.r
  !> @param[in,out] macr                 methacrolein m.m.r
  !> @param[in,out] macro2               lumped_mvk_macr_hydroperoxy_radical m.m.r
  !> @param[in,out] macrooh              lumped_mvk_macr_hydroperoxide m.m.r
  !> @param[in,out] mpan                 methacryolyl_peroxynitrate m.m.r
  !> @param[in,out] hacet                hydroxy_acetone m.m.r
  !> @param[in,out] mgly                 methylglyoxal m.m.r
  !> @param[in,out] nald                 nitroxy_acetaldehyde m.m.r
  !> @param[in,out] hcooh                formic_acid m.m.r
  !> @param[in,out] meco3h               peroxy_acetic_acid m.m.r
  !> @param[in,out] meco2h               acetic_acid m.m.r
  !> @param[in,out] h2                   hydrogen m.m.r
  !> @param[in,out] meoh                 methanol m.m.r
  !> @param[in,out] msa                  methyl_sulphonic_acid m.m.r
  !> @param[in,out] nh3                  ammonia m.m.r
  !> @param[in,out] cs2                  carbon_disulphide m.m.r
  !> @param[in,out] csul                 carbonyl_sulphide m.m.r
  !> @param[in,out] h2s                  hydrogen_sulphide m.m.r
  !> @param[in,out] so3                  sulphur_trioxide m.m.r
  !> @param[in,out] passive_o3           passive_ozone m.m.r
  !> @param[in,out] age_of_air           age_of_air m.m.r
  !> @param[in,out] dms                  Dimethyl sulfide m.m.r.
  !> @param[in,out] so2                  Sulfur dioxide m.m.r.
  !> @param[in,out] h2so4                Sulfuric acid m.m.r.
  !> @param[in,out] dmso                 Dimethyl sulfoxide m.m.r.
  !> @param[in,out] monoterpene          Monoterpene m.m.r.
  !> @param[in,out] secondary organic    Secondary organic m.m.r.
  !> @param[in,out] n_nuc_sol            Aerosol field: n.m.r. of soluble nucleation mode
  !> @param[in,out] nuc_sol_su           Aerosol field: m.m.r. of H2SO4 in soluble nucleation mode
  !> @param[in,out] nuc_sol_om           Aerosol field: m.m.r. of organic matter in soluble nucleation mode
  !> @param[in,out] n_ait_sol            Aerosol field: n.m.r. of soluble Aitken mode
  !> @param[in,out] ait_sol_su           Aerosol field: m.m.r. of H2SO4 in soluble Aitken mode
  !> @param[in,out] ait_sol_bc           Aerosol field: m.m.r. of black carbon in soluble Aitken mode
  !> @param[in,out] ait_sol_om           Aerosol field: m.m.r. of organic matter in soluble Aitken mode
  !> @param[in,out] n_acc_sol            Aerosol field: n.m.r. of soluble accumulation mode
  !> @param[in,out] acc_sol_su           Aerosol field: m.m.r. of H2SO4 in soluble accumulation mode
  !> @param[in,out] acc_sol_bc           Aerosol field: m.m.r. of black carbon in soluble accumulation mode
  !> @param[in,out] acc_sol_om           Aerosol field: m.m.r. of organic matter in soluble accumulation mode
  !> @param[in,out] acc_sol_ss           Aerosol field: m.m.r. of sea salt in soluble accumulation mode
  !> @param[in,out] n_cor_sol            Aerosol field: n.m.r. of soluble coarse mode
  !> @param[in,out] cor_sol_su           Aerosol field: m.m.r. of H2SO4 in soluble coarse mode
  !> @param[in,out] cor_sol_bc           Aerosol field: m.m.r. of black carbon in soluble coarse mode
  !> @param[in,out] cor_sol_om           Aerosol field: m.m.r. of organic matter in soluble coarse mode
  !> @param[in,out] cor_sol_ss           Aerosol field: m.m.r. of sea salt in soluble coarse mode
  !> @param[in,out] n_ait_ins            Aerosol field: n.m.r. of insoluble Aitken mode
  !> @param[in,out] ait_ins_bc           Aerosol field: m.m.r. of black carbon in insoluble Aitken mode
  !> @param[in,out] ait_ins_om           Aerosol field: m.m.r. of organic matter in insoluble Aitken mode
  !> @param[in,out] n_acc_ins            Aerosol field: n.m.r. of insoluble accumulation mode
  !> @param[in,out] acc_ins_du           Aerosol field: m.m.r. of dust in insoluble accumulation mode
  !> @param[in,out] n_cor_ins            Aerosol field: n.m.r. of insoluble coarse mode
  !> @param[in,out] cor_ins_du           Aerosol field: m.m.r. of dust in insoluble coarse mode
  !> @param[in]     tnuc                 Temperature of nucleation (K)
  !> @param[in,out] lowest_cv_base       Level number for start of convection in column
  !> @param[in,out] lowest_cv_top        Level number for end of lowest convection in column
  !> @param[in,out] cv_base              Level number of base of highest convection in column
  !> @param[in,out] cv_top               Level number for end of highest convection in column
  !> @param[in,out] pres_cv_base         Pressure at base of highest convection in column
  !> @param[in,out] pres_cv_top          Pressure at end of highest convection in column
  !> @param[in,out] pres_lowest_cv_base  Pressure at base of lowest convection in column
  !> @param[in,out] pres_lowest_cv_top   Pressure at end of lowest convection in column
  !> @param[in,out] lowest_cca_2d        2D convective cloud amount of lowest convecting layer
  !> @param[in,out] entrain_up           Convective upwards entrainment
  !> @param[in,out] entrain_down         Convective downwards entrainment
  !> @param[in,out] detrain_up           Convective upwards detrainment
  !> @param[in,out] detrain_down         Convective downwards detrainment
  !> @param[in,out] massflux_up_half     Convective upwards mass flux on half-levels (Pa/s)
  !> @param[in]     ndf_w3               Number of DOFs per cell for density space
  !> @param[in]     undf_w3              Number of unique DOFs  for density space
  !> @param[in]     map_w3               Dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_wth              Number of DOFs per cell for potential temperature space
  !> @param[in]     undf_wth             Number of unique DOFs for potential temperature space
  !> @param[in]     map_wth              Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_2d               Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d              Number of unique DOFs  for 2D fields
  !> @param[in]     map_2d               Dofmap for the cell at the base of the column for 2D fields
  !> @param[in]     ndf_bl               Number of DOFs per cell for BL types
  !> @param[in]     undf_bl              Number of total DOFs for BL types
  !> @param[in]     map_bl               Dofmap for cell for BL types
  !> @param[in]     ndf_surf             Number of DOFs per cell for surface variables
  !> @param[in]     undf_surf            Number of unique DOFs for surface variables
  !> @param[in]     map_surf             Dofmap for the cell at the base of the column for surface variables
  subroutine conv_comorph_code(nlayers, seg_len,             &
                          outer,                             &
                          rho_in_w3,                         &
                          rho_in_wth,                        &
                          wetrho_in_w3,                      &
                          wetrho_in_wth,                     &
                          exner_in_w3,                       &
                          exner_in_wth,                      &
                          delta,                             &
                          theta_n,                           &
                          u_in_w3,                           &
                          v_in_w3,                           &
                          w_in_wth,                          &
                          theta_star,                        &
                          u_in_w3_star,                      &
                          v_in_w3_star,                      &
                          height_w3,                         &
                          height_wth,                        &
                          dt_conv,                           &
                          dmv_conv,                          &
                          dmcl_conv,                         &
                          dms_conv,                          &
                          du_conv,                           &
                          dv_conv,                           &
                          conv_prog_dtheta,                  &
                          conv_prog_dmv,                     &
                          mv_n,                              &
                          mcl_n,                             &
                          mci_n,                             &
                          mr_n,                              &
                          mg_n,                              &
                          ms_n,                              &
                          m_v,                               &
                          m_cl,                              &
                          m_ci,                              &
                          m_r,                               &
                          m_g,                               &
                          m_s,                               &
                          cf_ice,                            &
                          cf_liq,                            &
                          cf_bulk,                           &
                          cca,                               &
                          ccw,                               &
                          precfrac,                          &
                          cf_liq_n,                          &
                          cf_fro_n,                          &
                          cf_bulk_n,                         &
                          cape_diluted,                      &
                          cca_2d,                            &
                          dd_mf_cb,                          &
                          massflux_up,                       &
                          massflux_down,                     &
                          tke_bl,                            &
                          pressure_inc_env,                  &
                          zh_2d,                             &
                          zh_nonloc,                         &
                          inv_depth,                         &
                          zhsc_2d,                           &
                          bl_type_ind,                       &
                          wvar,                              &
                          rhokm_bl,                          &
                          moist_flux,                        &
                          heat_flux,                         &
                          taux,                              &
                          tauy,                              &
                          surf_interp,                       &
                          ustar,                             &
                          ls_rain_2d,                        &
                          ls_snow_2d,                        &
                          dcfl_conv,                         &
                          dcff_conv,                         &
                          dbcf_conv,                         &
                          o3p,                               &
                          o1d,                               &
                          o3,                                &
                          nit,                               &
                          no,                                &
                          no3,                               &
                          lumped_n,                          &
                          n2o5,                              &
                          ho2no2,                            &
                          hono2,                             &
                          h2o2,                              &
                          ch4,                               &
                          co,                                &
                          hcho,                              &
                          meoo,                              &
                          meooh,                             &
                          h,                                 &
                          oh,                                &
                          ho2,                               &
                          cl,                                &
                          cl2o2,                             &
                          clo,                               &
                          oclo,                              &
                          br,                                &
                          lumped_br,                         &
                          brcl,                              &
                          brono2,                            &
                          n2o,                               &
                          lumped_cl,                         &
                          hocl,                              &
                          hbr,                               &
                          hobr,                              &
                          clono2,                            &
                          cfcl3,                             &
                          cf2cl2,                            &
                          mebr,                              &
                          hono,                              &
                          c2h6,                              &
                          etoo,                              &
                          etooh,                             &
                          mecho,                             &
                          meco3,                             &
                          pan,                               &
                          c3h8,                              &
                          n_proo,                            &
                          i_proo,                            &
                          n_prooh,                           &
                          i_prooh,                           &
                          etcho,                             &
                          etco3,                             &
                          me2co,                             &
                          mecoch2oo,                         &
                          mecoch2ooh,                        &
                          ppan,                              &
                          meono2,                            &
                          c5h8,                              &
                          iso2,                              &
                          isooh,                             &
                          ison,                              &
                          macr,                              &
                          macro2,                            &
                          macrooh,                           &
                          mpan,                              &
                          hacet,                             &
                          mgly,                              &
                          nald,                              &
                          hcooh,                             &
                          meco3h,                            &
                          meco2h,                            &
                          h2,                                &
                          meoh,                              &
                          msa,                               &
                          nh3,                               &
                          cs2,                               &
                          csul,                              &
                          h2s,                               &
                          so3,                               &
                          passive_o3,                        &
                          age_of_air,                        &
                          dms,                               &
                          so2,                               &
                          h2so4,                             &
                          dmso,                              &
                          monoterpene,                       &
                          secondary_organic,                 &
                          n_nuc_sol,                         &
                          nuc_sol_su,                        &
                          nuc_sol_om,                        &
                          n_ait_sol,                         &
                          ait_sol_su,                        &
                          ait_sol_bc,                        &
                          ait_sol_om,                        &
                          n_acc_sol,                         &
                          acc_sol_su,                        &
                          acc_sol_bc,                        &
                          acc_sol_om,                        &
                          acc_sol_ss,                        &
                          n_cor_sol,                         &
                          cor_sol_su,                        &
                          cor_sol_bc,                        &
                          cor_sol_om,                        &
                          cor_sol_ss,                        &
                          n_ait_ins,                         &
                          ait_ins_bc,                        &
                          ait_ins_om,                        &
                          n_acc_ins,                         &
                          acc_ins_du,                        &
                          n_cor_ins,                         &
                          cor_ins_du,                        &
                          tnuc,                              &
                          lowest_cv_base,                    &
                          lowest_cv_top,                     &
                          cv_base,                           &
                          cv_top,                            &
                          pres_cv_base,                      &
                          pres_cv_top,                       &
                          pres_lowest_cv_base,               &
                          pres_lowest_cv_top,                &
                          lowest_cca_2d,                     &
                          entrain_up,                        &
                          entrain_down,                      &
                          detrain_up,                        &
                          detrain_down,                      &
                          massflux_up_half,                  &
                          ndf_w3,                            &
                          undf_w3,                           &
                          map_w3,                            &
                          ndf_wth,                           &
                          undf_wth,                          &
                          map_wth,                           &
                          ndf_2d,                            &
                          undf_2d,                           &
                          map_2d,                            &
                          ndf_bl, undf_bl, map_bl,           &
                          ndf_surf, undf_surf, map_surf)

    !---------------------------------------
    ! LFRic modules
    !---------------------------------------
    use um_ukca_init_mod, only: fldname_o3p,                                   &
                                fldname_o1d,                                   &
                                fldname_o3,                                    &
                                fldname_n,                                     &
                                fldname_no,                                    &
                                fldname_no3,                                   &
                                fldname_lumped_n,                              &
                                fldname_n2o5,                                  &
                                fldname_ho2no2,                                &
                                fldname_hono2,                                 &
                                fldname_h2o2,                                  &
                                fldname_ch4,                                   &
                                fldname_co,                                    &
                                fldname_hcho,                                  &
                                fldname_meoo,                                  &
                                fldname_meooh,                                 &
                                fldname_h,                                     &
                                fldname_ch2o,                                  &
                                fldname_oh,                                    &
                                fldname_ho2,                                   &
                                fldname_cl,                                    &
                                fldname_cl2o2,                                 &
                                fldname_clo,                                   &
                                fldname_oclo,                                  &
                                fldname_br,                                    &
                                fldname_lumped_br,                             &
                                fldname_brcl,                                  &
                                fldname_brono2,                                &
                                fldname_n2o,                                   &
                                fldname_lumped_cl,                             &
                                fldname_hocl,                                  &
                                fldname_hbr,                                   &
                                fldname_hobr,                                  &
                                fldname_clono2,                                &
                                fldname_cfcl3,                                 &
                                fldname_cf2cl2,                                &
                                fldname_mebr,                                  &
                                fldname_hono,                                  &
                                fldname_c2h6,                                  &
                                fldname_etoo,                                  &
                                fldname_etooh,                                 &
                                fldname_mecho,                                 &
                                fldname_meco3,                                 &
                                fldname_pan,                                   &
                                fldname_c3h8,                                  &
                                fldname_n_proo,                                &
                                fldname_i_proo,                                &
                                fldname_n_prooh,                               &
                                fldname_i_prooh,                               &
                                fldname_etcho,                                 &
                                fldname_etco3,                                 &
                                fldname_me2co,                                 &
                                fldname_mecoch2oo,                             &
                                fldname_mecoch2ooh,                            &
                                fldname_ppan,                                  &
                                fldname_meono2,                                &
                                fldname_c5h8,                                  &
                                fldname_iso2,                                  &
                                fldname_isooh,                                 &
                                fldname_ison,                                  &
                                fldname_macr,                                  &
                                fldname_macro2,                                &
                                fldname_macrooh,                               &
                                fldname_mpan,                                  &
                                fldname_hacet,                                 &
                                fldname_mgly,                                  &
                                fldname_nald,                                  &
                                fldname_hcooh,                                 &
                                fldname_meco3h,                                &
                                fldname_meco2h,                                &
                                fldname_h2,                                    &
                                fldname_meoh,                                  &
                                fldname_msa,                                   &
                                fldname_nh3,                                   &
                                fldname_cs2,                                   &
                                fldname_csul,                                  &
                                fldname_h2s,                                   &
                                fldname_so3,                                   &
                                fldname_passive_o3,                            &
                                fldname_age_of_air,                            &
                                fldname_dms,                                   &
                                fldname_so2,                                   &
                                fldname_h2so4,                                 &
                                fldname_dmso,                                  &
                                fldname_monoterpene,                           &
                                fldname_secondary_organic,                     &
                                fldname_n_nuc_sol,                             &
                                fldname_nuc_sol_su,                            &
                                fldname_nuc_sol_om,                            &
                                fldname_n_ait_sol,                             &
                                fldname_ait_sol_su,                            &
                                fldname_ait_sol_bc,                            &
                                fldname_ait_sol_om,                            &
                                fldname_n_acc_sol,                             &
                                fldname_acc_sol_su,                            &
                                fldname_acc_sol_bc,                            &
                                fldname_acc_sol_om,                            &
                                fldname_acc_sol_ss,                            &
                                fldname_acc_sol_du,                            &
                                fldname_n_cor_sol,                             &
                                fldname_cor_sol_su,                            &
                                fldname_cor_sol_bc,                            &
                                fldname_cor_sol_om,                            &
                                fldname_cor_sol_ss,                            &
                                fldname_cor_sol_du,                            &
                                fldname_n_ait_ins,                             &
                                fldname_ait_ins_bc,                            &
                                fldname_ait_ins_om,                            &
                                fldname_n_acc_ins,                             &
                                fldname_acc_ins_du,                            &
                                fldname_n_cor_ins,                             &
                                fldname_cor_ins_du

    use aerosol_config_mod,        only: glomap_mode,               &
                                         glomap_mode_dust_and_clim, &
                                         glomap_mode_ukca

    use log_mod, only : log_event, log_scratch_space, LOG_LEVEL_ERROR
!$  use omp_lib, only : omp_get_max_threads

    !---------------------------------------
    ! UM modules containing switches or global constants
    !---------------------------------------
    use bl_option_mod, only: max_tke
    use cloud_inputs_mod, only: l_pc2_homog_conv_pressure
    use cv_run_mod, only: l_mom,                                          &
                          l_conv_prog_dtheta, l_conv_prog_dq,             &
                          tau_conv_prog_dtheta, tau_conv_prog_dq
    use jules_surface_mod, only: srf_ex_cnv_gust, IP_SrfExWithCnv
    use mphys_inputs_mod, only: l_mcr_qgraup, l_mcr_qrain, l_mcr_precfrac, &
         l_mcr_qcf2
    use nlsizes_namelist_mod, only: row_length, rows, bl_levels
    use planet_constants_mod, only: p_zero, kappa, planet_radius, g
    use timestep_mod, only: timestep
    use conversions_mod, only: zerodegc

    ! subroutines used
    use ukca_api_mod, only: ukca_get_tracer_varlist, ukca_maxlen_fieldname
    use comorph_ctl_mod, only: comorph_ctl
    use fields_type_mod, only: fields_type, fields_nullify
    use grid_type_mod, only: grid_type, grid_nullify
    use turb_type_mod, only: turb_type, turb_nullify
    use cloudfracs_type_mod, only: cloudfracs_type, cloudfracs_nullify
    use comorph_diags_type_mod, only: comorph_diags_type
    use set_constants_from_um_mod, only: set_constants_from_um
    use comorph_constants_mod, only: l_init_constants, l_turb_par_gen, &
         l_cv_rain, l_cv_cf, l_cv_snow, l_cv_graup
    use calc_conv_incs_mod, only: calc_conv_incs, i_call_save_before_conv, &
         i_call_diff_to_get_incs
    use fracs_consistency_mod, only: fracs_consistency
    use conv_update_precfrac_mod, only: conv_update_precfrac
    use interp_turb_mod, only: interp_turb
    use calc_turb_len_mod, only: calc_turb_len
    use limit_turb_perts_mod, only: limit_turb_perts
    use raise_error_mod, only: raise_fatal
    use umprintmgr, only: newline
    use assign_fields_mod, only: assign_fields
    use comorph_conv_cloud_extras_mod, only: comorph_conv_cloud_extras

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, seg_len
    integer(kind=i_def), intent(in) :: outer

    integer(kind=i_def), intent(in) :: ndf_wth, ndf_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: ndf_bl, undf_bl
    integer(kind=i_def), intent(in) :: ndf_surf, undf_surf
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3
    integer(kind=i_def), intent(in) :: map_wth(ndf_wth,seg_len)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3,seg_len)
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d,seg_len)
    integer(kind=i_def), intent(in) :: map_bl(ndf_bl,seg_len)
    integer(kind=i_def), intent(in) :: map_surf(ndf_surf,seg_len)

    real(kind=r_def), dimension(undf_w3), intent(in) :: rho_in_w3,          &
                                                        wetrho_in_w3,       &
                                                        exner_in_w3,        &
                                                        u_in_w3,            &
                                                        v_in_w3,            &
                                                        u_in_w3_star,       &
                                                        v_in_w3_star,       &
                                                        height_w3,          &
                                                        heat_flux,          &
                                                        moist_flux,         &
                                                        taux, tauy
    real(kind=r_def), dimension(undf_wth), intent(in) :: cf_ice,            &
                                                         cf_liq, cf_bulk,   &
                                                         mv_n, mcl_n, mci_n,&
                                                         mr_n, mg_n, ms_n,  &
                                                         m_v, m_cl, m_s,    &
                                                         rho_in_wth,        &
                                                         wetrho_in_wth,     &
                                                         exner_in_wth,      &
                                                         w_in_wth, delta,   &
                                                         theta_n,           &
                                                         theta_star,        &
                                                         height_wth,        &
                                                         rhokm_bl, wvar,    &
                                                         cf_liq_n, cf_fro_n,&
                                                         cf_bulk_n

    real(kind=r_def), dimension(undf_wth), intent(inout) :: dt_conv, dmv_conv, &
                                          dmcl_conv, dms_conv, cca, ccw,       &
                                          massflux_up, massflux_down,          &
                                          tke_bl, pressure_inc_env,            &
                                          conv_prog_dtheta, conv_prog_dmv,     &
                                          precfrac, m_r, m_g, m_ci

    real(kind=r_def), dimension(undf_w3), intent(inout) :: du_conv, dv_conv

    real(kind=r_def), dimension(undf_2d), intent(in) :: zh_2d,                &
                                                        zh_nonloc, inv_depth, &
                                                        zhsc_2d, ustar,       &
                                                        ls_rain_2d, ls_snow_2d
    integer(kind=i_def), dimension(undf_bl), intent(inout) :: bl_type_ind
    real(kind=r_def), dimension(undf_surf), intent(inout) :: surf_interp

    real(kind=r_def), dimension(undf_2d), intent(inout) :: cape_diluted,  &
                                                           cca_2d, dd_mf_cb
    real(kind=r_def), intent(in out), dimension(undf_wth) :: o3p
    real(kind=r_def), intent(in out), dimension(undf_wth) :: o1d
    real(kind=r_def), intent(in out), dimension(undf_wth) :: o3
    real(kind=r_def), intent(in out), dimension(undf_wth) :: nit
    real(kind=r_def), intent(in out), dimension(undf_wth) :: no
    real(kind=r_def), intent(in out), dimension(undf_wth) :: no3
    real(kind=r_def), intent(in out), dimension(undf_wth) :: lumped_n
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n2o5
    real(kind=r_def), intent(in out), dimension(undf_wth) :: ho2no2
    real(kind=r_def), intent(in out), dimension(undf_wth) :: hono2
    real(kind=r_def), intent(in out), dimension(undf_wth) :: h2o2
    real(kind=r_def), intent(in out), dimension(undf_wth) :: ch4
    real(kind=r_def), intent(in out), dimension(undf_wth) :: co
    real(kind=r_def), intent(in out), dimension(undf_wth) :: hcho
    real(kind=r_def), intent(in out), dimension(undf_wth) :: meoo
    real(kind=r_def), intent(in out), dimension(undf_wth) :: meooh
    real(kind=r_def), intent(in out), dimension(undf_wth) :: h
    real(kind=r_def), intent(in out), dimension(undf_wth) :: oh
    real(kind=r_def), intent(in out), dimension(undf_wth) :: ho2
    real(kind=r_def), intent(in out), dimension(undf_wth) :: cl
    real(kind=r_def), intent(in out), dimension(undf_wth) :: cl2o2
    real(kind=r_def), intent(in out), dimension(undf_wth) :: clo
    real(kind=r_def), intent(in out), dimension(undf_wth) :: oclo
    real(kind=r_def), intent(in out), dimension(undf_wth) :: br
    real(kind=r_def), intent(in out), dimension(undf_wth) :: lumped_br
    real(kind=r_def), intent(in out), dimension(undf_wth) :: brcl
    real(kind=r_def), intent(in out), dimension(undf_wth) :: brono2
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n2o
    real(kind=r_def), intent(in out), dimension(undf_wth) :: lumped_cl
    real(kind=r_def), intent(in out), dimension(undf_wth) :: hocl
    real(kind=r_def), intent(in out), dimension(undf_wth) :: hbr
    real(kind=r_def), intent(in out), dimension(undf_wth) :: hobr
    real(kind=r_def), intent(in out), dimension(undf_wth) :: clono2
    real(kind=r_def), intent(in out), dimension(undf_wth) :: cfcl3
    real(kind=r_def), intent(in out), dimension(undf_wth) :: cf2cl2
    real(kind=r_def), intent(in out), dimension(undf_wth) :: mebr
    real(kind=r_def), intent(in out), dimension(undf_wth) :: hono
    real(kind=r_def), intent(in out), dimension(undf_wth) :: c2h6
    real(kind=r_def), intent(in out), dimension(undf_wth) :: etoo
    real(kind=r_def), intent(in out), dimension(undf_wth) :: etooh
    real(kind=r_def), intent(in out), dimension(undf_wth) :: mecho
    real(kind=r_def), intent(in out), dimension(undf_wth) :: meco3
    real(kind=r_def), intent(in out), dimension(undf_wth) :: pan
    real(kind=r_def), intent(in out), dimension(undf_wth) :: c3h8
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n_proo
    real(kind=r_def), intent(in out), dimension(undf_wth) :: i_proo
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n_prooh
    real(kind=r_def), intent(in out), dimension(undf_wth) :: i_prooh
    real(kind=r_def), intent(in out), dimension(undf_wth) :: etcho
    real(kind=r_def), intent(in out), dimension(undf_wth) :: etco3
    real(kind=r_def), intent(in out), dimension(undf_wth) :: me2co
    real(kind=r_def), intent(in out), dimension(undf_wth) :: mecoch2oo
    real(kind=r_def), intent(in out), dimension(undf_wth) :: mecoch2ooh
    real(kind=r_def), intent(in out), dimension(undf_wth) :: ppan
    real(kind=r_def), intent(in out), dimension(undf_wth) :: meono2
    real(kind=r_def), intent(in out), dimension(undf_wth) :: c5h8
    real(kind=r_def), intent(in out), dimension(undf_wth) :: iso2
    real(kind=r_def), intent(in out), dimension(undf_wth) :: isooh
    real(kind=r_def), intent(in out), dimension(undf_wth) :: ison
    real(kind=r_def), intent(in out), dimension(undf_wth) :: macr
    real(kind=r_def), intent(in out), dimension(undf_wth) :: macro2
    real(kind=r_def), intent(in out), dimension(undf_wth) :: macrooh
    real(kind=r_def), intent(in out), dimension(undf_wth) :: mpan
    real(kind=r_def), intent(in out), dimension(undf_wth) :: hacet
    real(kind=r_def), intent(in out), dimension(undf_wth) :: mgly
    real(kind=r_def), intent(in out), dimension(undf_wth) :: nald
    real(kind=r_def), intent(in out), dimension(undf_wth) :: hcooh
    real(kind=r_def), intent(in out), dimension(undf_wth) :: meco3h
    real(kind=r_def), intent(in out), dimension(undf_wth) :: meco2h
    real(kind=r_def), intent(in out), dimension(undf_wth) :: h2
    real(kind=r_def), intent(in out), dimension(undf_wth) :: meoh
    real(kind=r_def), intent(in out), dimension(undf_wth) :: msa
    real(kind=r_def), intent(in out), dimension(undf_wth) :: nh3
    real(kind=r_def), intent(in out), dimension(undf_wth) :: cs2
    real(kind=r_def), intent(in out), dimension(undf_wth) :: csul
    real(kind=r_def), intent(in out), dimension(undf_wth) :: h2s
    real(kind=r_def), intent(in out), dimension(undf_wth) :: so3
    real(kind=r_def), intent(in out), dimension(undf_wth) :: passive_o3
    real(kind=r_def), intent(in out), dimension(undf_wth) :: age_of_air
    real(kind=r_def), intent(in out), dimension(undf_wth) :: dms
    real(kind=r_def), intent(in out), dimension(undf_wth) :: so2
    real(kind=r_def), intent(in out), dimension(undf_wth) :: h2so4
    real(kind=r_def), intent(in out), dimension(undf_wth) :: dmso
    real(kind=r_def), intent(in out), dimension(undf_wth) :: monoterpene
    real(kind=r_def), intent(in out), dimension(undf_wth) :: secondary_organic
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n_nuc_sol
    real(kind=r_def), intent(in out), dimension(undf_wth) :: nuc_sol_su
    real(kind=r_def), intent(in out), dimension(undf_wth) :: nuc_sol_om
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n_ait_sol
    real(kind=r_def), intent(in out), dimension(undf_wth) :: ait_sol_su
    real(kind=r_def), intent(in out), dimension(undf_wth) :: ait_sol_bc
    real(kind=r_def), intent(in out), dimension(undf_wth) :: ait_sol_om
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n_acc_sol
    real(kind=r_def), intent(in out), dimension(undf_wth) :: acc_sol_su
    real(kind=r_def), intent(in out), dimension(undf_wth) :: acc_sol_bc
    real(kind=r_def), intent(in out), dimension(undf_wth) :: acc_sol_om
    real(kind=r_def), intent(in out), dimension(undf_wth) :: acc_sol_ss
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n_cor_sol
    real(kind=r_def), intent(in out), dimension(undf_wth) :: cor_sol_su
    real(kind=r_def), intent(in out), dimension(undf_wth) :: cor_sol_bc
    real(kind=r_def), intent(in out), dimension(undf_wth) :: cor_sol_om
    real(kind=r_def), intent(in out), dimension(undf_wth) :: cor_sol_ss
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n_ait_ins
    real(kind=r_def), intent(in out), dimension(undf_wth) :: ait_ins_bc
    real(kind=r_def), intent(in out), dimension(undf_wth) :: ait_ins_om
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n_acc_ins
    real(kind=r_def), intent(in out), dimension(undf_wth) :: acc_ins_du
    real(kind=r_def), intent(in out), dimension(undf_wth) :: n_cor_ins
    real(kind=r_def), intent(in out), dimension(undf_wth) :: cor_ins_du
    real(kind=r_def), intent(in),     dimension(undf_wth) :: tnuc

    real(kind=r_def), pointer, intent(inout) :: lowest_cv_base(:),         &
                                                lowest_cv_top(:),          &
                                                cv_base(:),                &
                                                cv_top(:),                 &
                                                pres_cv_base(:),           &
                                                pres_cv_top(:),            &
                                                pres_lowest_cv_base(:),    &
                                                pres_lowest_cv_top(:),     &
                                                lowest_cca_2d(:)

    real(kind=r_def), pointer, intent(inout) :: entrain_up(:),       &
                                                entrain_down(:),     &
                                                detrain_up(:),       &
                                                detrain_down(:),     &
                                                massflux_up_half(:)

    real(kind=r_def), dimension(undf_wth), intent(inout) :: dcfl_conv
    real(kind=r_def), dimension(undf_wth), intent(inout) :: dcff_conv
    real(kind=r_def), dimension(undf_wth), intent(inout) :: dbcf_conv

    !-----------------------------------------------------------------------
    ! Local variables for the kernel
    !-----------------------------------------------------------------------
    ! loop counters etc
    integer(i_def) :: k, i, n

    ! local switches and scalars
    integer(i_um) :: segments, n_conv_levels, ntra_fld

    logical :: l_tracer

    real(r_um) :: decay_amount, interp

    ! profile fields from level 1 upwards
    real(r_um), dimension(row_length,rows,nlayers) ::                        &
         rho_wet, rho_dry, z_rho, z_theta, cca_3d, rho_wet_tq,               &
         rho_dry_tq, r_rho_levels,                                           &
         theta_conv, q_conv, qcl_conv, qcf_conv, dtheta_conv,                &
         qrain_conv, qcf2_conv, qgraup_conv, cf_liquid_conv, cf_frozen_conv, &
         bulk_cf_conv, u_conv, v_conv, ccw_3d, dubydt_p,                     &
         dvbydt_p, tnuc_new, cca_3d0, ccw_3d0

    ! profile fields from level 0 upwards
    real(r_um), dimension(row_length,rows,0:nlayers) ::                      &
         p_theta_levels, w, r_theta_levels, exner_theta_levels
    real(r_um), dimension(row_length,rows,nlayers+1) :: p_rho_levels

    ! single level real fields
    real(r_um), dimension(row_length,rows) :: zh

    ! single level integer fields
    integer(i_um), dimension(row_length,rows) :: lcbase, ccb, cct, lctop, &
         ccb0, cct0, lcbase0

    ! total tracer - has to be allocatable
    real(r_um), dimension(:,:,:,:), allocatable, target :: tot_tracer

    ! Variables for retrieving tracer names for a UKCA configuration
    integer :: ukca_errcode
    character(len=ukca_maxlen_fieldname), pointer :: ukca_tracer_names(:)
    character(len=ukca_maxlen_fieldname), target  ::                         &
                 local_dust_tracer_list(4) = [ 'Acc_INS_N ' , 'Acc_INS_DU' , &
                                               'Cor_INS_N ' , 'Cor_INS_DU' ]

    ! Heat and moisture fluxs from BL scheme
    real(r_um), dimension(row_length,rows,bl_levels) :: fqw, ftl

    ! Comorph variables
    ! Derived type structures storing pointers to the primary fields,
    ! at start-of-timestep and latest fields
    type(fields_type) :: fields_n
    type(fields_type) :: fields_np1

    ! Structure storing pointers to the model-grid fields and dry-rho
    type(grid_type) :: grid

    ! Structure containing pointers to turbulence fields passed
    ! in from the boundary-layer scheme
    type(turb_type) :: turb

    ! Structure containing diagnostic area fractions for
    ! cloud and rain
    type(cloudfracs_type) :: cloudfracs

    ! Structure storing meta-data and pointers to fields for all
    ! diagnostics calculated by CoMorph
    type(comorph_diags_type) :: comorph_diags

    ! Flag for whether to update w
    logical, parameter :: l_conv_inc_w = .false.

    ! Winds interpolated onto theta-levels
    real(kind=r_um), target :: u_th_n( row_length, rows, 1:nlayers-1)
    real(kind=r_um), target :: v_th_n( row_length, rows, 1:nlayers-1)
    real(kind=r_um), target :: u_th_np1( row_length, rows, 1:nlayers-1)
    real(kind=r_um), target :: v_th_np1( row_length, rows, 1:nlayers-1)
    ! Note: these do not exist at the top theta-level.

    real(kind=r_um), target, allocatable :: w_work(:,:,:)

    ! Wind velocity components / ms-1
    real(kind=r_um) :: u_p(row_length,rows,nlayers)
    real(kind=r_um) :: v_p(row_length,rows,nlayers)

    ! Vertical velocity increment:
    ! IN:  Increment so far this timestep
    ! OUT: Increment with contribution from convection added as well
    real(kind=r_um), target :: r_w(row_length, rows, 1:nlayers)

    real(kind=r_um) :: q_inc( row_length, rows, 1:nlayers )
    real(kind=r_um) :: qcl_inc( row_length, rows, 1:nlayers )
    real(kind=r_um) :: qcf_inc( row_length, rows, 1:nlayers )
    real(kind=r_um) :: qcf2_inc( row_length, rows, 1:nlayers )
    real(kind=r_um) :: qrain_inc( row_length, rows, 1:nlayers )
    real(kind=r_um) :: qgraup_inc( row_length, rows, 1:nlayers )
    real(KIND=r_um) :: cf_liquid_inc( row_length, rows, 1:nlayers )
    real(KIND=r_um) :: cf_frozen_inc( row_length, rows, 1:nlayers )
    real(KIND=r_um) :: bulk_cf_inc( row_length, rows, 1:nlayers )

    ! Array for start-of-timestep temperature
    ! (needs to be converted from theta to T)
    real(kind=r_um), target :: temperature_n( row_length, rows, nlayers )

    ! Rain fraction
    real(kind=r_um), target :: precfrac_star( row_length, rows, nlayers )

    ! Store to save qrain+qgraup before convection, for calculating precip
    ! fraction increment from convection
    real(kind=r_um), allocatable :: q_prec_b4(:,:,:)

    ! Separate convective bulk cloud fraction
    ! (doesn't yet have a separate field in the UM)
    real(kind=r_um), target :: frac_bulk_conv(row_length,rows,nlayers)

    ! "Effective" boundary-layer-top height
    real(kind=r_um), target :: zh_eff       ( row_length, rows )
    ! Resolved inversion thickness / m
    real(kind=r_um) :: dzh      ( row_length, rows )
    ! Surface-driven non-local BL-top height / m
    real(kind=r_um) :: zhnl     ( row_length, rows )
    ! Height of decoupled stratocumulus top / m
    real(kind=r_um) :: zhsc     ( row_length, rows )
    ! Indicator for shear-dominated boundary-layers
    real(kind=r_um) :: bl_type_7( row_length, rows )

    ! Arrays for turbulence fields interpolated onto rho-levels
    real(kind=r_um), target, allocatable :: w_var_rh(:,:,:)
    real(kind=r_um), target, allocatable :: fu_rh(:,:,:)
    real(kind=r_um), target, allocatable :: fv_rh(:,:,:)
    ! Turbulence lengthscale
    real(kind=r_um), target, allocatable :: turb_len(:,:,:)

    ! Sub-grid turbulent variance in vertical velocity
    real(kind=r_um) :: bl_w_var ( row_length, rows, 1:nlayers)
    ! Surface buoyancy flux / m2 s-3
    real(kind=r_um) :: fb_surf  ( row_length, rows)
    ! Surface friction velocity / m s-1
    real(kind=r_um) :: u_s      ( row_length, rows)
    ! Turbulent wind stress on p grid (on theta-levels) N m-2
    real(kind=r_um) :: taux_p ( row_length, rows, 0:bl_levels-1 )
    real(kind=r_um) :: tauy_p ( row_length, rows, 0:bl_levels-1 )
    ! rho * turbulent diffusivity for momentum, on theta-grid
    real(kind=r_um) :: rhokm ( row_length, rows, 0:bl_levels-1 )

    ! Surface precipitation rates from the "large-scale" microphysics scheme
    real(kind=r_um) :: ls_rain( row_length, rows)
    real(kind=r_um) :: ls_snow( row_length, rows)
    ! Grid size used in radius dependence on precip rate
    real(kind=r_um) :: delta_x( row_length, rows)

    ! Parcel radius amplification factor passed into comorph
    real(kind=r_um), target :: par_radius_amp_um(row_length,rows)

    character(len=*), parameter :: routinename = "conv_comorph_kernel"

    ! Flags for treating condensed water species diagnostically when
    ! they are off in the UM but on in CoMorph
    logical :: l_temporary_rain
    logical :: l_temporary_snow
    logical :: l_temporary_graup

    ! Work arrays for water species if they are on in CoMorph
    ! but off in the UM
    real(kind=r_um), target, allocatable :: q_rain_work(:,:,:)
    real(kind=r_um), target, allocatable :: q_snow_work(:,:,:)
    real(kind=r_um), target, allocatable :: q_graup_work(:,:,:)

    ! Cloud fraction fields
    real(kind=r_um), target :: cf_liquid_n(row_length,rows,0:nlayers)
    real(kind=r_um), target :: cf_frozen_n(row_length,rows,0:nlayers)
    real(kind=r_um), target :: bulk_cf_n(row_length,rows,0:nlayers)

    ! Start-of-timestep mixing ratios of water species
    ! Water vapour
    real(kind=r_um), target:: qv_n(row_length,rows,0:nlayers)
    ! Cloud liquid water
    real(kind=r_um), target :: qcl_n(row_length,rows,0:nlayers)
    ! Cloud ice and snow
    real(kind=r_um), target :: qcf_n(row_length,rows,0:nlayers)
    ! 2nd cloud ice category (optional)
    real(kind=r_um), target :: qcf2_n(row_length,rows,0:nlayers)
    ! Prognostic rain water (optional)
    real(kind=r_um), target :: qr_n(row_length,rows,0:nlayers)
    ! Prognostic graupel (optional)
    real(kind=r_um), target :: qgr_n(row_length,rows,0:nlayers)

    ! Counters for tracer fields
    integer :: i_field

    real(kind=r_um) :: cclwp  (row_length,rows)
    real(kind=r_um) :: cclwp0  (row_length,rows)
    real(kind=r_um) :: cca_2d_loc (row_length,rows)
    real(kind=r_um) :: lcca   (row_length,rows)

    ! Diagnostic fields
    real(kind=r_um), target :: cape_dil(row_length, rows)
    real(kind=r_um), target :: up_flux_half(row_length, rows, nlayers)
    real(kind=r_um), target :: down_flux_half(row_length, rows, nlayers)
    real(kind=r_um), target, allocatable, dimension(:,:,:) :: ent_up, ent_down,&
         det_up, det_down, pres_inc_env

    !-----------------------------------------------------------------------
    ! Mapping of LFRic fields into UM variables
    !-----------------------------------------------------------------------
    ! For the initial implementation we pass each individual column
    ! of data to an array sized (1,1,k) to match the UMs (i,j,k) data
    ! layout.
    ! assuming map_wth(1) points to level 0
    ! and map_w3(1) points to level 1
    !-----------------------------------------------------------------------
    do i = 1, row_length
      do k = 0, nlayers
        ! w wind on theta levels
        w(i,1,k) = w_in_wth(map_wth(1,i) + k)
        ! pressure on theta levels
        p_theta_levels(i,1,k) = p_zero*(exner_in_wth(map_wth(1,i) + k))**(1.0_r_def/kappa)
        exner_theta_levels(i,1,k) = exner_in_wth(map_wth(1,i) + k)
        ! height of theta levels from centre of planet
        r_theta_levels(i,1,k) = height_wth(map_wth(1,i) + k) + planet_radius
      end do
      do k = 1, nlayers
        ! wet density on theta and rho levels
        rho_wet_tq(i,1,k) = wetrho_in_wth(map_wth(1,i) + k)
        rho_wet(i,1,k) = wetrho_in_w3(map_w3(1,i) + k-1)
        ! dry density on theta and rho levels
        rho_dry_tq(i,1,k) = rho_in_wth(map_wth(1,i) + k)
        rho_dry(i,1,k) = rho_in_w3(map_w3(1,i) + k-1)
        ! pressure on rho levels
        p_rho_levels(i,1,k) = p_zero*(exner_in_w3(map_w3(1,i) + k-1))**(1.0_r_def/kappa)
        ! height of rho levels from centre of planet
        r_rho_levels(i,1,k) = height_w3(map_w3(1,i) + k-1) + planet_radius
        ! height of levels above surface
        z_rho(i,1,k) = r_rho_levels(i,1,k)-r_theta_levels(i,1,0)
        z_theta(i,1,k) = r_theta_levels(i,1,k)-r_theta_levels(i,1,0)
      end do
      ! adjusted from what is initialised above to match what Comorph expects
      p_rho_levels(i,1,1) = p_theta_levels(i,1,0)
      p_rho_levels(i,1,nlayers+1) = 0.0_r_um
    end do

    !-----------------------------------------------------------------------
    ! Things passed from other parametrization schemes on this timestep
    !-----------------------------------------------------------------------
    do i = 1, row_length
      zh(i,1) = zh_2d(map_2d(1,i))
      zhnl(i,1) = zh_nonloc(map_2d(1,i))
      dzh(i,1) = inv_depth(map_2d(1,i))
      zhsc(i,1) = zhsc_2d(map_2d(1,i))
      bl_type_7(i,1) = bl_type_ind(map_bl(1,i)+6)
    end do

    ! If this is the last solver outer loop then tracers may need convecting.
    ! Enable tracers for UKCA if a UKCA tracer list is available
    ! (This indicates that a UKCA configuration has been set up)
    ! Comorph needs tracer sizes setup on first call (which isn't an outer
    ! loop), hence establish these for all loops then set l_tracer
    ! to false unless outer loop
    if ( glomap_mode == glomap_mode_dust_and_clim ) then
      ukca_tracer_names => local_dust_tracer_list
      l_tracer = .true.
    else
      call ukca_get_tracer_varlist( ukca_tracer_names, ukca_errcode )
      l_tracer = ( ukca_errcode == 0 )
    end if
    if (l_tracer) then
      ntra_fld = size(ukca_tracer_names)
    else
      ntra_fld = 1
    end if
    if (outer /= outer_iterations) then
      l_tracer = .false.
    end if

    ! If not using OMP, just run with a single segment for all columns.
    segments = 1          ! i.e. one column or whole mpi rank
    ! Under OMP sentinel, set number of segments equal to number of threads
!$  segments = omp_get_max_threads()

    ! Set number of layers used by convection scheme
    n_conv_levels = nlayers - 1

    do i = 1, row_length
      do k = 1, nlayers
        theta_conv(i,1,k) = theta_star(map_wth(1,i) + k)

        q_conv(i,1,k)   = m_v(map_wth(1,i) + k)
        qcl_conv(i,1,k) = m_cl(map_wth(1,i) + k)

        qv_n(i,1,k)   = mv_n(map_wth(1,i) + k)
        qcl_n(i,1,k) = mcl_n(map_wth(1,i) + k)

        cf_liquid_conv(i,1,k) = cf_liq(map_wth(1,i) + k)
        cf_frozen_conv(i,1,k) = cf_ice(map_wth(1,i) + k)
        bulk_cf_conv(i,1,k)   = cf_bulk(map_wth(1,i) + k)

        cf_liquid_n(i,1,k) = cf_liq_n(map_wth(1,i) + k)
        cf_frozen_n(i,1,k) = cf_fro_n(map_wth(1,i) + k)
        bulk_cf_n(i,1,k)   = cf_bulk_n(map_wth(1,i) + k)
      end do
    end do

    do i = 1, row_length
      do k = 1, nlayers
        qcf_conv(i,1,k) = m_s(map_wth(1,i) + k)
        qcf_n(i,1,k) = ms_n(map_wth(1,i) + k)
        ! no 2nd ice prognostic without casim
      end do
    end do
    if (microphysics_casim) then
      do i = 1, row_length
        do k = 1, nlayers
          qcf2_conv(i,1,k) = m_ci(map_wth(1,i) + k)
          qcf2_n(i,1,k) = mci_n(map_wth(1,i) + k)
        end do
      end do
    end if
    if (l_mcr_qrain) then
      do i = 1, row_length
        do k = 1, nlayers
          qrain_conv(i,1,k) = m_r(map_wth(1,i) + k)
          qr_n(i,1,k) = mr_n(map_wth(1,i) + k)
        end do
      end do
    end if
    if (l_mcr_qgraup) then
      do i = 1, row_length
        do k = 1, nlayers
          qgraup_conv(i,1,k) = m_g(map_wth(1,i) + k)
          qgr_n(i,1,k) = mg_n(map_wth(1,i) + k)
        end do
      end do
    end if

    if (l_mom) then
      do i = 1, row_length
        do k = 1, nlayers
          u_p(i,1,k) = u_in_w3(map_w3(1,i) + k-1)
          v_p(i,1,k) = v_in_w3(map_w3(1,i) + k-1)
          u_conv(i,1,k) = u_in_w3_star(map_w3(1,i) + k-1)
          v_conv(i,1,k) = v_in_w3_star(map_w3(1,i) + k-1)
        end do ! k
      end do
    end if

    ! Map tracer fields to UM tracer array
    if ( outer == outer_iterations .and. l_tracer ) then

      allocate(tot_tracer( row_length, 1, nlayers, ntra_fld ))

      do n = 1, ntra_fld
        select case(ukca_tracer_names(n))
        case(fldname_o3p)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( o3p( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_o1d)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( o1d( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_o3)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( o3( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_n)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( nit( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_no)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( no( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_no3)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( no3( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_lumped_n)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( lumped_n( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_n2o5)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( n2o5( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_ho2no2)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( ho2no2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_hono2)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( hono2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_h2o2)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( h2o2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_ch4)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( ch4( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_co)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( co( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_hcho)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( hcho( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_meoo)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( meoo( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_meooh)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( meooh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_h)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( h( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_ch2o)
          ! Do nothing- H2O tracer from chemistry not transported.
        case(fldname_oh)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( oh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_ho2)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( ho2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_cl)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( cl( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_cl2o2)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( cl2o2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_clo)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( clo( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_oclo)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( oclo( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
            end do
        case(fldname_br)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( br( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_lumped_br)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( lumped_br( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_brcl)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( brcl( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_brono2)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( brono2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_n2o)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( n2o( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_lumped_cl)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( lumped_cl( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_hocl)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( hocl( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_hbr)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( hbr( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_hobr)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( hobr( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_clono2)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( clono2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_cfcl3)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( cfcl3( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_cf2cl2)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( cf2cl2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_mebr)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( mebr( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_hono)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( hono( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_c2h6)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( c2h6( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_etoo)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( etoo( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_etooh)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( etooh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_mecho)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( mecho( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_meco3)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( meco3( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_pan)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( pan( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
            end do
        case(fldname_c3h8)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( c3h8( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_n_proo)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( n_proo( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_i_proo)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( i_proo( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_n_prooh)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( n_prooh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_i_prooh)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( i_prooh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_etcho)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( etcho( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_etco3)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( etco3( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_me2co)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( me2co( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_mecoch2oo)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( mecoch2oo( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_mecoch2ooh)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( mecoch2ooh( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_ppan)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( ppan( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_meono2)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( meono2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_c5h8)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( c5h8( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_iso2)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( iso2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_isooh)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( isooh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_ison)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( ison( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_macr)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( macr( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_macro2)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( macro2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_macrooh)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( macrooh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_mpan)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( mpan( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_hacet)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( hacet( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_mgly)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( mgly( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_nald)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( nald( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_hcooh)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( hcooh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_meco3h)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( meco3h( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_meco2h)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( meco2h( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_h2)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( h2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_meoh)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( meoh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_msa)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( msa( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_nh3)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( nh3( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_cs2)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( cs2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_csul)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( csul( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_h2s)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( h2s( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_so3)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( so3( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_passive_o3)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( passive_o3( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_age_of_air)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( age_of_air( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_dms)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( dms( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_so2)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( so2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_h2so4)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( h2so4( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_dmso)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( dmso( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_monoterpene)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( monoterpene( map_wth(1,i)+1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_secondary_organic)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( secondary_organic( map_wth(1,i)+1:map_wth(1,i) + nlayers ),  &
                  r_um )
          end do
        case(fldname_n_nuc_sol)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( n_nuc_sol( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_nuc_sol_su)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( nuc_sol_su( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_nuc_sol_om)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( nuc_sol_om( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_n_ait_sol)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( n_ait_sol( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_ait_sol_su)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( ait_sol_su( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_ait_sol_bc)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( ait_sol_bc( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_ait_sol_om)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( ait_sol_om( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_n_acc_sol)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( n_acc_sol( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_acc_sol_su)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( acc_sol_su( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_acc_sol_bc)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( acc_sol_bc( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_acc_sol_om)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( acc_sol_om( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_acc_sol_ss)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( acc_sol_ss( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_acc_sol_du)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) = 0.0_r_um ! no prognostic, always zero
          end do
        case(fldname_n_cor_sol)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( n_cor_sol( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_cor_sol_su)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( cor_sol_su( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_cor_sol_bc)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( cor_sol_bc( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_cor_sol_om)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( cor_sol_om( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_cor_sol_ss)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( cor_sol_ss( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_cor_sol_du)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) = 0.0_r_um ! no prognostic, always zero
          end do
        case(fldname_n_ait_ins)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( n_ait_ins( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_ait_ins_bc)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( ait_ins_bc( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_ait_ins_om)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( ait_ins_om( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_n_acc_ins)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( n_acc_ins( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_acc_ins_du)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( acc_ins_du( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_n_cor_ins)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( n_cor_ins( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case(fldname_cor_ins_du)
          do i = 1, row_length
            tot_tracer( i, 1, :, n ) =                                         &
            real( cor_ins_du( map_wth(1,i) + 1:map_wth(1,i) + nlayers ), r_um )
          end do
        case default
          write( log_scratch_space, '(A,A)' )                                  &
                 'Missing required UKCA tracer field: ', ukca_tracer_names(n)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end select
      end do

    end if  ! outer == outer_iterations .AND. l_tracer

    ! Currently unused but needed for Comorph-B
    if (prog_tnuc) then
      ! Use tnuc from LFRic and map onto tnuc_new for UM to be passed to conv
      do i = 1, row_length
        do k = 1, nlayers
          tnuc_new(i,1,k) = real(tnuc(map_wth(1,i) + k),kind=r_um) + zerodegc
        end do ! k
      end do
    end if

    !----------------------------------------------------------------
    ! 1) Set timestepping, segmenting, array-dimension and moist
    !    thermodynamics constants for CoMorph consistent with the UM
    !----------------------------------------------------------------

    ! Set stuff in the CoMorph constants module, if not already set:
    if ( .not. l_init_constants ) then
      ! Comorph sets l_init_constants to true the first time it is called,
      ! so this block of code will only be entered before the first
      ! call to comorph, i.e. on the first timestep.

      ! Call routine to set comorph's switches and constants
      call set_constants_from_um( n_conv_levels, ntra_fld, 1 )

    end if

    !----------------------------------------------------------------
    ! 2) Save values of fields before convection, for use in
    !    calculating convective increments needed elsewhere in the UM
    !----------------------------------------------------------------

    if ( .not. l_conv_inc_w ) then
      ! Allocate temporary array for w passed into comorph if not updating w
      allocate( w_work( row_length, rows, 1:nlayers ) )
    else
      ! Minimal allocation if not used
      allocate( w_work(1,1,1) )
    end if

    ! set r_w to zero. The predictor for w is not good, so better to use
    ! start of timestep value for w
    r_w = 0.0_r_def

    call calc_conv_incs  ( i_call_save_before_conv,                            &
                           l_conv_inc_w, z_theta, z_rho,                       &
                           u_p, v_p, u_conv, v_conv, w, w_work,                &
                           u_th_n, v_th_n, u_th_np1, v_th_np1,                 &
                           theta_conv, q_conv, qcl_conv, qcf_conv,             &
                           qcf2_conv, qrain_conv, qgraup_conv,                 &
                           cf_liquid_conv, cf_frozen_conv, bulk_cf_conv,       &
                           dubydt_p, dvbydt_p, r_w,                            &
                           dtheta_conv, q_inc, qcl_inc, qcf_inc,               &
                           qcf2_inc, qrain_inc, qgraup_inc,                    &
                           cf_liquid_inc, cf_frozen_inc, bulk_cf_inc )

    !----------------------------------------------------------------
    ! 3) Conversions of input fields
    !----------------------------------------------------------------
    ! Convert potential temperature to actual temperature
    do k = 1, nlayers
      do i = 1, row_length
        temperature_n(i,1,k) = theta_n(map_wth(1,i)+k)                       &
                           * exner_theta_levels(i,1,k)
        theta_conv(i,1,k) = theta_conv(i,1,k)                                &
                        * exner_theta_levels(i,1,k)
      end do
    end do

    ! Values of ccw0 input are (as treated by the UM) in-cloud
    ! convective cloud-water = grid-mean q_cl_conv / cca.
    ! But CoMorph works using the grid-mean convective cloud-water.
    ! => Convert ccw0 to grid-mean
    do i = 1, row_length
      do k = 1, nlayers
        cca_3d(i,1,k) = cca(map_wth(1,i) + k)
        ccw_3d(i,1,k) = ccw(map_wth(1,i) + k) * cca(map_wth(1,i) + k)
      end do
    end do

    ! The "latest" (_star) fields used by convection are values
    ! interpolated to departure points by SL advection.
    ! The interpolation is not guaranteed to preserve consistency
    ! between the cloud fraction and cloud water fields.
    ! However, various things can go wrong within CoMorph if they
    ! are inconsistent; especially the routine calc_env_partitions,
    ! which attempts to calculate the internal properties of
    ! the in-cloud and clear sub-regions of the grid-box.
    ! E.g. a common problem is if qcf is large but CFF is small,
    ! the in-cloud qcf can get implausibly large, which can
    ! cause the phase-change calculations to yield nonsense.
    ! Therefore, apply safety checks to the cloud-fractions
    ! before passing them into convection...
    do i = 1, row_length
      do k = 1, nlayers
        precfrac_star(i,1,k) = precfrac(map_wth(1,i)+k)
      end do
    end do
    CALL fracs_consistency  ( qcl_conv, qcf_conv, qcf2_conv,                   &
                              qrain_conv, qgraup_conv,                         &
                              cf_liquid_conv, cf_frozen_conv, bulk_cf_conv,    &
                              precfrac_star )

    ! Note: if not using PC2, the cloud-fraction _star fields do exist
    ! but are just set to values after slow_physics
    ! Use the PC2 option l_cloud_call_b4_conv to ensure the latest
    ! cloud fractions are actually calculated before this point.

    ! If prognostic precip fraction is in use:
    if ( l_mcr_qrain .and. l_mcr_precfrac ) then
      ! Allocate array to store precip mixing ratio before convection
      allocate( q_prec_b4(row_length, rows, nlayers))
      ! Save precip mass before convection, for use in updating precfrac later
      call conv_update_precfrac( i_call_save_before_conv,                      &
                                 qrain_conv, qgraup_conv,                      &
                                 frac_bulk_conv, q_prec_b4, precfrac_star )
    end if

    ! For conservation purposes on theta-levels, the bottom rho-level
    ! is actually the surface (so that the bottom theta-level extends
    ! from the surface up to the next rho-level).
    ! Set the actual height of the bottom rho-level
    ! to zero for the remaining calculations:
    do i = 1, row_length
      z_rho(i,1,1) = 0.0_r_um
    end do
    ! Note: this is correct for the interpolation of turbulence
    ! fields below, since ftl(:,:,1) is actually at the surface,
    ! while the rest of the array is on rho-levels.

    ! Set an "effective" boundary-layer height, below-which comorph's
    ! increments are modified to avoid double-counting the boundary-layer
    ! scheme's non-local fluxes
    do i = 1, row_length
      ! Top height of surface-driven mixing is the surface-driven non-local
      ! BL-top height zhnl plus the inversion thickness dzh
      ! (note dzh defaults to a large negative number when not used
      !  due to the inversion being sub-grid; using MAX to remove these).
      ! Also sometimes the local mixing based BL-top height exceeds zhnl and
      ! replaces it, and is stored in zh in that case, so use zh if it is higher
      zh_eff(i,1) = max( zhnl(i,1) + max( dzh(i,1), 0.0_r_um ), zh(i,1) )
      if ( bl_type_7(i,1) == 1.0_r_um ) then
        ! For shear-dominated layers, the Sc-top is essentially coupled to the
        ! surface mixed-layer by shear-driven turbulence.  In this case, use
        ! the Sc-top height.  Comorph suppresses the increments from convection
        ! that only occurs within the layer below zh_eff, so this effectively
        ! disables convection within the shear-dominated layer.
        zh_eff(i,1) = max( zh_eff(i,1), zhsc(i,1) )
      end if
    end do

    ! If using turbulence-based parcel perturbations, need to convert
    ! the turbulence fields into the required units, and interpolate
    ! some of them onto rho-levels
    if ( l_turb_par_gen ) then

      ! Allocate arrays for momentum diffusivities and fluxes
      ! interpolated onto rho-levels, and turbulence lengthscale
      allocate( w_var_rh ( row_length, rows, 1:bl_levels ) )
      allocate( fu_rh    ( row_length, rows, 1:bl_levels ) )
      allocate( fv_rh    ( row_length, rows, 1:bl_levels ) )
      allocate( turb_len ( row_length, rows, 1:bl_levels ) )

      do i = 1, row_length
        do k = 1, bl_levels
          ! level indexing here is confusing but I think correct
          bl_w_var(i,1,k) = wvar(map_wth(1,i)+k)
          rhokm(i,1,k-1) = rhokm_bl(map_wth(1,i) + k-1)
          ftl(i,1,k) = heat_flux(map_w3(1,i) + k-1)
          fqw(i,1,k) = moist_flux(map_w3(1,i) + k-1)
          ! These fields are actually in wtheta, but stored as w3 fields to
          ! match the w2 versions which should really be in the fd1 space
          taux_p(i,1,k-1) = taux(map_w3(1,i)+k-1)
          tauy_p(i,1,k-1) = tauy(map_w3(1,i)+k-1)
        end do
      end do
      do i = 1, row_length
        ! Element 6 of this super array set in jules_exp_kernel contains the
        ! surface buoyancy flux
        fb_surf(i,1) = surf_interp(map_surf(1,i)+6)
        u_s(i,1) = ustar(map_2d(1,i))
        ls_rain(i,1) = ls_rain_2d(map_2d(1,i))
        ls_snow(i,1) = ls_snow_2d(map_2d(1,i))
        delta_x(i,1) = delta(map_wth(1,i))
      end do

      ! Interpolate momentum diffusivity and fluxes onto rho-levels
      call interp_turb( z_rho, z_theta,                                        &
                    bl_w_var, fb_surf, zh, u_s, taux_p, tauy_p,                &
                    w_var_rh, fu_rh, fv_rh )

      do k = 1, bl_levels
        do i = 1, row_length
          ! Convert wind stresses -rho <w'u'> / N m-2
          ! into just <w'u'> / m2 s-2
          fu_rh(i,1,k) = -fu_rh(i,1,k) / rho_wet(i,1,k)
          fv_rh(i,1,k) = -fv_rh(i,1,k) / rho_wet(i,1,k)
          ! Convert sensible heat flux rho_dry <w'Tl'> / K kg m-2
          ! into just <w'Tl'>
          ftl(i,1,k) = ftl(i,1,k) / rho_dry(i,1,k)
          ! Convert total water flux rho_dry <w'q'> / kg m-2
          ! into just <w'q'>
          ! (where q denotes water mixing-ratio)
          fqw(i,1,k) = fqw(i,1,k) / rho_dry(i,1,k)
          ! Note: dividing by dry-density here instead of wet, so
          ! that the flux is expressed in terms of mixing-ratio
          ! rather than specific humidity
        end do
      end do

      ! Calculate turb_len and a scaling factor applied to parcel initial radius
      call calc_turb_len( zh_eff, z_theta, z_rho, rho_wet_tq, qv_n,            &
                          rhokm, bl_w_var, ls_rain, ls_snow, w,                &
                          delta_x, delta_x,                                    &
                          turb_len, par_radius_amp_um )

      ! Check for instances of fluxes too big relative to the turbulent
      ! w-variance (causes excessive parcel perturbations);
      ! increase the w-variance where needed to avoid the problem
      call limit_turb_perts( z_theta, z_rho, p_theta_levels, theta_conv,       &
                             ftl, fqw, fu_rh, fv_rh, w_var_rh )

    end if  ! ( l_turb_par_gen )

    !----------------------------------------------------------------
    ! 4) Initialise temporary arrays needed by comorph
    !----------------------------------------------------------------

    ! Throw an error if any water species are switched on in the UM but
    ! switched off in comorph, since in this case comorph will not
    ! transport them consistently (e.g. qcf2 needs to be transported
    ! along with cf_frozen if it is used).
    if ( ( .not. l_cv_cf ) .or.                                                &
         ! Ice-cloud is always on in the UM, so must be on in comorph
         ( l_mcr_qcf2 .and. ( .not. ( l_cv_cf  .and. l_cv_snow ) ) ) .or.      &
         ! UM 2nd ice category requires both "ice" and "snow" in comorph
         ( l_mcr_qrain .and. ( .not. l_cv_rain ) ) .or.                        &
         ! UM prognostic rain requires rain to be on in comorph
         ( l_mcr_qgraup .and. ( .not. l_cv_graup ) ) ) then
      ! UM prognostic graupel requires graupel to be on in comorph
      call raise_fatal( routinename,                                           &
         "At least one condensed water species is switched "  //               &
         "on in the UM but switched off in CoMorph." //newline//               &
         "Code has not yet been implemented to handle this "  //               &
         "combination consistently." )
    end if

    ! Temporary fields needed for water species if they are switched
    ! on in comorph, but switched off in the UM...
    l_temporary_rain = l_cv_rain .and. (.not. l_mcr_qrain)
    l_temporary_snow = l_cv_snow .and. (.not. l_mcr_qcf2)
    l_temporary_graup = l_cv_graup .and. (.not. l_mcr_qgraup)

    ! In this case, CoMorph will operate using an initial rain / snow / graupel
    ! mixing ratio which is zero everywhere.  Any rain / snow / graupel
    ! produced by CoMorph will be either rained-out this timestep or added to
    ! the cloud-water, to let the large-scale microphysics rain it out.

    ! Allocate and initialise temporary water species to zero if needed...
    if ( l_temporary_rain ) then
      allocate( q_rain_work(row_length,rows,nlayers))
      q_rain_work = 0.0_r_um
    end if
    if ( l_temporary_snow ) then
      allocate( q_snow_work(row_length,rows,nlayers))
      q_snow_work = 0.0_r_um
    end if
    if ( l_temporary_graup ) then
      allocate( q_graup_work(row_length,rows,nlayers))
      q_graup_work = 0.0_r_um
    end if

    ! Initialise local array for bulk conv cloud fraction to zero
    frac_bulk_conv = 0.0_r_um

    !----------------------------------------------------------------
    ! 5) Assign pointers to fields to pass into comorph
    !----------------------------------------------------------------

    ! All the array arguments to this routine are passed into comorph_ctl,
    ! but we pass them in via pointers contained in the derived-type
    ! structures grid, turb, cloudfracs, fields_n, fields_np1.
    ! This routine assigns the pointers to the arrays.
    CALL assign_fields  ( z_theta, z_rho, p_theta_levels, p_rho_levels,        &
                          r_theta_levels,                                      &
                          rho_dry_tq, w_var_rh, ftl, fqw, fu_rh, fv_rh,        &
                          turb_len, par_radius_amp_um, zh_eff,                 &
                          cca_3d, ccw_3d, frac_bulk_conv,                      &
                          u_th_n, v_th_n, w, temperature_n,                    &
                          qv_n, qcl_n, qcf_n,                                  &
                          qcf2_n, qr_n, qgr_n,                                 &
                          cf_liquid_n, cf_frozen_n, bulk_cf_n,                 &
                          u_th_np1, v_th_np1, r_w, theta_conv,                 &
                          q_conv, qcl_conv, qcf_conv,                          &
                          qcf2_conv, qrain_conv, qgraup_conv,                  &
                          cf_liquid_conv, cf_frozen_conv, bulk_cf_conv,        &
                          precfrac_star, l_conv_inc_w,  l_temporary_snow,      &
                          l_temporary_rain, l_temporary_graup,                 &
                          w_work, q_snow_work, q_rain_work, q_graup_work,      &
                          grid, turb, cloudfracs, fields_n, fields_np1 )

    if (l_tracer) then
      ! Allocate list of pointers to tracer fields
      allocate( fields_np1 % tracers(ntra_fld) )

      do i_field = 1, ntra_fld
        fields_np1 % tracers(i_field)%pt(1:, 1:, 1:)          &
               => tot_tracer(:,:,:,i_field)
      end do

    end if

    !----------------------------------------------------------------
    ! 7) Diagnostic requests
    !----------------------------------------------------------------
    ! CoMorph has an internal diagnostic system; meta-data for
    ! each diagnostic, and pointers to the output arrays, are
    ! contained in the derived-type structure comorph_diags.

    ! Always request the CAPE diagnostic, as it can be used by
    ! other parts of the UM, regardless of whether it is
    ! requested as an output diagostic.
    ! Request the mass-flux-weighted CAPE, not straight CAPE,
    ! since CAPE itself can be noisy / not representative
    ! (e.g. when most of the mass-flux is shallow, but there
    !  is a tiny mass-flux going deep with high CAPE).
    comorph_diags % updraft_diags_2d % mfw_cape                                &
                                % request % x_y = .true.
    comorph_diags % updraft_diags_2d % mfw_cape                                &
                                % field_2d => cape_dil
    ! Updraft and downdraft mass-fluxes:
    ! CoMorph calculates the mass-fluxes on rho-levels.
    ! If they have been requested on theta-levels, we'll need to
    ! interpolate them onto theta-levels afterwards.
    comorph_diags % updraft % par % massflux_d                                 &
                                % request % x_y_z = .true.
    comorph_diags % updraft % par % massflux_d                                 &
                                % field_3d => up_flux_half
    comorph_diags % dndraft % par % massflux_d                                 &
                                % request % x_y_z = .true.
    comorph_diags % dndraft % par % massflux_d                                 &
                                % field_3d => down_flux_half
    ! Set requests and assign pointers for entrainment and
    ! detrainment diagnostics
    if (outer == outer_iterations) then
      if (.not. associated(entrain_up, empty_real_data) ) then
        allocate(ent_up(row_length,rows,nlayers))
        comorph_diags % updraft % plume_model % ent_mass_d                     &
                                % request % x_y_z = .true.
        comorph_diags % updraft % plume_model % ent_mass_d                     &
                                % field_3d => ent_up
      end if
      if (.not. associated(entrain_down, empty_real_data) ) then
        allocate(ent_down(row_length,rows,nlayers))
        comorph_diags % dndraft % plume_model % ent_mass_d                     &
                                % request % x_y_z = .true.
        comorph_diags % dndraft % plume_model % ent_mass_d                     &
                                % field_3d => ent_down
      end if
      if (.not. associated(detrain_up, empty_real_data) ) then
        allocate(det_up(row_length,rows,nlayers))
        comorph_diags % updraft % plume_model % det_mass_d                     &
                                % request % x_y_z = .true.
        comorph_diags % updraft % plume_model % det_mass_d                     &
                                % field_3d => det_up
      end if
      if (.not. associated(detrain_down, empty_real_data) ) then
        allocate(det_down(row_length,rows,nlayers))
        comorph_diags % dndraft % plume_model % det_mass_d                     &
                                % request % x_y_z = .true.
        comorph_diags % dndraft % plume_model % det_mass_d                     &
                                % field_3d => det_down
      end if
    end if
    if (l_pc2_homog_conv_pressure) then
      allocate(pres_inc_env(row_length,rows,nlayers))
      comorph_diags % pressure_incr_env                                        &
                                % request % x_y_z = .true.
      comorph_diags % pressure_incr_env                                        &
                                % field_3d => pres_inc_env
    end if

    !----------------------------------------------------------------
    ! 8) Call CoMorph
    !----------------------------------------------------------------

    ! The moment you've been waiting for!
    call comorph_ctl( l_tracer, segments,                                      &
                  grid, turb, cloudfracs, fields_n, fields_np1,                &
                  comorph_diags )

    ! Nullify pointers used to pass array arguments through CoMorph
    call grid_nullify( grid )
    call turb_nullify( turb )
    call cloudfracs_nullify( cloudfracs )
    call fields_nullify( fields_n )
    call fields_nullify( fields_np1 )

    ! Deallocate list of pointers to tracer fields
    if ( l_tracer )  deallocate( fields_np1 % tracers )

    !----------------------------------------------------------------
    ! 9) Conversions of output fields back into UM "format"
    !----------------------------------------------------------------

    ! Restore bottom level of z_rho, ready for whatever comes next...
    do i = 1, row_length
      z_rho(i,1,1) = r_rho_levels(i,1,1) - r_theta_levels(i,1,0)
    end do

    ! Convert final temperature back into potential temperature
    do k = 1, nlayers
      do i = 1, row_length
        theta_conv(i,1,k) = theta_conv(i,1,k)                                &
                        / exner_theta_levels(i,1,k)
      end do
    end do

    ! If any condensed water species are on in CoMorph but off
    ! in the UM, scatter any mixing-ratio of these species produced
    ! into the UM's equivalent fields
    if ( l_temporary_rain ) then
      ! Add rain to the liquid cloud if prognostic rain is off
      do k = 1, nlayers
        do i = 1, row_length
          qcl_conv(i,1,k) = qcl_conv(i,1,k) + q_rain_work(i,1,k)
        end do
      end do
    end if

    if ( l_temporary_snow ) then
      ! Add snow to single ice field qcf if 2nd ice category is off
      do k = 1, nlayers
        do i = 1, row_length
          qcf_conv(i,1,k) = qcf_conv(i,1,k) + q_snow_work(i,1,k)
        end do
      end do
    end if

    if ( l_temporary_graup ) then
      ! Add graupel to qcf field if prognostic graupel is off
      do k = 1, nlayers
        do i = 1, row_length
          qcf_conv(i,1,k) = qcf_conv(i,1,k) + q_graup_work(i,1,k)
        end do
      end do
    end if

    ! Note: it might be preferable to diagnostically rain-out these
    ! temporary fields instead of adding them to qcl or qcf, but
    ! that functionality hasn't been implemented yet.

    ! Deallocate temporary arrays for hydrometeor species
    if ( l_temporary_graup )  deallocate( q_graup_work )
    if ( l_temporary_snow )   deallocate( q_snow_work )
    if ( l_temporary_rain )   deallocate( q_rain_work )

    ! If prognostic precip fraction is in use:
    if ( l_mcr_qrain .and. l_mcr_precfrac ) then
      ! Update the precip fraction using the convective rain and graupel incs
      call conv_update_precfrac( i_call_diff_to_get_incs,                      &
                                 qrain_conv, qgraup_conv,                      &
                                 frac_bulk_conv, q_prec_b4, precfrac_star )
      ! Deallocate saved precip mass before convection, now we're done
      deallocate( q_prec_b4 )
    end if

    ! Note: currently CoMorph does not compute any cloud fraction
    ! increment associated with transfer of condensate by
    ! precipitation.  If the parcel contains ice, some of which
    ! falls out into the environment, but detrainment is zero,
    ! this will yield zero ice-cloud fraction but non-zero ice
    ! mixing-ratio, which shouldn't be allowed.
    ! Pending implementation of code in CoMorph to compute the
    ! appropriate cloud-fraction associated with precip fall,
    ! here there is a temporary fix to avoid the problem.
    ! A minimum limit is applied to the cloud-fractions after
    ! convection, proportional to the mixing-ratio.
    call fracs_consistency  ( qcl_conv, qcf_conv, qcf2_conv,                   &
                              qrain_conv, qgraup_conv,                         &
                              cf_liquid_conv, cf_frozen_conv, bulk_cf_conv,    &
                              precfrac_star )

    do i = 1, row_length
      do k = 1, nlayers
        precfrac(map_wth(1,i)+k) = precfrac_star(i,1,k)
      end do
      precfrac(map_wth(1,i)+0) = precfrac(map_wth(1,i)+1)
    end do

    ! CoMorph outputs grid-mean convective cloud-water, but the UM expects
    ! ccw to be the in-cloud water-content.
    ! => Convert to in-cloud value
    do k = 1, nlayers
      do i = 1, row_length
        if ( ccw_3d(i,1,k) > 0.0_r_um ) ccw_3d(i,1,k) = ccw_3d(i,1,k) / cca_3d(i,1,k)
      end do
    end do

    ! If using turbulence-based parcel perturbations
    if ( l_turb_par_gen ) then
      ! Deallocate the interpolated momentum diffusivity and fluxes
      ! on rho-levels
      deallocate( turb_len )
      deallocate( fv_rh )
      deallocate( fu_rh )
      deallocate( w_var_rh )
    end if  ! ( l_turb_par_gen )

    ! Calculate extra convective cloud fields from cca0, ccw0
    call comorph_conv_cloud_extras(                                            &
             n_conv_levels, rho_dry_tq, rho_wet_tq,                            &
             r_theta_levels, r_rho_levels,                                     &
             cca_3d, ccw_3d, frac_bulk_conv,                                   &
             cclwp0, ccb0, cct0, lcbase0, cca_3d0, ccw_3d0,                    &
             cclwp, cca_2d_loc, lcca, ccb, cct, lcbase, lctop )

    !----------------------------------------------------------------
    ! 10) Calculate convective increments required elsewhere in the UM
    !----------------------------------------------------------------
    call calc_conv_incs  ( i_call_diff_to_get_incs,                            &
                           l_conv_inc_w, z_theta, z_rho,                       &
                           u_p, v_p, u_conv, v_conv, w, w_work,                &
                           u_th_n, v_th_n, u_th_np1, v_th_np1,                 &
                           theta_conv, q_conv, qcl_conv, qcf_conv,             &
                           qcf2_conv, qrain_conv, qgraup_conv,                 &
                           cf_liquid_conv, cf_frozen_conv, bulk_cf_conv,       &
                           dubydt_p, dvbydt_p, r_w,                            &
                           dtheta_conv, q_inc, qcl_inc, qcf_inc,               &
                           qcf2_inc, qrain_inc, qgraup_inc,                    &
                           cf_liquid_inc, cf_frozen_inc, bulk_cf_inc )

    ! Deallocate temporary work array for w
    deallocate( w_work )

    ! single level convection diagnostics
    do i = 1, row_length
      cca_2d(map_2d(1,i)) = cca_2d_loc(i,1)
      cape_diluted(map_2d(1,i)) = cape_dil(i,1)
    end do

    if (outer == outer_iterations) then
      if (.not. associated(lowest_cca_2d, empty_real_data) ) then
        do i = 1, row_length
          lowest_cca_2d(map_2d(1,i)) = lcca(i,1)
        end do
      end if
    end if ! outer_iterations

    ! update input fields
    do k = 1, nlayers
      do i = 1, row_length
        dt_conv(map_wth(1,i) + k)   = dtheta_conv(i,1,k) &
                                    * exner_theta_levels(i,1,k)
        dmv_conv(map_wth(1,i) + k)  = q_inc(i,1,k)
        dmcl_conv(map_wth(1,i) + k) = qcl_inc(i,1,k)
        dms_conv(map_wth(1,i) + k) = qcf_inc(i,1,k)
      end do
    end do
    if (microphysics_casim) then
      do i = 1, row_length
        do k = 1, nlayers
          m_ci(map_wth(1,i) + k) = qcf2_conv(i,1,k)
        end do
        m_ci(map_wth(1,i) + 0) = m_ci(map_wth(1,i) + 1)
      end do
    end if
    if (l_mcr_qrain) then
      do i = 1, row_length
        do k = 1, nlayers
          m_r(map_wth(1,i) + k) = qrain_conv(i,1,k)
        end do
        m_r(map_wth(1,i) + 0) = m_r(map_wth(1,i) + 1)
      end do
    end if
    if (l_mcr_qgraup) then
      do i = 1, row_length
        do k = 1, nlayers
          m_g(map_wth(1,i) + k) = qgraup_conv(i,1,k)
        end do
        m_g(map_wth(1,i) + 0) = m_g(map_wth(1,i) + 1)
      end do
    end if

    ! Update diagnostics
    up_flux_half(:,1,nlayers) = 0.0_r_def
    if (.not. associated(massflux_up_half, empty_real_data) ) then
      do k = 1, nlayers
        do i = 1, row_length
          ! map_w3(1) points to level 1
          ! Scale output flux on rho-levels by g to convert to Pa s-1
          massflux_up_half(map_w3(1,i) + k-1) = up_flux_half(i,1,k) * g
          if (k <= lcbase(i,1) .or. lcbase(i,1) == 0) &
               massflux_up_half(map_w3(1,i) + k-1) = 0.0_r_def
        end do
      end do
    end if

    down_flux_half(:,1,nlayers) = 0.0_r_def
    do k = 1, n_conv_levels
      do i = 1, row_length
        ! Interpolate flux onto theta-levels
        ! Scale output flux on rho-levels by g to convert to Pa s-1
        interp = ( z_theta(i,1,k) - z_rho(i,1,k) )                             &
                 / ( z_rho(i,1,k+1) - z_rho(i,1,k) )
        massflux_up(map_wth(1,i) + k) = (1.0_r_def-interp) *                   &
                                       up_flux_half(i,1,k) * g                 &
                                    +  interp  * up_flux_half(i,1,k+1) * g
        if (k <= lcbase(i,1) .or. lcbase(i,1) == 0) &
             massflux_up(map_wth(1,i) + k) = 0.0_r_def
        massflux_down(map_wth(1,i) + k) = (1.0_r_def-interp) *                 &
                                       down_flux_half(i,1,k) * g               &
                                    +  interp  * down_flux_half(i,1,k+1) * g
      end do
    end do

    if (l_pc2_homog_conv_pressure) then
      do k = 1, n_conv_levels
        do i = 1, row_length
          pressure_inc_env(map_wth(1,i) + k) = pres_inc_env(i,1,k)
        end do
      end do
    end if

    ! Update optional diagnostics
    if (outer == outer_iterations) then
      if (.not. associated(entrain_up, empty_real_data) ) then
        do k = 1, n_conv_levels
          do i = 1, row_length
            ! Convert to Pa s-1
            entrain_up(map_wth(1,i) + k) = ent_up(i,1,k) * g
            if (k <= lcbase(i,1) .or. lcbase(i,1) == 0) &
                 entrain_up(map_wth(1,i) + k) = 0.0_r_def
          end do
        end do
        deallocate(ent_up)
      end if
      if (.not. associated(entrain_down, empty_real_data) ) then
        do k = 1, n_conv_levels
          do i = 1, row_length
            ! Convert to Pa s-1
            entrain_down(map_wth(1,i) + k) = ent_down(i,1,k) * g
          end do
        end do
        deallocate(ent_down)
      end if
      if (.not. associated(detrain_up, empty_real_data) ) then
        do k = 1, n_conv_levels
          do i = 1, row_length
            ! Convert to Pa s-1
            detrain_up(map_wth(1,i) + k) = det_up(i,1,k) * g
            if (k <= lcbase(i,1) .or. lcbase(i,1) == 0) &
                 detrain_up(map_wth(1,i) + k) = 0.0_r_def
          end do
        end do
        deallocate(det_up)
      end if
      if (.not. associated(detrain_down, empty_real_data) ) then
        do k = 1, n_conv_levels
          do i = 1, row_length
            ! Convert to Pa s-1
            detrain_down(map_wth(1,i) + k) = det_down(i,1,k) * g
          end do
        end do
        deallocate(det_down)
      end if
    end if ! outer_iterations

    if (l_mom) then
      do k = 1, n_conv_levels
        do i = 1, row_length
          ! total increments
          du_conv(map_w3(1,i) + k -1) = dubydt_p(i,1,k) * timestep
          dv_conv(map_w3(1,i) + k -1) = dvbydt_p(i,1,k) * timestep
        end do
      end do
    end if    !l_mom

    ! Update the time-smoothed convection prognostics
    if (l_conv_prog_dtheta) then
      decay_amount = timestep / tau_conv_prog_dtheta
      do k = 1, n_conv_levels
        do i = 1, row_length
          dt_conv(map_wth(1,i) + k)  = (decay_amount * dtheta_conv(i,1,k) &
              + (1.0_r_um - decay_amount) * conv_prog_dtheta(map_wth(1,i) + k))&
               * exner_in_wth(map_wth(1,i) + k)
        end do
      end do
      if (outer == outer_iterations) then
        do i = 1, row_length
          do k = 1, n_conv_levels
            conv_prog_dtheta(map_wth(1,i) + k) = dt_conv(map_wth(1,i) + k)     &
                                           / exner_in_wth(map_wth(1,i) + k)
          end do
          conv_prog_dtheta(map_wth(1,i)+0) = conv_prog_dtheta(map_wth(1,i) + 1)
        end do
      end if
    end if

    if (l_conv_prog_dq) then
      decay_amount = timestep / tau_conv_prog_dq
      do i = 1, row_length
        do k = 1, n_conv_levels
          dmv_conv(map_wth(1,i)+k) = decay_amount * dmv_conv(map_wth(1,i) + k) &
                 + (1.0_r_um - decay_amount) * conv_prog_dmv(map_wth(1,i) + k)
        end do
      end do
      if (outer == outer_iterations) then
        do i = 1, row_length
          do k = 1, n_conv_levels
            conv_prog_dmv(map_wth(1,i) + k) = dmv_conv(map_wth(1,i) + k)
          end do
        end do
        conv_prog_dmv(map_wth(1,i) + 0) = conv_prog_dmv(map_wth(1,i) + 1)
      end if
    end if

    ! Store cloud fraction increments for adding on later if using PC2
    do i = 1, row_length
      do k = 1, n_conv_levels
        dcfl_conv(map_wth(1,i) + k) = cf_liquid_conv(i,1,k) - cf_liq(map_wth(1,i) + k)
        dcff_conv(map_wth(1,i) + k) = cf_frozen_conv(i,1,k) - cf_ice(map_wth(1,i) + k)
        dbcf_conv(map_wth(1,i) + k) = bulk_cf_conv(i,1,k)   - cf_bulk(map_wth(1,i) + k)
      end do
      dcfl_conv(map_wth(1,i) + 0) = dcfl_conv(map_wth(1,i) + 1)
      dcff_conv(map_wth(1,i) + 0) = dcff_conv(map_wth(1,i) + 1)
      dbcf_conv(map_wth(1,i) + 0) = dbcf_conv(map_wth(1,i) + 1)
    end do

    ! Set level 0 increment such that theta increment will equal level 1
    do i = 1, row_length
      dt_conv (map_wth(1,i) + 0) = dt_conv  (map_wth(1,i) + 1)    &
                             * exner_in_wth(map_wth(1,i) + 0) &
                             / exner_in_wth(map_wth(1,i) + 1)
      dmv_conv (map_wth(1,i) + 0) = dmv_conv (map_wth(1,i) + 1)
      dmcl_conv(map_wth(1,i) + 0) = dmcl_conv(map_wth(1,i) + 1)
      dms_conv(map_wth(1,i) + 0) = dms_conv(map_wth(1,i) + 1)
    end do

    ! Store convective downdraught mass fluxes at cloud base
    ! if required for surface exchange.
    if (srf_ex_cnv_gust == ip_srfexwithcnv) then
      do i = 1, row_length
        if (ccb(i,1) > 0) then
          dd_mf_cb(map_2d(1,i))=massflux_down( map_wth(1,i) + ccb(i,1))
        else
          dd_mf_cb(map_2d(1,i))=0.0_r_def
        end if
      end do
    end if

    ! copy convective cloud fraction into prognostic array
    do k = 1, n_conv_levels
      do i = 1, row_length
        cca(map_wth(1,i) + k) =  min(cca_3d(i,1,k), 1.0_r_um)
        ccw(map_wth(1,i) + k) =  ccw_3d(i,1,k)
      end do
    end do

    if (outer == outer_iterations) then
     ! Copy integers into real diagnostic arrays
     if (.not. associated(cv_top, empty_real_data) ) then
       do i = 1, row_length
         cv_top(map_2d(1,i))        = real(cct(i,1))
       end do
      end if
      if (.not. associated(cv_base, empty_real_data) ) then
        do i = 1, row_length
          cv_base(map_2d(1,i))       = real(ccb(i,1))
        end do
      end if
      if (.not. associated(lowest_cv_top, empty_real_data) ) then
        do i = 1, row_length
          lowest_cv_top(map_2d(1,i)) = real(lctop(i,1))
        end do
      end if
      if (.not. associated(lowest_cv_base, empty_real_data) ) then
        do i = 1, row_length
          lowest_cv_base(map_2d(1,i)) = real(lcbase(i,1))
        end do
      end if

      ! pressure at cv top/base
      if (.not. associated(pres_cv_top, empty_real_data) ) then
        do i = 1, row_length
          if (cct(i,1) > 0) then
            pres_cv_top(map_2d(1,i)) = p_rho_levels(i,1,cct(i,1))
          else
            pres_cv_top(map_2d(1,i)) = 0.0_r_def
          end if
        end do
      end if
      if (.not. associated(pres_cv_base, empty_real_data) ) then
        do i = 1, row_length
          if (ccb(i,1) > 0) then
            pres_cv_base(map_2d(1,i)) = p_rho_levels(i,1,ccb(i,1))
          else
            pres_cv_base(map_2d(1,i))= 0.0_r_def
          end if
        end do
      end if

      ! pressure at lowest cv top/base
      if (.not. associated(pres_lowest_cv_top, empty_real_data) ) then
        do i = 1, row_length
          if (lctop(i,1) > 0) then
            pres_lowest_cv_top(map_2d(1,i)) = p_rho_levels(i,1,lctop(i,1))
          else
            pres_lowest_cv_top(map_2d(1,i)) = 0.0_r_def
        end if
      end do
      end if
      if (.not. associated(pres_lowest_cv_base, empty_real_data) ) then
        do i = 1, row_length
          if (lcbase(i,1) > 0) then
            pres_lowest_cv_base(map_2d(1,i)) = p_rho_levels(i,1,lcbase(i,1))
          else
            pres_lowest_cv_base(map_2d(1,i))= 0.0_r_def
          end if
        end do
      end if

      ! provide some estimate of TKE in convective plumes, based on
      ! the mass flux and convective cloud area
      do i = 1, row_length
        do k = 1, bl_levels
          tke_bl(map_wth(1,i)+k) = min(max_tke,max(tke_bl(map_wth(1,i)+k),     &
               ( massflux_up(map_wth(1,i)+k) / ( g*rho_wet_tq(i,1,k)*          &
                 min(0.5_r_um,max(0.05_r_um,cca_2d(map_2d(1,i)))) ) )**2))
                 ! 0.5 and 0.05 are used here as plausible max and min
                 ! values of CCA to prevent numerical problems
        end do
      end do

    end if ! outer_iterations

    ! Copy tracers back to LFRic fields
    if ( outer == outer_iterations .and. l_tracer ) then
      do n = 1, ntra_fld
        select case(ukca_tracer_names(n))
         case(fldname_o3p)
           do i = 1, row_length
             o3p( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =               &
                  real( tot_tracer( i, 1, :, n ), r_def )
             o3p( map_wth(1,i) + 0 ) = o3p( map_wth(1,i) + 1 )
           end do
         case(fldname_o1d)
           do i = 1, row_length
             o1d( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =               &
                  real( tot_tracer( i, 1, :, n ), r_def )
             o1d( map_wth(1,i) + 0 ) = o1d( map_wth(1,i) + 1 )
           end do
         case(fldname_o3)
           do i = 1, row_length
             o3( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =                &
                  real( tot_tracer( i, 1, :, n ), r_def )
             o3( map_wth(1,i) + 0 ) = o3( map_wth(1,i) + 1 )
           end do
         case(fldname_n)
           do i = 1, row_length
             nit( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =               &
                  real( tot_tracer( i, 1, :, n ), r_def )
             nit( map_wth(1,i) + 0 ) = nit( map_wth(1,i) + 1 )
           end do
         case(fldname_no)
           do i = 1, row_length
             no( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =                &
                  real( tot_tracer( i, 1, :, n ), r_def )
             no( map_wth(1,i) + 0 ) = no( map_wth(1,i) + 1 )
           end do
         case(fldname_no3)
           do i = 1, row_length
             no3( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =               &
                  real( tot_tracer( i, 1, :, n ), r_def )
             no3( map_wth(1,i) + 0 ) = no3( map_wth(1,i) + 1 )
           end do
         case(fldname_lumped_n)
           do i = 1, row_length
             lumped_n( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =          &
                  real( tot_tracer( i, 1, :, n ), r_def )
             lumped_n( map_wth(1,i) + 0 ) = lumped_n( map_wth(1,i) + 1 )
           end do
         case(fldname_n2o5)
           do i = 1, row_length
             n2o5( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             n2o5( map_wth(1,i) + 0 ) = n2o5( map_wth(1,i) + 1 )
           end do
         case(fldname_ho2no2)
           do i = 1, row_length
             ho2no2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =            &
                  real( tot_tracer( i, 1, :, n ), r_def )
             ho2no2( map_wth(1,i) + 0 ) = ho2no2( map_wth(1,i) + 1 )
           end do
         case(fldname_hono2)
           do i = 1, row_length
             hono2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =             &
                  real( tot_tracer( i, 1, :, n ), r_def )
             hono2( map_wth(1,i) + 0 ) = hono2( map_wth(1,i) + 1 )
           end do
        case(fldname_h2o2)
           do i = 1, row_length
             h2o2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             h2o2( map_wth(1,i) + 0 ) = h2o2( map_wth(1,i) + 1 )
           end do
         case(fldname_ch4)
            do i = 1, row_length
              ch4( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                   real( tot_tracer( i, 1, :, n ), r_def )
              ch4( map_wth(1,i) + 0 ) = ch4( map_wth(1,i) + 1 )
            end do
         case(fldname_co)
           do i = 1, row_length
             co( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =                &
                  real( tot_tracer( i, 1, :, n ), r_def )
             co( map_wth(1,i) + 0 ) = co( map_wth(1,i) + 1 )
           end do
         case(fldname_hcho)
           do i = 1, row_length
             hcho( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             hcho( map_wth(1,i) + 0 ) = hcho( map_wth(1,i) + 1 )
           end do
         case(fldname_meoo)
           do i = 1, row_length
             meoo( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             meoo( map_wth(1,i) + 0 ) = meoo( map_wth(1,i) + 1 )
           end do
         case(fldname_meooh)
           do i = 1, row_length
             meooh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =             &
                  real( tot_tracer( i, 1, :, n ), r_def )
             meooh( map_wth(1,i) + 0 ) = meooh( map_wth(1,i) + 1 )
           end do
         case(fldname_h)
           do i = 1, row_length
             h( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =                 &
                  real( tot_tracer( i, 1, :, n ), r_def )
             h( map_wth(1,i) + 0 ) = h( map_wth(1,i) + 1 )
           end do
         case(fldname_ch2o)
         case(fldname_oh)
           do i = 1, row_length
             oh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =                &
                  real( tot_tracer( i, 1, :, n ), r_def )
             oh( map_wth(1,i) + 0 ) = oh( map_wth(1,i) + 1 )
           end do
         case(fldname_ho2)
           do i = 1, row_length
             ho2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =               &
                  real( tot_tracer( i, 1, :, n ), r_def )
             ho2( map_wth(1,i) + 0 ) = ho2( map_wth(1,i) + 1 )
           end do
         case(fldname_cl)
           do i = 1, row_length
             cl( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =                &
                  real( tot_tracer( i, 1, :, n ), r_def )
             cl( map_wth(1,i) + 0 ) = cl( map_wth(1,i) + 1 )
           end do
         case(fldname_cl2o2)
           do i = 1, row_length
             cl2o2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =             &
                  real( tot_tracer( i, 1, :, n ), r_def )
             cl2o2( map_wth(1,i) + 0 ) = cl2o2( map_wth(1,i) + 1 )
           end do
         case(fldname_clo)
           do i = 1, row_length
             clo( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =               &
                  real( tot_tracer( i, 1, :, n ), r_def )
             clo( map_wth(1,i) + 0 ) = clo( map_wth(1,i) + 1 )
           end do
         case(fldname_oclo)
           do i = 1, row_length
             oclo( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             oclo( map_wth(1,i) + 0 ) = oclo( map_wth(1,i) + 1 )
           end do
         case(fldname_br)
           do i = 1, row_length
             br( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =                &
                  real( tot_tracer( i, 1, :, n ), r_def )
             br( map_wth(1,i) + 0 ) = br( map_wth(1,i) + 1 )
           end do
         case(fldname_lumped_br)
           do i = 1, row_length
             lumped_br( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                  real( tot_tracer( i, 1, :, n ), r_def )
             lumped_br( map_wth(1,i) + 0 ) = lumped_br( map_wth(1,i) + 1 )
           end do
         case(fldname_brcl)
           do i = 1, row_length
             brcl( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             brcl( map_wth(1,i) + 0 ) = brcl( map_wth(1,i) + 1 )
           end do
         case(fldname_brono2)
           do i = 1, row_length
             brono2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =            &
                  real( tot_tracer( i, 1, :, n ), r_def )
             brono2( map_wth(1,i) + 0 ) = brono2( map_wth(1,i) + 1 )
           end do
         case(fldname_n2o)
           do i = 1, row_length
             n2o( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =               &
                  real( tot_tracer( i, 1, :, n ), r_def )
             n2o( map_wth(1,i) + 0 ) = n2o( map_wth(1,i) + 1 )
           end do
         case(fldname_lumped_cl)
           do i = 1, row_length
             lumped_cl( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                  real( tot_tracer( i, 1, :, n ), r_def )
             lumped_cl( map_wth(1,i) + 0 ) = lumped_cl( map_wth(1,i) + 1 )
           end do
         case(fldname_hocl)
           do i = 1, row_length
             hocl( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             hocl( map_wth(1,i) + 0 ) = hocl( map_wth(1,i) + 1 )
           end do
         case(fldname_hbr)
           do i = 1, row_length
             hbr( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =               &
                  real( tot_tracer( i, 1, :, n ), r_def )
             hbr( map_wth(1,i) + 0 ) = hbr( map_wth(1,i) + 1 )
           end do
         case(fldname_hobr)
           do i = 1, row_length
             hobr( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             hobr( map_wth(1,i) + 0 ) = hobr( map_wth(1,i) + 1 )
           end do
         case(fldname_clono2)
           do i = 1, row_length
             clono2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =            &
                  real( tot_tracer( i, 1, :, n ), r_def )
             clono2( map_wth(1,i) + 0 ) = clono2( map_wth(1,i) + 1 )
           end do
         case(fldname_cfcl3)
           do i = 1, row_length
             cfcl3( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =             &
                  real( tot_tracer( i, 1, :, n ), r_def )
             cfcl3( map_wth(1,i) + 0 ) = cfcl3( map_wth(1,i) + 1 )
           end do
         case(fldname_cf2cl2)
           do i = 1, row_length
             cf2cl2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =            &
                  real( tot_tracer( i, 1, :, n ), r_def )
             cf2cl2( map_wth(1,i) + 0 ) = cf2cl2( map_wth(1,i) + 1 )
           end do
         case(fldname_mebr)
           do i = 1, row_length
             mebr( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             mebr( map_wth(1,i) + 0 ) = mebr( map_wth(1,i) + 1 )
           end do
         case(fldname_hono)
           do i = 1, row_length
             hono( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             hono( map_wth(1,i) + 0 ) = hono( map_wth(1,i) + 1 )
           end do
         case(fldname_c2h6)
           do i = 1, row_length
             c2h6( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             c2h6( map_wth(1,i) + 0 ) = c2h6( map_wth(1,i) + 1 )
           end do
         case(fldname_etoo)
           do i = 1, row_length
             etoo( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             etoo( map_wth(1,i) + 0 ) = etoo( map_wth(1,i) + 1 )
           end do
         case(fldname_etooh)
           do i = 1, row_length
             etooh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =             &
                  real( tot_tracer( i, 1, :, n ), r_def )
             etooh( map_wth(1,i) + 0 ) = etooh( map_wth(1,i) + 1 )
           end do
         case(fldname_mecho)
           do i = 1, row_length
             mecho( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =             &
                  real( tot_tracer( i, 1, :, n ), r_def )
             mecho( map_wth(1,i) + 0 ) = mecho( map_wth(1,i) + 1 )
           end do
         case(fldname_meco3)
           do i = 1, row_length
             meco3( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =             &
                  real( tot_tracer( i, 1, :, n ), r_def )
             meco3( map_wth(1,i) + 0 ) = meco3( map_wth(1,i) + 1 )
           end do
         case(fldname_pan)
           do i = 1, row_length
             pan( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =               &
                  real( tot_tracer( i, 1, :, n ), r_def )
             pan( map_wth(1,i) + 0 ) = pan( map_wth(1,i) + 1 )
           end do
         case(fldname_c3h8)
           do i = 1, row_length
             c3h8( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             c3h8( map_wth(1,i) + 0 ) = c3h8( map_wth(1,i) + 1 )
           end do
         case(fldname_n_proo)
           do i = 1, row_length
             n_proo( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =            &
                  real( tot_tracer( i, 1, :, n ), r_def )
             n_proo( map_wth(1,i) + 0 ) = n_proo( map_wth(1,i) + 1 )
           end do
         case(fldname_i_proo)
           do i = 1, row_length
             i_proo( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =            &
                  real( tot_tracer( i, 1, :, n ), r_def )
             i_proo( map_wth(1,i) + 0 ) = i_proo( map_wth(1,i) + 1 )
           end do
         case(fldname_n_prooh)
           do i = 1, row_length
             n_prooh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =           &
                  real( tot_tracer( i, 1, :, n ), r_def )
             n_prooh( map_wth(1,i) + 0 ) = n_prooh( map_wth(1,i) + 1 )
           end do
         case(fldname_i_prooh)
           do i = 1, row_length
             i_prooh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =           &
                  real( tot_tracer( i, 1, :, n ), r_def )
             i_prooh( map_wth(1,i) + 0 ) = i_prooh( map_wth(1,i) + 1 )
           end do
         case(fldname_etcho)
           do i = 1, row_length
             etcho( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =             &
                  real( tot_tracer( i, 1, :, n ), r_def )
             etcho( map_wth(1,i) + 0 ) = etcho( map_wth(1,i) + 1 )
           end do
         case(fldname_etco3)
           do i = 1, row_length
             etco3( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =             &
                  real( tot_tracer( i, 1, :, n ), r_def )
             etco3( map_wth(1,i) + 0 ) = etco3( map_wth(1,i) + 1 )
           end do
         case(fldname_me2co)
           do i = 1, row_length
             me2co( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =             &
                  real( tot_tracer( i, 1, :, n ), r_def )
             me2co( map_wth(1,i) + 0 ) = me2co( map_wth(1,i) + 1 )
           end do
         case(fldname_mecoch2oo)
           do i = 1, row_length
             mecoch2oo( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                  real( tot_tracer( i, 1, :, n ), r_def )
             mecoch2oo( map_wth(1,i) + 0 ) = mecoch2oo( map_wth(1,i) + 1 )
           end do
         case(fldname_mecoch2ooh)
           do i = 1, row_length
             mecoch2ooh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =        &
                  real( tot_tracer( i, 1, :, n ), r_def )
             mecoch2ooh( map_wth(1,i) + 0 ) = mecoch2ooh( map_wth(1,i) + 1 )
           end do
         case(fldname_ppan)
           do i = 1, row_length
             ppan( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             ppan( map_wth(1,i) + 0 ) = ppan( map_wth(1,i) + 1 )
           end do
         case(fldname_meono2)
           do i = 1, row_length
             meono2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =            &
                  real( tot_tracer( i, 1, :, n ), r_def )
             meono2( map_wth(1,i) + 0 ) = meono2( map_wth(1,i) + 1 )
           end do
         case(fldname_c5h8)
           do i = 1, row_length
             c5h8( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             c5h8( map_wth(1,i) + 0 ) = c5h8( map_wth(1,i) + 1 )
           end do
         case(fldname_iso2)
           do i = 1, row_length
             iso2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             iso2( map_wth(1,i) + 0 ) = iso2( map_wth(1,i) + 1 )
           end do
         case(fldname_isooh)
           do i = 1, row_length
             isooh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =             &
                  real( tot_tracer( i, 1, :, n ), r_def )
             isooh( map_wth(1,i) + 0 ) = isooh( map_wth(1,i) + 1 )
           end do
         case(fldname_ison)
           do i = 1, row_length
             ison( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             ison( map_wth(1,i) + 0 ) = ison( map_wth(1,i) + 1 )
           end do
         case(fldname_macr)
           do i = 1, row_length
             macr( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             macr( map_wth(1,i) + 0 ) = macr( map_wth(1,i) + 1 )
           end do
         case(fldname_macro2)
           do i = 1, row_length
             macro2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =            &
                  real( tot_tracer( i, 1, :, n ), r_def )
             macro2( map_wth(1,i) + 0 ) = macro2( map_wth(1,i) + 1 )
           end do
         case(fldname_macrooh)
           do i = 1, row_length
             macrooh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =           &
                  real( tot_tracer( i, 1, :, n ), r_def )
             macrooh( map_wth(1,i) + 0 ) = macrooh( map_wth(1,i) + 1 )
           end do
         case(fldname_mpan)
           do i = 1, row_length
             mpan( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             mpan( map_wth(1,i) + 0 ) = mpan( map_wth(1,i) + 1 )
           end do
         case(fldname_hacet)
           do i = 1, row_length
             hacet( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =             &
                  real( tot_tracer( i, 1, :, n ), r_def )
             hacet( map_wth(1,i) + 0 ) = hacet( map_wth(1,i) + 1 )
           end do
         case(fldname_mgly)
           do i = 1, row_length
             mgly( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             mgly( map_wth(1,i) + 0 ) = mgly( map_wth(1,i) + 1 )
           end do
         case(fldname_nald)
           do i = 1, row_length
             nald( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             nald( map_wth(1,i) + 0 ) = nald( map_wth(1,i) + 1 )
           end do
         case(fldname_hcooh)
           do i = 1, row_length
             hcooh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =             &
                  real( tot_tracer( i, 1, :, n ), r_def )
             hcooh( map_wth(1,i) + 0 ) = hcooh( map_wth(1,i) + 1 )
           end do
         case(fldname_meco3h)
           do i = 1, row_length
             meco3h( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =            &
                  real( tot_tracer( i, 1, :, n ), r_def )
             meco3h( map_wth(1,i) + 0 ) = meco3h( map_wth(1,i) + 1 )
           end do
         case(fldname_meco2h)
           do i = 1, row_length
             meco2h( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =            &
                  real( tot_tracer( i, 1, :, n ), r_def )
             meco2h( map_wth(1,i) + 0 ) = meco2h( map_wth(1,i) + 1 )
           end do
         case(fldname_h2)
           do i = 1, row_length
             h2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =                &
                  real( tot_tracer( i, 1, :, n ), r_def )
             h2( map_wth(1,i) + 0 ) = h2( map_wth(1,i) + 1 )
           end do
         case(fldname_meoh)
           do i = 1, row_length
             meoh( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             meoh( map_wth(1,i) + 0 ) = meoh( map_wth(1,i) + 1 )
           end do
         case(fldname_msa)
           do i = 1, row_length
             msa( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =               &
                  real( tot_tracer( i, 1, :, n ), r_def )
             msa( map_wth(1,i) + 0 ) = msa( map_wth(1,i) + 1 )
           end do
         case(fldname_nh3)
           do i = 1, row_length
             nh3( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =               &
                  real( tot_tracer( i, 1, :, n ), r_def )
             nh3( map_wth(1,i) + 0 ) = nh3( map_wth(1,i) + 1 )
           end do
         case(fldname_cs2)
           do i = 1, row_length
             cs2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =               &
                  real( tot_tracer( i, 1, :, n ), r_def )
             cs2( map_wth(1,i) + 0 ) = cs2( map_wth(1,i) + 1 )
           end do
         case(fldname_csul)
           do i = 1, row_length
             csul( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                  real( tot_tracer( i, 1, :, n ), r_def )
             csul( map_wth(1,i) + 0 ) = csul( map_wth(1,i) + 1 )
           end do
         case(fldname_h2s)
           do i = 1, row_length
             h2s( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =               &
                  real( tot_tracer( i, 1, :, n ), r_def )
             h2s( map_wth(1,i) + 0 ) = h2s( map_wth(1,i) + 1 )
           end do
         case(fldname_so3)
           do i = 1, row_length
             so3( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =               &
                  real( tot_tracer( i, 1, :, n ), r_def )
             so3( map_wth(1,i) + 0 ) = so3( map_wth(1,i) + 1 )
           end do
         case(fldname_passive_o3)
           do i = 1, row_length
             passive_o3( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =        &
                  real( tot_tracer( i, 1, :, n ), r_def )
             passive_o3( map_wth(1,i) + 0 ) = passive_o3( map_wth(1,i) + 1 )
           end do
         case(fldname_age_of_air)
           do i = 1, row_length
             age_of_air( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =        &
                  real( tot_tracer( i, 1, :, n ), r_def )
             age_of_air( map_wth(1,i) + 0 ) = age_of_air( map_wth(1,i) + 1 )
           end do
        case(fldname_dms)
          do i = 1, row_length
            dms( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =                &
                 real( tot_tracer( i, 1, :, n ), r_def )
            dms( map_wth(1,i) + 0 ) = dms( map_wth(1,i) + 1 )
          end do
        case(fldname_so2)
          do i = 1, row_length
            so2( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =                &
                 real( tot_tracer( i, 1, :, n ), r_def )
            so2( map_wth(1,i) + 0 ) = so2( map_wth(1,i) + 1 )
          end do
        case(fldname_h2so4)
          do i = 1, row_length
            h2so4( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =              &
                 real( tot_tracer( i, 1, :, n ), r_def )
            h2so4( map_wth(1,i) + 0 ) = h2so4( map_wth(1,i) + 1 )
          end do
        case(fldname_dmso)
          do i = 1, row_length
            dmso( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =               &
                 real( tot_tracer( i, 1, :, n ), r_def )
            dmso( map_wth(1,i) + 0 ) = dmso( map_wth(1,i) + 1 )
          end do
        case(fldname_monoterpene)
          do i = 1, row_length
            monoterpene( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =        &
                 real( tot_tracer( i, 1, :, n ), r_def )
            monoterpene( map_wth(1,i) + 0 ) = monoterpene( map_wth(1,i) + 1 )
          end do
        case(fldname_secondary_organic)
          do i = 1, row_length
            secondary_organic( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =  &
                 real( tot_tracer( i, 1, :, n ), r_def )
            secondary_organic( map_wth(1,i) + 0 ) =                           &
                 secondary_organic( map_wth(1,i) + 1 )
          end do
        case(fldname_n_nuc_sol)
          do i = 1, row_length
            n_nuc_sol( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =          &
                 real( tot_tracer( i, 1, :, n ), r_def )
            n_nuc_sol( map_wth(1,i) + 0 ) = n_nuc_sol( map_wth(1,i) + 1 )
          end do
        case(fldname_nuc_sol_su)
          do i = 1, row_length
            nuc_sol_su( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                 real( tot_tracer( i, 1, :, n ), r_def )
            nuc_sol_su( map_wth(1,i) + 0 ) = nuc_sol_su( map_wth(1,i) + 1 )
          end do
        case(fldname_nuc_sol_om)
          do i = 1, row_length
            nuc_sol_om( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                 real( tot_tracer( i, 1, :, n ), r_def )
            nuc_sol_om( map_wth(1,i) + 0 ) = nuc_sol_om( map_wth(1,i) + 1 )
          end do
        case(fldname_n_ait_sol)
          do i = 1, row_length
            n_ait_sol( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =          &
                 real( tot_tracer( i, 1, :, n ), r_def )
            n_ait_sol( map_wth(1,i) + 0 ) = n_ait_sol( map_wth(1,i) + 1 )
          end do
        case(fldname_ait_sol_su)
          do i = 1, row_length
            ait_sol_su( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                 real( tot_tracer( i, 1, :, n ), r_def )
            ait_sol_su( map_wth(1,i) + 0 ) = ait_sol_su( map_wth(1,i) + 1 )
          end do
        case(fldname_ait_sol_bc)
          do i = 1, row_length
            ait_sol_bc( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                 real( tot_tracer( i, 1, :, n ), r_def )
            ait_sol_bc( map_wth(1,i) + 0 ) = ait_sol_bc( map_wth(1,i) + 1 )
          end do
        case(fldname_ait_sol_om)
          do i = 1, row_length
            ait_sol_om( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                 real( tot_tracer( i, 1, :, n ), r_def )
            ait_sol_om( map_wth(1,i) + 0 ) = ait_sol_om( map_wth(1,i) + 1 )
          end do
        case(fldname_n_acc_sol)
          do i = 1, row_length
            n_acc_sol( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =          &
                 real( tot_tracer( i, 1, :, n ), r_def )
            n_acc_sol( map_wth(1,i) + 0 ) = n_acc_sol( map_wth(1,i) + 1 )
          end do
        case(fldname_acc_sol_su)
          do i = 1, row_length
            acc_sol_su( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                 real( tot_tracer( i, 1, :, n ), r_def )
            acc_sol_su( map_wth(1,i) + 0 ) = acc_sol_su( map_wth(1,i) + 1 )
          end do
        case(fldname_acc_sol_bc)
          do i = 1, row_length
            acc_sol_bc( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                 real( tot_tracer( i, 1, :, n ), r_def )
            acc_sol_bc( map_wth(1,i) + 0 ) = acc_sol_bc( map_wth(1,i) + 1 )
          end do
        case(fldname_acc_sol_om)
          do i = 1, row_length
            acc_sol_om( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                 real( tot_tracer( i, 1, :, n ), r_def )
            acc_sol_om( map_wth(1,i) + 0 ) = acc_sol_om( map_wth(1,i) + 1 )
          end do
        case(fldname_acc_sol_ss)
          do i = 1, row_length
            acc_sol_ss( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                 real( tot_tracer( i, 1, :, n ), r_def )
            acc_sol_ss( map_wth(1,i) + 0 ) = acc_sol_ss( map_wth(1,i) + 1 )
          end do
        case(fldname_acc_sol_du)
          ! No field to update
        case(fldname_n_cor_sol)
          do i = 1, row_length
            n_cor_sol( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =          &
                 real( tot_tracer( i, 1, :, n ), r_def )
            n_cor_sol( map_wth(1,i) + 0 ) = n_cor_sol( map_wth(1,i) + 1 )
          end do
        case(fldname_cor_sol_su)
          do i = 1, row_length
            cor_sol_su( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                 real( tot_tracer( i, 1, :, n ), r_def )
            cor_sol_su( map_wth(1,i) + 0 ) = cor_sol_su( map_wth(1,i) + 1 )
          end do
        case(fldname_cor_sol_bc)
          do i = 1, row_length
            cor_sol_bc( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                 real( tot_tracer( i, 1, :, n ), r_def )
            cor_sol_bc( map_wth(1,i) + 0 ) = cor_sol_bc( map_wth(1,i) + 1 )
          end do
        case(fldname_cor_sol_om)
          do i = 1, row_length
            cor_sol_om( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                 real( tot_tracer( i, 1, :, n ), r_def )
            cor_sol_om( map_wth(1,i) + 0 ) = cor_sol_om( map_wth(1,i) + 1 )
          end do
        case(fldname_cor_sol_ss)
          do i = 1, row_length
            cor_sol_ss( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                 real( tot_tracer( i, 1, :, n ), r_def )
            cor_sol_ss( map_wth(1,i) + 0 ) = cor_sol_ss( map_wth(1,i) + 1 )
          end do
        case(fldname_cor_sol_du)
          ! No field to update
        case(fldname_n_ait_ins)
          do i = 1, row_length
            n_ait_ins( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =          &
                 real( tot_tracer( i, 1, :, n ), r_def )
            n_ait_ins( map_wth(1,i) + 0 ) = n_ait_ins( map_wth(1,i) + 1 )
          end do
        case(fldname_ait_ins_bc)
          do i = 1, row_length
            ait_ins_bc( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                 real( tot_tracer( i, 1, :, n ), r_def )
            ait_ins_bc( map_wth(1,i) + 0 ) = ait_ins_bc( map_wth(1,i) + 1 )
          end do
        case(fldname_ait_ins_om)
          do i = 1, row_length
            ait_ins_om( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                 real( tot_tracer( i, 1, :, n ), r_def )
            ait_ins_om( map_wth(1,i) + 0 ) = ait_ins_om( map_wth(1,i) + 1 )
          end do
        case(fldname_n_acc_ins)
          do i = 1, row_length
            n_acc_ins( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =          &
                 real( tot_tracer( i, 1, :, n ), r_def )
            n_acc_ins( map_wth(1,i) + 0 ) = n_acc_ins( map_wth(1,i) + 1 )
          end do
        case(fldname_acc_ins_du)
          do i = 1, row_length
            acc_ins_du( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                 real( tot_tracer( i, 1, :, n ), r_def )
            acc_ins_du( map_wth(1,i) + 0 ) = acc_ins_du( map_wth(1,i) + 1 )
          end do
        case(fldname_n_cor_ins)
          do i = 1, row_length
            n_cor_ins( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =          &
                 real( tot_tracer( i, 1, :, n ), r_def )
            n_cor_ins( map_wth(1,i) + 0 ) = n_cor_ins( map_wth(1,i) + 1 )
          end do
        case(fldname_cor_ins_du)
          do i = 1, row_length
            cor_ins_du( map_wth(1,i) + 1 : map_wth(1,i) + nlayers ) =         &
                 real( tot_tracer( i, 1, :, n ), r_def )
            cor_ins_du( map_wth(1,i) + 0 ) = cor_ins_du( map_wth(1,i) + 1 )
          end do
        end select
      end do
      deallocate(tot_tracer)
    end if  ! outer == outer_iterations .AND. l_tracer

  end subroutine conv_comorph_code

end module conv_comorph_kernel_mod
