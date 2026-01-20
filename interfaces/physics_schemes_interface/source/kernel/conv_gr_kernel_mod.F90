!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the UM Gregory Rowntree Convection scheme.
!>
module conv_gr_kernel_mod

  use argument_mod,            only : arg_type,                  &
                                      GH_FIELD, GH_SCALAR,       &
                                      GH_INTEGER, GH_REAL,       &
                                      GH_READ, GH_WRITE,         &
                                      GH_READWRITE, DOMAIN,      &
                                      ANY_DISCONTINUOUS_SPACE_1, &
                                      ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,           only : i_def, i_um, r_def, r_um
  use tuning_segments_mod,     only : conv_gr_segment_size
  use empty_data_mod,          only : empty_real_data
  use fs_continuity_mod,       only : W3, Wtheta
  use kernel_mod,              only : kernel_type
  use mixing_config_mod,       only : smagorinsky
  use timestepping_config_mod, only : outer_iterations
  use microphysics_config_mod, only : prog_tnuc

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: conv_gr_kernel_type
    private
    type(arg_type) :: meta_args(210) = (/                                         &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                &! outer
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! rho_in_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! rho_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! wetrho_in_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! wetrho_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! exner_in_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! exner_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! w_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! theta_star
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! u_in_w3_star
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! v_in_w3_star
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! height_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! height_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! delta
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! ntml_2d
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! cumulus_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! tile_fraction
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dt_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dmv_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dmcl_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dms_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! du_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! dv_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! conv_prog_precip
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! conv_prog_dtheta
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! conv_prog_dmv
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_v
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_cl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_cf
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! cf_ice
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! cf_liq
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! cf_bulk
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cca
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! ccw
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! cape_diluted
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! conv_rain
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! conv_snow
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! conv_rain_3d
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! conv_snow_3d
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! cca_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! dd_mf_cb
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! massflux_up
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! massflux_down
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! tke_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! zh_2d
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! shallow_flag
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! uw0_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! vw0_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! lcl_height
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! parcel_top
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! level_parcel_top
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! wstar_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! thv_flux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! parcel_buoyancy
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! qsat_at_lcl
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
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! tnuc_nlcl
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! deep_in_col
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! shallow_in_col
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! mid_in_col
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! freeze_level
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! deep_prec
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! shallow_prec
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! mid_prec
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! deep_term
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! cape_timescale
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! lowest_cv_base
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! lowest_cv_top
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! cv_base
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! cv_top
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! pres_cv_base
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! pres_cv_top
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! pres_lowest_cv_base
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! pres_lowest_cv_top
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! lowest_cca_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! deep_cfl_limited
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),&! mid_cfl_limited
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! entrain_up
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! entrain_down
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! detrain_up
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! detrain_down
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dd_dt
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! dd_dq
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! deep_dt
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! deep_dq
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! deep_massflux
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! shallow_dt
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! shallow_dq
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! shallow_massflux
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! mid_dt
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! mid_dq
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! mid_massflux
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! deep_tops
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! massflux_up_half
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! massflux_up_cmpta
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTHETA),                   &! cca_unadjusted
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     WTHETA),                   &! dth_conv_noshal
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     WTHETA)                    &! dmv_conv_noshal
        /)
    integer :: operates_on = DOMAIN
  contains
    procedure, nopass :: conv_gr_code
  end type

  public :: conv_gr_code

contains

  !> @brief Interface to the UM Gregory Rowntree convection scheme
  !> @details The UM Gregory Rowntree convection scheme does:
  !>             vertical mixing of heat, momentum and moisture,
  !>             as documented in UMDP27
  !> @param[in]     nlayers              Number of layers
  !> @param[in]     ncells               Number of cells
  !> @param[in]     outer                Outer loop counter
  !> @param[in]     rho_in_w3            Density field in density space
  !> @param[in]     rho_in_wth           Density field in wth space
  !> @param[in]     wetrho_in_w3         Wet density field in density space
  !> @param[in]     wetrho_in_wth        Wet density field in wth space
  !> @param[in]     exner_in_w3          Exner pressure field in density space
  !> @param[in]     exner_in_wth         Exner pressure field in wth space
  !> @param[in]     w_in_wth             'Vertical' wind in theta space
  !> @param[in]     theta_star           Potential temperature after advection
  !> @param[in]     u_in_w3_star         'Zonal' wind after advection
  !> @param[in]     v_in_w3_star         'Meridional' wind after advection
  !> @param[in]     height_w3            Height of density space above surface
  !> @param[in]     height_wth           Height of theta space above surface
  !> @param[in]     delta                Edge length on wtheta points
  !> @param[in]     ntml_2d              Number of turbulently mixed levels
  !> @param[in]     cumulus_2d           Cumulus flag (true/false)
  !> @param[in]     tile_fraction        Surface tile fractions
  !> @param[in,out] dt_conv              Convection temperature increment
  !> @param[in,out] dmv_conv             Convection vapour increment
  !> @param[in,out] dmcl_conv            Convection cloud liquid increment
  !> @param[in,out] dms_conv             Convection cloud snow increment
  !> @param[in,out] du_conv              Convection 'zonal' wind increment
  !> @param[in,out] dv_conv              Convection 'meridional' wind increment
  !> @param[in,out] conv_prog_precip     Surface precipitation based 3d convective prognostic (kg/m2/s)
  !> @param[in,out] conv_prog_dtheta     Time smoothed convective theta increment
  !> @param[in,out] conv_prog_dmv        Time smoothed convective humidity increment
  !> @param[in]     m_v                  Vapour mixing ratio after advection
  !> @param[in]     m_cl                 Cloud liq mixing ratio after advection
  !> @param[in]     m_cf                 Cloud fro mixing ratio after advection
  !> @param[in]     cf_ice               Ice cloud fraction
  !> @param[in]     cf_liq               Liquid cloud fraction
  !> @param[in]     cf_bulk              Bulk cloud fraction
  !> @param[in,out] cca                  Convective cloud amount (fraction)
  !> @param[in,out] ccw                  Convective cloud water (kg/kg) (can be ice or liquid)
  !> @param[in,out] cape_diluted         CAPE value
  !> @param[in,out] conv_rain            Surface rainfall rate from convection (kg/m2/s)
  !> @param[in,out] conv_snow            Surface snowfall rate from convection (kg/m2/s)
  !> @param[in,out] conv_rain_3d         Rainfall rate from convection (kg/m2/s)
  !> @param[in,out] conv_snow_3d         Snowfall rate from convection (kg/m2/s)
  !> @param[in,out] cca_2d               Convective cloud amout (2d) with no anvil
  !> @param[in,out] dd_mf_cb             Downdraft massflux at cloud base (Pa/s)
  !> @param[in,out] massflux_up          Convective upwards mass flux (Pa/s)
  !> @param[in,out] massflux_down        Convective downwards mass flux (Pa/s)
  !> @param[in,out] tke_bl               Turbulent kinetic energy (m2/s2)
  !> @param[in]     zh_2d                Boundary layer depth
  !> @param[in]     shallow_flag         Indicator of shallow convection
  !> @param[in]     uw0_flux             'Zonal' surface momentum flux
  !> @param[in]     vw0 flux             'Meridional' surface momentum flux
  !> @param[in]     lcl_height           Height of lifting condensation level
  !> @param[in]     parcel_top           Height of surface based parcel ascent
  !> @param[in]     level_parcel_top     Model level of parcel_top
  !> @param[in]     wstar_2d             BL velocity scale
  !> @param[in]     thv_flux             Surface flux of theta_v
  !> @param[in]     parcel_buoyancy      Integral of parcel buoyancy
  !> @param[in]     qsat_at_lcl          Saturation specific hum at LCL
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
  !> @param[in]     tnuc_nlcl            Temperature of nucleation (K) (2D)
  !> @param[in,out] deep_in_col          Indicator of deep in column
  !> @param[in,out] shallow_in_col       Indicator of shallow in column
  !> @param[in,out] mid_in_col           Indicator of mid in column
  !> @param[in,out] freeze_level         Level number of freezing level
  !> @param[in,out] deep_prec            Precipitation rate from deep convection(kg/m2/s)
  !> @param[in,out] shallow_prec         Precipitation rate from shallow convection(kg/m2/s)
  !> @param[in,out] mid_prec             Precipitation rate from mid convection(kg/m2/s)
  !> @param[in,out] deep_term            Termination level number of deep convection
  !> @param[in,out] cape_timescale       CAPE timescale (s)
  !> @param[in,out] lowest_cv_base       Level number for start of convection in column
  !> @param[in,out] lowest_cv_top        Level number for end of lowest convection in column
  !> @param[in,out] cv_base              Level number of base of highest convection in column
  !> @param[in,out] cv_top               Level number for end of highest convection in column
  !> @param[in,out] pres_cv_base         Pressure at base of highest convection in column
  !> @param[in,out] pres_cv_top          Pressure at end of highest convection in column
  !> @param[in,out] pres_lowest_cv_base  Pressure at base of lowest convection in column
  !> @param[in,out] pres_lowest_cv_top   Pressure at end of lowest convection in column
  !> @param[in,out] lowest_cca_2d        2D convective cloud amount of lowest convecting layer
  !> @param[in,out] deep_cfl_limited     Deep convection CFL limited
  !> @param[in,out] mid_cfl_limited      Mid convection CFL limited
  !> @param[in,out] entrain_up           Convective upwards entrainment
  !> @param[in,out] entrain_down         Convective downwards entrainment
  !> @param[in,out] detrain_up           Convective upwards detrainment
  !> @param[in,out] detrain_down         Convective downwards detrainment
  !> @param[in,out] dd_dt                Temperature increment from downdraughts per timestep
  !> @param[in,out] dd_dq                Vapour increment from downdraughts per timestep
  !> @param[in,out] deep_dt              Temperature increment from deep convection per timestep
  !> @param[in,out] deep_dq              Vapour increment from deep convection per timestep
  !> @param[in,out] deep_massflux        Upward mass flux from deep convection
  !> @param[in,out] shallow_dt           Temperature increment from shallow convection per timestep
  !> @param[in,out] shallow_dq           Vapour increment from shallow convection per timestep
  !> @param[in,out] shallow_massflux     Upward mass flux from shallow convection
  !> @param[in,out] mid_dt               Temperature increment from mid convection per timestep
  !> @param[in,out] mid_dq               Vapour increment from mid convection per timestep
  !> @param[in,out] mid_massflux         Upward mass flux from mid convection
  !> @param[in,out] deep_tops            Set to 1.0 if deep stops at this model level
  !> @param[in,out] massflux_up_half     Convective upwards mass flux on half-levels (Pa/s)
  !> @param[in,out] massflux_up_cmpta    Convective upwards mass flux component A (Pa/s)
  !> @param[in,out] cca_unadjusted       Convective cloud amout unadjusted for radiation calc.
  !> @param[out]    dth_conv_noshal      Convection theta increment from non-shallow regions
  !> @param[out]    dmv_conv_noshal      Convection vapour increment from non-shallow regions
  !> @param[in]     ndf_w3               Number of DOFs per cell for density space
  !> @param[in]     undf_w3              Number of unique DOFs  for density space
  !> @param[in]     map_w3               Dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_wth              Number of DOFs per cell for potential temperature space
  !> @param[in]     undf_wth             Number of unique DOFs for potential temperature space
  !> @param[in]     map_wth              Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_2d               Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d              Number of unique DOFs  for 2D fields
  !> @param[in]     map_2d               Dofmap for the cell at the base of the column for 2D fields
  !> @param[in]     ndf_tile             Number of DOFs per cell for tiles
  !> @param[in]     undf_tile            Number of total DOFs for tiles
  !> @param[in]     map_tile             Dofmap for cell for surface tiles
  subroutine conv_gr_code(nlayers,                           &
                          ncells,                            &
                          outer,                             &
                          rho_in_w3,                         &
                          rho_in_wth,                        &
                          wetrho_in_w3,                      &
                          wetrho_in_wth,                     &
                          exner_in_w3,                       &
                          exner_in_wth,                      &
                          w_in_wth,                          &
                          theta_star,                        &
                          u_in_w3_star,                      &
                          v_in_w3_star,                      &
                          height_w3,                         &
                          height_wth,                        &
                          delta,                             &
                          ntml_2d,                           &
                          cumulus_2d,                        &
                          tile_fraction,                     &
                          dt_conv,                           &
                          dmv_conv,                          &
                          dmcl_conv,                         &
                          dms_conv,                          &
                          du_conv,                           &
                          dv_conv,                           &
                          conv_prog_precip,                  &
                          conv_prog_dtheta,                  &
                          conv_prog_dmv,                     &
                          m_v,                               &
                          m_cl,                              &
                          m_cf,                              &
                          cf_ice,                            &
                          cf_liq,                            &
                          cf_bulk,                           &
                          cca,                               &
                          ccw,                               &
                          cape_diluted,                      &
                          conv_rain,                         &
                          conv_snow,                         &
                          conv_rain_3d,                      &
                          conv_snow_3d,                      &
                          cca_2d,                            &
                          dd_mf_cb,                          &
                          massflux_up,                       &
                          massflux_down,                     &
                          tke_bl,                            &
                          zh_2d,                             &
                          shallow_flag,                      &
                          uw0_flux,                          &
                          vw0_flux,                          &
                          lcl_height,                        &
                          parcel_top,                        &
                          level_parcel_top,                  &
                          wstar_2d,                          &
                          thv_flux,                          &
                          parcel_buoyancy,                   &
                          qsat_at_lcl,                       &
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
                          tnuc_nlcl,                         &
                          deep_in_col,                       &
                          shallow_in_col,                    &
                          mid_in_col,                        &
                          freeze_level,                      &
                          deep_prec,                         &
                          shallow_prec,                      &
                          mid_prec,                          &
                          deep_term,                         &
                          cape_timescale,                    &
                          lowest_cv_base,                    &
                          lowest_cv_top,                     &
                          cv_base,                           &
                          cv_top,                            &
                          pres_cv_base,                      &
                          pres_cv_top,                       &
                          pres_lowest_cv_base,               &
                          pres_lowest_cv_top,                &
                          lowest_cca_2d,                     &
                          deep_cfl_limited,                  &
                          mid_cfl_limited,                   &
                          entrain_up,                        &
                          entrain_down,                      &
                          detrain_up,                        &
                          detrain_down,                      &
                          dd_dt,                             &
                          dd_dq,                             &
                          deep_dt,                           &
                          deep_dq,                           &
                          deep_massflux,                     &
                          shallow_dt,                        &
                          shallow_dq,                        &
                          shallow_massflux,                  &
                          mid_dt,                            &
                          mid_dq,                            &
                          mid_massflux,                      &
                          deep_tops,                         &
                          massflux_up_half,                  &
                          massflux_up_cmpta,                 &
                          cca_unadjusted,                    &
                          dth_conv_noshal,                   &
                          dmv_conv_noshal,                   &
                          ndf_w3,                            &
                          undf_w3,                           &
                          map_w3,                            &
                          ndf_wth,                           &
                          undf_wth,                          &
                          map_wth,                           &
                          ndf_2d,                            &
                          undf_2d,                           &
                          map_2d,                            &
                          ndf_tile,                          &
                          undf_tile,                         &
                          map_tile )

    !---------------------------------------
    ! LFRic modules
    !---------------------------------------

    use jules_control_init_mod, only: n_land_tile

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

    !---------------------------------------
    ! UM modules containing switches or global constants
    !---------------------------------------
    use bl_option_mod, only: max_tke
    use cloud_inputs_mod, only: i_cld_vn
    use cv_run_mod, only: n_conv_calls, iconv_deep, iconv_shallow, l_mom, &
                          qmin_conv, l_safe_conv,                         &
                          l_conv_prog_dtheta, l_conv_prog_dq,             &
                          l_conv_prog_precip,                             &
                          tau_conv_prog_dtheta, tau_conv_prog_dq,         &
                          tau_conv_prog_precip
    use cv_param_mod, only: max_mf_fall, dthetadt_conv_active_threshold, &
                            conv_prog_precip_min_threshold
    use jules_surface_mod, only: srf_ex_cnv_gust, IP_SrfExWithCnv
    use nlsizes_namelist_mod, only: bl_levels, n_cca_lev
    use pc2_constants_mod, only: i_cld_pc2
    use planet_constants_mod, only: p_zero, kappa, planet_radius, g
    use scm_convss_dg_mod, only: scm_convss_dg_type
    use timestep_mod, only: timestep

    use free_tracers_inputs_mod, only: n_wtrac
    use wtrac_conv_mod, only: l_wtrac_conv

    ! subroutines used
    use glue_conv_6a_mod, only: glue_conv_6a
    use ukca_api_mod, only: ukca_get_tracer_varlist, ukca_maxlen_fieldname

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, ncells
    integer(kind=i_def), intent(in) :: outer

    integer(kind=i_def), intent(in) :: ndf_wth, ndf_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3
    integer(kind=i_def), intent(in) :: map_wth(ndf_wth, ncells)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3, ncells)
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d, ncells)

    integer(kind=i_def), intent(in) :: ndf_tile, undf_tile
    integer(kind=i_def), intent(in) :: map_tile(ndf_tile, ncells)

    real(kind=r_def), dimension(undf_w3), intent(in) :: rho_in_w3,          &
                                                        wetrho_in_w3,       &
                                                        exner_in_w3,        &
                                                        u_in_w3_star,       &
                                                        v_in_w3_star,       &
                                                        height_w3
    real(kind=r_def), dimension(undf_wth), intent(in) :: cf_ice,            &
                                                         cf_liq, cf_bulk,   &
                                                         m_v, m_cl, m_cf,   &
                                                         rho_in_wth,        &
                                                         wetrho_in_wth,     &
                                                         exner_in_wth,      &
                                                         w_in_wth,          &
                                                         theta_star,        &
                                                         height_wth,        &
                                                         delta

    real(kind=r_def), dimension(undf_wth), intent(inout) :: dt_conv, dmv_conv, &
                                          dmcl_conv, dms_conv, cca, ccw,       &
                                          massflux_up, massflux_down,          &
                                          conv_rain_3d, conv_snow_3d, tke_bl,  &
                                          conv_prog_precip,                    &
                                          conv_prog_dtheta, conv_prog_dmv

    real(kind=r_def), dimension(undf_w3), intent(inout) :: du_conv, dv_conv

    integer(kind=i_def), dimension(undf_2d), intent(in) :: ntml_2d,            &
                                                           cumulus_2d,         &
                                                           shallow_flag,       &
                                                           level_parcel_top

    real(kind=r_def), dimension(undf_2d), intent(in) :: zh_2d,                &
                                                        uw0_flux,             &
                                                        vw0_flux,             &
                                                        lcl_height,           &
                                                        parcel_top,           &
                                                        wstar_2d,             &
                                                        thv_flux,             &
                                                        parcel_buoyancy,      &
                                                        qsat_at_lcl,          &
                                                        tnuc_nlcl

    real(kind=r_def), dimension(undf_2d), intent(inout) :: cape_diluted,  &
                                   conv_rain, conv_snow, cca_2d, dd_mf_cb
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

    real(kind=r_def), pointer, intent(inout) :: deep_in_col(:),            &
                                                shallow_in_col(:),         &
                                                mid_in_col(:),             &
                                                freeze_level(:),           &
                                                deep_prec(:),              &
                                                shallow_prec(:),           &
                                                mid_prec(:),               &
                                                deep_term(:),              &
                                                cape_timescale(:),         &
                                                lowest_cv_base(:),         &
                                                lowest_cv_top(:),          &
                                                cv_base(:),                &
                                                cv_top(:),                 &
                                                pres_cv_base(:),           &
                                                pres_cv_top(:),            &
                                                pres_lowest_cv_base(:),    &
                                                pres_lowest_cv_top(:),     &
                                                lowest_cca_2d(:),          &
                                                deep_cfl_limited(:),       &
                                                mid_cfl_limited(:)

    real(kind=r_def), pointer, intent(inout) :: entrain_up(:),       &
                                                entrain_down(:),     &
                                                detrain_up(:),       &
                                                detrain_down(:),     &
                                                dd_dt(:),            &
                                                dd_dq(:),            &
                                                deep_dt(:),          &
                                                deep_dq(:),          &
                                                deep_massflux(:),    &
                                                shallow_dt(:),       &
                                                shallow_dq(:),       &
                                                shallow_massflux(:), &
                                                mid_dt(:),           &
                                                mid_dq(:),           &
                                                mid_massflux(:),     &
                                                deep_tops(:),        &
                                                massflux_up_half(:), &
                                                massflux_up_cmpta(:),&
                                                cca_unadjusted(:),   &
                                                dth_conv_noshal(:),  &
                                                dmv_conv_noshal(:)

    real(kind=r_def), dimension(undf_wth), intent(inout) :: dcfl_conv
    real(kind=r_def), dimension(undf_wth), intent(inout) :: dcff_conv
    real(kind=r_def), dimension(undf_wth), intent(inout) :: dbcf_conv

    real(kind=r_def), intent(in) :: tile_fraction(undf_tile)

    !-----------------------------------------------------------------------
    ! Local variables for the kernel
    !-----------------------------------------------------------------------
    ! loop counters etc
    integer(i_def) :: k, i, n, ii, j

    ! local switches and scalars
    integer(i_um) :: n_deep, n_shallow, n_congestus, n_mid,                  &
                     ntra_lev, nml_segment_size, seg_size, n_conv_levels,    &
                     call_number, ntra_fld, seg_num, num_seg

    logical :: l_tracer, l_calc_dxek, l_q_interact, l_scm_convss_dg

    real(r_um) :: timestep_conv, one_over_conv_calls, orig_value,            &
                  decay_amount, conv_active, theta_inc_threshold

    ! profile fields from level 1 upwards
    real(r_um), dimension(ncells,1,nlayers) ::                               &
         p_rho_levels, rho_wet, rho_dry, z_rho, z_theta, cca_3d, rho_wet_tq, &
         rho_dry_theta, exner_rho_levels, r_rho_levels,                      &
         theta_conv, q_conv, qcl_conv, qcf_conv,                             &
         dtheta_conv,                                                        &
         cf_liquid_conv, cf_frozen_conv,                                     &
         bulk_cf_conv, u_conv, v_conv, dq_add, ccw_3d, dthbydt, dqbydt,      &
         dqclbydt, dqcfbydt, dcflbydt, dcffbydt, dbcfbydt, dubydt_p,         &
         dvbydt_p, it_ccw, it_ccw0, it_conv_rain_3d, it_conv_snow_3d, it_cca,&
         it_cca0, it_w2p, it_cca0_dp, it_cca0_md, it_cca0_sh, it_up_flux,    &
         it_dwn_flux, it_entrain_up, it_detrain_up, it_entrain_dwn,          &
         it_detrain_dwn, it_mf_deep, it_mf_shall, it_mf_midlev,              &
         it_up_flux_half,                                                    &
         it_dt_deep, it_dt_shall, it_dt_midlev,                              &
         it_dq_deep, it_dq_shall, it_dq_midlev,                              &
         it_du_deep, it_du_shall, it_du_midlev,                              &
         it_dv_deep, it_dv_shall, it_dv_midlev,                              &
         it_dt_dd, it_dq_dd, tnuc_new

    ! profile fields from level 0 upwards
    real(r_um), dimension(ncells,1,0:nlayers) ::                             &
         p_theta_levels, p_rho_minus_one, w, r_theta_levels,                 &
         exner_rho_minus_one, exner_theta_levels

    ! single level real fields
    real(r_um), dimension(ncells,1) ::                                       &
         p_star, zhpar, zh, wstar, wthvs, zlcl_uv, entrain_coef,             &
         qsat_lcl, delthvu, flandg, uw0, vw0, it_lcca, it_cca_2d, it_cclwp,  &
         it_cclwp0, it_conv_rain, it_conv_snow, it_precip_dp, it_precip_sh,  &
         it_precip_md, it_cape_diluted, it_dp_cfl_limited,                   &
         it_md_cfl_limited, cape_ts_used, it_ind_deep, it_ind_shall,         &
         it_precip_cg, it_wstar_up, it_mb1, it_mb2, tot_conv_precip_2d,      &
         delta_smag, tnuc_nlcl_um

    ! single level integer fields
    integer(i_um), dimension(ncells,1) :: ntml, ntpar, lcbase,               &
         it_lcbase,it_lctop, it_ccb, it_cct, it_ccb0, it_cct0, it_kterm_deep,&
         it_kterm_shall, it_cg_term, it_lcbase0, freeze_lev, ccb, cct, lctop

    ! single level logical fields
    logical, dimension(ncells,1) :: land_sea_mask, cumulus,                  &
                                           l_shallow, l_congestus,           &
                                           it_mid_level, l_mid

    ! total tracer - has to be allocatable
    real(r_um), dimension(:,:,:,:), allocatable :: tot_tracer

    ! Variables for retrieving tracer names for a UKCA configuration
    integer :: ukca_errcode
    character(len=ukca_maxlen_fieldname), pointer :: ukca_tracer_names(:)
    character(len=ukca_maxlen_fieldname), target  ::                         &
                 local_dust_tracer_list(4) = [ 'Acc_INS_N ' , 'Acc_INS_DU' , &
                                               'Cor_INS_N ' , 'Cor_INS_DU' ]

    ! Fields which are not used and only required for subroutine argument list,
    ! hence are unset in the kernel
    ! if they become set, please move up to be with other variables
    type(scm_convss_dg_type), allocatable :: scm_convss_dg(:)

    real(r_um), dimension(ncells,1,nlayers) ::                               &
         it_mf_congest, it_dt_congest, it_dq_congest, it_du_congest,         &
         it_dv_congest, it_du_dd, it_dv_dd, it_area_ud, it_area_dd,          &
         it_uw_dp, it_vw_dp, it_uw_shall, it_vw_shall,                       &
         it_uw_mid, it_vw_mid, it_wqt_flux, it_wthetal_flux, it_wthetav_flux,&
         it_wql_flux, conv_prog_flx

    real(r_um), dimension(ncells,1) :: fqw, ftl

    real(r_um), dimension(ncells,1,nlayers) :: conv_prog_precip_conv

    real(r_um), dimension(ncells,1) :: zlcl, t1_sd, q1_sd, w_max,            &
         deep_flag, past_conv_ht, ql_ad, ind_cape_reduced,                   &
         it_wstar_dn, g_ccp, h_ccp, ccp_strength

    integer(i_um), dimension(ncells,1) :: conv_type

    ! Water tracer fields which are not currently used but are required by
    ! UM routine
    real(r_um), dimension(ncells,nlayers,n_wtrac) :: q_wtrac, qcl_wtrac,     &
          qcf_wtrac, dqbydt_wtrac, dqclbydt_wtrac, dqcfbydt_wtrac
    real(r_um), dimension(ncells,n_wtrac) :: rain_wtrac, snow_wtrac

    ! Current assumptions about setup based on GA9

    ! If this is the last solver outer loop then tracers may need convecting.
    ! Enable tracers for UKCA if a UKCA tracer list is available
    ! (This indicates that a UKCA configuration has been set up)
    if (outer == outer_iterations) then
      if ( glomap_mode == glomap_mode_dust_and_clim ) then
        ukca_tracer_names => local_dust_tracer_list
        l_tracer = .true.
      else
        call ukca_get_tracer_varlist( ukca_tracer_names, ukca_errcode )
        l_tracer = ( ukca_errcode == 0 )
      end if
    else
      l_tracer = .false.
    end if

    if ( outer == outer_iterations .AND. l_tracer ) then

      ntra_fld = size(ukca_tracer_names)
      ntra_lev = nlayers
      allocate(tot_tracer( ncells, 1, ntra_lev, ntra_fld ))

    else

      ! No tracers: set up a dummy tracer array
      ntra_fld = 1
      ntra_lev = 1
      allocate(tot_tracer( ncells, 1, ntra_lev, ntra_fld ))

    end if  ! outer == outer_iterations .AND. l_tracer

    !-----------------------------------------------------------------------
    ! Mapping of LFRic fields into UM variables
    !-----------------------------------------------------------------------

    ! Land sea mask
    do i = 1, ncells
      flandg(i,1) = 0.0_r_um
      do n = 1, n_land_tile
        flandg(i,1) = flandg(i,1) + real(tile_fraction(map_tile(1,i)+n-1), r_um)
      end do
    end do
    do i = 1, ncells
      ! Jules requires fractions with respect to the land area
      if (flandg(i,1) > 0.0_r_um) then
        land_sea_mask(i,1) = .true.
      else
        land_sea_mask(i,1) = .false.
      end if
    end do

    !-----------------------------------------------------------------------
    ! Build vectors of columns in the horizontal by copying data to
    ! arrays sized (i,1,k) to match the UMs (i,j,k) data layout.
    ! Assuming map_wth(1,i) points to level 0
    ! and map_w3(1,i) points to level 1
    !-----------------------------------------------------------------------
    do k = 0, nlayers
      do i = 1, ncells
        ! pressure on theta levels
        p_theta_levels(i,1,k) = p_zero*(exner_in_wth(map_wth(1,i) + k))**(1.0_r_def/kappa)
        ! exner pressure on theta levels
        exner_theta_levels(i,1,k) = exner_in_wth(map_wth(1,i) + k)
        ! height of theta levels from centre of planet
        r_theta_levels(i,1,k) = height_wth(map_wth(1,i) + k) + planet_radius
        ! w wind on theta levels
        w(i,1,k) = w_in_wth(map_wth(1,i) + k)
      end do
    end do
    do k = 1, nlayers
      do i = 1, ncells
        ! wet density on theta and rho levels
        rho_wet_tq(i,1,k) = wetrho_in_wth(map_wth(1,i) + k)
        rho_wet(i,1,k) = wetrho_in_w3(map_w3(1,i) + k-1)
        ! dry density on theta and rho levels
        rho_dry_theta(i,1,k) = rho_in_wth(map_wth(1,i) + k)
        rho_dry(i,1,k) = rho_in_w3(map_w3(1,i) + k-1)
        ! pressure on rho levels
        p_rho_levels(i,1,k) = p_zero*(exner_in_w3(map_w3(1,i) + k-1))**(1.0_r_def/kappa)
        ! exner pressure on rho levels
        exner_rho_levels(i,1,k) = exner_in_w3(map_w3(1,i) + k-1)
        ! height of rho levels from centre of planet
        r_rho_levels(i,1,k) = height_w3(map_w3(1,i) + k-1) + planet_radius
      end do
    end do

    if ( smagorinsky ) then
      do i = 1, ncells
        delta_smag(i,1) = delta(map_wth(1,i))
      end do
    end if

    do i = 1, ncells
      ! surface pressure
      p_star(i,1) = p_theta_levels(i,1,0)
      ! setup odd array which is on rho levels but without level 1
      p_rho_minus_one(i,1,0) = p_theta_levels(i,1,0)
      do k = 1, nlayers-1
        p_rho_minus_one(i,1,k) = p_rho_levels(i,1,k+1)
      end do
      p_rho_minus_one(i,1,nlayers) = 0.0_r_um
      ! and similar array for exner
      exner_rho_minus_one(i,1,0)   = exner_theta_levels(i,1,0)
      do k = 1, nlayers-1
        exner_rho_minus_one(i,1,k) = exner_rho_levels(i,1,k+1)
      end do
      exner_rho_minus_one(i,1,nlayers) = 0.0_r_um
      ! height of levels above surface
      do k = 1, nlayers
        z_rho(i,1,k) = r_rho_levels(i,1,k)-r_theta_levels(i,1,0)
        z_theta(i,1,k) = r_theta_levels(i,1,k)-r_theta_levels(i,1,0)
      end do
    end do

    !-----------------------------------------------------------------------
    ! Things passed from other parametrization schemes on this timestep
    !----------------------------------------------------------------------
    do i = 1, ncells
      cumulus(i,1) = (cumulus_2d(map_2d(1,i)) == 1_i_def)
      ntml(i,1) = ntml_2d(map_2d(1,i))
  
      zh(i,1) = zh_2d(map_2d(1,i))
      l_shallow(i,1) = (shallow_flag(map_2d(1,i)) == 1_i_def)
      uw0(i,1) = uw0_flux(map_2d(1,i))
      vw0(i,1) = vw0_flux(map_2d(1,i))
      zlcl_uv(i,1) = lcl_height(map_2d(1,i))
      zhpar(i,1) = parcel_top(map_2d(1,i))
      ntpar(i,1) = level_parcel_top(map_2d(1,i))
      wstar(i,1) = wstar_2d(map_2d(1,i))
      wthvs(i,1) = thv_flux(map_2d(1,i))
      delthvu(i,1) = parcel_buoyancy(map_2d(1,i))
      qsat_lcl(i,1) = qsat_at_lcl(map_2d(1,i))
    end do

    !========================================================================
    ! Call to 6A Gregory-Rowntree convection scheme
    !========================================================================


    l_calc_dxek  = ( i_cld_vn == i_cld_pc2 )
    l_q_interact = l_calc_dxek
    l_congestus = .false. ! never used in GA9
    entrain_coef = -99.0_r_um  ! unused default value

    n_conv_levels = nlayers
    if (l_mom) then
      ! Limit convection calling levels to maximum of model_levels - 1
      ! This is because CMT increments to u and v exist on n_levels + 1
      if (n_conv_levels  >   nlayers - 1 ) then
        n_conv_levels = nlayers - 1
      end if
    end if

    timestep_conv=timestep/real(n_conv_calls)

    ! Sub-timestep scheme
    one_over_conv_calls = 1.0_r_um/real(n_conv_calls)

    do i = 1, ncells
      l_mid(i,1) = .true.
    end do

    do k=1,nlayers
      do i = 1, ncells
        ! Pointing to _star values
        theta_conv(i,1,k) = theta_star(map_wth(1,i) + k)
        q_conv(i,1,k)   = m_v(map_wth(1,i) + k)
        qcl_conv(i,1,k) = m_cl(map_wth(1,i) + k)
        qcf_conv(i,1,k) = m_cf(map_wth(1,i) + k)
  
        cf_liquid_conv(i,1,k) = cf_liq(map_wth(1,i) + k)
        cf_frozen_conv(i,1,k) = cf_ice(map_wth(1,i) + k)
        bulk_cf_conv(i,1,k)   = cf_bulk(map_wth(1,i) + k)
        ! Total theta increment
        dtheta_conv(i,1,k) = 0.0_r_um
        cca_3d(i,1,k) = 0.0_r_um
        ccw_3d(i,1,k) = 0.0_r_um
      end do
    end do
    if (l_conv_prog_precip) then
      do k = 1, nlayers
        do i = 1, ncells
          conv_prog_precip_conv(i,1,k) = conv_prog_precip(map_wth(1,i) + k)
        end do
      end do
    end if

    do i = 1, ncells
      ccb(i,1) = 0_i_um
      cct(i,1) = 0_i_um
      lctop(i,1) = 0_i_um
      lcbase(i,1) = 0_i_um
      w_max(i,1) = 0.0_r_um
    end do

    ! Turn water tracers off in convection and initialise dummy fields
    l_wtrac_conv = .FALSE.

    do k = 1,n_wtrac
      do j = 1,nlayers
        do i = 1,ncells
          q_wtrac(i,j,k) = 0.0_r_um
          qcl_wtrac(i,j,k) = 0.0_r_um
          qcf_wtrac(i,j,k) = 0.0_r_um
          dqbydt_wtrac(i,j,k) = 0.0_r_um
          dqclbydt_wtrac(i,j,k) = 0.0_r_um
          dqcfbydt_wtrac(i,j,k) = 0.0_r_um
        end do
      end do
    end do
    do j = 1,n_wtrac
      do i = 1,ncells
        rain_wtrac(i,j) = 0.0_r_um
        snow_wtrac(i,j) = 0.0_r_um
      end do
    end do

    ! Check for negative (less than a minimum) q being passed to convection
    if (l_safe_conv) then
      do k=1,nlayers
        do i = 1, ncells
          if (q_conv(i,1,k) < qmin_conv) then
            ! Should really be warning print statements if this is happening
            ! but log_event calls not allowed.
            ! store q added to ensure sensible profile
            dq_add(i,1,k) = qmin_conv - q_conv(i,1,k)
            ! Reset to qmin in non-conservative way
            q_conv(i,1,k) = qmin_conv
          else
            dq_add(i,1,k) = 0.0_r_um
          end if
        end do ! i
      end do ! k
    end if
    if (l_mom) then
      do k=1,nlayers
        do i = 1, ncells
          u_conv(i,1,k) = u_in_w3_star(map_w3(1,i) + k-1)
          v_conv(i,1,k) = v_in_w3_star(map_w3(1,i) + k-1)
        end do ! i
      end do ! k
    end if

    ! Map tracer fields to UM tracer array

    if ( outer == outer_iterations .AND. l_tracer ) then
      do n = 1, ntra_fld
        select case(ukca_tracer_names(n))
        case(fldname_o3p)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( o3p( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_o1d)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( o1d( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_o3)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( o3( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_n)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( nit( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_no)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( no( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_no3)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( no3( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_lumped_n)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( lumped_n( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_n2o5)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( n2o5( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_ho2no2)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( ho2no2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_hono2)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( hono2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_h2o2)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( h2o2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_ch4)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( ch4( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_co)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( co( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_hcho)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( hcho( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
            end do
        case(fldname_meoo)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( meoo( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_meooh)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( meooh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_h)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( h( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_ch2o)
          do i = 1, ncells
            ! H2O tracer from chemistry ns not transported
            tot_tracer( i, 1, :, n ) = 0.0_r_um
          end do
        case(fldname_oh)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( oh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_ho2)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( ho2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_cl)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( cl( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_cl2o2)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( cl2o2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_clo)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( clo( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_oclo)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( oclo( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_br)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( br( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_lumped_br)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( lumped_br( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_brcl)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( brcl( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_brono2)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( brono2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_n2o)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( n2o( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_lumped_cl)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( lumped_cl( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_hocl)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( hocl( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_hbr)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( hbr( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_hobr)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( hobr( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_clono2)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( clono2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_cfcl3)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( cfcl3( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_cf2cl2)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( cf2cl2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_mebr)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( mebr( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_hono)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( hono( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_c2h6)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( c2h6( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_etoo)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( etoo( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_etooh)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( etooh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_mecho)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( mecho( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_meco3)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( meco3( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_pan)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( pan( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_c3h8)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( c3h8( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_n_proo)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( n_proo( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_i_proo)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( i_proo( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_n_prooh)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( n_prooh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_i_prooh)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( i_prooh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_etcho)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( etcho( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_etco3)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( etco3( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_me2co)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( me2co( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_mecoch2oo)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( mecoch2oo( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_mecoch2ooh)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( mecoch2ooh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_ppan)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( ppan( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_meono2)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( meono2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_c5h8)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( c5h8( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_iso2)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( iso2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_isooh)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( isooh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_ison)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( ison( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_macr)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( macr( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_macro2)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( macro2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_macrooh)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( macrooh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_mpan)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( mpan( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_hacet)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( hacet( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_mgly)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( mgly( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_nald)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( nald( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_hcooh)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( hcooh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_meco3h)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( meco3h( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_meco2h)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( meco2h( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_h2)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( h2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_meoh)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( meoh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_msa)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( msa( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_nh3)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( nh3( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_cs2)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( cs2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_csul)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( csul( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_h2s)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( h2s( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_so3)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( so3( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_passive_o3)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( passive_o3( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_age_of_air)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( age_of_air( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_dms)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( dms( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_so2)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( so2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_h2so4)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( h2so4( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_dmso)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( dmso( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_monoterpene)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( monoterpene( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_secondary_organic)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( secondary_organic( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), &
                    r_um )
          end do
        case(fldname_n_nuc_sol)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( n_nuc_sol( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_nuc_sol_su)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( nuc_sol_su( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_nuc_sol_om)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( nuc_sol_om( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_n_ait_sol)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( n_ait_sol( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_ait_sol_su)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( ait_sol_su( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_ait_sol_bc)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( ait_sol_bc( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_ait_sol_om)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( ait_sol_om( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_n_acc_sol)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( n_acc_sol( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_acc_sol_su)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( acc_sol_su( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_acc_sol_bc)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( acc_sol_bc( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_acc_sol_om)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( acc_sol_om( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_acc_sol_ss)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( acc_sol_ss( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_acc_sol_du)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) = 0.0_r_um ! no prognostic, always zero
          end do
        case(fldname_n_cor_sol)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( n_cor_sol( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_cor_sol_su)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( cor_sol_su( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_cor_sol_bc)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( cor_sol_bc( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_cor_sol_om)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( cor_sol_om( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_cor_sol_ss)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( cor_sol_ss( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_cor_sol_du)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) = 0.0_r_um ! no prognostic, always zero
          end do
        case(fldname_n_ait_ins)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( n_ait_ins( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_ait_ins_bc)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( ait_ins_bc( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_ait_ins_om)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( ait_ins_om( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_n_acc_ins)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( n_acc_ins( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_acc_ins_du)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( acc_ins_du( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_n_cor_ins)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( n_cor_ins( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case(fldname_cor_ins_du)
          do i = 1, ncells
            tot_tracer( i, 1, :, n ) =                                           &
              real( cor_ins_du( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ), r_um )
          end do
        case default
          write( log_scratch_space, '(A,A)' )                                  &
                 'Missing required UKCA tracer field: ', ukca_tracer_names(n)
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end select
      end do
    end if  ! outer == outer_iterations .AND. l_tracer

    ! We do not want any sub-timestep SCM diagnostics but we still have
    ! to allocate some amount of memory to it as something is still trying
    ! to access it. Needless to say this will be caught by run-time checking.
    !
    l_scm_convss_dg = .false.
    allocate( scm_convss_dg(0) )

    ! Loop over convection calls per model time step
    do call_number = 1, n_conv_calls

      ! Initialise convection work arrays holding information from each
      ! iteration (sub-step)
      do k = 1, n_cca_lev
        do i = 1, ncells
          it_cca(i,1,k)  = 0.0_r_um
          it_cca0(i,1,k) = 0.0_r_um
          it_cca0_dp(i,1,k) = 0.0_r_um
          it_cca0_sh(i,1,k) = 0.0_r_um
          it_cca0_md(i,1,k) = 0.0_r_um
        end do
      end do

      do k = 1, nlayers
        do i = 1, ncells
          it_ccw(i,1,k)  = 0.0_r_um
          it_ccw0(i,1,k) = 0.0_r_um
          it_conv_rain_3d(i,1,k) = 0.0_r_um
          it_conv_snow_3d(i,1,k) = 0.0_r_um
          it_w2p(i,1,k) = 0.0_r_um
          it_up_flux(i,1,k) = 0.0_r_um
        end do
      end do

      do i = 1, ncells
        it_lcca(i,1)   = 0.0_r_um
        it_lcbase(i,1) = 0
        it_lctop(i,1)  = 0

        it_ccb(i,1)    = 0
        it_cct(i,1)    = 0
        it_cca_2d(i,1) = 0.0_r_um
        it_cclwp(i,1)  = 0.0_r_um

        it_ccb0(i,1)   = 0
        it_cct0(i,1)   = 0
        it_cclwp0(i,1) = 0.0_r_um

        it_conv_rain(i,1) = 0.0_r_um
        it_conv_snow(i,1) = 0.0_r_um
        it_precip_dp(i,1) = 0.0_r_um
        it_precip_sh(i,1) = 0.0_r_um
        it_precip_md(i,1) = 0.0_r_um
        it_cape_diluted(i,1)  = 0.0_r_um
        it_kterm_deep(i,1)  = 0
        it_kterm_shall(i,1) = 0
        it_mid_level(i,1) = .FALSE.
        it_dp_cfl_limited(i,1) = 0.0_r_um
        it_md_cfl_limited(i,1) = 0.0_r_um

        it_precip_cg(i,1) = 0.0_r_um
        it_wstar_up(i,1)  = 0.0_r_um
        it_mb1(i,1) = 0.0_r_um
        it_mb2(i,1) = 0.0_r_um
        it_cg_term(i,1) = 0
      end do

      if (prog_tnuc) then
        ! Use tnuc from LFRic and map onto tnuc_new for UM to be passed to glue_conv_6a
        do k = 1, nlayers
          do i = 1, ncells
            tnuc_new(i,1,k) = real(tnuc(map_wth(1,i) + k),kind=r_um)
          end do ! i
        end do ! k

        ! Use tnuc_nlcl from LFRic and map onto tnuc_nlcl_um for UM to the be passed to glue_conv_6a
        do i = 1, ncells
          tnuc_nlcl_um(i,1) = real(tnuc_nlcl(map_2d(1,i)),kind=r_um)
        end do
      end if

      nml_segment_size = min(conv_gr_segment_size,ncells)
      num_seg = ceiling(real(ncells,r_um)/real(nml_segment_size,r_um))

      do i = 1, num_seg
        ii = (i-1)*nml_segment_size + 1
        seg_size = min(ncells, ii+nml_segment_size-1) - ii + 1

        seg_num = map_wth(1,ii)  ! Only used by debugging error message from UM
                            ! This is probably most useful to indicate
                            ! which column has a problem

        n_deep = 0
        if (iconv_deep > 0 ) then
          do j = ii, ii+seg_size-1
            if (cumulus(j,1) .AND. .NOT. l_shallow(j,1) ) then
              n_deep = n_deep + 1
            endif
          end do
        end if

        n_shallow = 0
        if (iconv_shallow > 0 ) then
          do j = ii, ii+seg_size-1
            if (cumulus(j,1)  .AND. l_shallow(j,1) ) then
              n_shallow = n_shallow + 1
            endif
          end do
        end if

        n_congestus = 0
        do j = ii, ii+seg_size-1
            if (cumulus(j,1)) then
              n_congestus = n_congestus + 1   ! as UM though not actually using scheme
          end if
        end do

        n_mid = 0
        do j = ii, ii+seg_size-1
          if (l_mid(j,1)) then
            n_mid = n_mid + 1
          end if
        end do

        call glue_conv_6a                                                     &
          ( ncells, seg_size, n_conv_levels, n_wtrac, bl_levels               &
          , call_number, seg_num, theta_conv(ii,1,1), q_conv(ii,1,1)          &
          , qcl_conv(ii,1,1), qcf_conv(ii,1,1)                                &
          , q_wtrac(ii,1,1), qcl_wtrac(ii,1,1), qcf_wtrac(ii,1,1)             &
          , cf_liquid_conv(ii,1,1), cf_frozen_conv(ii,1,1)                    &
          , bulk_cf_conv(ii,1,1)                                              &
          , p_star(ii,1), land_sea_mask(ii,1)                                 &
          , u_conv(ii,1,1), v_conv(ii,1,1), w(ii,1,1)                         &
          , tot_tracer(ii,1,1,1), dthbydt(ii,1,1), dqbydt(ii,1,1)             &
          , dqclbydt(ii,1,1), dqcfbydt(ii,1,1)                                &
          , dcflbydt(ii,1,1), dcffbydt(ii,1,1), dbcfbydt(ii,1,1)              &
          , dubydt_p(ii,1,1), dvbydt_p(ii,1,1)                                &
          , dqbydt_wtrac(ii,1,1), dqclbydt_wtrac(ii,1,1)                      &
          , dqcfbydt_wtrac(ii,1,1)                                            &
          , it_conv_rain(ii,1), it_conv_snow(ii,1)                            &
          , it_conv_rain_3d(ii,1,1), it_conv_snow_3d(ii,1,1)                  &
          , rain_wtrac(ii,1), snow_wtrac(ii,1)                                &
          , it_cca0_dp(ii,1,1), it_cca0_md(ii,1,1), it_cca0_sh(ii,1,1)        &
          , it_cca0(ii,1,1),  it_ccb0(ii,1), it_cct0(ii,1)                    &
          , it_cclwp0(ii,1), it_ccw0(ii,1,1)                                  &
          , it_lcbase0(ii,1), it_lctop(ii,1),  it_lcca(ii,1)                  &
          , it_cca(ii,1,1),   it_ccb(ii,1),  it_cct(ii,1),  it_cclwp(ii,1)    &
          , it_ccw(ii,1,1),  it_lcbase(ii,1)                                  &
          , it_cca_2d(ii,1), freeze_lev(ii,1)                                 &
          , it_dp_cfl_limited(ii,1), it_md_cfl_limited(ii,1)                  &
          , it_mid_level(ii,1), it_kterm_deep(ii,1), it_kterm_shall(ii,1)     &
          , it_precip_dp(ii,1), it_precip_sh(ii,1)                            &
          , it_precip_md(ii,1), it_precip_cg(ii,1)                            &
          , it_wstar_dn(ii,1),  it_wstar_up(ii,1)                             &
          , it_mb1(ii,1), it_mb2(ii,1), it_cg_term(ii,1)                      &
          , uw0(ii,1), vw0(ii,1), w_max(ii,1)                                 &
          , zlcl(ii,1), zlcl_uv(ii,1), tnuc_new(ii,1,1), tnuc_nlcl_um(ii,1)   &
          , zhpar(ii,1), entrain_coef(ii,1)                                   &
          , conv_prog_precip_conv(ii,1,1), conv_prog_flx(ii,1,1)              &
          , deep_flag(ii,1)                                                   &
          , past_conv_ht(ii,1), it_cape_diluted(ii,1)                         &
          , n_deep, n_congestus, n_shallow, n_mid                             &
          , r_rho_levels(ii,1,1), r_theta_levels(ii,1,0)                      &
          , rho_wet(ii,1,1), rho_wet_tq(ii,1,1), rho_dry(ii,1,1)              &
          , rho_dry_theta(ii,1,1), delta_smag(ii,1)                           &
          , exner_rho_levels(ii,1,1), exner_rho_minus_one(ii,1,0)             &
          , exner_theta_levels(ii,1,0)                                        &
          , p_rho_minus_one(ii,1,0), p_theta_levels(ii,1,0)                   &
          , z_theta(ii,1,1), z_rho(ii,1,1), timestep_conv                     &
          , t1_sd(ii,1), q1_sd(ii,1), ntml(ii,1), ntpar(ii,1)                 &
          , conv_type(ii,1), l_shallow(ii,1)                                  &
          , l_congestus(ii,1), l_mid(ii,1), cumulus(ii,1)                     &
          , wstar(ii,1), wthvs(ii,1), delthvu(ii,1), ql_ad(ii,1)              &
          , qsat_lcl(ii,1), ftl(ii,1), fqw(ii,1)                              &
          , l_tracer, ntra_fld, ntra_lev, n_cca_lev                           &
          , l_calc_dxek , l_q_interact                                        &
          , it_up_flux_half(ii,1,1), it_up_flux(ii,1,1), it_dwn_flux(ii,1,1)  &
          , it_entrain_up(ii,1,1),   it_detrain_up(ii,1,1)                    &
          , it_entrain_dwn(ii,1,1),  it_detrain_dwn(ii,1,1)                   &
          , it_uw_dp(ii,1,1),        it_vw_dp(ii,1,1)                         &
          , it_uw_shall(ii,1,1),     it_vw_shall(ii,1,1)                      &
          , it_uw_mid(ii,1,1),  it_vw_mid(ii,1,1)                             &
          , it_wqt_flux(ii,1,1),  it_wthetal_flux(ii,1,1)                     &
          , it_wthetav_flux(ii,1,1), it_wql_flux(ii,1,1)                      &
          , it_mf_deep(ii,1,1),      it_mf_congest(ii,1,1)                    &
          , it_mf_shall(ii,1,1),  it_mf_midlev(ii,1,1)                        &
          , it_dt_deep(ii,1,1),      it_dt_congest(ii,1,1)                    &
          , it_dt_shall(ii,1,1),  it_dt_midlev(ii,1,1)                        &
          , it_dq_deep(ii,1,1),      it_dq_congest(ii,1,1)                    &
          , it_dq_shall(ii,1,1),  it_dq_midlev(ii,1,1)                        &
          , it_du_deep(ii,1,1),      it_du_congest(ii,1,1)                    &
          , it_du_shall(ii,1,1),  it_du_midlev(ii,1,1)                        &
          , it_dv_deep(ii,1,1),      it_dv_congest(ii,1,1)                    &
          , it_dv_shall(ii,1,1),  it_dv_midlev(ii,1,1)                        &
          , ind_cape_reduced(ii,1),  cape_ts_used(ii,1)                       &
          , it_ind_deep(ii,1), it_ind_shall(ii,1)                             &
          , it_w2p(ii,1,1), it_dt_dd(ii,1,1), it_dq_dd(ii,1,1)                &
          , it_du_dd(ii,1,1), it_dv_dd(ii,1,1)                                &
          , it_area_ud(ii,1,1), it_area_dd(ii,1,1)                            &
          , scm_convss_dg, l_scm_convss_dg                                    &
          , g_ccp(ii,1), h_ccp(ii,1), ccp_strength(ii,1)                      &
          )
      end do ! i -> num_seg

      ! Mid-level convection only possible on subsequent sub-steps if
      ! occurs on first step or column is diagnosed as cumulus
      do i = 1, ncells
        l_mid(i,1) = it_mid_level(i,1) .or. cumulus(i,1)
      end do

      ! Update cloud info from substep
      ! Highest convective layer properties - diagnostic
      !------------------------------------
      ! max cct across total number of calls to convection
      ! Note that diagnostic is a real not an integer
      do i = 1, ncells
        cct(i,1) = max( cct(i,1),it_cct(i,1))
        ! min ccb across total number of calls to convection
        ! excluding ccb=0
        if (ccb(i,1) > 0 .AND. it_ccb(i,1) > 0) then
          ccb(i,1) = min(ccb(i,1),it_ccb(i,1))
        else
          ccb(i,1) = max(ccb(i,1),it_ccb(i,1))
        end if
      end do

      ! Lowest convective layer properties
      !------------------------------------
      ! max lctop across total number of calls to convection
      do i = 1, ncells
        lctop(i,1) = max(lctop(i,1),it_lctop(i,1))

        ! min lcbase across total number of calls to convection
        ! excluding lcbase=0
        if (lcbase(i,1) > 0 .AND. it_lcbase(i,1) > 0) then
          lcbase(i,1) = min(lcbase(i,1), it_lcbase(i,1))
        else
          lcbase(i,1) = max(lcbase(i,1), it_lcbase(i,1))
        end if
      end do

      do k=1, nlayers
        do i = 1, ncells
          ccw_3d(i,1,k) = ccw_3d(i,1,k) + one_over_conv_calls*it_ccw0(i,1,k)
          cca_3d(i,1,k) = cca_3d(i,1,k) + one_over_conv_calls*it_cca0(i,1,k)
          ! Assuming lccrad = .true.
          cca_3d(i,1,k) = min(cca_3d(i,1,k), 1.0_r_um)
        end do
      end do

      ! single level convection diagnostics
      do i = 1, ncells
        conv_rain(map_2d(1,i)) = conv_rain(map_2d(1,i))+                          &
                                 it_conv_rain(i,1) *one_over_conv_calls
        conv_snow(map_2d(1,i)) = conv_snow(map_2d(1,i))+                          &
                                 it_conv_snow(i,1) *one_over_conv_calls
        cca_2d(map_2d(1,i)) = cca_2d(map_2d(1,i))+                                &
                                 it_cca_2d(i,1) *one_over_conv_calls
        cape_diluted(map_2d(1,i)) = cape_diluted(map_2d(1,i)) +                   &
                                it_cape_diluted(i,1)*one_over_conv_calls
      end do

      if (outer == outer_iterations) then
        if (.not. associated(lowest_cca_2d, empty_real_data) ) then
          do i = 1, ncells
            lowest_cca_2d(map_2d(1,i)) = lowest_cca_2d(map_2d(1,i)) +               &
                                       it_lcca(i,1)*one_over_conv_calls
          end do
        end if
        if (.not. associated(deep_in_col, empty_real_data) ) then
          do i = 1, ncells
            deep_in_col(map_2d(1,i)) = deep_in_col(map_2d(1,i)) +                   &
                                     it_ind_deep(i,1)*one_over_conv_calls
          end do
        end if
        if (.not. associated(shallow_in_col, empty_real_data) ) then
          do i = 1, ncells
            shallow_in_col(map_2d(1,i)) = shallow_in_col(map_2d(1,i)) +             &
                                        it_ind_shall(i,1)*one_over_conv_calls
          end do
        end if
        if (.not. associated(mid_in_col, empty_real_data) ) then
          do i = 1, ncells
            if (it_mid_level(i,1)) then
              mid_in_col(map_2d(1,i)) = mid_in_col(map_2d(1,i)) + one_over_conv_calls
            end if
          end do
        end if
        if (.not. associated(freeze_level, empty_real_data) ) then
          do i = 1, ncells
            freeze_level(map_2d(1,i)) = freeze_level(map_2d(1,i)) +                 &
                                    real(freeze_lev(i,1)) *one_over_conv_calls
          end do
        end if
        if (.not. associated(deep_prec, empty_real_data) ) then
          do i = 1, ncells
            deep_prec(map_2d(1,i)) = deep_prec(map_2d(1,i)) +                       &
                                   it_precip_dp(i,1) *one_over_conv_calls
          end do
        end if
        if (.not. associated(shallow_prec, empty_real_data) ) then
          do i = 1, ncells
            shallow_prec(map_2d(1,i)) = shallow_prec(map_2d(1,i)) +                 &
                                      it_precip_sh(i,1) *one_over_conv_calls
          end do
        end if
        if (.not. associated(mid_prec, empty_real_data) ) then
          do i = 1, ncells
            mid_prec(map_2d(1,i)) = mid_prec(map_2d(1,i)) +                         &
                                it_precip_md(i,1) *one_over_conv_calls
          end do
        end if
        if (.not. associated(deep_term, empty_real_data) ) then
          do i = 1, ncells
            deep_term(map_2d(1,i)) = deep_term(map_2d(1,i)) +                       &
                                 real(it_kterm_deep(i,1)) *one_over_conv_calls
          end do
        end if
        if (.not. associated(cape_timescale, empty_real_data) ) then
          do i = 1, ncells
            cape_timescale(map_2d(1,i)) =cape_timescale(map_2d(1,i)) +              &
                                     cape_ts_used(i,1) *one_over_conv_calls
          end do
        end if
        if (.not. associated(deep_cfl_limited, empty_real_data) ) then
          do i = 1, ncells
            deep_cfl_limited(map_2d(1,i)) = deep_cfl_limited(map_2d(1,i)) +         &
                                  it_dp_cfl_limited(i,1) *one_over_conv_calls
          end do
        end if
        if (.not. associated(mid_cfl_limited, empty_real_data) ) then
          do i = 1, ncells
          mid_cfl_limited(map_2d(1,i)) = mid_cfl_limited(map_2d(1,i)) +             &
                                it_md_cfl_limited(i,1) *one_over_conv_calls
          end do
        end if

        ! Frequency of deep convection terminating on level k
        if (.not. associated(deep_tops, empty_real_data) ) then
          k = 0  ! dummy initialisation for psyclone firstprivate bug and cray compiler
          do i = 1, ncells
            if (it_ind_deep(i,1) == 1.0_r_um) then
              k = it_kterm_deep(i,1)
              if (k > 0) then  ! in case still get a zero value
                deep_tops(map_wth(1,i)+k) = deep_tops(map_wth(1,i)+k) + one_over_conv_calls
              end if
            end if
          end do
        end if
      end if ! outer_iterations

      ! update input fields *_conv for next substep
      do k = 1, n_conv_levels
        do i = 1, ncells
          theta_conv(i,1,k) = theta_conv(i,1,k)                         &
                              + dthbydt(i,1,k) * timestep_conv
          q_conv(i,1,k)     = q_conv(i,1,k)                             &
                              + dqbydt(i,1,k) * timestep_conv
          qcl_conv(i,1,k)   = qcl_conv(i,1,k)                           &
                              +(dqclbydt(i,1,k) * timestep_conv)
          qcf_conv(i,1,k)   = qcf_conv(i,1,k)                           &
                              +(dqcfbydt(i,1,k) * timestep_conv)
          cf_liquid_conv(i,1,k) = cf_liquid_conv(i,1,k)                 &
                                  +(dcflbydt(i,1,k) * timestep_conv)
          cf_frozen_conv(i,1,k) = cf_frozen_conv(i,1,k)                 &
                                  + (dcffbydt(i,1,k) * timestep_conv)
          bulk_cf_conv(i,1,k)   = bulk_cf_conv(i,1,k)                   &
                                  +(dbcfbydt(i,1,k) * timestep_conv)
          dtheta_conv(i,1,k)   = dtheta_conv(i,1,k)                     &
                                      + dthbydt(i,1,k) * timestep_conv
          dt_conv(map_wth(1,i) + k)   = dt_conv(map_wth(1,i) + k)       &
                                      + dthbydt(i,1,k) * timestep_conv  &
                                      * exner_theta_levels(i,1,k)
          dmv_conv(map_wth(1,i) + k)  =  dmv_conv(map_wth(1,i) + k)     &
                                       + dqbydt(i,1,k) * timestep_conv
          dmcl_conv(map_wth(1,i) + k) =  dmcl_conv(map_wth(1,i) + k)    &
                                       + dqclbydt(i,1,k) * timestep_conv
          dms_conv(map_wth(1,i) + k) =  dms_conv(map_wth(1,i) + k)      &
                                       + dqcfbydt(i,1,k) * timestep_conv

          ! Update diagnostics
          massflux_up(map_wth(1,i) + k) = massflux_up(map_wth(1,i) + k) +          &
                                         it_up_flux(i,1,k)*one_over_conv_calls
          massflux_down(map_wth(1,i) + k) = massflux_down(map_wth(1,i) + k)+       &
                                         it_dwn_flux(i,1,k)*one_over_conv_calls

          conv_rain_3d(map_wth(1,i) + k) = conv_rain_3d(map_wth(1,i) + k) +        &
                                         it_conv_rain_3d(i,1,k) *                  &
                                         one_over_conv_calls
          conv_snow_3d(map_wth(1,i) + k) = conv_snow_3d(map_wth(1,i) + k) +        &
                                         it_conv_snow_3d(i,1,k) *                  &
                                         one_over_conv_calls
        end do
      end do

      ! Update optional diagnostics
      if (outer == outer_iterations) then
        if (.not. associated(entrain_up, empty_real_data) ) then
          do k = 1, n_conv_levels
            do i = 1, ncells
              entrain_up(map_wth(1,i) + k) = entrain_up(map_wth(1,i) + k) +          &
                                           it_entrain_up(i,1,k)*one_over_conv_calls
            end do
          end do
        end if
        if (.not. associated(entrain_down, empty_real_data) ) then
          do k = 1, n_conv_levels
            do i = 1, ncells
              entrain_down(map_wth(1,i) + k) = entrain_down(map_wth(1,i) + k) +      &
                                           it_entrain_dwn(i,1,k)*one_over_conv_calls
            end do
          end do
        end if
        if (.not. associated(detrain_up, empty_real_data) ) then
          do k = 1, n_conv_levels
            do i = 1, ncells
              detrain_up(map_wth(1,i) + k) = detrain_up(map_wth(1,i) + k) +          &
                                           it_detrain_up(i,1,k)*one_over_conv_calls
            end do
          end do
        end if
        if (.not. associated(detrain_down, empty_real_data) ) then
          do k = 1, n_conv_levels
            do i = 1, ncells
              detrain_down(map_wth(1,i) + k) = detrain_down(map_wth(1,i) + k) +      &
                                           it_detrain_dwn(i,1,k)*one_over_conv_calls
            end do
          end do
        end if
        if (.not. associated(dd_dt, empty_real_data) ) then
          do k = 1, n_conv_levels
            do i = 1, ncells
              dd_dt(map_wth(1,i) + k) = dd_dt(map_wth(1,i) + k) +                    &
                                           it_dt_dd(i,1,k)*one_over_conv_calls
            end do
          end do
        end if
        if (.not. associated(dd_dq, empty_real_data) ) then
          do k = 1, n_conv_levels
            do i = 1, ncells
              dd_dq(map_wth(1,i) + k) = dd_dq(map_wth(1,i) + k) +                    &
                                           it_dq_dd(i,1,k)*one_over_conv_calls
            end do
          end do
        end if
        if (.not. associated(deep_massflux, empty_real_data) ) then
          do k = 1, n_conv_levels
            do i = 1, ncells
              deep_massflux(map_wth(1,i) + k) = deep_massflux(map_wth(1,i) + k) +    &
                                           it_mf_deep(i,1,k)*one_over_conv_calls
            end do
          end do
        end if
        if (.not. associated(deep_dt, empty_real_data) ) then
          do k = 1, n_conv_levels
            do i = 1, ncells
              deep_dt(map_wth(1,i) + k) = deep_dt(map_wth(1,i) + k) +                &
                                           it_dt_deep(i,1,k)*one_over_conv_calls
            end do
          end do
        end if
        if (.not. associated(deep_dq, empty_real_data) ) then
          do k = 1, n_conv_levels
            do i = 1, ncells
               deep_dq(map_wth(1,i) + k) = deep_dq(map_wth(1,i) + k) +                &
                                           it_dq_deep(i,1,k)*one_over_conv_calls
           end do
          end do
        end if
        if (.not. associated(shallow_massflux, empty_real_data) ) then
          do k = 1, n_conv_levels
            do i = 1, ncells
              shallow_massflux(map_wth(1,i) + k) = shallow_massflux(map_wth(1,i) + k) + &
                                           it_mf_shall(i,1,k)*one_over_conv_calls
            end do
          end do
        end if
        if (.not. associated(shallow_dt, empty_real_data) ) then
          do k = 1, n_conv_levels
            do i = 1, ncells
              shallow_dt(map_wth(1,i) + k) = shallow_dt(map_wth(1,i) + k) +          &
                                           it_dt_shall(i,1,k)*one_over_conv_calls
            end do
          end do
        end if
        if (.not. associated(shallow_dq, empty_real_data) ) then
          do k = 1, n_conv_levels
            do i = 1, ncells
              shallow_dq(map_wth(1,i) + k) = shallow_dq(map_wth(1,i) + k) +          &
                                           it_dq_shall(i,1,k)*one_over_conv_calls
            end do
          end do
        end if
        if (.not. associated(mid_massflux, empty_real_data) ) then
          do k = 1, n_conv_levels
            do i = 1, ncells
              mid_massflux(map_wth(1,i) + k) = mid_massflux(map_wth(1,i) + k) +      &
                                           it_mf_midlev(i,1,k)*one_over_conv_calls
            end do
          end do
        end if
        if (.not. associated(mid_dt, empty_real_data) ) then
          do k = 1, n_conv_levels
            do i = 1, ncells
              mid_dt(map_wth(1,i) + k) = mid_dt(map_wth(1,i) + k) +                  &
                                           it_dt_midlev(i,1,k)*one_over_conv_calls
            end do
          end do
        end if
        if (.not. associated(mid_dq, empty_real_data) ) then
          do k = 1, n_conv_levels
            do i = 1, ncells
              mid_dq(map_wth(1,i) + k) = mid_dq(map_wth(1,i) + k) +                  &
                                           it_dq_midlev(i,1,k)*one_over_conv_calls
            end do
          end do
        end if
        if (.not. associated(cca_unadjusted, empty_real_data) ) then
          do k = 1, n_conv_levels
            do i = 1, ncells
              cca_unadjusted(map_wth(1,i) + k) = cca_unadjusted(map_wth(1,i) + k) +  &
                                           it_cca(i,1,k)*one_over_conv_calls
            end do
          end do
        end if
        if (.not. associated(massflux_up_half, empty_real_data) ) then
          do k = 1, n_conv_levels
            do i = 1, ncells
            massflux_up_half(map_w3(1,i) + k-1) = massflux_up_half(map_w3(1,i) + k-1) +&
                                         it_up_flux_half(i,1,k)*one_over_conv_calls
            end do
          end do
        end if
      end if ! outer_iterations

      if (l_mom) then
        do k = 1, n_conv_levels
          do i = 1, ncells
            u_conv(i,1,k) = u_conv(i,1,k)   + dubydt_p(i,1,k) * timestep_conv
            v_conv(i,1,k) = v_conv(i,1,k)   + dvbydt_p(i,1,k) * timestep_conv
            ! total increments
            du_conv(map_w3(1,i) + k -1) = du_conv(map_w3(1,i) + k -1) + dubydt_p(i,1,k) * timestep_conv
            dv_conv(map_w3(1,i) + k -1) = dv_conv(map_w3(1,i) + k -1) + dvbydt_p(i,1,k) * timestep_conv
          end do
        end do
      end if    !l_mom

      ! Would have PC2 checks here

    end do    ! loop over calls to convection

    if (l_safe_conv) then
      do k = 1, n_conv_levels
        do i = 1, ncells
          dmv_conv(map_wth(1,i) + k) = dmv_conv(map_wth(1,i) + k) + dq_add(i,1,k)
        end do
      end do
    end if

    ! Update the time-smoothed convection prognostics
    if (l_conv_prog_precip .and. outer == outer_iterations) then
      decay_amount            = timestep / tau_conv_prog_precip
      theta_inc_threshold     = dthetadt_conv_active_threshold * timestep_conv
      do i = 1, ncells
        tot_conv_precip_2d(i,1) = MAX(conv_rain(map_2d(1,i)) +        &
                                  conv_snow(map_2d(1,i)),             &
                                  conv_prog_precip_min_threshold )
      end do
      conv_active = 0.0_r_um  ! dummy initialisation for psyclone firstprivate bug and cray compiler
      do k = 1, n_conv_levels
        do i = 1, ncells
          if (abs(dtheta_conv(i,1,k)) > theta_inc_threshold) then
            conv_active = 1.0_r_um
          else
            conv_active = 0.0_r_um
          end if
          conv_prog_precip(map_wth(1,i) + k)                                      &
              = decay_amount * tot_conv_precip_2d(i,1) * conv_active              &
              + (1.0_r_def - decay_amount) * conv_prog_precip(map_wth(1,i) + k)
        end do
      end do
      do i = 1, ncells
          conv_prog_precip(map_wth(1,i) + 0) = conv_prog_precip(map_wth(1,i) + 1)
      end do
    end if

    if (l_conv_prog_dtheta) then
      decay_amount = timestep / tau_conv_prog_dtheta
      do k = 1, n_conv_levels
        do i = 1, ncells
          dt_conv(map_wth(1,i) + k)  = (decay_amount * dtheta_conv(i,1,k)                &
                      + (1.0_r_def - decay_amount) * conv_prog_dtheta(map_wth(1,i) + k)) &
                      * exner_in_wth(map_wth(1,i) + k)
        end do
      end do

      if (outer == outer_iterations) then
        do k = 1, n_conv_levels
          do i = 1, ncells
            conv_prog_dtheta(map_wth(1,i) + k) = dt_conv(map_wth(1,i) + k)         &
                                             / exner_in_wth(map_wth(1,i) + k)
          end do
        end do
        do i = 1, ncells
          conv_prog_dtheta(map_wth(1,i) + 0) = conv_prog_dtheta(map_wth(1,i) + 1)
        end do
      end if
    end if

    if (l_conv_prog_dq) then
      decay_amount = timestep / tau_conv_prog_dq
      do k = 1, n_conv_levels
        do i = 1, ncells
          dmv_conv(map_wth(1,i) + k) = decay_amount * dmv_conv(map_wth(1,i) + k)   &
                       + (1.0_r_def - decay_amount) * conv_prog_dmv(map_wth(1,i) + k)
        end do
      end do

      if (outer == outer_iterations) then
        do k = 1, n_conv_levels
          do i = 1, ncells
            conv_prog_dmv(map_wth(1,i) + k) = dmv_conv(map_wth(1,i) + k)
          end do
        end do
        do i = 1, ncells
          conv_prog_dmv(map_wth(1,i) + 0) = conv_prog_dmv(map_wth(1,i) + 1)
        end do
      end if
    end if

    ! Convection/PC2 checks
    ! Protect against generation of inconsistently low cloud
    ! fraction implying very high in-cloud condensate amounts.
    ! In-cloud condensate amounts above 2.0e-3 lead to
    ! cloud fraction being increased (up to a value of 1.0)

    orig_value = 0.0_r_um ! dummy initialisation for psyclone firstprivate bug and cray compiler
    do k = 1, n_conv_levels
      do i = 1, ncells
        ! Liquid cloud fraction
        if (cf_liquid_conv(i,1,k) > 0.0_r_um) then
          if ( (qcl_conv(i,1,k)/cf_liquid_conv(i,1,k) ) > 2.0e-3_r_um ) then
            orig_value = cf_liquid_conv(i,1,k)
            cf_liquid_conv(i,1,k) = min(1.0_r_um,qcl_conv(i,1,k)/2.0e-3_r_um)
            bulk_cf_conv(i,1,k) = bulk_cf_conv(i,1,k)                      &
                                  + cf_liquid_conv(i,1,k) - orig_value
          end if
        end if

        ! Ice cloud fraction
        if (cf_frozen_conv(i,1,k) > 0.0_r_um) then
          if ( (qcf_conv(i,1,k)/cf_frozen_conv(i,1,k)) > 2.0e-3_r_um ) then
            orig_value = cf_frozen_conv(i,1,k)
            cf_frozen_conv(i,1,k) = min(1.0_r_um,qcf_conv(i,1,k)/2.0e-3_r_um)
            bulk_cf_conv(i,1,k) = bulk_cf_conv(i,1,k)                     &
                                  + cf_frozen_conv(i,1,k) - orig_value
          end if
        end if
      end do
    end do

    ! Store cloud fraction increments for adding on later if using PC2
    do k = 1, n_conv_levels
      do i = 1, ncells
        dcfl_conv(map_wth(1,i) + k) = cf_liquid_conv(i,1,k) - cf_liq(map_wth(1,i) + k)
        dcff_conv(map_wth(1,i) + k) = cf_frozen_conv(i,1,k) - cf_ice(map_wth(1,i) + k)
        dbcf_conv(map_wth(1,i) + k) = bulk_cf_conv(i,1,k)   - cf_bulk(map_wth(1,i) + k)
      end do
    end do
    do i = 1, ncells
      dcfl_conv(map_wth(1,i) + 0) = dcfl_conv(map_wth(1,i) + 1)
      dcff_conv(map_wth(1,i) + 0) = dcff_conv(map_wth(1,i) + 1)
      dbcf_conv(map_wth(1,i) + 0) = dbcf_conv(map_wth(1,i) + 1)
    end do

    ! Set level 0 increment such that theta increment will equal level 1
    do i = 1, ncells
      dt_conv (map_wth(1,i) + 0) = dt_conv  (map_wth(1,i) + 1)    &
                               * exner_in_wth(map_wth(1,i) + 0)   &
                               / exner_in_wth(map_wth(1,i) + 1)
      dmv_conv (map_wth(1,i) + 0) = dmv_conv (map_wth(1,i) + 1)
      dmcl_conv(map_wth(1,i) + 0) = dmcl_conv(map_wth(1,i) + 1)
      dms_conv(map_wth(1,i) + 0) = dms_conv(map_wth(1,i) + 1)
    end do

    ! Store convective downdraught mass fluxes at cloud base
    ! if required for surface exchange.
    if (srf_ex_cnv_gust == ip_srfexwithcnv) then
      do i = 1, ncells
        if (ccb(i,1) > 0) then
          dd_mf_cb(map_2d(1,i)) = massflux_down( map_wth(1,i) + ccb(i,1))
        else
          dd_mf_cb(map_2d(1,i)) = 0.0_r_def
        end if
      end do
    end if

    ! copy convective cloud fraction into prognostic array
    do k = 1, n_conv_levels
      do i = 1, ncells
        cca(map_wth(1,i) + k) =  min(cca_3d(i,1,k), 1.0_r_um)
        ccw(map_wth(1,i) + k) =  ccw_3d(i,1,k)
      end do
    end do

    if (outer == outer_iterations) then
     ! Copy integers into real diagnostic arrays
      if (.not. associated(cv_top, empty_real_data) ) then
        do i = 1, ncells
          cv_top(map_2d(1,i))        = real(cct(i,1))
        end do
      end if
      if (.not. associated(cv_base, empty_real_data) ) then
        do i = 1, ncells
          cv_base(map_2d(1,i))       = real(ccb(i,1))
        end do
      end if
      if (.not. associated(lowest_cv_top, empty_real_data) ) then
        do i = 1, ncells
          lowest_cv_top(map_2d(1,i)) = real(lctop(i,1))
        end do
      end if
      if (.not. associated(lowest_cv_base, empty_real_data) ) then
        do i = 1, ncells
          lowest_cv_base(map_2d(1,i)) = real(lcbase(i,1))
        end do
      end if

      ! pressure at cv top/base
      if (.not. associated(pres_cv_top, empty_real_data) ) then
        do i = 1, ncells
          if (cct(i,1) > 0) then
            pres_cv_top(map_2d(1,i)) = p_rho_levels(i,1,cct(i,1))
          else
            pres_cv_top(map_2d(1,i)) = 0.0_r_def
          end if
        end do
      end if
      if (.not. associated(pres_cv_base, empty_real_data) ) then
        do i = 1, ncells
          if (ccb(i,1) > 0) then
            pres_cv_base(map_2d(1,i)) = p_rho_levels(i,1,ccb(i,1))
          else
            pres_cv_base(map_2d(1,i))= 0.0_r_def
          end if
        end do
      end if

      ! pressure at lowest cv top/base
      if (.not. associated(pres_lowest_cv_top, empty_real_data) ) then
        do i = 1, ncells
          if (lctop(i,1) > 0) then
            pres_lowest_cv_top(map_2d(1,i)) = p_rho_levels(i,1,lctop(i,1))
          else
            pres_lowest_cv_top(map_2d(1,i)) = 0.0_r_def
          end if
        end do
      end if
      if (.not. associated(pres_lowest_cv_base, empty_real_data) ) then
        do i = 1, ncells
          if (lcbase(i,1) > 0) then
            pres_lowest_cv_base(map_2d(1,i)) = p_rho_levels(i,1,lcbase(i,1))
          else
            pres_lowest_cv_base(map_2d(1,i))= 0.0_r_def
          end if
        end do
      end if

      ! component A of upward mass flux
      if (.not. associated(massflux_up_cmpta, empty_real_data) ) then
        do k = 1, n_conv_levels - 1
          do i = 1, ncells
            if ( (1.0_r_def - max_mf_fall) * massflux_up_half(map_w3(1,i) + k-1) < &
                                       massflux_up_half(map_w3(1,i) + k) ) then
              massflux_up_cmpta(map_w3(1,i) + k-1) = massflux_up_half(map_w3(1,i) + k-1) &
                                               * (1.0_r_def - shallow_in_col(map_2d(1,i)))
            else
              massflux_up_cmpta(map_w3(1,i) + k-1) = 0.0_r_def
            end if
          end do ! i
        end do ! k
      end if

      ! Convection theta increment without shallow for VAR
      if (.not. associated(dth_conv_noshal, empty_real_data) ) then
        do k = 0, nlayers
          do i = 1, ncells
            dth_conv_noshal(map_wth(1,i) + k) = dt_conv(map_wth(1,i) + k) / &
                                            exner_in_wth(map_wth(1,i) + k) * &
                                           (1.0_r_def - shallow_in_col(map_2d(1,i)))
          end do
        end do
      end if

      ! Convection mixing ratio increment without shallow for VAR
      if (.not. associated(dmv_conv_noshal, empty_real_data) ) then
        do k = 0, nlayers
          do i = 1, ncells
            dmv_conv_noshal(map_wth(1,i) + k) = dmv_conv(map_wth(1,i) + k) * &
                                           (1.0_r_def - shallow_in_col(map_2d(1,i)))
          end do
        end do
      end if

      ! provide some estimate of TKE in convective plumes, based on
      ! the mass flux and convective cloud area
      do k = 1, bl_levels
        do i = 1, ncells
          tke_bl(map_wth(1,i)+k) = MIN(max_tke,MAX(tke_bl(map_wth(1,i)+k),         &
                 ( massflux_up(map_wth(1,i)+k) / ( g*rho_wet_tq(i,1,k)*            &
                   MIN(0.5_r_um,MAX(0.05_r_um,cca_2d(map_2d(1,i)))) ) )**2))
                   ! 0.5 and 0.05 are used here as plausible max and min
                   ! values of CCA to prevent numerical problems
        end do
      end do

    end if ! outer_iterations

    ! Copy tracers back to LFRic fields
    if ( outer == outer_iterations .AND. l_tracer ) then
      do n = 1, ntra_fld
        select case(ukca_tracer_names(n))
         case(fldname_o3p)
          do i = 1, ncells
            o3p( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            o3p( map_wth(1,i) + 0 ) = o3p( map_wth(1,i) + 1 )
          end do
         case(fldname_o1d)
          do i = 1, ncells
            o1d( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            o1d( map_wth(1,i) + 0 ) = o1d( map_wth(1,i) + 1 )
          end do
         case(fldname_o3)
          do i = 1, ncells
            o3( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            o3( map_wth(1,i) + 0 ) = o3( map_wth(1,i) + 1 )
          end do
         case(fldname_n)
          do i = 1, ncells
            nit( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            nit( map_wth(1,i) + 0 ) = nit( map_wth(1,i) + 1 )
          end do
         case(fldname_no)
          do i = 1, ncells
            no( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            no( map_wth(1,i) + 0 ) = no( map_wth(1,i) + 1 )
          end do
         case(fldname_no3)
          do i = 1, ncells
            no3( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            no3( map_wth(1,i) + 0 ) = no3( map_wth(1,i) + 1 )
          end do
         case(fldname_lumped_n)
          do i = 1, ncells
            lumped_n( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            lumped_n( map_wth(1,i) + 0 ) = lumped_n( map_wth(1,i) + 1 )
          end do
         case(fldname_n2o5)
          do i = 1, ncells
            n2o5( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            n2o5( map_wth(1,i) + 0 ) = n2o5( map_wth(1,i) + 1 )
          end do
         case(fldname_ho2no2)
          do i = 1, ncells
            ho2no2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            ho2no2( map_wth(1,i) + 0 ) = ho2no2( map_wth(1,i) + 1 )
          end do
         case(fldname_hono2)
          do i = 1, ncells
            hono2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            hono2( map_wth(1,i) + 0 ) = hono2( map_wth(1,i) + 1 )
          end do
        case(fldname_h2o2)
          do i = 1, ncells
            h2o2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            h2o2( map_wth(1,i) + 0 ) = h2o2( map_wth(1,i) + 1 )
          end do
         case(fldname_ch4)
          do i = 1, ncells
            ch4( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            ch4( map_wth(1,i) + 0 ) = ch4( map_wth(1,i) + 1 )
          end do
         case(fldname_co)
          do i = 1, ncells
            co( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            co( map_wth(1,i) + 0 ) = co( map_wth(1,i) + 1 )
          end do
         case(fldname_hcho)
          do i = 1, ncells
            hcho( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            hcho( map_wth(1,i) + 0 ) = hcho( map_wth(1,i) + 1 )
          end do
         case(fldname_meoo)
          do i = 1, ncells
            meoo( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            meoo( map_wth(1,i) + 0 ) = meoo( map_wth(1,i) + 1 )
          end do
         case(fldname_meooh)
          do i = 1, ncells
            meooh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                    &
              real( tot_tracer( i, 1, :, n ), r_def )
            meooh( map_wth(1,i) + 0 ) = meooh( map_wth(1,i) + 1 )
          end do
         case(fldname_h)
          do i = 1, ncells
            h( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            h( map_wth(1,i) + 0 ) = h( map_wth(1,i) + 1 )
          end do
         case(fldname_ch2o)
          ! H2O tracer from chemistry ns not transported
         case(fldname_oh)
          do i = 1, ncells
            oh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            oh( map_wth(1,i) + 0 ) = oh( map_wth(1,i) + 1 )
          end do
         case(fldname_ho2)
          do i = 1, ncells
            ho2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            ho2( map_wth(1,i) + 0 ) = ho2( map_wth(1,i) + 1 )
          end do
         case(fldname_cl)
          do i = 1, ncells
            cl( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            cl( map_wth(1,i) + 0 ) = cl( map_wth(1,i) + 1 )
          end do
         case(fldname_cl2o2)
          do i = 1, ncells
            cl2o2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                   &
              real( tot_tracer( i, 1, :, n ), r_def )
            cl2o2( map_wth(1,i) + 0 ) = cl2o2( map_wth(1,i) + 1 )
          end do
         case(fldname_clo)
          do i = 1, ncells
            clo( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            clo( map_wth(1,i) + 0 ) = clo( map_wth(1,i) + 1 )
          end do
         case(fldname_oclo)
          do i = 1, ncells
            oclo( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            oclo( map_wth(1,i) + 0 ) = oclo( map_wth(1,i) + 1 )
          end do
         case(fldname_br)
          do i = 1, ncells
            br( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            br( map_wth(1,i) + 0 ) = br( map_wth(1,i) + 1 )
          end do
         case(fldname_lumped_br)
          do i = 1, ncells
            lumped_br( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            lumped_br( map_wth(1,i) + 0 ) = lumped_br( map_wth(1,i) + 1 )
          end do
         case(fldname_brcl)
          do i = 1, ncells
            brcl( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            brcl( map_wth(1,i) + 0 ) = brcl( map_wth(1,i) + 1 )
          end do
         case(fldname_brono2)
          do i = 1, ncells
            brono2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            brono2( map_wth(1,i) + 0 ) = brono2( map_wth(1,i) + 1 )
          end do
         case(fldname_n2o)
          do i = 1, ncells
            n2o( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            n2o( map_wth(1,i) + 0 ) = n2o( map_wth(1,i) + 1 )
          end do
         case(fldname_lumped_cl)
          do i = 1, ncells
            lumped_cl( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            lumped_cl( map_wth(1,i) + 0 ) = lumped_cl( map_wth(1,i) + 1 )
          end do
         case(fldname_hocl)
          do i = 1, ncells
            hocl( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            hocl( map_wth(1,i) + 0 ) = hocl( map_wth(1,i) + 1 )
          end do
         case(fldname_hbr)
          do i = 1, ncells
            hbr( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            hbr( map_wth(1,i) + 0 ) = hbr( map_wth(1,i) + 1 )
          end do
         case(fldname_hobr)
          do i = 1, ncells
            hobr( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            hobr( map_wth(1,i) + 0 ) = hobr( map_wth(1,i) + 1 )
          end do
         case(fldname_clono2)
          do i = 1, ncells
            clono2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            clono2( map_wth(1,i) + 0 ) = clono2( map_wth(1,i) + 1 )
          end do
         case(fldname_cfcl3)
          do i = 1, ncells
            cfcl3( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            cfcl3( map_wth(1,i) + 0 ) = cfcl3( map_wth(1,i) + 1 )
          end do
         case(fldname_cf2cl2)
          do i = 1, ncells
            cf2cl2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            cf2cl2( map_wth(1,i) + 0 ) = cf2cl2( map_wth(1,i) + 1 )
          end do
         case(fldname_mebr)
          do i = 1, ncells
            mebr( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            mebr( map_wth(1,i) + 0 ) = mebr( map_wth(1,i) + 1 )
          end do
         case(fldname_hono)
          do i = 1, ncells
            hono( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            hono( map_wth(1,i) + 0 ) = hono( map_wth(1,i) + 1 )
          end do
         case(fldname_c2h6)
          do i = 1, ncells
            c2h6( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            c2h6( map_wth(1,i) + 0 ) = c2h6( map_wth(1,i) + 1 )
          end do
         case(fldname_etoo)
          do i = 1, ncells
            etoo( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            etoo( map_wth(1,i) + 0 ) = etoo( map_wth(1,i) + 1 )
          end do
         case(fldname_etooh)
          do i = 1, ncells
            etooh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            etooh( map_wth(1,i) + 0 ) = etooh( map_wth(1,i) + 1 )
          end do
         case(fldname_mecho)
          do i = 1, ncells
            mecho( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            mecho( map_wth(1,i) + 0 ) = mecho( map_wth(1,i) + 1 )
          end do
         case(fldname_meco3)
          do i = 1, ncells
            meco3( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            meco3( map_wth(1,i) + 0 ) = meco3( map_wth(1,i) + 1 )
          end do
         case(fldname_pan)
          do i = 1, ncells
            pan( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            pan( map_wth(1,i) + 0 ) = pan( map_wth(1,i) + 1 )
          end do
         case(fldname_c3h8)
          do i = 1, ncells
            c3h8( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            c3h8( map_wth(1,i) + 0 ) = c3h8( map_wth(1,i) + 1 )
          end do
         case(fldname_n_proo)
          do i = 1, ncells
            n_proo( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            n_proo( map_wth(1,i) + 0 ) = n_proo( map_wth(1,i) + 1 )
          end do
         case(fldname_i_proo)
          do i = 1, ncells
            i_proo( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            i_proo( map_wth(1,i) + 0 ) = i_proo( map_wth(1,i) + 1 )
          end do
         case(fldname_n_prooh)
          do i = 1, ncells
            n_prooh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            n_prooh( map_wth(1,i) + 0 ) = n_prooh( map_wth(1,i) + 1 )
          end do
         case(fldname_i_prooh)
          do i = 1, ncells
            i_prooh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            i_prooh( map_wth(1,i) + 0 ) = i_prooh( map_wth(1,i) + 1 )
          end do
         case(fldname_etcho)
          do i = 1, ncells
            etcho( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            etcho( map_wth(1,i) + 0 ) = etcho( map_wth(1,i) + 1 )
          end do
         case(fldname_etco3)
          do i = 1, ncells
            etco3( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            etco3( map_wth(1,i) + 0 ) = etco3( map_wth(1,i) + 1 )
         end do
         case(fldname_me2co)
          do i = 1, ncells
            me2co( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            me2co( map_wth(1,i) + 0 ) = me2co( map_wth(1,i) + 1 )
          end do
         case(fldname_mecoch2oo)
          do i = 1, ncells
            mecoch2oo( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            mecoch2oo( map_wth(1,i) + 0 ) = mecoch2oo( map_wth(1,i) + 1 )
          end do
         case(fldname_mecoch2ooh)
          do i = 1, ncells
            mecoch2ooh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            mecoch2ooh( map_wth(1,i) + 0 ) = mecoch2ooh( map_wth(1,i) + 1 )
          end do
         case(fldname_ppan)
          do i = 1, ncells
            ppan( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            ppan( map_wth(1,i) + 0 ) = ppan( map_wth(1,i) + 1 )
          end do
         case(fldname_meono2)
          do i = 1, ncells
            meono2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            meono2( map_wth(1,i) + 0 ) = meono2( map_wth(1,i) + 1 )
          end do
         case(fldname_c5h8)
          do i = 1, ncells
            c5h8( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            c5h8( map_wth(1,i) + 0 ) = c5h8( map_wth(1,i) + 1 )
          end do
         case(fldname_iso2)
          do i = 1, ncells
            iso2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            iso2( map_wth(1,i) + 0 ) = iso2( map_wth(1,i) + 1 )
          end do
         case(fldname_isooh)
          do i = 1, ncells
            isooh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            isooh( map_wth(1,i) + 0 ) = isooh( map_wth(1,i) + 1 )
          end do
         case(fldname_ison)
          do i = 1, ncells
            ison( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            ison( map_wth(1,i) + 0 ) = ison( map_wth(1,i) + 1 )
          end do
         case(fldname_macr)
          do i = 1, ncells
            macr( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            macr( map_wth(1,i) + 0 ) = macr( map_wth(1,i) + 1 )
          end do
         case(fldname_macro2)
          do i = 1, ncells
            macro2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            macro2( map_wth(1,i) + 0 ) = macro2( map_wth(1,i) + 1 )
          end do
         case(fldname_macrooh)
          do i = 1, ncells
            macrooh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            macrooh( map_wth(1,i) + 0 ) = macrooh( map_wth(1,i) + 1 )
          end do
         case(fldname_mpan)
          do i = 1, ncells
            mpan( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            mpan( map_wth(1,i) + 0 ) = mpan( map_wth(1,i) + 1 )
          end do
         case(fldname_hacet)
          do i = 1, ncells
            hacet( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            hacet( map_wth(1,i) + 0 ) = hacet( map_wth(1,i) + 1 )
          end do
         case(fldname_mgly)
          do i = 1, ncells
            mgly( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            mgly( map_wth(1,i) + 0 ) = mgly( map_wth(1,i) + 1 )
          end do
         case(fldname_nald)
          do i = 1, ncells
            nald( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            nald( map_wth(1,i) + 0 ) = nald( map_wth(1,i) + 1 )
          end do
         case(fldname_hcooh)
          do i = 1, ncells
            hcooh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            hcooh( map_wth(1,i) + 0 ) = hcooh( map_wth(1,i) + 1 )
          end do
         case(fldname_meco3h)
          do i = 1, ncells
            meco3h( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            meco3h( map_wth(1,i) + 0 ) = meco3h( map_wth(1,i) + 1 )
          end do
         case(fldname_meco2h)
          do i = 1, ncells
            meco2h( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            meco2h( map_wth(1,i) + 0 ) = meco2h( map_wth(1,i) + 1 )
          end do
         case(fldname_h2)
          do i = 1, ncells
            h2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            h2( map_wth(1,i) + 0 ) = h2( map_wth(1,i) + 1 )
          end do
         case(fldname_meoh)
          do i = 1, ncells
            meoh( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            meoh( map_wth(1,i) + 0 ) = meoh( map_wth(1,i) + 1 )
          end do
         case(fldname_msa)
          do i = 1, ncells
            msa( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            msa( map_wth(1,i) + 0 ) = msa( map_wth(1,i) + 1 )
          end do
         case(fldname_nh3)
          do i = 1, ncells
            nh3( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            nh3( map_wth(1,i) + 0 ) = nh3( map_wth(1,i) + 1 )
          end do
         case(fldname_cs2)
          do i = 1, ncells
            cs2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            cs2( map_wth(1,i) + 0 ) = cs2( map_wth(1,i) + 1 )
          end do
         case(fldname_csul)
          do i = 1, ncells
            csul( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            csul( map_wth(1,i) + 0 ) = csul( map_wth(1,i) + 1 )
          end do
         case(fldname_h2s)
          do i = 1, ncells
            h2s( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            h2s( map_wth(1,i) + 0 ) = h2s( map_wth(1,i) + 1 )
          end do
         case(fldname_so3)
          do i = 1, ncells
            so3( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            so3( map_wth(1,i) + 0 ) = so3( map_wth(1,i) + 1 )
          end do
         case(fldname_passive_o3)
          do i = 1, ncells
            passive_o3( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            passive_o3( map_wth(1,i) + 0 ) = passive_o3( map_wth(1,i) + 1 )
          end do
         case(fldname_age_of_air)
          do i = 1, ncells
            age_of_air( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            age_of_air( map_wth(1,i) + 0 ) = age_of_air( map_wth(1,i) + 1 )
          end do
        case(fldname_dms)
          do i = 1, ncells
            dms( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                      &
              real( tot_tracer( i, 1, :, n ), r_def )
            dms( map_wth(1,i) + 0 ) = dms( map_wth(1,i) + 1 )
          end do
        case(fldname_so2)
          do i = 1, ncells
            so2( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                      &
             real( tot_tracer( i, 1, :, n ), r_def )
            so2( map_wth(1,i) + 0 ) = so2( map_wth(1,i) + 1 )
          end do
        case(fldname_h2so4)
          do i = 1, ncells
            h2so4( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                    &
              real( tot_tracer( i, 1, :, n ), r_def )
            h2so4( map_wth(1,i) + 0 ) = h2so4( map_wth(1,i) + 1 )
          end do
        case(fldname_dmso)
          do i = 1, ncells
            dmso( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                     &
              real( tot_tracer( i, 1, :, n ), r_def )
            dmso( map_wth(1,i) + 0 ) = dmso( map_wth(1,i) + 1 )
          end do
        case(fldname_monoterpene)
          do i = 1, ncells
            monoterpene( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =              &
              real( tot_tracer( i, 1, :, n ), r_def )
            monoterpene( map_wth(1,i) + 0 ) = monoterpene( map_wth(1,i) + 1 )
          end do
        case(fldname_secondary_organic)
          do i = 1, ncells
            secondary_organic( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =        &
              real( tot_tracer( i, 1, :, n ), r_def )
            secondary_organic( map_wth(1,i) + 0 ) =                                &
              secondary_organic( map_wth(1,i) + 1 )
          end do
        case(fldname_n_nuc_sol)
          do i = 1, ncells
            n_nuc_sol( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                &
              real( tot_tracer( i, 1, :, n ), r_def )
            n_nuc_sol( map_wth(1,i) + 0 ) = n_nuc_sol( map_wth(1,i) + 1 )
          end do
        case(fldname_nuc_sol_su)
          do i = 1, ncells
            nuc_sol_su( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =               &
              real( tot_tracer( i, 1, :, n ), r_def )
            nuc_sol_su( map_wth(1,i) + 0 ) = nuc_sol_su( map_wth(1,i) + 1 )
          end do
        case(fldname_nuc_sol_om)
          do i = 1, ncells
            nuc_sol_om( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =               &
              real( tot_tracer( i, 1, :, n ), r_def )
            nuc_sol_om( map_wth(1,i) + 0 ) = nuc_sol_om( map_wth(1,i) + 1 )
          end do
        case(fldname_n_ait_sol)
          do i = 1, ncells
            n_ait_sol( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                &
              real( tot_tracer( i, 1, :, n ), r_def )
            n_ait_sol( map_wth(1,i) + 0 ) = n_ait_sol( map_wth(1,i) + 1 )
          end do
        case(fldname_ait_sol_su)
          do i = 1, ncells
            ait_sol_su( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =               &
              real( tot_tracer( i, 1, :, n ), r_def )
            ait_sol_su( map_wth(1,i) + 0 ) = ait_sol_su( map_wth(1,i) + 1 )
          end do
        case(fldname_ait_sol_bc)
          do i = 1, ncells
            ait_sol_bc( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =               &
              real( tot_tracer( i, 1, :, n ), r_def )
            ait_sol_bc( map_wth(1,i) + 0 ) = ait_sol_bc( map_wth(1,i) + 1 )
          end do
        case(fldname_ait_sol_om)
          do i = 1, ncells
            ait_sol_om( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =               &
              real( tot_tracer( i, 1, :, n ), r_def )
            ait_sol_om( map_wth(1,i) + 0 ) = ait_sol_om( map_wth(1,i) + 1 )
          end do
        case(fldname_n_acc_sol)
          do i = 1, ncells
            n_acc_sol( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                &
              real( tot_tracer( i, 1, :, n ), r_def )
            n_acc_sol( map_wth(1,i) + 0 ) = n_acc_sol( map_wth(1,i) + 1 )
          end do
        case(fldname_acc_sol_su)
          do i = 1, ncells
            acc_sol_su( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =               &
              real( tot_tracer( i, 1, :, n ), r_def )
            acc_sol_su( map_wth(1,i) + 0 ) = acc_sol_su( map_wth(1,i) + 1 )
          end do
        case(fldname_acc_sol_bc)
          do i = 1, ncells
            acc_sol_bc( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =               &
              real( tot_tracer( i, 1, :, n ), r_def )
            acc_sol_bc( map_wth(1,i) + 0 ) = acc_sol_bc( map_wth(1,i) + 1 )
          end do
        case(fldname_acc_sol_om)
          do i = 1, ncells
            acc_sol_om( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =               &
              real( tot_tracer( i, 1, :, n ), r_def )
            acc_sol_om( map_wth(1,i) + 0 ) = acc_sol_om( map_wth(1,i) + 1 )
          end do
        case(fldname_acc_sol_ss)
          do i = 1, ncells
            acc_sol_ss( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =               &
              real( tot_tracer( i, 1, :, n ), r_def )
            acc_sol_ss( map_wth(1,i) + 0 ) = acc_sol_ss( map_wth(1,i) + 1 )
          end do
        case(fldname_acc_sol_du)
          ! No field to update (always zero)
        case(fldname_n_cor_sol)
          do i = 1, ncells
            n_cor_sol( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                &
              real( tot_tracer( i, 1, :, n ), r_def )
            n_cor_sol( map_wth(1,i) + 0 ) = n_cor_sol( map_wth(1,i) + 1 )
          end do
        case(fldname_cor_sol_su)
          do i = 1, ncells
            cor_sol_su( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =               &
              real( tot_tracer( i, 1, :, n ), r_def )
            cor_sol_su( map_wth(1,i) + 0 ) = cor_sol_su( map_wth(1,i) + 1 )
          end do
        case(fldname_cor_sol_bc)
          do i = 1, ncells
            cor_sol_bc( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =               &
              real( tot_tracer( i, 1, :, n ), r_def )
            cor_sol_bc( map_wth(1,i) + 0 ) = cor_sol_bc( map_wth(1,i) + 1 )
          end do
        case(fldname_cor_sol_om)
          do i = 1, ncells
            cor_sol_om( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =               &
              real( tot_tracer( i, 1, :, n ), r_def )
            cor_sol_om( map_wth(1,i) + 0 ) = cor_sol_om( map_wth(1,i) + 1 )
          end do
        case(fldname_cor_sol_ss)
          do i = 1, ncells
            cor_sol_ss( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =               &
              real( tot_tracer( i, 1, :, n ), r_def )
            cor_sol_ss( map_wth(1,i) + 0 ) = cor_sol_ss( map_wth(1,i) + 1 )
          end do
        case(fldname_cor_sol_du)
          ! No field to update (always zero)
        case(fldname_n_ait_ins)
          do i = 1, ncells
            n_ait_ins( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                &
              real( tot_tracer( i, 1, :, n ), r_def )
            n_ait_ins( map_wth(1,i) + 0 ) = n_ait_ins( map_wth(1,i) + 1 )
          end do
        case(fldname_ait_ins_bc)
          do i = 1, ncells
            ait_ins_bc( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =               &
              real( tot_tracer( i, 1, :, n ), r_def )
            ait_ins_bc( map_wth(1,i) + 0 ) = ait_ins_bc( map_wth(1,i) + 1 )
          end do
        case(fldname_ait_ins_om)
          do i = 1, ncells
            ait_ins_om( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =               &
              real( tot_tracer( i, 1, :, n ), r_def )
            ait_ins_om( map_wth(1,i) + 0 ) = ait_ins_om( map_wth(1,i) + 1 )
          end do
        case(fldname_n_acc_ins)
          do i = 1, ncells
            n_acc_ins( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                &
              real( tot_tracer( i, 1, :, n ), r_def )
            n_acc_ins( map_wth(1,i) + 0 ) = n_acc_ins( map_wth(1,i) + 1 )
          end do
        case(fldname_acc_ins_du)
          do i = 1, ncells
            acc_ins_du( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =               &
              real( tot_tracer( i, 1, :, n ), r_def )
            acc_ins_du( map_wth(1,i) + 0 ) = acc_ins_du( map_wth(1,i) + 1 )
          end do
        case(fldname_n_cor_ins)
          do i = 1, ncells
            n_cor_ins( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =                &
              real( tot_tracer( i, 1, :, n ), r_def )
            n_cor_ins( map_wth(1,i) + 0 ) = n_cor_ins( map_wth(1,i) + 1 )
          end do
        case(fldname_cor_ins_du)
          do i = 1, ncells
            cor_ins_du( map_wth(1,i) + 1 : map_wth(1,i) + ntra_lev ) =               &
              real( tot_tracer( i, 1, :, n ), r_def )
            cor_ins_du( map_wth(1,i) + 0 ) = cor_ins_du( map_wth(1,i) + 1 )
          end do
        end select
      end do
    end if  ! outer == outer_iterations .AND. l_tracer
    deallocate(tot_tracer)

  end subroutine conv_gr_code

end module conv_gr_kernel_mod
