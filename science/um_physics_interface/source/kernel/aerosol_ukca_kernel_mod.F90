!-------------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Interface to UKCA GLOMAP-mode aerosol scheme

module aerosol_ukca_kernel_mod

use argument_mod,      only: arg_type, GH_FIELD, GH_SCALAR, GH_INTEGER,        &
                             GH_REAL, GH_READWRITE, GH_READ,                   &
                             ANY_DISCONTINUOUS_SPACE_1,                        &
                             ANY_DISCONTINUOUS_SPACE_2,                        &
                             ANY_DISCONTINUOUS_SPACE_3,                        &
                             ANY_DISCONTINUOUS_SPACE_4,                        &
                             ANY_DISCONTINUOUS_SPACE_5,                        &
                             ANY_DISCONTINUOUS_SPACE_6,                        &
                             CELL_COLUMN

use fs_continuity_mod, only: WTHETA, W3
use kernel_mod,        only: kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: aerosol_ukca_kernel_type
  private
  type(arg_type) :: meta_args(247) = (/            &
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! o3p
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! o1d
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! o3
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! n
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! no
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! no3
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! lumped_n
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! no2
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! n2o5
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! ho2no2
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! hono2
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! h2o2
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! ch4
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! co
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! hcho
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! meoo
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! meooh
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! h
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! oh
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! ho2
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! cl
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! cl2o2
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! clo
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! oclo
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! br
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! lumped_br
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! bro
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! brcl
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! brono2
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! n2o
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! lumped_cl
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! hcl
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! hocl
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! hbr
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! hobr
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! clono2
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! cfcl3
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! cf2cl2
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! mebr
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! hono
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! c2h6
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! etoo
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! etooh
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! mecho
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! meco3
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pan
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! c3h8
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! n_proo
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! i_proo
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! n_prooh
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! i_prooh
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! etcho
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! etco3
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! me2co
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! mecoch2oo
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! mecoch2ooh
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! ppan
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! meono2
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! c5h8
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! iso2
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! isooh
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! ison
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! macr
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! macro2
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! macrooh
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! mpan
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! hacet
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! mgly
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! nald
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! hcooh
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! meco3h
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! meco2h
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! h2
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! meoh
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! msa
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! nh3
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! cs2
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! csul
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! h2s
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! so3
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! passive_o3
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! age_of_air
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! dms
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! so2
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! h2so4
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! dmso
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! monoterpene
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! secondary_organic
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! n_nuc_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! nuc_sol_su
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! nuc_sol_om
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! n_ait_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! ait_sol_su
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! ait_sol_bc
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! ait_sol_om
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! n_acc_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! acc_sol_su
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! acc_sol_bc
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! acc_sol_om
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! acc_sol_ss
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! acc_sol_du
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! n_cor_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! cor_sol_su
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! cor_sol_bc
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! cor_sol_om
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! cor_sol_ss
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! cor_sol_du
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! n_ait_ins
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! ait_ins_bc
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! ait_ins_om
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! n_acc_ins
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! acc_ins_du
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! n_cor_ins
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! cor_ins_du
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! cloud_drop_no_conc
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! drydp_ait_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! drydp_acc_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! drydp_cor_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! drydp_ait_ins
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! drydp_acc_ins
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! drydp_cor_ins
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! wetdp_ait_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! wetdp_acc_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! wetdp_cor_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! rhopar_ait_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! rhopar_acc_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! rhopar_cor_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! rhopar_ait_ins
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! rhopar_acc_ins
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! rhopar_cor_ins
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_wat_ait_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_wat_acc_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_wat_cor_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_su_ait_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_bc_ait_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_om_ait_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_su_acc_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_bc_acc_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_om_acc_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_ss_acc_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_du_acc_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_su_cor_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_bc_cor_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_om_cor_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_ss_cor_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_du_cor_sol
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_bc_ait_ins
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_om_ait_ins
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_du_acc_ins
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, WTHETA ), & ! pvol_du_cor_ins
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! timestep_number
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! current_time_year
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! current_time_month
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! current_time_day
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! current_time_hour
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! current_time_minute
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! current_time_second
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! current_time_daynum
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! previous_time_year
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! previous_time_month
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! previous_time_day
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! previous_time_hour
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! previous_time_minute
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! previous_time_second
       arg_type( GH_SCALAR, GH_INTEGER, GH_READ ),          & ! previous_time_daynum
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! theta_wth
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! exner_in_w3
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! exner_in_wth
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! u3_in_wth
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! m_v_n
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! m_cl_n
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! m_ci_n
       arg_type( GH_FIELD, GH_REAL, GH_READ, W3 ),          & ! height_w3
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! height_wth
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! ls_rain_3d
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! ls_snow_3d
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! autoconv
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! accretion
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! rim_cry
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! rim_agg
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! tke_bl
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! conv_rain_3d
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! conv_snow_3d
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! cf_bulk
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! cf_liq
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1 ), & ! tile_fraction
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1 ), & ! tile_temperature
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1 ), & ! tile_heat_flux
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1 ), & ! gc_tile
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_2 ), & ! leaf_area_index
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_2 ), & ! canopy_height
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3 ), & ! soil_moisture
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! latitude
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! longitude
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! sin_stellar_declination_rts
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! stellar_eqn_of_time_rts
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! soil_roughness
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! z0m
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1 ), & ! urbztm
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! ustar
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! wspd10m
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! chloro_sea
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! zh
       arg_type( GH_FIELD, GH_REAL, GH_READ, W3 ),          & ! wetrho_in_w3
       arg_type( GH_FIELD, GH_REAL, GH_READ, W3 ),          & ! rhokh_bl
       arg_type( GH_FIELD, GH_REAL, GH_READ, W3 ),          & ! rdz_tq_bl
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! dtrdz_tq_bl
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! zhsc
       arg_type( GH_FIELD, GH_INTEGER, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! level_ent
       arg_type( GH_FIELD, GH_INTEGER, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! level_ent_dsc
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_5 ), & ! ent_we_lim
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_5 ), & ! ent_t_frac
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_5 ), & ! ent_zrzi
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_5 ), & ! ent_we_lim_dsc
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_5 ), & ! ent_t_frac_dsc
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_5 ), & ! ent_zrzi_dsc
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! h2o2_limit
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! dms_conc_ocean
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_6 ), & ! dust_flux
       arg_type( GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_4 ), & ! surf_wetness
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_c2h6
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_c3h8
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_c5h8
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_ch4
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_co
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_hcho
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_me2co
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_mecho
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_meoh
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_nh3
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_no
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_bc_biofuel
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_bc_biomass_high
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_bc_biomass_low
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_bc_fossil
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_dms_land
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_monoterp
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_om_biofuel
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_om_biomass_high
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_om_biomass_low
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_om_fossil
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_so2_high
       arg_type( GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_4 ), & ! emiss_so2_low
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! emiss_no_aircrft
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! emiss_bc_biomass
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA ),      & ! emiss_om_biomass
       arg_type( GH_FIELD, GH_REAL, GH_READ, WTHETA )       & ! emiss_so2_nat
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: aerosol_ukca_code
end type

public aerosol_ukca_code
contains

!> @brief UKCA GLOMAP-mode aerosol scheme time step.
!> @details Copy the required tracer and non-transported prognostic fields
!>          to the corresponding UM fields, likewise copy or create the
!>          required environmental driver fields, do the UKCA time step
!>          and finally copy the updated UM tracer and non-transported
!>          prognostics back to the LFRic fields.
!>          Note that where the expected emission input fields are indicated
!>          to be 'expressed as C' or 'expressed as S' in the argument list
!>          this assumes that they are registered as such when UKCA is
!>          initialised.
!>          Mass mixing ratios are in kg per kg of air, number mixing ratio
!>          for the aerosol modes are in particles per molecule of air.
!> @param[in]     nlayers             Number of layers
!> @param[in,out] o3p                 Oxygen_ground_state m.m.r
!> @param[in,out] o1d                 Oxygen_excited_state m.m.r
!> @param[in,out] o3                  Ozone m.m.r.
!> @param[in,out] n                   Nitrogen_radical m.m.r
!> @param[in,out] no                  Nitric_oxide m.m.r
!> @param[in,out] no3                 Nitrate m.m.r.
!> @param[in,out] lumped_n            Lumped N species as Nitrogen_dioxide m.m.r
!> @param[in,out] no2                 Nitrogen_dioxide m.m.r
!> @param[in,out] n2o5                Dinitrogen_pentoxide m.m.r
!> @param[in,out] ho2no2              Peroxynitric_acid m.m.r
!> @param[in,out] hono2               Nitric_acid m.m.r
!> @param[in,out] h2o2                Hydrogen peroxide m.m.r.
!> @param[in,out] ch4                 Methane m.m.r
!> @param[in,out] co                  Carbon_monoxide m.m.r
!> @param[in,out] hcho                Formaldehyde m.m.r
!> @param[in,out] meoo                Methyl_peroxy_radical m.m.r
!> @param[in,out] meooh               Methyl_hydroperoxide m.m.r
!> @param[in,out] oh                  Hydroxyl radical m.m.r.
!> @param[in,out] ho2                 Hydroperoxyl radical m.m.r.
!> @param[in,out] h                   Hydrogen_radical m.m.r
!> @param[in,out] cl                  Chloride_radical m.m.r
!> @param[in,out] cl2o2               Chlorine_monoxide_dimer m.m.r
!> @param[in,out] clo                 Chlorine_monoxide m.m.r
!> @param[in,out] oclo                Chlorine_dioxide m.m.r
!> @param[in,out] br                  Bromine m.m.r
!> @param[in,out] lumped_br           Lumped Br species as Bromine_monoxide m.m.r
!> @param[in,out] bro                 Bromine_monoxide m.m.r
!> @param[in,out] brcl                Bromine_chloride m.m.r
!> @param[in,out] brono2              Bromine_nitrate m.m.r
!> @param[in,out] n2o                 Nitrous_oxide m.m.r
!> @param[in,out] lumped_cl           Lumped Cl spc as Hydrogen_chloride m.m.r
!> @param[in,out] hcl                 Hydrogen_chloride m.m.r
!> @param[in,out] hocl                Hydrochlorous_acid m.m.r
!> @param[in,out] hbr                 Hydrogen_bromide m.m.r
!> @param[in,out] hobr                Hydrobromous_acid m.m.r
!> @param[in,out] clono2              Chlorine_nitrate m.m.r
!> @param[in,out] cfcl3               CFC11 m.m.r
!> @param[in,out] cf2cl2              CFC12 m.m.r
!> @param[in,out] mebr                Methyl_bromide m.m.r
!> @param[in,out] hono                Nitrous_acid m.m.r
!> @param[in,out] c2h6                Ethane m.m.r
!> @param[in,out] etoo                Ethyl_peroxy_radical m.m.r
!> @param[in,out] etooh               Ethyl_hydroperoxide m.m.r
!> @param[in,out] mecho               Acetaldehyde m.m.r
!> @param[in,out] meco3               Peroxy_acetic_acide m.m.r
!> @param[in,out] pan                 Peroxy_acetyl_nitrate m.m.r
!> @param[in,out] c3h8                Propane m.m.r
!> @param[in,out] n_proo              N-propyl_peroxy_radical m.m.r
!> @param[in,out] i_proo              I-propyl_peroxy_radical m.m.r
!> @param[in,out] n_prooh             N-propyl_alcohol m.m.r
!> @param[in,out] i_prooh             I-propyl_alcohol m.m.r
!> @param[in,out] etcho               Propionaldehyde m.m.r
!> @param[in,out] etco3               Ethyl_carbonate m.m.r
!> @param[in,out] me2co               Acetone m.m.r
!> @param[in,out] mecoch2oo           Acetone_peroxy_radical m.m.r
!> @param[in,out] mecoch2ooh          Acetone_hydroperoxide m.m.r
!> @param[in,out] ppan                Peroxypropionyl_nitrate m.m.r
!> @param[in,out] meono2              Methyl_nitrate m.m.r
!> @param[in,out] c5h8                Isoprene m.m.r
!> @param[in,out] iso2                Isoprene_peroxy_radical m.m.r
!> @param[in,out] isooh               Isoprene_hydroperoxide m.m.r
!> @param[in,out] ison                Isoprene_peroxy_acetyl_nitrate m.m.r
!> @param[in,out] macr                Methacrolein m.m.r
!> @param[in,out] macro2              Lumped_mvk_macr_hydroperoxy_radical m.m.r
!> @param[in,out] macrooh             Lumped_mvk_macr_hydroperoxide m.m.r
!> @param[in,out] mpan                Methacryolyl_peroxynitrate m.m.r
!> @param[in,out] hacet               Hydroxy_acetone m.m.r
!> @param[in,out] mgly                Methylglyoxal m.m.r
!> @param[in,out] nald                Nitroxy_acetaldehyde m.m.r
!> @param[in,out] hcooh               Formic_acid m.m.r
!> @param[in,out] meco3h              Peroxy_acetic_acid m.m.r
!> @param[in,out] meco2h              Acetic_acid m.m.r
!> @param[in,out] h2                  Hydrogen m.m.r
!> @param[in,out] meoh                Methanol m.m.r
!> @param[in,out] msa                 Methyl_sulphonic_acid m.m.r
!> @param[in,out] nh3                 Ammonia m.m.r
!> @param[in,out] cs2                 Carbon_disulphide m.m.r
!> @param[in,out] csul                Carbonyl_sulphide m.m.r
!> @param[in,out] h2s                 Hydrogen_sulphide m.m.r
!> @param[in,out] so3                 Sulphur_trioxide m.m.r
!> @param[in,out] passive_o3          Passive ozone tracer m.m.r
!> @param[in,out] age_of_air          Age of air tracer (s)
!> @param[in,out] dms                 Dimethyl sulfide m.m.r.
!> @param[in,out] so2                 Sulfur dioxide m.m.r.
!> @param[in,out] h2so4               Sulfuric acid m.m.r.
!> @param[in,out] dmso                Dimethyl sulfoxide m.m.r.
!> @param[in,out] monoterpene         Monoterpene m.m.r.
!> @param[in,out] secondary organic   Secondary organic m.m.r.
!> @param[in,out] n_nuc_sol           Aerosol field: n.m.r. of soluble nucleation mode
!> @param[in,out] nuc_sol_su          Aerosol field: m.m.r. of H2SO4 in soluble nucleation mode
!> @param[in,out] nuc_sol_om          Aerosol field: m.m.r. of organic matter in soluble nucleation mode
!> @param[in,out] n_ait_sol           Aerosol field: n.m.r. of soluble Aitken mode
!> @param[in,out] ait_sol_su          Aerosol field: m.m.r. of H2SO4 in soluble Aitken mode
!> @param[in,out] ait_sol_bc          Aerosol field: m.m.r. of black carbon in soluble Aitken mode
!> @param[in,out] ait_sol_om          Aerosol field: m.m.r. of organic matter in soluble Aitken mode
!> @param[in,out] n_acc_sol           Aerosol field: n.m.r. of soluble accumulation mode
!> @param[in,out] acc_sol_su          Aerosol field: m.m.r. of H2SO4 in soluble accumulation mode
!> @param[in,out] acc_sol_bc          Aerosol field: m.m.r. of black carbon in soluble accumulation mode
!> @param[in,out] acc_sol_om          Aerosol field: m.m.r. of organic matter in soluble accumulation mode
!> @param[in,out] acc_sol_ss          Aerosol field: m.m.r. of sea salt in soluble accumulation mode
!> @param[in,out] acc_sol_du          Aerosol field: m.m.r. of dust in soluble accumulation mode
!> @param[in,out] n_cor_sol           Aerosol field: n.m.r. of soluble coarse mode
!> @param[in,out] cor_sol_su          Aerosol field: m.m.r. of H2SO4 in soluble coarse mode
!> @param[in,out] cor_sol_bc          Aerosol field: m.m.r. of black carbon in soluble coarse mode
!> @param[in,out] cor_sol_om          Aerosol field: m.m.r. of organic matter in soluble coarse mode
!> @param[in,out] cor_sol_ss          Aerosol field: m.m.r. of sea salt in soluble coarse mode
!> @param[in,out] cor_sol_du          Aerosol field: m.m.r. of dust in soluble coarse mode
!> @param[in,out] n_ait_ins           Aerosol field: n.m.r. of insoluble Aitken mode
!> @param[in,out] ait_ins_bc          Aerosol field: m.m.r. of black carbon in insoluble Aitken mode
!> @param[in,out] ait_ins_om          Aerosol field: m.m.r. of organic matter in insoluble Aitken mode
!> @param[in,out] n_acc_ins           Aerosol field: n.m.r. of insoluble accumulation mode
!> @param[in,out] acc_ins_du          Aerosol field: m.m.r. of dust in insoluble accumulation mode
!> @param[in,out] n_cor_ins           Aerosol field: n.m.r. of insoluble coarse mode
!> @param[in,out] cor_ins_du          Aerosol field: m.m.r. of dust in insoluble coarse mode
!> @param[in,out] cloud_drop_no_conc  Cloud Droplet Number Concentration
!> @param[in,out] drydp_ait_sol       Dry modal diam for aitken soluble mode (m)
!> @param[in,out] drydp_acc_sol       Dry modal diam for accum soluble mode (m)
!> @param[in,out] drydp_cor_sol       Dry modal diam for coarse soluble mode (m)
!> @param[in,out] drydp_ait_ins       Dry modal diam for aitken insoluble mode (m)
!> @param[in,out] drydp_acc_ins       Dry modal diam for accum insoluble mode (m)
!> @param[in,out] drydp_cor_ins       Dry modal diam for coarse insoluble mode (m)
!> @param[in,out] wetdp_ait_sol       Wet modal diam for aitken soluble mode (m)
!> @param[in,out] wetdp_acc_sol       Wet modal diam for accum soluble mode (m)
!> @param[in,out] wetdp_cor_sol       Wet modal diam for coarse soluble mode (m)
!> @param[in,out] rhopar_ait_sol      Particle density - aitken soluble mode (kg m-3)
!> @param[in,out] rhopar_acc_sol      Particle density - accum soluble mode (kg m-3)
!> @param[in,out] rhopar_cor_sol      Particle density - coarse soluble mode (kg m-3)
!> @param[in,out] rhopar_ait_ins      Particle density - aitken insoluble mode (kg m-3)
!> @param[in,out] rhopar_acc_ins      Particle density - accum insoluble mode (kg m-3)
!> @param[in,out] rhopar_cor_ins      Particle density - coarse soluble mode (kg m-3)
!> @param[in,out] pvol_wat_ait_sol    Partial vol. of water in aitken soluble mode (m3)
!> @param[in,out] pvol_wat_acc_sol    Partial vol. of water in accum soluble mode (m3)
!> @param[in,out] pvol_wat_cor_sol    Partial vol. of water in coarse soluble mode (m3)
!> @param[in,out] pvol_su_ait_sol     Partial vol. of h2so4 in aitken soluble mode (m3)
!> @param[in,out] pvol_bc_ait_sol     Partial vol. of black carbon in aitken soluble mode (m3)
!> @param[in,out] pvol_om_ait_sol     Partial vol. of organic matter in aitken soluble mode (m3)
!> @param[in,out] pvol_su_acc_sol     Partial vol. of h2so4 in accum soluble mode (m3)
!> @param[in,out] pvol_bc_acc_sol     Partial vol. of black carbon in accum soluble mode (m3)
!> @param[in,out] pvol_om_acc_sol     Partial vol. of organic matter in accum soluble mode (m3)
!> @param[in,out] pvol_ss_acc_sol     Partial vol. of seasalt in accum soluble mode (m3)
!> @param[in,out] pvol_du_acc_sol     Partial vol. of dust in accum soluble mode (m3)
!> @param[in,out] pvol_su_cor_sol     Partial vol. of h2so4 in coarse soluble mode (m3)
!> @param[in,out] pvol_bc_cor_sol     Partial vol. of black carbon in coarse soluble mode (m3)
!> @param[in,out] pvol_om_cor_sol     Partial vol. of organic matter in coarse soluble mode (m3)
!> @param[in,out] pvol_ss_cor_sol     Partial vol. of seasalt in coarse soluble mode (m3)
!> @param[in,out] pvol_du_cor_sol     Partial vol. of dust in coarse soluble mode (m3)
!> @param[in,out] pvol_bc_ait_ins     Partial vol. of black carbon in aitken insoluble mode (m3)
!> @param[in,out] pvol_om_ait_ins     Partial vol. of organic matter in aitken insoluble mode (m3)
!> @param[in,out] pvol_du_acc_ins     Partial vol. of dust in accum insoluble mode (m3)
!> @param[in,out] pvol_du_cor_ins     Partial vol. of dust in coarse insoluble mode (m3)
!> @param[in]     timestep_number     Time step number
!> @param[in]     current_time_year   Current model year
!> @param[in]     current_time_month  Current model month
!> @param[in]     current_time_day    Current model day
!> @param[in]     current_time_hour   Current model hour
!> @param[in]     current_time_minute Current model minute
!> @param[in]     current_time_second Current model second
!> @param[in]     current_time_daynum Current model day number
!> @param[in]     previous_time_year  Model year at previous time step
!> @param[in]     previous_time_month Model month at previous time step
!> @param[in]     previous_time_day   Model day of month at previous time step
!> @param[in]     previous_time_hour  Model hour at previous time step
!> @param[in]     previous_time_minute Model minute at previous time step
!> @param[in]     previous_time_second Model second at previous time step
!> @param[in]     previous_time_daynum Model day number at previous time step
!> @param[in]     theta_wth           Potential temperature field (K)
!> @param[in]     exner_in_w3         Exner pressure in w3 space
!> @param[in]     exner_in_wth        Exner pressure in theta space
!> @param[in]     u3_in_wth           'Vertical' wind in theta space
!> @param[in]     m_v_n               Vapour mixing ratio at time level n
!> @param[in]     m_cl_n              Cloud liq mixing ratio at time level n
!> @param[in]     m_ci_n              Cloud ice mixing ratio at time level n
!> @param[in]     height_w3           Height of density space above surface (m)
!> @param[in]     height_wth          Height of theta space above surface (m)
!> @param[in]     ls_rain_3d          Large scale rainfall flux (kg m-2 s-1)
!> @param[in]     ls_snow_3d          Large scale snowfall flux (kg m-2 s-1)
!> @param[in]     autoconv            Rain autoconversion rate (kg kg-1 s-1)
!> @param[in]     accretion           Rain accretion rate (kg kg-1 s-1)
!> @param[in]     rim_cry             Riming rate for ice crystals (kg kg-1 s-1)
!> @param[in]     rim_agg             Riming rate for ice aggregates (kg kg-1 s-1)
!> @param[in]     tke_bl              Turbulent_kinetic_energy (m2 s-2)
!> @param[in]     conv_rain_3d        Convective rainfall flux (kg m-2 s-1)
!> @param[in]     conv_snow_3d        Convective snowfall flux (kg m-2 s-1)
!> @param[in]     cf_bulk             Bulk cloud fraction
!> @param[in]     cf_liq              Liquid cloud fraction
!> @param[in]     tile_fraction       Surface tile fractions
!> @param[in]     tile_temperature    Surface tile temperatures (K)
!> @param[in]     tile_heat_flux      Surface tile heat fluxes (W m-2)
!> @param[in]     gc_tile             Stomatal conductance (m s-1)
!> @param[in]     leaf_area_index     Leaf Area Index
!> @param[in]     canopy_height       Canopy height (m)
!> @param[in]     soil_moisture       Soil moisture content (kg m-2)
!> @param[in]     latitude            Latitude field (radians)
!> @param[in]     longitude           Longitude field (radians)
!> @param[in]     sin_stellar_declination_rts Sin of stellar declination
!> @param[in]     stellar_eqn_of_time Stellar equation of time (radians)
!> @param[in]     soil_roughness      Bare soil surface roughness length (m)
!> @param[in]     z0m                 Cell surface roughness length (m)
!> @param[in]     urbztm              Urban effective roughness length
!> @param[in]     ustar               Friction velocity (m s-1)
!> @param[in]     wspd10m             Wind speed at 10 m (m s-1)
!> @param[in]     chloro_sea          Sea surface chlorophyll content
!> @param[in]     zh                  Boundary layer depth (m)
!> @param[in]     wetrho_in_w3        Wet density (kg m-3)
!> @param[in]     rhokh_bl            Scalar eddy diffusivity * rho (kg m-1 s-1)
!> @param[in]     rdz_tq_bl           1/dz (m-1)
!> @param[in]     dtrdz_tq_bl         dt/(rho*r*r*dz) (s kg-1)
!> @param[in]     zhsc                Height at top of decoupled stratocumulus layer (m)
!> @param[in]     level_ent           Level of surface mixed layer inversion
!> @param[in]     level_ent_dsc       Level of decoupled stratocumulus inversion
!> @param[in]     ent_we_lim          Rho * entrainment rate at surface ML inversion (kg m-2 s-1)
!> @param[in]     ent_t_frac          Fraction of time surface ML inversion is above level
!> @param[in]     ent_zrzi            Level height as fraction of SML inversion height above ML base
!> @param[in]     ent_we_lim_dsc      Rho * entrainment rate at DSC inversion (kg m-2 s-1)
!> @param[in]     ent_t_frac_dsc      Fraction of time DSC inversion is above level
!> @param[in]     ent_zrzi_dsc        Level height as fraction of DSC inversion height above DSC ML base
!> @param[in]     h2o2_limit          Hydrogen peroxide m.m.r. upper limit
!> @param[in]     dms_conc_ocean      DMS concentration in seawater (nmol l-1)
!> @param[in]     dust_flux           Dust emission fluxes in CLASSIC size divisions (kg m-2 s-1)
!> @param[in,out] surf_wetness        Surface wetness prognostic (dimensionless)
!> @param[in]     emiss_c2h6          C2H6 emissions (kg m-2 s-1)
!> @param[in]     emiss_c3h8          C3H8 emissions (kg m-2 s-1)
!> @param[in]     emiss_c5h8          Biogenic C5H8 emissions (kg m-2 s-1)
!> @param[in]     emiss_ch4           CH4 emissions (kg m-2 s-1)
!> @param[in]     emiss_co            CO emissions (kg m-2 s-1)
!> @param[in]     emiss_hcho          HCHO emissions (kg m-2 s-1)
!> @param[in]     emiss_me2co         Me2CO emissions (kg m-2 s-1)
!> @param[in]     emiss_mecho         MeCHO emissions (kg m-2 s-1)
!> @param[in]     emiss_meoh          Biogenic methanol (CH3OH) emissions (kg m-2 s-1)
!> @param[in]     emiss_nh3           Ammonia gas emissions (kg m-2 s-1)
!> @param[in]     emiss_no            NOx emissions (kg m-2 s-1)
!> @param[in]     emiss_bc_biofuel    Black C biofuel emissions (kg m-2 s-1)
!> @param[in]     emiss_bc_biomass_high Black C emissions from biomass burning (kg m-2 s-1)
!> @param[in]     emiss_bc_biomass_low  Black C emissions from biomass burning (kg m-2 s-1)
!> @param[in]     emiss_bc_fossil     Black C fossil fuel emissions (kg m-2 s-1)
!> @param[in]     emiss_dms_land      DMS emissions from land surface (kg m-2 s-1)
!> @param[in]     emiss_monoterp      Monoterpene emissions expressed as C
!> @param[in]     emiss_om_biofuel    Organic matter biofuel emissions expressed as C (kg m-2 s-1)
!> @param[in]     emiss_om_biomass_high Organic matter emissions from biomass burning expressed as C (kg m-2 s-1)
!> @param[in]     emiss_om_biomass_low  Organic matter emissions from biomass burning expressed as C (kg m-2 s-1)
!> @param[in]     emiss_om_fossil     Organic matter fossil fuel emissions expressed as C (kg m-2 s-1)
!> @param[in]     emiss_so2_high      High-level SO2 emissions expressed as S (kg m-2 s-1)
!> @param[in]     emiss_so2_low       Low-level SO2 emissions expressed as S (kg m-2 s-1)
!> @param[in]     emiss_no_aircrft    NOx aircraft emissions (kg m-2 s-1)
!> @param[in]     emiss_bc_biomass    Black C emissions from biomass burning (kg m-2 s-1)
!> @param[in]     emiss_om_biomass    Organic matter emissions from biomass burning expressed as C (kg m-2 s-1)
!> @param[in]     emiss_so2_nat       SO2 natural emissions expressed as S (kg m-2 s-1)
!> @param[in]     ndf_wth             Number of DOFs per cell for potential temperature space
!> @param[in]     undf_wth            Number of unique DOFs for potential temperature space
!> @param[in]     map_wth             Dofmap for the cell at the base of the column for potential temperature space
!> @param[in]     ndf_w3              Number of DOFs per cell for density space
!> @param[in]     undf_w3             Number of unique DOFs for density space
!> @param[in]     map_w3              Dofmap for the cell at the base of the column for density space
!> @param[in]     ndf_tile            Number of DOFs per cell for surface tiles
!> @param[in]     undf_tile           Number of unique DOFs for surface tiles
!> @param[in]     map_tile            Dofmap for cell for surface tiles
!> @param[in]     ndf_pft             Number of DOFs per cell for PFTs
!> @param[in]     undf_pft            Number of unique DOFs for PFTs
!> @param[in]     map_pft             Dofmap for cell for PFTs
!> @param[in]     ndf_soil            Number of DOFs per cell for soil levels
!> @param[in]     undf_soil           Number of unique DOFs for soil levels
!> @param[in]     map_soil            Dofmap for cell for soil levels
!> @param[in]     ndf_2d              Number of DOFs per cell for 2D fields
!> @param[in]     undf_2d             Number of unique DOFs for 2D fields
!> @param[in]     map_2d              Dofmap for cell for 2D fields
!> @param[in]     ndf_ent             Number of DOFs per cell for entrainment levels
!> @param[in]     undf_ent            Number of unique DOFs for entrainment levels
!> @param[in]     map_ent             Dofmap for cell for entrainment levels
!> @param[in]     ndf_dust            Number of DOFs per cell for dust divisions
!> @param[in]     undf_dust           Number of unique DOFs for dust divisions
!> @param[in]     map_dust            Dofmap for cell for dust divisions

subroutine aerosol_ukca_code( nlayers,                                         &
                              o3p,                                             &
                              o1d,                                             &
                              o3,                                              &
                              n,                                               &
                              no,                                              &
                              no3,                                             &
                              lumped_n,                                        &
                              no2,                                             &
                              n2o5,                                            &
                              ho2no2,                                          &
                              hono2,                                           &
                              h2o2,                                            &
                              ch4,                                             &
                              co,                                              &
                              hcho,                                            &
                              meoo,                                            &
                              meooh,                                           &
                              h,                                               &
                              oh,                                              &
                              ho2,                                             &
                              cl,                                              &
                              cl2o2,                                           &
                              clo,                                             &
                              oclo,                                            &
                              br,                                              &
                              lumped_br,                                       &
                              bro,                                             &
                              brcl,                                            &
                              brono2,                                          &
                              n2o,                                             &
                              lumped_cl,                                       &
                              hcl,                                             &
                              hocl,                                            &
                              hbr,                                             &
                              hobr,                                            &
                              clono2,                                          &
                              cfcl3,                                           &
                              cf2cl2,                                          &
                              mebr,                                            &
                              hono,                                            &
                              c2h6,                                            &
                              etoo,                                            &
                              etooh,                                           &
                              mecho,                                           &
                              meco3,                                           &
                              pan,                                             &
                              c3h8,                                            &
                              n_proo,                                          &
                              i_proo,                                          &
                              n_prooh,                                         &
                              i_prooh,                                         &
                              etcho,                                           &
                              etco3,                                           &
                              me2co,                                           &
                              mecoch2oo,                                       &
                              mecoch2ooh,                                      &
                              ppan,                                            &
                              meono2,                                          &
                              c5h8,                                            &
                              iso2,                                            &
                              isooh,                                           &
                              ison,                                            &
                              macr,                                            &
                              macro2,                                          &
                              macrooh,                                         &
                              mpan,                                            &
                              hacet,                                           &
                              mgly,                                            &
                              nald,                                            &
                              hcooh,                                           &
                              meco3h,                                          &
                              meco2h,                                          &
                              h2,                                              &
                              meoh,                                            &
                              msa,                                             &
                              nh3,                                             &
                              cs2,                                             &
                              csul,                                            &
                              h2s,                                             &
                              so3,                                             &
                              passive_o3,                                      &
                              age_of_air,                                      &
                              dms,                                             &
                              so2,                                             &
                              h2so4,                                           &
                              dmso,                                            &
                              monoterpene,                                     &
                              secondary_organic,                               &
                              n_nuc_sol,                                       &
                              nuc_sol_su,                                      &
                              nuc_sol_om,                                      &
                              n_ait_sol,                                       &
                              ait_sol_su,                                      &
                              ait_sol_bc,                                      &
                              ait_sol_om,                                      &
                              n_acc_sol,                                       &
                              acc_sol_su,                                      &
                              acc_sol_bc,                                      &
                              acc_sol_om,                                      &
                              acc_sol_ss,                                      &
                              acc_sol_du,                                      &
                              n_cor_sol,                                       &
                              cor_sol_su,                                      &
                              cor_sol_bc,                                      &
                              cor_sol_om,                                      &
                              cor_sol_ss,                                      &
                              cor_sol_du,                                      &
                              n_ait_ins,                                       &
                              ait_ins_bc,                                      &
                              ait_ins_om,                                      &
                              n_acc_ins,                                       &
                              acc_ins_du,                                      &
                              n_cor_ins,                                       &
                              cor_ins_du,                                      &
                              cloud_drop_no_conc,                              &
                              drydp_ait_sol,                                   &
                              drydp_acc_sol,                                   &
                              drydp_cor_sol,                                   &
                              drydp_ait_ins,                                   &
                              drydp_acc_ins,                                   &
                              drydp_cor_ins,                                   &
                              wetdp_ait_sol,                                   &
                              wetdp_acc_sol,                                   &
                              wetdp_cor_sol,                                   &
                              rhopar_ait_sol,                                  &
                              rhopar_acc_sol,                                  &
                              rhopar_cor_sol,                                  &
                              rhopar_ait_ins,                                  &
                              rhopar_acc_ins,                                  &
                              rhopar_cor_ins,                                  &
                              pvol_wat_ait_sol,                                &
                              pvol_wat_acc_sol,                                &
                              pvol_wat_cor_sol,                                &
                              pvol_su_ait_sol,                                 &
                              pvol_bc_ait_sol,                                 &
                              pvol_om_ait_sol,                                 &
                              pvol_su_acc_sol,                                 &
                              pvol_bc_acc_sol,                                 &
                              pvol_om_acc_sol,                                 &
                              pvol_ss_acc_sol,                                 &
                              pvol_du_acc_sol,                                 &
                              pvol_su_cor_sol,                                 &
                              pvol_bc_cor_sol,                                 &
                              pvol_om_cor_sol,                                 &
                              pvol_ss_cor_sol,                                 &
                              pvol_du_cor_sol,                                 &
                              pvol_bc_ait_ins,                                 &
                              pvol_om_ait_ins,                                 &
                              pvol_du_acc_ins,                                 &
                              pvol_du_cor_ins,                                 &
                              timestep_number,                                 &
                              current_time_year,                               &
                              current_time_month,                              &
                              current_time_day,                                &
                              current_time_hour,                               &
                              current_time_minute,                             &
                              current_time_second,                             &
                              current_time_daynum,                             &
                              previous_time_year,                              &
                              previous_time_month,                             &
                              previous_time_day,                               &
                              previous_time_hour,                              &
                              previous_time_minute,                            &
                              previous_time_second,                            &
                              previous_time_daynum,                            &
                              theta_wth,                                       &
                              exner_in_w3,                                     &
                              exner_in_wth,                                    &
                              u3_in_wth,                                       &
                              m_v_n,                                           &
                              m_cl_n,                                          &
                              m_ci_n,                                          &
                              height_w3,                                       &
                              height_wth,                                      &
                              ls_rain_3d,                                      &
                              ls_snow_3d,                                      &
                              autoconv,                                        &
                              accretion,                                       &
                              rim_cry,                                         &
                              rim_agg,                                         &
                              tke_bl,                                          &
                              conv_rain_3d,                                    &
                              conv_snow_3d,                                    &
                              cf_bulk,                                         &
                              cf_liq,                                          &
                              tile_fraction,                                   &
                              tile_temperature,                                &
                              tile_heat_flux,                                  &
                              gc_tile,                                         &
                              leaf_area_index,                                 &
                              canopy_height,                                   &
                              soil_moisture,                                   &
                              latitude,                                        &
                              longitude,                                       &
                              sin_stellar_declination_rts,                     &
                              stellar_eqn_of_time_rts,                         &
                              soil_roughness,                                  &
                              z0m,                                             &
                              urbztm,                                          &
                              ustar,                                           &
                              wspd10m,                                         &
                              chloro_sea,                                      &
                              zh,                                              &
                              wetrho_in_w3,                                    &
                              rhokh_bl,                                        &
                              rdz_tq_bl,                                       &
                              dtrdz_tq_bl,                                     &
                              zhsc,                                            &
                              level_ent,                                       &
                              level_ent_dsc,                                   &
                              ent_we_lim,                                      &
                              ent_t_frac,                                      &
                              ent_zrzi,                                        &
                              ent_we_lim_dsc,                                  &
                              ent_t_frac_dsc,                                  &
                              ent_zrzi_dsc,                                    &
                              h2o2_limit,                                      &
                              dms_conc_ocean,                                  &
                              dust_flux,                                       &
                              surf_wetness,                                    &
                              emiss_c2h6,                                      &
                              emiss_c3h8,                                      &
                              emiss_c5h8,                                      &
                              emiss_ch4,                                       &
                              emiss_co,                                        &
                              emiss_hcho,                                      &
                              emiss_me2co,                                     &
                              emiss_mecho,                                     &
                              emiss_meoh,                                      &
                              emiss_nh3,                                       &
                              emiss_no,                                        &
                              emiss_bc_biofuel,                                &
                              emiss_bc_biomass_high,                           &
                              emiss_bc_biomass_low,                            &
                              emiss_bc_fossil,                                 &
                              emiss_dms_land,                                  &
                              emiss_monoterp,                                  &
                              emiss_om_biofuel,                                &
                              emiss_om_biomass_high,                           &
                              emiss_om_biomass_low,                            &
                              emiss_om_fossil,                                 &
                              emiss_so2_high,                                  &
                              emiss_so2_low,                                   &
                              emiss_no_aircrft,                                &
                              emiss_bc_biomass,                                &
                              emiss_om_biomass,                                &
                              emiss_so2_nat,                                   &
                              ndf_wth, undf_wth, map_wth,                      &
                              ndf_w3, undf_w3, map_w3,                         &
                              ndf_tile, undf_tile, map_tile,                   &
                              ndf_pft, undf_pft, map_pft,                      &
                              ndf_soil, undf_soil, map_soil,                   &
                              ndf_2d, undf_2d, map_2d,                         &
                              ndf_ent, undf_ent, map_ent,                      &
                              ndf_dust, undf_dust, map_dust)

  use constants_mod,    only: r_def, i_def, r_um, i_um, i_timestep,            &
                              radians_to_degrees
  use conversions_mod,  only: rsec_per_day, rsec_per_hour
  use jules_control_init_mod, only: n_surf_tile, n_land_tile, n_sea_ice_tile,  &
                                    first_sea_ice_tile
  use um_ukca_init_mod, only: tracer_names, ntp_names,                         &
                              env_names_scalar_real,                           &
                              env_names_flat_integer,                          &
                              env_names_flat_real,                             &
                              env_names_flat_logical,                          &
                              env_names_flatpft_real,                          &
                              env_names_fullht_real,                           &
                              env_names_fullht0_real,                          &
                              env_names_fullhtp1_real,                         &
                              env_names_bllev_real,                            &
                              env_names_entlev_real,                           &
                              env_names_land_real,                             &
                              env_names_landtile_real,                         &
                              env_names_landtile_logical,                      &
                              env_names_landpft_real,                          &
                              emiss_names_flat,                                &
                              emiss_names_fullht,                              &
                              fldname_o3p,                                     &
                              fldname_o3,                                      &
                              fldname_n,                                       &
                              fldname_no,                                      &
                              fldname_no3,                                     &
                              fldname_lumped_n,                                &
                              fldname_n2o5,                                    &
                              fldname_ho2no2,                                  &
                              fldname_hono2,                                   &
                              fldname_ch4,                                     &
                              fldname_co,                                      &
                              fldname_hcho,                                    &
                              fldname_meooh,                                   &
                              fldname_h,                                       &
                              fldname_ch2o,                                    &
                              fldname_oh,                                      &
                              fldname_ho2,                                     &
                              fldname_cl,                                      &
                              fldname_cl2o2,                                   &
                              fldname_clo,                                     &
                              fldname_oclo,                                    &
                              fldname_br,                                      &
                              fldname_lumped_br,                               &
                              fldname_brcl,                                    &
                              fldname_brono2,                                  &
                              fldname_n2o,                                     &
                              fldname_lumped_cl,                               &
                              fldname_hocl,                                    &
                              fldname_hbr,                                     &
                              fldname_hobr,                                    &
                              fldname_clono2,                                  &
                              fldname_cfcl3,                                   &
                              fldname_cf2cl2,                                  &
                              fldname_mebr,                                    &
                              fldname_hono,                                    &
                              fldname_c2h6,                                    &
                              fldname_etooh,                                   &
                              fldname_mecho,                                   &
                              fldname_pan,                                     &
                              fldname_c3h8,                                    &
                              fldname_n_prooh,                                 &
                              fldname_i_prooh,                                 &
                              fldname_etcho,                                   &
                              fldname_me2co,                                   &
                              fldname_mecoch2ooh,                              &
                              fldname_ppan,                                    &
                              fldname_meono2,                                  &
                              fldname_c5h8,                                    &
                              fldname_isooh,                                   &
                              fldname_ison,                                    &
                              fldname_macr,                                    &
                              fldname_macrooh,                                 &
                              fldname_mpan,                                    &
                              fldname_hacet,                                   &
                              fldname_mgly,                                    &
                              fldname_nald,                                    &
                              fldname_hcooh,                                   &
                              fldname_meco3h,                                  &
                              fldname_meco2h,                                  &
                              fldname_h2,                                      &
                              fldname_meoh,                                    &
                              fldname_msa,                                     &
                              fldname_nh3,                                     &
                              fldname_cs2,                                     &
                              fldname_csul,                                    &
                              fldname_h2s,                                     &
                              fldname_so3,                                     &
                              fldname_passive_o3,                              &
                              fldname_age_of_air,                              &
                              fldname_no2,                                     &
                              fldname_hcl,                                     &
                              fldname_bro,                                     &
                              fldname_o1d,                                     &
                              fldname_meoo,                                    &
                              fldname_meco3,                                   &
                              fldname_etoo,                                    &
                              fldname_iso2,                                    &
                              fldname_n_proo,                                  &
                              fldname_i_proo,                                  &
                              fldname_etco3,                                   &
                              fldname_mecoch2oo,                               &
                              fldname_macro2,                                  &
                              fldname_h2o2,                                    &
                              fldname_dms,                                     &
                              fldname_so2,                                     &
                              fldname_h2so4,                                   &
                              fldname_dmso,                                    &
                              fldname_monoterpene,                             &
                              fldname_secondary_organic,                       &
                              fldname_n_nuc_sol,                               &
                              fldname_nuc_sol_su,                              &
                              fldname_nuc_sol_om,                              &
                              fldname_n_ait_sol,                               &
                              fldname_ait_sol_su,                              &
                              fldname_ait_sol_bc,                              &
                              fldname_ait_sol_om,                              &
                              fldname_n_acc_sol,                               &
                              fldname_acc_sol_su,                              &
                              fldname_acc_sol_bc,                              &
                              fldname_acc_sol_om,                              &
                              fldname_acc_sol_ss,                              &
                              fldname_acc_sol_du,                              &
                              fldname_n_cor_sol,                               &
                              fldname_cor_sol_su,                              &
                              fldname_cor_sol_bc,                              &
                              fldname_cor_sol_om,                              &
                              fldname_cor_sol_ss,                              &
                              fldname_cor_sol_du,                              &
                              fldname_n_ait_ins,                               &
                              fldname_ait_ins_bc,                              &
                              fldname_ait_ins_om,                              &
                              fldname_n_acc_ins,                               &
                              fldname_acc_ins_du,                              &
                              fldname_n_cor_ins,                               &
                              fldname_cor_ins_du,                              &
                              fldname_cloud_drop_no_conc,                      &
                              fldname_surfarea,                                &
                              fldname_drydp_ait_sol,                           &
                              fldname_drydp_acc_sol,                           &
                              fldname_drydp_cor_sol,                           &
                              fldname_drydp_ait_ins,                           &
                              fldname_drydp_acc_ins,                           &
                              fldname_drydp_cor_ins,                           &
                              fldname_wetdp_ait_sol,                           &
                              fldname_wetdp_acc_sol,                           &
                              fldname_wetdp_cor_sol,                           &
                              fldname_rhopar_ait_sol,                          &
                              fldname_rhopar_acc_sol,                          &
                              fldname_rhopar_cor_sol,                          &
                              fldname_rhopar_ait_ins,                          &
                              fldname_rhopar_acc_ins,                          &
                              fldname_rhopar_cor_ins,                          &
                              fldname_pvol_su_ait_sol,                         &
                              fldname_pvol_bc_ait_sol,                         &
                              fldname_pvol_om_ait_sol,                         &
                              fldname_pvol_wat_ait_sol,                        &
                              fldname_pvol_su_acc_sol,                         &
                              fldname_pvol_bc_acc_sol,                         &
                              fldname_pvol_om_acc_sol,                         &
                              fldname_pvol_ss_acc_sol,                         &
                              fldname_pvol_du_acc_sol,                         &
                              fldname_pvol_wat_acc_sol,                        &
                              fldname_pvol_su_cor_sol,                         &
                              fldname_pvol_bc_cor_sol,                         &
                              fldname_pvol_om_cor_sol,                         &
                              fldname_pvol_ss_cor_sol,                         &
                              fldname_pvol_du_cor_sol,                         &
                              fldname_pvol_wat_cor_sol,                        &
                              fldname_pvol_bc_ait_ins,                         &
                              fldname_pvol_om_ait_ins,                         &
                              fldname_pvol_du_acc_ins,                         &
                              fldname_pvol_du_cor_ins,                         &
                              fldname_sin_declination,                         &
                              fldname_equation_of_time,                        &
                              fldname_kent,                                    &
                              fldname_kent_dsc,                                &
                              fldname_latitude,                                &
                              fldname_longitude,                               &
                              fldname_sin_latitude,                            &
                              fldname_cos_latitude,                            &
                              fldname_tan_latitude,                            &
                              fldname_frac_seaice,                             &
                              fldname_tstar,                                   &
                              fldname_pstar,                                   &
                              fldname_rough_length,                            &
                              fldname_ustar,                                   &
                              fldname_surf_hf,                                 &
                              fldname_zbl,                                     &
                              fldname_u_scalar_10m,                            &
                              fldname_chloro_sea,                              &
                              fldname_dms_sea_conc,                            &
                              fldname_dust_flux_div1,                          &
                              fldname_dust_flux_div2,                          &
                              fldname_dust_flux_div3,                          &
                              fldname_dust_flux_div4,                          &
                              fldname_dust_flux_div5,                          &
                              fldname_dust_flux_div6,                          &
                              fldname_zhsc,                                    &
                              fldname_surf_wetness,                            &
                              fldname_l_land,                                  &
                              fldname_stcon,                                   &
                              fldname_theta,                                   &
                              fldname_exner_theta_lev,                         &
                              fldname_p_theta_lev,                             &
                              fldname_p_rho_lev,                               &
                              fldname_q,                                       &
                              fldname_qcl,                                     &
                              fldname_qcf,                                     &
                              fldname_bulk_cloud_frac,                         &
                              fldname_ls_rain3d,                               &
                              fldname_ls_snow3d,                               &
                              fldname_conv_rain3d,                             &
                              fldname_conv_snow3d,                             &
                              fldname_rho_r2,                                  &
                              fldname_o3_offline,                              &
                              fldname_no3_offline,                             &
                              fldname_oh_offline,                              &
                              fldname_ho2_offline,                             &
                              fldname_h2o2_limit,                              &
                              fldname_liq_cloud_frac,                          &
                              fldname_autoconv,                                &
                              fldname_accretion,                               &
                              fldname_rim_cry,                                 &
                              fldname_rim_agg,                                 &
                              fldname_vertvel,                                 &
                              fldname_interf_z,                                &
                              fldname_exner_rho_lev,                           &
                              fldname_rhokh_rdz,                               &
                              fldname_dtrdz,                                   &
                              fldname_bl_tke,                                  &
                              fldname_we_lim,                                  &
                              fldname_t_frac,                                  &
                              fldname_zrzi,                                    &
                              fldname_we_lim_dsc,                              &
                              fldname_t_frac_dsc,                              &
                              fldname_zrzi_dsc,                                &
                              fldname_frac_land,                               &
                              fldname_soil_moisture_layer1,                    &
                              fldname_frac_surft,                              &
                              fldname_tstar_surft,                             &
                              fldname_z0_surft,                                &
                              fldname_l_active_surft,                          &
                              fldname_lai_pft,                                 &
                              fldname_canht_pft,                               &
                              nlev_ent_tr_mix

  use log_mod,          only: log_event, log_scratch_space, LOG_LEVEL_ERROR
  use chemistry_config_mod, only: chem_scheme, chem_scheme_strattrop

  ! UM modules
  use nlsizes_namelist_mod, only: bl_levels
  use planet_constants_mod,  only: p_zero, kappa, planet_radius
  use timestep_mod, only: timestep

  ! JULES modules
  use jules_surface_types_mod, only: npft, ntype
  use jules_surface_mod,        only: l_urban2t
  use jules_sea_seaice_mod, only: nice
  use jules_urban_mod,           only: l_moruses
  use sparm_mod, only: sparm
  use tilepts_mod, only: tilepts
  use ancil_info, only: ainfo_data_type, ainfo_type, ancil_info_alloc,         &
          ancil_info_dealloc, ancil_info_nullify, ancil_info_assoc, nsoilt
  use urban_param_mod, only: urban_param_data_type, urban_param_type,          &
          urban_param_alloc, urban_param_assoc, urban_param_nullify,           &
          urban_param_dealloc
  use theta_field_sizes,        only: t_i_length, t_j_length

  ! UKCA API module
  use ukca_api_mod, only: ukca_step_control, ukca_maxlen_message,              &
                          ukca_maxlen_procname

  implicit none

  ! Arguments

  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_w3
  integer(kind=i_def), intent(in) :: undf_w3
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), intent(in) :: ndf_wth
  integer(kind=i_def), intent(in) :: undf_wth
  integer(kind=i_def), dimension(ndf_wth), intent(in) :: map_wth
  integer(kind=i_def), intent(in) :: ndf_tile
  integer(kind=i_def), intent(in) :: undf_tile
  integer(kind=i_def), dimension(ndf_tile), intent(in) :: map_tile
  integer(kind=i_def), intent(in) :: ndf_pft
  integer(kind=i_def), intent(in) :: undf_pft
  integer(kind=i_def), dimension(ndf_pft), intent(in) :: map_pft
  integer(kind=i_def), intent(in) :: ndf_soil
  integer(kind=i_def), intent(in) :: undf_soil
  integer(kind=i_def), dimension(ndf_soil), intent(in) :: map_soil
  integer(kind=i_def), intent(in) :: ndf_2d
  integer(kind=i_def), intent(in) :: undf_2d
  integer(kind=i_def), dimension(ndf_2d), intent(in) :: map_2d
  integer(kind=i_def), intent(in) :: ndf_ent
  integer(kind=i_def), intent(in) :: undf_ent
  integer(kind=i_def), dimension(ndf_ent), intent(in) :: map_ent
  integer(kind=i_def), intent(in) :: ndf_dust
  integer(kind=i_def), intent(in) :: undf_dust
  integer(kind=i_def), dimension(ndf_dust), intent(in) :: map_dust

  real(kind=r_def), intent(in out), dimension(undf_wth) :: o3p
  real(kind=r_def), intent(in out), dimension(undf_wth) :: o1d
  real(kind=r_def), intent(in out), dimension(undf_wth) :: o3
  real(kind=r_def), intent(in out), dimension(undf_wth) :: n
  real(kind=r_def), intent(in out), dimension(undf_wth) :: no
  real(kind=r_def), intent(in out), dimension(undf_wth) :: no3
  real(kind=r_def), intent(in out), dimension(undf_wth) :: no2
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
  real(kind=r_def), intent(in out), dimension(undf_wth) :: bro
  real(kind=r_def), intent(in out), dimension(undf_wth) :: brcl
  real(kind=r_def), intent(in out), dimension(undf_wth) :: brono2
  real(kind=r_def), intent(in out), dimension(undf_wth) :: n2o
  real(kind=r_def), intent(in out), dimension(undf_wth) :: lumped_cl
  real(kind=r_def), intent(in out), dimension(undf_wth) :: hcl
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
  real(kind=r_def), intent(in out), dimension(undf_wth) :: acc_sol_du
  real(kind=r_def), intent(in out), dimension(undf_wth) :: n_cor_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: cor_sol_su
  real(kind=r_def), intent(in out), dimension(undf_wth) :: cor_sol_bc
  real(kind=r_def), intent(in out), dimension(undf_wth) :: cor_sol_om
  real(kind=r_def), intent(in out), dimension(undf_wth) :: cor_sol_ss
  real(kind=r_def), intent(in out), dimension(undf_wth) :: cor_sol_du
  real(kind=r_def), intent(in out), dimension(undf_wth) :: n_ait_ins
  real(kind=r_def), intent(in out), dimension(undf_wth) :: ait_ins_bc
  real(kind=r_def), intent(in out), dimension(undf_wth) :: ait_ins_om
  real(kind=r_def), intent(in out), dimension(undf_wth) :: n_acc_ins
  real(kind=r_def), intent(in out), dimension(undf_wth) :: acc_ins_du
  real(kind=r_def), intent(in out), dimension(undf_wth) :: n_cor_ins
  real(kind=r_def), intent(in out), dimension(undf_wth) :: cor_ins_du
  real(kind=r_def), intent(in out), dimension(undf_wth) :: cloud_drop_no_conc
  real(kind=r_def), intent(in out), dimension(undf_wth) :: drydp_ait_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: drydp_acc_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: drydp_cor_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: drydp_ait_ins
  real(kind=r_def), intent(in out), dimension(undf_wth) :: drydp_acc_ins
  real(kind=r_def), intent(in out), dimension(undf_wth) :: drydp_cor_ins
  real(kind=r_def), intent(in out), dimension(undf_wth) :: wetdp_ait_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: wetdp_acc_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: wetdp_cor_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: rhopar_ait_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: rhopar_acc_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: rhopar_cor_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: rhopar_ait_ins
  real(kind=r_def), intent(in out), dimension(undf_wth) :: rhopar_acc_ins
  real(kind=r_def), intent(in out), dimension(undf_wth) :: rhopar_cor_ins
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_wat_ait_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_wat_acc_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_wat_cor_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_su_ait_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_bc_ait_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_om_ait_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_su_acc_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_bc_acc_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_om_acc_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_ss_acc_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_du_acc_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_su_cor_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_bc_cor_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_om_cor_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_ss_cor_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_du_cor_sol
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_bc_ait_ins
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_om_ait_ins
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_du_acc_ins
  real(kind=r_def), intent(in out), dimension(undf_wth) :: pvol_du_cor_ins

  integer(kind=i_timestep), intent(in) :: timestep_number
  integer(kind=i_def), intent(in) :: current_time_year
  integer(kind=i_def), intent(in) :: current_time_month
  integer(kind=i_def), intent(in) :: current_time_day
  integer(kind=i_def), intent(in) :: current_time_hour
  integer(kind=i_def), intent(in) :: current_time_minute
  integer(kind=i_def), intent(in) :: current_time_second
  integer(kind=i_def), intent(in) :: current_time_daynum
  integer(kind=i_def), intent(in) :: previous_time_year
  integer(kind=i_def), intent(in) :: previous_time_month
  integer(kind=i_def), intent(in) :: previous_time_day
  integer(kind=i_def), intent(in) :: previous_time_hour
  integer(kind=i_def), intent(in) :: previous_time_minute
  integer(kind=i_def), intent(in) :: previous_time_second
  integer(kind=i_def), intent(in) :: previous_time_daynum

  real(kind=r_def), intent(in), dimension(undf_wth) :: theta_wth
  real(kind=r_def), intent(in), dimension(undf_w3) :: exner_in_w3
  real(kind=r_def), intent(in), dimension(undf_wth) :: exner_in_wth
  real(kind=r_def), intent(in), dimension(undf_wth) :: u3_in_wth
  real(kind=r_def), intent(in), dimension(undf_wth) :: m_v_n
  real(kind=r_def), intent(in), dimension(undf_wth) :: m_cl_n
  real(kind=r_def), intent(in), dimension(undf_wth) :: m_ci_n
  real(kind=r_def), intent(in), dimension(undf_w3) :: height_w3
  real(kind=r_def), intent(in), dimension(undf_wth) :: height_wth
  real(kind=r_def), intent(in), dimension(undf_wth) :: ls_rain_3d
  real(kind=r_def), intent(in), dimension(undf_wth) :: ls_snow_3d
  real(kind=r_def), intent(in), dimension(undf_wth) :: autoconv
  real(kind=r_def), intent(in), dimension(undf_wth) :: accretion
  real(kind=r_def), intent(in), dimension(undf_wth) :: rim_cry
  real(kind=r_def), intent(in), dimension(undf_wth) :: rim_agg
  real(kind=r_def), intent(in), dimension(undf_wth) :: tke_bl
  real(kind=r_def), intent(in), dimension(undf_wth) :: conv_rain_3d
  real(kind=r_def), intent(in), dimension(undf_wth) :: conv_snow_3d
  real(kind=r_def), intent(in), dimension(undf_wth) :: cf_bulk
  real(kind=r_def), intent(in), dimension(undf_wth) :: cf_liq
  real(kind=r_def), intent(in), dimension(undf_tile) :: tile_fraction
  real(kind=r_def), intent(in), dimension(undf_tile) :: tile_temperature
  real(kind=r_def), intent(in), dimension(undf_tile) :: tile_heat_flux
  real(kind=r_def), intent(in), dimension(undf_tile) :: gc_tile
  real(kind=r_def), intent(in), dimension(undf_pft) :: leaf_area_index
  real(kind=r_def), intent(in), dimension(undf_pft) :: canopy_height
  real(kind=r_def), intent(in), dimension(undf_soil) :: soil_moisture
  real(kind=r_def), intent(in), dimension(undf_2d) :: latitude
  real(kind=r_def), intent(in), dimension(undf_2d) :: longitude
  real(kind=r_def), intent(in), dimension(undf_2d) ::                          &
    sin_stellar_declination_rts
  real(kind=r_def), intent(in), dimension(undf_2d) :: stellar_eqn_of_time_rts
  real(kind=r_def), intent(in), dimension(undf_2d) :: soil_roughness
  real(kind=r_def), intent(in), dimension(undf_2d) :: z0m
  real(kind=r_def), intent(in), dimension(undf_2d) :: urbztm
  real(kind=r_def), intent(in), dimension(undf_2d) :: ustar
  real(kind=r_def), intent(in), dimension(undf_2d) :: wspd10m
  real(kind=r_def), intent(in), dimension(undf_2d) :: chloro_sea
  real(kind=r_def), intent(in), dimension(undf_2d) :: zh
  real(kind=r_def), intent(in), dimension(undf_w3) :: wetrho_in_w3
  real(kind=r_def), intent(in), dimension(undf_w3) :: rhokh_bl
  real(kind=r_def), intent(in), dimension(undf_w3) :: rdz_tq_bl
  real(kind=r_def), intent(in), dimension(undf_wth) :: dtrdz_tq_bl
  real(kind=r_def), intent(in), dimension(undf_2d) :: zhsc
  integer(kind=i_def), intent(in), dimension(undf_2d) :: level_ent
  integer(kind=i_def), intent(in), dimension(undf_2d) :: level_ent_dsc
  real(kind=r_def), intent(in), dimension(undf_ent) :: ent_we_lim
  real(kind=r_def), intent(in), dimension(undf_ent) :: ent_t_frac
  real(kind=r_def), intent(in), dimension(undf_ent) :: ent_zrzi
  real(kind=r_def), intent(in), dimension(undf_ent) :: ent_we_lim_dsc
  real(kind=r_def), intent(in), dimension(undf_ent) :: ent_t_frac_dsc
  real(kind=r_def), intent(in), dimension(undf_ent) :: ent_zrzi_dsc
  real(kind=r_def), intent(in), dimension(undf_wth) :: h2o2_limit
  real(kind=r_def), intent(in), dimension(undf_2d) :: dms_conc_ocean
  real(kind=r_def), intent(in), dimension(undf_dust) :: dust_flux
  real(kind=r_def), intent(inout), dimension(undf_2d) :: surf_wetness
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_c2h6
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_c3h8
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_c5h8
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_ch4
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_co
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_hcho
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_me2co
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_mecho
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_meoh
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_nh3
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_no
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_bc_biofuel
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_bc_biomass_high
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_bc_biomass_low
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_bc_fossil
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_dms_land
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_monoterp
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_om_biofuel
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_om_biomass_high
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_om_biomass_low
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_om_fossil
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_so2_high
  real(kind=r_def), intent(in), dimension(undf_2d) :: emiss_so2_low
  real(kind=r_def), intent(in), dimension(undf_wth) :: emiss_no_aircrft
  real(kind=r_def), intent(in), dimension(undf_wth) :: emiss_bc_biomass
  real(kind=r_def), intent(in), dimension(undf_wth) :: emiss_om_biomass
  real(kind=r_def), intent(in), dimension(undf_wth) :: emiss_so2_nat


  ! Local variables for the kernel

  ! Current model time (year, month, day, hour, minute, second, day number)
  integer(i_um) :: current_time(7)
  ! Model time at previous time step
  integer(i_um) :: previous_time(7)

  ! Prognostics to be updated in the time step
  real(r_um), allocatable :: tracer(:,:) ! Tracer fields
  real(r_um), allocatable :: ntp(:,:)    ! NTP fields
                                         ! (for microphysics & RADAER)

  real(r_um) :: r_theta_levels(1,1,0:nlayers)
  real(r_um) :: r_rho_levels(1,1,nlayers)
  ! Environmental driver fields (including emissions)

  integer(i_um), allocatable :: environ_flat_integer(:)

  real(r_um), allocatable :: environ_scalar_real(:)
  real(r_um), allocatable :: environ_flat_real(:)
  real(r_um), allocatable :: environ_flatpft_real(:,:)
  real(r_um), allocatable :: environ_fullht_real(:,:)
  real(r_um), allocatable :: environ_fullht0_real(:,:)
  real(r_um), allocatable :: environ_fullhtp1_real(:,:)
  real(r_um), allocatable :: environ_bllev_real(:,:)
  real(r_um), allocatable :: environ_entlev_real(:,:)
  real(r_um), allocatable :: environ_land_real(:,:)
  real(r_um), allocatable :: environ_landtile_real(:,:,:)
  real(r_um), allocatable :: environ_landpft_real(:,:,:)
  real(r_um), allocatable :: emissions_flat(:)
  real(r_um), allocatable :: emissions_fullht(:,:)

  logical, allocatable :: environ_flat_logical(:)
  logical, allocatable :: environ_landtile_logical(:,:,:)

  ! Working variables

  integer(i_um) :: n_fields   ! Number of fields in a group
  integer(i_um) :: n_land_pts ! Number of land points
  integer(i_um) :: i
  integer(i_um) :: j
  integer(i_um) :: k
  integer(i_um) :: error_code

  integer(i_um) :: surft_index(n_land_tile)
                              ! Indices of land points that include surface type
  real(r_um) :: frac_land     ! Land fraction of cell
  real(r_um) :: frac_sea      ! Sea fraction of cell
  real(r_um) :: frac_seaice   ! Sea fraction with respect to sea area
  real(r_um) :: meanval       ! Arbitrary mean value

  real(r_um) :: z0m_soil_gb(1) ! Soil roughness length (m)

  real(r_um) :: frac_surft( 1, n_land_tile )
                              ! Fraction of surf. type with respect to land area
  real(r_um) :: z0_surft( 1, n_land_tile )
                              ! Surface roughness length on tiles (m)
  real(r_um) :: lai_pft( 1, n_land_tile )
                              ! Leaf area index of plant functional type
  real(r_um) :: canht_pft( 1, n_land_tile )
                              ! Canopy_height of plant functional type (m)
  real(r_um) :: rho_r2        ! Density * r * r
  real(r_um) :: exner_rho_top ! Exner pressure at top rho level
  real(r_um) :: exner_theta_top ! Exner pressure at top theta level

  logical :: l_land           ! Land/sea indicator (True for land point)

  logical :: l_tile_active( 1, n_land_tile )
                              ! active tile indicator (True if tile is in use)

  ! Unused output fields from JULES sparm routine
  real(r_um) :: catch_snow_surft( 1, n_land_tile )
  real(r_um) :: catch_surft( 1, n_land_tile )
  real(r_um) :: z0h_bare_surft( 1, n_land_tile )

  type(ainfo_data_type) :: ainfo_data
  type(ainfo_type) :: ainfo
  type(urban_param_data_type) :: urban_param_data
  type(urban_param_type) :: urban_param

  ! surface wetness calculation
  real(r_def) :: tot_precip
  real(r_def), parameter :: tol = 1.0e-10_r_def
  real(r_def), parameter :: raincrit = 0.5_r_def/rsec_per_day
  real(r_def), parameter :: ztodry = 3.0_r_def * rsec_per_hour

  ! UKCA error reporting variables
  character(len=ukca_maxlen_message)  :: ukca_errmsg  ! Error return message
  character(len=ukca_maxlen_procname) :: ukca_errproc ! Routine in which error
                                                      ! was trapped

  !-----------------------------------------------------------------------
  ! Map LFRic fields into UKCA arrays
  !-----------------------------------------------------------------------

  ! -- Tracers --

  n_fields = size(tracer_names)
  allocate(tracer( nlayers, n_fields ))
  tracer = 0.0_r_um

  do i = 1, n_fields
    select case(tracer_names(i))
    case(fldname_o3p)
      tracer( :, i ) =                                                         &
        real( o3p( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_o3)
      tracer( :, i ) =                                                         &
        real( o3( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_n)
      tracer( :, i ) =                                                         &
        real( n( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_no)
      tracer( :, i ) =                                                         &
        real( no( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_no3)
      tracer( :, i ) =                                                         &
        real( no3( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_lumped_n)
      tracer( :, i ) =                                                         &
        real( lumped_n( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_n2o5)
      tracer( :, i ) =                                                         &
        real( n2o5( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_ho2no2)
      tracer( :, i ) =                                                         &
        real( ho2no2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_hono2)
      tracer( :, i ) =                                                         &
        real( hono2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_h2o2)
      tracer( :, i ) =                                                         &
        real( h2o2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_ch4)
      tracer( :, i ) =                                                         &
        real( ch4( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_co)
      tracer( :, i ) =                                                         &
        real( co( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_hcho)
      tracer( :, i ) =                                                         &
        real( hcho( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_meoo)
      tracer( :, i ) =                                                         &
        real( meoo( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_meooh)
      tracer( :, i ) =                                                         &
        real( meooh( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_h)
      tracer( :, i ) =                                                         &
        real( h( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_ch2o)
      ! Copy water mixing ratio in H2O tracer in chemistry
      tracer( :, i ) =                                                         &
        real( m_v_n( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_oh)
      tracer( :, i ) =                                                         &
        real( oh( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_ho2)
      tracer( :, i ) =                                                         &
        real( ho2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_cl)
      tracer( :, i ) =                                                         &
        real( cl( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_cl2o2)
      tracer( :, i ) =                                                         &
        real( cl2o2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_clo)
      tracer( :, i ) =                                                         &
        real( clo( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_oclo)
      tracer( :, i ) =                                                         &
        real( oclo( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_br)
      tracer( :, i ) =                                                         &
        real( br( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_lumped_br)
      tracer( :, i ) =                                                         &
        real( lumped_br( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_brcl)
      tracer( :, i ) =                                                         &
        real( brcl( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_brono2)
      tracer( :, i ) =                                                         &
        real( brono2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_n2o)
      tracer( :, i ) =                                                         &
        real( n2o( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_lumped_cl)
      tracer( :, i ) =                                                         &
        real( lumped_cl( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_hocl)
      tracer( :, i ) =                                                         &
        real( hocl( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_hbr)
      tracer( :, i ) =                                                         &
        real( hbr( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_hobr)
      tracer( :, i ) =                                                         &
        real( hobr( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_clono2)
      tracer( :, i ) =                                                         &
        real( clono2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_cfcl3)
      tracer( :, i ) =                                                         &
        real( cfcl3( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_cf2cl2)
      tracer( :, i ) =                                                         &
        real( cf2cl2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_mebr)
      tracer( :, i ) =                                                         &
        real( mebr( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_hono)
      tracer( :, i ) =                                                         &
        real( hono( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_c2h6)
      tracer( :, i ) =                                                         &
        real( c2h6( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_etoo)
      tracer( :, i ) =                                                         &
        real( etoo( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_etooh)
      tracer( :, i ) =                                                         &
        real( etooh( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_mecho)
      tracer( :, i ) =                                                         &
        real( mecho( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_meco3)
      tracer( :, i ) =                                                         &
        real( meco3( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pan)
      tracer( :, i ) =                                                         &
        real( pan( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_c3h8)
      tracer( :, i ) =                                                         &
        real( c3h8( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_n_proo)
      tracer( :, i ) =                                                         &
        real( n_proo( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_i_proo)
      tracer( :, i ) =                                                         &
        real( i_proo( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_n_prooh)
      tracer( :, i ) =                                                         &
        real( n_prooh( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_i_prooh)
      tracer( :, i ) =                                                         &
        real( i_prooh( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_etcho)
      tracer( :, i ) =                                                         &
        real( etcho( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_etco3)
      tracer( :, i ) =                                                         &
        real( etco3( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_me2co)
      tracer( :, i ) =                                                         &
        real( me2co( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_mecoch2oo)
      tracer( :, i ) =                                                         &
        real( mecoch2oo( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_mecoch2ooh)
      tracer( :, i ) =                                                         &
        real( mecoch2ooh( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_ppan)
      tracer( :, i ) =                                                         &
        real( ppan( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_meono2)
      tracer( :, i ) =                                                         &
        real( meono2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_c5h8)
      tracer( :, i ) =                                                         &
        real( c5h8( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_iso2)
      tracer( :, i ) =                                                         &
        real( iso2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_isooh)
      tracer( :, i ) =                                                         &
        real( isooh( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_ison)
      tracer( :, i ) =                                                         &
        real( ison( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_macr)
      tracer( :, i ) =                                                         &
        real( macr( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_macrooh)
      tracer( :, i ) =                                                         &
        real( macrooh( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_macro2)
      tracer( :, i ) =                                                         &
        real( macro2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_mpan)
      tracer( :, i ) =                                                         &
        real( mpan( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_hacet)
      tracer( :, i ) =                                                         &
        real( hacet( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_mgly)
      tracer( :, i ) =                                                         &
        real( mgly( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_nald)
      tracer( :, i ) =                                                         &
        real( nald( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_hcooh)
      tracer( :, i ) =                                                         &
        real( hcooh( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_meco3h)
      tracer( :, i ) =                                                         &
        real( meco3h( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_meco2h)
      tracer( :, i ) =                                                         &
        real( meco2h( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_h2)
      tracer( :, i ) =                                                         &
        real( h2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_meoh)
      tracer( :, i ) =                                                         &
        real( meoh( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_msa)
      tracer( :, i ) =                                                         &
        real( msa( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_nh3)
      tracer( :, i ) =                                                         &
        real( nh3( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_cs2)
      tracer( :, i ) =                                                         &
        real( cs2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_csul)
      tracer( :, i ) =                                                         &
        real( csul( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_h2s)
      tracer( :, i ) =                                                         &
        real( h2s( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_so3)
      tracer( :, i ) =                                                         &
        real( so3( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_passive_o3)
      tracer( :, i ) =                                                         &
        real( passive_o3( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_age_of_air)
      tracer( :, i ) =                                                         &
        real( age_of_air( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_dms)
      tracer( :, i ) =                                                         &
        real( dms( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_so2)
      tracer( :, i ) =                                                         &
        real( so2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_h2so4)
      tracer( :, i ) =                                                         &
        real( h2so4( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_dmso)
      tracer( :, i ) =                                                         &
        real( dmso( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_monoterpene)
      tracer( :, i ) =                                                         &
        real( monoterpene( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_secondary_organic)
      tracer( :, i ) =                                                         &
        real( secondary_organic( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_n_nuc_sol)
      tracer( :, i ) =                                                         &
        real( n_nuc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_nuc_sol_su)
      tracer( :, i ) =                                                         &
        real( nuc_sol_su( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_nuc_sol_om)
      tracer( :, i ) =                                                         &
        real( nuc_sol_om( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_n_ait_sol)
      tracer( :, i ) =                                                         &
        real( n_ait_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_ait_sol_su)
      tracer( :, i ) =                                                         &
        real( ait_sol_su( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_ait_sol_bc)
      tracer( :, i ) =                                                         &
        real( ait_sol_bc( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_ait_sol_om)
      tracer( :, i ) =                                                         &
        real( ait_sol_om( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_n_acc_sol)
      tracer( :, i ) =                                                         &
        real( n_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_acc_sol_su)
      tracer( :, i ) =                                                         &
        real( acc_sol_su( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_acc_sol_bc)
      tracer( :, i ) =                                                         &
        real( acc_sol_bc( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_acc_sol_om)
      tracer( :, i ) =                                                         &
        real( acc_sol_om( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_acc_sol_ss)
      tracer( :, i ) =                                                         &
        real( acc_sol_ss( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_acc_sol_du)
      tracer( :, i ) =                                                         &
        real( acc_sol_du( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_n_cor_sol)
      tracer( :, i ) =                                                         &
        real( n_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_cor_sol_su)
      tracer( :, i ) =                                                         &
        real( cor_sol_su( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_cor_sol_bc)
      tracer( :, i ) =                                                         &
        real( cor_sol_bc( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_cor_sol_om)
      tracer( :, i ) =                                                         &
        real( cor_sol_om( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_cor_sol_ss)
      tracer( :, i ) =                                                         &
        real( cor_sol_ss( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_cor_sol_du)
      tracer( :, i ) =                                                         &
        real( cor_sol_du( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_n_ait_ins)
      tracer( :, i ) =                                                         &
        real( n_ait_ins( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_ait_ins_bc)
      tracer( :, i ) =                                                         &
        real( ait_ins_bc( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_ait_ins_om)
      tracer( :, i ) =                                                         &
        real( ait_ins_om( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_n_acc_ins)
      tracer( :, i ) =                                                         &
        real( n_acc_ins( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_acc_ins_du)
      tracer( :, i ) =                                                         &
        real( acc_ins_du( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_n_cor_ins)
      tracer( :, i ) =                                                         &
        real( n_cor_ins( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_cor_ins_du)
      tracer( :, i ) =                                                         &
        real( cor_ins_du( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case default
      write( log_scratch_space, '(A,A)' )                                      &
        'Missing required UKCA tracer field: ', tracer_names(i)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  ! -- Non-transported prognostics --
  n_fields = size(ntp_names)
  allocate(ntp(nlayers, n_fields))
  ntp = 0.0_r_um

  do i = 1, n_fields
    select case(ntp_names(i))
    case(fldname_cloud_drop_no_conc)
      ntp( :, i ) =                                                            &
        real( cloud_drop_no_conc( map_wth(1) + 1 : map_wth(1) + nlayers ),     &
              r_um )
    case(fldname_surfarea)
      ! Aerosol surface area is not used in this configuration but UKCA
      ! requires it to be present for potential use in diagnostic calculations.
      ! Pass zero for now pending UKCA API extension to support diagnostics.
    case(fldname_drydp_ait_sol)
      ntp( :, i ) =                                                            &
        real( drydp_ait_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_drydp_acc_sol)
      ntp( :, i ) =                                                            &
        real( drydp_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_drydp_cor_sol)
      ntp( :, i ) =                                                            &
        real( drydp_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_drydp_ait_ins)
      ntp( :, i ) =                                                            &
        real( drydp_ait_ins( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_drydp_acc_ins)
      ntp( :, i ) =                                                            &
        real( drydp_acc_ins( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_drydp_cor_ins)
      ntp( :, i ) =                                                            &
        real( drydp_cor_ins( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_wetdp_ait_sol)
      ntp( :, i ) =                                                            &
        real( wetdp_ait_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_wetdp_acc_sol)
      ntp( :, i ) =                                                            &
        real( wetdp_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_wetdp_cor_sol)
      ntp( :, i ) =                                                            &
        real( wetdp_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_rhopar_ait_sol)
      ntp( :, i ) =                                                            &
        real( rhopar_ait_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_rhopar_acc_sol)
      ntp( :, i ) =                                                            &
        real( rhopar_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_rhopar_cor_sol)
      ntp( :, i ) =                                                            &
        real( rhopar_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_rhopar_ait_ins)
      ntp( :, i ) =                                                            &
        real( rhopar_ait_ins( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_rhopar_acc_ins)
      ntp( :, i ) =                                                            &
        real( rhopar_acc_ins( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_rhopar_cor_ins)
      ntp( :, i ) =                                                            &
        real( rhopar_cor_ins( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_su_ait_sol)
      ntp( :, i ) =                                                            &
        real( pvol_su_ait_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_bc_ait_sol)
      ntp( :, i ) =                                                            &
        real( pvol_bc_ait_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_om_ait_sol)
      ntp( :, i ) =                                                            &
        real( pvol_om_ait_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_wat_ait_sol)
      ntp( :, i ) =                                                            &
        real( pvol_wat_ait_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_su_acc_sol)
      ntp( :, i ) =                                                            &
        real( pvol_su_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_bc_acc_sol)
      ntp( :, i ) =                                                            &
        real( pvol_bc_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_om_acc_sol)
      ntp( :, i ) =                                                            &
        real( pvol_om_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_ss_acc_sol)
      ntp( :, i ) =                                                            &
        real( pvol_ss_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_du_acc_sol)
      ntp( :, i ) =                                                            &
        real( pvol_du_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_wat_acc_sol)
      ntp( :, i ) =                                                            &
        real( pvol_wat_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_su_cor_sol)
      ntp( :, i ) =                                                            &
        real( pvol_su_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_bc_cor_sol)
      ntp( :, i ) =                                                            &
        real( pvol_bc_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_om_cor_sol)
      ntp( :, i ) =                                                            &
        real( pvol_om_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_ss_cor_sol)
      ntp( :, i ) =                                                            &
        real( pvol_ss_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_du_cor_sol)
      ntp( :, i ) =                                                            &
        real( pvol_du_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_wat_cor_sol)
      ntp( :, i ) =                                                            &
        real( pvol_wat_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_bc_ait_ins)
      ntp( :, i ) =                                                            &
        real( pvol_bc_ait_ins( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_om_ait_ins)
      ntp( :, i ) =                                                            &
        real( pvol_om_ait_ins( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_du_acc_ins)
      ntp( :, i ) =                                                            &
        real( pvol_du_acc_ins( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_pvol_du_cor_ins)
      ntp( :, i ) =                                                            &
        real( pvol_du_cor_ins( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_no2)
      ntp( :, i ) =                                                            &
        real( no2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_hcl)
      ntp( :, i ) =                                                            &
        real( hcl( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_bro)
      ntp( :, i ) =                                                            &
        real( bro( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_o1d)
      ntp( :, i ) =                                                            &
        real( o1d( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    ! The next set can be tracers or ntp depending on l_ukca_ro2_ntp
    case(fldname_meoo)
      ntp( :, i ) =                                                            &
        real( meoo( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_meco3)
      ntp( :, i ) =                                                            &
        real( meco3( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_etoo)
      ntp( :, i ) =                                                            &
        real( etoo( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_iso2)
      ntp( :, i ) =                                                            &
        real( iso2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_n_proo)
      ntp( :, i ) =                                                            &
        real( n_proo( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_i_proo)
      ntp( :, i ) =                                                            &
        real( i_proo( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_etco3)
      ntp( :, i ) =                                                            &
        real( etco3( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_mecoch2oo)
      ntp( :, i) =                                                             &
        real( mecoch2oo( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_macro2)
      ntp( :, i ) =                                                            &
        real( macro2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case default
      write( log_scratch_space, '(A,A)' )                                      &
        'Missing required UKCA NTP: ', ntp_names(i)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  ! -- Environmental Drivers --

  ! Set up UM level_height_mod driver fields. These are not yet
  ! supported as environmental drivers in the UKCA API.
  ! N.B. These are required for interf_z calculation below

  do k = 0, nlayers
    ! height of theta levels from centre of planet
    r_theta_levels( 1, 1, k ) = real( height_wth( map_wth(1) + k ), r_um ) +   &
                                planet_radius
  end do
  do k = 1, nlayers
    ! height of rho levels from centre of planet
    r_rho_levels( 1, 1, k ) = real( height_w3( map_w3(1) + k - 1 ), r_um ) +   &
                              planet_radius
  end do

  ! Determine no. of land points (0 or 1) and set land/sea indicator

  frac_land = 0.0_r_um
  do i = 1, n_land_tile
    frac_land = frac_land + real( tile_fraction( map_tile(1) + i - 1 ), r_um )
  end do

  if (frac_land > 0.0_r_um) then
    l_land = .true.
    n_land_pts = 1
  else
    l_land = .false.
    n_land_pts = 0
  end if

!  ! Set up JULES fields needed
  call ancil_info_alloc(n_land_pts, t_i_length, t_j_length, nice, nsoilt,      &
                        ntype, ainfo_data)
  call ancil_info_assoc(ainfo, ainfo_data)
  call urban_param_alloc(n_land_pts, l_urban2t, l_moruses, urban_param_data)
  call urban_param_assoc(urban_param, urban_param_data)

  if (l_land) then
    ! Tile fractions with respect to land area
    do i = 1, n_land_tile
      frac_surft( 1, i ) =                                                     &
        real( tile_fraction( map_tile(1) + i - 1 ), r_um ) / frac_land
    end do
    ! Call JULES subroutine to determine number of active points (0 or 1)
    ! for each tile type
    call tilepts( n_land_pts, frac_surft, ainfo%surft_pts, surft_index,        &
                  ainfo%l_lice_point, ainfo%l_lice_surft )
    ! Fields on plant functional type tiles: leaf area index & canopy height
    do i = 1, npft
      lai_pft( 1, i ) = real( leaf_area_index(map_pft(1) + i - 1), r_um )
      canht_pft( 1, i ) = real( canopy_height(map_pft(1) + i - 1), r_um )
    end do
    ! Roughness length on tiles (z0_surft) from JULES
    z0m_soil_gb(1)     = real( soil_roughness(map_2d(1)), r_um )
    urban_param%ztm_gb = real(urbztm(map_2d(1)), r_um)
    call sparm( n_land_pts, n_land_tile, ainfo%surft_pts, surft_index,         &
                frac_surft, canht_pft, lai_pft, z0m_soil_gb,                   &
                catch_snow_surft, catch_surft, z0_surft, z0h_bare_surft,       &
                urban_param%ztm_gb )
  end if

  call urban_param_nullify(urban_param)
  call urban_param_dealloc(urban_param_data)

  ! Drivers in scalar group

  n_fields = size(env_names_scalar_real)
  allocate(environ_scalar_real(n_fields))
  environ_scalar_real = 0

  do i = 1, n_fields
    select case(env_names_scalar_real(i))
    case(fldname_sin_declination)
      environ_scalar_real(i) =                                                 &
        real( sin_stellar_declination_rts(map_2d(1)), r_um )
    case(fldname_equation_of_time)
      environ_scalar_real(i) = real( stellar_eqn_of_time_rts(map_2d(1)), r_um )
    case default
      write( log_scratch_space, '(A,A)' )                                      &
        'Missing required UKCA environment field: ', env_names_scalar_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  ! Drivers in flat grid groups (scalar values)

  n_fields = size(env_names_flat_integer)
  allocate(environ_flat_integer(n_fields))
  environ_flat_integer = 0

  do i = 1, n_fields
    select case(env_names_flat_integer(i))
    case(fldname_kent)
      environ_flat_integer(i) = int( level_ent(map_2d(1)), i_um )
    case(fldname_kent_dsc)
      environ_flat_integer(i) = int( level_ent_dsc(map_2d(1)), i_um )
    case default
      write( log_scratch_space, '(A,A)' )                                      &
        'Missing required UKCA environment field: ', env_names_flat_integer(i)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  n_fields = size(env_names_flat_real)
  allocate(environ_flat_real(n_fields))
  environ_flat_real = 0.0_r_um

  do i = 1, n_fields
    select case(env_names_flat_real(i))
    case(fldname_latitude)
      ! Latitude (degrees N)
      environ_flat_real(i) =                                                   &
        real( radians_to_degrees * latitude(map_2d(1)), r_um )
    case(fldname_longitude)
      ! Longitude (degrees E, >=0, <360)
      environ_flat_real(i) =                                                   &
        real( radians_to_degrees * longitude(map_2d(1)), r_um )
      if (environ_flat_real(i) < 0.0_r_um)                                     &
        environ_flat_real(i) = environ_flat_real(i) + 360.0_r_um
    case(fldname_sin_latitude)
      ! Sin latitude
      environ_flat_real(i) = real( sin(latitude(map_2d(1))), r_um )
    case(fldname_cos_latitude)
      ! Cos latitude
      environ_flat_real(i) = real( cos(latitude(map_2d(1))), r_um )
    case(fldname_tan_latitude)
      ! Tan latitude
      environ_flat_real(i) = real( tan(latitude(map_2d(1))), r_um )
    case(fldname_frac_seaice)
      ! Sea-ice fraction with respect to sea area
      frac_seaice = 0.0_r_um
      frac_sea = 1.0_r_um - frac_land
      if (frac_sea > 0.0_r_um) then
        do j = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
          frac_seaice = frac_seaice +                                          &
            real( tile_fraction( map_tile(1) + j - 1 ), r_um ) / frac_sea
        end do
      end if
      environ_flat_real(i) = frac_seaice
    case(fldname_tstar)
      ! Cell mean surface temperature
      meanval = 0.0_r_um
      do j = 1, n_surf_tile
        meanval = meanval +                                                    &
          real( tile_fraction( map_tile(1) + j - 1 ), r_um ) *                 &
          real( tile_temperature( map_tile(1) + j - 1 ), r_um )
      end do
      environ_flat_real(i) = meanval
    case(fldname_pstar)
      ! Surface air pressure
      environ_flat_real(i) = p_zero *                                          &
        real( exner_in_wth(map_wth(1) + 0), r_um ) ** ( 1.0_r_um / kappa )
    case(fldname_rough_length)
      ! Cell mean surface roughness
      environ_flat_real(i) = real( z0m(map_2d(1)) , r_um)
    case(fldname_ustar)
      ! Friction velocity
      environ_flat_real(i) = real( ustar(map_2d(1)), r_um )
    case(fldname_surf_hf)
      ! Cell mean surface heat flux
      meanval = 0.0_r_um
      do j = 1, n_surf_tile
        meanval = meanval +                                                    &
          real( tile_fraction( map_tile(1) + j - 1 ), r_um ) *                 &
          real( tile_heat_flux( map_tile(1) + j - 1 ), r_um )
      end do
      environ_flat_real(i) = meanval
    case(fldname_zbl)
      ! Boundary layer height
      environ_flat_real(i) = real( zh(map_2d(1)), r_um )
    case(fldname_u_scalar_10m)
      ! Wind speed at 10 m
      environ_flat_real(i) = real( wspd10m(map_2d(1)), r_um )
    case(fldname_chloro_sea)
      environ_flat_real(i) = real( chloro_sea(map_2d(1)), r_um )
    case(fldname_dms_sea_conc)
      ! Sea surface DMS concentration
      ! (Replace fill value with zeros to avoid potential unsafe multiplication
      ! of fill values by zero over land in UKCA DMS flux calculation and/or
      ! use of fill values if present at coastal grid points)
      environ_flat_real(i) = max( real( dms_conc_ocean(map_2d(1)), r_um ),     &
                                  0.0_r_um )
    case(fldname_dust_flux_div1)
      ! Dust emission flux in division 1
      environ_flat_real(i) = real( dust_flux(map_dust(1) + 0), r_um )
    case(fldname_dust_flux_div2)
      ! Dust emission flux in division 2
      environ_flat_real(i) = real( dust_flux(map_dust(1) + 1), r_um )
    case(fldname_dust_flux_div3)
      ! Dust emission flux in division 3
      environ_flat_real(i) = real( dust_flux(map_dust(1) + 2), r_um )
    case(fldname_dust_flux_div4)
      ! Dust emission flux in division 4
      environ_flat_real(i) = real( dust_flux(map_dust(1) + 3), r_um )
    case(fldname_dust_flux_div5)
      ! Dust emission flux in division 5
      environ_flat_real(i) = real( dust_flux(map_dust(1) + 4), r_um )
    case(fldname_dust_flux_div6)
      ! Dust emission flux in division 6
      environ_flat_real(i) = real( dust_flux(map_dust(1) + 5), r_um )
    case(fldname_zhsc)
      ! Height at top of decoupled stratocumulus layer
      environ_flat_real(i) = real( zhsc(map_2d(1)), r_um )
    case(fldname_surf_wetness)
      tot_precip = max(0.0_r_def, ls_rain_3d(map_wth(1)+1)   &
                                + ls_snow_3d(map_wth(1)+1)   &
                                + conv_rain_3d(map_wth(1)+1) &
                                + conv_snow_3d(map_wth(1)+1) )
      if ( surf_wetness(map_2d(1)) < tol) then
        ! Surface was dry last timestep, check if it has rained since
        if (tot_precip > raincrit) then
          surf_wetness(map_2d(1)) = 1.0_r_def
        end if
      else
        ! Surface was wet last timestep, check if more rain since
        if (tot_precip > tol) then
          ! Check if there is enough rain to reset surf_wet to 1.0.
          ! Otherwise, it is raining but less than raincrit so surface
          ! does not dry but is not reset to 1 it just keeps the same
          ! value of surf_wet
          if (tot_precip > raincrit) then
            surf_wetness(map_2d(1)) = 1.0_r_def
          end if
        else
          ! Not raining, so the surface is gradually drying.
          surf_wetness(map_2d(1)) = max(0.0_r_def, surf_wetness(map_2d(1)) &
                                                 - timestep/ztodry)
        end if
      end if
      environ_flat_real(i) = real(surf_wetness(map_2d(1)), r_um)
    case default
      write( log_scratch_space, '(A,A)' )                                      &
        'Missing required UKCA environment field: ', env_names_flat_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  n_fields = size(env_names_flat_logical)
  allocate(environ_flat_logical(n_fields))
  environ_flat_logical = .false.

  do i = 1, n_fields
    select case(env_names_flat_logical(i))
    ! Land-sea mask
    case(fldname_l_land)
      environ_flat_logical(i) = l_land
    case default
      write( log_scratch_space, '(A,A)' )                                      &
        'Missing required UKCA environment field: ', env_names_flat_logical(i)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  ! Drivers in flat grid plant functional type tile group (1D fields)

  n_fields = size(env_names_flatpft_real)
  allocate(environ_flatpft_real( npft, n_fields ))
  environ_flatpft_real = 0.0_r_um

  do i = 1, n_fields
    select case(env_names_flatpft_real(i))
    case(fldname_stcon)
      ! Stomatal conductance
      if (l_land) then
        do j = 1, npft
          environ_flatpft_real( j, i ) =                                       &
            real( gc_tile( map_tile(1) + j - 1 ), r_um )
        end do
      end if
    case default
      write( log_scratch_space, '(A,A)' )                                      &
        'Missing required UKCA environment field: ', env_names_flatpft_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  ! Drivers in full-height grid group (1D fields)

  n_fields = size(env_names_fullht_real)
  allocate(environ_fullht_real( nlayers, n_fields ))
  environ_fullht_real = 0.0_r_um

  do i = 1, n_fields
    select case(env_names_fullht_real(i))
    case(fldname_o3_offline)
      ! Ozone mmr
      environ_fullht_real( :, i ) =                                            &
        real( o3( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_no3_offline)
      ! Nitrate mmr
      environ_fullht_real( :, i ) =                                            &
        real( no3( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_oh_offline)
      ! Hydroxyl radical mmr
      environ_fullht_real( :, i ) =                                            &
        real( oh( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_ho2_offline)
      ! Hydroperoxyl radical mmr
      environ_fullht_real( :, i ) =                                            &
        real( ho2( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_h2o2_limit)
      ! Hydrogen peroxide mmr upper limit
      environ_fullht_real( :, i ) =                                            &
        real( h2o2_limit( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_theta)
      ! Potential temperature
      environ_fullht_real( :, i ) =                                            &
        real( theta_wth( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_exner_theta_lev)
      ! Dimensionless Exner function
      environ_fullht_real( :, i ) =                                            &
        real( exner_in_wth( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_p_theta_lev)
      ! Air pressure on theta levels
      environ_fullht_real( :, i ) = p_zero *                                   &
        real( exner_in_wth( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )    &
        ** ( 1.0_r_um / kappa )
    case(fldname_p_rho_lev)
      ! Air pressure on rho levels
      environ_fullht_real( :, i ) = p_zero *                                   &
        real( exner_in_w3( map_w3(1) + 0 : map_w3(1) + nlayers - 1 ), r_um )   &
        ** ( 1.0_r_um / kappa )
    case(fldname_q)
      ! Water vapour mixing ratio
      environ_fullht_real( :, i ) =                                            &
        real( m_v_n( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_qcl)
      ! Cloud liquid water mixing ratio
      environ_fullht_real( :, i ) =                                            &
        real( m_cl_n( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_qcf)
      ! Cloud ice mixing ratio
      environ_fullht_real( :, i ) =                                            &
        real( m_ci_n( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_bulk_cloud_frac)
      ! Bulk cloud fraction
      environ_fullht_real( :, i ) =                                            &
        real( cf_bulk( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_ls_rain3d)
      ! Large scale rainfall flux
      environ_fullht_real( :, i ) =                                            &
        real( ls_rain_3d( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_ls_snow3d)
      ! Large scale snowfall flux
      environ_fullht_real( :, i ) =                                            &
        real( ls_snow_3d( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_conv_rain3d)
      ! Convective rainfall flux
      environ_fullht_real( :, i ) =                                            &
        real( conv_rain_3d( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_conv_snow3d)
      ! Convective snowfall flux
      environ_fullht_real( :, i ) =                                            &
        real( conv_snow_3d( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_rho_r2)
      do k = 1, nlayers
        ! Density * r * r
        rho_r2 = real( wetrho_in_w3( map_w3(1) + k - 1 ) *                     &
                       (r_rho_levels( 1, 1, k )**2 ), r_um )
        environ_fullht_real( k, i) = rho_r2
      end do
    case(fldname_liq_cloud_frac)
      ! Liquid cloud fraction
      environ_fullht_real( :, i ) =                                            &
        real( cf_liq( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_autoconv)
      ! Rain autoconversion rate
      environ_fullht_real( :, i ) =                                            &
        real( autoconv( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_accretion)
      ! Rain accretion rate
      environ_fullht_real( :, i ) =                                            &
        real( accretion( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_rim_cry)
      ! Riming rate for ice crystals
      environ_fullht_real( :, i ) =                                            &
        real( rim_cry( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_rim_agg)
      ! Riming rate for ice aggregates
      environ_fullht_real( :, i ) =                                            &
        real( rim_agg( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case(fldname_vertvel)
      ! Vertical velocity
      environ_fullht_real( :, i ) =                                            &
        real( u3_in_wth( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case default
      write( log_scratch_space, '(A,A)' )                                      &
        'Missing required UKCA environment field: ', env_names_fullht_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  ! Drivers in full-height plus level zero grid group (1D fields)

  n_fields = size(env_names_fullht0_real)
  allocate(environ_fullht0_real( 0:nlayers, n_fields ))
  environ_fullht0_real = 0.0_r_um

  do i = 1, n_fields
    select case(env_names_fullht0_real(i))
    case(fldname_interf_z)
      environ_fullht0_real( 0, i ) = 0.0_r_um
      do k = 1, nlayers - 1
        environ_fullht0_real( k, i ) =                                         &
          r_rho_levels( 1, 1, k + 1 ) - r_theta_levels( 1, 1, 0 )
      end do
      environ_fullht0_real( nlayers, i ) =                                     &
        r_theta_levels( 1, 1, nlayers ) - r_theta_levels( 1, 1, 0 )
    case default
      write( log_scratch_space, '(A,A)' )                                      &
        'Missing required UKCA environment field: ', env_names_fullht0_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  ! Drivers in full-height plus one grid group (1D fields)

  n_fields = size(env_names_fullhtp1_real)
  allocate(environ_fullhtp1_real( nlayers + 1, n_fields ))
  environ_fullhtp1_real = 0.0_r_um

  do i = 1, n_fields
    select case(env_names_fullhtp1_real(i))
    case(fldname_exner_rho_lev)
      environ_fullhtp1_real( 1 : nlayers, i ) =                                &
        real( exner_in_w3( map_w3(1) + 0 : map_w3(1) + nlayers - 1 ), r_um )
      ! exner_rho_lev is required above the level at which exner_in_w3 is
      ! defined so is approximated here by extrapolation.
      ! This is a temporary solution pending changes to UKCA to replace the
      ! UM-specific routine for adding emissions to tracers trsrce.
      exner_rho_top = real( exner_in_w3( map_w3(1) + nlayers - 1 ), r_um )
      exner_theta_top = real( exner_in_wth( map_wth(1) + nlayers ), r_um )
      environ_fullhtp1_real( nlayers + 1, i ) =                                &
        exner_rho_top + 2.0_r_um * (exner_theta_top - exner_rho_top)
    case default
      write( log_scratch_space, '(A,A)' )                                      &
        'Missing required UKCA environment field: ', env_names_fullhtp1_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  ! Drivers in boundary layer levels group (1D fields)

  n_fields = size(env_names_bllev_real)
  allocate(environ_bllev_real( bl_levels, n_fields ))
  environ_bllev_real = 0.0_r_um

  do i = 1, n_fields
    select case(env_names_bllev_real(i))

    case(fldname_rhokh_rdz)
      ! Mixing coefficient: eddy diffusivity * rho / dz
      do k = 2, bl_levels
        environ_bllev_real( k, i ) = real( rhokh_bl( map_w3(1) + k - 1 ) *     &
                                           rdz_tq_bl( map_w3(1) + k - 1 ), r_um)
      end do
    case(fldname_dtrdz)
      ! dt / (rho * r * r * dz) for scalar flux divergence
      environ_bllev_real( :, i ) =                                             &
        real( dtrdz_tq_bl( map_wth(1) + 1 : map_wth(1) + bl_levels ), r_um )
    case(fldname_bl_tke)
      ! Turbulent kinetic energy
      environ_bllev_real( 1 : bl_levels - 1, i ) =                             &
        real( tke_bl( map_wth(1) + 1 : map_wth(1) + bl_levels - 1 ), r_um )
    case default
      write( log_scratch_space, '(A,A)' )                                      &
        'Missing required UKCA environment field: ', env_names_bllev_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  ! Drivers in entrainment layers group (1D fields)

  n_fields = size(env_names_entlev_real)
  allocate(environ_entlev_real( nlev_ent_tr_mix, n_fields ))
  environ_entlev_real = 0.0_r_um

  do i = 1, n_fields
    select case(env_names_entlev_real(i))
    case(fldname_we_lim)
      ! rho * entrainment rate at surface mixed layer inversion
      environ_entlev_real( :, i ) =                                            &
        real( ent_we_lim( map_ent(1) + 0 : map_ent(1) + nlev_ent_tr_mix - 1 ), &
              r_um )
    case(fldname_t_frac)
      ! Fraction of time surface mixed layer inversion is above level
      environ_entlev_real( :, i ) =                                            &
        real( ent_t_frac( map_ent(1) + 0 : map_ent(1) + nlev_ent_tr_mix - 1 ), &
              r_um )
    case(fldname_zrzi)
      ! Level height as fraction of surface mixed layer inversion height above
      ! ML base
      environ_entlev_real( :, i ) =                                            &
        real( ent_zrzi( map_ent(1) + 0 : map_ent(1) + nlev_ent_tr_mix - 1 ),   &
              r_um )
    case(fldname_we_lim_dsc)
      ! rho * entrainment rate at decoupled stratocumulus inversion
      environ_entlev_real( :, i ) = real(                                      &
        ent_we_lim_dsc( map_ent(1) + 0 : map_ent(1) + nlev_ent_tr_mix - 1 ),   &
        r_um )
    case(fldname_t_frac_dsc)
      ! Fraction of time decoupled stratocumulus inversion is above level
      environ_entlev_real( :, i ) = real(                                      &
        ent_t_frac_dsc( map_ent(1) + 0 : map_ent(1) + nlev_ent_tr_mix - 1 ),   &
        r_um )
    case(fldname_zrzi_dsc)
      ! Level height as fraction of decoupled stratocumulus inversion height
      ! above DSC ML base
      environ_entlev_real( :, i ) = real(                                      &
        ent_zrzi_dsc( map_ent(1) + 0 : map_ent(1) + nlev_ent_tr_mix - 1 ),     &
        r_um )
    case default
      write( log_scratch_space, '(A,A)' )                                      &
        'Missing required UKCA environment field: ', env_names_entlev_real(i)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  ! Drivers in land point group (1D fields)

  n_fields = size(env_names_land_real)
  allocate(environ_land_real( n_land_pts, n_fields ))
  environ_land_real = 0.0_r_um

  if (l_land) then
    do i = 1, n_fields
      select case(env_names_land_real(i))
      case(fldname_frac_land)
        ! Land_fraction
        environ_land_real( 1, i ) = frac_land
      case(fldname_soil_moisture_layer1)
        ! Soil moisture content of top layer
        environ_land_real( 1, i ) = real( soil_moisture(map_soil(1)) + 0, r_um )
      case default
        write( log_scratch_space, '(A,A)' )                                    &
          'Missing required UKCA environment field: ', env_names_land_real(i)
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end select
    end do
  endif

  ! Drivers in land-point tile groups (2D fields)

  n_fields = size(env_names_landtile_real)
  allocate(environ_landtile_real( n_land_pts, n_land_tile, n_fields ))
  environ_landtile_real = 0.0_r_um

  if (l_land) then
    do i = 1, n_fields
      select case(env_names_landtile_real(i))
      case(fldname_frac_surft)
        do j = 1, n_land_tile
          environ_landtile_real( 1, j, i ) = frac_surft( 1, j )
        end do
      case(fldname_tstar_surft)
        do j = 1, n_land_tile
          environ_landtile_real( 1, j, i ) =                                   &
            real( tile_temperature( map_tile(1) + j - 1 ), r_um )
        end do
      case(fldname_z0_surft)
        do j = 1, n_land_tile
          environ_landtile_real( 1, j, i ) = z0_surft( 1, j )
        end do
      case default
        write( log_scratch_space, '(A,A)' )                                    &
          'Missing required UKCA environment field: ',                         &
          env_names_landtile_real(i)
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end select
    end do
  endif

  n_fields = size(env_names_landtile_logical)
  allocate(environ_landtile_logical( n_land_pts, n_land_tile, n_fields) )
  environ_landtile_logical = .false.

  if (l_land) then
    do i = 1, n_fields
      select case(env_names_landtile_logical(i))
      case(fldname_l_active_surft)
        l_tile_active( 1, : ) = .false.
        do j = 1, n_land_tile
          do k = 1, ainfo%surft_pts(j)   ! i.e 0 or 1 times
            l_tile_active( 1, j ) = .true.
          end do
          environ_landtile_logical( 1, j, i ) = l_tile_active( 1, j )
        end do
      case default
        write( log_scratch_space, '(A,A)' )                                    &
          'Missing required UKCA environment field: ',                         &
          env_names_landtile_logical(i)
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end select
    end do
  endif

  call ancil_info_nullify(ainfo)
  call ancil_info_dealloc(ainfo_data)

  ! Drivers in land-point plant functional type tile group (2D fields)

  n_fields = size(env_names_landpft_real)
  allocate(environ_landpft_real( n_land_pts, npft, n_fields ))
  environ_landpft_real = 0.0_r_um

  if (l_land) then
    do i = 1, n_fields
      select case(env_names_landpft_real(i))
      case(fldname_lai_pft)
        ! Leaf area index
        do j = 1, npft
          environ_landpft_real( 1, j, i ) = lai_pft( 1, j )
        end do
      case(fldname_canht_pft)
        ! Canopy height
        do j = 1, npft
          environ_landpft_real( 1, j, i ) = canht_pft( 1, j )
        end do
      case default
        write( log_scratch_space, '(A,A)' )                                    &
          'Missing required UKCA environment field: ', env_names_landpft_real(i)
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end select
    end do
  endif

  ! Emissions in flat group

  n_fields = size(emiss_names_flat)
  allocate(emissions_flat( n_fields ))
  emissions_flat = 0.0_r_um

  do i = 1, size(emiss_names_flat)
    select case(emiss_names_flat(i))
    case('emissions_C2H6')
      emissions_flat(i) = real( emiss_c2h6(map_2d(1)), r_um )
    case('emissions_C3H8')
      emissions_flat(i) = real( emiss_c3h8(map_2d(1)), r_um )
    case('emissions_C5H8')
      emissions_flat(i) = real( emiss_c5h8(map_2d(1)), r_um )
    case('emissions_CH4')
      emissions_flat(i) = real( emiss_ch4(map_2d(1)), r_um )
    case('emissions_CO')
      emissions_flat(i) = real( emiss_co(map_2d(1)), r_um )
    case('emissions_HCHO')
      emissions_flat(i) = real( emiss_hcho(map_2d(1)), r_um )
    case('emissions_Me2CO')
      emissions_flat(i) = real( emiss_me2co(map_2d(1)), r_um )
    case('emissions_MeCHO')
      emissions_flat(i) = real( emiss_mecho(map_2d(1)), r_um )
    case('emissions_MeOH')
      emissions_flat(i) = real( emiss_meoh(map_2d(1)), r_um )
    case('emissions_NH3')
      emissions_flat(i) = real( emiss_nh3(map_2d(1)), r_um )
    case('emissions_NO')
      emissions_flat(i) = real( emiss_no(map_2d(1)), r_um )
    case('emissions_BC_biofuel')
      emissions_flat(i) = real( emiss_bc_biofuel(map_2d(1)), r_um )
    case('emissions_BC_fossil')
      emissions_flat(i) = real( emiss_bc_fossil(map_2d(1)), r_um )
    case('emissions_BC_biomass_high')
      emissions_flat(i) = real( emiss_bc_biomass_high(map_2d(1)), r_um )
    case('emissions_BC_biomass_low')
      emissions_flat(i) = real( emiss_bc_biomass_low(map_2d(1)), r_um )
    case('emissions_DMS')
      emissions_flat(i) = real( emiss_dms_land(map_2d(1)), r_um )
    case('emissions_Monoterp')
      emissions_flat(i) = real( emiss_monoterp(map_2d(1)), r_um )
    case('emissions_OM_biofuel')
      emissions_flat(i) = real( emiss_om_biofuel(map_2d(1)), r_um )
    case('emissions_OM_biomass_high')
      emissions_flat(i) = real( emiss_om_biomass_high(map_2d(1)), r_um )
    case('emissions_OM_biomass_low')
      emissions_flat(i) = real( emiss_om_biomass_low(map_2d(1)), r_um )
    case('emissions_OM_fossil')
      emissions_flat(i) = real( emiss_om_fossil(map_2d(1)), r_um )
    case('emissions_SO2_high')
      emissions_flat(i) = real( emiss_so2_high(map_2d(1)), r_um )
    case('emissions_SO2_low')
      emissions_flat(i) = real( emiss_so2_low(map_2d(1)), r_um )
    case default
      write( log_scratch_space, '(A,A)' )                                      &
        'Missing required UKCA emission field: ', emiss_names_flat(i)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  n_fields = size(emiss_names_fullht)
  allocate(emissions_fullht( nlayers, n_fields ))
  emissions_fullht = 0.0_r_um

  ! Emissions in full height group

  do i = 1, size(emiss_names_fullht)
    select case(emiss_names_fullht(i))
    case('emissions_NO_aircrft')
      emissions_fullht( :, i ) =                                               &
        real( emiss_no_aircrft( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case('emissions_BC_biomass')
      emissions_fullht( :, i ) =                                               &
        real( emiss_bc_biomass( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case('emissions_OM_biomass')
      emissions_fullht( :, i ) =                                               &
        real( emiss_om_biomass( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case('emissions_SO2_nat')
      emissions_fullht( :, i ) =                                               &
        real( emiss_so2_nat( map_wth(1) + 1 : map_wth(1) + nlayers ), r_um )
    case default
      write( log_scratch_space, '(A,A)' )                                      &
        'Missing required UKCA emission field: ', emiss_names_fullht(i)
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )
    end select
  end do

  !-----------------------------------------------------------------------
  ! Do UKCA time step
  !-----------------------------------------------------------------------

  ! Collate time variables
  current_time(1) = int( current_time_year, i_um )
  current_time(2) = int( current_time_month, i_um )
  current_time(3) = int( current_time_day, i_um )
  current_time(4) = int( current_time_hour, i_um )
  current_time(5) = int( current_time_minute, i_um )
  current_time(6) = int( current_time_second, i_um )
  current_time(7) = int( current_time_daynum, i_um )
  previous_time(1) = int( previous_time_year, i_um )
  previous_time(2) = int( previous_time_month, i_um )
  previous_time(3) = int( previous_time_day, i_um )
  previous_time(4) = int( previous_time_hour, i_um )
  previous_time(5) = int( previous_time_minute, i_um )
  previous_time(6) = int( previous_time_second, i_um )
  previous_time(7) = int( previous_time_daynum, i_um )

  call ukca_step_control( int( timestep_number, i_um ), current_time,          &
                          tracer, ntp, r_theta_levels, r_rho_levels,           &
                          error_code, previous_time=previous_time,             &
                          envgroup_scalar_real=environ_scalar_real,            &
                          envgroup_flat_integer=environ_flat_integer,          &
                          envgroup_flat_real=environ_flat_real,                &
                          envgroup_flat_logical=environ_flat_logical,          &
                          envgroup_flatpft_real=environ_flatpft_real,          &
                          envgroup_fullht_real=environ_fullht_real,            &
                          envgroup_fullht0_real=environ_fullht0_real,          &
                          envgroup_fullhtp1_real=environ_fullhtp1_real,        &
                          envgroup_bllev_real=environ_bllev_real,              &
                          envgroup_entlev_real=environ_entlev_real,            &
                          envgroup_land_real=environ_land_real,                &
                          envgroup_landtile_real=environ_landtile_real,        &
                          envgroup_landtile_logical=environ_landtile_logical,  &
                          envgroup_landpft_real=environ_landpft_real,          &
                          emissions_flat=emissions_flat,                       &
                          emissions_fullht=emissions_fullht,                   &
                          error_message=ukca_errmsg,                           &
                          error_routine=ukca_errproc )

  if (error_code > 0) then
    write( log_scratch_space, '(A)' )                                          &
      trim(ukca_errmsg) // ' in UKCA procedure ' // trim(ukca_errproc)
    call log_event( log_scratch_space, LOG_LEVEL_ERROR )
  end if

  ! Clear emissions fields
  deallocate(emissions_fullht)
  deallocate(emissions_flat)

  ! Clear environmental driver fields
  deallocate(environ_landpft_real)
  deallocate(environ_landtile_logical)
  deallocate(environ_landtile_real)
  deallocate(environ_land_real)
  deallocate(environ_entlev_real)
  deallocate(environ_bllev_real)
  deallocate(environ_fullhtp1_real)
  deallocate(environ_fullht0_real)
  deallocate(environ_fullht_real)
  deallocate(environ_flatpft_real)
  deallocate(environ_flat_logical)
  deallocate(environ_flat_real)
  deallocate(environ_flat_integer)
  deallocate(environ_scalar_real)

  !-----------------------------------------------------------------------
  ! Copy updated UKCA arrays back to LFRic fields
  !-----------------------------------------------------------------------

  ! Tracers

  do i = 1, size(tracer_names)
    select case(tracer_names(i))
    case(fldname_o3p)
      o3p( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( tracer( :, i ), r_def )
      o3p( map_wth(1) + 0 ) = o3p( map_wth(1) + 1 )
    case(fldname_o3)
      o3( map_wth(1) + 1 : map_wth(1) + nlayers ) =                            &
        real( tracer( :, i ), r_def )
      o3( map_wth(1) + 0 ) = o3( map_wth(1) + 1 )
    case(fldname_n)
      n( map_wth(1) + 1 : map_wth(1) + nlayers ) =                             &
        real( tracer( :, i ), r_def )
      n( map_wth(1) + 0 ) = n( map_wth(1) + 1 )
    case(fldname_no)
      no( map_wth(1) + 1 : map_wth(1) + nlayers ) =                            &
        real( tracer( :, i ), r_def )
      no( map_wth(1) + 0 ) = no( map_wth(1) + 1 )
    case(fldname_no3)
      no3( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( tracer( :, i ), r_def )
      no3( map_wth(1) + 0 ) = no3( map_wth(1) + 1 )
    case(fldname_lumped_n)
      lumped_n( map_wth(1) + 1 : map_wth(1) + nlayers ) =                      &
        real( tracer( :, i ), r_def )
      lumped_n( map_wth(1) + 0 ) = lumped_n( map_wth(1) + 1 )
    case(fldname_n2o5)
      n2o5( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      n2o5( map_wth(1) + 0 ) = n2o5( map_wth(1) + 1 )
    case(fldname_ho2no2)
      ho2no2( map_wth(1) + 1 : map_wth(1) + nlayers ) =                        &
        real( tracer( :, i ), r_def )
      ho2no2( map_wth(1) + 0 ) = ho2no2( map_wth(1) + 1 )
    case(fldname_hono2)
      hono2( map_wth(1) + 1 : map_wth(1) + nlayers ) =                         &
        real( tracer( :, i ), r_def )
      hono2( map_wth(1) + 0 ) = hono2( map_wth(1) + 1 )
    case(fldname_h2o2)
      h2o2( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      h2o2( map_wth(1) + 0 ) = h2o2( map_wth(1) + 1 )
    case(fldname_ch4)
      ch4( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( tracer( :, i ), r_def )
      ch4( map_wth(1) + 0 ) = ch4( map_wth(1) + 1 )
    case(fldname_co)
      co( map_wth(1) + 1 : map_wth(1) + nlayers ) =                            &
        real( tracer( :, i ), r_def )
      co( map_wth(1) + 0 ) = co( map_wth(1) + 1 )
    case(fldname_hcho)
      hcho( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      hcho( map_wth(1) + 0 ) = hcho( map_wth(1) + 1 )
    case(fldname_meoo)
      meoo( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      meoo( map_wth(1) + 0 ) = meoo( map_wth(1) + 1 )
    case(fldname_meooh)
      meooh( map_wth(1) + 1 : map_wth(1) + nlayers ) =                         &
        real( tracer( :, i ), r_def )
      meooh( map_wth(1) + 0 ) = meooh( map_wth(1) + 1 )
    case(fldname_h)
      h( map_wth(1) + 1 : map_wth(1) + nlayers ) =                             &
        real( tracer( :, i ), r_def )
      h( map_wth(1) + 0 ) = h( map_wth(1) + 1 )
    case(fldname_ch2o)
      ! DO nothing - UKCA water tracer not directly fed back to Q/m_v_n
    case(fldname_oh)
      oh( map_wth(1) + 1 : map_wth(1) + nlayers ) =                            &
        real( tracer( :, i ), r_def )
      oh( map_wth(1) + 0 ) = oh( map_wth(1) + 1 )
    case(fldname_ho2)
      ho2( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( tracer( :, i ), r_def )
      ho2( map_wth(1) + 0 ) = ho2( map_wth(1) + 1 )
    case(fldname_cl)
      cl( map_wth(1) + 1 : map_wth(1) + nlayers ) =                            &
        real( tracer( :, i ), r_def )
      cl( map_wth(1) + 0 ) = cl( map_wth(1) + 1 )
    case(fldname_cl2o2)
      cl2o2( map_wth(1) + 1 : map_wth(1) + nlayers ) =                         &
        real( tracer( :, i ), r_def )
      cl2o2( map_wth(1) + 0 ) = cl2o2( map_wth(1) + 1 )
    case(fldname_clo)
      clo( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( tracer( :, i ), r_def )
      clo( map_wth(1) + 0 ) = clo( map_wth(1) + 1 )
    case(fldname_oclo)
      oclo( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      oclo( map_wth(1) + 0 ) = oclo( map_wth(1) + 1 )
    case(fldname_br)
      br( map_wth(1) + 1 : map_wth(1) + nlayers ) =                            &
        real( tracer( :, i ), r_def )
      br( map_wth(1) + 0 ) = br( map_wth(1) + 1 )
    case(fldname_lumped_br)
      lumped_br( map_wth(1) + 1 : map_wth(1) + nlayers ) =                     &
        real( tracer( :, i ), r_def )
      lumped_br( map_wth(1) + 0 ) = lumped_br( map_wth(1) + 1 )
    case(fldname_brcl)
      brcl( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      brcl( map_wth(1) + 0 ) = brcl( map_wth(1) + 1 )
    case(fldname_brono2)
      brono2( map_wth(1) + 1 : map_wth(1) + nlayers ) =                        &
        real( tracer( :, i ), r_def )
      brono2( map_wth(1) + 0 ) = brono2( map_wth(1) + 1 )
    case(fldname_n2o)
      n2o( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( tracer( :, i ), r_def )
      n2o( map_wth(1) + 0 ) = n2o( map_wth(1) + 1 )
    case(fldname_lumped_cl)
      lumped_cl( map_wth(1) + 1 : map_wth(1) + nlayers ) =                     &
        real( tracer( :, i ), r_def )
      lumped_cl( map_wth(1) + 0 ) = lumped_cl( map_wth(1) + 1 )
    case(fldname_hocl)
      hocl( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      hocl( map_wth(1) + 0 ) = hocl( map_wth(1) + 1 )
    case(fldname_hbr)
      hbr( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( tracer( :, i ), r_def )
      hbr( map_wth(1) + 0 ) = hbr( map_wth(1) + 1 )
    case(fldname_hobr)
      hobr( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      hobr( map_wth(1) + 0 ) = hobr( map_wth(1) + 1 )
    case(fldname_clono2)
      clono2( map_wth(1) + 1 : map_wth(1) + nlayers ) =                        &
        real( tracer( :, i ), r_def )
      clono2( map_wth(1) + 0 ) = clono2( map_wth(1) + 1 )
    case(fldname_cfcl3)
      cfcl3( map_wth(1) + 1 : map_wth(1) + nlayers ) =                         &
        real( tracer( :, i ), r_def )
      cfcl3( map_wth(1) + 0 ) = cfcl3( map_wth(1) + 1 )
    case(fldname_cf2cl2)
      cf2cl2( map_wth(1) + 1 : map_wth(1) + nlayers ) =                        &
        real( tracer( :, i ), r_def )
      cf2cl2( map_wth(1) + 0 ) = cf2cl2( map_wth(1) + 1 )
    case(fldname_mebr)
      mebr( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      mebr( map_wth(1) + 0 ) = mebr( map_wth(1) + 1 )
    case(fldname_hono)
      hono( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      hono( map_wth(1) + 0 ) = hono( map_wth(1) + 1 )
    case(fldname_c2h6)
      c2h6( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      c2h6( map_wth(1) + 0 ) = c2h6( map_wth(1) + 1 )
    case(fldname_etoo)
      etoo( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      etoo( map_wth(1) + 0 ) = etoo( map_wth(1) + 1 )
    case(fldname_etooh)
      etooh( map_wth(1) + 1 : map_wth(1) + nlayers ) =                         &
        real( tracer( :, i ), r_def )
      etooh( map_wth(1) + 0 ) = etooh( map_wth(1) + 1 )
    case(fldname_mecho)
      mecho( map_wth(1) + 1 : map_wth(1) + nlayers ) =                         &
        real( tracer( :, i ), r_def )
      mecho( map_wth(1) + 0 ) = mecho( map_wth(1) + 1 )
    case(fldname_meco3)
      meco3( map_wth(1) + 1 : map_wth(1) + nlayers ) =                         &
        real( tracer( :, i ), r_def )
      meco3( map_wth(1) + 0 ) = meco3( map_wth(1) + 1 )
    case(fldname_pan)
      pan( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( tracer( :, i ), r_def )
      pan( map_wth(1) + 0 ) = pan( map_wth(1) + 1 )
    case(fldname_c3h8)
      c3h8( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      c3h8( map_wth(1) + 0 ) = c3h8( map_wth(1) + 1 )
    case(fldname_n_proo)
      n_proo( map_wth(1) + 1 : map_wth(1) + nlayers ) =                        &
        real( tracer( :, i ), r_def )
      n_proo( map_wth(1) + 0 ) = n_proo( map_wth(1) + 1 )
    case(fldname_i_proo)
      i_proo( map_wth(1) + 1 : map_wth(1) + nlayers ) =                        &
        real( tracer( :, i ), r_def )
      i_proo( map_wth(1) + 0 ) = i_proo( map_wth(1) + 1 )
    case(fldname_n_prooh)
      n_prooh( map_wth(1) + 1 : map_wth(1) + nlayers ) =                       &
        real( tracer( :, i ), r_def )
      n_prooh( map_wth(1) + 0 ) = n_prooh( map_wth(1) + 1 )
    case(fldname_i_prooh)
      i_prooh( map_wth(1) + 1 : map_wth(1) + nlayers ) =                       &
        real( tracer( :, i ), r_def )
      i_prooh( map_wth(1) + 0 ) = i_prooh( map_wth(1) + 1 )
    case(fldname_etcho)
      etcho( map_wth(1) + 1 : map_wth(1) + nlayers ) =                         &
        real( tracer( :, i ), r_def )
      etcho( map_wth(1) + 0 ) = etcho( map_wth(1) + 1 )
    case(fldname_etco3)
      etco3( map_wth(1) + 1 : map_wth(1) + nlayers ) =                         &
        real( tracer( :, i ), r_def )
      etco3( map_wth(1) + 0 ) = etco3( map_wth(1) + 1 )
    case(fldname_me2co)
      me2co( map_wth(1) + 1 : map_wth(1) + nlayers ) =                         &
        real( tracer( :, i ), r_def )
      me2co( map_wth(1) + 0 ) = me2co( map_wth(1) + 1 )
    case(fldname_mecoch2oo)
      mecoch2oo( map_wth(1) + 1 : map_wth(1) + nlayers ) =                     &
        real( tracer( :, i ), r_def )
      mecoch2oo( map_wth(1) + 0 ) = mecoch2oo( map_wth(1) + 1 )
    case(fldname_mecoch2ooh)
      mecoch2ooh( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      mecoch2ooh( map_wth(1) + 0 ) = mecoch2ooh( map_wth(1) + 1 )
    case(fldname_ppan)
      ppan( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      ppan( map_wth(1) + 0 ) = ppan( map_wth(1) + 1 )
    case(fldname_meono2)
      meono2( map_wth(1) + 1 : map_wth(1) + nlayers ) =                        &
        real( tracer( :, i ), r_def )
      meono2( map_wth(1) + 0 ) = meono2( map_wth(1) + 1 )
    case(fldname_c5h8)
      c5h8( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      c5h8( map_wth(1) + 0 ) = c5h8( map_wth(1) + 1 )
    case(fldname_iso2)
      iso2( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      iso2( map_wth(1) + 0 ) = iso2( map_wth(1) + 1 )
    case(fldname_isooh)
      isooh( map_wth(1) + 1 : map_wth(1) + nlayers ) =                         &
        real( tracer( :, i ), r_def )
      isooh( map_wth(1) + 0 ) = isooh( map_wth(1) + 1 )
    case(fldname_ison)
      ison( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      ison( map_wth(1) + 0 ) = ison( map_wth(1) + 1 )
    case(fldname_macr)
      macr( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      macr( map_wth(1) + 0 ) = macr( map_wth(1) + 1 )
    case(fldname_macrooh)
      macrooh( map_wth(1) + 1 : map_wth(1) + nlayers ) =                       &
        real( tracer( :, i ), r_def )
      macrooh( map_wth(1) + 0 ) = macrooh( map_wth(1) + 1 )
    case(fldname_macro2)
      macro2( map_wth(1) + 1 : map_wth(1) + nlayers ) =                        &
        real( tracer( :, i ), r_def )
      macro2( map_wth(1) + 0 ) = macro2( map_wth(1) + 1 )
    case(fldname_mpan)
      mpan( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      mpan( map_wth(1) + 0 ) = mpan( map_wth(1) + 1 )
    case(fldname_hacet)
      hacet( map_wth(1) + 1 : map_wth(1) + nlayers ) =                         &
        real( tracer( :, i ), r_def )
      hacet( map_wth(1) + 0 ) = hacet( map_wth(1) + 1 )
    case(fldname_mgly)
      mgly( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      mgly( map_wth(1) + 0 ) = mgly( map_wth(1) + 1 )
    case(fldname_nald)
      nald( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      nald( map_wth(1) + 0 ) = nald( map_wth(1) + 1 )
    case(fldname_hcooh)
      hcooh( map_wth(1) + 1 : map_wth(1) + nlayers ) =                         &
        real( tracer( :, i ), r_def )
      hcooh( map_wth(1) + 0 ) = hcooh( map_wth(1) + 1 )
    case(fldname_meco3h)
      meco3h( map_wth(1) + 1 : map_wth(1) + nlayers ) =                        &
        real( tracer( :, i ), r_def )
      meco3h( map_wth(1) + 0 ) = meco3h( map_wth(1) + 1 )
    case(fldname_meco2h)
      meco2h( map_wth(1) + 1 : map_wth(1) + nlayers ) =                        &
        real( tracer( :, i ), r_def )
      meco2h( map_wth(1) + 0 ) = meco2h( map_wth(1) + 1 )
    case(fldname_h2)
      h2( map_wth(1) + 1 : map_wth(1) + nlayers ) =                            &
        real( tracer( :, i ), r_def )
      h2( map_wth(1) + 0 ) = h2( map_wth(1) + 1 )
    case(fldname_meoh)
      meoh( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      meoh( map_wth(1) + 0 ) = meoh( map_wth(1) + 1 )
    case(fldname_msa)
      msa( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( tracer( :, i ), r_def )
      msa( map_wth(1) + 0 ) = msa( map_wth(1) + 1 )
    case(fldname_nh3)
      nh3( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( tracer( :, i ), r_def )
      nh3( map_wth(1) + 0 ) = nh3( map_wth(1) + 1 )
    case(fldname_cs2)
      cs2( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( tracer( :, i ), r_def )
      cs2( map_wth(1) + 0 ) = cs2( map_wth(1) + 1 )
    case(fldname_csul)
      csul( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( tracer( :, i ), r_def )
      csul( map_wth(1) + 0 ) = csul( map_wth(1) + 1 )
    case(fldname_h2s)
      h2s( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( tracer( :, i ), r_def )
      h2s( map_wth(1) + 0 ) = h2s( map_wth(1) + 1 )
    case(fldname_so3)
      so3( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( tracer( :, i ), r_def )
      so3( map_wth(1) + 0 ) = so3( map_wth(1) + 1 )
    case(fldname_passive_o3)
      passive_o3( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      passive_o3( map_wth(1) + 0 ) = passive_o3( map_wth(1) + 1 )
    case(fldname_age_of_air)
      age_of_air( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      age_of_air( map_wth(1) + 0 ) = age_of_air( map_wth(1) + 1 )
    case(fldname_dms)
      dms( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( tracer( :, i ), r_def )
      dms( map_wth(1) + 0 ) = dms( map_wth(1) + 1 )
    case(fldname_so2)
      so2( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( tracer( :, i ), r_def )
      so2( map_wth(1) + 0 ) = so2( map_wth(1) + 1 )
    case(fldname_h2so4)
      h2so4( map_wth(1) + 1 : map_wth(1) + nlayers ) =                         &
        real( tracer( :, i ), r_def )
      h2so4( map_wth(1) + 0 ) = h2so4( map_wth(1) + 1 )
    case(fldname_dmso)
      dmso( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( tracer( :, i ), r_def )
      dmso( map_wth(1) + 0 ) = dmso( map_wth(1) + 1 )
    case(fldname_monoterpene)
      monoterpene( map_wth(1) + 1 : map_wth(1) + nlayers ) =                   &
        real( tracer( :, i ), r_def )
      monoterpene( map_wth(1) + 0 ) = monoterpene( map_wth(1) + 1 )
    case(fldname_secondary_organic)
      secondary_organic( map_wth(1) + 1 : map_wth(1) + nlayers ) =             &
        real( tracer( :, i ), r_def )
      secondary_organic( map_wth(1) + 0 ) = secondary_organic( map_wth(1) + 1 )
    case(fldname_n_nuc_sol)
      n_nuc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =                     &
        real( tracer( :, i ), r_def )
      n_nuc_sol( map_wth(1) + 0 ) = n_nuc_sol( map_wth(1) + 1 )
    case(fldname_nuc_sol_su)
      nuc_sol_su( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      nuc_sol_su( map_wth(1) + 0 ) = nuc_sol_su( map_wth(1) + 1 )
    case(fldname_nuc_sol_om)
      nuc_sol_om( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      nuc_sol_om( map_wth(1) + 0 ) = nuc_sol_om( map_wth(1) + 1 )
    case(fldname_n_ait_sol)
      n_ait_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =                     &
        real( tracer( :, i ), r_def )
      n_ait_sol( map_wth(1) + 0 ) = n_ait_sol( map_wth(1) + 1 )
    case(fldname_ait_sol_su)
      ait_sol_su( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      ait_sol_su( map_wth(1) + 0 ) = ait_sol_su( map_wth(1) + 1 )
    case(fldname_ait_sol_bc)
      ait_sol_bc( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      ait_sol_bc( map_wth(1) + 0 ) = ait_sol_bc( map_wth(1) + 1 )
    case(fldname_ait_sol_om)
      ait_sol_om( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      ait_sol_om( map_wth(1) + 0 ) = ait_sol_om( map_wth(1) + 1 )
    case(fldname_n_acc_sol)
      n_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =                     &
        real( tracer( :, i ), r_def )
      n_acc_sol( map_wth(1) + 0 ) = n_acc_sol( map_wth(1) + 1 )
    case(fldname_acc_sol_su)
      acc_sol_su( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      acc_sol_su( map_wth(1) + 0 ) = acc_sol_su( map_wth(1) + 1 )
    case(fldname_acc_sol_bc)
      acc_sol_bc( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      acc_sol_bc( map_wth(1) + 0 ) = acc_sol_bc( map_wth(1) + 1 )
    case(fldname_acc_sol_om)
      acc_sol_om( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      acc_sol_om( map_wth(1) + 0 ) = acc_sol_om( map_wth(1) + 1 )
    case(fldname_acc_sol_ss)
      acc_sol_ss( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      acc_sol_ss( map_wth(1) + 0 ) = acc_sol_ss( map_wth(1) + 1 )
    case(fldname_acc_sol_du)
      acc_sol_du( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      acc_sol_du( map_wth(1) + 0 ) = acc_sol_du( map_wth(1) + 1 )
    case(fldname_n_cor_sol)
      n_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =                     &
        real( tracer( :, i ), r_def )
      n_cor_sol( map_wth(1) + 0 ) = n_cor_sol( map_wth(1) + 1 )
    case(fldname_cor_sol_su)
      cor_sol_su( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      cor_sol_su( map_wth(1) + 0 ) = cor_sol_su( map_wth(1) + 1 )
    case(fldname_cor_sol_bc)
      cor_sol_bc( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      cor_sol_bc( map_wth(1) + 0 ) = cor_sol_bc( map_wth(1) + 1 )
    case(fldname_cor_sol_om)
      cor_sol_om( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      cor_sol_om( map_wth(1) + 0 ) = cor_sol_om( map_wth(1) + 1 )
    case(fldname_cor_sol_ss)
      cor_sol_ss( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      cor_sol_ss( map_wth(1) + 0 ) = cor_sol_ss( map_wth(1) + 1 )
    case(fldname_cor_sol_du)
      cor_sol_du( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      cor_sol_du( map_wth(1) + 0 ) = cor_sol_du( map_wth(1) + 1 )
    case(fldname_n_ait_ins)
      n_ait_ins( map_wth(1) + 1 : map_wth(1) + nlayers ) =                     &
        real( tracer( :, i ), r_def )
      n_ait_ins( map_wth(1) + 0 ) = n_ait_ins( map_wth(1) + 1 )
    case(fldname_ait_ins_bc)
      ait_ins_bc( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      ait_ins_bc( map_wth(1) + 0 ) = ait_ins_bc( map_wth(1) + 1 )
    case(fldname_ait_ins_om)
      ait_ins_om( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      ait_ins_om( map_wth(1) + 0 ) = ait_ins_om( map_wth(1) + 1 )
    case(fldname_n_acc_ins)
      n_acc_ins( map_wth(1) + 1 : map_wth(1) + nlayers ) =                     &
        real( tracer( :, i ), r_def )
      n_acc_ins( map_wth(1) + 0 ) = n_acc_ins( map_wth(1) + 1 )
    case(fldname_acc_ins_du)
      acc_ins_du( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      acc_ins_du( map_wth(1) + 0 ) = acc_ins_du( map_wth(1) + 1 )
    case(fldname_n_cor_ins)
      n_cor_ins( map_wth(1) + 1 : map_wth(1) + nlayers ) =                     &
        real( tracer( :, i ), r_def )
      n_cor_ins( map_wth(1) + 0 ) = n_cor_ins( map_wth(1) + 1 )
    case(fldname_cor_ins_du)
      cor_ins_du( map_wth(1) + 1 : map_wth(1) + nlayers ) =                    &
        real( tracer( :, i ), r_def )
      cor_ins_du( map_wth(1) + 0 ) = cor_ins_du( map_wth(1) + 1 )
    end select
  end do

  ! Non-transported prognostics

  do i = 1, size(ntp_names)
    select case(ntp_names(i))
    case(fldname_cloud_drop_no_conc)
      ! Impose lower limit on CDNC (5 cm-3)
      cloud_drop_no_conc( map_wth(1) + 1 : map_wth(1) + nlayers ) =            &
        real( max( 5.0e6_r_um, ntp( :, i ) ), r_def )
      cloud_drop_no_conc( map_wth(1) + 0 ) =                                   &
        cloud_drop_no_conc( map_wth(1) + 1 )
    case(fldname_surfarea)
       ! UKCA output not currently used
    case(fldname_drydp_ait_sol)
      drydp_ait_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      drydp_ait_sol( map_wth(1) + 0 ) = drydp_ait_sol( map_wth(1) + 1 )
    case(fldname_drydp_acc_sol)
      drydp_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      drydp_acc_sol( map_wth(1) + 0 ) = drydp_acc_sol( map_wth(1) + 1 )
    case(fldname_drydp_cor_sol)
      drydp_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      drydp_cor_sol( map_wth(1) + 0 ) = drydp_cor_sol( map_wth(1) + 1 )
    case(fldname_drydp_ait_ins)
      drydp_ait_ins( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      drydp_ait_ins( map_wth(1) + 0 ) = drydp_ait_ins( map_wth(1) + 1 )
    case(fldname_drydp_acc_ins)
      drydp_acc_ins( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      drydp_acc_ins( map_wth(1) + 0 ) = drydp_acc_ins( map_wth(1) + 1 )
    case(fldname_drydp_cor_ins)
      drydp_cor_ins( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      drydp_cor_ins( map_wth(1) + 0 ) = drydp_cor_ins( map_wth(1) + 1 )
    case(fldname_wetdp_ait_sol)
      wetdp_ait_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      wetdp_ait_sol( map_wth(1) + 0 ) = wetdp_ait_sol( map_wth(1) + 1 )
    case(fldname_wetdp_acc_sol)
      wetdp_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      wetdp_acc_sol( map_wth(1) + 0 ) = wetdp_acc_sol( map_wth(1) + 1 )
    case(fldname_wetdp_cor_sol)
      wetdp_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      wetdp_cor_sol( map_wth(1) + 0 ) = wetdp_cor_sol( map_wth(1) + 1 )
    case(fldname_rhopar_ait_sol)
      rhopar_ait_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =                &
        real( ntp( :, i ), r_def )
      rhopar_ait_sol( map_wth(1) + 0 ) = rhopar_ait_sol( map_wth(1) + 1 )
    case(fldname_rhopar_acc_sol)
      rhopar_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =                &
        real( ntp( :, i ), r_def )
      rhopar_acc_sol( map_wth(1) + 0 ) = rhopar_acc_sol( map_wth(1) + 1 )
    case(fldname_rhopar_cor_sol)
      rhopar_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =                &
        real( ntp( :, i ), r_def )
      rhopar_cor_sol( map_wth(1) + 0 ) = rhopar_cor_sol( map_wth(1) + 1 )
    case(fldname_rhopar_ait_ins)
      rhopar_ait_ins( map_wth(1) + 1 : map_wth(1) + nlayers ) =                &
        real( ntp( :, i ), r_def )
      rhopar_ait_ins( map_wth(1) + 0 ) = rhopar_ait_ins( map_wth(1) + 1 )
    case(fldname_rhopar_acc_ins)
      rhopar_acc_ins( map_wth(1) + 1 : map_wth(1) + nlayers ) =                &
        real( ntp( :, i ), r_def )
      rhopar_acc_ins( map_wth(1) + 0 ) = rhopar_acc_ins( map_wth(1) + 1 )
    case(fldname_rhopar_cor_ins)
      rhopar_cor_ins( map_wth(1) + 1 : map_wth(1) + nlayers ) =                &
        real( ntp( :, i ), r_def )
      rhopar_cor_ins( map_wth(1) + 0 ) = rhopar_cor_ins( map_wth(1) + 1 )
    case(fldname_pvol_su_ait_sol)
      pvol_su_ait_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      pvol_su_ait_sol( map_wth(1) + 0 ) = pvol_su_ait_sol( map_wth(1) + 1 )
    case(fldname_pvol_bc_ait_sol)
      pvol_bc_ait_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      pvol_bc_ait_sol( map_wth(1) + 0 ) = pvol_bc_ait_sol( map_wth(1) + 1 )
    case(fldname_pvol_om_ait_sol)
      pvol_om_ait_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      pvol_om_ait_sol( map_wth(1) + 0 ) = pvol_om_ait_sol( map_wth(1) + 1 )
    case(fldname_pvol_wat_ait_sol)
      pvol_wat_ait_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =              &
        real( ntp( :, i ), r_def )
      pvol_wat_ait_sol( map_wth(1) + 0 ) = pvol_wat_ait_sol( map_wth(1) + 1 )
    case(fldname_pvol_su_acc_sol)
      pvol_su_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      pvol_su_acc_sol( map_wth(1) + 0 ) = pvol_su_acc_sol( map_wth(1) + 1 )
    case(fldname_pvol_bc_acc_sol)
      pvol_bc_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      pvol_bc_acc_sol( map_wth(1) + 0 ) = pvol_bc_acc_sol( map_wth(1) + 1 )
    case(fldname_pvol_om_acc_sol)
      pvol_om_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      pvol_om_acc_sol( map_wth(1) + 0 ) = pvol_om_acc_sol( map_wth(1) + 1 )
    case(fldname_pvol_ss_acc_sol)
      pvol_ss_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      pvol_ss_acc_sol( map_wth(1) + 0 ) = pvol_ss_acc_sol( map_wth(1) + 1 )
    case(fldname_pvol_du_acc_sol)
      pvol_du_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      pvol_du_acc_sol( map_wth(1) + 0 ) = pvol_du_acc_sol( map_wth(1) + 1 )
    case(fldname_pvol_wat_acc_sol)
      pvol_wat_acc_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =              &
        real( ntp( :, i ), r_def )
      pvol_wat_acc_sol( map_wth(1) + 0 ) = pvol_wat_acc_sol( map_wth(1) + 1 )
    case(fldname_pvol_su_cor_sol)
      pvol_su_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      pvol_su_cor_sol( map_wth(1) + 0 ) = pvol_su_cor_sol( map_wth(1) + 1 )
    case(fldname_pvol_bc_cor_sol)
      pvol_bc_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      pvol_bc_cor_sol( map_wth(1) + 0 ) = pvol_bc_cor_sol( map_wth(1) + 1 )
    case(fldname_pvol_om_cor_sol)
      pvol_om_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      pvol_om_cor_sol( map_wth(1) + 0 ) = pvol_om_cor_sol( map_wth(1) + 1 )
    case(fldname_pvol_ss_cor_sol)
      pvol_ss_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      pvol_ss_cor_sol( map_wth(1) + 0 ) = pvol_ss_cor_sol( map_wth(1) + 1 )
    case(fldname_pvol_du_cor_sol)
      pvol_du_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      pvol_du_cor_sol( map_wth(1) + 0 ) = pvol_du_cor_sol( map_wth(1) + 1 )
    case(fldname_pvol_wat_cor_sol)
      pvol_wat_cor_sol( map_wth(1) + 1 : map_wth(1) + nlayers ) =              &
        real( ntp( :, i ), r_def )
      pvol_wat_cor_sol( map_wth(1) + 0 ) = pvol_wat_cor_sol( map_wth(1) + 1 )
    case(fldname_pvol_bc_ait_ins)
      pvol_bc_ait_ins( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      pvol_bc_ait_ins( map_wth(1) + 0 ) = pvol_bc_ait_ins( map_wth(1) + 1 )
    case(fldname_pvol_om_ait_ins)
      pvol_om_ait_ins( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      pvol_om_ait_ins( map_wth(1) + 0 ) = pvol_om_ait_ins( map_wth(1) + 1 )
    case(fldname_pvol_du_acc_ins)
      pvol_du_acc_ins( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      pvol_du_acc_ins( map_wth(1) + 0 ) = pvol_du_acc_ins( map_wth(1) + 1 )
    case(fldname_pvol_du_cor_ins)
      pvol_du_cor_ins( map_wth(1) + 1 : map_wth(1) + nlayers ) =               &
        real( ntp( :, i ), r_def )
      pvol_du_cor_ins( map_wth(1) + 0 ) = pvol_du_cor_ins( map_wth(1) + 1 )
    case(fldname_no2)
      no2( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( ntp( :, i ), r_def )
      no2( map_wth(1) + 0 ) = no2( map_wth(1) + 1 )
    case(fldname_bro)
      bro( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( ntp( :, i ), r_def )
      bro( map_wth(1) + 0 ) = bro( map_wth(1) + 1 )
    case(fldname_hcl)
      hcl( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( ntp( :, i ), r_def )
      hcl( map_wth(1) + 0 ) = hcl( map_wth(1) + 1 )
    case(fldname_o1d)
      o1d( map_wth(1) + 1 : map_wth(1) + nlayers ) =                           &
        real( ntp( :, i ), r_def )
      o1d( map_wth(1) + 0 ) = o1d( map_wth(1) + 1 )
     ! These fields are ntp if l_ukca_ro2_ntp=true, else tracers.
    case(fldname_meoo)
      meoo( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( ntp( :, i ), r_def )
      meoo( map_wth(1) + 0 ) = meoo( map_wth(1) + 1 )
    case(fldname_etoo)
      etoo( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( ntp( :, i ), r_def )
      etoo( map_wth(1) + 0 ) = etoo( map_wth(1) + 1 )
    case(fldname_meco3)
      meco3( map_wth(1) + 1 : map_wth(1) + nlayers ) =                         &
        real( ntp( :, i ), r_def )
      meco3( map_wth(1) + 0 ) = meco3( map_wth(1) + 1 )
    case(fldname_n_proo)
      n_proo( map_wth(1) + 1 : map_wth(1) + nlayers ) =                        &
        real( ntp( :, i ), r_def )
      n_proo( map_wth(1) + 0 ) = n_proo( map_wth(1) + 1 )
    case(fldname_i_proo)
      i_proo( map_wth(1) + 1 : map_wth(1) + nlayers ) =                        &
        real( ntp( :, i ), r_def )
      i_proo( map_wth(1) + 0 ) = i_proo( map_wth(1) + 1 )
    case(fldname_etco3)
      etco3( map_wth(1) + 1 : map_wth(1) + nlayers ) =                         &
        real( ntp( :, i ), r_def )
      etco3( map_wth(1) + 0 ) = etco3( map_wth(1) + 1 )
    case(fldname_mecoch2oo)
      mecoch2oo( map_wth(1) + 1 : map_wth(1) + nlayers ) =                     &
        real( ntp( :, i ), r_def )
      mecoch2oo( map_wth(1) + 0 ) = mecoch2oo( map_wth(1) + 1 )
    case(fldname_macro2)
      macro2( map_wth(1) + 1 : map_wth(1) + nlayers ) =                        &
        real( ntp( :, i ), r_def )
      macro2( map_wth(1) + 0 ) = macro2( map_wth(1) + 1 )
    case(fldname_iso2)
      iso2( map_wth(1) + 1 : map_wth(1) + nlayers ) =                          &
        real( ntp( :, i ), r_def )
      iso2( map_wth(1) + 0 ) = iso2( map_wth(1) + 1 )
    end select
  end do

  deallocate(ntp)
  deallocate(tracer)

end subroutine aerosol_ukca_code

end module aerosol_ukca_kernel_mod
