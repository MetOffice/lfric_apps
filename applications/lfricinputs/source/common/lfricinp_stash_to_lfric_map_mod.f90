! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_stash_to_lfric_map_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int64
! lfric modules
use constants_mod, only: str_def, imdi
use log_mod,  only : log_event, log_scratch_space, LOG_LEVEL_ERROR

implicit none

private

public :: lfricinp_init_stash_to_lfric_map, get_field_name,                    &
          get_lfric_field_kind, w2h_field, w3_field, w3_field_2d,              &
          w3_soil_field, wtheta_field

integer(kind=int64), parameter :: max_lfric_field_names = 500
character(len=str_def) :: field_name(max_lfric_field_names)
integer(kind=int64), save :: field_counter = 0

! The get_index array is dimensioned to the maximum number of stashcodes
! possible in the UM. The array index will correspond to the
! stashcode. First two digits are the section code, remaining 3 digits
! are the items code. E.g. Section code 0, item code 3 will be at
! get_index(3), section code 2, item code 21 would be get_index(2021)
! The contents of the array will an index to the character array of lfric
! field names. This intermediate step is used to avoid having a hugely
! oversized character array
integer(kind=int64) :: get_index(99999) = int(imdi, int64)

integer(kind=int64), parameter :: w2h_field     = 1
integer(kind=int64), parameter :: w3_field      = 2
integer(kind=int64), parameter :: w3_field_2d   = 3
integer(kind=int64), parameter :: w3_soil_field = 4
integer(kind=int64), parameter :: wtheta_field  = 5

contains

subroutine lfricinp_init_stash_to_lfric_map()

! Description:
!  Map stashcode to lfric field name

use lfricinp_stashmaster_mod, only: &
    stashcode_u, stashcode_v, stashcode_theta, stashcode_soil_moist,          &
    stashcode_soil_temp, stashcode_tstar, stashcode_bl_depth, stashcode_orog, &
    stashcode_ozone, stashcode_w,                                             &
    stashcode_cpl_sw_rad_sea,stashcode_cpl_sw_surf_sea,                       &
    stashcode_cpl_lw_rad_sea,stashcode_cpl_lw_surf_sea,                       &
    stashcode_cpl_xcomp_windstr,stashcode_cpl_ycomp_windstr,                  &
    stashcode_can_conduct,                                                    &
    stashcode_unfrozen_soil, stashcode_frozen_soil, stashcode_snow_soot_tile, &
    stashcode_can_water_tile, stashcode_rgrain, stashcode_tstar_tile,         &
    stashcode_snow_tile, stashcode_snow_grnd, stashcode_area_cf,              &
    stashcode_bulk_cf, stashcode_liquid_cf, stashcode_frozen_cf, stashcode_zw,&
    stashcode_fsat, stashcode_fwetl, stashcode_sthzw,                         &
    stashcode_snowdep_grd_tile, stashcode_snowpack_bk_dens,                   &
    stashcode_nsnow_layrs_tiles, stashcode_snow_laythk_tiles,                 &
    stashcode_snow_ice_tile, stashcode_snow_liq_tile, stashcode_snow_T_tile,  &
    stashcode_snow_grnsiz_tiles, stashcode_dry_rho, stashcode_mv,             &
    stashcode_ice_conc_cat,stashcode_ice_thick_cat,                           &
    stashcode_ice_temp_cat,stashcode_ice_snow_depth_cat,                      &
    stashcode_mcl, stashcode_mcf, stashcode_mr,                               &
    stashcode_ice_surf_cond_cat, stashcode_ice_surf_temp_cat,                 &
    stashcode_ddmfx,                                                          &
    stashcode_tstar_sea, stashcode_tstar_sice, stashcode_sea_ice_temp,        &
    stashcode_z0, stashcode_q, stashcode_qcf, stashcode_qcl, stashcode_qrain, &
    stashcode_rhor2, stashcode_lsm, stashcode_icefrac, stashcode_icethick,    &
    stashcode_total_aero, stashcode_z0h_tile, stashcode_dust1_mmr,            &
    stashcode_dust2_mmr, stashcode_ls_snow_rate, stashcode_conv_rain_rate,    &
    stashcode_qt , stashcode_p, stashcode_exner,                              &
    stashcode_o3p, stashcode_o1d, stashcode_o3, stashcode_n, stashcode_no,    &
    stashcode_no3, stashcode_no2, stashcode_n2o5, stashcode_ho2no2,           &
    stashcode_hono2, stashcode_h2o2, stashcode_ch4, stashcode_co,             &
    stashcode_hcho, stashcode_meoo, stashcode_meooh, stashcode_h,             &
    stashcode_oh, stashcode_ho2, stashcode_cl, stashcode_cl2o2,               &
    stashcode_clo, stashcode_oclo, stashcode_br, stashcode_bro,               &
    stashcode_brcl, stashcode_brono2, stashcode_n2o, stashcode_hcl,           &
    stashcode_hocl, stashcode_hbr, stashcode_hobr, stashcode_clono2,          &
    stashcode_cfcl3, stashcode_cf2cl2, stashcode_mebr, stashcode_hono,        &
    stashcode_c2h6, stashcode_etoo, stashcode_etooh, stashcode_mecho,         &
    stashcode_meco3, stashcode_pan, stashcode_c3h8, stashcode_n_proo,         &
    stashcode_i_proo, stashcode_n_prooh, stashcode_i_prooh, stashcode_etcho,  &
    stashcode_etco3, stashcode_me2co, stashcode_mecoch2oo,                    &
    stashcode_mecoch2ooh, stashcode_ppan, stashcode_meono2, stashcode_c5h8,   &
    stashcode_isooh, stashcode_ison, stashcode_macr, stashcode_macrooh,       &
    stashcode_mpan, stashcode_hacet, stashcode_mgly, stashcode_nald,          &
    stashcode_hcooh, stashcode_meco3h, stashcode_meco2h, stashcode_h2,        &
    stashcode_dms_mmr, stashcode_so2_mmr, stashcode_h2so4, stashcode_meoh,    &
    stashcode_msa, stashcode_dmso, stashcode_nh3, stashcode_cs2,              &
    stashcode_csul, stashcode_h2s, stashcode_so3, stashcode_passive_o3,       &
    stashcode_age_of_air, stashcode_lumped_n, stashcode_lumped_br,            &
    stashcode_lumped_cl, stashcode_monoterpene, stashcode_sec_org,            &
    stashcode_n_nuc_sol, stashcode_nuc_sol_su, stashcode_n_ait_sol,           &
    stashcode_ait_sol_su, stashcode_ait_sol_bc, stashcode_ait_sol_om,         &
    stashcode_n_acc_sol, stashcode_acc_sol_su, stashcode_acc_sol_bc,          &
    stashcode_acc_sol_om, stashcode_acc_sol_ss, stashcode_acc_sol_du,         &
    stashcode_n_cor_sol, stashcode_cor_sol_su, stashcode_cor_sol_bc,          &
    stashcode_cor_sol_om, stashcode_cor_sol_ss, stashcode_cor_sol_du,         &
    stashcode_n_ait_ins, stashcode_ait_ins_bc, stashcode_ait_ins_om,          &
    stashcode_n_acc_ins, stashcode_acc_ins_du, stashcode_n_cor_ins,           &
    stashcode_cor_ins_du, stashcode_nuc_sol_om

use lfricinp_regrid_options_mod, only: winds_on_w3

implicit none

field_name(:) = trim('unset')
field_counter = 0

! PLEASE KEEP THIS LIST OF SUPPORTED STASHCODES in NUMERICAL ORDER
if (winds_on_w3) then
  call map_field_name(stashcode_u, 'ew_wind')                        ! stash 2
  call map_field_name(stashcode_v, 'ns_wind')                        ! stash 3
else
  call map_field_name(stashcode_u, 'h_wind')                         ! stash 2&3
                                                 !combined into single W2H field
  call map_field_name(stashcode_v, 'h_wind')
ENDIF
call map_field_name(stashcode_theta, 'theta')                        ! stash 4
call map_field_name(stashcode_soil_moist, 'soil_moisture')           ! stash 9
call map_field_name(stashcode_q, 'q')                                ! stash 10
call map_field_name(stashcode_qcf, 'qcf')                            ! stash 12
call map_field_name(stashcode_soil_temp, 'soil_temperature')         ! stash 20
call map_field_name(stashcode_tstar, 'tstar')                        ! stash 24
call map_field_name(stashcode_bl_depth, 'zh')                        ! stash 25
call map_field_name(stashcode_z0, 'z0msea')                          ! stash 26
call map_field_name(stashcode_lsm, 'land_mask')                      ! stash 30
call map_field_name(stashcode_icefrac, 'icefrac')                    ! stash 31
call map_field_name(stashcode_icethick, 'icethick')                  ! stash 32
call map_field_name(stashcode_orog, 'surface_altitude')              ! stash 33
call map_field_name(stashcode_sea_ice_temp, 'sea_ice_temperature')   ! stash 49
call map_field_name(stashcode_ozone, 'ozone')                        ! stash 60
call map_field_name(stashcode_total_aero, 'total_aero')              ! stash 90
call map_field_name(stashcode_w, 'upward_wind')                      ! stash 150
call map_field_name(stashcode_cpl_sw_rad_sea,'cpl_sw_rad')           ! stash 171
call map_field_name(stashcode_cpl_sw_surf_sea,'cpl_sw_surf')         ! stash 172
call map_field_name(stashcode_cpl_lw_rad_sea,'cpl_lw_rad')           ! stash 173
call map_field_name(stashcode_cpl_lw_surf_sea,'cpl_lw_surf')         ! stash 174
call map_field_name(stashcode_cpl_xcomp_windstr,'cpl_xcomp_windstr') ! stash 175
call map_field_name(stashcode_cpl_ycomp_windstr,'cpl_ycomp_windstr') ! stash 176
call map_field_name(stashcode_ls_snow_rate, 'ls_snow_rate')          ! stash 187
call map_field_name(stashcode_conv_rain_rate, 'conv_rain_rate')      ! stash 188
call map_field_name(stashcode_can_conduct, 'surface_conductance')    ! stash 213
call map_field_name(stashcode_unfrozen_soil, 'unfrozen_soil_moisture')
                                                                     ! stash 214
call map_field_name(stashcode_frozen_soil, 'frozen_soil_moisture')   ! stash 215
call map_field_name(stashcode_snow_soot_tile, 'snow_soot')           ! stash 221
call map_field_name(stashcode_can_water_tile, 'tile_canopy_water')   ! stash 229
call map_field_name(stashcode_rgrain, 'tile_snow_rgrain')            ! stash 231
call map_field_name(stashcode_tstar_tile, 'tile_temperature')        ! stash 233
call map_field_name(stashcode_snow_tile , 'tile_snow_mass')          ! stash 240
call map_field_name(stashcode_snow_grnd, 'tile_snow_under_canopy')   ! stash 242
call map_field_name(stashcode_z0h_tile, 'z0h_tile')                  ! stash 246
call map_field_name(stashcode_rhor2, 'rho_r2')                       ! stash 253
call map_field_name(stashcode_qcl, 'qcl')                            ! stash 254
call map_field_name(stashcode_exner, 'exner')                        ! stash 255
call map_field_name(stashcode_area_cf, 'area_fraction')              ! stash 265
call map_field_name(stashcode_bulk_cf, 'bulk_fraction')              ! stash 266
call map_field_name(stashcode_liquid_cf, 'liquid_fraction')          ! stash 267
call map_field_name(stashcode_frozen_cf, 'frozen_fraction')          ! stash 268
call map_field_name(stashcode_qrain, 'qrain')                        ! stash 272
call map_field_name(stashcode_zw, 'water_table')                     ! stash 278
call map_field_name(stashcode_fsat, 'soil_sat_frac')                 ! stash 279
call map_field_name(stashcode_fwetl, 'soil_wet_frac')                ! stash 280
call map_field_name(stashcode_sthzw, 'wetness_under_soil')           ! stash 281
call map_field_name(stashcode_snowdep_grd_tile, 'tile_snow_depth')   ! stash 376
call map_field_name(stashcode_snowpack_bk_dens, 'tile_snowpack_density')
                                                                     ! stash 377
call map_field_name(stashcode_nsnow_layrs_tiles, 'tile_n_snow_layers')
                                                                     ! stash 380
call map_field_name(stashcode_snow_laythk_tiles,                    &! stash 381
     'tile_snow_layer_thickness')
call map_field_name(stashcode_snow_liq_tile, 'tile_snow_layer_liq_mass')
                                                                     ! stash 383
call map_field_name(stashcode_snow_ice_tile, 'tile_snow_layer_ice_mass')
                                                                     ! stash 382
call map_field_name(stashcode_snow_T_tile, 'tile_snow_layer_temp')   ! stash 384
call map_field_name(stashcode_snow_grnsiz_tiles,                    &! stash 386
     'tile_snow_layer_rgrain')
call map_field_name(stashcode_dry_rho, 'rho')                        ! stash 389
call map_field_name(stashcode_mv, 'm_v')                             ! stash 391
call map_field_name(stashcode_mcl, 'm_cl')                           ! stash 392
call map_field_name(stashcode_mcf, 'm_cf')                           ! stash 393
call map_field_name(stashcode_mr, 'm_r')                             ! stash 394
call map_field_name(stashcode_p, 'pressure')                         ! stash 407
call map_field_name(stashcode_ice_conc_cat, 'ice_conc_cat')              ! stash 413
call map_field_name(stashcode_ice_thick_cat, 'ice_thick_cat')            ! stash 414
call map_field_name(stashcode_ice_temp_cat, 'ice_temp_cat')              ! stash 415
call map_field_name(stashcode_ice_snow_depth_cat, 'ice_snow_depth_cat')  ! stash 416
call map_field_name(stashcode_dust1_mmr, 'dust1_mmr')                ! stash 431
call map_field_name(stashcode_dust2_mmr, 'dust2_mmr')                ! stash 432
call map_field_name(stashcode_ice_surf_cond_cat, 'ice_surf_cond_cat')    ! stash 440
call map_field_name(stashcode_ice_surf_temp_cat, 'ice_surf_temp_cat')    ! stash 441
call map_field_name(stashcode_ddmfx , 'dd_mf_cb')                    ! stash 493
call map_field_name(stashcode_tstar_sea, 'tstar_sea')                ! stash 507
call map_field_name(stashcode_tstar_sice, 'tstar_sea_ice')           ! stash 508
call map_field_name(stashcode_qt, 'qt')                              ! stash 16207

call map_field_name(stashcode_o3, 'o3')                              ! stash 34001
call map_field_name(stashcode_no, 'no')                              ! stash 34002
call map_field_name(stashcode_no3, 'no3')                            ! stash 34003
call map_field_name(stashcode_n2o5, 'n2o5')                          ! stash 34005
call map_field_name(stashcode_ho2no2, 'ho2no2')                      ! stash 34006
call map_field_name(stashcode_hono2, 'hono2')                        ! stash 34007
call map_field_name(stashcode_h2o2, 'h2o2')                          ! stash 34008
call map_field_name(stashcode_ch4, 'ch4')                            ! stash 34009
call map_field_name(stashcode_co, 'co')                              ! stash 34010
call map_field_name(stashcode_hcho, 'hcho')                          ! stash 34011
call map_field_name(stashcode_meooh, 'meooh')                        ! stash 34012
call map_field_name(stashcode_hono, 'hono')                          ! stash 34013
call map_field_name(stashcode_c2h6, 'c2h6')                          ! stash 34014
call map_field_name(stashcode_etooh, 'etooh')                        ! stash 34015
call map_field_name(stashcode_mecho, 'mecho')                        ! stash 34016
call map_field_name(stashcode_pan, 'pan')                            ! stash 34017
call map_field_name(stashcode_c3h8, 'c3h8')                          ! stash 34018
call map_field_name(stashcode_n_prooh, 'n_prooh')                    ! stash 34019
call map_field_name(stashcode_i_prooh, 'i_prooh')                    ! stash 34020
call map_field_name(stashcode_etcho, 'etcho')                        ! stash 34021
call map_field_name(stashcode_me2co, 'me2co')                        ! stash 34022
call map_field_name(stashcode_mecoch2ooh, 'mecoch2ooh')              ! stash 34023
call map_field_name(stashcode_ppan, 'ppan')                          ! stash 34024
call map_field_name(stashcode_meono2, 'meono2')                      ! stash 34025
call map_field_name(stashcode_c5h8, 'c5h8')                          ! stash 34027
call map_field_name(stashcode_isooh, 'isooh')                        ! stash 34028
call map_field_name(stashcode_ison, 'ison')                          ! stash 34029
call map_field_name(stashcode_macr, 'macr')                          ! stash 34030
call map_field_name(stashcode_macrooh, 'macrooh')                    ! stash 34031
call map_field_name(stashcode_mpan, 'mpan')                          ! stash 34032
call map_field_name(stashcode_hacet, 'hacet')                        ! stash 34033
call map_field_name(stashcode_mgly, 'mgly')                          ! stash 34034
call map_field_name(stashcode_nald, 'nald')                          ! stash 34035
call map_field_name(stashcode_hcooh, 'hcooh')                        ! stash 34036
call map_field_name(stashcode_meco3h, 'meco3h')                      ! stash 34037
call map_field_name(stashcode_meco2h, 'meco2h')                      ! stash 34038
call map_field_name(stashcode_cl, 'cl')                              ! stash 34041
call map_field_name(stashcode_clo, 'clo')                            ! stash 34042
call map_field_name(stashcode_cl2o2, 'cl2o2')                        ! stash 34043
call map_field_name(stashcode_oclo, 'oclo')                          ! stash 34044
call map_field_name(stashcode_br, 'br')                              ! stash 34045
call map_field_name(stashcode_brcl, 'brcl')                          ! stash 34047
call map_field_name(stashcode_brono2, 'brono2')                      ! stash 34048
call map_field_name(stashcode_n2o, 'n2o')                            ! stash 34049
call map_field_name(stashcode_hocl, 'hocl')                          ! stash 34051
call map_field_name(stashcode_hbr, 'hbr')                            ! stash 34052
call map_field_name(stashcode_hobr, 'hobr')                          ! stash 34053
call map_field_name(stashcode_clono2, 'clono2')                      ! stash 34054
call map_field_name(stashcode_cfcl3, 'cfcl3')                        ! stash 34055
call map_field_name(stashcode_cf2cl2, 'cf2cl2')                      ! stash 34056
call map_field_name(stashcode_mebr, 'mebr')                          ! stash 34057
call map_field_name(stashcode_n, 'n')                                ! stash 34058
call map_field_name(stashcode_o3p, 'o3p')                            ! stash 34059
call map_field_name(stashcode_h2, 'h2')                              ! stash 34070
call map_field_name(stashcode_dms_mmr, 'dms_mmr')                    ! stash 34071
call map_field_name(stashcode_so2_mmr, 'so2_mmr')                    ! stash 34072
call map_field_name(stashcode_h2so4, 'h2so4')                        ! stash 34073
call map_field_name(stashcode_msa, 'msa')                            ! stash 34074
call map_field_name(stashcode_dmso, 'dmso')                          ! stash 34075
call map_field_name(stashcode_nh3, 'nh3')                            ! stash 34076
call map_field_name(stashcode_cs2, 'cs2')                            ! stash 34077
call map_field_name(stashcode_csul, 'csul')                          ! stash 34078
call map_field_name(stashcode_h2s, 'h2s')                            ! stash 34079
call map_field_name(stashcode_h, 'h')                                ! stash 34080
call map_field_name(stashcode_oh, 'oh')                              ! stash 34081
call map_field_name(stashcode_ho2, 'ho2')                            ! stash 34082
call map_field_name(stashcode_meoo, 'meoo')                          ! stash 34083
call map_field_name(stashcode_etoo, 'etoo')                          ! stash 34084
call map_field_name(stashcode_meco3, 'meco3')                        ! stash 34085
call map_field_name(stashcode_n_proo, 'n_proo')                      ! stash 34086
call map_field_name(stashcode_i_proo, 'i_proo')                      ! stash 34087
call map_field_name(stashcode_etco3, 'etco3')                        ! stash 34088
call map_field_name(stashcode_mecoch2oo, 'mecoch2oo')                ! stash 34089
call map_field_name(stashcode_meoh, 'meoh')                          ! stash 34090
call map_field_name(stashcode_monoterpene, 'monoterpene')            ! stash 34091
call map_field_name(stashcode_sec_org, 'sec_org')                    ! stash 34092
call map_field_name(stashcode_so3, 'so3')                            ! stash 34094
call map_field_name(stashcode_lumped_n, 'lumped_n')                  ! stash 34098
call map_field_name(stashcode_lumped_br, 'lumped_br')                ! stash 34099
call map_field_name(stashcode_lumped_cl, 'lumped_cl')                ! stash 34100
call map_field_name(stashcode_n_nuc_sol, 'n_nuc_sol')                ! stash 34101
call map_field_name(stashcode_nuc_sol_su, 'nuc_sol_su')              ! stash 34102
call map_field_name(stashcode_n_ait_sol, 'n_ait_sol')                ! stash 34103
call map_field_name(stashcode_ait_sol_su, 'ait_sol_su')              ! stash 34104
call map_field_name(stashcode_ait_sol_bc, 'ait_sol_bc')              ! stash 34105
call map_field_name(stashcode_ait_sol_om, 'ait_sol_om')              ! stash 34106
call map_field_name(stashcode_n_acc_sol, 'n_acc_sol')                ! stash 34107
call map_field_name(stashcode_acc_sol_su, 'acc_sol_su')              ! stash 34108
call map_field_name(stashcode_acc_sol_bc, 'acc_sol_bc')              ! stash 34109
call map_field_name(stashcode_acc_sol_om, 'acc_sol_om')              ! stash 34110
call map_field_name(stashcode_acc_sol_ss, 'acc_sol_ss')              ! stash 34111
call map_field_name(stashcode_acc_sol_du, 'acc_sol_du')              ! stash 34112
call map_field_name(stashcode_n_cor_sol, 'n_cor_sol')                ! stash 34113
call map_field_name(stashcode_cor_sol_su, 'cor_sol_su')              ! stash 34114
call map_field_name(stashcode_cor_sol_bc, 'cor_sol_bc')              ! stash 34115
call map_field_name(stashcode_cor_sol_om, 'cor_sol_om')              ! stash 34116
call map_field_name(stashcode_cor_sol_ss, 'cor_sol_ss')              ! stash 34117
call map_field_name(stashcode_cor_sol_du, 'cor_sol_du')              ! stash 34118
call map_field_name(stashcode_n_ait_ins, 'n_ait_ins')                ! stash 34119
call map_field_name(stashcode_ait_ins_bc, 'ait_ins_bc')              ! stash 34120
call map_field_name(stashcode_ait_ins_om, 'ait_ins_om')              ! stash 34121
call map_field_name(stashcode_n_acc_ins, 'n_acc_ins')                ! stash 34122
call map_field_name(stashcode_acc_ins_du, 'acc_ins_du')              ! stash 34123
call map_field_name(stashcode_n_cor_ins, 'n_cor_ins')                ! stash 34124
call map_field_name(stashcode_cor_ins_du, 'cor_ins_du')              ! stash 34125
call map_field_name(stashcode_nuc_sol_om, 'nuc_sol_om')              ! stash 34126
call map_field_name(stashcode_passive_o3, 'passive_o3')              ! stash 34149
call map_field_name(stashcode_age_of_air, 'age_of_air')              ! stash 34150

call map_field_name(stashcode_hcl, 'hcl')                            ! stash 34992
call map_field_name(stashcode_bro, 'bro')                            ! stash 34994
call map_field_name(stashcode_no2, 'no2')                            ! stash 34996
call map_field_name(stashcode_o1d, 'o1d')                            ! stash 34997

! PLEASE KEEP THIS LIST OF SUPPORTED STASHCODES in NUMERICAL ORDER

end subroutine lfricinp_init_stash_to_lfric_map


subroutine map_field_name(stashcode, lfric_field_name)

implicit none

integer(kind=int64), intent(in) :: stashcode
character(len=*), intent(in) :: lfric_field_name

field_counter = field_counter + 1

if (field_counter >  max_lfric_field_names) then
  write(log_scratch_space, '(A)') "field_counter is greater than" // &
       " max_lfric_field_names. Recompile to increase."
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
else
  field_name(field_counter) = lfric_field_name
  get_index(stashcode) = field_counter
end if

end subroutine map_field_name


function get_field_name(stashcode) result(name)

implicit none

integer(kind=int64), intent(in) :: stashcode
character(len=str_def) :: name ! Result

if (get_index(stashcode) == imdi) then
  write(log_scratch_space, '(A,I0,A)') "Stashcode ", stashcode, &
       " has not been mapped to lfric field name"
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
else
  name = trim(field_name(get_index(stashcode)))
end if

end function get_field_name


function get_lfric_field_kind(stashcode) result(lfric_field_kind)

use lfricinp_regrid_options_mod, only: winds_on_w3
use lfricinp_stashmaster_mod,    only: get_stashmaster_item, levelt,       &
                                       rho_levels, theta_levels,           &
                                       single_level, soil_levels,          &
                                       stashcode_u, stashcode_v
implicit none

integer(kind=int64), intent(in) :: stashcode
integer(kind=int64) :: lfric_field_kind, level_code

level_code = get_stashmaster_item(stashcode, levelt)

if ((stashcode == stashcode_u .or. stashcode == stashcode_v) .and.     &
    (.not. winds_on_w3)) then
  lfric_field_kind = w2h_field

else if (level_code == rho_levels) then
  lfric_field_kind = w3_field

else if (level_code == theta_levels) then
  lfric_field_kind = wtheta_field

else if (level_code == single_level) then
  lfric_field_kind = w3_field_2d

else if (level_code == soil_levels) then
  lfric_field_kind = w3_soil_field

else

  write(log_scratch_space, '(A,I0,A)') "Stashcode ", stashcode, &
       " is not mapped to a lfric field type"
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)

end if

end function get_lfric_field_kind

end module lfricinp_stash_to_lfric_map_mod
