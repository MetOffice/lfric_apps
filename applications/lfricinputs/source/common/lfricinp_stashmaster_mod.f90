! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_stashmaster_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int64
! Shumlib modules
use f_shum_stashmaster_mod, only: shum_STASHmaster, f_shum_read_stashmaster
! LFRic modules
use log_mod, only: log_event, log_scratch_space, LOG_LEVEL_ERROR
use constants_mod, only: imdi

implicit none

! Define STASHmaster array. The size of which must match the size
! of the array in Shumlib. The array index will correspond to the
! stashcode. First two digits are the section code, remaining 3 digits
! are the items code. E.g. Section code 0, item code 3 will be at
! stashmaster(3), section code 2, item code 21 would be stashmaster(2021)
type(shum_STASHmaster), public :: stashmaster(99999)

! Parameters for horizontal grid code
! Atmospheric theta / p points
integer(kind=int64), public, parameter :: p_points = 1
! Atmospheric theta / p points - values over sea only
integer(kind=int64), public, parameter :: p_points_values_over_sea = 3
! Atmospheric u points on the 'c' grid
integer(kind=int64), public, parameter :: u_points = 18
! Atmospheric v points on the 'c' grid
integer(kind=int64), public, parameter :: v_points = 19
! Land compressed point
integer(kind=int64), public, parameter :: land_compressed = 21
! Ozone points
integer(kind=int64), public, parameter :: ozone_points = 22

! Parameters for level code
integer(kind=int64), public, parameter :: rho_levels = 1
integer(kind=int64), public, parameter :: theta_levels = 2
integer(kind=int64), public, parameter :: single_level = 5
integer(kind=int64), public, parameter :: soil_levels = 6

! Parameters for pseudo level type
integer(kind=int64), public, parameter :: snow_layers_and_tiles = 11
integer(KIND=int64), public, parameter :: ice_cats = 10

! STASHcodes
integer(kind=int64), public, parameter :: stashcode_u                 =   2
integer(kind=int64), public, parameter :: stashcode_v                 =   3
integer(kind=int64), public, parameter :: stashcode_theta             =   4
integer(kind=int64), public, parameter :: stashcode_orog_x_grad       =   5
integer(kind=int64), public, parameter :: stashcode_orog_y_grad       =   6
integer(kind=int64), public, parameter :: stashcode_unfilt_orog       =   7
integer(kind=int64), public, parameter :: stashcode_soil_moist        =   9

integer(kind=int64), public, parameter :: stashcode_q                 =  10
integer(kind=int64), public, parameter :: stashcode_qcf               =  12
integer(kind=int64), public, parameter :: stashcode_cca               =  13
integer(kind=int64), public, parameter :: stashcode_ccb               =  14
integer(kind=int64), public, parameter :: stashcode_cct               =  15
integer(kind=int64), public, parameter :: stashcode_cc_lwp            =  16
integer(kind=int64), public, parameter :: stashcode_sil_orog_rough    =  17
integer(kind=int64), public, parameter :: stashcode_hlf_pk_to_trf_ht  =  18

integer(kind=int64), public, parameter :: stashcode_soil_temp         =  20
integer(kind=int64), public, parameter :: stashcode_lcbase            =  21
integer(kind=int64), public, parameter :: stashcode_mean_canopyw      =  22
integer(kind=int64), public, parameter :: stashcode_snow_amount       =  23
integer(kind=int64), public, parameter :: stashcode_mean_snow         =  23
integer(kind=int64), public, parameter :: stashcode_surftemp          =  24
integer(kind=int64), public, parameter :: stashcode_tstar             =  24
integer(kind=int64), public, parameter :: stashcode_bl_depth          =  25
integer(kind=int64), public, parameter :: stashcode_rough_length      =  26
integer(kind=int64), public, parameter :: stashcode_z0                =  26
integer(kind=int64), public, parameter :: stashcode_snow_edge         =  27
integer(kind=int64), public, parameter :: stashcode_surf_z_curr       =  28
integer(kind=int64), public, parameter :: stashcode_surf_m_curr       =  29

integer(kind=int64), public, parameter :: stashcode_lsm               =  30
integer(kind=int64), public, parameter :: stashcode_icefrac           =  31
integer(kind=int64), public, parameter :: stashcode_icethick          =  32
integer(kind=int64), public, parameter :: stashcode_orog              =  33
integer(kind=int64), public, parameter :: stashcode_orog_var          =  34
integer(kind=int64), public, parameter :: stashcode_orog_gdxx         =  35
integer(kind=int64), public, parameter :: stashcode_orog_gdxy         =  36
integer(kind=int64), public, parameter :: stashcode_orog_gdyy         =  37
integer(kind=int64), public, parameter :: stashcode_ice_edge_ancil    =  38
integer(kind=int64), public, parameter :: stashcode_ice_edge_inancil  =  38
integer(kind=int64), public, parameter :: stashcode_tstar_anom        =  39

integer(kind=int64), public, parameter :: stashcode_vol_smc_wilt      =  40
integer(kind=int64), public, parameter :: stashcode_vol_smc_cri       =  41
integer(kind=int64), public, parameter :: stashcode_vol_smc_sat       =  43
integer(kind=int64), public, parameter :: stashcode_Ksat              =  44
integer(kind=int64), public, parameter :: stashcode_thermal_capacity  =  46
integer(kind=int64), public, parameter :: stashcode_thermal_conduct   =  47
integer(kind=int64), public, parameter :: stashcode_soil_suction      =  48
integer(kind=int64), public, parameter :: stashcode_sea_ice_temp      =  49

integer(kind=int64), public, parameter :: stashcode_veg_frac          =  50
integer(kind=int64), public, parameter :: stashcode_total_aero_emiss  =  57
integer(kind=int64), public, parameter :: stashcode_SO2_emiss         =  58
integer(kind=int64), public, parameter :: stashcode_dimethyl_sul_emiss=  59

integer(kind=int64), public, parameter :: stashcode_ozone             =  60
integer(kind=int64), public, parameter :: stashcode_e_trb             =  70
integer(kind=int64), public, parameter :: stashcode_tsq_trb           =  71
integer(kind=int64), public, parameter :: stashcode_qsq_trb           =  72
integer(kind=int64), public, parameter :: stashcode_cov_trb           =  73
integer(kind=int64), public, parameter :: stashcode_zhpar_shcu        =  74

integer(kind=int64), public, parameter :: stashcode_cloud_number      =  75
integer(kind=int64), public, parameter :: stashcode_rain_number       =  76
integer(kind=int64), public, parameter :: stashcode_rain_3mom         =  77
integer(kind=int64), public, parameter :: stashcode_ice_number        =  78
integer(kind=int64), public, parameter :: stashcode_snow_number       =  79
integer(kind=int64), public, parameter :: stashcode_snow_3mom         =  80
integer(kind=int64), public, parameter :: stashcode_graup_number      =  81
integer(kind=int64), public, parameter :: stashcode_graup_3mom        =  82

integer(kind=int64), public, parameter :: stashcode_activesol_liquid  =  83
integer(kind=int64), public, parameter :: stashcode_activesol_rain    =  84
integer(kind=int64), public, parameter :: stashcode_active_insol_ice  =  85
integer(kind=int64), public, parameter :: stashcode_active_sol_ice    =  86
integer(kind=int64), public, parameter :: stashcode_active_insol_liq  =  87
integer(kind=int64), public, parameter :: stashcode_active_sol_num    =  88
integer(kind=int64), public, parameter :: stashcode_active_insol_num  =  89

integer(kind=int64), public, parameter :: stashcode_total_aero        =  90
integer(kind=int64), public, parameter :: stashcode_flash_pot         =  91
integer(kind=int64), public, parameter :: stashcode_runoff_coast_out  =  93
integer(kind=int64), public, parameter :: stashcode_snow_on_ice       =  95
integer(kind=int64), public, parameter :: stashcode_ocnsrf_chlorophyll=  96
integer(kind=int64), public, parameter :: stashcode_chlorophyll       =  96
integer(kind=int64), public, parameter :: stashcode_z0m_soil          =  97

integer(kind=int64), public, parameter :: stashcode_blwvariance       =  99

integer(kind=int64), public, parameter :: stashcode_so2               = 101
integer(kind=int64), public, parameter :: stashcode_dms               = 102
integer(kind=int64), public, parameter :: stashcode_mmr_so4_aitken    = 103
integer(kind=int64), public, parameter :: stashcode_mmr_so4_accum     = 104
integer(kind=int64), public, parameter :: stashcode_mmr_so4_diss      = 105
integer(kind=int64), public, parameter :: stashcode_mmr_nh3           = 107
integer(kind=int64), public, parameter :: stashcode_mmr_bc_fr         = 108
integer(kind=int64), public, parameter :: stashcode_mmr_bc_ag         = 109
integer(kind=int64), public, parameter :: stashcode_mmr_bc_cl         = 110
integer(kind=int64), public, parameter :: stashcode_mmr_smoke_fr      = 111
integer(kind=int64), public, parameter :: stashcode_mmr_smoke_ag      = 112
integer(kind=int64), public, parameter :: stashcode_mmr_smoke_cl      = 113
integer(kind=int64), public, parameter :: stashcode_mmr_ocff_fr       = 114
integer(kind=int64), public, parameter :: stashcode_mmr_ocff_ag       = 115
integer(kind=int64), public, parameter :: stashcode_mmr_ocff_cl       = 116
integer(kind=int64), public, parameter :: stashcode_mmr_nitr_acc      = 117
integer(kind=int64), public, parameter :: stashcode_mmr_nitr_diss     = 118

integer(kind=int64), public, parameter :: stashcode_biom_elev_em_h1   = 119
integer(kind=int64), public, parameter :: stashcode_biom_elev_em_h2   = 120
integer(kind=int64), public, parameter :: stashcode_3d_nat_so2_em     = 121
integer(kind=int64), public, parameter :: stashcode_3d_oh_conc        = 122
integer(kind=int64), public, parameter :: stashcode_3d_ho2_conc       = 123
integer(kind=int64), public, parameter :: stashcode_3dh2o2_mixrat     = 124
integer(kind=int64), public, parameter :: stashcode_3d_ozone_mixrat   = 125
integer(kind=int64), public, parameter :: stashcode_hi_SO2_emiss_emiss= 126
integer(kind=int64), public, parameter :: stashcode_ammonia_gas_emiss = 127
integer(kind=int64), public, parameter :: stashcode_soot_surf         = 128
integer(kind=int64), public, parameter :: stashcode_soot_hi_lev       = 129

integer(kind=int64), public, parameter :: stashcode_biom_surf_em      = 130
integer(kind=int64), public, parameter :: stashcode_biom_elev_em      = 131
integer(kind=int64), public, parameter :: stashcode_dms_conc_sea      = 132
integer(kind=int64), public, parameter :: stashcode_dms_conc_sw       = 132
integer(kind=int64), public, parameter :: stashcode_ocff_surf_emiss   = 134
integer(kind=int64), public, parameter :: stashcode_ocff_hilev_emiss  = 135

integer(kind=int64), public, parameter :: stashcode_w                 = 150
integer(kind=int64), public, parameter :: stashcode_riv_sequence      = 151
integer(kind=int64), public, parameter :: stashcode_riv_direction     = 152
integer(kind=int64), public, parameter :: stashcode_riv_storage       = 153
integer(kind=int64), public, parameter :: stashcode_riv_number        = 154
integer(kind=int64), public, parameter :: stashcode_riv_iarea         = 160
integer(kind=int64), public, parameter :: stashcode_riv_slope         = 161
integer(kind=int64), public, parameter :: stashcode_riv_flowobs1      = 162
integer(kind=int64), public, parameter :: stashcode_riv_inext         = 163
integer(kind=int64), public, parameter :: stashcode_riv_jnext         = 164
integer(kind=int64), public, parameter :: stashcode_riv_land          = 165
integer(kind=int64), public, parameter :: stashcode_riv_sfcstorage    = 166
integer(kind=int64), public, parameter :: stashcode_riv_substorage    = 167
integer(kind=int64), public, parameter :: stashcode_riv_flowin        = 168
integer(kind=int64), public, parameter :: stashcode_riv_bflowin       = 169

integer(KIND=int64), public, parameter :: stashcode_cpl_sw_rad_sea    = 171
integer(KIND=int64), public, parameter :: stashcode_cpl_sw_surf_sea   = 172
integer(KIND=int64), public, parameter :: stashcode_cpl_lw_rad_sea    = 173
integer(KIND=int64), public, parameter :: stashcode_cpl_lw_surf_sea   = 174
integer(KIND=int64), public, parameter :: stashcode_cpl_xcomp_windstr = 175
integer(KIND=int64), public, parameter :: stashcode_cpl_ycomp_windstr = 176
integer(kind=int64), public, parameter :: stashcode_ice_subl_cat      = 182
integer(kind=int64), public, parameter :: stashcode_ls_snow_rate      = 187
integer(kind=int64), public, parameter :: stashcode_conv_rain_rate    = 188
integer(kind=int64), public, parameter :: stashcode_iceberg_calving   = 190
integer(kind=int64), public, parameter :: stashcode_sstfrz            = 194
integer(kind=int64), public, parameter :: stashcode_tstar_ice_cat_cpl = 195
integer(kind=int64), public, parameter :: stashcode_riv_outflow_cpl   = 198

integer(kind=int64), public, parameter :: stashcode_u_compnt_pert     = 202
integer(kind=int64), public, parameter :: stashcode_v_compnt_pert     = 203
integer(kind=int64), public, parameter :: stashcode_clapp_hb          = 207
integer(kind=int64), public, parameter :: stashcode_3d_cca            = 211
integer(kind=int64), public, parameter :: stashcode_3d_ccw            = 212
integer(kind=int64), public, parameter :: stashcode_can_conduct       = 213
integer(kind=int64), public, parameter :: stashcode_unfrozen_soil     = 214
integer(kind=int64), public, parameter :: stashcode_frozen_soil       = 215
integer(kind=int64), public, parameter :: stashcode_frac_surf_type    = 216
integer(kind=int64), public, parameter :: stashcode_lai               = 217
integer(kind=int64), public, parameter :: stashcode_canopy_height     = 218
integer(kind=int64), public, parameter :: stashcode_disturb_frac_veg  = 219

integer(kind=int64), public, parameter :: stashcode_snw_free_alb_bs   = 220
integer(kind=int64), public, parameter :: stashcode_snow_soot_tile    = 221
integer(kind=int64), public, parameter :: stashcode_soil_carbon_cont  = 223
integer(kind=int64), public, parameter :: stashcode_npp_pft_acc       = 224
integer(kind=int64), public, parameter :: stashcode_g_lf_pft_acc      = 225
integer(kind=int64), public, parameter :: stashcode_g_ph_lf_pft_acc   = 226
integer(kind=int64), public, parameter :: stashcode_rsp_w_pft_acc     = 227
integer(kind=int64), public, parameter :: stashcode_rsp_s_acc         = 228
integer(kind=int64), public, parameter :: stashcode_can_water_tile    = 229

integer(kind=int64), public, parameter :: stashcode_catch_tile        = 230
integer(kind=int64), public, parameter :: stashcode_rgrain            = 231
integer(kind=int64), public, parameter :: stashcode_tstar_tile        = 233
integer(kind=int64), public, parameter :: stashcode_z0_tile           = 234
integer(kind=int64), public, parameter :: stashcode_infil_max_tile    = 236
integer(kind=int64), public, parameter :: stashcode_sw_down_tile      = 237
integer(kind=int64), public, parameter :: stashcode_sw_down           = 238
integer(kind=int64), public, parameter :: stashcode_lw_up_diff        = 239

integer(kind=int64), public, parameter :: stashcode_snow_tile         = 240
integer(kind=int64), public, parameter :: stashcode_catch_snow        = 241
integer(kind=int64), public, parameter :: stashcode_snow_grnd         = 242
integer(kind=int64), public, parameter :: stashcode_surf_sw_alb       = 243
integer(kind=int64), public, parameter :: stashcode_surf_vis_alb      = 244
integer(kind=int64), public, parameter :: stashcode_surf_nir_alb      = 245
integer(kind=int64), public, parameter :: stashcode_z0h_tile          = 246

integer(kind=int64), public, parameter :: stashcode_CO2_surf_emiss    = 251
integer(kind=int64), public, parameter :: stashcode_rhor2             = 253
integer(kind=int64), public, parameter :: stashcode_qcl               = 254
integer(kind=int64), public, parameter :: stashcode_exner             = 255
integer(kind=int64), public, parameter :: stashcode_u_adv             = 256
integer(kind=int64), public, parameter :: stashcode_v_adv             = 257
integer(kind=int64), public, parameter :: stashcode_w_adv             = 258
integer(kind=int64), public, parameter :: stashcode_n_turb_mixlvs     = 259

integer(kind=int64), public, parameter :: stashcode_lvl_bse_dp_sc     = 260
integer(kind=int64), public, parameter :: stashcode_lvl_top_dp_sc     = 261
integer(kind=int64), public, parameter :: stashcode_bl_conv_flag      = 262
integer(kind=int64), public, parameter :: stashcode_turb_temp         = 263
integer(kind=int64), public, parameter :: stashcode_turb_humid        = 264
integer(kind=int64), public, parameter :: stashcode_area_cf           = 265
integer(kind=int64), public, parameter :: stashcode_bulk_cf           = 266
integer(kind=int64), public, parameter :: stashcode_liquid_cf         = 267
integer(kind=int64), public, parameter :: stashcode_frozen_cf         = 268
integer(kind=int64), public, parameter :: stashcode_sfc_zonal_cur     = 269

integer(kind=int64), public, parameter :: stashcode_sfc_merid_cur     = 270
integer(kind=int64), public, parameter :: stashcode_qcf2              = 271
integer(kind=int64), public, parameter :: stashcode_qrain             = 272
integer(kind=int64), public, parameter :: stashcode_qgraup            = 273
integer(kind=int64), public, parameter :: stashcode_top_ind_mean      = 274
integer(kind=int64), public, parameter :: stashcode_Ti_Mean           = 274
integer(kind=int64), public, parameter :: stashcode_top_ind_stddev    = 275
integer(kind=int64), public, parameter :: stashcode_Ti_Sig            = 275
integer(kind=int64), public, parameter :: stashcode_fexp              = 276
integer(kind=int64), public, parameter :: stashcode_gamtot            = 277
integer(kind=int64), public, parameter :: stashcode_zw                = 278
integer(kind=int64), public, parameter :: stashcode_fsat              = 279

integer(kind=int64), public, parameter :: stashcode_fwetl             = 280
integer(kind=int64), public, parameter :: stashcode_sthzw             = 281
integer(kind=int64), public, parameter :: stashcode_a_fsat            = 282
integer(kind=int64), public, parameter :: stashcode_c_fsat            = 283
integer(kind=int64), public, parameter :: stashcode_a_fwet            = 284
integer(kind=int64), public, parameter :: stashcode_c_fwet            = 285

integer(kind=int64), public, parameter :: stashcode_disturb_frac_veg_prev = 286
integer(kind=int64), public, parameter :: stashcode_wood_prod_fast    = 287
integer(kind=int64), public, parameter :: stashcode_wood_prod_med     = 288
integer(kind=int64), public, parameter :: stashcode_wood_prod_slow    = 289

integer(kind=int64), public, parameter :: stashcode_flake_depth       = 291
integer(kind=int64), public, parameter :: stashcode_flake_fetch       = 292
integer(kind=int64), public, parameter :: stashcode_flake_t_mean      = 293
integer(kind=int64), public, parameter :: stashcode_flake_t_mxl       = 294
integer(kind=int64), public, parameter :: stashcode_flake_t_ice       = 295
integer(kind=int64), public, parameter :: stashcode_flake_h_mxl       = 296
integer(kind=int64), public, parameter :: stashcode_flake_h_ice       = 297
integer(kind=int64), public, parameter :: stashcode_flake_shape       = 298
integer(kind=int64), public, parameter :: stashcode_flake_g_over_dt   = 299

integer(kind=int64), public, parameter :: stashcode_user_anc_sing1    = 301
integer(kind=int64), public, parameter :: stashcode_user_anc_sing20   = 320
integer(kind=int64), public, parameter :: stashcode_user_anc_mult1    = 321

integer(kind=int64), public, parameter :: stashcode_user_anc_mult20   = 340
integer(kind=int64), public, parameter :: stashcode_tppsozone         = 341
integer(kind=int64), public, parameter :: stashcode_deep_conv_flag    = 342
integer(kind=int64), public, parameter :: stashcode_past_conv_precip  = 343
integer(kind=int64), public, parameter :: stashcode_past_conv_depth   = 344
integer(kind=int64), public, parameter :: stashcode_cca_dp            = 345
integer(kind=int64), public, parameter :: stashcode_cca_md            = 346
integer(kind=int64), public, parameter :: stashcode_cca_sh            = 347
integer(kind=int64), public, parameter :: stashcode_total_precip      = 348

integer(kind=int64), public, parameter :: stashcode_clim_biogenic_aero= 351
integer(kind=int64), public, parameter :: stashcode_clim_delta_aero   = 371
integer(kind=int64), public, parameter :: stashcode_snowdep_grd_tile  = 376
integer(kind=int64), public, parameter :: stashcode_snowpack_bk_dens  = 377

integer(kind=int64), public, parameter :: stashcode_nsnow_layrs_tiles = 380
integer(kind=int64), public, parameter :: stashcode_snow_laythk_tiles = 381
integer(kind=int64), public, parameter :: stashcode_snow_ice_tile     = 382
integer(kind=int64), public, parameter :: stashcode_snow_liq_tile     = 383
integer(kind=int64), public, parameter :: stashcode_snow_T_tile       = 384
integer(kind=int64), public, parameter :: stashcode_snow_laydns_tiles = 385
integer(kind=int64), public, parameter :: stashcode_snow_grnsiz_tiles = 386
integer(kind=int64), public, parameter :: stashcode_etadot            = 387
integer(kind=int64), public, parameter :: stashcode_thetavd           = 388
integer(kind=int64), public, parameter :: stashcode_dry_rho           = 389

integer(kind=int64), public, parameter :: stashcode_psiw_surf         = 390
integer(kind=int64), public, parameter :: stashcode_psiw_lid          = 397
integer(kind=int64), public, parameter :: stashcode_mv                = 391
integer(kind=int64), public, parameter :: stashcode_mcl               = 392
integer(kind=int64), public, parameter :: stashcode_mcf               = 393
integer(kind=int64), public, parameter :: stashcode_mr                = 394
integer(kind=int64), public, parameter :: stashcode_mgr               = 395
integer(kind=int64), public, parameter :: stashcode_mcf2              = 396
integer(kind=int64), public, parameter :: stashcode_exner_surf        = 398

integer(kind=int64), public, parameter :: stashcode_p                 = 407
integer(kind=int64), public, parameter :: stashcode_pstar             = 409
integer(kind=int64), public, parameter :: stashcode_ice_conc_cat      = 413
integer(kind=int64), public, parameter :: stashcode_ice_thick_cat     = 414
integer(kind=int64), public, parameter :: stashcode_ice_temp_cat      = 415
integer(kind=int64), public, parameter :: stashcode_ice_snow_depth_cat= 416
integer(kind=int64), public, parameter :: stashcode_dust_parent_clay  = 418
integer(kind=int64), public, parameter :: stashcode_dust_parent_silt  = 419

integer(kind=int64), public, parameter :: stashcode_dust_parent_sand  = 420
integer(kind=int64), public, parameter :: stashcode_dust_soil_mf1     = 421
integer(kind=int64), public, parameter :: stashcode_dust_soil_mf2     = 422
integer(kind=int64), public, parameter :: stashcode_dust_soil_mf3     = 423
integer(kind=int64), public, parameter :: stashcode_dust_soil_mf4     = 424
integer(kind=int64), public, parameter :: stashcode_dust_soil_mf5     = 425
integer(kind=int64), public, parameter :: stashcode_dust_soil_mf6     = 426
integer(kind=int64), public, parameter :: stashcode_soil_massfrac6    = 426
integer(kind=int64), public, parameter :: stashcode_pond_frac_cat     = 428
integer(kind=int64), public, parameter :: stashcode_pond_depth_cat    = 429

integer(kind=int64), public, parameter :: stashcode_dust1_mmr         = 431
integer(kind=int64), public, parameter :: stashcode_dust2_mmr         = 432
integer(kind=int64), public, parameter :: stashcode_dust3_mmr         = 433
integer(kind=int64), public, parameter :: stashcode_dust4_mmr         = 434
integer(kind=int64), public, parameter :: stashcode_dust5_mmr         = 435
integer(kind=int64), public, parameter :: stashcode_dust6_mmr         = 436

integer(kind=int64), public, parameter :: stashcode_ice_surf_cond_cat = 440
integer(kind=int64), public, parameter :: stashcode_ice_surf_temp_cat = 441

integer(kind=int64), public, parameter :: stashcode_soilnitro_dpm     = 442
integer(kind=int64), public, parameter :: stashcode_soilnitro_rpm     = 443
integer(kind=int64), public, parameter :: stashcode_soilnitro_bio     = 444
integer(kind=int64), public, parameter :: stashcode_soilnitro_hum     = 445
integer(kind=int64), public, parameter :: stashcode_soil_inorgnit     = 446
integer(kind=int64), public, parameter :: stashcode_nitrogen_deposition = 447

integer(kind=int64), public, parameter :: stashcode_crop_frac         = 448
integer(kind=int64), public, parameter :: stashcode_pasture_frac      = 458

integer(kind=int64), public, parameter :: stashcode_soilcarb_dpm      = 466
integer(kind=int64), public, parameter :: stashcode_soilcarb_rpm      = 467
integer(kind=int64), public, parameter :: stashcode_soilcarb_bio      = 468
integer(kind=int64), public, parameter :: stashcode_soilcarb_hum      = 469

integer(kind=int64), public, parameter :: stashcode_ozone_tracer      = 480
integer(kind=int64), public, parameter :: stashcode_o3_prod_loss      = 481
integer(kind=int64), public, parameter :: stashcode_o3_p_l_vmr        = 482
integer(kind=int64), public, parameter :: stashcode_o3_vmr            = 483
integer(kind=int64), public, parameter :: stashcode_o3_p_l_temp       = 484
integer(kind=int64), public, parameter :: stashcode_o3_temp           = 485
integer(kind=int64), public, parameter :: stashcode_o3_p_l_colo3      = 486
integer(kind=int64), public, parameter :: stashcode_o3_colo3          = 487

integer(kind=int64), public, parameter :: stashcode_dctemp_tile       = 490
integer(kind=int64), public, parameter :: stashcode_dctemp_ssi        = 491
integer(kind=int64), public, parameter :: stashcode_tm_trans          = 492
integer(kind=int64), public, parameter :: stashcode_ddmfx             = 493
integer(kind=int64), public, parameter :: stashcode_urbhgt            = 494
integer(kind=int64), public, parameter :: stashcode_urbhwr            = 495
integer(kind=int64), public, parameter :: stashcode_urbwrr            = 496
integer(kind=int64), public, parameter :: stashcode_urbdisp           = 497
integer(kind=int64), public, parameter :: stashcode_urbztm            = 498
integer(kind=int64), public, parameter :: stashcode_urbalbwl          = 499

integer(kind=int64), public, parameter :: stashcode_urbalbrd          = 500
integer(kind=int64), public, parameter :: stashcode_urbemisw          = 501
integer(kind=int64), public, parameter :: stashcode_urbemisr          = 502
integer(kind=int64), public, parameter :: stashcode_land_frac         = 505
integer(kind=int64), public, parameter :: stashcode_tstar_land        = 506
integer(kind=int64), public, parameter :: stashcode_tstar_sea         = 507
integer(kind=int64), public, parameter :: stashcode_tstar_sice        = 508
integer(kind=int64), public, parameter :: stashcode_albedo_sice       = 509
integer(kind=int64), public, parameter :: stashcode_albedo_land       = 510

integer(kind=int64), public, parameter :: stashcode_u10_cpl           = 515
integer(kind=int64), public, parameter :: stashcode_v10_cpl           = 516
integer(kind=int64), public, parameter :: stashcode_charnock_cpl      = 517

integer(kind=int64), public, parameter :: stashcode_ux_ccp            = 569
integer(kind=int64), public, parameter :: stashcode_uy_ccp            = 570
integer(kind=int64), public, parameter :: stashcode_um_ccp            = 571
integer(kind=int64), public, parameter :: stashcode_g_ccp             = 572
integer(kind=int64), public, parameter :: stashcode_h_ccp             = 573
integer(kind=int64), public, parameter :: stashcode_riso_ccp          = 574
integer(kind=int64), public, parameter :: stashcode_rdir_ccp          = 575
integer(kind=int64), public, parameter :: stashcode_tsurf_elev_surft  = 576

! PV-tracers
integer(kind=int64), public, parameter :: stashcode_dPV_rad           = 577
integer(kind=int64), public, parameter :: stashcode_dPV_sw            = 578
integer(kind=int64), public, parameter :: stashcode_dPV_lw            = 579
integer(kind=int64), public, parameter :: stashcode_dPV_mic           = 580
integer(kind=int64), public, parameter :: stashcode_dPV_gwd           = 581
integer(kind=int64), public, parameter :: stashcode_dPV_ph1           = 582
integer(kind=int64), public, parameter :: stashcode_dPV_conv          = 583
integer(kind=int64), public, parameter :: stashcode_dPV_bl            = 584
integer(kind=int64), public, parameter :: stashcode_dPV_stph          = 585
integer(kind=int64), public, parameter :: stashcode_dPV_cld           = 586
integer(kind=int64), public, parameter :: stashcode_dPV_iau           = 587
integer(kind=int64), public, parameter :: stashcode_dPV_nud           = 588
integer(kind=int64), public, parameter :: stashcode_dPV_tot           = 589
integer(kind=int64), public, parameter :: stashcode_dEPS_I            = 590
integer(kind=int64), public, parameter :: stashcode_dPV_sol           = 591
integer(kind=int64), public, parameter :: stashcode_dPV_mass          = 592
integer(kind=int64), public, parameter :: stashcode_dPV_0             = 593

! Stochastic Physics
integer(kind=int64), public, parameter :: stashcode_bl_pert_rand_fld  = 595
integer(kind=int64), public, parameter :: stashcode_bl_pert_flag      = 596

! More PV-tracers, diab-friction split
integer(kind=int64), public, parameter :: stashcode_dPV_conv_d        = 620
integer(kind=int64), public, parameter :: stashcode_dPV_conv_f        = 621
integer(kind=int64), public, parameter :: stashcode_dPV_bl_d          = 622
integer(kind=int64), public, parameter :: stashcode_dPV_bl_f          = 623
integer(kind=int64), public, parameter :: stashcode_dPV_PC2c          = 624
! Theta tracers
integer(kind=int64), public, parameter :: stashcode_dtheta_0          = 600
integer(kind=int64), public, parameter :: stashcode_dtheta_bl         = 601
integer(kind=int64), public, parameter :: stashcode_dtheta_bl_mix     = 602
integer(kind=int64), public, parameter :: stashcode_dtheta_bl_LH      = 603
integer(kind=int64), public, parameter :: stashcode_dtheta_conv       = 604
integer(kind=int64), public, parameter :: stashcode_dtheta_mic        = 605
integer(kind=int64), public, parameter :: stashcode_dtheta_rad        = 606
integer(kind=int64), public, parameter :: stashcode_dtheta_SW         = 607
integer(kind=int64), public, parameter :: stashcode_dtheta_LW         = 608
integer(kind=int64), public, parameter :: stashcode_dtheta_slow       = 609
integer(kind=int64), public, parameter :: stashcode_dtheta_cld        = 610
integer(kind=int64), public, parameter :: stashcode_dtheta_PC2c       = 611

! INFERNO Ignition Ancillaries
integer(kind=int64), public, parameter :: stashcode_flash_rate_ancil  = 626
integer(kind=int64), public, parameter :: stashcode_pop_den_ancil     = 627
integer(kind=int64), public, parameter :: stashcode_wealth_index_ancil= 628

! Irrigation
integer(kind=int64), public, parameter :: stashcode_sthu_irr          = 630
integer(kind=int64), public, parameter :: stashcode_frac_irr          = 631

! New
integer(kind=int64), public, parameter :: stashcode_qt                = 16207

! Chemistry fields
integer(kind=int64), public, parameter :: stashcode_o3                = 34001
integer(kind=int64), public, parameter :: stashcode_no                = 34002
integer(kind=int64), public, parameter :: stashcode_no3               = 34003
integer(kind=int64), public, parameter :: stashcode_n2o5              = 34005
integer(kind=int64), public, parameter :: stashcode_ho2no2            = 34006
integer(kind=int64), public, parameter :: stashcode_hono2             = 34007
integer(kind=int64), public, parameter :: stashcode_h2o2              = 34008
integer(kind=int64), public, parameter :: stashcode_ch4               = 34009
integer(kind=int64), public, parameter :: stashcode_co                = 34010
integer(kind=int64), public, parameter :: stashcode_hcho              = 34011
integer(kind=int64), public, parameter :: stashcode_meooh             = 34012
integer(kind=int64), public, parameter :: stashcode_hono              = 34013
integer(kind=int64), public, parameter :: stashcode_c2h6              = 34014
integer(kind=int64), public, parameter :: stashcode_etooh             = 34015
integer(kind=int64), public, parameter :: stashcode_mecho             = 34016
integer(kind=int64), public, parameter :: stashcode_pan               = 34017
integer(kind=int64), public, parameter :: stashcode_c3h8              = 34018
integer(kind=int64), public, parameter :: stashcode_n_prooh           = 34019
integer(kind=int64), public, parameter :: stashcode_i_prooh           = 34020
integer(kind=int64), public, parameter :: stashcode_etcho             = 34021
integer(kind=int64), public, parameter :: stashcode_me2co             = 34022
integer(kind=int64), public, parameter :: stashcode_mecoch2ooh        = 34023
integer(kind=int64), public, parameter :: stashcode_ppan              = 34024
integer(kind=int64), public, parameter :: stashcode_meono2            = 34025
integer(kind=int64), public, parameter :: stashcode_c5h8              = 34027
integer(kind=int64), public, parameter :: stashcode_isooh             = 34028
integer(kind=int64), public, parameter :: stashcode_ison              = 34029
integer(kind=int64), public, parameter :: stashcode_macr              = 34030
integer(kind=int64), public, parameter :: stashcode_macrooh           = 34031
integer(kind=int64), public, parameter :: stashcode_mpan              = 34032
integer(kind=int64), public, parameter :: stashcode_hacet             = 34033
integer(kind=int64), public, parameter :: stashcode_mgly              = 34034
integer(kind=int64), public, parameter :: stashcode_nald              = 34035
integer(kind=int64), public, parameter :: stashcode_hcooh             = 34036
integer(kind=int64), public, parameter :: stashcode_meco3h            = 34037
integer(kind=int64), public, parameter :: stashcode_meco2h            = 34038
integer(kind=int64), public, parameter :: stashcode_cl                = 34041
integer(kind=int64), public, parameter :: stashcode_clo               = 34042
integer(kind=int64), public, parameter :: stashcode_cl2o2             = 34043
integer(kind=int64), public, parameter :: stashcode_oclo              = 34044
integer(kind=int64), public, parameter :: stashcode_br                = 34045
integer(kind=int64), public, parameter :: stashcode_brcl              = 34047
integer(kind=int64), public, parameter :: stashcode_brono2            = 34048
integer(kind=int64), public, parameter :: stashcode_n2o               = 34049
integer(kind=int64), public, parameter :: stashcode_hocl              = 34051
integer(kind=int64), public, parameter :: stashcode_hbr               = 34052
integer(kind=int64), public, parameter :: stashcode_hobr              = 34053
integer(kind=int64), public, parameter :: stashcode_clono2            = 34054
integer(kind=int64), public, parameter :: stashcode_cfcl3             = 34055
integer(kind=int64), public, parameter :: stashcode_cf2cl2            = 34056
integer(kind=int64), public, parameter :: stashcode_mebr              = 34057
integer(kind=int64), public, parameter :: stashcode_n                 = 34058
integer(kind=int64), public, parameter :: stashcode_o3p               = 34059
integer(kind=int64), public, parameter :: stashcode_h2                = 34070
integer(kind=int64), public, parameter :: stashcode_dms_mmr           = 34071
integer(kind=int64), public, parameter :: stashcode_so2_mmr           = 34072
integer(kind=int64), public, parameter :: stashcode_h2so4             = 34073
integer(kind=int64), public, parameter :: stashcode_msa               = 34074
integer(kind=int64), public, parameter :: stashcode_dmso              = 34075
integer(kind=int64), public, parameter :: stashcode_nh3               = 34076
integer(kind=int64), public, parameter :: stashcode_cs2               = 34077
integer(kind=int64), public, parameter :: stashcode_csul              = 34078
integer(kind=int64), public, parameter :: stashcode_h2s               = 34079
integer(kind=int64), public, parameter :: stashcode_h                 = 34080
integer(kind=int64), public, parameter :: stashcode_oh                = 34081
integer(kind=int64), public, parameter :: stashcode_ho2               = 34082
integer(kind=int64), public, parameter :: stashcode_meoo              = 34083
integer(kind=int64), public, parameter :: stashcode_etoo              = 34084
integer(kind=int64), public, parameter :: stashcode_meco3             = 34085
integer(kind=int64), public, parameter :: stashcode_n_proo            = 34086
integer(kind=int64), public, parameter :: stashcode_i_proo            = 34087
integer(kind=int64), public, parameter :: stashcode_etco3             = 34088
integer(kind=int64), public, parameter :: stashcode_mecoch2oo         = 34089
integer(kind=int64), public, parameter :: stashcode_meoh              = 34090
integer(kind=int64), public, parameter :: stashcode_monoterpene       = 34091
integer(kind=int64), public, parameter :: stashcode_sec_org           = 34092
integer(kind=int64), public, parameter :: stashcode_so3               = 34094
integer(kind=int64), public, parameter :: stashcode_lumped_n          = 34098
integer(kind=int64), public, parameter :: stashcode_lumped_br         = 34099
integer(kind=int64), public, parameter :: stashcode_lumped_cl         = 34100
integer(kind=int64), public, parameter :: stashcode_n_nuc_sol         = 34101
integer(kind=int64), public, parameter :: stashcode_nuc_sol_su        = 34102
integer(kind=int64), public, parameter :: stashcode_n_ait_sol         = 34103
integer(kind=int64), public, parameter :: stashcode_ait_sol_su        = 34104
integer(kind=int64), public, parameter :: stashcode_ait_sol_bc        = 34105
integer(kind=int64), public, parameter :: stashcode_ait_sol_om        = 34106
integer(kind=int64), public, parameter :: stashcode_n_acc_sol         = 34107
integer(kind=int64), public, parameter :: stashcode_acc_sol_su        = 34108
integer(kind=int64), public, parameter :: stashcode_acc_sol_bc        = 34109
integer(kind=int64), public, parameter :: stashcode_acc_sol_om        = 34110
integer(kind=int64), public, parameter :: stashcode_acc_sol_ss        = 34111
integer(kind=int64), public, parameter :: stashcode_acc_sol_du        = 34112
integer(kind=int64), public, parameter :: stashcode_n_cor_sol         = 34113
integer(kind=int64), public, parameter :: stashcode_cor_sol_su        = 34114
integer(kind=int64), public, parameter :: stashcode_cor_sol_bc        = 34115
integer(kind=int64), public, parameter :: stashcode_cor_sol_om        = 34116
integer(kind=int64), public, parameter :: stashcode_cor_sol_ss        = 34117
integer(kind=int64), public, parameter :: stashcode_cor_sol_du        = 34118
integer(kind=int64), public, parameter :: stashcode_n_ait_ins         = 34119
integer(kind=int64), public, parameter :: stashcode_ait_ins_bc        = 34120
integer(kind=int64), public, parameter :: stashcode_ait_ins_om        = 34121
integer(kind=int64), public, parameter :: stashcode_n_acc_ins         = 34122
integer(kind=int64), public, parameter :: stashcode_acc_ins_du        = 34123
integer(kind=int64), public, parameter :: stashcode_n_cor_ins         = 34124
integer(kind=int64), public, parameter :: stashcode_cor_ins_du        = 34125
integer(kind=int64), public, parameter :: stashcode_nuc_sol_om        = 34126
integer(kind=int64), public, parameter :: stashcode_passive_o3        = 34149
integer(kind=int64), public, parameter :: stashcode_age_of_air        = 34150

integer(kind=int64), public, parameter :: stashcode_hcl               = 34992
integer(kind=int64), public, parameter :: stashcode_bro               = 34994
integer(kind=int64), public, parameter :: stashcode_no2               = 34996
integer(kind=int64), public, parameter :: stashcode_o1d               = 34997


! STASHmaster element codes
integer(kind=int64), public, parameter :: model = 1
integer(kind=int64), public, parameter :: section = 2
integer(kind=int64), public, parameter :: item = 3
integer(kind=int64), public, parameter :: name = 4
integer(kind=int64), public, parameter :: space = 5
integer(kind=int64), public, parameter :: point = 6
integer(kind=int64), public, parameter :: time = 7
integer(kind=int64), public, parameter :: grid = 8
integer(kind=int64), public, parameter :: levelt = 9
integer(kind=int64), public, parameter :: levelf = 10
integer(kind=int64), public, parameter :: levell = 11
integer(kind=int64), public, parameter :: pseudt = 12
integer(kind=int64), public, parameter :: pseudf = 13
integer(kind=int64), public, parameter :: pseudl = 14
integer(kind=int64), public, parameter :: levcom = 15
integer(kind=int64), public, parameter :: option = 16
integer(kind=int64), public, parameter :: version_mask = 17
integer(kind=int64), public, parameter :: halo = 18
integer(kind=int64), public, parameter :: datat = 19
integer(kind=int64), public, parameter :: dumpp = 20
integer(kind=int64), public, parameter :: packing_codes = 21
integer(kind=int64), public, parameter :: rotate = 22
integer(kind=int64), public, parameter :: ppfc = 23
integer(kind=int64), public, parameter :: user = 24
integer(kind=int64), public, parameter :: lbvc = 25
integer(kind=int64), public, parameter :: blev = 26
integer(kind=int64), public, parameter :: tlev = 27
integer(kind=int64), public, parameter :: rblevv = 28
integer(kind=int64), public, parameter :: cfll = 29
integer(kind=int64), public, parameter :: cfff = 30

private

public :: lfricinp_read_stashmaster, get_stashmaster_item

contains

!-------------------------------------------------------------------

subroutine lfricinp_read_stashmaster(stashmaster_path)
! Description:
!  Read in STASHmaster using shumlib and check return code
!  of shumlib function.

implicit none
! Arguments
character(len=*), intent(in) :: stashmaster_path

! Error handling
integer(kind=int64) :: status
character(len=1024) :: message

status = f_shum_read_stashmaster(stashmaster_path, stashmaster,  &
                                 message)
if (status /= 0_int64) then
  call log_event(message, LOG_LEVEL_ERROR)
end if

end subroutine lfricinp_read_stashmaster

!-------------------------------------------------------------------

function get_stashmaster_item(stashcode, item_id) result(item)
! Currently only coded for 64-bit scalar integers
implicit none

integer(kind=int64), intent(in) :: stashcode
integer(kind=int64), intent(in) :: item_id
! Result
integer(kind=int64) :: item

item = imdi

! Check that the STASHmaster record exists.
if (.not. associated(stashmaster(stashcode) % record)) then
  write(log_scratch_space, '(A,I0)') &
       "Unassociated STASHmaster record for stashcode ", stashcode
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

select case(item_id)
case(grid)
  item = stashmaster(stashcode) % record % grid
case(levelt)
  item = stashmaster(stashcode) % record % levelt
case(levelf)
  item = stashmaster(stashcode) % record % levelf
case(levell)
  item = stashmaster(stashcode) % record % levell
case(pseudt)
  item = stashmaster(stashcode) % record % pseudt
case(pseudf)
  item = stashmaster(stashcode) % record % pseudf
case(pseudl)
  item = stashmaster(stashcode) % record % pseudl
case(ppfc)
  item = stashmaster(stashcode) % record % ppfc
case(lbvc)
  item = stashmaster(stashcode) % record % lbvc
case(datat)
  item = stashmaster(stashcode) % record % datat
case(cfff)
  item = stashmaster(stashcode) % record % cfff
case(cfll)
  item = stashmaster(stashcode) % record % cfll
case DEFAULT
  write(log_scratch_space, '(A,I0,A)') &
       "STASHmaster item id ", item_id, " not supported"
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end select

end function get_stashmaster_item

end module lfricinp_stashmaster_mod
