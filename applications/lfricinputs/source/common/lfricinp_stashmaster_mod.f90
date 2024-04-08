! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_stashmaster_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int64
! Shumlib modules
USE f_shum_stashmaster_mod, ONLY: shum_STASHmaster, f_shum_read_stashmaster
! LFRic modules
USE log_mod, ONLY: log_event, log_scratch_space, LOG_LEVEL_ERROR
USE constants_mod, ONLY: imdi

IMPLICIT NONE

! Define STASHmaster array. The size of which must match the size
! of the array in Shumlib. The array index will correspond to the
! stashcode. First two digits are the section code, remaining 3 digits
! are the items code. E.g. Section code 0, item code 3 will be at
! stashmaster(3), section code 2, item code 21 would be stashmaster(2021)
TYPE(shum_STASHmaster), PUBLIC :: stashmaster(99999)

! Parameters for horizontal grid code
! Atmospheric theta / p points
INTEGER(KIND=int64), PUBLIC, PARAMETER :: p_points = 1
! Atmospheric theta / p points - values over sea only
INTEGER(KIND=int64), PUBLIC, PARAMETER :: p_points_values_over_sea = 3
! Atmospheric u points on the 'c' grid
INTEGER(KIND=int64), PUBLIC, PARAMETER :: u_points = 18
! Atmospheric v points on the 'c' grid
INTEGER(KIND=int64), PUBLIC, PARAMETER :: v_points = 19
! Land compressed point
INTEGER(KIND=int64), PUBLIC, PARAMETER :: land_compressed = 21
! Ozone points
INTEGER(KIND=int64), PUBLIC, PARAMETER :: ozone_points = 22

! Parameters for level code
INTEGER(KIND=int64), PUBLIC, PARAMETER :: rho_levels = 1
INTEGER(KIND=int64), PUBLIC, PARAMETER :: theta_levels = 2
INTEGER(KIND=int64), PUBLIC, PARAMETER :: single_level = 5
INTEGER(KIND=int64), PUBLIC, PARAMETER :: soil_levels = 6

! Parameters for pseudo level type
INTEGER(KIND=int64), PUBLIC, PARAMETER :: snow_layers_and_tiles = 11

! STASHcodes
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_u                 =   2
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_v                 =   3
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_theta             =   4
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_orog_x_grad       =   5
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_orog_y_grad       =   6
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_unfilt_orog       =   7
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_soil_moist        =   9

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_q                 =  10
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_qcf               =  12
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_cca               =  13
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ccb               =  14
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_cct               =  15
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_cc_lwp            =  16
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_sil_orog_rough    =  17
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_hlf_pk_to_trf_ht  =  18

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_soil_temp         =  20
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_lcbase            =  21
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mean_canopyw      =  22
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_snow_amount       =  23
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mean_snow         =  23
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_surftemp          =  24
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_tstar             =  24
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_bl_depth          =  25
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_rough_length      =  26
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_z0                =  26
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_snow_edge         =  27
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_surf_z_curr       =  28
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_surf_m_curr       =  29

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_lsm               =  30
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_icefrac           =  31
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_icethick          =  32
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_orog              =  33
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_orog_var          =  34
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_orog_gdxx         =  35
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_orog_gdxy         =  36
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_orog_gdyy         =  37
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ice_edge_ancil    =  38
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ice_edge_inancil  =  38
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_tstar_anom        =  39

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_vol_smc_wilt      =  40
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_vol_smc_cri       =  41
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_vol_smc_sat       =  43
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_Ksat              =  44
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_thermal_capacity  =  46
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_thermal_conduct   =  47
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_soil_suction      =  48
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_sea_ice_temp      =  49

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_veg_frac          =  50
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_total_aero_emiss  =  57
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_SO2_emiss         =  58
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dimethyl_sul_emiss=  59

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ozone             =  60
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_e_trb             =  70
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_tsq_trb           =  71
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_qsq_trb           =  72
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_cov_trb           =  73
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_zhpar_shcu        =  74

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_cloud_number      =  75
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_rain_number       =  76
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_rain_3mom         =  77
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ice_number        =  78
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_snow_number       =  79
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_snow_3mom         =  80
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_graup_number      =  81
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_graup_3mom        =  82

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_activesol_liquid  =  83
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_activesol_rain    =  84
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_active_insol_ice  =  85
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_active_sol_ice    =  86
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_active_insol_liq  =  87
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_active_sol_num    =  88
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_active_insol_num  =  89

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_total_aero        =  90
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_flash_pot         =  91
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_runoff_coast_out  =  93
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_snow_on_ice       =  95
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ocnsrf_chlorophyll=  96
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_chlorophyll       =  96
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_z0m_soil          =  97

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_blwvariance       =  99

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_so2               = 101
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dms               = 102
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mmr_so4_aitken    = 103
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mmr_so4_accum     = 104
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mmr_so4_diss      = 105
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mmr_nh3           = 107
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mmr_bc_fr         = 108
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mmr_bc_ag         = 109
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mmr_bc_cl         = 110
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mmr_smoke_fr      = 111
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mmr_smoke_ag      = 112
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mmr_smoke_cl      = 113
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mmr_ocff_fr       = 114
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mmr_ocff_ag       = 115
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mmr_ocff_cl       = 116
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mmr_nitr_acc      = 117
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mmr_nitr_diss     = 118

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_biom_elev_em_h1   = 119
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_biom_elev_em_h2   = 120
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_3d_nat_so2_em     = 121
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_3d_oh_conc        = 122
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_3d_ho2_conc       = 123
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_3dh2o2_mixrat     = 124
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_3d_ozone_mixrat   = 125
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_hi_SO2_emiss_emiss= 126
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ammonia_gas_emiss = 127
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_soot_surf         = 128
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_soot_hi_lev       = 129

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_biom_surf_em      = 130
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_biom_elev_em      = 131
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dms_conc_sea      = 132
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dms_conc_sw       = 132
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ocff_surf_emiss   = 134
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ocff_hilev_emiss  = 135

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_w                 = 150
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_riv_sequence      = 151
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_riv_direction     = 152
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_riv_storage       = 153
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_riv_number        = 154
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_riv_iarea         = 160
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_riv_slope         = 161
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_riv_flowobs1      = 162
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_riv_inext         = 163
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_riv_jnext         = 164
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_riv_land          = 165
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_riv_sfcstorage    = 166
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_riv_substorage    = 167
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_riv_flowin        = 168
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_riv_bflowin       = 169
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ice_subl_cat      = 182
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ls_snow_rate      = 187
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_conv_rain_rate    = 188
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_iceberg_calving   = 190
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_sstfrz            = 194
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_tstar_ice_cat_cpl = 195
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_riv_outflow_cpl   = 198

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_u_compnt_pert     = 202
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_v_compnt_pert     = 203
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_clapp_hb          = 207
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_3d_cca            = 211
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_3d_ccw            = 212
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_can_conduct       = 213
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_unfrozen_soil     = 214
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_frozen_soil       = 215
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_frac_surf_type    = 216
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_lai               = 217
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_canopy_height     = 218
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_disturb_frac_veg  = 219

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_snw_free_alb_bs   = 220
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_snow_soot_tile    = 221
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_soil_carbon_cont  = 223
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_npp_pft_acc       = 224
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_g_lf_pft_acc      = 225
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_g_ph_lf_pft_acc   = 226
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_rsp_w_pft_acc     = 227
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_rsp_s_acc         = 228
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_can_water_tile    = 229

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_catch_tile        = 230
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_rgrain            = 231
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_tstar_tile        = 233
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_z0_tile           = 234
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_infil_max_tile    = 236
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_sw_down_tile      = 237
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_sw_down           = 238
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_lw_up_diff        = 239

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_snow_tile         = 240
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_catch_snow        = 241
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_snow_grnd         = 242
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_surf_sw_alb       = 243
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_surf_vis_alb      = 244
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_surf_nir_alb      = 245
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_z0h_tile          = 246

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_CO2_surf_emiss    = 251
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_rhor2             = 253
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_qcl               = 254
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_exner             = 255
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_u_adv             = 256
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_v_adv             = 257
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_w_adv             = 258
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_n_turb_mixlvs     = 259

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_lvl_bse_dp_sc     = 260
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_lvl_top_dp_sc     = 261
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_bl_conv_flag      = 262
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_turb_temp         = 263
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_turb_humid        = 264
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_area_cf           = 265
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_bulk_cf           = 266
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_liquid_cf         = 267
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_frozen_cf         = 268
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_sfc_zonal_cur     = 269

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_sfc_merid_cur     = 270
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_qcf2              = 271
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_qrain             = 272
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_qgraup            = 273
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_top_ind_mean      = 274
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_Ti_Mean           = 274
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_top_ind_stddev    = 275
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_Ti_Sig            = 275
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_fexp              = 276
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_gamtot            = 277
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_zw                = 278
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_fsat              = 279

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_fwetl             = 280
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_sthzw             = 281
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_a_fsat            = 282
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_c_fsat            = 283
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_a_fwet            = 284
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_c_fwet            = 285

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_disturb_frac_veg_prev = 286
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_wood_prod_fast    = 287
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_wood_prod_med     = 288
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_wood_prod_slow    = 289

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_flake_depth       = 291
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_flake_fetch       = 292
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_flake_t_mean      = 293
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_flake_t_mxl       = 294
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_flake_t_ice       = 295
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_flake_h_mxl       = 296
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_flake_h_ice       = 297
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_flake_shape       = 298
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_flake_g_over_dt   = 299

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_user_anc_sing1    = 301
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_user_anc_sing20   = 320
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_user_anc_mult1    = 321

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_user_anc_mult20   = 340
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_tppsozone         = 341
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_deep_conv_flag    = 342
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_past_conv_precip  = 343
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_past_conv_depth   = 344
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_cca_dp            = 345
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_cca_md            = 346
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_cca_sh            = 347
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_total_precip      = 348

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_clim_biogenic_aero= 351
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_clim_delta_aero   = 371
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_snowdep_grd_tile  = 376
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_snowpack_bk_dens  = 377

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_nsnow_layrs_tiles = 380
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_snow_laythk_tiles = 381
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_snow_ice_tile     = 382
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_snow_liq_tile     = 383
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_snow_T_tile       = 384
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_snow_laydns_tiles = 385
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_snow_grnsiz_tiles = 386
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_etadot            = 387
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_thetavd           = 388
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dry_rho           = 389

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_psiw_surf         = 390
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_psiw_lid          = 397
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mv                = 391
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mcl               = 392
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mcf               = 393
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mr                = 394
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mgr               = 395
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mcf2              = 396
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_exner_surf        = 398

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_p                 = 407
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_pstar             = 409
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ice_conc_cat      = 413
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ice_thick_cat     = 414
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ice_temp_cat      = 415
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ice_snow_depth    = 416
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dust_parent_clay  = 418
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dust_parent_silt  = 419

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dust_parent_sand  = 420
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dust_soil_mf1     = 421
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dust_soil_mf2     = 422
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dust_soil_mf3     = 423
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dust_soil_mf4     = 424
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dust_soil_mf5     = 425
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dust_soil_mf6     = 426
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_soil_massfrac6    = 426
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_pond_frac_cat     = 428
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_pond_depth_cat    = 429

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dust1_mmr         = 431
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dust2_mmr         = 432
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dust3_mmr         = 433
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dust4_mmr         = 434
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dust5_mmr         = 435
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dust6_mmr         = 436

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ice_surf_cond_cat = 440
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ice_surf_temp_cat = 441

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_soilnitro_dpm     = 442
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_soilnitro_rpm     = 443
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_soilnitro_bio     = 444
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_soilnitro_hum     = 445
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_soil_inorgnit     = 446
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_nitrogen_deposition = 447

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_crop_frac         = 448
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_pasture_frac      = 458

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_soilcarb_dpm      = 466
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_soilcarb_rpm      = 467
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_soilcarb_bio      = 468
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_soilcarb_hum      = 469

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ozone_tracer      = 480
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_o3_prod_loss      = 481
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_o3_p_l_vmr        = 482
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_o3_vmr            = 483
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_o3_p_l_temp       = 484
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_o3_temp           = 485
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_o3_p_l_colo3      = 486
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_o3_colo3          = 487

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dctemp_tile       = 490
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dctemp_ssi        = 491
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_tm_trans          = 492
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ddmfx             = 493
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_urbhgt            = 494
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_urbhwr            = 495
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_urbwrr            = 496
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_urbdisp           = 497
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_urbztm            = 498
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_urbalbwl          = 499

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_urbalbrd          = 500
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_urbemisw          = 501
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_urbemisr          = 502
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_land_frac         = 505
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_tstar_land        = 506
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_tstar_sea         = 507
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_tstar_sice        = 508
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_albedo_sice       = 509
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_albedo_land       = 510

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_u10_cpl           = 515
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_v10_cpl           = 516
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_charnock_cpl      = 517

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ux_ccp            = 569
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_uy_ccp            = 570
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_um_ccp            = 571
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_g_ccp             = 572
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_h_ccp             = 573
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_riso_ccp          = 574
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_rdir_ccp          = 575
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_tsurf_elev_surft  = 576

! PV-tracers
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_rad           = 577
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_sw            = 578
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_lw            = 579
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_mic           = 580
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_gwd           = 581
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_ph1           = 582
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_conv          = 583
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_bl            = 584
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_stph          = 585
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_cld           = 586
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_iau           = 587
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_nud           = 588
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_tot           = 589
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dEPS_I            = 590
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_sol           = 591
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_mass          = 592
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_0             = 593

! Stochastic Physics
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_bl_pert_rand_fld  = 595
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_bl_pert_flag      = 596

! More PV-tracers, diab-friction split
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_conv_d        = 620
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_conv_f        = 621
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_bl_d          = 622
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_bl_f          = 623
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dPV_PC2c          = 624
! Theta tracers
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dtheta_0          = 600
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dtheta_bl         = 601
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dtheta_bl_mix     = 602
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dtheta_bl_LH      = 603
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dtheta_conv       = 604
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dtheta_mic        = 605
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dtheta_rad        = 606
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dtheta_SW         = 607
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dtheta_LW         = 608
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dtheta_slow       = 609
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dtheta_cld        = 610
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_dtheta_PC2c       = 611

! INFERNO Ignition Ancillaries
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_flash_rate_ancil  = 626
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_pop_den_ancil     = 627
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_wealth_index_ancil= 628

! Irrigation
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_sthu_irr          = 630
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_frac_irr          = 631

! New
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_qt                = 16207

! Chemistry fields
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_o3                = 34001
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_no                = 34002
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_no3               = 34003
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_n2o5              = 34005
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ho2no2            = 34006
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_hono2             = 34007
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_h2o2              = 34008
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ch4               = 34009
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_co                = 34010
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_hcho              = 34011
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_meooh             = 34012
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_hono              = 34013
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_c2h6              = 34014
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_etooh             = 34015
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mecho             = 34016
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_pan               = 34017
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_c3h8              = 34018
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_n_prooh           = 34019
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_i_prooh           = 34020
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_etcho             = 34021
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_me2co             = 34022
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mecoch2ooh        = 34023
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ppan              = 34024
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_meono2            = 34025
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_c5h8              = 34027
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_isooh             = 34028
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ison              = 34029
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_macr              = 34030
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_macrooh           = 34031
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mpan              = 34032
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_hacet             = 34033
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mgly              = 34034
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_nald              = 34035
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_hcooh             = 34036
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_meco3h            = 34037
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_meco2h            = 34038
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_cl                = 34041
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_clo               = 34042
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_cl2o2             = 34043
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_oclo              = 34044
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_br                = 34045
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_brcl              = 34047
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_brono2            = 34048
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_n2o               = 34049
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_hocl              = 34051
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_hbr               = 34052
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_hobr              = 34053
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_clono2            = 34054
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_cfcl3             = 34055
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_cf2cl2            = 34056
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mebr              = 34057
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_n                 = 34058
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_o3p               = 34059
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_h2                = 34070
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_msa               = 34074
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_nh3               = 34076
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_cs2               = 34077
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_csul              = 34078
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_h2s               = 34079
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_h                 = 34080
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_oh                = 34081
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_ho2               = 34082
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_meoo              = 34083
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_etoo              = 34084
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_meco3             = 34085
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_n_proo            = 34086
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_i_proo            = 34087
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_etco3             = 34088
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_mecoch2oo         = 34089
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_meoh              = 34090
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_so3               = 34094
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_lumped_n          = 34098
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_lumped_br         = 34099
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_lumped_cl         = 34100
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_passive_o3        = 34149
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_age_of_air        = 34150

INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_hcl               = 34992
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_bro               = 34994
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_no2               = 34996
INTEGER(KIND=int64), PUBLIC, PARAMETER :: stashcode_o1d               = 34997


! STASHmaster element codes
INTEGER(KIND=int64), PUBLIC, PARAMETER :: model = 1
INTEGER(KIND=int64), PUBLIC, PARAMETER :: section = 2
INTEGER(KIND=int64), PUBLIC, PARAMETER :: item = 3
INTEGER(KIND=int64), PUBLIC, PARAMETER :: name = 4
INTEGER(KIND=int64), PUBLIC, PARAMETER :: space = 5
INTEGER(KIND=int64), PUBLIC, PARAMETER :: point = 6
INTEGER(KIND=int64), PUBLIC, PARAMETER :: time = 7
INTEGER(KIND=int64), PUBLIC, PARAMETER :: grid = 8
INTEGER(KIND=int64), PUBLIC, PARAMETER :: levelt = 9
INTEGER(KIND=int64), PUBLIC, PARAMETER :: levelf = 10
INTEGER(KIND=int64), PUBLIC, PARAMETER :: levell = 11
INTEGER(KIND=int64), PUBLIC, PARAMETER :: pseudt = 12
INTEGER(KIND=int64), PUBLIC, PARAMETER :: pseudf = 13
INTEGER(KIND=int64), PUBLIC, PARAMETER :: pseudl = 14
INTEGER(KIND=int64), PUBLIC, PARAMETER :: levcom = 15
INTEGER(KIND=int64), PUBLIC, PARAMETER :: option = 16
INTEGER(KIND=int64), PUBLIC, PARAMETER :: version_mask = 17
INTEGER(KIND=int64), PUBLIC, PARAMETER :: halo = 18
INTEGER(KIND=int64), PUBLIC, PARAMETER :: datat = 19
INTEGER(KIND=int64), PUBLIC, PARAMETER :: dumpp = 20
INTEGER(KIND=int64), PUBLIC, PARAMETER :: packing_codes = 21
INTEGER(KIND=int64), PUBLIC, PARAMETER :: rotate = 22
INTEGER(KIND=int64), PUBLIC, PARAMETER :: ppfc = 23
INTEGER(KIND=int64), PUBLIC, PARAMETER :: user = 24
INTEGER(KIND=int64), PUBLIC, PARAMETER :: lbvc = 25
INTEGER(KIND=int64), PUBLIC, PARAMETER :: blev = 26
INTEGER(KIND=int64), PUBLIC, PARAMETER :: tlev = 27
INTEGER(KIND=int64), PUBLIC, PARAMETER :: rblevv = 28
INTEGER(KIND=int64), PUBLIC, PARAMETER :: cfll = 29
INTEGER(KIND=int64), PUBLIC, PARAMETER :: cfff = 30

PRIVATE

PUBLIC :: lfricinp_read_stashmaster, get_stashmaster_item

CONTAINS

!-------------------------------------------------------------------

SUBROUTINE lfricinp_read_stashmaster(stashmaster_path)
! Description:
!  Read in STASHmaster using shumlib and check return code
!  of shumlib function.

IMPLICIT NONE
! Arguments
CHARACTER(LEN=*), INTENT(IN) :: stashmaster_path

! Error handling
INTEGER(KIND=int64) :: status
CHARACTER(LEN=1024) :: message

status = f_shum_read_stashmaster(stashmaster_path, stashmaster,  &
                                 message)
IF (status /= 0_int64) THEN
  CALL log_event(message, LOG_LEVEL_ERROR)
END IF

END SUBROUTINE lfricinp_read_stashmaster

!-------------------------------------------------------------------

FUNCTION get_stashmaster_item(stashcode, item_id) RESULT(item)
! Currently only coded for 64-bit scalar integers
IMPLICIT NONE

INTEGER(KIND=int64), INTENT(IN) :: stashcode
INTEGER(KIND=int64), INTENT(IN) :: item_id
! Result
INTEGER(KIND=int64) :: item

item = imdi

! Check that the STASHmaster record exists.
IF (.NOT. ASSOCIATED(stashmaster(stashcode) % record)) THEN
  WRITE(log_scratch_space, '(A,I0)') &
       "Unassociated STASHmaster record for stashcode ", stashcode
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

SELECT CASE(item_id)
CASE(grid)
  item = stashmaster(stashcode) % record % grid
CASE(levelt)
  item = stashmaster(stashcode) % record % levelt
CASE(levelf)
  item = stashmaster(stashcode) % record % levelf
CASE(levell)
  item = stashmaster(stashcode) % record % levell
CASE(pseudt)
  item = stashmaster(stashcode) % record % pseudt
CASE(pseudf)
  item = stashmaster(stashcode) % record % pseudf
CASE(pseudl)
  item = stashmaster(stashcode) % record % pseudl
CASE(ppfc)
  item = stashmaster(stashcode) % record % ppfc
CASE(lbvc)
  item = stashmaster(stashcode) % record % lbvc
CASE(datat)
  item = stashmaster(stashcode) % record % datat
CASE(cfff)
  item = stashmaster(stashcode) % record % cfff
CASE(cfll)
  item = stashmaster(stashcode) % record % cfll
CASE DEFAULT
  WRITE(log_scratch_space, '(A,I0,A)') &
       "STASHmaster item id ", item_id, " not supported"
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END SELECT

END FUNCTION get_stashmaster_item

END MODULE lfricinp_stashmaster_mod
