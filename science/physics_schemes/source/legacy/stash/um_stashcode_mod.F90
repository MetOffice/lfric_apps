!
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  List of stashcode magic numbers

module um_stashcode_mod

! Description:
!   Stash code definitions used in the RCF, ancillary update and PWS codes
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Stash
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 v6 programming standards.

implicit none

! Section numbers
integer, parameter :: stashcode_prog_sec          =   0
integer, parameter :: stashcode_bl_sec            =   3
integer, parameter :: stashcode_ls_sec            =   4
integer, parameter :: stashcode_conv_sec          =   5
integer, parameter :: stashcode_gwd_sec           =   6
integer, parameter :: stashcode_proc_dyn_sec      =  15
integer, parameter :: stashcode_proc_phys_sec     =  16
integer, parameter :: stashcode_aerosol_sec       =  17
integer, parameter :: stashcode_pws_sec           =  20
integer, parameter :: stashcode_lbc_input_sec     =  31
integer, parameter :: stashcode_lbc_output_sec    =  32
integer, parameter :: stashcode_tracer_sec        =  33
integer, parameter :: stashcode_ukca_sec          =  34
integer, parameter :: stashcode_tracer_lbc_sec    =  36
integer, parameter :: stashcode_ukca_lbc_sec      =  37
integer, parameter :: stashcode_glomap_sec        =  38
integer, parameter :: stashcode_ukca_chem_diag    =  50
integer, parameter :: stashcode_glomap_clim_sec   =  54

! Stashcodes used in the reconfiguration and/or ancillary updating

integer, parameter :: stashcode_u                 =   2
integer, parameter :: stashcode_v                 =   3
integer, parameter :: stashcode_theta             =   4
integer, parameter :: stashcode_orog_x_grad       =   5
integer, parameter :: stashcode_orog_y_grad       =   6
integer, parameter :: stashcode_unfilt_orog       =   7
integer, parameter :: stashcode_soil_moist        =   9

integer, parameter :: stashcode_q                 =  10
integer, parameter :: stashcode_qcf               =  12
integer, parameter :: stashcode_cca               =  13
integer, parameter :: stashcode_ccb               =  14
integer, parameter :: stashcode_cct               =  15
integer, parameter :: stashcode_cc_lwp            =  16
integer, parameter :: stashcode_sil_orog_rough    =  17
integer, parameter :: stashcode_hlf_pk_to_trf_ht  =  18

integer, parameter :: stashcode_soil_temp         =  20
integer, parameter :: stashcode_lcbase            =  21
integer, parameter :: stashcode_mean_canopyw      =  22
integer, parameter :: stashcode_snow_amount       =  23
integer, parameter :: stashcode_mean_snow         =  23
integer, parameter :: stashcode_surftemp          =  24
integer, parameter :: stashcode_tstar             =  24
integer, parameter :: stashcode_bl_depth          =  25
integer, parameter :: stashcode_rough_length      =  26
integer, parameter :: stashcode_z0                =  26
integer, parameter :: stashcode_snow_edge         =  27
integer, parameter :: stashcode_surf_z_curr       =  28
integer, parameter :: stashcode_surf_m_curr       =  29

integer, parameter :: stashcode_lsm               =  30
integer, parameter :: stashcode_icefrac           =  31
integer, parameter :: stashcode_icethick          =  32
integer, parameter :: stashcode_orog              =  33
integer, parameter :: stashcode_orog_var          =  34
integer, parameter :: stashcode_orog_gdxx         =  35
integer, parameter :: stashcode_orog_gdxy         =  36
integer, parameter :: stashcode_orog_gdyy         =  37
integer, parameter :: stashcode_ice_edge_ancil    =  38
integer, parameter :: stashcode_ice_edge_inancil  =  38
integer, parameter :: stashcode_tstar_anom        =  39

integer, parameter :: stashcode_vol_smc_wilt      =  40
integer, parameter :: stashcode_vol_smc_cri       =  41
integer, parameter :: stashcode_vol_smc_sat       =  43
integer, parameter :: stashcode_Ksat              =  44
integer, parameter :: stashcode_thermal_capacity  =  46
integer, parameter :: stashcode_thermal_conduct   =  47
integer, parameter :: stashcode_soil_suction      =  48
integer, parameter :: stashcode_sea_ice_temp      =  49

integer, parameter :: stashcode_veg_frac          =  50
integer, parameter :: stashcode_total_aero_emiss  =  57
integer, parameter :: stashcode_SO2_emiss         =  58
integer, parameter :: stashcode_dimethyl_sul_emiss=  59

integer, parameter :: stashcode_ozone             =  60

integer, parameter :: stashcode_orog_f1           =  61
integer, parameter :: stashcode_orog_f2           =  62
integer, parameter :: stashcode_orog_f3           =  63
integer, parameter :: stashcode_orog_amp          =  64

integer, parameter :: stashcode_e_trb             =  70
integer, parameter :: stashcode_tsq_trb           =  71
integer, parameter :: stashcode_qsq_trb           =  72
integer, parameter :: stashcode_cov_trb           =  73
integer, parameter :: stashcode_zhpar_shcu        =  74

integer, parameter :: stashcode_cloud_number      =  75
integer, parameter :: stashcode_rain_number       =  76
integer, parameter :: stashcode_rain_3mom         =  77
integer, parameter :: stashcode_ice_number        =  78
integer, parameter :: stashcode_snow_number       =  79
integer, parameter :: stashcode_snow_3mom         =  80
integer, parameter :: stashcode_graup_number      =  81
integer, parameter :: stashcode_graup_3mom        =  82

integer, parameter :: stashcode_activesol_liquid  =  83
integer, parameter :: stashcode_activesol_rain    =  84
integer, parameter :: stashcode_active_insol_ice  =  85
integer, parameter :: stashcode_active_sol_ice    =  86
integer, parameter :: stashcode_active_insol_liq  =  87
integer, parameter :: stashcode_active_sol_num    =  88
integer, parameter :: stashcode_active_insol_num  =  89

integer, parameter :: stashcode_total_aero        =  90
integer, parameter :: stashcode_flash_pot         =  91
integer, parameter :: stashcode_runoff_coast_out  =  93
integer, parameter :: stashcode_snow_on_ice       =  95
integer, parameter :: stashcode_ocnsrf_chlorophyll=  96
integer, parameter :: stashcode_chlorophyll       =  96
integer, parameter :: stashcode_z0m_soil          =  97

integer, parameter :: stashcode_blwvariance       =  99

integer, parameter :: stashcode_so2               = 101
integer, parameter :: stashcode_dms               = 102
integer, parameter :: stashcode_mmr_so4_aitken    = 103
integer, parameter :: stashcode_mmr_so4_accum     = 104
integer, parameter :: stashcode_mmr_so4_diss      = 105
integer, parameter :: stashcode_mmr_nh3           = 107
integer, parameter :: stashcode_mmr_bc_fr         = 108
integer, parameter :: stashcode_mmr_bc_ag         = 109
integer, parameter :: stashcode_mmr_bc_cl         = 110
integer, parameter :: stashcode_mmr_smoke_fr      = 111
integer, parameter :: stashcode_mmr_smoke_ag      = 112
integer, parameter :: stashcode_mmr_smoke_cl      = 113
integer, parameter :: stashcode_mmr_ocff_fr       = 114
integer, parameter :: stashcode_mmr_ocff_ag       = 115
integer, parameter :: stashcode_mmr_ocff_cl       = 116
integer, parameter :: stashcode_mmr_nitr_acc      = 117
integer, parameter :: stashcode_mmr_nitr_diss     = 118

integer, parameter :: stashcode_biom_elev_em_h1   = 119
integer, parameter :: stashcode_biom_elev_em_h2   = 120
integer, parameter :: stashcode_3d_nat_so2_em     = 121
integer, parameter :: stashcode_3d_oh_conc        = 122
integer, parameter :: stashcode_3d_ho2_conc       = 123
integer, parameter :: stashcode_3dh2o2_mixrat     = 124
integer, parameter :: stashcode_3d_ozone_mixrat   = 125
integer, parameter :: stashcode_hi_SO2_emiss_emiss= 126
integer, parameter :: stashcode_ammonia_gas_emiss = 127
integer, parameter :: stashcode_soot_surf         = 128
integer, parameter :: stashcode_soot_hi_lev       = 129

integer, parameter :: stashcode_biom_surf_em      = 130
integer, parameter :: stashcode_biom_elev_em      = 131
integer, parameter :: stashcode_dms_conc_sea      = 132
integer, parameter :: stashcode_dms_conc_sw       = 132
integer, parameter :: stashcode_ocff_surf_emiss   = 134
integer, parameter :: stashcode_ocff_hilev_emiss  = 135

integer, parameter :: stashcode_w                 = 150
integer, parameter :: stashcode_riv_sequence      = 151
integer, parameter :: stashcode_riv_direction     = 152
integer, parameter :: stashcode_riv_storage       = 153
integer, parameter :: stashcode_riv_number        = 154
integer, parameter :: stashcode_riv_tot_surfroff  = 155
integer, parameter :: stashcode_riv_tot_subroff   = 156
integer, parameter :: stashcode_riv_iarea         = 160
integer, parameter :: stashcode_riv_slope         = 161
integer, parameter :: stashcode_riv_flowobs1      = 162
integer, parameter :: stashcode_riv_inext         = 163
integer, parameter :: stashcode_riv_jnext         = 164
integer, parameter :: stashcode_riv_land          = 165
integer, parameter :: stashcode_riv_sfcstorage    = 166
integer, parameter :: stashcode_riv_substorage    = 167
integer, parameter :: stashcode_riv_flowin        = 168
integer, parameter :: stashcode_riv_bflowin       = 169
integer, parameter :: stashcode_ice_subl_cat      = 182
integer, parameter :: stashcode_iceberg_calving   = 190
integer, parameter :: stashcode_sstfrz            = 194
integer, parameter :: stashcode_tstar_ice_cat_cpl = 195
integer, parameter :: stashcode_riv_outflow_cpl   = 198
integer, parameter :: stashcode_penabs_cpl        = 200

integer, parameter :: stashcode_u_compnt_pert     = 202
integer, parameter :: stashcode_v_compnt_pert     = 203
integer, parameter :: stashcode_clapp_hb          = 207
integer, parameter :: stashcode_3d_cca            = 211
integer, parameter :: stashcode_3d_ccw            = 212
integer, parameter :: stashcode_can_conduct       = 213
integer, parameter :: stashcode_unfrozen_soil     = 214
integer, parameter :: stashcode_frozen_soil       = 215
integer, parameter :: stashcode_frac_surf_type    = 216
integer, parameter :: stashcode_lai               = 217
integer, parameter :: stashcode_canopy_height     = 218
integer, parameter :: stashcode_disturb_frac_veg  = 219

integer, parameter :: stashcode_snw_free_alb_bs   = 220
integer, parameter :: stashcode_soil_carbon_cont  = 223
integer, parameter :: stashcode_npp_pft_acc       = 224
integer, parameter :: stashcode_g_lf_pft_acc      = 225
integer, parameter :: stashcode_g_ph_lf_pft_acc   = 226
integer, parameter :: stashcode_rsp_w_pft_acc     = 227
integer, parameter :: stashcode_resp_s_acc        = 228
integer, parameter :: stashcode_can_water_tile    = 229

integer, parameter :: stashcode_catch_tile        = 230
integer, parameter :: stashcode_rgrain            = 231
integer, parameter :: stashcode_tstar_tile        = 233
integer, parameter :: stashcode_tsurf_elev_surft  = 576
integer, parameter :: stashcode_z0_tile           = 234
integer, parameter :: stashcode_infil_max_tile    = 236
integer, parameter :: stashcode_sw_down_tile      = 237
integer, parameter :: stashcode_sw_down           = 238
integer, parameter :: stashcode_lw_up_diff        = 239

integer, parameter :: stashcode_snow_tile         = 240
integer, parameter :: stashcode_catch_snow        = 241
integer, parameter :: stashcode_snow_grnd         = 242
integer, parameter :: stashcode_surf_sw_alb       = 243
integer, parameter :: stashcode_surf_vis_alb      = 244
integer, parameter :: stashcode_surf_nir_alb      = 245
integer, parameter :: stashcode_z0h_tile          = 246

integer, parameter :: stashcode_CO2_surf_emiss    = 251
integer, parameter :: stashcode_rho               = 253
integer, parameter :: stashcode_qcl               = 254
integer, parameter :: stashcode_exner             = 255
integer, parameter :: stashcode_u_adv             = 256
integer, parameter :: stashcode_v_adv             = 257
integer, parameter :: stashcode_w_adv             = 258
integer, parameter :: stashcode_n_turb_mixlvs     = 259

integer, parameter :: stashcode_lvl_bse_dp_sc     = 260
integer, parameter :: stashcode_lvl_top_dp_sc     = 261
integer, parameter :: stashcode_bl_conv_flag      = 262
integer, parameter :: stashcode_turb_temp         = 263
integer, parameter :: stashcode_turb_humid        = 264
integer, parameter :: stashcode_area_cf           = 265
integer, parameter :: stashcode_bulk_cf           = 266
integer, parameter :: stashcode_liquid_cf         = 267
integer, parameter :: stashcode_frozen_cf         = 268
integer, parameter :: stashcode_sfc_zonal_cur     = 269

integer, parameter :: stashcode_sfc_merid_cur     = 270
integer, parameter :: stashcode_qcf2              = 271
integer, parameter :: stashcode_qrain             = 272
integer, parameter :: stashcode_qgraup            = 273
integer, parameter :: stashcode_top_ind_mean      = 274
integer, parameter :: stashcode_Ti_Mean           = 274
integer, parameter :: stashcode_top_ind_stddev    = 275
integer, parameter :: stashcode_Ti_Sig            = 275
integer, parameter :: stashcode_fexp              = 276
integer, parameter :: stashcode_gamtot            = 277
integer, parameter :: stashcode_zw                = 278
integer, parameter :: stashcode_fsat              = 279

integer, parameter :: stashcode_fwetl             = 280
integer, parameter :: stashcode_sthzw             = 281
integer, parameter :: stashcode_a_fsat            = 282
integer, parameter :: stashcode_c_fsat            = 283
integer, parameter :: stashcode_a_fwet            = 284
integer, parameter :: stashcode_c_fwet            = 285

integer, parameter :: stashcode_disturb_frac_veg_prev = 286
integer, parameter :: stashcode_wood_prod_fast    = 287
integer, parameter :: stashcode_wood_prod_med     = 288
integer, parameter :: stashcode_wood_prod_slow    = 289
integer, parameter :: stashcode_acc_lake_evap     = 290
integer, parameter :: stashcode_flake_depth       = 291
integer, parameter :: stashcode_flake_fetch       = 292
integer, parameter :: stashcode_flake_t_mean      = 293
integer, parameter :: stashcode_flake_t_mxl       = 294
integer, parameter :: stashcode_flake_t_ice       = 295
integer, parameter :: stashcode_flake_h_mxl       = 296
integer, parameter :: stashcode_flake_h_ice       = 297
integer, parameter :: stashcode_flake_shape       = 298
integer, parameter :: stashcode_flake_g_over_dt   = 299

integer, parameter :: stashcode_user_anc_sing1    = 301
integer, parameter :: stashcode_user_anc_sing20   = 320
integer, parameter :: stashcode_user_anc_mult1    = 321

integer, parameter :: stashcode_user_anc_mult20   = 340
integer, parameter :: stashcode_tppsozone         = 341
integer, parameter :: stashcode_deep_conv_flag    = 342
integer, parameter :: stashcode_past_conv_precip  = 343
integer, parameter :: stashcode_past_conv_depth   = 344
integer, parameter :: stashcode_cca_dp            = 345
integer, parameter :: stashcode_cca_md            = 346
integer, parameter :: stashcode_cca_sh            = 347
integer, parameter :: stashcode_total_precip      = 348

integer, parameter :: stashcode_clim_biogenic_aero= 351
integer, parameter :: stashcode_clim_delta_aero   = 371
integer, parameter :: stashcode_snowdep_grd_tile  = 376
integer, parameter :: stashcode_snowpack_bk_dens  = 377

integer, parameter :: stashcode_nsnow_layrs_tiles = 380
integer, parameter :: stashcode_snow_laythk_tiles = 381
integer, parameter :: stashcode_snow_ice_tile     = 382
integer, parameter :: stashcode_snow_liq_tile     = 383
integer, parameter :: stashcode_snow_T_tile       = 384
integer, parameter :: stashcode_snow_laydns_tiles = 385
integer, parameter :: stashcode_snow_grnsiz_tiles = 386

integer, parameter :: stashcode_p                 = 407
integer, parameter :: stashcode_pstar             = 409
integer, parameter :: stashcode_ice_conc_cat      = 413
integer, parameter :: stashcode_ice_thick_cat     = 414
integer, parameter :: stashcode_ice_temp_cat      = 415
integer, parameter :: stashcode_ice_snow_depth    = 416
integer, parameter :: stashcode_dust_parent_clay  = 418
integer, parameter :: stashcode_dust_parent_silt  = 419

integer, parameter :: stashcode_dust_parent_sand  = 420
integer, parameter :: stashcode_dust_soil_mf1     = 421
integer, parameter :: stashcode_dust_soil_mf2     = 422
integer, parameter :: stashcode_dust_soil_mf3     = 423
integer, parameter :: stashcode_dust_soil_mf4     = 424
integer, parameter :: stashcode_dust_soil_mf5     = 425
integer, parameter :: stashcode_dust_soil_mf6     = 426
integer, parameter :: stashcode_soil_massfrac6    = 426
integer, parameter :: stashcode_dust_psti         = 427
integer, parameter :: stashcode_pond_frac_cat     = 428
integer, parameter :: stashcode_pond_depth_cat    = 429

integer, parameter :: stashcode_dust1_mmr         = 431
integer, parameter :: stashcode_dust2_mmr         = 432
integer, parameter :: stashcode_dust3_mmr         = 433
integer, parameter :: stashcode_dust4_mmr         = 434
integer, parameter :: stashcode_dust5_mmr         = 435
integer, parameter :: stashcode_dust6_mmr         = 436

integer, parameter :: stashcode_ice_surf_cond_cat = 440
integer, parameter :: stashcode_ice_surf_temp_cat = 441

integer, parameter :: stashcode_soilnitro_dpm     = 442
integer, parameter :: stashcode_soilnitro_rpm     = 443
integer, parameter :: stashcode_soilnitro_bio     = 444
integer, parameter :: stashcode_soilnitro_hum     = 445
integer, parameter :: stashcode_soil_inorgnit     = 446
integer, parameter :: stashcode_nitrogen_deposition = 447

integer, parameter :: stashcode_crop_frac         = 448
integer, parameter :: stashcode_pasture_frac      = 458

integer, parameter :: stashcode_soilcarb_dpm      = 466
integer, parameter :: stashcode_soilcarb_rpm      = 467
integer, parameter :: stashcode_soilcarb_bio      = 468
integer, parameter :: stashcode_soilcarb_hum      = 469

integer, parameter :: stashcode_soilcarblyr_dpm   = 730
integer, parameter :: stashcode_soilcarblyr_rpm   = 731
integer, parameter :: stashcode_soilcarblyr_bio   = 732
integer, parameter :: stashcode_soilcarblyr_hum   = 733

integer, parameter :: stashcode_soilnitrolyr_dpm  = 734
integer, parameter :: stashcode_soilnitrolyr_rpm  = 735
integer, parameter :: stashcode_soilnitrolyr_bio  = 736
integer, parameter :: stashcode_soilnitrolyr_hum  = 737
integer, parameter :: stashcode_soil_inorgnit_lyr = 738
integer, parameter :: stashcode_soil_inorgnit_availpft = 743

integer, parameter :: stashcode_resplyr_dpm       = 739
integer, parameter :: stashcode_resplyr_rpm       = 740
integer, parameter :: stashcode_resplyr_bio       = 741
integer, parameter :: stashcode_resplyr_hum       = 742

integer, parameter :: stashcode_ozone_tracer      = 480
integer, parameter :: stashcode_o3_prod_loss      = 481
integer, parameter :: stashcode_o3_p_l_vmr        = 482
integer, parameter :: stashcode_o3_vmr            = 483
integer, parameter :: stashcode_o3_p_l_temp       = 484
integer, parameter :: stashcode_o3_temp           = 485
integer, parameter :: stashcode_o3_p_l_colo3      = 486
integer, parameter :: stashcode_o3_colo3          = 487

integer, parameter :: stashcode_dctemp_tile       = 490
integer, parameter :: stashcode_dctemp_ssi        = 491
integer, parameter :: stashcode_tm_trans          = 492
integer, parameter :: stashcode_ddmfx             = 493
integer, parameter :: stashcode_urbhgt            = 494
integer, parameter :: stashcode_urbhwr            = 495
integer, parameter :: stashcode_urbwrr            = 496
integer, parameter :: stashcode_urbdisp           = 497
integer, parameter :: stashcode_urbztm            = 498
integer, parameter :: stashcode_urbalbwl          = 499

integer, parameter :: stashcode_urbalbrd          = 500
integer, parameter :: stashcode_urbemisw          = 501
integer, parameter :: stashcode_urbemisr          = 502
integer, parameter :: stashcode_land_frac         = 505
integer, parameter :: stashcode_tstar_land        = 506
integer, parameter :: stashcode_tstar_sea         = 507
integer, parameter :: stashcode_tstar_sice        = 508
integer, parameter :: stashcode_albedo_sice       = 509
integer, parameter :: stashcode_albedo_land       = 510
integer, parameter :: stashcode_inlandatm         = 511
integer, parameter :: stashcode_z0m_sice_fmd      = 513
integer, parameter :: stashcode_z0m_sice_skin     = 514
integer, parameter :: stashcode_blendht           = 680

integer, parameter :: stashcode_u10_cpl           = 515
integer, parameter :: stashcode_v10_cpl           = 516
integer, parameter :: stashcode_charnock_cpl      = 517
integer, parameter :: stashcode_rhoair_cpl        = 518

integer, parameter :: stashcode_ux_ccp            = 569
integer, parameter :: stashcode_uy_ccp            = 570
integer, parameter :: stashcode_um_ccp            = 571
integer, parameter :: stashcode_g_ccp             = 572
integer, parameter :: stashcode_h_ccp             = 573
integer, parameter :: stashcode_riso_ccp          = 574
integer, parameter :: stashcode_rdir_ccp          = 575

! PV-tracers
integer, parameter :: stashcode_dPV_rad           = 577
integer, parameter :: stashcode_dPV_sw            = 578
integer, parameter :: stashcode_dPV_lw            = 579
integer, parameter :: stashcode_dPV_mic           = 580
integer, parameter :: stashcode_dPV_gwd           = 581
integer, parameter :: stashcode_dPV_ph1           = 582
integer, parameter :: stashcode_dPV_conv          = 583
integer, parameter :: stashcode_dPV_bl            = 584
integer, parameter :: stashcode_dPV_stph          = 585
integer, parameter :: stashcode_dPV_cld           = 586
integer, parameter :: stashcode_dPV_iau           = 587
integer, parameter :: stashcode_dPV_nud           = 588
integer, parameter :: stashcode_dPV_tot           = 589
integer, parameter :: stashcode_dEPS_I            = 590
integer, parameter :: stashcode_dPV_sol           = 591
integer, parameter :: stashcode_dPV_mass          = 592
integer, parameter :: stashcode_dPV_0             = 593
! More PV-tracers, diab-friction split
integer, parameter :: stashcode_dPV_conv_d        = 620
integer, parameter :: stashcode_dPV_conv_f        = 621
integer, parameter :: stashcode_dPV_bl_d          = 622
integer, parameter :: stashcode_dPV_bl_f          = 623
integer, parameter :: stashcode_dPV_PC2c          = 624
! Theta tracers
integer, parameter :: stashcode_dtheta_0          = 600
integer, parameter :: stashcode_dtheta_bl         = 601
integer, parameter :: stashcode_dtheta_bl_mix     = 602
integer, parameter :: stashcode_dtheta_bl_LH      = 603
integer, parameter :: stashcode_dtheta_conv       = 604
integer, parameter :: stashcode_dtheta_mic        = 605
integer, parameter :: stashcode_dtheta_rad        = 606
integer, parameter :: stashcode_dtheta_SW         = 607
integer, parameter :: stashcode_dtheta_LW         = 608
integer, parameter :: stashcode_dtheta_slow       = 609
integer, parameter :: stashcode_dtheta_cld        = 610
integer, parameter :: stashcode_dtheta_PC2c       = 611

! Stochastic Physics
integer, parameter :: stashcode_bl_pert_rand_fld  = 595
integer, parameter :: stashcode_bl_pert_flag      = 596

! INFERNO Ignition Ancillaries
integer, parameter :: stashcode_flash_rate_ancil  = 626
integer, parameter :: stashcode_pop_den_ancil     = 627
integer, parameter :: stashcode_wealth_index_ancil= 628

! Irrigation
integer, parameter :: stashcode_sthu_irr          = 630
integer, parameter :: stashcode_frac_irr          = 631

! Wetness of surface
integer, parameter :: stashcode_surf_wetness      = 634

! Hybrid resolution model
integer, parameter :: stashcode_filter_land_hyb   = 635
integer, parameter :: stashcode_tke_activ_hyb     = 636

! Thermal acclimation of photosynthesis
integer, parameter :: stashcode_t_growth          = 637
integer, parameter :: stashcode_t_home            = 638

!----------------------------------------------------------
! Section 16 fields that may be reconfigured for VAR
!----------------------------------------------------------
integer, parameter :: stashcode_t                 = 16004
integer, parameter :: stashcode_qc                = 16206
integer, parameter :: stashcode_qT                = 16207

!----------------------------------------------------------
! Water Tracers - section 33
!----------------------------------------------------------
integer, parameter :: stashcode_q_wtrac                = 33151
integer, parameter :: stashcode_qcl_wtrac              = 33152
integer, parameter :: stashcode_qcf_wtrac              = 33153
integer, parameter :: stashcode_qcf2_wtrac             = 33154
integer, parameter :: stashcode_qr_wtrac               = 33155
integer, parameter :: stashcode_qgr_wtrac              = 33156
integer, parameter :: stashcode_mv_wtrac               = 33157
integer, parameter :: stashcode_mcl_wtrac              = 33158
integer, parameter :: stashcode_mcf_wtrac              = 33159
integer, parameter :: stashcode_mcf2_wtrac             = 33160
integer, parameter :: stashcode_mr_wtrac               = 33161
integer, parameter :: stashcode_mgr_wtrac              = 33162
integer, parameter :: stashcode_conv_prog_dq_wtrac     = 33163
integer, parameter :: stashcode_snow_tile_wtrac        = 33178
integer, parameter :: stashcode_snow_grnd_wtrac        = 33179
integer, parameter :: stashcode_snow_ice_tile_wtrac    = 33180
integer, parameter :: stashcode_snow_liq_tile_wtrac    = 33181
integer, parameter :: stashcode_can_water_tile_wtrac   = 33182
integer, parameter :: stashcode_unfrozen_soil_wtrac    = 33183
integer, parameter :: stashcode_frozen_soil_wtrac      = 33184
integer, parameter :: stashcode_soil_moist_wtrac       = 33185
integer, parameter :: stashcode_sthzw_wtrac            = 33186
integer, parameter :: stashcode_riv_tot_surfroff_wtrac = 33187
integer, parameter :: stashcode_riv_tot_subroff_wtrac  = 33188
integer, parameter :: stashcode_acc_lake_evap_wtrac    = 33189
integer, parameter :: stashcode_riv_storage_wtrac      = 33190
integer, parameter :: stashcode_inlandatm_wtrac        = 33191

!----------------------------------------------------------
! UKCA stashcodes - section 34
!----------------------------------------------------------
integer, parameter :: stashcode_NO2               = 34004
integer, parameter :: stashcode_CH4               = 34009
integer, parameter :: stashcode_CO                = 34010
integer, parameter :: stashcode_HCHO              = 34011
integer, parameter :: stashcode_O3                = 34001
integer, parameter :: stashcode_NO                = 34002
integer, parameter :: stashcode_HNO3              = 34007
integer, parameter :: stashcode_PAN               = 34017
integer, parameter :: stashcode_C2H6              = 34014
integer, parameter :: stashcode_C3H8              = 34018
! Aerosols
integer, parameter :: stashcode_nucsol_no         = 34101
integer, parameter :: stashcode_nucsol_so4        = 34102
integer, parameter :: stashcode_Aitsol_no         = 34103
integer, parameter :: stashcode_Aitsol_so4        = 34104
integer, parameter :: stashcode_Aitsol_bc         = 34105
integer, parameter :: stashcode_Aitsol_oc         = 34106
integer, parameter :: stashcode_accsol_no         = 34107
integer, parameter :: stashcode_accsol_so4        = 34108
integer, parameter :: stashcode_accsol_bc         = 34109
integer, parameter :: stashcode_accsol_oc         = 34110
integer, parameter :: stashcode_accsol_ss         = 34111
integer, parameter :: stashcode_accsol_du         = 34112
integer, parameter :: stashcode_corsol_no         = 34113
integer, parameter :: stashcode_corsol_so4        = 34114
integer, parameter :: stashcode_corsol_bc         = 34115
integer, parameter :: stashcode_corsol_oc         = 34116
integer, parameter :: stashcode_corsol_ss         = 34117
integer, parameter :: stashcode_corsol_du         = 34118
integer, parameter :: stashcode_Aitinsol_no       = 34119
integer, parameter :: stashcode_Aitinsol_bc       = 34120
integer, parameter :: stashcode_Aitinsol_oc       = 34121
integer, parameter :: stashcode_accinsol_no       = 34122
integer, parameter :: stashcode_accinsol_du       = 34123
integer, parameter :: stashcode_corinsol_no       = 34124
integer, parameter :: stashcode_corinsol_du       = 34125
integer, parameter :: stashcode_nucsol_oc         = 34126
integer, parameter :: stashcode_nucsol_so         = 34128
integer, parameter :: stashcode_Aitsol_so         = 34129
integer, parameter :: stashcode_accsol_so         = 34130
integer, parameter :: stashcode_corsol_so         = 34131
integer, parameter :: stashcode_nucsol_nh4        = 34132
integer, parameter :: stashcode_aitsol_nh4        = 34133
integer, parameter :: stashcode_accsol_nh4        = 34134
integer, parameter :: stashcode_corsol_nh4        = 34135
integer, parameter :: stashcode_nucsol_no3        = 34136
integer, parameter :: stashcode_aitsol_no3        = 34137
integer, parameter :: stashcode_accsol_no3        = 34138
integer, parameter :: stashcode_corsol_no3        = 34139
integer, parameter :: stashcode_accsol_nn         = 34250
integer, parameter :: stashcode_corsol_nn         = 34251
integer, parameter :: stashcode_supinsol_no       = 34262
integer, parameter :: stashcode_supinsol_du       = 34263
integer, parameter :: stashcode_aitsol_mp         = 34264
integer, parameter :: stashcode_accsol_mp         = 34265
integer, parameter :: stashcode_corsol_mp         = 34266
integer, parameter :: stashcode_aitinsol_mp       = 34267
integer, parameter :: stashcode_accinsol_mp       = 34268
integer, parameter :: stashcode_corinsol_mp       = 34269
integer, parameter :: stashcode_supinsol_mp       = 34270

! RADAER prognostics
integer, parameter :: stashcode_dryd_ait_sol      = 34921
integer, parameter :: stashcode_dryd_acc_sol      = 34922
integer, parameter :: stashcode_dryd_cor_sol      = 34923
integer, parameter :: stashcode_dryd_ait_insol    = 34924
integer, parameter :: stashcode_dryd_acc_insol    = 34925
integer, parameter :: stashcode_dryd_cor_insol    = 34926
integer, parameter :: stashcode_wetd_ait_sol      = 34927
integer, parameter :: stashcode_wetd_acc_sol      = 34928
integer, parameter :: stashcode_wetd_cor_sol      = 34929
integer, parameter :: stashcode_rho_ait_sol       = 34930
integer, parameter :: stashcode_rho_acc_sol       = 34931
integer, parameter :: stashcode_rho_cor_sol       = 34932
integer, parameter :: stashcode_rho_ait_insol     = 34933
integer, parameter :: stashcode_rho_acc_insol     = 34934
integer, parameter :: stashcode_rho_cor_insol     = 34935
integer, parameter :: stashcode_pvol_ait_su_sol   = 34936
integer, parameter :: stashcode_pvol_ait_bc_sol   = 34937
integer, parameter :: stashcode_pvol_ait_oc_sol   = 34938
integer, parameter :: stashcode_pvol_ait_no3_sol  = 34939
integer, parameter :: stashcode_pvol_ait_so_sol   = 34940
integer, parameter :: stashcode_pvol_ait_h2o_sol  = 34941
integer, parameter :: stashcode_pvol_acc_su_sol   = 34942
integer, parameter :: stashcode_pvol_acc_bc_sol   = 34943
integer, parameter :: stashcode_pvol_acc_oc_sol   = 34944
integer, parameter :: stashcode_pvol_acc_ss_sol   = 34945
integer, parameter :: stashcode_pvol_acc_no3_sol  = 34946
integer, parameter :: stashcode_pvol_acc_du_sol   = 34947
integer, parameter :: stashcode_pvol_acc_so_sol   = 34948
integer, parameter :: stashcode_pvol_acc_h2o_sol  = 34951
integer, parameter :: stashcode_pvol_cor_su_sol   = 34952
integer, parameter :: stashcode_pvol_cor_bc_sol   = 34953
integer, parameter :: stashcode_pvol_cor_oc_sol   = 34954
integer, parameter :: stashcode_pvol_cor_ss_sol   = 34955
integer, parameter :: stashcode_pvol_cor_no3_sol  = 34956
integer, parameter :: stashcode_pvol_cor_du_sol   = 34957
integer, parameter :: stashcode_pvol_cor_so_sol   = 34958
integer, parameter :: stashcode_pvol_cor_h2o_sol  = 34961
integer, parameter :: stashcode_pvol_ait_bc_insol = 34962
integer, parameter :: stashcode_pvol_ait_oc_insol = 34963
integer, parameter :: stashcode_pvol_acc_du_insol = 34964
integer, parameter :: stashcode_pvol_cor_du_insol = 34965
integer, parameter :: stashcode_pvol_ait_nh4_sol  = 34866
integer, parameter :: stashcode_pvol_acc_nh4_sol  = 34867
integer, parameter :: stashcode_pvol_cor_nh4_sol  = 34868
integer, parameter :: stashcode_pvol_acc_nn_sol   = 34869
integer, parameter :: stashcode_pvol_cor_nn_sol   = 34870
integer, parameter :: stashcode_dryd_sup_insol    = 34859
integer, parameter :: stashcode_rho_sup_insol     = 34860
integer, parameter :: stashcode_pvol_sup_du_insol = 34861
integer, parameter :: stashcode_pvol_ait_mp_sol   = 34852
integer, parameter :: stashcode_pvol_acc_mp_sol   = 34853
integer, parameter :: stashcode_pvol_cor_mp_sol   = 34854
integer, parameter :: stashcode_pvol_ait_mp_insol = 34855
integer, parameter :: stashcode_pvol_acc_mp_insol = 34856
integer, parameter :: stashcode_pvol_cor_mp_insol = 34857
integer, parameter :: stashcode_pvol_sup_mp_insol = 34858

! Needed in hybrid resolution model when running ACTIVATE in senior UM.
integer, parameter :: stashcode_dryd_nuc_sol    = 34862

! Non transported prognostics
integer, parameter :: stashcode_delta_q       = 34871
integer, parameter :: stashcode_n_activ_sum   = 34914
integer, parameter :: stashcode_surfarea      = 34966
integer, parameter :: stashcode_cdnc          = 34968
integer, parameter :: stashcode_ho2s          = 34969
integer, parameter :: stashcode_ohs           = 34970
integer, parameter :: stashcode_o1ds          = 34971
integer, parameter :: stashcode_o3ps          = 34972
integer, parameter :: stashcode_het_ho2       = 34973
integer, parameter :: stashcode_het_n2o5      = 34974
integer, parameter :: stashcode_tolp1         = 34975
integer, parameter :: stashcode_hoipo2        = 34976
integer, parameter :: stashcode_homvko2       = 34977
integer, parameter :: stashcode_memald1       = 34978
integer, parameter :: stashcode_oxyl1         = 34979
integer, parameter :: stashcode_hoc3h6o2      = 34980
integer, parameter :: stashcode_hoc2h4o2      = 34981
integer, parameter :: stashcode_meko2         = 34982
integer, parameter :: stashcode_mecoch2oo     = 34983
integer, parameter :: stashcode_mecoc2oo      = 34984
integer, parameter :: stashcode_etco3         = 34985
integer, parameter :: stashcode_iproo         = 34986
integer, parameter :: stashcode_sbuoo         = 34987
integer, parameter :: stashcode_nproo         = 34988
integer, parameter :: stashcode_meco3         = 34989
integer, parameter :: stashcode_etoo          = 34990
integer, parameter :: stashcode_meoo          = 34991
integer, parameter :: stashcode_hcl_unlmp     = 34992
integer, parameter :: stashcode_ho2_ntp       = 34993
integer, parameter :: stashcode_bro_unlmp     = 34994
integer, parameter :: stashcode_oh_ntp        = 34995
integer, parameter :: stashcode_no2_unlmp     = 34996
integer, parameter :: stashcode_o1d_ntp       = 34997
integer, parameter :: stashcode_o3p_ntp       = 34998

!----------------------------------------------------------
! LBC input stashcodes
!----------------------------------------------------------

integer, parameter :: stashcode_lbc_orog           = 31001
integer, parameter :: stashcode_lbc_u              = 31002
integer, parameter :: stashcode_lbc_v              = 31003
integer, parameter :: stashcode_lbc_w              = 31004
integer, parameter :: stashcode_lbc_density        = 31005
integer, parameter :: stashcode_lbc_theta          = 31006
integer, parameter :: stashcode_lbc_q              = 31007
integer, parameter :: stashcode_lbc_qcl            = 31008
integer, parameter :: stashcode_lbc_qcf            = 31009
integer, parameter :: stashcode_lbc_exner          = 31010

integer, parameter :: stashcode_lbc_u_adv          = 31011
integer, parameter :: stashcode_lbc_v_adv          = 31012
integer, parameter :: stashcode_lbc_w_adv          = 31013
integer, parameter :: stashcode_lbc_qcf2           = 31014
integer, parameter :: stashcode_lbc_qrain          = 31015
integer, parameter :: stashcode_lbc_qgraup         = 31016
integer, parameter :: stashcode_lbc_cf_bulk        = 31017
integer, parameter :: stashcode_lbc_cf_liquid      = 31018
integer, parameter :: stashcode_lbc_cf_frozen      = 31019
integer, parameter :: stashcode_lbc_murk           = 31020

integer, parameter :: stashcode_lbc_dust1_mmr      = 31023
integer, parameter :: stashcode_lbc_dust2_mmr      = 31024
integer, parameter :: stashcode_lbc_dust3_mmr      = 31025
integer, parameter :: stashcode_lbc_dust4_mmr      = 31026
integer, parameter :: stashcode_lbc_dust5_mmr      = 31027
integer, parameter :: stashcode_lbc_dust6_mmr      = 31028
integer, parameter :: stashcode_lbc_so2            = 31029
integer, parameter :: stashcode_lbc_dms            = 31030

integer, parameter :: stashcode_lbc_so4_aitken     = 31031
integer, parameter :: stashcode_lbc_so4_accu       = 31032
integer, parameter :: stashcode_lbc_so4_diss       = 31033
integer, parameter :: stashcode_lbc_nh3            = 31035
integer, parameter :: stashcode_lbc_soot_new       = 31036
integer, parameter :: stashcode_lbc_soot_agd       = 31037
integer, parameter :: stashcode_lbc_soot_cld       = 31038
integer, parameter :: stashcode_lbc_bmass_new      = 31039
integer, parameter :: stashcode_lbc_bmass_agd      = 31040

integer, parameter :: stashcode_lbc_bmass_cld      = 31041
integer, parameter :: stashcode_lbc_ocff_new       = 31042
integer, parameter :: stashcode_lbc_ocff_agd       = 31043
integer, parameter :: stashcode_lbc_ocff_cld       = 31044
integer, parameter :: stashcode_lbc_nitr_acc       = 31045
integer, parameter :: stashcode_lbc_nitr_diss      = 31046

!----------------------------------------------------------
! LBC input stashcodes - Free Tracers
!----------------------------------------------------------

integer, parameter :: stashcode_lbc_free_tracer_1   = 36001
integer, parameter :: stashcode_lbc_free_tracer_150 = 36150

!----------------------------------------------------------
! LBC input stashcodes - UKCA Tracers
!----------------------------------------------------------

integer, parameter :: stashcode_lbc_ukca_1         = 37001
integer, parameter :: stashcode_lbc_ukca_499       = 37499

!----------------------------------------------------------
! UKCA/GLOMAP diagnostic stashcodes - section 38
!----------------------------------------------------------

! CMIP6 diagnostics for component and number densities
integer, parameter :: stashcode_so4_nuc_sol   = 38485
integer, parameter :: stashcode_so4_ait_sol   = 38486
integer, parameter :: stashcode_so4_acc_sol   = 38487
integer, parameter :: stashcode_so4_cor_sol   = 38488
integer, parameter :: stashcode_bc_ait_sol    = 38489
integer, parameter :: stashcode_bc_acc_sol    = 38490
integer, parameter :: stashcode_bc_cor_sol    = 38491
integer, parameter :: stashcode_bc_ait_insol  = 38492
integer, parameter :: stashcode_oc_nuc_sol    = 38493
integer, parameter :: stashcode_oc_ait_sol    = 38494
integer, parameter :: stashcode_oc_acc_sol    = 38495
integer, parameter :: stashcode_oc_cor_sol    = 38496
integer, parameter :: stashcode_oc_ait_insol  = 38497
integer, parameter :: stashcode_ss_acc_sol    = 38498
integer, parameter :: stashcode_ss_cor_sol    = 38499
integer, parameter :: stashcode_du_acc_sol    = 38500
integer, parameter :: stashcode_du_cor_sol    = 38501
integer, parameter :: stashcode_du_acc_insol  = 38502
integer, parameter :: stashcode_du_cor_insol  = 38503
integer, parameter :: stashcode_n_nuc_sol     = 38504
integer, parameter :: stashcode_n_ait_sol     = 38505
integer, parameter :: stashcode_n_acc_sol     = 38506
integer, parameter :: stashcode_n_cor_sol     = 38507
integer, parameter :: stashcode_n_ait_insol   = 38508
integer, parameter :: stashcode_n_acc_insol   = 38509
integer, parameter :: stashcode_n_cor_insol   = 38510
integer, parameter :: stashcode_h2o_nuc_sol   = 38511
integer, parameter :: stashcode_h2o_ait_sol   = 38512
integer, parameter :: stashcode_h2o_acc_sol   = 38513
integer, parameter :: stashcode_h2o_cor_sol   = 38514
integer, parameter :: stashcode_h2o_total     = 38515
integer, parameter :: stashcode_so4_nuc_sol_load  = 38516
integer, parameter :: stashcode_so4_ait_sol_load  = 38517
integer, parameter :: stashcode_so4_acc_sol_load  = 38518
integer, parameter :: stashcode_so4_cor_sol_load  = 38519
integer, parameter :: stashcode_so4_total_load    = 38520
integer, parameter :: stashcode_bc_ait_sol_load   = 38521
integer, parameter :: stashcode_bc_acc_sol_load   = 38522
integer, parameter :: stashcode_bc_cor_sol_load   = 38523
integer, parameter :: stashcode_bc_ait_insol_load = 38524
integer, parameter :: stashcode_bc_total_load     = 38525
integer, parameter :: stashcode_oc_nuc_sol_load   = 38526
integer, parameter :: stashcode_oc_ait_sol_load   = 38527
integer, parameter :: stashcode_oc_acc_sol_load   = 38528
integer, parameter :: stashcode_oc_cor_sol_load   = 38529
integer, parameter :: stashcode_oc_ait_insol_load = 38530
integer, parameter :: stashcode_oc_total_load     = 38531
integer, parameter :: stashcode_du_acc_sol_load   = 38532
integer, parameter :: stashcode_du_cor_sol_load   = 38533
integer, parameter :: stashcode_du_acc_insol_load = 38534
integer, parameter :: stashcode_du_cor_insol_load = 38535
integer, parameter :: stashcode_du_total_load     = 38536
integer, parameter :: stashcode_ss_acc_sol_load   = 38537
integer, parameter :: stashcode_ss_cor_sol_load   = 38538
integer, parameter :: stashcode_ss_total_load     = 38539
integer, parameter :: stashcode_h2o_nuc_sol_load  = 38540
integer, parameter :: stashcode_h2o_ait_sol_load  = 38541
integer, parameter :: stashcode_h2o_acc_sol_load  = 38542
integer, parameter :: stashcode_h2o_cor_sol_load  = 38543
integer, parameter :: stashcode_h2o_total_load    = 38544
integer, parameter :: stashcode_h2o_mmr           = 38545
integer, parameter :: stashcode_so4_nuc_sol_ps    = 38900
integer, parameter :: stashcode_oc_nuc_sol_ps     = 38901
integer, parameter :: stashcode_so_nuc_sol_ps     = 38902
integer, parameter :: stashcode_so4_ait_sol_ps    = 38905
integer, parameter :: stashcode_bc_ait_sol_ps     = 38906
integer, parameter :: stashcode_oc_ait_sol_ps     = 38907
integer, parameter :: stashcode_so_ait_sol_ps     = 38908
integer, parameter :: stashcode_no3_ait_sol_ps    = 38909
integer, parameter :: stashcode_nh4_ait_sol_ps    = 38910
integer, parameter :: stashcode_so4_acc_sol_ps    = 38911
integer, parameter :: stashcode_bc_acc_sol_ps     = 38912
integer, parameter :: stashcode_oc_acc_sol_ps     = 38913
integer, parameter :: stashcode_ss_acc_sol_ps     = 38914
integer, parameter :: stashcode_du_acc_sol_ps     = 38916
integer, parameter :: stashcode_so_acc_sol_ps     = 38917
integer, parameter :: stashcode_no3_acc_sol_ps    = 38918
integer, parameter :: stashcode_nh4_acc_sol_ps    = 38919
integer, parameter :: stashcode_so4_cor_sol_ps    = 38920
integer, parameter :: stashcode_bc_cor_sol_ps     = 38921
integer, parameter :: stashcode_oc_cor_sol_ps     = 38922
integer, parameter :: stashcode_ss_cor_sol_ps     = 38923
integer, parameter :: stashcode_du_cor_sol_ps     = 38925
integer, parameter :: stashcode_so_cor_sol_ps     = 38926
integer, parameter :: stashcode_no3_cor_sol_ps    = 38927
integer, parameter :: stashcode_nh4_cor_sol_ps    = 38928
integer, parameter :: stashcode_bc_ait_insol_ps   = 38929
integer, parameter :: stashcode_oc_ait_insol_ps   = 38930
integer, parameter :: stashcode_du_acc_insol_ps   = 38931
integer, parameter :: stashcode_du_cor_insol_ps   = 38932
integer, parameter :: stashcode_nn_acc_sol_ps     = 38933
integer, parameter :: stashcode_nn_cor_sol_ps     = 38934
!UKCA dust super-coarse insoluble diagnostics
integer, parameter :: stashcode_du_sup_insol      = 38696
integer, parameter :: stashcode_n_sup_insol       = 38697
integer, parameter :: stashcode_du_sup_insol_load = 38698
integer, parameter :: stashcode_du_sup_insol_ps   = 38935

!PM10 and PM2.5 diagnostics
integer, parameter :: stashcode_pm10_wet  = 38560
integer, parameter :: stashcode_pm2p5_wet = 38561
integer, parameter :: stashcode_pm10_dry  = 38562
integer, parameter :: stashcode_pm2p5_dry = 38563
integer, parameter :: stashcode_pm10_so4  = 38564
integer, parameter :: stashcode_pm2p5_so4 = 38565
integer, parameter :: stashcode_pm10_bc   = 38566
integer, parameter :: stashcode_pm2p5_bc  = 38567
integer, parameter :: stashcode_pm10_oc   = 38568
integer, parameter :: stashcode_pm2p5_oc  = 38569
integer, parameter :: stashcode_pm10_ss   = 38570
integer, parameter :: stashcode_pm2p5_ss  = 38571
integer, parameter :: stashcode_pm10_du   = 38572
integer, parameter :: stashcode_pm2p5_du  = 38573
!First and last PM diagnostics
integer, parameter :: stashcode_pm_first  = 38560
integer, parameter :: stashcode_pm_last   = 38573

!UKCA Nitrate diagnostics
integer, parameter :: stashcode_nh4_nuc_sol       = 38646
integer, parameter :: stashcode_nh4_ait_sol       = 38647
integer, parameter :: stashcode_nh4_acc_sol       = 38648
integer, parameter :: stashcode_nh4_cor_sol       = 38649
integer, parameter :: stashcode_no3_nuc_sol       = 38650
integer, parameter :: stashcode_no3_ait_sol       = 38651
integer, parameter :: stashcode_no3_acc_sol       = 38652
integer, parameter :: stashcode_no3_cor_sol       = 38653
integer, parameter :: stashcode_nn_acc_sol        = 38654
integer, parameter :: stashcode_nn_cor_sol        = 38655
integer, parameter :: stashcode_nh4_nuc_sol_load  = 38656
integer, parameter :: stashcode_nh4_ait_sol_load  = 38657
integer, parameter :: stashcode_nh4_acc_sol_load  = 38658
integer, parameter :: stashcode_nh4_cor_sol_load  = 38659
integer, parameter :: stashcode_nh4_total_load    = 38660
integer, parameter :: stashcode_no3_nuc_sol_load  = 38661
integer, parameter :: stashcode_no3_ait_sol_load  = 38662
integer, parameter :: stashcode_no3_acc_sol_load  = 38663
integer, parameter :: stashcode_no3_cor_sol_load  = 38664
integer, parameter :: stashcode_no3_total_load    = 38665
integer, parameter :: stashcode_nn_acc_sol_load   = 38666
integer, parameter :: stashcode_nn_cor_sol_load   = 38667
integer, parameter :: stashcode_nn_total_load     = 38668
integer, parameter :: stashcode_pm10_nh4          = 38669
integer, parameter :: stashcode_pm2p5_nh4         = 38670
integer, parameter :: stashcode_pm10_no3          = 38671
integer, parameter :: stashcode_pm2p5_no3         = 38672
integer, parameter :: stashcode_pm10_nn           = 38673
integer, parameter :: stashcode_pm2p5_nn          = 38674

!UKCA Microplastic diagnostics
integer, parameter :: stashcode_mp_ait_sol        = 38747
integer, parameter :: stashcode_mp_acc_sol        = 38748
integer, parameter :: stashcode_mp_cor_sol        = 38749
integer, parameter :: stashcode_mp_ait_insol      = 38750
integer, parameter :: stashcode_mp_acc_insol      = 38751
integer, parameter :: stashcode_mp_cor_insol      = 38752
integer, parameter :: stashcode_mp_sup_insol      = 38753
integer, parameter :: stashcode_mp_ait_sol_load   = 38754
integer, parameter :: stashcode_mp_acc_sol_load   = 38755
integer, parameter :: stashcode_mp_cor_sol_load   = 38756
integer, parameter :: stashcode_mp_ait_insol_load = 38757
integer, parameter :: stashcode_mp_acc_insol_load = 38758
integer, parameter :: stashcode_mp_cor_insol_load = 38759
integer, parameter :: stashcode_mp_sup_insol_load = 38760
integer, parameter :: stashcode_mp_total_load     = 38761
integer, parameter :: stashcode_pm10_mp           = 38762
integer, parameter :: stashcode_pm2p5_mp          = 38763
integer, parameter :: stashcode_mp_ait_sol_ps     = 38936
integer, parameter :: stashcode_mp_acc_sol_ps     = 38937
integer, parameter :: stashcode_mp_cor_sol_ps     = 38938
integer, parameter :: stashcode_mp_ait_insol_ps   = 38939
integer, parameter :: stashcode_mp_acc_insol_ps   = 38940
integer, parameter :: stashcode_mp_cor_insol_ps   = 38941
integer, parameter :: stashcode_mp_sup_insol_ps   = 38942

!----------------------------------------------------------
! UKCA chemical diagnostics - section 50
!----------------------------------------------------------

integer, parameter :: stashcode_ukca_rxnflux_oh_ch4_trop = 50041
integer, parameter :: stashcode_ukca_nat                 = 50218
integer, parameter :: stashcode_ukca_o3_column_du        = 50219
integer, parameter :: stashcode_ukca_trop_ch4            = 50220
integer, parameter :: stashcode_ukca_trop_o3             = 50221
integer, parameter :: stashcode_ukca_trop_oh             = 50222
integer, parameter :: stashcode_ukca_strat_ch4           = 50223
integer, parameter :: stashcode_ukca_p_tropopause        = 50224
integer, parameter :: stashcode_ukca_strt_ch4_lss        = 50226
integer, parameter :: stashcode_ukca_jo1d                = 50228
integer, parameter :: stashcode_ukca_jno2                = 50229
integer, parameter :: stashcode_ukca_atmos_ch4           = 50231
integer, parameter :: stashcode_ukca_atmos_co            = 50232
integer, parameter :: stashcode_ukca_atmos_n2o           = 50233
integer, parameter :: stashcode_ukca_atmos_cfc12         = 50234
integer, parameter :: stashcode_ukca_atmos_cfc11         = 50235
integer, parameter :: stashcode_ukca_atmos_ch3br         = 50236
integer, parameter :: stashcode_ukca_atmos_h2            = 50237
integer, parameter :: stashcode_ukca_h2o_incr            = 50240
integer, parameter :: stashcode_ukca_jo2                 = 50245
integer, parameter :: stashcode_ukca_jo3p                = 50246
integer, parameter :: stashcode_ukca_so4_sad             = 50256
integer, parameter :: stashcode_ukca_h_plus              = 50442
integer, parameter :: stashcode_ukca_plumeria_height     = 50450

!----------------------------------------------------------
! GLOMAP_CLIM stashcodes - section 54
!----------------------------------------------------------

! Number density and component mass mixing ratios
integer, parameter :: stashcode_gc_nd_nuc_sol     = 54101
integer, parameter :: stashcode_gc_nuc_sol_su     = 54102
integer, parameter :: stashcode_gc_nd_ait_sol     = 54103
integer, parameter :: stashcode_gc_ait_sol_su     = 54104
integer, parameter :: stashcode_gc_ait_sol_bc     = 54105
integer, parameter :: stashcode_gc_ait_sol_oc     = 54106
integer, parameter :: stashcode_gc_nd_acc_sol     = 54107
integer, parameter :: stashcode_gc_acc_sol_su     = 54108
integer, parameter :: stashcode_gc_acc_sol_bc     = 54109
integer, parameter :: stashcode_gc_acc_sol_oc     = 54110
integer, parameter :: stashcode_gc_acc_sol_ss     = 54111
integer, parameter :: stashcode_gc_acc_sol_du     = 54112
integer, parameter :: stashcode_gc_nd_cor_sol     = 54113
integer, parameter :: stashcode_gc_cor_sol_su     = 54114
integer, parameter :: stashcode_gc_cor_sol_bc     = 54115
integer, parameter :: stashcode_gc_cor_sol_oc     = 54116
integer, parameter :: stashcode_gc_cor_sol_ss     = 54117
integer, parameter :: stashcode_gc_cor_sol_du     = 54118
integer, parameter :: stashcode_gc_nd_ait_ins     = 54119
integer, parameter :: stashcode_gc_ait_ins_bc     = 54120
integer, parameter :: stashcode_gc_ait_ins_oc     = 54121
integer, parameter :: stashcode_gc_nd_acc_ins     = 54122
integer, parameter :: stashcode_gc_acc_ins_du     = 54123
integer, parameter :: stashcode_gc_nd_cor_ins     = 54124
integer, parameter :: stashcode_gc_cor_ins_du     = 54125
integer, parameter :: stashcode_gc_nuc_sol_oc     = 54126
integer, parameter :: stashcode_gc_ait_sol_ss     = 54127
integer, parameter :: stashcode_gc_nuc_sol_so     = 54128
integer, parameter :: stashcode_gc_ait_sol_so     = 54129
integer, parameter :: stashcode_gc_acc_sol_so     = 54130
integer, parameter :: stashcode_gc_cor_sol_so     = 54131
integer, parameter :: stashcode_gc_acc_sol_nh4    = 54134
integer, parameter :: stashcode_gc_cor_sol_nh4    = 54135
integer, parameter :: stashcode_gc_acc_sol_no3    = 54138
integer, parameter :: stashcode_gc_cor_sol_no3    = 54139

! Fields required by RADAER
integer, parameter :: stashcode_gc_dryd_ait_sol      = 54921
integer, parameter :: stashcode_gc_dryd_acc_sol      = 54922
integer, parameter :: stashcode_gc_dryd_cor_sol      = 54923
integer, parameter :: stashcode_gc_dryd_ait_ins      = 54924
integer, parameter :: stashcode_gc_dryd_acc_ins      = 54925
integer, parameter :: stashcode_gc_dryd_cor_ins      = 54926
integer, parameter :: stashcode_gc_wetd_ait_sol      = 54927
integer, parameter :: stashcode_gc_wetd_acc_sol      = 54928
integer, parameter :: stashcode_gc_wetd_cor_sol      = 54929
integer, parameter :: stashcode_gc_rho_ait_sol       = 54930
integer, parameter :: stashcode_gc_rho_acc_sol       = 54931
integer, parameter :: stashcode_gc_rho_cor_sol       = 54932
integer, parameter :: stashcode_gc_rho_ait_ins       = 54933
integer, parameter :: stashcode_gc_rho_acc_ins       = 54934
integer, parameter :: stashcode_gc_rho_cor_ins       = 54935
integer, parameter :: stashcode_gc_pvol_ait_su_sol   = 54936
integer, parameter :: stashcode_gc_pvol_ait_bc_sol   = 54937
integer, parameter :: stashcode_gc_pvol_ait_oc_sol   = 54938
integer, parameter :: stashcode_gc_pvol_ait_no_sol   = 54939
integer, parameter :: stashcode_gc_pvol_ait_so_sol   = 54940
integer, parameter :: stashcode_gc_pvol_ait_h2o_sol  = 54941
integer, parameter :: stashcode_gc_pvol_acc_su_sol   = 54942
integer, parameter :: stashcode_gc_pvol_acc_bc_sol   = 54943
integer, parameter :: stashcode_gc_pvol_acc_oc_sol   = 54944
integer, parameter :: stashcode_gc_pvol_acc_ss_sol   = 54945
integer, parameter :: stashcode_gc_pvol_acc_no3_sol  = 54946
integer, parameter :: stashcode_gc_pvol_acc_du_sol   = 54947
integer, parameter :: stashcode_gc_pvol_acc_so_sol   = 54948
integer, parameter :: stashcode_gc_pvol_acc_h2o_sol  = 54951
integer, parameter :: stashcode_gc_pvol_cor_su_sol   = 54952
integer, parameter :: stashcode_gc_pvol_cor_bc_sol   = 54953
integer, parameter :: stashcode_gc_pvol_cor_oc_sol   = 54954
integer, parameter :: stashcode_gc_pvol_cor_ss_sol   = 54955
integer, parameter :: stashcode_gc_pvol_cor_no3_sol  = 54956
integer, parameter :: stashcode_gc_pvol_cor_du_sol   = 54957
integer, parameter :: stashcode_gc_pvol_cor_so_sol   = 54958
integer, parameter :: stashcode_gc_pvol_cor_h2o_sol  = 54961
integer, parameter :: stashcode_gc_pvol_ait_bc_ins   = 54962
integer, parameter :: stashcode_gc_pvol_ait_oc_ins   = 54963
integer, parameter :: stashcode_gc_pvol_acc_du_ins   = 54964
integer, parameter :: stashcode_gc_pvol_cor_du_ins   = 54965

! Fields required by ACTIVATE

integer, parameter :: stashcode_gc_cf_liquid         = 54476
integer, parameter :: stashcode_gc_cdncwt            = 54477
integer, parameter :: stashcode_gc_cdnc              = 54968

!----------------------------------------------------------
! ENDGame stashcodes
!----------------------------------------------------------

integer, parameter :: stashcode_etadot            = 387
integer, parameter :: stashcode_thetavd           = 388
integer, parameter :: stashcode_dry_rho           = 389
integer, parameter :: stashcode_exner_surf        = 398
integer, parameter :: stashcode_psiw_surf         = 390
integer, parameter :: stashcode_psiw_lid          = 397
integer, parameter :: stashcode_mv                = 391
integer, parameter :: stashcode_mcl               = 392
integer, parameter :: stashcode_mcf               = 393
integer, parameter :: stashcode_mr                = 394
integer, parameter :: stashcode_mgr               = 395
integer, parameter :: stashcode_mcf2              = 396

! ---------------------------------------------------------
! PWS stashcodes - section 20
! ---------------------------------------------------------
integer, parameter :: stashcode_pws_thickness500   =  1
integer, parameter :: stashcode_pws_thickness850   =  2
integer, parameter :: stashcode_pws_windspeed10m   =  3
integer, parameter :: stashcode_pws_windspeedplev  =  4
integer, parameter :: stashcode_pws_divergence     =  5
integer, parameter :: stashcode_pws_rel_vorticity  =  6
integer, parameter :: stashcode_pws_mtn_wave_turb  =  7
integer, parameter :: stashcode_pws_conv_cld_dep   =  12
integer, parameter :: stashcode_pws_precip_sym     =  14
integer, parameter :: stashcode_pws_cat_turb       =  16
integer, parameter :: stashcode_pws_max_cat        =  17
integer, parameter :: stashcode_pws_max_cat_press  =  18
integer, parameter :: stashcode_pws_max_wind_ub    =  20
integer, parameter :: stashcode_pws_max_wind_vb    =  21
integer, parameter :: stashcode_pws_max_wind_pb    =  22
integer, parameter :: stashcode_pws_max_wind_icao  =  23
integer, parameter :: stashcode_pws_snow_prob      =  28
integer, parameter :: stashcode_pws_contrail_bot   =  29
integer, parameter :: stashcode_pws_contrail_top   =  30
integer, parameter :: stashcode_pws_freezing_ht    =  33
integer, parameter :: stashcode_pws_freezing_press =  34
integer, parameter :: stashcode_pws_freezing_icao  =  35
integer, parameter :: stashcode_pws_isotherm_ms20_ht    =  36
integer, parameter :: stashcode_pws_isotherm_ms20_press =  37
integer, parameter :: stashcode_pws_isotherm_ms20_icao  =  38
integer, parameter :: stashcode_pws_conv_icao_base =  39
integer, parameter :: stashcode_pws_conv_icao_top  =  40
integer, parameter :: stashcode_pws_max_wind_base  =  41
integer, parameter :: stashcode_pws_max_wind_top   =  42
integer, parameter :: stashcode_pws_icing_pot_diag =  43
integer, parameter :: stashcode_pws_cloudturb_pot_diag = 45
integer, parameter :: stashcode_pws_wafc_caturb    =  47
integer, parameter :: stashcode_pws_dustconc_surf  =  58
integer, parameter :: stashcode_pws_dustconc_5000  =  59
integer, parameter :: stashcode_pws_zenithdelay    =  60
integer, parameter :: stashcode_pws_isotherm_ms70_ht    =  61
integer, parameter :: stashcode_pws_isotherm_ms70_press =  62
integer, parameter :: stashcode_pws_isotherm_ms70_icao  =  63
integer, parameter :: stashcode_pws_1p5m_vis_tot        =  70
integer, parameter :: stashcode_pws_1p5m_vis_dust       =  71
integer, parameter :: stashcode_pws_inv_richardson      =  76
integer, parameter :: stashcode_pws_1p5m_vis_land       =  78
integer, parameter :: stashcode_pws_1p5m_vis_ssi        =  79
integer, parameter :: stashcode_pws_upd_helicity_5k     =  80
integer, parameter :: stashcode_pws_tropopause_press    =  84
integer, parameter :: stashcode_pws_tropopause_temp     =  85
integer, parameter :: stashcode_pws_tropopause_ht       =  86
integer, parameter :: stashcode_pws_tropopause_icao     =  87

integer, parameter :: stashcode_pws_mlcape     =  110
integer, parameter :: stashcode_pws_mlcin      =  111
integer, parameter :: stashcode_pws_mucape     =  112
integer, parameter :: stashcode_pws_mucin      =  113
integer, parameter :: stashcode_pws_sbcape     =  114
integer, parameter :: stashcode_pws_sbcin      =  115
integer, parameter :: stashcode_pws_mlelz      =  116
integer, parameter :: stashcode_pws_mulplz     =  117
integer, parameter :: stashcode_pws_muelz      =  118
integer, parameter :: stashcode_pws_eibasez    =  119
integer, parameter :: stashcode_pws_eitopz     =  120

! ---------------------------------------------------------
! Stashcodes which PWS diags depend on
! ---------------------------------------------------------
! Boundary layer, section 3
integer, parameter :: stashcode_bl_windspeed_10mb = 227
integer, parameter :: stashcode_bl_1p5m_temp      = 236
integer, parameter :: stashcode_bl_1p5m_temp_land = 547
integer, parameter :: stashcode_bl_1p5m_temp_ssi  = 545
integer, parameter :: stashcode_bl_1p5m_q         = 237
integer, parameter :: stashcode_bl_1p5m_vis_tot   = 281
integer, parameter :: stashcode_bl_tke            = 473
integer, parameter :: stashcode_bl_1p5m_vis_land  = 576
integer, parameter :: stashcode_bl_1p5m_vis_ssi   = 577
! Convection, section 5
integer, parameter :: stashcode_conv_icao_base    = 210
integer, parameter :: stashcode_conv_icao_top     = 211
! Gravity wave drag, section 6
integer, parameter :: stashcode_gwd_stress_lev_u  = 201
integer, parameter :: stashcode_gwd_stress_lev_v  = 202
! Processed dynamics, section 15
integer, parameter :: stashcode_dyn_wind_ub       = 201
integer, parameter :: stashcode_dyn_wind_vb       = 202
! Processed physics, section 16
integer, parameter :: stashcode_phy_geopht        = 202
! Aerosols
integer, parameter :: stashcode_aero_total_dust   = 257

end module um_stashcode_mod
