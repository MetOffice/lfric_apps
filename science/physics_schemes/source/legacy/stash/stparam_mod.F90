! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module stparam_mod

implicit none

!  Purpose: Meaningful parameter names for STASH processing routines.
!           Both a long name and short name have been declared, to
!           reduce the existence of "magic" numbers in STASH.
!           Format is that the address of the item is declare in
!           both long and short form. example is;
!           !Item number
!           integer,parameter:: st_item_code=3,s_item=3
!
!  Programming standard: UM Doc Paper 3
!
!  External documentation:
!    UMDP C04 - Storage handling and diagnostic system (STASH)

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH

! Internal model number address
integer,parameter:: st_model_code = 28
integer,parameter:: s_modl        = 28

! Section Number address
integer,parameter:: st_sect_no_code = 2
integer,parameter:: s_sect          = 2

! Item number address
integer,parameter:: st_item_code=1,s_item=1

! Processing Code address
integer,parameter:: st_proc_no_code=3,s_proc=3

! subsidiary codes for st_proc_no_code now
integer,parameter:: st_replace_code=1
integer,parameter:: st_accum_code=2
integer,parameter:: st_time_mean_code=3
integer,parameter:: st_time_series_code=4
integer,parameter:: st_max_code=5
integer,parameter:: st_min_code=6
integer,parameter:: st_append_traj_code=7
integer,parameter:: st_time_series_mean=8
integer,parameter:: st_variance_code=9

! Frequency (Input & output) address
integer,parameter:: st_freq_code=4,s_freq=4

! Offset for sampling
integer,parameter:: st_offset_code=30,s_offs=30

! start timestep address
integer,parameter:: st_start_time_code=5,s_times=5

! timestep zero output logical address
integer,parameter:: st_output_ts0_code=34

! end timestep address
integer,parameter:: st_end_time_code=6,s_timee=6

! period in timesteps address
integer,parameter:: st_period_code=7,s_period=7

! infinite end/period value
integer,parameter:: st_infinite_time=-1

integer,parameter:: st_end_of_list=-1 !end-of-list marker in times

! grid point stuff
! gridpoint info address
integer,parameter:: st_gridpoint_code=8,s_grid=8

! now subsid grid point stuff
! no masking done
integer,parameter:: stash_null_mask_code=1,s_nomask=1

! land mask conds
integer,parameter:: stash_land_mask_code=2,s_lndms=2

! sea mask code
integer,parameter:: stash_sea_mask_code=3,s_seams =3

! minimal mdi masking done
integer,parameter:: stash_nmdi_mask_code=4,s_nmdims=4

! processing options

! size of block for gridpoint code
integer,parameter:: block_size=10

! max code for vertical mean subroutine
integer,parameter:: extract_base=block_size*0

! base codes for vertical mean subroutine
integer,parameter:: extract_top=block_size*1

! max code for vertical mean subroutine
integer,parameter:: vert_mean_base=block_size*1

! base codes for vertical mean subroutine
integer,parameter:: vert_mean_top=block_size*2

! max code for zonal mean subroutine
integer,parameter:: zonal_mean_base=block_size*2

! base codes for zonal mean subroutine
integer,parameter:: zonal_mean_top=block_size*3

! max code for meridional mean subroutine
integer,parameter:: merid_mean_base=block_size*3

! base codes for meridional mean subroutine
integer,parameter:: merid_mean_top=block_size*4

! max code for field mean subroutine
integer,parameter:: field_mean_base=block_size*4

! base codes for field mean subroutine
integer,parameter:: field_mean_top=block_size*5

! max code for global mean subroutine
integer,parameter:: global_mean_base=block_size*5

! base codes for global mean subroutine
integer,parameter:: global_mean_top=block_size*6

! Weighting

! weighting info address
integer,parameter:: st_weight_code=9,s_weight=9

integer,parameter:: stash_weight_null_code  =0,s_noweight  =0
integer,parameter:: stash_weight_area_code  =1,s_areaweight=1
integer,parameter:: stash_weight_volume_code=2,s_volweight =2
integer,parameter:: stash_weight_mass_code  =3,s_massweight=3

! Domain definition

! row addresses
integer,parameter:: st_north_code=12,s_north=12
integer,parameter:: st_south_code=13,s_south=13
integer,parameter:: st_west_code =14,s_west =14
integer,parameter:: st_east_code =15,s_east =15

! Levels

! input bottom level address
integer,parameter:: st_input_bottom=10,s_bottom =10

! special code
integer,parameter:: st_special_code=100,s_special=100

! input top level address
integer,parameter:: st_input_top=11,s_top=11

! output bottom level address
integer,parameter:: st_output_bottom=21,s_outbot=21

! output top level address
integer,parameter:: st_output_top=22,s_outtop=22

! input code addres
integer,parameter:: st_input_code=16,s_input=16

! input length of diagnostic address
integer,parameter:: st_input_length=17,s_length=17

! output code address
integer,parameter:: st_output_code=18,s_output=18

! Pointer to D1 addressing information
! Pos of item in D1 for relevant submodel
integer,parameter:: st_position_in_d1=29,st_d1pos=29

! Output destination options

integer,parameter:: st_dump=1
integer,parameter:: st_secondary=2

! output length of diagnostic address
integer,parameter:: st_output_length=19,s_outlen=19
! output length on dump
integer,parameter:: st_dump_output_length=32,s_doutlen=32
! output length of a single level on dump
integer,parameter:: st_dump_level_output_length=33,s_dlevoutlen=33
! start locn of diag after stash output address
integer,parameter:: st_output_addr=20,s_outadd=20
! output address on dump
integer,parameter:: st_dump_output_addr=31,s_doutadd=31

! ptr to dump lookup header address
integer,parameter:: st_lookup_ptr=23

! ptr into stash_series where control data address
integer,parameter:: st_series_ptr=24

! subsid stuff for time series
integer,parameter:: series_grid_type=1
integer,parameter:: series_grid_code=0
integer,parameter:: series_long_code=1
integer,parameter:: series_size=2
integer,parameter:: series_proc_code=3
integer,parameter:: series_north=4
integer,parameter:: series_south=5
integer,parameter:: series_west=6
integer,parameter:: series_east=7
integer,parameter:: series_list_start=8
integer,parameter:: series_list_end=9
integer,parameter:: record_size=9

! Miscellaneous parameters

! system/user tag field in stlist address
integer,parameter:: st_macrotag=25

! Pseudo-level list pointers

! pseudo-levels input list address
integer,parameter:: st_pseudo_in=26

! pseudo-levels output list address
integer,parameter:: st_pseudo_out=27

! Internal horizontal gridtype codes common to all diagnostics

integer,parameter:: st_tp_grid =1 ! T-p grid
integer,parameter:: st_uv_grid =2 ! u-v grid
integer,parameter:: st_cu_grid =3 ! C-grid u point
integer,parameter:: st_cv_grid =4 ! C-grid v point
integer,parameter:: st_zt_grid =5 ! Zonal T-grid
integer,parameter:: st_zu_grid =6 ! Zonal u-grid
integer,parameter:: st_mt_grid =7 ! Meridional T-grid
integer,parameter:: st_mu_grid =8 ! Meridional u-grid
integer,parameter:: st_riv_grid= 23    ! river_routing grid
integer,parameter:: st_scalar  =9 ! Scalar (ie. single value)

! Pointer to original diagnostic address in STASH
integer,parameter:: st_diag_address=35

! Pointer to required domain area
integer,parameter:: st_domain_code=36

! Domain area possible values (IOPA, see UMDP-C04)
integer,parameter :: st_domain_global        = 1
integer,parameter :: st_domain_n_hemisphere  = 2
integer,parameter :: st_domain_s_hemisphere  = 3
integer,parameter :: st_domain_30_to_90_N    = 4
integer,parameter :: st_domain_30_to_90_S    = 5
integer,parameter :: st_domain_0_to_30_N     = 6
integer,parameter :: st_domain_0_to_30_S     = 7
integer,parameter :: st_domain_30_S_to_30_N  = 8
integer,parameter :: st_domain_whole_degrees = 9
integer,parameter :: st_domain_gridpoints    = 10

! Domain level possible values (IOPL, see UMDP-C04)
integer,parameter :: st_levels_model_rho        = 1
integer,parameter :: st_levels_model_theta      = 2
integer,parameter :: st_levels_pressure         = 3
integer,parameter :: st_levels_geometric_height = 4
integer,parameter :: st_levels_single           = 5
integer,parameter :: st_levels_deep_soil        = 6
integer,parameter :: st_levels_potential_temp   = 7
integer,parameter :: st_levels_potential_vort   = 8
integer,parameter :: st_levels_cloud_thresh     = 9

! Pointer to units for sampling period (frequency)
integer,parameter:: st_time_unit2_code=37

! Pointer to units for output times of diagnostic
integer,parameter:: st_time_unit3_code=38

! Pointer to output file type
integer,parameter:: st_output_type=39

! Output file types
integer,parameter:: st_dump_store=0
integer,parameter:: st_fieldsfile=1
integer,parameter:: st_netcdf=2

end module stparam_mod
