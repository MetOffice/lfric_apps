!----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Provides an implementation of the Psy layer for physics

!> @details Contains hand-rolled versions of the Psy layer that can be used for
!> simple testing and development of the scientific code

module psykal_lite_phys_mod

  use constants_mod,         only : i_def, r_def
  use field_mod,             only : field_type, field_proxy_type
  use integer_field_mod,     only : integer_field_type, integer_field_proxy_type
  use mesh_mod,              only : mesh_type
  
  implicit none
  public

contains
  !---------------------------------------------------------------------
  !> LFRic and PSyclone currently do not have a mechanism to loop over a subset
  !> of cells in a horizontal domain. Furthermore this is a psykal_lite file which is not PSycloned.
  !> PSyclone development to this aim can be tracked at https://github.com/stfc/PSyclone/issues/3193
  !> Future refactoring may seek to push these loops into a kernel layer and then to apply PSyclone.
  !>
  !> The orographic drag kernel only needs to be applied to a subset of land
  !> points where the standard deviation of subgrid orography is more than zero.
  !>
  !> invoke_orographic_drag_kernel: Invokes the kernel which calls the UM
  !> orographic drag scheme only on points where the standard deviation
  !> of subgrid orography is more than zero.
  subroutine invoke_orographic_drag_kernel(                              &
                      du_orog_blk, dv_orog_blk, du_orog_gwd, dv_orog_gwd,&
                      dtemp_orog_blk, dtemp_orog_gwd, u_in_w3, v_in_w3,&
                      wetrho_in_w3, theta_in_wth, exner_in_w3, sd_orog,  &
                      grad_xx_orog, grad_xy_orog, grad_yy_orog,          &
                      mr_v, mr_cl, mr_ci,                                &
                      height_w3, height_wth,                             &
                      taux_orog_blk, tauy_orog_blk,                      &
                      taux_orog_gwd, tauy_orog_gwd )

    use orographic_drag_kernel_mod, only: orographic_drag_kernel_code
    use mesh_mod, only: mesh_type
    use physics_config_mod,  only : gw_segment
    implicit none

    ! Increments from orographic drag
    type(field_type), intent(inout) :: du_orog_blk, dv_orog_blk, & ! Winds
                                       du_orog_gwd, dv_orog_gwd, &
                                       dtemp_orog_blk, dtemp_orog_gwd ! Temperature

    ! Inputs to orographic drag scheme
    type(field_type), intent(in) :: u_in_w3, v_in_w3,           & ! Winds
                                    wetrho_in_w3, theta_in_wth, & ! Density, Temperature
                                    exner_in_w3,                & ! Exner pressure
                                    sd_orog, grad_xx_orog,      & ! Orography ancils
                                    grad_xy_orog, grad_yy_orog, & !
                                    mr_v, mr_cl, mr_ci,         & ! mixing ratios
                                    height_w3, height_wth         ! Heights
    ! Diagnostics from orographic drag
    ! Stress from ...
    type(field_type), intent(inout) :: taux_orog_blk, tauy_orog_blk, & ! ... orographic flow blocking drag
                                       taux_orog_gwd, tauy_orog_gwd    ! ... orographic gravity wave drag

    integer :: cell

    ! Number of degrees of freedom
    integer :: ndf_w3, undf_w3, ndf_wtheta, undf_wtheta

    ! Integers for segmentation
    integer :: applicable_points, nlayers, loop_upper_bound, segment, seg_len, &
               l_bound, u_bound, n_segments , seg_target
    
    ! These are in ANY_DISCONTINUOUS_SPACE_1
    integer :: ndf_adspc1_sd_orog, undf_adspc1_sd_orog

    integer, allocatable ::    cell_index(:) !dhc record the points with orography

    
    type(field_proxy_type) :: du_blk_proxy, dv_blk_proxy,             &
                              du_orog_gwd_proxy, dv_orog_gwd_proxy,   &
                              dtemp_blk_proxy, dtemp_orog_gwd_proxy,  &
                              u_in_w3_proxy, v_in_w3_proxy,           &
                              wetrho_in_w3_proxy, theta_in_wth_proxy, &
                              exner_in_w3_proxy,                      &
                              sd_orog_proxy, grad_xx_orog_proxy,      &
                              grad_xy_orog_proxy, grad_yy_orog_proxy, &
                              mr_v_proxy, mr_cl_proxy, mr_ci_proxy,   &
                              height_w3_proxy, height_wth_proxy,      &
                              taux_blk_proxy, tauy_blk_proxy,         &
                              taux_gwd_proxy, tauy_gwd_proxy

    integer, pointer :: map_adspc1_sd_orog(:,:) => null(), &
                        map_w3(:,:) => null(),                  &
                        map_wtheta(:,:) => null()

    TYPE(mesh_type), pointer :: mesh => null()

    ! Initialise field and/or operator proxies
    du_blk_proxy = du_orog_blk%get_proxy()
    dv_blk_proxy = dv_orog_blk%get_proxy()
    du_orog_gwd_proxy = du_orog_gwd%get_proxy()
    dv_orog_gwd_proxy = dv_orog_gwd%get_proxy()
    dtemp_blk_proxy = dtemp_orog_blk%get_proxy()
    dtemp_orog_gwd_proxy = dtemp_orog_gwd%get_proxy()
    u_in_w3_proxy = u_in_w3%get_proxy()
    v_in_w3_proxy = v_in_w3%get_proxy()
    wetrho_in_w3_proxy = wetrho_in_w3%get_proxy()
    theta_in_wth_proxy = theta_in_wth%get_proxy()
    exner_in_w3_proxy = exner_in_w3%get_proxy()
    sd_orog_proxy = sd_orog%get_proxy()
    grad_xx_orog_proxy = grad_xx_orog%get_proxy()
    grad_xy_orog_proxy = grad_xy_orog%get_proxy()
    grad_yy_orog_proxy = grad_yy_orog%get_proxy()
    mr_v_proxy = mr_v%get_proxy()
    mr_cl_proxy = mr_cl%get_proxy()
    mr_ci_proxy = mr_ci%get_proxy()
    height_w3_proxy = height_w3%get_proxy()
    height_wth_proxy = height_wth%get_proxy()

    taux_blk_proxy = taux_orog_blk%get_proxy()
    tauy_blk_proxy = tauy_orog_blk%get_proxy()
    taux_gwd_proxy = taux_orog_gwd%get_proxy()
    tauy_gwd_proxy = tauy_orog_gwd%get_proxy()

    ! Initialise number of layers
    nlayers = du_blk_proxy%vspace%get_nlayers()

    ! Create a mesh object
    mesh => du_orog_blk%get_mesh()

    ! Look-up dofmaps for each function space
    map_w3 => du_blk_proxy%vspace%get_whole_dofmap()
    map_wtheta => dtemp_blk_proxy%vspace%get_whole_dofmap()
    map_adspc1_sd_orog => sd_orog_proxy%vspace%get_whole_dofmap()

    ! Initialise number of DoFs for w3
    ndf_w3 = du_blk_proxy%vspace%get_ndf()
    undf_w3 = du_blk_proxy%vspace%get_undf()

    ! Initialise number of DoFs for wtheta
    ndf_wtheta = dtemp_blk_proxy%vspace%get_ndf()
    undf_wtheta = dtemp_blk_proxy%vspace%get_undf()

    ! Initialise number of DoFs for adspc1_sd_orog
    ndf_adspc1_sd_orog = sd_orog_proxy%vspace%get_ndf()
    undf_adspc1_sd_orog = sd_orog_proxy%vspace%get_undf()

    ! Temporary variable for loop bound helps OpenMP
    loop_upper_bound = mesh%get_last_edge_cell()    
  
    applicable_points = 0   
    !cell index padded with 0s beyond applicable_points
    allocate(cell_index(loop_upper_bound))
    cell_index = 0
    do cell=1, loop_upper_bound
      ! Only call orographic_drag_kernel_code at points where the
      ! standard deviation of the subgrid orography is more than zero.
      if ( sd_orog_proxy%data(map_adspc1_sd_orog(1, cell)) > 0.0_r_def ) then

         ! We want to select the cells which have orographic drag
         applicable_points = applicable_points+1
         cell_index(applicable_points) = cell

      end if ! sd_orog_proxy%data(map_adspc1_sd_orog(1, cell)) > 0.0_r_def

    end do

    ! Check that the value of gw_segment is acceptable as a target length, otherwise use 1 for safety
    if (gw_segment .le. 0 .or. gw_segment .gt. applicable_points) then
      seg_target = 1
    else
      seg_target = gw_segment
    end if

    n_segments = (applicable_points + seg_target - 1) / seg_target

    ! Call orographic_drag_kernel_code if you have applicable points
    ! seg_len will be seg_target unless you have too few points at the end - l and u bound are of the segment iterated over
    if (applicable_points .gt. 0) then
      !$omp parallel do default(shared), private(segment, seg_len, l_bound, u_bound), schedule(dynamic)
      do segment=1, n_segments
        seg_len = min(seg_target, applicable_points - (segment - 1) * seg_target)
        l_bound = (segment - 1) * seg_target + 1
        u_bound = l_bound + seg_len - 1
        call orographic_drag_kernel_code(                                                &
                        loop_upper_bound, nlayers, du_blk_proxy%data, dv_blk_proxy%data, &
                        du_orog_gwd_proxy%data, dv_orog_gwd_proxy%data,                  &
                        dtemp_blk_proxy%data, dtemp_orog_gwd_proxy%data,                 &
                        u_in_w3_proxy%data, v_in_w3_proxy%data,                          &
                        wetrho_in_w3_proxy%data, theta_in_wth_proxy%data,                &
                        exner_in_w3_proxy%data,                                          &
                        sd_orog_proxy%data, grad_xx_orog_proxy%data,                     &
                        grad_xy_orog_proxy%data, grad_yy_orog_proxy%data,                &
                        mr_v_proxy%data, mr_cl_proxy%data, mr_ci_proxy%data,             &
                        height_w3_proxy%data, height_wth_proxy%data,                     &
                        taux_blk_proxy%data, tauy_blk_proxy%data,                        &
                        taux_gwd_proxy%data, tauy_gwd_proxy%data,                        &
                        ndf_w3, undf_w3, map_w3,                                         &
                        ndf_wtheta, undf_wtheta, map_wtheta,                             &
                        ndf_adspc1_sd_orog, undf_adspc1_sd_orog,                         &
                        map_adspc1_sd_orog, seg_len, cell_index(l_bound:u_bound))
      end do
      !$omp end parallel do
    end if
    deallocate(cell_index)

    ! Set halos dirty/clean for fields modified in the above loop
    call du_blk_proxy%set_dirty()
    call dv_blk_proxy%set_dirty()
    call du_orog_gwd_proxy%set_dirty()
    call dv_orog_gwd_proxy%set_dirty()
    call dtemp_blk_proxy%set_dirty()
    call dtemp_orog_gwd_proxy%set_dirty()
    call taux_blk_proxy%set_dirty()
    call tauy_blk_proxy%set_dirty()
    call taux_gwd_proxy%set_dirty()
    call tauy_gwd_proxy%set_dirty()

  end subroutine invoke_orographic_drag_kernel
  !---------------------------------------------------------------------
  !> Contains the PSy-layer to build the stochastic physics
  !> forcing pattern. At the moment it requires to pass the arrays
  !> my_coeff_rad and my_phi_stph via the kernel argument.
  !> PSyclone does not recognize arrays in the argument yet
  !> this functionality is being developed in PSyclone ticket 1312
  !> at https://github.com/stfc/PSyclone/issues/1312
  !> Hence this module could be removed once the PSyclone ticket is
  !> completed
  subroutine invoke_spectral_2_cs_kernel_type(fp, longitude, pnm_star,         &
                                              coeffc_phase, coeffs_phase,      &
                                              stph_level_bottom,               &
                                              stph_level_top,                  &
                                              stph_n_min, stph_n_max)

  use spectral_2_cs_kernel_mod, ONLY: spectral_2_cs_code
  use mesh_mod, ONLY: mesh_type

  implicit none

  integer(KIND=i_def), intent(in) :: stph_level_bottom, stph_level_top,stph_n_min, stph_n_max
  type(field_type), intent(in) :: fp, longitude, pnm_star
  integer(KIND=i_def) cell
  integer(KIND=i_def) nlayers
  type(field_proxy_type) fp_proxy, longitude_proxy, pnm_star_proxy
  integer(KIND=i_def), pointer :: map_adspc1_longitude(:,:) => null(), map_adspc2_pnm_star(:,:) => null(), &
  &map_wspace(:,:) => null()
  integer(KIND=i_def) ndf_wspace, undf_wspace, ndf_adspc1_longitude, undf_adspc1_longitude, ndf_adspc2_pnm_star, &
  &undf_adspc2_pnm_star
  type(mesh_type), pointer :: mesh => null()

  ! Add arrays my_coeff_rad & my_phi_stph by hand
  real(kind=r_def), intent(in), dimension(:,:) :: coeffc_phase
  real(kind=r_def), intent(in), dimension(:,:) :: coeffs_phase
  integer(kind=i_def), parameter :: nranks_array = 2
  integer(kind=i_def), dimension(nranks_array) :: dims_array

  ! Get the upper bound for each rank of each scalar array
  ! Do dims_array for coeffc and coeffs?
  dims_array = shape(coeffc_phase)

  !
  ! Initialise field and/or operator proxies
  !
  fp_proxy     = fp%get_proxy()
  longitude_proxy  = longitude%get_proxy()
  pnm_star_proxy   = pnm_star%get_proxy()
  !
  ! Initialise number of layers
  !
  nlayers = fp_proxy%vspace%get_nlayers()
  !
  ! Create a mesh object
  !
  mesh => fp_proxy%vspace%get_mesh()
  !
  ! Look-up dofmaps for each function space
  !
  map_wspace           => fp_proxy%vspace%get_whole_dofmap()
  map_adspc1_longitude => longitude_proxy%vspace%get_whole_dofmap()
  map_adspc2_pnm_star  => pnm_star_proxy%vspace%get_whole_dofmap()
  !
  ! Initialise number of DoFs for wtheta
  !
  ndf_wspace  = fp_proxy%vspace%get_ndf()
  undf_wspace = fp_proxy%vspace%get_undf()
  !
  ! Initialise number of DoFs for adspc1_longitude
  !
  ndf_adspc1_longitude  = longitude_proxy%vspace%get_ndf()
  undf_adspc1_longitude = longitude_proxy%vspace%get_undf()
  !
  ! Initialise number of DoFs for adspc2_pnm_star
  !
  ndf_adspc2_pnm_star  = pnm_star_proxy%vspace%get_ndf()
  undf_adspc2_pnm_star = pnm_star_proxy%vspace%get_undf()
  !
  ! Call kernels and communication routines
  !
  do cell=1,mesh%get_last_edge_cell()
      call spectral_2_cs_code(nlayers, &
                              ! Add fields
                              fp_proxy%data, longitude_proxy%data, pnm_star_proxy%data, &
                              ! Add arrays
                              nranks_array, dims_array, coeffc_phase, coeffs_phase, &
                              ! Add SPT scalars
                              stph_level_bottom, stph_level_top, stph_n_min, stph_n_max, &
                              ! Add fields' assoc. space variables
                              ndf_wspace, undf_wspace, map_wspace(:,cell), &
                              ndf_adspc1_longitude, undf_adspc1_longitude, map_adspc1_longitude(:,cell), &
                              ndf_adspc2_pnm_star, undf_adspc2_pnm_star, map_adspc2_pnm_star(:,cell))
  end do
  !
  ! Set halos dirty/clean for fields modified in the above loop
  !
  call fp_proxy%set_dirty()
  !
  !
  end subroutine invoke_spectral_2_cs_kernel_type
  !---------------------------------------------------------------------
  ! PSyclone currently doesn't work for DOMAIN kernels with stencil fields
  ! See PSyclone #1948
    SUBROUTINE invoke_jules_exp_kernel_type(ncells, ncells_halo, theta, exner_in_wth, u_in_w3, v_in_w3, mr_n, mr_n_1, mr_n_2, height_w3, height_wth, &
&zh, z0msea, z0m, tile_fraction, leaf_area_index, canopy_height, peak_to_trough_orog, silhouette_area_orog, soil_albedo, &
&soil_roughness, soil_moist_wilt, soil_moist_crit, soil_moist_sat, soil_thermal_cond, soil_suction_sat, clapp_horn_b, &
&soil_respiration, thermal_cond_wet_soil, sea_u_current_ptr, sea_v_current_ptr, sea_ice_temperature, sea_ice_conductivity, sea_ice_pensolar, &
&sea_ice_pensolar_frac_direct, sea_ice_pensolar_frac_diffuse, tile_temperature, tile_snow_mass, n_snow_layers, snow_depth, &
&snow_layer_thickness, snow_layer_ice_mass, snow_layer_liq_mass, snow_layer_temp, surface_conductance, canopy_water, &
&soil_temperature, soil_moisture, unfrozen_soil_moisture, frozen_soil_moisture, tile_heat_flux, tile_moisture_flux, net_prim_prod, &
&cos_zenith_angle, skyview, &
sw_up_tile, tile_lw_grey_albedo, &
sw_down_surf, lw_down_surf, sw_down_blue_surf, sw_direct_blue_surf, dd_mf_cb, ozone, cf_bulk, cf_liquid, rhokm_bl, &
&surf_interp, rhokh_bl, moist_flux_bl, heat_flux_bl, gradrinr, &
&alpha1_tile, ashtf_prime_tile, dtstar_tile, fracaero_t_tile, fracaero_s_tile, z0h_tile, &
&z0m_tile, rhokh_tile, chr1p5m_tile, resfs_tile, gc_tile, canhc_tile, tile_water_extract, blend_height_tq, z0m_eff, ustar, &
&soil_moist_avail, snow_unload_rate, albedo_obs_scaling, soil_clay, soil_sand, dust_mrel, dust_flux, day_of_year, second_of_day, &
flux_e, flux_h, urbwrr, urbhwr, urbhgt, urbztm, urbdisp, &
&rhostar, recip_l_mo_sea, &
&t1_sd_2d, q1_sd_2d, gross_prim_prod, z0h_eff, ocn_cpl_point, stencil_depth)
      USE jules_exp_kernel_mod, ONLY: jules_exp_code
      USE mesh_mod, ONLY: mesh_type
      USE stencil_dofmap_mod, ONLY: STENCIL_REGION
      USE stencil_dofmap_mod, ONLY: stencil_dofmap_type
      implicit none
      TYPE(field_type), intent(in) :: theta, exner_in_wth, u_in_w3, v_in_w3, mr_n, mr_n_1, mr_n_2, height_w3, height_wth, zh, &
&z0msea, z0m, tile_fraction, leaf_area_index, canopy_height, peak_to_trough_orog, silhouette_area_orog, soil_albedo, &
&soil_roughness, soil_moist_wilt, soil_moist_crit, soil_moist_sat, soil_thermal_cond, soil_suction_sat, clapp_horn_b, &
&soil_respiration, thermal_cond_wet_soil, sea_u_current_ptr, sea_v_current_ptr, sea_ice_temperature, sea_ice_conductivity, &
sea_ice_pensolar, sea_ice_pensolar_frac_direct, sea_ice_pensolar_frac_diffuse, tile_temperature, tile_snow_mass, snow_depth, snow_layer_thickness, &
&snow_layer_ice_mass, snow_layer_liq_mass, snow_layer_temp, surface_conductance, canopy_water, soil_temperature, soil_moisture, &
&unfrozen_soil_moisture, frozen_soil_moisture, tile_heat_flux, tile_moisture_flux, net_prim_prod, cos_zenith_angle, &
skyview, sw_up_tile, tile_lw_grey_albedo,&
&sw_down_surf, lw_down_surf, sw_down_blue_surf, sw_direct_blue_surf, dd_mf_cb, ozone, cf_bulk, cf_liquid, rhokm_bl, surf_interp, rhokh_bl, &
&moist_flux_bl, heat_flux_bl, gradrinr, &
&alpha1_tile, ashtf_prime_tile, dtstar_tile, fracaero_t_tile, fracaero_s_tile, z0h_tile, z0m_tile, rhokh_tile, &
&chr1p5m_tile, resfs_tile, gc_tile, canhc_tile, tile_water_extract, z0m_eff, ustar, soil_moist_avail, snow_unload_rate, &
&albedo_obs_scaling, soil_clay, soil_sand, dust_mrel, dust_flux, &
urbwrr, urbhwr, urbhgt, urbztm, urbdisp, &
rhostar, recip_l_mo_sea, t1_sd_2d, q1_sd_2d, &
&gross_prim_prod, z0h_eff
      TYPE(integer_field_type), intent(in) :: n_snow_layers, blend_height_tq, ocn_cpl_point
      INTEGER(KIND=i_def), intent(in) :: stencil_depth, ncells, ncells_halo, day_of_year, second_of_day
      REAL(KIND=r_def), intent(in) :: flux_e, flux_h
      INTEGER(KIND=i_def) nlayers
      TYPE(integer_field_proxy_type) n_snow_layers_proxy, blend_height_tq_proxy, ocn_cpl_point_proxy
      TYPE(field_proxy_type) theta_proxy, exner_in_wth_proxy, u_in_w3_proxy, v_in_w3_proxy, mr_n_proxy, mr_n_1_proxy, &
&mr_n_2_proxy, height_w3_proxy, height_wth_proxy, zh_proxy, z0msea_proxy, z0m_proxy, tile_fraction_proxy, leaf_area_index_proxy, &
&canopy_height_proxy, peak_to_trough_orog_proxy, silhouette_area_orog_proxy, soil_albedo_proxy, soil_roughness_proxy, &
&soil_moist_wilt_proxy, soil_moist_crit_proxy, soil_moist_sat_proxy, soil_thermal_cond_proxy, soil_suction_sat_proxy, &
&clapp_horn_b_proxy, soil_respiration_proxy, thermal_cond_wet_soil_proxy, &
sea_u_current_ptr_proxy, sea_v_current_ptr_proxy, sea_ice_temperature_proxy, sea_ice_conductivity_proxy, &
sea_ice_pensolar_proxy, sea_ice_pensolar_frac_direct_proxy, sea_ice_pensolar_frac_diffuse_proxy, tile_temperature_proxy, &
&tile_snow_mass_proxy, snow_depth_proxy, snow_layer_thickness_proxy, snow_layer_ice_mass_proxy, snow_layer_liq_mass_proxy, &
&snow_layer_temp_proxy, surface_conductance_proxy, canopy_water_proxy, soil_temperature_proxy, soil_moisture_proxy, &
&unfrozen_soil_moisture_proxy, frozen_soil_moisture_proxy, tile_heat_flux_proxy, tile_moisture_flux_proxy, net_prim_prod_proxy, &
&cos_zenith_angle_proxy, skyview_proxy, &
sw_up_tile_proxy, tile_lw_grey_albedo_proxy, &
sw_down_surf_proxy, lw_down_surf_proxy, sw_down_blue_surf_proxy, sw_direct_blue_surf_proxy, dd_mf_cb_proxy, &
&ozone_proxy, cf_bulk_proxy, cf_liquid_proxy, rhokm_bl_proxy, surf_interp_proxy, rhokh_bl_proxy, moist_flux_bl_proxy, &
&heat_flux_bl_proxy, gradrinr_proxy, alpha1_tile_proxy, ashtf_prime_tile_proxy, dtstar_tile_proxy, &
&fracaero_t_tile_proxy, fracaero_s_tile_proxy, &
&z0h_tile_proxy, z0m_tile_proxy, rhokh_tile_proxy, chr1p5m_tile_proxy, resfs_tile_proxy, gc_tile_proxy, canhc_tile_proxy, &
&tile_water_extract_proxy, z0m_eff_proxy, ustar_proxy, soil_moist_avail_proxy, snow_unload_rate_proxy, albedo_obs_scaling_proxy, &
&soil_clay_proxy, soil_sand_proxy, dust_mrel_proxy, dust_flux_proxy, &
urbwrr_proxy, urbhwr_proxy, urbhgt_proxy, urbztm_proxy, urbdisp_proxy, &
rhostar_proxy, recip_l_mo_sea_proxy, &
&t1_sd_2d_proxy, q1_sd_2d_proxy, gross_prim_prod_proxy, z0h_eff_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc10_dust_mrel(:,:) => null(), map_adspc1_zh(:,:) => null(), &
&map_adspc2_tile_fraction(:,:) => null(), map_adspc3_leaf_area_index(:,:) => null(), &
&map_adspc4_sea_ice_temperature(:,:) => null(), map_adspc5_snow_layer_thickness(:,:) => null(), &
&map_adspc6_soil_temperature(:,:) => null(), map_adspc7_surf_interp(:,:) => null(), map_adspc8_tile_water_extract(:,:) => null(), &
&map_adspc9_albedo_obs_scaling(:,:) => null(), map_w3(:,:) => null(), map_wtheta(:,:) => null()
      INTEGER(KIND=i_def) ndf_wtheta, undf_wtheta, ndf_w3, undf_w3, ndf_adspc1_zh, undf_adspc1_zh, ndf_adspc2_tile_fraction, &
&undf_adspc2_tile_fraction, ndf_adspc3_leaf_area_index, undf_adspc3_leaf_area_index, ndf_adspc4_sea_ice_temperature, &
&undf_adspc4_sea_ice_temperature, ndf_adspc5_snow_layer_thickness, undf_adspc5_snow_layer_thickness, ndf_adspc6_soil_temperature, &
&undf_adspc6_soil_temperature, ndf_adspc7_surf_interp, undf_adspc7_surf_interp, ndf_adspc8_tile_water_extract, &
&undf_adspc8_tile_water_extract, ndf_adspc9_albedo_obs_scaling, undf_adspc9_albedo_obs_scaling, ndf_adspc10_dust_mrel, &
&undf_adspc10_dust_mrel
      INTEGER(KIND=i_def) max_halo_depth_mesh
      TYPE(mesh_type), pointer :: mesh => null()
      INTEGER(KIND=i_def), pointer :: tile_fraction_stencil_size(:) => null()
      INTEGER(KIND=i_def), pointer :: tile_fraction_stencil_dofmap(:,:,:) => null()
      TYPE(stencil_dofmap_type), pointer :: tile_fraction_stencil_map => null()
      INTEGER(KIND=i_def), pointer :: u_in_w3_stencil_size(:) => null()
      INTEGER(KIND=i_def), pointer :: u_in_w3_stencil_dofmap(:,:,:) => null()
      TYPE(stencil_dofmap_type), pointer :: u_in_w3_stencil_map => null()

      INTEGER(KIND=i_def), pointer :: sea_u_current_ptr_stencil_size(:) => null()
      INTEGER(KIND=i_def), pointer :: sea_u_current_ptr_stencil_dofmap(:,:,:) => null()
      TYPE(stencil_dofmap_type), pointer :: sea_u_current_ptr_stencil_map => null()
      INTEGER(KIND=i_def), pointer :: sea_v_current_ptr_stencil_size(:) => null()
      INTEGER(KIND=i_def), pointer :: sea_v_current_ptr_stencil_dofmap(:,:,:) => null()
      TYPE(stencil_dofmap_type), pointer :: sea_v_current_ptr_stencil_map => null()
      !
      ! Initialise field and/or operator proxies
      !
      theta_proxy = theta%get_proxy()
      exner_in_wth_proxy = exner_in_wth%get_proxy()
      u_in_w3_proxy = u_in_w3%get_proxy()
      v_in_w3_proxy = v_in_w3%get_proxy()
      mr_n_proxy = mr_n%get_proxy()
      mr_n_1_proxy = mr_n_1%get_proxy()
      mr_n_2_proxy = mr_n_2%get_proxy()
      height_w3_proxy = height_w3%get_proxy()
      height_wth_proxy = height_wth%get_proxy()
      zh_proxy = zh%get_proxy()
      z0msea_proxy = z0msea%get_proxy()
      z0m_proxy = z0m%get_proxy()
      tile_fraction_proxy = tile_fraction%get_proxy()
      leaf_area_index_proxy = leaf_area_index%get_proxy()
      canopy_height_proxy = canopy_height%get_proxy()
      peak_to_trough_orog_proxy = peak_to_trough_orog%get_proxy()
      silhouette_area_orog_proxy = silhouette_area_orog%get_proxy()
      soil_albedo_proxy = soil_albedo%get_proxy()
      soil_roughness_proxy = soil_roughness%get_proxy()
      soil_moist_wilt_proxy = soil_moist_wilt%get_proxy()
      soil_moist_crit_proxy = soil_moist_crit%get_proxy()
      soil_moist_sat_proxy = soil_moist_sat%get_proxy()
      soil_thermal_cond_proxy = soil_thermal_cond%get_proxy()
      soil_suction_sat_proxy = soil_suction_sat%get_proxy()
      clapp_horn_b_proxy = clapp_horn_b%get_proxy()
      soil_respiration_proxy = soil_respiration%get_proxy()
      thermal_cond_wet_soil_proxy = thermal_cond_wet_soil%get_proxy()
      sea_u_current_ptr_proxy = sea_u_current_ptr%get_proxy()
      sea_v_current_ptr_proxy = sea_v_current_ptr%get_proxy()
      sea_ice_temperature_proxy = sea_ice_temperature%get_proxy()
      sea_ice_conductivity_proxy = sea_ice_conductivity%get_proxy()
      sea_ice_pensolar_proxy = sea_ice_pensolar%get_proxy()
      sea_ice_pensolar_frac_direct_proxy = sea_ice_pensolar_frac_direct%get_proxy()
      sea_ice_pensolar_frac_diffuse_proxy = sea_ice_pensolar_frac_diffuse%get_proxy()
      tile_temperature_proxy = tile_temperature%get_proxy()
      tile_snow_mass_proxy = tile_snow_mass%get_proxy()
      n_snow_layers_proxy = n_snow_layers%get_proxy()
      snow_depth_proxy = snow_depth%get_proxy()
      snow_layer_thickness_proxy = snow_layer_thickness%get_proxy()
      snow_layer_ice_mass_proxy = snow_layer_ice_mass%get_proxy()
      snow_layer_liq_mass_proxy = snow_layer_liq_mass%get_proxy()
      snow_layer_temp_proxy = snow_layer_temp%get_proxy()
      surface_conductance_proxy = surface_conductance%get_proxy()
      canopy_water_proxy = canopy_water%get_proxy()
      soil_temperature_proxy = soil_temperature%get_proxy()
      soil_moisture_proxy = soil_moisture%get_proxy()
      unfrozen_soil_moisture_proxy = unfrozen_soil_moisture%get_proxy()
      frozen_soil_moisture_proxy = frozen_soil_moisture%get_proxy()
      tile_heat_flux_proxy = tile_heat_flux%get_proxy()
      tile_moisture_flux_proxy = tile_moisture_flux%get_proxy()
      net_prim_prod_proxy = net_prim_prod%get_proxy()
      cos_zenith_angle_proxy = cos_zenith_angle%get_proxy()
      skyview_proxy = skyview%get_proxy()
      sw_up_tile_proxy = sw_up_tile%get_proxy()
      tile_lw_grey_albedo_proxy = tile_lw_grey_albedo%get_proxy()
      sw_down_surf_proxy = sw_down_surf%get_proxy()
      lw_down_surf_proxy = lw_down_surf%get_proxy()
      sw_down_blue_surf_proxy = sw_down_blue_surf%get_proxy()
      sw_direct_blue_surf_proxy = sw_direct_blue_surf%get_proxy()
      dd_mf_cb_proxy = dd_mf_cb%get_proxy()
      ozone_proxy = ozone%get_proxy()
      cf_bulk_proxy = cf_bulk%get_proxy()
      cf_liquid_proxy = cf_liquid%get_proxy()
      rhokm_bl_proxy = rhokm_bl%get_proxy()
      surf_interp_proxy = surf_interp%get_proxy()
      rhokh_bl_proxy = rhokh_bl%get_proxy()
      moist_flux_bl_proxy = moist_flux_bl%get_proxy()
      heat_flux_bl_proxy = heat_flux_bl%get_proxy()
      gradrinr_proxy = gradrinr%get_proxy()
      alpha1_tile_proxy = alpha1_tile%get_proxy()
      ashtf_prime_tile_proxy = ashtf_prime_tile%get_proxy()
      dtstar_tile_proxy = dtstar_tile%get_proxy()
      fracaero_t_tile_proxy = fracaero_t_tile%get_proxy()
      fracaero_s_tile_proxy = fracaero_s_tile%get_proxy()
      z0h_tile_proxy = z0h_tile%get_proxy()
      z0m_tile_proxy = z0m_tile%get_proxy()
      rhokh_tile_proxy = rhokh_tile%get_proxy()
      chr1p5m_tile_proxy = chr1p5m_tile%get_proxy()
      resfs_tile_proxy = resfs_tile%get_proxy()
      gc_tile_proxy = gc_tile%get_proxy()
      canhc_tile_proxy = canhc_tile%get_proxy()
      tile_water_extract_proxy = tile_water_extract%get_proxy()
      blend_height_tq_proxy = blend_height_tq%get_proxy()
      z0m_eff_proxy = z0m_eff%get_proxy()
      ustar_proxy = ustar%get_proxy()
      soil_moist_avail_proxy = soil_moist_avail%get_proxy()
      snow_unload_rate_proxy = snow_unload_rate%get_proxy()
      albedo_obs_scaling_proxy = albedo_obs_scaling%get_proxy()
      soil_clay_proxy = soil_clay%get_proxy()
      soil_sand_proxy = soil_sand%get_proxy()
      dust_mrel_proxy = dust_mrel%get_proxy()
      dust_flux_proxy = dust_flux%get_proxy()
      urbwrr_proxy = urbwrr%get_proxy()
      urbhwr_proxy = urbhwr%get_proxy()
      urbhgt_proxy = urbhgt%get_proxy()
      urbztm_proxy = urbztm%get_proxy()
      urbdisp_proxy = urbdisp%get_proxy()
      rhostar_proxy = rhostar%get_proxy()
      recip_l_mo_sea_proxy = recip_l_mo_sea%get_proxy()
      t1_sd_2d_proxy = t1_sd_2d%get_proxy()
      q1_sd_2d_proxy = q1_sd_2d%get_proxy()
      gross_prim_prod_proxy = gross_prim_prod%get_proxy()
      z0h_eff_proxy = z0h_eff%get_proxy()
      ocn_cpl_point_proxy = ocn_cpl_point%get_proxy()
      !
      ! Initialise number of layers
      !
      nlayers = theta_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => theta_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Initialise stencil dofmaps
      !
      u_in_w3_stencil_map => u_in_w3_proxy%vspace%get_stencil_dofmap(STENCIL_REGION,stencil_depth)
      u_in_w3_stencil_dofmap => u_in_w3_stencil_map%get_whole_dofmap()
      u_in_w3_stencil_size => u_in_w3_stencil_map%get_stencil_sizes()
      tile_fraction_stencil_map => tile_fraction_proxy%vspace%get_stencil_dofmap(STENCIL_REGION,stencil_depth)
      tile_fraction_stencil_dofmap => tile_fraction_stencil_map%get_whole_dofmap()
      tile_fraction_stencil_size => tile_fraction_stencil_map%get_stencil_sizes()
      sea_u_current_ptr_stencil_map => sea_u_current_ptr_proxy%vspace%get_stencil_dofmap(STENCIL_REGION,stencil_depth)
      sea_u_current_ptr_stencil_dofmap => sea_u_current_ptr_stencil_map%get_whole_dofmap()
      sea_u_current_ptr_stencil_size => sea_u_current_ptr_stencil_map%get_stencil_sizes()
      sea_v_current_ptr_stencil_map => sea_v_current_ptr_proxy%vspace%get_stencil_dofmap(STENCIL_REGION,stencil_depth)
      sea_v_current_ptr_stencil_dofmap => sea_v_current_ptr_stencil_map%get_whole_dofmap()
      sea_v_current_ptr_stencil_size => sea_v_current_ptr_stencil_map%get_stencil_sizes()
      !
      ! Look-up dofmaps for each function space
      !
      map_wtheta => theta_proxy%vspace%get_whole_dofmap()
      map_w3 => u_in_w3_proxy%vspace%get_whole_dofmap()
      map_adspc1_zh => zh_proxy%vspace%get_whole_dofmap()
      map_adspc2_tile_fraction => tile_fraction_proxy%vspace%get_whole_dofmap()
      map_adspc3_leaf_area_index => leaf_area_index_proxy%vspace%get_whole_dofmap()
      map_adspc4_sea_ice_temperature => sea_ice_temperature_proxy%vspace%get_whole_dofmap()
      map_adspc5_snow_layer_thickness => snow_layer_thickness_proxy%vspace%get_whole_dofmap()
      map_adspc6_soil_temperature => soil_temperature_proxy%vspace%get_whole_dofmap()
      map_adspc7_surf_interp => surf_interp_proxy%vspace%get_whole_dofmap()
      map_adspc8_tile_water_extract => tile_water_extract_proxy%vspace%get_whole_dofmap()
      map_adspc9_albedo_obs_scaling => albedo_obs_scaling_proxy%vspace%get_whole_dofmap()
      map_adspc10_dust_mrel => dust_mrel_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for wtheta
      !
      ndf_wtheta = theta_proxy%vspace%get_ndf()
      undf_wtheta = theta_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for w3
      !
      ndf_w3 = u_in_w3_proxy%vspace%get_ndf()
      undf_w3 = u_in_w3_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc1_zh
      !
      ndf_adspc1_zh = zh_proxy%vspace%get_ndf()
      undf_adspc1_zh = zh_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc2_tile_fraction
      !
      ndf_adspc2_tile_fraction = tile_fraction_proxy%vspace%get_ndf()
      undf_adspc2_tile_fraction = tile_fraction_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc3_leaf_area_index
      !
      ndf_adspc3_leaf_area_index = leaf_area_index_proxy%vspace%get_ndf()
      undf_adspc3_leaf_area_index = leaf_area_index_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc4_sea_ice_temperature
      !
      ndf_adspc4_sea_ice_temperature = sea_ice_temperature_proxy%vspace%get_ndf()
      undf_adspc4_sea_ice_temperature = sea_ice_temperature_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc5_snow_layer_thickness
      !
      ndf_adspc5_snow_layer_thickness = snow_layer_thickness_proxy%vspace%get_ndf()
      undf_adspc5_snow_layer_thickness = snow_layer_thickness_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc6_soil_temperature
      !
      ndf_adspc6_soil_temperature = soil_temperature_proxy%vspace%get_ndf()
      undf_adspc6_soil_temperature = soil_temperature_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc7_surf_interp
      !
      ndf_adspc7_surf_interp = surf_interp_proxy%vspace%get_ndf()
      undf_adspc7_surf_interp = surf_interp_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc8_tile_water_extract
      !
      ndf_adspc8_tile_water_extract = tile_water_extract_proxy%vspace%get_ndf()
      undf_adspc8_tile_water_extract = tile_water_extract_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc9_albedo_obs_scaling
      !
      ndf_adspc9_albedo_obs_scaling = albedo_obs_scaling_proxy%vspace%get_ndf()
      undf_adspc9_albedo_obs_scaling = albedo_obs_scaling_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc10_dust_mrel
      !
      ndf_adspc10_dust_mrel = dust_mrel_proxy%vspace%get_ndf()
      undf_adspc10_dust_mrel = dust_mrel_proxy%vspace%get_undf()
      !
      ! Call kernels and communication routines
      !
      IF (u_in_w3_proxy%is_dirty(depth=stencil_depth)) THEN
        CALL u_in_w3_proxy%halo_exchange(depth=stencil_depth)
      END IF
      !
      IF (v_in_w3_proxy%is_dirty(depth=stencil_depth)) THEN
        CALL v_in_w3_proxy%halo_exchange(depth=stencil_depth)
      END IF
      !
      IF (tile_fraction_proxy%is_dirty(depth=stencil_depth)) THEN
        CALL tile_fraction_proxy%halo_exchange(depth=stencil_depth)
      END IF
      !
      IF (sea_u_current_ptr_proxy%is_dirty(depth=stencil_depth)) THEN
        CALL sea_u_current_ptr_proxy%halo_exchange(depth=stencil_depth)
      END IF
      !
      IF (sea_v_current_ptr_proxy%is_dirty(depth=stencil_depth)) THEN
        CALL sea_v_current_ptr_proxy%halo_exchange(depth=stencil_depth)
      END IF
      !
      CALL jules_exp_code(nlayers, ncells, ncells_halo, theta_proxy%data, exner_in_wth_proxy%data, u_in_w3_proxy%data, u_in_w3_stencil_size, &
&u_in_w3_stencil_dofmap, v_in_w3_proxy%data, u_in_w3_stencil_size, u_in_w3_stencil_dofmap, &
&mr_n_proxy%data, mr_n_1_proxy%data, mr_n_2_proxy%data, height_w3_proxy%data, height_wth_proxy%data, zh_proxy%data, &
&z0msea_proxy%data, z0m_proxy%data, tile_fraction_proxy%data, tile_fraction_stencil_size, &
&tile_fraction_stencil_dofmap, leaf_area_index_proxy%data, canopy_height_proxy%data, peak_to_trough_orog_proxy%data, &
&silhouette_area_orog_proxy%data, soil_albedo_proxy%data, soil_roughness_proxy%data, soil_moist_wilt_proxy%data, &
&soil_moist_crit_proxy%data, soil_moist_sat_proxy%data, soil_thermal_cond_proxy%data, soil_suction_sat_proxy%data, &
&clapp_horn_b_proxy%data, soil_respiration_proxy%data, thermal_cond_wet_soil_proxy%data,  &
sea_u_current_ptr_proxy%data,sea_u_current_ptr_stencil_size,sea_u_current_ptr_stencil_dofmap, &
sea_v_current_ptr_proxy%data,sea_v_current_ptr_stencil_size,sea_v_current_ptr_stencil_dofmap, &
sea_ice_temperature_proxy%data, sea_ice_conductivity_proxy%data, sea_ice_pensolar_proxy%data, sea_ice_pensolar_frac_direct_proxy%data, &
sea_ice_pensolar_frac_diffuse_proxy%data, &
&tile_temperature_proxy%data, tile_snow_mass_proxy%data, n_snow_layers_proxy%data, snow_depth_proxy%data, &
&snow_layer_thickness_proxy%data, snow_layer_ice_mass_proxy%data, snow_layer_liq_mass_proxy%data, snow_layer_temp_proxy%data, &
&surface_conductance_proxy%data, canopy_water_proxy%data, soil_temperature_proxy%data, soil_moisture_proxy%data, &
&unfrozen_soil_moisture_proxy%data, frozen_soil_moisture_proxy%data, tile_heat_flux_proxy%data, tile_moisture_flux_proxy%data, &
&net_prim_prod_proxy%data, cos_zenith_angle_proxy%data, skyview_proxy%data, &
sw_up_tile_proxy%data, tile_lw_grey_albedo_proxy%data, sw_down_surf_proxy%data, lw_down_surf_proxy%data, &
&sw_down_blue_surf_proxy%data, sw_direct_blue_surf_proxy%data, dd_mf_cb_proxy%data, ozone_proxy%data, cf_bulk_proxy%data, cf_liquid_proxy%data, &
&rhokm_bl_proxy%data, surf_interp_proxy%data, rhokh_bl_proxy%data, moist_flux_bl_proxy%data, heat_flux_bl_proxy%data, &
&gradrinr_proxy%data, alpha1_tile_proxy%data, ashtf_prime_tile_proxy%data, dtstar_tile_proxy%data, &
&fracaero_t_tile_proxy%data, fracaero_s_tile_proxy%data, &
&z0h_tile_proxy%data, z0m_tile_proxy%data, rhokh_tile_proxy%data, chr1p5m_tile_proxy%data, resfs_tile_proxy%data, &
&gc_tile_proxy%data, canhc_tile_proxy%data, tile_water_extract_proxy%data, blend_height_tq_proxy%data, z0m_eff_proxy%data, &
&ustar_proxy%data, soil_moist_avail_proxy%data, snow_unload_rate_proxy%data, albedo_obs_scaling_proxy%data, soil_clay_proxy%data, &
&soil_sand_proxy%data, dust_mrel_proxy%data, dust_flux_proxy%data, day_of_year, second_of_day, &
flux_e, flux_h, &
urbwrr_proxy%data, urbhwr_proxy%data, urbhgt_proxy%data, urbztm_proxy%data, &
urbdisp_proxy%data, &
rhostar_proxy%data, recip_l_mo_sea_proxy%data, &
&t1_sd_2d_proxy%data, q1_sd_2d_proxy%data, gross_prim_prod_proxy%data, &
z0h_eff_proxy%data, ocn_cpl_point_proxy%data, ndf_wtheta, &
&undf_wtheta, map_wtheta, ndf_w3, undf_w3, map_w3, ndf_adspc1_zh, undf_adspc1_zh, map_adspc1_zh, &
&ndf_adspc2_tile_fraction, undf_adspc2_tile_fraction, map_adspc2_tile_fraction, ndf_adspc3_leaf_area_index, &
&undf_adspc3_leaf_area_index, map_adspc3_leaf_area_index, ndf_adspc4_sea_ice_temperature, undf_adspc4_sea_ice_temperature, &
&map_adspc4_sea_ice_temperature, ndf_adspc5_snow_layer_thickness, undf_adspc5_snow_layer_thickness, &
&map_adspc5_snow_layer_thickness, ndf_adspc6_soil_temperature, undf_adspc6_soil_temperature, &
&map_adspc6_soil_temperature, ndf_adspc7_surf_interp, undf_adspc7_surf_interp, map_adspc7_surf_interp, &
&ndf_adspc8_tile_water_extract, undf_adspc8_tile_water_extract, map_adspc8_tile_water_extract, &
&ndf_adspc9_albedo_obs_scaling, undf_adspc9_albedo_obs_scaling, map_adspc9_albedo_obs_scaling, ndf_adspc10_dust_mrel, &
&undf_adspc10_dust_mrel, map_adspc10_dust_mrel)
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL z0msea_proxy%set_dirty()
      CALL z0m_proxy%set_dirty()
      CALL soil_respiration_proxy%set_dirty()
      CALL thermal_cond_wet_soil_proxy%set_dirty()
      CALL tile_temperature_proxy%set_dirty()
      CALL surface_conductance_proxy%set_dirty()
      CALL tile_heat_flux_proxy%set_dirty()
      CALL tile_moisture_flux_proxy%set_dirty()
      CALL net_prim_prod_proxy%set_dirty()
      CALL rhokm_bl_proxy%set_dirty()
      CALL surf_interp_proxy%set_dirty()
      CALL rhokh_bl_proxy%set_dirty()
      CALL moist_flux_bl_proxy%set_dirty()
      CALL heat_flux_bl_proxy%set_dirty()
      CALL gradrinr_proxy%set_dirty()
      CALL alpha1_tile_proxy%set_dirty()
      CALL ashtf_prime_tile_proxy%set_dirty()
      CALL dtstar_tile_proxy%set_dirty()
      CALL fracaero_t_tile_proxy%set_dirty()
      CALL fracaero_s_tile_proxy%set_dirty()
      CALL z0h_tile_proxy%set_dirty()
      CALL z0m_tile_proxy%set_dirty()
      CALL rhokh_tile_proxy%set_dirty()
      CALL chr1p5m_tile_proxy%set_dirty()
      CALL resfs_tile_proxy%set_dirty()
      CALL gc_tile_proxy%set_dirty()
      CALL canhc_tile_proxy%set_dirty()
      CALL tile_water_extract_proxy%set_dirty()
      CALL blend_height_tq_proxy%set_dirty()
      CALL z0m_eff_proxy%set_dirty()
      CALL ustar_proxy%set_dirty()
      CALL soil_moist_avail_proxy%set_dirty()
      CALL snow_unload_rate_proxy%set_dirty()
      CALL dust_flux_proxy%set_dirty()
      CALL rhostar_proxy%set_dirty()
      CALL recip_l_mo_sea_proxy%set_dirty()
      CALL t1_sd_2d_proxy%set_dirty()
      CALL q1_sd_2d_proxy%set_dirty()
      CALL gross_prim_prod_proxy%set_dirty()
      CALL z0h_eff_proxy%set_dirty()
      !
      !
    END SUBROUTINE invoke_jules_exp_kernel_type

  !---------------------------------------------------------------------
  !> Contains the PSy-layer to build the pressure level diagnostics
  !> These require passing an array "plevs" into each kernel
  !> which is currently unsupported by PSyclone
  !> see https://github.com/stfc/PSyclone/issues/1312
  !> Hence this module could be removed once the PSyclone ticket is
  !> completed
    SUBROUTINE invoke_heaviside_kernel_type(exner_wth, nplev, plevs, plev_heaviside, p_zero, kappa)
      USE heaviside_kernel_mod, ONLY: heaviside_code
      USE mesh_mod, ONLY: mesh_type
      REAL(KIND=r_def), intent(in) :: p_zero, kappa
      INTEGER(KIND=i_def), intent(in) :: nplev
      REAL(KIND=r_def), intent(in) :: plevs(nplev)
      TYPE(field_type), intent(in) :: exner_wth, plev_heaviside
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers
      REAL(KIND=r_def), pointer, dimension(:) :: plev_heaviside_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: exner_wth_data => null()
      TYPE(field_proxy_type) exner_wth_proxy, plev_heaviside_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_exner_wth(:,:) => null(), map_adspc2_plev_heaviside(:,:) => null()
      INTEGER(KIND=i_def) ndf_adspc1_exner_wth, undf_adspc1_exner_wth, ndf_adspc2_plev_heaviside, undf_adspc2_plev_heaviside
      INTEGER(KIND=i_def) max_halo_depth_mesh
      TYPE(mesh_type), pointer :: mesh => null()
      !
      ! Initialise field and/or operator proxies
      !
      exner_wth_proxy = exner_wth%get_proxy()
      exner_wth_data => exner_wth_proxy%data
      plev_heaviside_proxy = plev_heaviside%get_proxy()
      plev_heaviside_data => plev_heaviside_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers = exner_wth_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => exner_wth_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Look-up dofmaps for each function space
      !
      map_adspc1_exner_wth => exner_wth_proxy%vspace%get_whole_dofmap()
      map_adspc2_plev_heaviside => plev_heaviside_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for adspc1_exner_wth
      !
      ndf_adspc1_exner_wth = exner_wth_proxy%vspace%get_ndf()
      undf_adspc1_exner_wth = exner_wth_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc2_plev_heaviside
      !
      ndf_adspc2_plev_heaviside = plev_heaviside_proxy%vspace%get_ndf()
      undf_adspc2_plev_heaviside = plev_heaviside_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = mesh%get_last_edge_cell()
      !
      ! Call kernels and communication routines
      !
      DO cell=loop0_start,loop0_stop
        !
        CALL heaviside_code(nlayers, exner_wth_data, nplev, plevs, plev_heaviside_data, p_zero, kappa, ndf_adspc1_exner_wth, &
&undf_adspc1_exner_wth, map_adspc1_exner_wth(:,cell), ndf_adspc2_plev_heaviside, undf_adspc2_plev_heaviside, &
&map_adspc2_plev_heaviside(:,cell))
      END DO
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL plev_heaviside_proxy%set_dirty()
      !
      !
    END SUBROUTINE invoke_heaviside_kernel_type

  !---------------------------------------------------------------------
  !> Contains the PSy-layer to build the pressure level diagnostics
  !> These require passing an array "plevs" into each kernel
  !> which is currently unsupported by PSyclone
  !> see https://github.com/stfc/PSyclone/issues/1312
  !> Hence this module could be removed once the PSyclone ticket is
  !> completed
    SUBROUTINE invoke_temp_on_pres_kernel_type(temp, exner_wth, &
         height_wth, nplev, plevs, plev_temp, p_zero, kappa, ex_power)
      USE temp_on_pres_kernel_mod, ONLY: temp_on_pres_code
      USE mesh_mod, ONLY: mesh_type
      REAL(KIND=r_def), intent(in) :: p_zero, kappa, ex_power
      INTEGER(KIND=i_def), intent(in) :: nplev
      REAL(KIND=r_def), intent(in) :: plevs(nplev)
      TYPE(field_type), intent(in) :: temp, exner_wth, height_wth, plev_temp
      INTEGER(KIND=i_def) cell
      INTEGER df
      INTEGER(KIND=i_def) loop1_start, loop1_stop
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers
      REAL(KIND=r_def), pointer, dimension(:) :: plev_temp_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: height_wth_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: exner_wth_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: temp_data => null()
      TYPE(field_proxy_type) temp_proxy, exner_wth_proxy, height_wth_proxy, plev_temp_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_temp(:,:) => null(), map_adspc2_plev_temp(:,:) => null(), map_wtheta(:,:) => null()
      INTEGER(KIND=i_def) ndf_aspc1_temp, undf_aspc1_temp, ndf_adspc1_temp, undf_adspc1_temp, ndf_wtheta, undf_wtheta, &
&ndf_adspc2_plev_temp, undf_adspc2_plev_temp
      INTEGER(KIND=i_def) max_halo_depth_mesh
      TYPE(mesh_type), pointer :: mesh => null()
      !
      ! Initialise field and/or operator proxies
      !
      temp_proxy = temp%get_proxy()
      temp_data => temp_proxy%data
      exner_wth_proxy = exner_wth%get_proxy()
      exner_wth_data => exner_wth_proxy%data
      height_wth_proxy = height_wth%get_proxy()
      height_wth_data => height_wth_proxy%data
      plev_temp_proxy = plev_temp%get_proxy()
      plev_temp_data => plev_temp_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers = temp_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => temp_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Look-up dofmaps for each function space
      !
      map_adspc1_temp => temp_proxy%vspace%get_whole_dofmap()
      map_wtheta => height_wth_proxy%vspace%get_whole_dofmap()
      map_adspc2_plev_temp => plev_temp_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for aspc1_temp
      !
      ndf_aspc1_temp = temp_proxy%vspace%get_ndf()
      undf_aspc1_temp = temp_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc1_temp
      !
      ndf_adspc1_temp = temp_proxy%vspace%get_ndf()
      undf_adspc1_temp = temp_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for wtheta
      !
      ndf_wtheta = height_wth_proxy%vspace%get_ndf()
      undf_wtheta = height_wth_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc2_plev_temp
      !
      ndf_adspc2_plev_temp = plev_temp_proxy%vspace%get_ndf()
      undf_adspc2_plev_temp = plev_temp_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop1_start = 1
      loop1_stop = mesh%get_last_edge_cell()
      !
      ! Call kernels and communication routines
      !
      DO cell=loop1_start,loop1_stop
        !
        CALL temp_on_pres_code(nlayers, temp_data, exner_wth_data, height_wth_data, nplev, plevs, plev_temp_data, p_zero, kappa, &
&ex_power, ndf_adspc1_temp, undf_adspc1_temp, map_adspc1_temp(:,cell), ndf_wtheta, undf_wtheta, map_wtheta(:,cell), &
&ndf_adspc2_plev_temp, undf_adspc2_plev_temp, map_adspc2_plev_temp(:,cell))
      END DO
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL plev_temp_proxy%set_dirty()
      !
      !
    END SUBROUTINE invoke_temp_on_pres_kernel_type

  !---------------------------------------------------------------------
  !> Contains the PSy-layer to build the pressure level diagnostics
  !> These require passing an array "plevs" into each kernel
  !> which is currently unsupported by PSyclone
  !> see https://github.com/stfc/PSyclone/issues/1312
  !> Hence this module could be removed once the PSyclone ticket is
  !> completed
    SUBROUTINE invoke_pres_interp_kernel_type(u_in_w3, exner_w3, nplev, plevs, plev_u, p_zero, kappa)
      USE pres_interp_kernel_mod, ONLY: pres_interp_code
      USE mesh_mod, ONLY: mesh_type
      REAL(KIND=r_def), intent(in) :: p_zero, kappa
      INTEGER(KIND=i_def), intent(in) :: nplev
      REAL(KIND=r_def), intent(in) :: plevs(nplev)
      TYPE(field_type), intent(in) :: u_in_w3, exner_w3, plev_u
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers
      REAL(KIND=r_def), pointer, dimension(:) :: plev_u_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: exner_w3_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: u_in_w3_data => null()
      TYPE(field_proxy_type) u_in_w3_proxy, exner_w3_proxy, plev_u_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_u_in_w3(:,:) => null(), map_adspc2_plev_u(:,:) => null()
      INTEGER(KIND=i_def) ndf_adspc1_u_in_w3, undf_adspc1_u_in_w3, ndf_adspc2_plev_u, undf_adspc2_plev_u
      INTEGER(KIND=i_def) max_halo_depth_mesh
      TYPE(mesh_type), pointer :: mesh => null()
      !
      ! Initialise field and/or operator proxies
      !
      u_in_w3_proxy = u_in_w3%get_proxy()
      u_in_w3_data => u_in_w3_proxy%data
      exner_w3_proxy = exner_w3%get_proxy()
      exner_w3_data => exner_w3_proxy%data
      plev_u_proxy = plev_u%get_proxy()
      plev_u_data => plev_u_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers = u_in_w3_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => u_in_w3_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Look-up dofmaps for each function space
      !
      map_adspc1_u_in_w3 => u_in_w3_proxy%vspace%get_whole_dofmap()
      map_adspc2_plev_u => plev_u_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for adspc1_u_in_w3
      !
      ndf_adspc1_u_in_w3 = u_in_w3_proxy%vspace%get_ndf()
      undf_adspc1_u_in_w3 = u_in_w3_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc2_plev_u
      !
      ndf_adspc2_plev_u = plev_u_proxy%vspace%get_ndf()
      undf_adspc2_plev_u = plev_u_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = mesh%get_last_edge_cell()
      !
      ! Call kernels and communication routines
      !
      DO cell=loop0_start,loop0_stop
        !
        CALL pres_interp_code(nlayers, u_in_w3_data, exner_w3_data, nplev, plevs, plev_u_data, p_zero, kappa, ndf_adspc1_u_in_w3, &
&undf_adspc1_u_in_w3, map_adspc1_u_in_w3(:,cell), ndf_adspc2_plev_u, undf_adspc2_plev_u, map_adspc2_plev_u(:,cell))
      END DO
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL plev_u_proxy%set_dirty()
      !
      !
    END SUBROUTINE invoke_pres_interp_kernel_type

  !---------------------------------------------------------------------
  !> Contains the PSy-layer to build the pressure level diagnostics
  !> These require passing an array "plevs" into each kernel
  !> which is currently unsupported by PSyclone
  !> see https://github.com/stfc/PSyclone/issues/1312
  !> Hence this module could be removed once the PSyclone ticket is
  !> completed
    SUBROUTINE invoke_geo_on_pres_kernel_type(height_w3, exner_w3, theta_wth, height_wth, exner_wth, nplev, plevs, plev_geopot, &
&p_zero, kappa, cp, gravity, ex_power)
      USE geo_on_pres_kernel_mod, ONLY: geo_on_pres_code
      USE mesh_mod, ONLY: mesh_type
      REAL(KIND=r_def), intent(in) :: p_zero, kappa, cp, gravity, ex_power
      INTEGER(KIND=i_def), intent(in) :: nplev
      REAL(KIND=r_def), intent(in) :: plevs(nplev)
      TYPE(field_type), intent(in) :: height_w3, exner_w3, theta_wth, height_wth, exner_wth, plev_geopot
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers
      REAL(KIND=r_def), pointer, dimension(:) :: plev_geopot_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: exner_wth_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: height_wth_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: theta_wth_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: exner_w3_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: height_w3_data => null()
      TYPE(field_proxy_type) height_w3_proxy, exner_w3_proxy, theta_wth_proxy, height_wth_proxy, exner_wth_proxy, plev_geopot_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_height_w3(:,:) => null(), map_adspc2_plev_geopot(:,:) => null(), &
&map_wtheta(:,:) => null()
      INTEGER(KIND=i_def) ndf_adspc1_height_w3, undf_adspc1_height_w3, ndf_wtheta, undf_wtheta, ndf_adspc2_plev_geopot, &
&undf_adspc2_plev_geopot
      INTEGER(KIND=i_def) max_halo_depth_mesh
      TYPE(mesh_type), pointer :: mesh => null()
      !
      ! Initialise field and/or operator proxies
      !
      height_w3_proxy = height_w3%get_proxy()
      height_w3_data => height_w3_proxy%data
      exner_w3_proxy = exner_w3%get_proxy()
      exner_w3_data => exner_w3_proxy%data
      theta_wth_proxy = theta_wth%get_proxy()
      theta_wth_data => theta_wth_proxy%data
      height_wth_proxy = height_wth%get_proxy()
      height_wth_data => height_wth_proxy%data
      exner_wth_proxy = exner_wth%get_proxy()
      exner_wth_data => exner_wth_proxy%data
      plev_geopot_proxy = plev_geopot%get_proxy()
      plev_geopot_data => plev_geopot_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers = height_w3_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => height_w3_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Look-up dofmaps for each function space
      !
      map_adspc1_height_w3 => height_w3_proxy%vspace%get_whole_dofmap()
      map_wtheta => theta_wth_proxy%vspace%get_whole_dofmap()
      map_adspc2_plev_geopot => plev_geopot_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for adspc1_height_w3
      !
      ndf_adspc1_height_w3 = height_w3_proxy%vspace%get_ndf()
      undf_adspc1_height_w3 = height_w3_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for wtheta
      !
      ndf_wtheta = theta_wth_proxy%vspace%get_ndf()
      undf_wtheta = theta_wth_proxy%vspace%get_undf()
      !
      ! Initialise number of DoFs for adspc2_plev_geopot
      !
      ndf_adspc2_plev_geopot = plev_geopot_proxy%vspace%get_ndf()
      undf_adspc2_plev_geopot = plev_geopot_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = mesh%get_last_edge_cell()
      !
      ! Call kernels and communication routines
      !
      DO cell=loop0_start,loop0_stop
        !
        CALL geo_on_pres_code(nlayers, height_w3_data, exner_w3_data, theta_wth_data, height_wth_data, exner_wth_data, nplev, &
&plevs, plev_geopot_data, p_zero, kappa, cp, gravity, ex_power, ndf_adspc1_height_w3, undf_adspc1_height_w3, &
&map_adspc1_height_w3(:,cell), ndf_wtheta, undf_wtheta, map_wtheta(:,cell), ndf_adspc2_plev_geopot, undf_adspc2_plev_geopot, &
&map_adspc2_plev_geopot(:,cell))
      END DO
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL plev_geopot_proxy%set_dirty()
      !
      !
    END SUBROUTINE invoke_geo_on_pres_kernel_type
  !---------------------------------------------------------------------
  !> Contains the PSy-layer to build the pressure level diagnostics
  !> These require passing an array "plevs" into each kernel
  !> which is currently unsupported by PSyclone
  !> see https://github.com/stfc/PSyclone/issues/1312
  !> Hence this module could be removed once the PSyclone ticket is
  !> completed
    SUBROUTINE invoke_thetaw_kernel_type(plev_thetaw, plev_temp, plev_qv, nplev, plevs)
      USE thetaw_kernel_mod, ONLY: thetaw_code
      USE mesh_mod, ONLY: mesh_type
      INTEGER(KIND=i_def), intent(in) :: nplev
      TYPE(field_type), intent(in) :: plev_temp, plev_qv, plev_thetaw
      REAL(KIND=r_def), intent(in) :: plevs(nplev)
      INTEGER(KIND=i_def) cell
      INTEGER(KIND=i_def) loop0_start, loop0_stop
      INTEGER(KIND=i_def) nlayers_plev_thetaw
      REAL(KIND=r_def), pointer, dimension(:) :: plev_thetaw_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: plev_qv_data => null()
      REAL(KIND=r_def), pointer, dimension(:) :: plev_temp_data => null()
      TYPE(field_proxy_type) plev_qv_proxy, plev_thetaw_proxy, plev_temp_proxy
      INTEGER(KIND=i_def), pointer :: map_adspc1_plev_thetaw(:,:) => null()
      INTEGER(KIND=i_def) ndf_adspc1_plev_thetaw, undf_adspc1_plev_thetaw
      INTEGER(KIND=i_def) max_halo_depth_mesh
      TYPE(mesh_type), pointer :: mesh => null()
      !
      ! Initialise field and/or operator proxies
      !
      plev_qv_proxy = plev_qv%get_proxy()
      plev_qv_data => plev_qv_proxy%data
      plev_temp_proxy = plev_temp%get_proxy()
      plev_temp_data => plev_temp_proxy%data
      plev_thetaw_proxy = plev_thetaw%get_proxy()
      plev_thetaw_data => plev_thetaw_proxy%data
      !
      ! Initialise number of layers
      !
      nlayers_plev_thetaw = plev_thetaw_proxy%vspace%get_nlayers()
      !
      ! Create a mesh object
      !
      mesh => plev_thetaw_proxy%vspace%get_mesh()
      max_halo_depth_mesh = mesh%get_halo_depth()
      !
      ! Look-up dofmaps for each function space
      !
      map_adspc1_plev_thetaw => plev_thetaw_proxy%vspace%get_whole_dofmap()
      !
      ! Initialise number of DoFs for adspc1_plev_thetaw
      !
      ndf_adspc1_plev_thetaw = plev_thetaw_proxy%vspace%get_ndf()
      undf_adspc1_plev_thetaw = plev_thetaw_proxy%vspace%get_undf()
      !
      ! Set-up all of the loop bounds
      !
      loop0_start = 1
      loop0_stop = mesh%get_last_edge_cell()
      !
      ! Call kernels and communication routines
      !
      DO cell = loop0_start, loop0_stop, 1
        CALL thetaw_code(nlayers_plev_thetaw, plev_thetaw_data, plev_temp_data, plev_qv_data, nplev, plevs, ndf_adspc1_plev_thetaw, &
&undf_adspc1_plev_thetaw, map_adspc1_plev_thetaw(:,cell))
      END DO
      !
      ! Set halos dirty/clean for fields modified in the above loop
      !
      CALL plev_thetaw_proxy%set_dirty()
      !
      !
    END SUBROUTINE invoke_thetaw_kernel_type
end module psykal_lite_phys_mod
