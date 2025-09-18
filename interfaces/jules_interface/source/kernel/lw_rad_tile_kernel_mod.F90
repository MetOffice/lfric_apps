!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Interface to Jules for LW surface tile radiative properties

module lw_rad_tile_kernel_mod

use argument_mod,      only : arg_type,                  &
                              GH_FIELD, GH_SCALAR,       &
                              GH_REAL, GH_INTEGER,       &
                              GH_READ, GH_WRITE,         &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              ANY_DISCONTINUOUS_SPACE_3, &
                              ANY_DISCONTINUOUS_SPACE_4, &
                              DOMAIN
use constants_mod,     only : r_def, i_def, r_um, i_um, l_def
use kernel_mod,        only : kernel_type

implicit none

private

public :: lw_rad_tile_kernel_type
public :: lw_rad_tile_code

!------------------------------------------------------------------------------
! Public types
!------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, extends(kernel_type) :: lw_rad_tile_kernel_type
  private
  type(arg_type) :: meta_args(11) = (/                                    &
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! tile_lw_albedo
    arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), & ! tile_lw_grey_albedo
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! tile_fraction
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! tile_temperature
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! snow_tile
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! leaf_area_index
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! canopy_height
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! soil_roughness
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! urbztm
    arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! urbemisc
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ                            )  & ! n_band
    /)
  integer :: operates_on = DOMAIN
contains
  procedure, nopass :: lw_rad_tile_code
end type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @param[in]     nlayers                Number of layers
!> @param[in,out] tile_lw_albedo         LW tile albedos
!> @param[in,out] tile_lw_grey_albedo    LW tile grey albedos
!> @param[in]     tile_fraction          Surface tile fractions
!> @param[in]     tile_temperature       Surface tile temperature
!> @param[in]     tile_snow_mass         Snow mass on tiles (kg/m2)
!> @param[in]     leaf_area_index        Leaf area index on plant types
!> @param[in]     canopy_height          Canopy height of plant types (m)
!> @param[in]     soil_roughness         Bare soil roughness length (m)
!> @param[in]     urbztm                 Urban effective roughness length (m)
!> @param[in]     urbemisc               Urban canyon emissivity
!> @param[in]     n_band                 Number of spectral bands
!> @param[in]     ndf_lw_tile            DOFs per cell for tiles and lw bands
!> @param[in]     undf_lw_tile           Total DOFs for tiles and lw bands
!> @param[in]     map_lw_tile            Dofmap for cell at the base of the column
!> @param[in]     ndf_tile               Number of DOFs per cell for tiles
!> @param[in]     undf_tile              Number of total DOFs for tiles
!> @param[in]     map_tile               Dofmap for cell at the base of the column
!> @param[in]     ndf_pft                Number of DOFs per cell for pfts
!> @param[in]     undf_pft               Number of total DOFs for pfts
!> @param[in]     map_pft                Dofmap for cell at the base of the column
!> @param[in]     ndf_2d                 Number of DOFs per cell for 2d
!> @param[in]     undf_2d                Number of total DOFs for 2d
!> @param[in]     map_2d                 Dofmap for cell at the base of the column
subroutine lw_rad_tile_code(nlayers, seg_len,                       &
                            tile_lw_albedo,                         &
                            tile_lw_grey_albedo,                    &
                            tile_fraction,                          &
                            tile_temperature,                       &
                            tile_snow_mass,                         &
                            leaf_area_index,                        &
                            canopy_height,                          &
                            soil_roughness,                         &
                            urbztm,                                 &
                            urbemisc,                               &
                            n_band,                                 &
                            ndf_lw_tile, undf_lw_tile, map_lw_tile, &
                            ndf_tile, undf_tile, map_tile,          &
                            ndf_pft, undf_pft, map_pft,             &
                            ndf_2d, undf_2d, map_2d)

  use socrates_init_mod, only: n_band_exclude, index_exclude, &
                               wavelength_short, wavelength_long
  use jules_control_init_mod, only: n_surf_tile, n_land_tile, n_sea_tile, &
                                    n_sea_ice_tile, first_sea_tile,       &
                                    first_sea_ice_tile
  use jules_surface_types_mod, only: soil, npft, ntype, nnpft, ice, urban_canyon
  use specemis_mod, only: specemis
  use surface_config_mod, only: emis_method_sea, emis_method_sea_fixed,   &
                                emis_method_sea_feldman,                  &
                                emis_method_sea_iremis,                   &
                                emis_method_soil, emis_method_soil_fixed, &
                                emis_method_soil_feldman_desert
  use ancil_info,               only: ainfo_type, ainfo_data_type,           &
                                      ancil_info_assoc, ancil_info_alloc,    &
                                      ancil_info_dealloc, ancil_info_nullify
  use prognostics,              only: progs_data_type, progs_type,           &
                                      prognostics_alloc, prognostics_assoc,  &
                                      prognostics_dealloc, prognostics_nullify
  use p_s_parms,                only: psparms_type, psparms_data_type,       &
                                      psparms_alloc, psparms_assoc,          &
                                      psparms_dealloc, psparms_nullify
  use urban_param_mod,          only: emiss, urban_param_type,               &
                                      urban_param_data_type,                 &
                                      urban_param_assoc, urban_param_alloc,  &
                                      urban_param_dealloc, urban_param_nullify
  use ancil_info, only: nsoilt, dim_cslayer, nmasst
  use atm_step_local, only: dim_cs1
  use jules_sea_seaice_mod, only: nice_use
  use jules_soil_mod, only: ns_deep, l_bedrock
  use jules_snow_mod, only: nsmax, rho_snow_const
  use jules_soil_biogeochem_mod, only: dim_ch4layer, soil_bgc_model, &
       soil_model_ecosse, l_layeredc
  use jules_vegetation_mod, only: l_triffid, l_phenol, l_acclim, &
       l_sugar, l_use_pft_psi, l_red
  use jules_surface_mod, only: l_urban2t
  use jules_urban_mod, only: l_moruses, l_moruses_emissivity
  use nlsizes_namelist_mod, only: sm_levels
  use nvegparm, only: emis_nvg
  use pftparm,  only: emis_pft
  use sparm_mod, only: sparm
  use tilepts_mod, only: tilepts

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers, n_band, seg_len
  integer(i_def), intent(in) :: ndf_lw_tile, undf_lw_tile
  integer(i_def), intent(in) :: map_lw_tile(ndf_lw_tile, seg_len)
  integer(i_def), intent(in) :: ndf_tile, undf_tile
  integer(i_def), intent(in) :: map_tile(ndf_tile, seg_len)
  integer(i_def), intent(in) :: ndf_pft, undf_pft
  integer(i_def), intent(in) :: map_pft(ndf_pft, seg_len)
  integer(i_def), intent(in) :: ndf_2d, undf_2d
  integer(i_def), intent(in) :: map_2d(ndf_2d, seg_len)

  real(r_def), intent(inout) :: tile_lw_albedo(undf_lw_tile)

  real(r_def), intent(inout) :: tile_lw_grey_albedo(undf_tile)

  real(r_def), intent(in) :: tile_fraction(undf_tile)

  real(r_def), intent(in) :: tile_temperature(undf_tile)
  real(r_def), intent(in) :: tile_snow_mass(undf_tile)

  real(kind=r_def), intent(in) :: leaf_area_index(undf_pft)
  real(kind=r_def), intent(in) :: canopy_height(undf_pft)

  real(kind=r_def), intent(in) :: soil_roughness(undf_2d)
  real(kind=r_def), intent(in) :: urbztm(undf_2d)
  real(kind=r_def), intent(in) :: urbemisc(undf_2d)

  ! Local variables for the kernel
  integer(i_def) :: i_tile, i_band, k, l, n, i
  integer(i_def) :: df_rtile

  real(r_def) :: specemis_sea(n_band, seg_len)
  real(r_def) :: specemis_soil(n_band, seg_len)

  real(r_def) :: greyemis_sea
  real(r_def) :: greyemis_soil

  logical(l_def) :: update_required
  real(r_def) :: flandg(seg_len)
  real(r_def) :: snow_fac(seg_len, n_land_tile)
  integer(i_def) :: land_field
  type(ainfo_type) :: ainfo
  type(ainfo_data_type) :: ainfo_data
  type(progs_type) :: progs
  type(progs_data_type) :: progs_data
  type(psparms_type) :: psparms
  type(psparms_data_type) :: psparms_data
  type(urban_param_type) :: urban_param
  type(urban_param_data_type) :: urban_param_data

  ! Include effect of snow on emissivity
  ! First check if this is required
  update_required = .false.
  iloop: do i = 1, seg_len
    nloop: do n = 1, n_land_tile
      if (tile_snow_mass(map_tile(1,i)+n-1) > 0.0_r_def) then
        update_required = .true.
        exit iloop
      end if
    end do nloop
  end do iloop

  do n = 1, n_land_tile
    do i = 1, seg_len
      snow_fac(i,n) = 0.0_r_def
    end do
  end do
  if (update_required) then

    ! Calculate Jules information needed
    land_field = 0
    do i = 1, seg_len
      flandg(i) = 0.0_r_def
      do n = 1, n_land_tile
        flandg(i) = flandg(i) + tile_fraction(map_tile(1,i)+n-1)
      end do
      if (flandg(i) > 0.0_r_def) then
        land_field = land_field + 1
      end if
    end do

    call ancil_info_alloc(land_field, seg_len, 1, nice_use, nsoilt, ntype,     &
                          ainfo_data)
    call ancil_info_assoc(ainfo, ainfo_data)
    call prognostics_alloc(land_field, seg_len, 1, n_land_tile, npft, nsoilt,  &
                           sm_levels, ns_deep, nsmax, dim_cslayer, dim_cs1,    &
                           dim_ch4layer, nice_use, nice_use, soil_bgc_model,   &
                           soil_model_ecosse, l_layeredc, l_triffid, l_phenol, &
                           l_bedrock, l_red, nmasst, nnpft, l_acclim,          &
                           l_sugar, progs_data)
    call prognostics_assoc(progs,progs_data)
    call psparms_alloc(land_field, seg_len, 1, nsoilt, sm_levels, dim_cslayer, &
                       n_land_tile, npft, soil_bgc_model, soil_model_ecosse,   &
                       l_use_pft_psi, psparms_data)
    call psparms_assoc(psparms, psparms_data)
    call urban_param_alloc(land_field, l_urban2t, l_moruses, urban_param_data)
    call urban_param_assoc(urban_param, urban_param_data)

    l = 0
    do i = 1, seg_len
      if (flandg(i) > 0.0_r_def) then
        l = l+1
        ainfo%land_index(l) = i
      end if
    end do

    do l = 1, land_field
      i = ainfo%land_index(l)
      do n = 1, n_land_tile
        ainfo%frac_surft(l,n) = tile_fraction(map_tile(1,i)+n-1) / flandg(i)
      end do
    end do

    ! Set type_pts and type_index
    call tilepts(land_field, ainfo%frac_surft, ainfo%surft_pts,                &
                 ainfo%surft_index, ainfo%l_lice_point, ainfo%l_lice_surft)

    do l = 1, land_field
      i = ainfo%land_index(l)
      do n = 1, n_land_tile
        progs%snow_surft(l,n) = tile_snow_mass(map_tile(1,i)+n-1)
      end do
      do n = 1, npft
        ! Leaf area index
        progs%lai_pft(l,n) = leaf_area_index(map_pft(1,i)+n-1)
        ! Canopy height
        progs%canht_pft(l,n) = canopy_height(map_pft(1,i)+n-1)
      end do
      ! Roughness length (z0_tile)
      psparms%z0m_soil_gb(l) = soil_roughness(map_2d(1,i))
    end do
    ! Urban ancillaries
    if ( l_urban2t ) then
      do l = 1, land_field
        i = ainfo%land_index(l)
        urban_param%ztm_gb(l) = urbztm(map_2d(1,i))
      end do
    end if

    call sparm(land_field, n_land_tile, ainfo%surft_pts, ainfo%surft_index,   &
               ainfo%frac_surft, progs%canht_pft, progs%lai_pft,              &
               psparms%z0m_soil_gb, psparms%catch_snow_surft,                 &
               psparms%catch_surft, psparms%z0_surft, psparms%z0h_bare_surft, &
               urban_param%ztm_gb)

    do n = 1, n_land_tile
      do k = 1, ainfo%surft_pts(n)
        l = ainfo%surft_index(k,n)
        i = ainfo%land_index(l)
        snow_fac(i,n) = ( max(0.0_r_def, progs%snow_surft(l,n)) /              &
                        ( max(0.0_r_def, progs%snow_surft(l,n)) +              &
                        10.0_r_def * psparms%z0_surft(l,n) * rho_snow_const ) )
      end do
    end do

    call ancil_info_nullify(ainfo)
    call ancil_info_dealloc(ainfo_data)
    call urban_param_nullify(urban_param)
    call urban_param_dealloc(urban_param_data)
    call psparms_nullify(psparms)
    call psparms_dealloc(psparms_data)
    call prognostics_nullify(progs)
    call prognostics_dealloc(progs_data)

  end if ! snow update

  ! Recompute grey albedos for land tiles (snow could have melted)
  do n = 1, npft
    do i = 1, seg_len
      tile_lw_grey_albedo(map_tile(1,i)+n-1) = 1.0_r_def - &
           (emis_pft(n)+(emis_nvg(ice-npft)-emis_pft(n))*snow_fac(i,n))
    end do
  end do
  do n = npft+1, n_land_tile
    if (l_moruses_emissivity .and. n == urban_canyon) then
      do i = 1, seg_len
        tile_lw_grey_albedo(map_tile(1,i)+n-1) = 1.0_r_def - &
             urbemisc(map_2d(1,i))
      end do
    else
      do i = 1, seg_len
        tile_lw_grey_albedo(map_tile(1,i)+n-1) = 1.0_r_def - &
          (emis_nvg(n-npft)+(emis_nvg(ice-npft)-emis_nvg(n-npft))*snow_fac(i,n))
      end do
    end if
  end do

  ! If spectrally varying albedos/emissivities are used for desert tiles
  ! call specemis and set tile_lw_grey_albedo and tile_lw_albedo in the
  ! later part of this routine accordingly
  if (emis_method_soil == emis_method_soil_feldman_desert) then
    do i = 1, seg_len
      if ( (tile_fraction(map_tile(1,i)+soil-1) > 0.0_r_def)) then

        call specemis('desert_feldman', n_band, wavelength_short,              &
                      wavelength_long, tile_temperature(map_tile(1,i)+soil-1), &
                      n_band_exclude, index_exclude, specemis_soil(:,i),       &
                      greyemis_soil)
        tile_lw_grey_albedo(map_tile(1,i)+soil-1) = 1.0_r_def - &
             (greyemis_soil+(emis_nvg(ice-npft)-greyemis_soil)*snow_fac(i,soil))
      end if
    end do
  end if

  ! Perform the same steps if spectrally varying albedos/emissivities
  ! are used for sea tiles (from feldman or iremis)
  if (emis_method_sea == emis_method_sea_feldman .or. &
       emis_method_sea == emis_method_sea_iremis) then
    do i = 1, seg_len
      if (tile_fraction(map_tile(1,i)+first_sea_tile-1) > 0.0_r_def) then

        if (emis_method_sea == emis_method_sea_feldman) then
          call specemis('sea_feldman', n_band, wavelength_short,               &
                        wavelength_long,                                       &
                        tile_temperature(map_tile(1,i)+first_sea_tile-1),      &
                        n_band_exclude, index_exclude, specemis_sea(:,i),      &
                        greyemis_sea)
        else
          call specemis('sea_iremis', n_band, wavelength_short,                &
                        wavelength_long,                                       &
                        tile_temperature(map_tile(1,i)+first_sea_tile-1),      &
                        n_band_exclude, index_exclude, specemis_sea(:,i),      &
                        greyemis_sea)
        end if
        tile_lw_grey_albedo(map_tile(1,i)+first_sea_tile-1) &
             = 1.0_r_def - greyemis_sea
      end if
    end do
  end if

  ! Now set the tile_lw_albedo for land, sea and sea-ice tiles
  ! If constant albedos/emissivities are used these can are copied from tile_lw_grey_albedo
  ! If spectrally varying emissivities are used these have been calculated
  ! per band in the calls to specemis above and are applied below to set tile_lw_albedo
  do i = 1, seg_len
    do i_band = 1, n_band

      ! Land tile albedos
      df_rtile = n_surf_tile * (i_band-1)
      do i_tile = 1, n_land_tile
        df_rtile = df_rtile + 1
        if ( (tile_fraction(map_tile(1,i)+i_tile-1) > 0.0_r_def) .and. &
             ((i_tile == soil) .and. &
             (emis_method_soil /= emis_method_soil_fixed)) ) then
          tile_lw_albedo(map_lw_tile(1,i)+df_rtile-1) = 1.0_r_def - &
               ( specemis_soil(i_band,i) + &
              (emis_nvg(ice-npft)-specemis_soil(i_band,i)) * snow_fac(i,i_tile))
        else
          tile_lw_albedo(map_lw_tile(1,i)+df_rtile-1) &
               = tile_lw_grey_albedo(map_tile(1,i)+i_tile-1)
        end if
      end do

      ! Sea tile albedos
      df_rtile = first_sea_tile-1 + n_surf_tile * (i_band-1)
      do i_tile = first_sea_tile, first_sea_tile + n_sea_tile - 1
        df_rtile = df_rtile + 1
        if ( (tile_fraction(map_tile(1,i)+i_tile-1) > 0.0_r_def).and. &
             (emis_method_sea /= emis_method_sea_fixed)) then
          tile_lw_albedo(map_lw_tile(1,i)+df_rtile-1) &
               = 1.0_r_def - specemis_sea(i_band,i)
        else
          tile_lw_albedo(map_lw_tile(1,i)+df_rtile-1) &
               = tile_lw_grey_albedo(map_tile(1,i)+i_tile-1)
        end if
      end do

      ! Sea-ice tile albedos
      df_rtile = first_sea_ice_tile - 1 + n_surf_tile * (i_band-1)
      do i_tile = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
        df_rtile = df_rtile + 1
        tile_lw_albedo(map_lw_tile(1,i)+df_rtile-1) &
             = tile_lw_grey_albedo(map_tile(1,i)+i_tile-1)
      end do

    end do ! bands
  end do ! seg_len

end subroutine lw_rad_tile_code

end module lw_rad_tile_kernel_mod
