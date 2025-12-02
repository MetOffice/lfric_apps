!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Interface to Jules for SW surface tile radiative properties
module sw_rad_tile_kernel_mod

use argument_mod,      only : arg_type,                  &
                              GH_FIELD, GH_SCALAR,       &
                              GH_REAL, GH_INTEGER,       &
                              GH_READ, GH_WRITE,         &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              ANY_DISCONTINUOUS_SPACE_2, &
                              ANY_DISCONTINUOUS_SPACE_3, &
                              ANY_DISCONTINUOUS_SPACE_4, &
                              ANY_DISCONTINUOUS_SPACE_5, &
                              ANY_DISCONTINUOUS_SPACE_6, &
                              DOMAIN
use fs_continuity_mod, only:  W3, WTheta
use constants_mod,     only : r_def, i_def, l_def, r_um, i_um
use kernel_mod,        only : kernel_type

implicit none

private

public :: sw_rad_tile_kernel_type
public :: sw_rad_tile_code
CHARACTER(LEN=*), PARAMETER, PRIVATE :: ModuleName='SW_RAD_TILE_KERNEL_MOD'

!------------------------------------------------------------------------------
! Public types
!------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, extends(kernel_type) :: sw_rad_tile_kernel_type
  private
  type(arg_type) :: meta_args(33) = (/                                &
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! tile_sw_direct_albedo
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! tile_sw_diffuse_albedo
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! tile_fraction
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! leaf_area_index
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), & ! canopy_height
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! sd_orog_2d
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! soil_albedo
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! soil_roughness
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! albedo_obs_vis
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! albedo_obs_nir
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_5), & ! albedo_obs_scaling
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! tile_temperature
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! tile_snow_mass
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! tile_snow_rgrain
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! snow_depth
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! snowpack_density
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! snow_soot
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! chloro_sea
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_6), & ! sea_ice_thickness
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_6), & ! melt_pond_fraction
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_6), & ! melt_pond_depth
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_6), & ! sea_ice_pensolar_frac_direct
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_6), & ! sea_ice_pensolar_frac_diffuse
    arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        & ! u_in_w3
    arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        & ! v_in_w3
    arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! dz_wth
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! z0msea
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_4), & ! cos_zenith_angle_rts
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! urbhwr
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! urbztm
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! urbalbwl
    arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! urbalbrd
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                          & ! n_band
    /)
  integer :: operates_on = DOMAIN
contains
  procedure, nopass :: sw_rad_tile_code
end type

!------------------------------------------------------------------------------
! Contained functions/subroutines
!------------------------------------------------------------------------------
contains

!> @param[in]     nlayers                Number of layers
!> @param[in]     seg_len                Number of horizontal cells in segment
!> @param[in,out] tile_sw_direct_albedo  SW direct tile albedos
!> @param[in,out] tile_sw_diffuse_albedo SW diffuse tile albedos
!> @param[in]     tile_fraction          Surface tile fractions
!> @param[in]     leaf_area_index        Leaf Area Index
!> @param[in]     canopy_height          Canopy height
!> @param[in]     sd_orog_2d             Standard deviation of orography
!> @param[in]     soil_albedo            Snow-free soil albedo
!> @param[in]     soil_roughness         Bare soil surface roughness length
!> @param[in]     albedo_obs_vis         Observed snow-free visible albedo
!> @param[in]     albedo_obs_nir         Observed snow-free near-IR albedo
!> @param[in,out] albedo_obs_scaling     Scaling factor to adjust albedos by
!> @param[in]     tile_temperature       Surface tile temperatures
!> @param[in]     tile_snow_mass         Snow mass on tiles (kg/m2)
!> @param[in]     tile_snow_rgrain       Snow grain size on tiles (microns)
!> @param[in]     snow_depth             Snow depth on tiles (m)
!> @param[in]     snowpack_density       Density of snow on ground (kg m-3)
!> @param[in]     snow_soot              Snow soot content (kg/kg)
!> @param[in]     chloro_sea             Chlorophyll content of the sea
!> @param[in]     sea_ice_thickness      Sea ice thickness (m)
!> @param[in]     melt_pond_fraction     Sea ice melt pond fraction
!> @param[in]     melt_pond_depth        Sea ice melt pond depth (m)
!> @param[in,out] sea_ice_pensolar_frac_direct  Fraction of penetrating solar (direct visible) into sea ice
!> @param[in,out] sea_ice_pensolar_frac_diffuse  Fraction of penetrating solar (diffuse visible) into sea ice
!> @param[in]     u_in_w3                'Zonal' wind in density space
!> @param[in]     v_in_w3                'Meridional' wind in density space
!> @param[in]     dz_wth                 Delta z at wtheta levels
!> @param[in]     z0msea                 Roughness length of sea
!> @param[in]     cos_zenith_angle_rts   Cosine of the stellar zenith angle
!> @param[in]     urbhwr                 Urban height-to-width ratio
!> @param[in]     urbztm                 Urban effective roughness length
!> @param[in]     urbalbwl               Urban wall albedo
!> @param[in]     urbalbrd               Urban road albedo
!> @param[in]     n_band                 Number of spectral bands
!> @param[in]     ndf_sw_tile            DOFs per cell for tiles and sw bands
!> @param[in]     undf_sw_tile           Total DOFs for tiles and sw bands
!> @param[in]     map_sw_tile            Dofmap for a segment of sw_tiles
!> @param[in]     ndf_tile               Number of DOFs per cell for tiles
!> @param[in]     undf_tile              Number of total DOFs for tiles
!> @param[in]     map_tile               Dofmap for a segment of tiles
!> @param[in]     ndf_pft                Number of DOFs per cell for PFTs
!> @param[in]     undf_pft               Number of total DOFs for PFTs
!> @param[in]     map_pft                Dofmap for a segment of PFTs
!> @param[in]     ndf_2d                 Number of DOFs per cell for 2D fields
!> @param[in]     undf_2d                Number of total DOFs for 2D fields
!> @param[in]     map_2d                 Dofmap for a segment of 2D fields
!> @param[in]     ndf_scal               Number of DOFs per cell for albedo scaling
!> @param[in]     undf_scal              Number of total DOFs for albedo scaling
!> @param[in]     map_scal               Dofmap for a segment of scaling
!> @param[in]     ndf_sice               Number of DOFs per cell for sea ice tiles
!> @param[in]     undf_sice              Number of total DOFs for sea ice tiles
!> @param[in]     map_sice               Dofmap for a segment of sea ice tiles
!> @param[in]     ndf_w3                 Number of DOFs per cell for density space
!> @param[in]     undf_w3                Number of unique DOFs for density space
!> @param[in]     map_w3                 Dofmap for a segment of density space
!> @param[in]     ndf_wth                Number of DOFs per cell for theta space
!> @param[in]     undf_wth               Number of unique DOFs for theta space
!> @param[in]     map_wth                Dofmap for a segment of potential temperatures
subroutine sw_rad_tile_code(nlayers, seg_len,                       &
                            tile_sw_direct_albedo,                  &
                            tile_sw_diffuse_albedo,                 &
                            tile_fraction,                          &
                            leaf_area_index,                        &
                            canopy_height,                          &
                            sd_orog_2d,                             &
                            soil_albedo,                            &
                            soil_roughness,                         &
                            albedo_obs_vis,                         &
                            albedo_obs_nir,                         &
                            albedo_obs_scaling,                     &
                            tile_temperature,                       &
                            tile_snow_mass,                         &
                            tile_snow_rgrain,                       &
                            snow_depth,                             &
                            snowpack_density,                       &
                            snow_soot,                              &
                            chloro_sea,                             &
                            sea_ice_thickness,                      &
                            melt_pond_fraction,                     &
                            melt_pond_depth,                        &
                            sea_ice_pensolar_frac_direct,           &
                            sea_ice_pensolar_frac_diffuse,          &
                            u_in_w3,                                &
                            v_in_w3,                                &
                            dz_wth,                                 &
                            z0msea,                                 &
                            cos_zenith_angle_rts,                   &
                            urbhwr,                                 &
                            urbztm,                                 &
                            urbalbwl,                               &
                            urbalbrd,                               &
                            n_band,                                 &
                            ndf_sw_tile, undf_sw_tile, map_sw_tile, &
                            ndf_tile, undf_tile, map_tile,          &
                            ndf_pft, undf_pft, map_pft,             &
                            ndf_2d, undf_2d, map_2d,                &
                            ndf_scal, undf_scal, map_scal,          &
                            ndf_sice, undf_sice, map_sice,          &
                            ndf_w3, undf_w3, map_w3,                &
                            ndf_wth, undf_wth, map_wth)

  !-----------------------------------------------------------------------------
  ! LFRic modules
  !-----------------------------------------------------------------------------
  use jules_control_init_mod,    only: n_surf_tile, n_land_tile, n_sea_tile,   &
                                       n_sea_ice_tile, first_sea_tile,         &
                                       first_sea_ice_tile

  !-----------------------------------------------------------------------------
  ! UM/Jules modules containing switches or global constants
  !-----------------------------------------------------------------------------
  use socrates_init_mod,         only: wavelength_short, wavelength_long,      &
                                       weight_blue
  use atm_step_local,            only: dim_cs1
  use atm_fields_bounds_mod,     only: pdims_s, pdims
  use nlsizes_namelist_mod,      only: ntiles, sm_levels, bl_levels
  use jules_surface_types_mod,   only: ntype, npft, ice, nnpft
  use jules_sea_seaice_mod,      only: nice, nice_use
  use ancil_info,                only: rad_nband, dim_cslayer, nsoilt, nmasst
  use jules_vegetation_mod,      only: l_triffid, l_phenol, l_use_pft_psi,     &
                                       can_rad_mod, l_acclim, l_sugar, l_red
  use jules_soil_mod,            only: ns_deep, l_bedrock
  use jules_soil_biogeochem_mod, only: dim_ch4layer, soil_bgc_model,           &
                                       soil_model_ecosse, l_layeredc
  use jules_radiation_mod,       only: l_albedo_obs
  use jules_snow_mod,            only: nsmax, cansnowtile
  use jules_deposition_mod,      only: l_deposition
  use jules_surface_mod,         only: l_urban2t, l_flake_model
  use jules_urban_mod,           only: l_moruses

  !-----------------------------------------------------------------------------
  ! JULES modules
  !-----------------------------------------------------------------------------
  use prognostics,               only: progs_data_type, progs_type,            &
                                       prognostics_alloc, prognostics_assoc,   &
                                       prognostics_nullify, prognostics_dealloc
  use jules_vars_mod,            only: jules_vars_type, jules_vars_data_type,  &
                                       jules_vars_alloc, jules_vars_assoc,     &
                                       jules_vars_nullify, jules_vars_dealloc
  use p_s_parms,                 only: psparms_type, psparms_data_type,        &
                                       psparms_alloc, psparms_assoc,           &
                                       psparms_nullify, psparms_dealloc
  use ancil_info,                only: ainfo_data_type, ainfo_type,            &
                                       ancil_info_alloc, ancil_info_assoc,     &
                                       ancil_info_nullify, ancil_info_dealloc
  use coastal,                   only: coastal_type, coastal_data_type,        &
                                       coastal_assoc, coastal_alloc,           &
                                       coastal_nullify, coastal_dealloc
  use urban_param_mod,           only: urban_param_data_type, urban_param_type, &
                                       urban_param_alloc, urban_param_assoc,   &
                                       urban_param_nullify, urban_param_dealloc
  use lake_mod,                  only: lake_data_type, lake_type,              &
                                       lake_assoc, lake_alloc,                 &
                                       lake_nullify, lake_dealloc
  use fluxes_mod,                only: fluxes_type, fluxes_data_type,          &
                                       fluxes_alloc, fluxes_assoc,             &
                                       fluxes_nullify, fluxes_dealloc
  use cable_fields_mod,          only: progs_cbl_vars

  use tilepts_mod,               only: tilepts
  use sparm_mod,                 only: sparm
  use surf_couple_radiation_mod, only: surf_couple_radiation

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers, seg_len, n_band
  integer(i_def), intent(in) :: ndf_sw_tile, undf_sw_tile
  integer(i_def), intent(in) :: map_sw_tile(ndf_sw_tile, seg_len)
  integer(i_def), intent(in) :: ndf_tile, undf_tile
  integer(i_def), intent(in) :: map_tile(ndf_tile, seg_len)
  integer(i_def), intent(in) :: ndf_pft, undf_pft
  integer(i_def), intent(in) :: map_pft(ndf_pft, seg_len)
  integer(i_def), intent(in) :: ndf_2d, undf_2d
  integer(i_def), intent(in) :: map_2d(ndf_2d, seg_len)
  integer(i_def), intent(in) :: ndf_scal, undf_scal
  integer(i_def), intent(in) :: map_scal(ndf_scal, seg_len)
  integer(i_def), intent(in) :: ndf_sice, undf_sice
  integer(i_def), intent(in) :: map_sice(ndf_sice, seg_len)
  integer(i_def), intent(in) :: ndf_w3, undf_w3
  integer(i_def), intent(in) :: map_w3(ndf_w3, seg_len)
  integer(i_def), intent(in) :: ndf_wth, undf_wth
  integer(i_def), intent(in) :: map_wth(ndf_wth, seg_len)
  real(r_def), intent(inout) :: tile_sw_direct_albedo(undf_sw_tile)
  real(r_def), intent(inout) :: tile_sw_diffuse_albedo(undf_sw_tile)

  real(r_def), intent(in)    :: tile_fraction(undf_tile)
  real(r_def), intent(in)    :: tile_temperature(undf_tile)
  real(r_def), intent(in)    :: tile_snow_mass(undf_tile)
  real(r_def), intent(in)    :: tile_snow_rgrain(undf_tile)
  real(r_def), intent(in)    :: snow_depth(undf_tile)
  real(r_def), intent(in)    :: snowpack_density(undf_tile)

  real(r_def), intent(in)    :: leaf_area_index(undf_pft)
  real(r_def), intent(in)    :: canopy_height(undf_pft)

  real(r_def), intent(in)    :: sd_orog_2d(undf_2d)
  real(r_def), intent(in)    :: soil_albedo(undf_2d)
  real(r_def), intent(in)    :: soil_roughness(undf_2d)
  real(r_def), intent(in)    :: albedo_obs_vis(undf_2d)
  real(r_def), intent(in)    :: albedo_obs_nir(undf_2d)
  real(r_def), intent(inout) :: albedo_obs_scaling(undf_scal)
  real(r_def), intent(in)    :: snow_soot(undf_2d)
  real(r_def), intent(in)    :: chloro_sea(undf_2d)
  real(r_def), intent(in)    :: z0msea(undf_2d)
  real(r_def), intent(in)    :: cos_zenith_angle_rts(undf_2d)
  real(r_def), intent(in)    :: urbhwr(undf_2d)
  real(r_def), intent(in)    :: urbztm(undf_2d)
  real(r_def), intent(in)    :: urbalbwl(undf_2d)
  real(r_def), intent(in)    :: urbalbrd(undf_2d)

  real(r_def), intent(in)    :: sea_ice_thickness(undf_sice)
  real(r_def), intent(in)    :: melt_pond_fraction(undf_sice)
  real(r_def), intent(in)    :: melt_pond_depth(undf_sice)
  real(r_def), intent(inout) :: sea_ice_pensolar_frac_direct(undf_sice)
  real(r_def), intent(inout) :: sea_ice_pensolar_frac_diffuse(undf_sice)

  real(r_def), intent(in)    :: u_in_w3(undf_w3)
  real(r_def), intent(in)    :: v_in_w3(undf_w3)
  real(r_def), intent(in)    :: dz_wth(undf_wth)

  !-----------------------------------------------------------------------
  ! Local variables for the kernel
  !-----------------------------------------------------------------------
  integer(i_def) :: i, i_sice, i_band, l, n, r1, r2, r3
  integer(i_def) :: df_rtile, land_field, sea_pts

  ! fields on land points
  real(r_um), dimension(:), allocatable          :: fland
  real(r_um), dimension(:,:), allocatable        :: snow_surft

  ! single level logical fields
  logical(l_def), dimension(seg_len,1)                  :: land_sea_mask

  ! Reverse index of land points
  integer(i_def), dimension(seg_len,1)                  :: reverse_land_index

  ! Inputs to surf_couple_radiation
  real(r_um), dimension(seg_len, 1)              :: flandg, ws_10m_sea, chloro
  integer, dimension(ntype)                      :: type_pts

  ! Outputs from surf_couple_radiation
  real(r_um), dimension(seg_len, 1, 4)           :: sea_ice_albedo
  real(r_um), dimension(seg_len, 1, ntiles, 2)   :: albobs_sc
  real(r_um), dimension(seg_len, 1,  2, n_band)  :: open_sea_albedo

  !-----------------------------------------------------------------------
  ! JULES Types
  !-----------------------------------------------------------------------
  type(progs_type)            :: progs
  type(progs_data_type)       :: progs_data
  type(jules_vars_type)       :: jules_vars
  type(jules_vars_data_type)  :: jules_vars_data
  type(psparms_type)          :: psparms
  type(psparms_data_type)     :: psparms_data
  type(ainfo_type)            :: ainfo
  type(ainfo_data_type)       :: ainfo_data
  type(urban_param_type)      :: urban_param
  type(urban_param_data_type) :: urban_param_data
  type(coastal_type)          :: coast
  type(coastal_data_type)     :: coastal_data
  type(lake_type)             :: lake_vars
  type(lake_data_type)        :: lake_data
  type(fluxes_type)           :: fluxes
  type(fluxes_data_type)      :: fluxes_data

  !-----------------------------------------------------------------------
  ! Initialisation of JULES data and pointer types
  !-----------------------------------------------------------------------

  ! Land tile fractions
  land_field = 0
  do i = 1, seg_len
    flandg(i,1) = 0.0_r_um
    do n = 1, n_land_tile
      flandg(i,1) = flandg(i,1) + real(tile_fraction(map_tile(1,i)+n-1), r_um)
    end do
    if ( flandg(i,1) > 0.0_r_um ) then
      land_field = land_field + 1
      land_sea_mask(i,1) = .true.
    else
      land_sea_mask(i,1) = .false.
    end if
  end do ! i

  call prognostics_alloc(land_field, seg_len, 1,                              &
                      n_land_tile, npft, nsoilt, sm_levels, ns_deep, nsmax,   &
                      dim_cslayer, dim_cs1, dim_ch4layer,                     &
                      nice, nice_use, soil_bgc_model, soil_model_ecosse,      &
                      l_layeredc, l_triffid, l_phenol, l_bedrock, l_red,      &
                      nmasst, nnpft, l_acclim, l_sugar, progs_data)
  call prognostics_assoc(progs,progs_data)

  call jules_vars_alloc(land_field,ntype,n_land_tile,rad_nband,nsoilt,        &
                   sm_levels, seg_len, 1, npft, bl_levels, pdims_s, pdims,    &
                   l_albedo_obs, cansnowtile, l_deposition, jules_vars_data)
  call jules_vars_assoc(jules_vars,jules_vars_data)

  if ( can_rad_mod == 6 ) then
    jules_vars%diff_frac = 0.4_r_um
  end if

  call psparms_alloc(land_field, seg_len, 1, nsoilt,sm_levels,dim_cslayer,    &
                   n_land_tile,npft, soil_bgc_model,soil_model_ecosse,        &
                   l_use_pft_psi, psparms_data)
  call psparms_assoc(psparms, psparms_data)

  call ancil_info_alloc(land_field, seg_len, 1, nice,nsoilt,ntype, ainfo_data)
  call ancil_info_assoc(ainfo, ainfo_data)

  call urban_param_alloc(land_field, l_urban2t, l_moruses, urban_param_data)
  call urban_param_assoc(urban_param, urban_param_data)

  call coastal_alloc(land_field, seg_len, 1, seg_len, 1, seg_len, 1,          &
                   nice_use,nice,coastal_data)
  call coastal_assoc(coast, coastal_data)

  call lake_alloc(land_field, l_flake_model, lake_data)
  call lake_assoc(lake_vars, lake_data)

  call fluxes_alloc(land_field, seg_len, 1, n_land_tile, npft, nsoilt,        &
                   sm_levels, nice, nice_use, fluxes_data)
  call fluxes_assoc(fluxes, fluxes_data)

  ! ---------------------------------------------------------------------------
  ! SW tile albedos
  ! ---------------------------------------------------------------------------

  allocate(snow_surft(land_field, ntiles))

  ! Land_index
  reverse_land_index = 0
  l = 0
  do i = 1, seg_len
    if ( flandg(i,1) > 0.0_r_um ) then
      l = l+1
      ainfo%land_index(l) = i
      reverse_land_index(i,1) = l
    end if
  end do ! i

  ! Sea-ice fraction
  do i = 1, seg_len
    ainfo%ice_fract_ij(i,1) = 0.0_r_um
    do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = n - first_sea_ice_tile + 1
      ainfo%ice_fract_ij(i,1) = ainfo%ice_fract_ij(i,1) + &
                           real(tile_fraction(map_tile(1,i)+n-1), r_um)
      ainfo%ice_fract_ncat_sicat(i,1, i_sice) = &
                           real(tile_fraction(map_tile(1,i)+n-1), r_um)
    end do ! n

    ! Because Jules tests on flandg < 1, we need to ensure this is exactly
    ! 1 when no sea or sea-ice is present
    if ( tile_fraction(map_tile(1,i)+first_sea_tile-1) == 0.0_r_def .and. &
         ainfo%ice_fract_ij(i,1) == 0.0_r_um ) then
      flandg(i,1) = 1.0_r_um
    end if

    ! Jules requires sea-ice fractions with respect to the sea area
    if ( ainfo%ice_fract_ij(i,1) > 0.0_r_um ) then
      ainfo%ice_fract_ij(i,1) = ainfo%ice_fract_ij(i,1) / (1.0_r_um - flandg(i,1))
      ainfo%ice_fract_ncat_sicat(i,1,1:n_sea_ice_tile) = &
           ainfo%ice_fract_ncat_sicat(i,1,1:n_sea_ice_tile) / (1.0_r_um - flandg(i,1))
    end if

  end do ! i

  ! Allocate land_field
  allocate(fland(land_field))
  do l = 1, land_field
    fland(l) = flandg(ainfo%land_index(l),1)
    do n = 1, n_land_tile
      ! Jules requires fractions with respect to the land area
      ainfo%frac_surft(l, n) = &
          real(tile_fraction(map_tile(1,ainfo%land_index(l))+n-1), r_um) / fland(l)
    end do ! n
    if ( ainfo%frac_surft(l,ice) > 0.0_r_def ) then
      ainfo%l_lice_point(l) = .true.
    end if
  end do ! l

  ! Set type_pts and type_index
  call tilepts(land_field, ainfo%frac_surft, type_pts, ainfo%surft_index, &
                ainfo%l_lice_point, ainfo%l_lice_surft)

  ! Land tile temperatures
  do l = 1, land_field
    do n = 1, n_land_tile
      progs%tstar_surft(l,n) = &
            real(tile_temperature(map_tile(1,ainfo%land_index(l))+n-1), r_um)
    end do ! n
  end do ! l

  ! Sea temperature
  do i = 1, seg_len
    fluxes%tstar_ij(i,1) = &
            real(tile_temperature(map_tile(1,i)+first_sea_tile-1), r_um)
  end do ! i

  ! Sea-ice temperatures
  do i = 1, seg_len
    if ( ainfo%ice_fract_ij(i,1) > 0.0_r_um ) then
      do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
        i_sice = n - first_sea_ice_tile + 1
        coast%tstar_sice_sicat(i,1,i_sice) = &
            real(tile_temperature(map_tile(1,i)+n-1), r_um)
      end do ! n
    end if
  end do ! i

  ! Combined sea and sea-ice index
  do i = 1, seg_len
    if ( flandg(i,1) < 1.0_r_um ) then
      ainfo%ssi_index(i) = i
    end if
  end do ! i

  ! Individual sea and sea-ice indices
  ! first set defaults
  sea_pts = 0
  ! Then adjust based on state
  do i = 1, seg_len
    if ( ainfo%ssi_index(i) > 0 ) then
      if ( ainfo%ice_fract_ij(i,1) < 1.0_r_um ) then
        sea_pts = sea_pts + 1
        ainfo%sea_index(sea_pts) = i
      end if
    end if
  end do ! i

  ! Multi-category sea-ice index
  do i = 1, seg_len
    do n = 1, nice_use
      if ( ainfo%ssi_index(i) > 0 .and. ainfo%ice_fract_ncat_sicat(i,1,n) > 0.0_r_um ) then
        ainfo%sice_pts_ncat(n) = ainfo%sice_pts_ncat(n)+1
        ainfo%sice_index_ncat(ainfo%sice_pts_ncat(n),n) = i
        ainfo%sice_frac_ncat(i,n) = ainfo%ice_fract_ncat_sicat(i,1,n)
      end if
    end do ! n
  end do ! i

  ! Other sea ice fields
  do i = 1, seg_len
    do n = 1, n_sea_ice_tile
      ! Sea-ice thickness
      progs%di_ncat_sicat(i,1,n) = real(sea_ice_thickness(map_sice(1,i)+n-1),r_um)
      ! Sea ice melt ponds
      ainfo%pond_frac_cat_sicat(i,1,n) = real(melt_pond_fraction(map_sice(1,i)+n-1),r_um)
      ainfo%pond_depth_cat_sicat(i,1,n) = real(melt_pond_depth(map_sice(1,i)+n-1),r_um)
    end do ! n
  end do ! i

  do l = 1, land_field
    do n = 1, npft
      ! Leaf area index
      progs%lai_pft(l,n) = real(leaf_area_index(map_pft(1,ainfo%land_index(l))+n-1), r_um)
      ! Canopy height
      progs%canht_pft(l,n) = real(canopy_height(map_pft(1,ainfo%land_index(l))+n-1), r_um)
    end do ! n
    ! Roughness length (z0_surft)
    psparms%z0m_soil_gb(l) = real(soil_roughness(map_2d(1,ainfo%land_index(l))), r_um)
  end do ! l

  if ( l_urban2t ) then
    ! Urban fields (including morphology required for surf_couple_radiation).
    ! Allocated with size(1) if not l_urban2t
    do l = 1, land_field
      urban_param%ztm_gb(l)   = real(urbztm(map_2d(1,ainfo%land_index(l))), r_um)
      urban_param%hwr_gb(l)   = real(urbhwr(map_2d(1,ainfo%land_index(l))), r_um)
      urban_param%albwl_gb(l) = real(urbalbwl(map_2d(1,ainfo%land_index(l))), r_um)
      urban_param%albrd_gb(l) = real(urbalbrd(map_2d(1,ainfo%land_index(l))), r_um)
    end do
  end if

  call sparm(land_field, n_land_tile, type_pts, ainfo%surft_index,         &
    ainfo%frac_surft, progs%canht_pft, progs%lai_pft, psparms%z0m_soil_gb, &
    psparms%catch_snow_surft, psparms%catch_surft, psparms%z0_surft,       &
    psparms%z0h_bare_surft, urban_param%ztm_gb)

  ! Snow-free soil albedo
  do l = 1, land_field
    psparms%albsoil_soilt(l,1) = real(soil_albedo(map_2d(1,ainfo%land_index(l))), r_um)
  end do ! l

  ! Cosine of the solar zenith angle
  do i = 1, seg_len
    psparms%cosz_ij(i,1) = real(cos_zenith_angle_rts(map_2d(1,i)), r_um)
  end do ! i

  ! Standard deviation of orography - note that the variables names here
  ! appear to mismatch; this mirrors what is done in the UM; it's possible
  ! that the variable is misnamed in JULES
  do l = 1, land_field
    jules_vars%ho2r2_orog_gb(l) = real(sd_orog_2d(map_2d(1,ainfo%land_index(l))), r_um)
  end do ! l

  do i = 1, seg_len
    ! 10m wind speed over the sea
    ws_10m_sea(i,1) = sqrt(u_in_w3(map_w3(1,i))**2 + v_in_w3(map_w3(1,i))**2) &
         * log(10.0_r_def / z0msea(map_2d(1,i))) &
         / log((dz_wth(map_wth(1,i))) / z0msea(map_2d(1,i)))
    ! Chlorophyll content of the sea
    chloro(i,1) = real(chloro_sea(map_2d(1,i)), r_um)
  end do ! i

  ! Observed albedo
  if (l_albedo_obs) then
    do l = 1, land_field
      psparms%albobs_vis_gb(l) = real(albedo_obs_vis(map_2d(1,ainfo%land_index(l))), r_um)
      psparms%albobs_nir_gb(l) = real(albedo_obs_nir(map_2d(1,ainfo%land_index(l))), r_um)
    end do ! l
  end if

  ! Lying snow mass on land tiles
  do l = 1, land_field
    do n = 1, n_land_tile
      snow_surft(l,n) = real(tile_snow_mass(map_tile(1,ainfo%land_index(l))+n-1), r_um)
      progs%rgrain_surft(l,n) = &
                           real(tile_snow_rgrain(map_tile(1,ainfo%land_index(l))+n-1), r_um)
      progs%snowdepth_surft(l,n) = &
                           real(snow_depth(map_tile(1,ainfo%land_index(l))+n-1), r_um)
      progs%rho_snow_grnd_surft(l,n) = &
                           real(snowpack_density(map_tile(1,ainfo%land_index(l))+n-1), r_um)
    end do ! n
  end do ! l

  do i = 1, seg_len
    ! Lying snow mass on sea ice categories
    do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      i_sice = n - first_sea_ice_tile + 1
      progs%snow_mass_sea_sicat(i,1,i_sice) = &
                          real(tile_snow_mass(map_tile(1,i)+n-1), r_um)
    end do ! n
    ! Snow soot content
    progs%soot_ij(i,1) = real(snow_soot(map_2d(1,i)), r_um)
  end do ! i

  !-----------------------------------------------------------------------
  ! External science code called
  !-----------------------------------------------------------------------

  call surf_couple_radiation( &
    ! Misc INTENT(IN)
    ws_10m_sea, chloro, &
    n_band, n_band, real(wavelength_short,r_um), real(wavelength_long,r_um), &
    ! Misc INTENT(OUT)
    sea_ice_albedo, &
    ! (ancil_info mod)
    ntiles, land_field, sea_pts, type_pts, seg_len, 1, &
    ! (coastal mod)
    flandg, &
    ! (prognostics mod)
    snow_surft, &
    ! UM-only args: INTENT(OUT)
    albobs_sc, open_sea_albedo, &
    ! JULES types
    psparms, ainfo, urban_param, progs, coast, jules_vars, &
    fluxes, lake_vars, &
    progs_cbl_vars  &
    )

  do i_band = 1, n_band
    ! Land tile albedos
    do n = 1, n_land_tile
      df_rtile = n_surf_tile*(i_band-1) + n
      do i = 1, seg_len
        r1 = map_tile(1,i)+n-1
        r2 = map_sw_tile(1,i)+df_rtile-1
        l = reverse_land_index(i,1)
        if ( tile_fraction(r1) > 0.0_r_def ) then
          tile_sw_direct_albedo(r2) &
               = weight_blue(i_band) &
               * real(fluxes%alb_surft(l, n, 1), r_def) & ! visible direct albedo
               + (1.0_r_def - weight_blue(i_band)) &
               * real(fluxes%alb_surft(l, n, 3), r_def)   ! near-ir direct albedo
          tile_sw_diffuse_albedo(r2) &
               = weight_blue(i_band) &
               * real(fluxes%alb_surft(l, n, 2), r_def) & ! visible diffuse albedo
               + (1.0_r_def - weight_blue(i_band)) &
               * real(fluxes%alb_surft(l, n, 4), r_def)   ! near-ir diffuse albedo
        else
          tile_sw_direct_albedo(r2)  = 0.0_r_def
          tile_sw_diffuse_albedo(r2) = 0.0_r_def
        end if
      end do ! i
    end do ! n

    ! Sea tile albedos
    do n = first_sea_tile, first_sea_tile + n_sea_tile - 1
      df_rtile = n_surf_tile*(i_band-1) + n
      do i = 1, seg_len
        r1 = map_tile(1,i)+n-1
        r2 = map_sw_tile(1,i)+df_rtile-1
        if ( tile_fraction(r1) > 0.0_r_def ) then
          tile_sw_direct_albedo(r2)  = real(open_sea_albedo(i,1,1,i_band), r_def)
          tile_sw_diffuse_albedo(r2) = real(open_sea_albedo(i,1,2,i_band), r_def)
        else
          tile_sw_direct_albedo(r2)  = 0.0_r_def
          tile_sw_diffuse_albedo(r2) = 0.0_r_def
        end if
      end do ! i
    end do ! n

    ! Sea-ice tile albedos
    do n = first_sea_ice_tile, first_sea_ice_tile + n_sea_ice_tile - 1
      df_rtile = n_surf_tile*(i_band-1) + n
      i_sice = n - first_sea_ice_tile + 1
      do i = 1, seg_len
        r1 = map_tile(1,i)+n-1
        r2 = map_sw_tile(1,i)+df_rtile-1
        r3 = map_sice(1,i)+i_sice-1
        if ( tile_fraction(r1) > 0.0_r_def ) then
          tile_sw_direct_albedo(r2) &
               = weight_blue(i_band) &
               * real(fluxes%alb_sicat(i,i_sice,1), r_def) &
               + (1.0_r_def - weight_blue(i_band)) &
               * real(fluxes%alb_sicat(i,i_sice,3), r_def)
          tile_sw_diffuse_albedo(r2) &
               = weight_blue(i_band) &
               * real(fluxes%alb_sicat(i,i_sice,2), r_def) &
               + (1.0_r_def - weight_blue(i_band)) &
               * real(fluxes%alb_sicat(i,i_sice,4), r_def)
          sea_ice_pensolar_frac_direct(r3) &
               = real(fluxes%penabs_rad_frac(i,i_sice,1), r_def)
          sea_ice_pensolar_frac_diffuse(r3) &
               = real(fluxes%penabs_rad_frac(i,i_sice,2), r_def)
        else
          tile_sw_direct_albedo(r2)  = 0.0_r_def
          tile_sw_diffuse_albedo(r2) = 0.0_r_def
          sea_ice_pensolar_frac_direct(r3) = 0.0_r_def
          sea_ice_pensolar_frac_diffuse(r3) = 0.0_r_def
        end if
      end do ! i
    end do ! n
  end do ! i_band

  ! Scaling factors needed for use in surface exchange code
  if ( l_albedo_obs ) then
    do i_band = 1, rad_nband
      do n = 1, n_land_tile
         df_rtile=n_land_tile*(i_band-1) + n - 1
         do l = 1, land_field
          if ( fland(l) > 0.0_r_um ) then
             albedo_obs_scaling(map_scal(1,ainfo%land_index(l))+df_rtile) = &
                          jules_vars%albobs_scaling_surft(l,n,i_band)
          end if
        end do ! l
      end do ! n
    end do ! i_band
  end if

  deallocate(fland)
  deallocate(snow_surft)

  call ancil_info_nullify(ainfo)
  call ancil_info_dealloc(ainfo_data)

  call lake_nullify(lake_vars)
  call lake_dealloc(lake_data)

  call coastal_nullify(coast)
  call coastal_dealloc(coastal_data)

  call urban_param_nullify(urban_param)
  call urban_param_dealloc(urban_param_data)

  call psparms_nullify(psparms)
  call psparms_dealloc(psparms_data)

  call jules_vars_dealloc(jules_vars_data)
  call jules_vars_nullify(jules_vars)

  call prognostics_nullify(progs)
  call prognostics_dealloc(progs_data)

  call fluxes_nullify(fluxes)
  call fluxes_dealloc(fluxes_data)

end subroutine sw_rad_tile_code

end module sw_rad_tile_kernel_mod
