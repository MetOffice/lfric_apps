!----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Controls the setting of Jules high level variables which are either
!>         fixed in LFRic or derived from LFRic inputs.

module jules_control_init_mod

  ! LFRic namelists which have been read
  use surface_config_mod,          only : n_sea_ice_tile_in => n_sea_ice_tile

  ! Other LFRic modules used
  use constants_mod,        only : i_def

  implicit none

  integer(kind=i_def), parameter :: n_sea_tile  = 1

  ! This is the size of array surf_interp, 2D variables that need interpolating
  ! to cell faces in um_physics/source/kernel/bl_exp_kernel_mod.F90
  ! Has nothing to do with number of surface types
  integer(kind=i_def), parameter :: n_surf_interp = 11

  integer(kind=i_def), protected :: n_land_tile
  integer(kind=i_def), protected :: n_sea_ice_tile
  integer(kind=i_def), protected :: n_surf_tile
  integer(kind=i_def), protected :: first_sea_tile
  integer(kind=i_def), protected :: first_sea_ice_tile
  integer(kind=i_def), protected :: soil_lev_tile

  private
  public :: n_land_tile, n_sea_tile, n_sea_ice_tile, n_surf_tile, &
       first_sea_tile, first_sea_ice_tile, jules_control_init,    &
       soil_lev_tile, n_surf_interp

contains

  !>@brief Initialise JULES high levels variables with are either fixed in LFRic
  !>        or derived from LFRic inputs.
  !>@details Most variables in this file need to be set consistent with the
  !>          ancillary files which are provided. As we do not yet have access
  !>          to these, we set the "variables" as parameters here. Hopefully
  !>          they can be read directly from the ancillary file header in
  !>          due course, but if not, they will need to be promoted to the
  !>          LFRic namelists. We then derive other JULES information
  !>          from these parameters.
  subroutine jules_control_init()

    ! LFRic namelists which have been read
    use jules_surface_types_config_mod, only : npft_in => npft, &
                                               nnvg_in => nnvg,     &
                                               brd_leaf_in => brd_leaf, &
                                               ndl_leaf_in => ndl_leaf, &
                                               c3_grass_in => c3_grass, &
                                               c4_grass_in => c4_grass, &
                                               shrub_in => shrub,    &
                                               urban_in => urban,    &
                                               urban_canyon_in => urban_canyon,&
                                               urban_roof_in => urban_roof,    &
                                               lake_in => lake,     &
                                               soil_in => soil,     &
                                               ice_in => ice
    use section_choice_config_mod, only : surface, surface_jules

    ! UM/JULES modules containing things that need setting
    use ancil_info, only: jules_dim_cs1 => dim_cs1, nsurft
    use atm_step_local, only: co2_dim_len, co2_dim_row, &
        dim_cs1
    use jules_soil_mod, only: jules_sm_levels => sm_levels
    use jules_surface_types_mod, only: npft, nnvg, ntype, brd_leaf, ndl_leaf,  &
        c3_grass, c4_grass, shrub, urban, urban_canyon, urban_roof, lake,      &
        soil, ice
    use jules_vegetation_mod, only: l_triffid
    use jules_model_environment_mod, only: lsm_id, jules
    use nlsizes_namelist_mod, only: ntiles, sm_levels
    use jules_surface_types_mod, only:                                         &
        set_derived_variables_jules_surface_types,                             &
        print_nlist_jules_surface_types, check_jules_surface_types
    use land_tile_ids_mod, only: set_surface_type_ids

    use log_mod, only : log_event, log_scratch_space, LOG_LEVEL_INFO

    implicit none

    call log_event( 'jules_control_init', LOG_LEVEL_INFO )

    ! ----------------------------------------------------------------
    ! Surface tile information - set in JULES module jules_surface_types
    ! ----------------------------------------------------------------

    if (surface == surface_jules) then
      write( log_scratch_space, '(A)' ) 'JULES surface scheme is being used.'
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
      npft         = npft_in
      nnvg         = nnvg_in
      brd_leaf     = brd_leaf_in
      ndl_leaf     = ndl_leaf_in
      c3_grass     = c3_grass_in
      c4_grass     = c4_grass_in
      shrub        = shrub_in
      urban        = urban_in
      urban_canyon = urban_canyon_in
      urban_roof   = urban_roof_in
      lake         = lake_in
      soil         = soil_in
      ice          = ice_in
      ! Calculate ntype and nnpft
      call set_derived_variables_jules_surface_types()
      call print_nlist_jules_surface_types()
      call set_surface_type_ids()
      call check_jules_surface_types()
      n_sea_ice_tile = n_sea_ice_tile_in
    else
      write( log_scratch_space, '(A)' ) 'No surface scheme is being used.'
      call log_event( log_scratch_space, LOG_LEVEL_INFO )
      ! create_physics_prognostics needs something to be allocated for
      ! pft_space, surft_space and sice_space
      npft           = 1
      ntype          = 1
      n_sea_ice_tile = 1
    end if

    ! Number of land tiles; currently called different things in different
    ! places. Should probably be consolidated to only use one in LFRic
    ! interface; n_land_tile probably as ntiles now exists only in the UM, while
    ! nsurft in used in JULES. nsurft would be prefered, although it could be
    ! confused with n_surf_tile, which includes the sea and sea-ice tiles in
    ! LFRic.
    n_land_tile = ntype
    nsurft      = ntype
    ntiles      = ntype

    ! Total number of surface tiles, used to dimension LFRic
    ! multidata fields
    n_surf_tile = n_land_tile + n_sea_tile + n_sea_ice_tile

    ! Indices of the first sea and sea-ice tiles. By convection the tile
    ! order is always land, sea, sea-ice
    first_sea_tile     = n_land_tile + 1
    first_sea_ice_tile = n_land_tile + n_sea_tile + 1

    ! ----------------------------------------------------------------
    ! Model dimensions - in each case, the first variable is
    !  contained in UM module nlsizes_namelist_mod. It must then be
    !  copied across into variables which live in JULES modules so that
    !  allocate_jules_arrays can access via modules. Ultimately the
    !  UM variables should be removed and only the JULES ones will exist.
    ! ----------------------------------------------------------------
    ! Number of soil levels - set to a constant as this rarely changes
    ! but may migrate to namelist or read from ancillary in due course.
    sm_levels       = 4
    jules_sm_levels = sm_levels

    ! Product of soil levels and land tiles for water extraction
    soil_lev_tile = sm_levels * n_land_tile

    ! ----------------------------------------------------------------
    ! More model dimensions, this time from atm_step_local
    ! ----------------------------------------------------------------
    ! Dimensions for triffid - this is not yet implemented, but we set
    ! the variables correctly based on the use or not of triffid here
    ! so that hopefully it works when it is implemented.
    if (l_triffid) then
      dim_cs1       = 4
    else
      dim_cs1       = 1
    end if

    ! Now pass into JULES module variables
    jules_dim_cs1 = dim_cs1

    ! Dimensions of co2 array - set to 1 to match kernel size, but may change
    ! if multiple cells are passed to kernels.
    co2_dim_len = 1
    co2_dim_row = 1

    ! Initialise LSM to be JULES (other options do exist; CABLE, Rivers-only)
    lsm_id = jules

  end subroutine jules_control_init

end module jules_control_init_mod
