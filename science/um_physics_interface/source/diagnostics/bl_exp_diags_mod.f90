!-------------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Processes diagnostics for bl_exp_alg

module bl_exp_diags_mod

  use constants_mod,       only: l_def
  use field_mod,           only: field_type
  use integer_field_mod,   only: integer_field_type
  use io_config_mod,       only: subroutine_timers
  use timer_mod,           only: timer
  use initialise_diagnostics_mod,     only : init_diag => init_diagnostic_field

  implicit none

  private

  ! Logical indicating whether diagnostics are requested
  logical( l_def ) :: gross_prim_prod_flag, zht_flag, z0h_eff_flag, oblen_flag,  &
                      soil_respiration_flag

  public :: initialise_diags_for_bl_exp
  public :: output_diags_for_bl_exp

contains

  !> @brief Initialise fields for locally-computed diagnostics
  !> @param[in]    zh                 Surface-based boundary layer depth
  !> @param[inout] zht                Turbulent mixing height
  !> @param[inout] z0h_eff            Gridbox mean effective roughness length for scalars
  !> @param[inout] gross_prim_prod    Gross Primary Productivity
  !> @param[inout] oblen              Obukhov length
  !> @param[inout] soil_respiration   Soil heterotrophic respiration
  subroutine initialise_diags_for_bl_exp(zh, zht, z0h_eff, gross_prim_prod, &
                                         oblen, soil_respiration)

    implicit none

    type( field_type ), intent(in)    :: zh
    type( field_type ), intent(inout) :: zht
    type( field_type ), intent(inout) :: z0h_eff
    type( field_type ), intent(inout) :: gross_prim_prod
    type( field_type ), intent(inout) :: oblen
    type( field_type ), intent(inout) :: soil_respiration

    if ( subroutine_timers ) call timer("bl_exp_diags")

    zht_flag = init_diag(zht, 'turbulence__zht')
    z0h_eff_flag = init_diag(z0h_eff, 'surface__z0h_eff')
    oblen_flag = init_diag(oblen, 'turbulence__oblen')
    gross_prim_prod_flag = init_diag(gross_prim_prod, 'surface__gross_prim_prod')
    soil_respiration_flag = init_diag(soil_respiration, 'surface__soil_respiration')

    if ( subroutine_timers ) call timer("bl_exp_diags")

  end subroutine initialise_diags_for_bl_exp

  !> @brief Output diagnostics from bl_exp_alg
  !> @param[in] ntml              Number of turbulently mixed levels
  !> @param[in] cumulus           Cumulus flag (true/false)
  !> @param[in] bl_type_ind       Diagnosed BL types
  !> @param[in] wvar              Vertical velocity variance in wth
  !> @param[in] dsldzm            Liquid potential temperature gradient in wth
  !> @param[in] gradrinr          Gradient Richardson number in wth
  !> @param[in] rhokh_bl          Heat eddy diffusivity on BL levels
  !> @param[in] tke_bl            Turbulent kinetic energy (m2 s-2)
  !> @param[in] dtrdz_tq_bl       dt/(rho*r*r*dz) in wth
  !> @param[in] rdz_tq_bl         1/dz in w3
  !> @param[in] zhsc              Height of decoupled layer top
  !> @param[in] level_ent         Level of surface mixed layer inversion
  !> @param[in] level_ent_dsc     Level of decoupled stratocumulus inversion
  !> @param[in] ent_we_lim        Rho * entrainment rate at surface ML inversion (kg m-2 s-1)
  !> @param[in] ent_t_frac        Fraction of time surface ML inversion is above level
  !> @param[in] ent_zrzi          Level height as fraction of DSC inversion height above DSC ML base
  !> @param[in] ent_we_lim_dsc    Rho * entrainment rate at DSC inversion (kg m-2 s-1)
  !> @param[in] ent_t_frac_dsc    Fraction of time DSC inversion is above level
  !> @param[in] ent_zrzi_dsc      Level height as fraction of DSC inversion height above DSC ML base
  !> @param[in] zht               Turbulent mixing height
  !> @param[in] z0h_eff           Gridbox mean effective roughness length for scalars
  !> @param[in] oblen             Obukhov length
  !> @param[in] tile_fraction     Surface tile fractions
  !> @param[in] z0m_tile          tile roughness length for momentum
  !> @param[in] z0m               Cell roughness length
  !> @param[in] canopy_height     Height of the plant canopy
  !> @param[in] leaf_area_index   Leaf area index of vegetated tiles
  !> @param[in] gross_prim_prod   Gross Primary Productivity
  !> @param[in] net_prim_prod     Net Primary Productivity
  !> @param[in] gc_tile           Stomatal conductance on tiles (m s-1)
  !> @param[in] soil_respiration  Soil respiration  (kg m-2 s-1)
  !> @param[in] ustar             surface friction velocity
  !> @param[in] z0m_eff           Gridbox mean effective roughness length for momentum
  !> @param[in] bl_weight_1dbl    Blending weight to 1D BL scheme in the BL
  !> @param[in] dust_flux         Flux of mineral dust by division
  subroutine output_diags_for_bl_exp(ntml, cumulus, bl_type_ind, wvar, dsldzm, &
                                     gradrinr, rhokh_bl, tke_bl, dtrdz_tq_bl,  &
                                     rdz_tq_bl, zhsc, level_ent, level_ent_dsc,&
                                     ent_we_lim, ent_t_frac, ent_zrzi,         &
                                     ent_we_lim_dsc, ent_t_frac_dsc,           &
                                     ent_zrzi_dsc, zht, z0h_eff, oblen,        &
                                     tile_fraction, z0m_tile, z0m,             &
                                     canopy_height, leaf_area_index,           &
                                     gross_prim_prod, net_prim_prod, gc_tile,  &
                                     soil_respiration, ustar, z0m_eff,         &
                                     bl_weight_1dbl, dust_flux)

    implicit none

    ! Prognostic fields to output
    type( field_type ), intent(in)    :: wvar,                                 &
                                         dsldzm, gradrinr, rhokh_bl, tke_bl,   &
                                         dtrdz_tq_bl, rdz_tq_bl, zhsc,         &
                                         zht, z0h_eff,                         &
                                         ent_we_lim, ent_t_frac, ent_zrzi,     &
                                         ent_we_lim_dsc, ent_t_frac_dsc,       &
                                         ent_zrzi_dsc, oblen, bl_weight_1dbl
    type(integer_field_type), intent(in) :: level_ent, level_ent_dsc, ntml,    &
                                            cumulus, bl_type_ind
    type( field_type ), intent(in)    :: tile_fraction, z0m_tile, z0m,         &
                                         gross_prim_prod, net_prim_prod,       &
                                         gc_tile, soil_respiration, ustar,     &
                                         z0m_eff, canopy_height, leaf_area_index
    type( field_type ),  intent(in)   :: dust_flux


    if ( subroutine_timers ) call timer("bl_exp_diags")

    ! Prognostic fields from turbulence collection
    call ntml%write_field('turbulence__ntml')
    call cumulus%write_field('turbulence__cumulus')
    call bl_type_ind%write_field('turbulence__bl_type_ind')
    call wvar%write_field('turbulence__wvar')
    call dsldzm%write_field('turbulence__dsldzm')
    call gradrinr%write_field('turbulence__gradrinr')
    call rhokh_bl%write_field('turbulence__rhokh')
    call tke_bl%write_field('turbulence__tke')
    call dtrdz_tq_bl%write_field('turbulence__dtrdz_tq')
    call rdz_tq_bl%write_field('turbulence__rdz_tq')
    call zhsc%write_field('turbulence__zhsc')
    call level_ent%write_field('turbulence__level_ent')
    call level_ent_dsc%write_field('turbulence__level_ent_dsc')
    call ent_we_lim%write_field('turbulence__ent_we_lim')
    call ent_t_frac%write_field('turbulence__ent_t_frac')
    call ent_zrzi%write_field('turbulence__ent_zrzi')
    call ent_we_lim_dsc%write_field('turbulence__ent_we_lim_dsc')
    call ent_t_frac_dsc%write_field('turbulence__ent_t_frac_dsc')
    call ent_zrzi_dsc%write_field('turbulence__ent_zrzi_dsc')
    call bl_weight_1dbl%write_field('turbulence__bl_weight_1dbl')
    ! Prognostic fields from surface collection
    call tile_fraction%write_field('surface__tile_fraction')
    call z0m_tile%write_field('surface__z0m_tile')
    call z0m%write_field('surface__z0m')
    call z0m_eff%write_field('surface__z0m_eff')
    call canopy_height%write_field('surface__canopy_height')
    call leaf_area_index%write_field('surface__leaf_area_index')
    call net_prim_prod%write_field('surface__net_prim_prod')
    call gc_tile%write_field('surface__gc_tile')
    call ustar%write_field('surface__ustar')
    ! Prognostic fields from aerosol collection
    call dust_flux%write_field('aerosol__dust_flux')

    ! Diagnostics computed in the kernel
    if (zht_flag) call zht%write_field()
    if (oblen_flag) call oblen%write_field()
    if (z0h_eff_flag) call z0h_eff%write_field()
    if (gross_prim_prod_flag) call gross_prim_prod%write_field()
    if (soil_respiration_flag) call soil_respiration%write_field()

    if ( subroutine_timers ) call timer("bl_exp_diags")

  end subroutine output_diags_for_bl_exp
end module bl_exp_diags_mod
