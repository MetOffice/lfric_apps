!-------------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Processes diagnostics for jules_exp_alg

module jules_exp_diags_mod

  use constants_mod,       only: l_def
  use field_mod,           only: field_type
  use integer_field_mod,   only: integer_field_type
  use timing_mod,          only: start_timing, stop_timing, &
                                 tik, LPROF
  use initialise_diagnostics_mod, only : init_diag => init_diagnostic_field, &
                                         samp_diag => diagnostic_to_be_sampled
  use empty_data_mod,             only : empty_real_data
  use extrusion_mod,              only: TWOD
  use mesh_collection_mod,        only: mesh_collection
  use mesh_mod,                   only: mesh_type
  use function_space_collection_mod, only: function_space_collection
  use function_space_mod,            only: function_space_type
  use fs_continuity_mod,             only: W3

  implicit none

  private

  ! Logical indicating whether diagnostics are requested
  logical( l_def ) :: gross_prim_prod_flag, z0h_eff_flag,  &
                      soil_respiration_flag

  public :: initialise_diags_for_jules_exp
  public :: output_diags_for_jules_exp

contains

  !> @brief Initialise fields for locally-computed diagnostics
  !> @param[inout] z0h_eff            Gridbox mean effective roughness length for scalars
  !> @param[inout] gross_prim_prod    Gross Primary Productivity
  !> @param[inout] soil_respiration   Soil heterotrophic respiration
  !> @param[inout] chr10m             10m transfer coefficient
  subroutine initialise_diags_for_jules_exp(z0h_eff, gross_prim_prod, &
                                            soil_respiration, chr10m, mesh)

    implicit none

    type( field_type ), intent(inout) :: z0h_eff
    type( field_type ), intent(inout) :: gross_prim_prod
    type( field_type ), intent(inout) :: soil_respiration
    type( field_type ), intent(inout) :: chr10m

    type( mesh_type ), intent(in), pointer :: mesh
    type(mesh_type), pointer :: twod_mesh
    type(function_space_type), pointer :: vector_space

    integer( tik )                    :: id

    if ( LPROF ) call start_timing( id, 'diags.jules_exp' )

    z0h_eff_flag = init_diag(z0h_eff, 'surface__z0h_eff')
    gross_prim_prod_flag = init_diag(gross_prim_prod, 'surface__gross_prim_prod')
    soil_respiration_flag = init_diag(soil_respiration, 'surface__soil_respiration')
    ! chr10m is required if either t10m_ssi or q10m_ssi are requested
    ! but needed as an empty field otherwise
    twod_mesh => mesh_collection%get_mesh(mesh, TWOD)
    vector_space => function_space_collection%get_fs(twod_mesh, 0, 0, W3)
    if (samp_diag('surface__t10m_ssi') .or. samp_diag('surface__q10m_ssi')) then
      call chr10m%initialise(vector_space, name='chr10m')
    else
      call chr10m%initialise(vector_space, name='chr10m', &
                             override_data = empty_real_data)
    end if

    if ( LPROF ) call stop_timing( id, 'diags.jules_exp' )

  end subroutine initialise_diags_for_jules_exp

  !> @brief Output diagnostics from jules_exp_alg
  !> @param[in] z0h_eff           Gridbox mean effective roughness length for scalars
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
  !> @param[in] dust_flux         Flux of mineral dust by division
  subroutine output_diags_for_jules_exp(z0h_eff, tile_fraction, z0m_tile, z0m, &
                                     canopy_height, leaf_area_index,           &
                                     gross_prim_prod, net_prim_prod, gc_tile,  &
                                     soil_respiration, ustar, z0m_eff,         &
                                     dust_flux)

    implicit none

    ! Prognostic fields to output
    type( field_type ), intent(in)    :: z0h_eff
    type( field_type ), intent(in)    :: tile_fraction, z0m_tile, z0m,         &
                                         gross_prim_prod, net_prim_prod,       &
                                         gc_tile, soil_respiration, ustar,     &
                                         z0m_eff, canopy_height, leaf_area_index
    type( field_type ),  intent(in)   :: dust_flux
    integer( tik )                    :: id

    if ( LPROF ) call start_timing( id, 'diags.jules_exp' )

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
    if (z0h_eff_flag) call z0h_eff%write_field()
    if (gross_prim_prod_flag) call gross_prim_prod%write_field()
    if (soil_respiration_flag) call soil_respiration%write_field()

    if ( LPROF ) call stop_timing( id, 'diags.jules_exp' )

  end subroutine output_diags_for_jules_exp
end module jules_exp_diags_mod
