!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Create fields used for spectral nudging
module create_nudging_fields_mod

  use clock_mod,                   only : clock_type
  use constants_mod,               only : i_def, l_def, str_def
  use field_mod,                   only : field_type
  use field_collection_mod,        only : field_collection_type
  use field_mapper_mod,            only : field_mapper_type
  use field_maker_mod,             only : field_maker_type
  use fs_continuity_mod,           only : W3
  use function_space_mod,          only : function_space_type
  use gungho_time_axes_mod,        only : gungho_time_axes_type
  use log_mod,                     only : log_event, log_scratch_space,        &
                                          LOG_LEVEL_INFO, LOG_LEVEL_ERROR
  use mesh_mod,                    only : mesh_type

  use external_forcing_config_mod, only : theta_forcing_nudging,               &
                                          theta_forcing,                       &
                                          wind_forcing_nudging,                &
                                          wind_forcing
  implicit none

  public :: create_nudging_fields
  public :: process_nudging_fields

  contains

  !> @brief   Create and add nudging fields.
  !> @details Create reference fields to for nudging in the derived field
  !!          collection. On every timestep these fields will be updated.
  !> @param[in]    mesh       The current 3d mesh
  !> @param[in]    twod_mesh  The current 2d mesh (not used here)
  !> @param[in]    mapper     Provides access to the field collections
  !> @param[in]    clock      The model clock
  subroutine create_nudging_fields(mesh, twod_mesh, mapper, clock)
    implicit none
    type(mesh_type), intent(in), pointer    :: mesh
    type(mesh_type), intent(in), pointer    :: twod_mesh
    type(field_mapper_type), intent(in)     :: mapper
    class(clock_type), intent(in)           :: clock

    type(gungho_time_axes_type), pointer    :: gungho_axes

    type(field_maker_type) :: creator

    call log_event('GungHo: Creating nudging fields...', LOG_LEVEL_INFO)

    gungho_axes => mapper%get_gungho_axes()
    call gungho_axes%make_nudging_time_axis()

    call creator%init(mesh, twod_mesh, mapper, clock)

    call process_nudging_fields(creator)

    call gungho_axes%save_nudging_time_axis()

  end subroutine create_nudging_fields

  !> @brief Iterate over active nudging fields and apply an arbitrary
  !! processor to the field specifiers.
  subroutine process_nudging_fields(processor)
    use field_spec_mod,                only : main => main_coll_dict,          &
                                              axis => time_axis_dict,          &
                                              processor_type,                  &
                                              make_spec
    use multires_coupling_config_mod,  only : coarse_nudging,                  &
                                              nudging_mesh_name

    implicit none

    class(processor_type) :: processor

    !------ Fields updated directly from nudging file-----------------
    call processor%apply(make_spec(                                            &
            'surface_pressure_nudging_ext_ref', main%derived, W3,              &
            coarse=coarse_nudging,                                             &
            coarse_mesh_name=nudging_mesh_name,                                &
            twod=.true., time_axis=axis%nudging                                &
    ))

    if ( theta_forcing == theta_forcing_nudging ) then
      call processor%apply(make_spec(                                          &
              'temperature_nudging_ext_ref', main%derived, W3, twod=.true.,    &
              coarse=coarse_nudging,                                           &
              coarse_mesh_name=nudging_mesh_name,                              &
              mult='ecmwf_levels', time_axis=axis%nudging                      &
      ))
    end if

    if ( wind_forcing == wind_forcing_nudging ) then
      call processor%apply(make_spec(                                          &
              'u_nudging_ext_ref', main%derived, W3, twod=.true.,              &
              coarse=coarse_nudging,                                           &
              coarse_mesh_name=nudging_mesh_name,                              &
              mult='ecmwf_levels', time_axis=axis%nudging                      &
      ))
      call processor%apply(make_spec(                                          &
              'v_nudging_ext_ref', main%derived, W3, twod=.true.,              &
              coarse=coarse_nudging,                                           &
              coarse_mesh_name=nudging_mesh_name,                              &
              mult='ecmwf_levels', time_axis=axis%nudging                      &
      ))
    end if

  end subroutine process_nudging_fields

end module create_nudging_fields_mod
