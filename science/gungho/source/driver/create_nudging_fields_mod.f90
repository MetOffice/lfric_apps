!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Create fields used for spectral nudging
module create_nudging_fields_mod

  use clock_mod,                   only : clock_type
  use config_mod,                  only : config_type
  use constants_mod,               only : i_def, l_def, str_def, EMDI
  use driver_modeldb_mod,          only : modeldb_type
  use field_mod,                   only : field_type
  use field_collection_mod,        only : field_collection_type
  use field_mapper_mod,            only : field_mapper_type
  use field_maker_mod,             only : field_maker_type
  use fs_continuity_mod,           only : W3
  use function_space_mod,          only : function_space_type
  use gungho_time_axes_mod,        only : gungho_time_axes_type
  use log_mod,                     only : log_event, log_scratch_space,        &
                                          LOG_LEVEL_TRACE, LOG_LEVEL_ERROR
  use mesh_mod,                    only : mesh_type

  use external_forcing_config_mod, only : theta_forcing_nudging,               &
                                          wind_forcing_nudging,                &
                                          external_forcing_is_loaded

  implicit none

  public :: create_nudging_fields
  public :: process_nudging_fields

  contains

  !> @brief   Create and add nudging fields.
  !> @details Create reference fields to for nudging in the derived field
  !!          collection. On every timestep these fields will be updated.
  !> @param[in]    modeldb    The dataset for this run
  !> @param[in]    mesh       The current 3d mesh
  !> @param[in]    twod_mesh  The current 2d mesh
  !> @param[in]    mapper     Provides access to the field collections
  subroutine create_nudging_fields(modeldb, mesh, twod_mesh, mapper)
    implicit none
    type(modeldb_type), intent(in)          :: modeldb
    type(mesh_type), intent(in), pointer    :: mesh
    type(mesh_type), intent(in), pointer    :: twod_mesh
    type(field_mapper_type), intent(in)     :: mapper

    type(gungho_time_axes_type), pointer    :: gungho_axes

    type(field_maker_type) :: creator

    call log_event('GungHo: Creating nudging fields...', LOG_LEVEL_TRACE)

    gungho_axes => mapper%get_gungho_axes()
    call gungho_axes%make_nudging_time_axis()

    call creator%init(mesh, twod_mesh, mapper, modeldb%clock)

    call process_nudging_fields(modeldb%config, creator)

    call gungho_axes%save_nudging_time_axis()

  end subroutine create_nudging_fields

  !> @brief Iterate over active nudging fields and apply an arbitrary
  !! processor to the field specifiers.
  !> @param        config      Argument holding model configuration
  !> @param        processor   Processor to be applied to field specifiers
  subroutine process_nudging_fields(config, processor)
    use field_spec_mod,                only : main => main_coll_dict,          &
                                              axis => time_axis_dict,          &
                                              processor_type,                  &
                                              make_spec
    implicit none

    type(config_type),  intent(in)   :: config
    class(processor_type)            :: processor

    logical(l_def)                   :: coarse_nudging
    character(str_def)               :: nudging_mesh_name
    integer(i_def)                   :: theta_forcing
    integer(i_def)                   :: wind_forcing

    if ( config%formulation%use_multires_coupling() ) then
      coarse_nudging = config%multires_coupling%coarse_nudging()
      nudging_mesh_name = config%multires_coupling%nudging_mesh_name()
    else
      coarse_nudging = .false.
    end if

    if ( external_forcing_is_loaded() ) then
      theta_forcing = config%external_forcing%theta_forcing()
      wind_forcing = config%external_forcing%wind_forcing()
    else
      theta_forcing = EMDI
      wind_forcing = EMDI
    end if

    !------ Fields updated directly from nudging file-----------------

    ! The surface pressure field from data file is used to re-create the
    !  vertical grid (sigma or pressure-levels) of the source data while
    !  remapping from nudging to model grid. Hence this field is needed
    !  when any one of theta or wind is forced.
    call processor%apply(make_spec(                                            &
            'surface_pressure_nudging_ext_ref', main%none, W3,                 &
            coarse=coarse_nudging,                                             &
            coarse_mesh_name=nudging_mesh_name,                                &
            twod=.true., time_axis=axis%nudging                                &
    ))

    if ( theta_forcing == theta_forcing_nudging ) then
      call processor%apply(make_spec(                                          &
              'temperature_nudging_ext_ref', main%none, W3, twod=.true.,       &
              coarse=coarse_nudging,                                           &
              coarse_mesh_name=nudging_mesh_name,                              &
              mult='nudging_levels', time_axis=axis%nudging                    &
      ))
    end if

    if ( wind_forcing == wind_forcing_nudging ) then
      call processor%apply(make_spec(                                          &
              'u_nudging_ext_ref', main%none, W3, twod=.true.,                 &
              coarse=coarse_nudging,                                           &
              coarse_mesh_name=nudging_mesh_name,                              &
              mult='nudging_levels', time_axis=axis%nudging                    &
      ))
      call processor%apply(make_spec(                                          &
              'v_nudging_ext_ref', main%none, W3, twod=.true.,                 &
              coarse=coarse_nudging,                                           &
              coarse_mesh_name=nudging_mesh_name,                              &
              mult='nudging_levels', time_axis=axis%nudging                    &
      ))
    end if

  end subroutine process_nudging_fields

end module create_nudging_fields_mod
