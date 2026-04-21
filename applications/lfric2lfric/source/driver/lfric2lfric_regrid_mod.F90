!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Performs regridding of a field collection
!>
module lfric2lfric_regrid_mod

  use base_mesh_config_mod,     only: geometry_spherical,      &
                                      geometry_planar,         &
                                      topology_fully_periodic, &
                                      topology_non_periodic
  use constants_mod,            only: str_def, i_def
  use driver_modeldb_mod,       only: modeldb_type
  use field_parent_mod,         only: field_parent_type
  use field_mod,                only: field_type
  use field_collection_mod,     only: field_collection_type
  use field_collection_iterator_mod, only: &
                                      field_collection_iterator_type
  use fs_continuity_mod,        only: W2, W3, Wtheta
  use function_space_collection_mod, only: function_space_collection
  use function_space_mod,       only: function_space_type
  use interpolation_alg_mod,    only: interp_w2_to_w3wth_alg, &
                                      interp_w3wth_to_w2_alg
  use log_mod,                  only: log_event, &
                                      log_level_info, &
                                      log_scratch_space
  use mesh_collection_mod,      only: mesh_collection
  use mesh_mod,                 only: mesh_type
  use model_clock_mod,          only: model_clock_type
  use namelist_mod,             only: namelist_type

  !------------------------------------
  ! lfric2lfric modules
  !------------------------------------
  use lfric2lfric_config_mod,         only: regrid_method_map,         &
                                            regrid_method_lfric2lfric, &
                                            regrid_method_oasis
  use lfric2lfric_map_regrid_mod,     only: lfric2lfric_map_regrid
  use lfric2lfric_no_regrid_mod,      only: lfric2lfric_no_regrid
  use lfric2lfric_oasis_regrid_mod,   only: lfric2lfric_oasis_regrid

  implicit none

  private
  public lfric2lfric_regrid

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief    Performs regridding of a field collection
  !> @details  Over one time step, all regridding is performed by algorithm
  !!           modules specific to the regrid method defined in the
  !!           configuration.
  !!           Fields to be regridded are extracted from the source field
  !!           collection. Corresponding source and target field pairs
  !!           are passed to the regridding algorithm, and written to fields
  !!           in the destination field collection.
  !> @param [in,out] modeldb                 The structure that holds model
  !!                                         state
  !> @param [in,out] oasis_clock             Clock for OASIS exchanges
  !> @param [in]     source_fields           Collection of fields to be regridded
  !> @param [in]     target_fields           Collection of regridded fields
  !> @param [in]     regrid_method           Method for regridding between the
  !>                                         source and destination meshes
  subroutine lfric2lfric_regrid( modeldb, oasis_clock,            &
                  source_fields, target_fields, regrid_method )

    implicit none

    type(modeldb_type),                   intent(inout) :: modeldb
    type(model_clock_type), allocatable,  intent(inout) :: oasis_clock
    type(field_collection_type), pointer, intent(in)    :: source_fields
    type(field_collection_type), pointer, intent(inout) :: target_fields
    integer(kind=i_def),                  intent(in)    :: regrid_method


    type(field_collection_iterator_type) :: iter

    class(field_parent_type),  pointer :: field
    type(field_type),          pointer :: field_src
    type(field_type),          pointer :: field_dst

    type(field_type),          target  :: u_in_w3_src,  u_in_w3_dst
    type(field_type),          target  :: v_in_w3_src,  v_in_w3_dst
    type(field_type),          target  :: w_in_wth_src, w_in_wth_dst

    integer(kind=i_def)                :: fs_id
    type(function_space_type), pointer :: fs_w3_src, fs_w3_dst
    type(function_space_type), pointer :: fs_wth_src, fs_wth_dst
    type(mesh_type),           pointer :: mesh_dst

    character(len=str_def)   :: mesh_names(2)
    integer(kind=i_def)      :: element_order_h
    integer(kind=i_def)      :: element_order_v
    integer(kind=i_def)      :: geometry
    integer(kind=i_def)      :: topology

    character(len=str_def)   :: field_name

    integer(kind=i_def), parameter :: dst = 1
    integer(kind=i_def), parameter :: src = 2


    ! Obtain namelist parameters
    mesh_names(dst) = modeldb%config%lfric2lfric%destination_mesh_name()
    mesh_names(src) = modeldb%config%lfric2lfric%source_mesh_name()
    element_order_h = modeldb%config%finite_element%element_order_h()
    element_order_v = modeldb%config%finite_element%element_order_v()

    ! Function spaces for creating temporary fields used in W2 interpolation
    fs_w3_src => function_space_collection%get_fs(                   &
                          mesh_collection%get_mesh(mesh_names(src)), &
                          element_order_h, element_order_v, W3)
    fs_wth_src => function_space_collection%get_fs(                  &
                          mesh_collection%get_mesh(mesh_names(src)), &
                          element_order_h, element_order_v, Wtheta)
    fs_w3_dst => function_space_collection%get_fs(                   &
                          mesh_collection%get_mesh(mesh_names(dst)), &
                          element_order_h, element_order_v, W3)
    fs_wth_dst => function_space_collection%get_fs(                  &
                          mesh_collection%get_mesh(mesh_names(dst)), &
                          element_order_h, element_order_v, Wtheta)

    ! Get the geometry and topology of the destination mesh
    mesh_dst => mesh_collection%get_mesh(trim(mesh_names(dst)))
    if (mesh_dst%is_geometry_spherical()) then
      geometry = geometry_spherical
    else if (mesh_dst%is_geometry_planar()) then
      geometry = geometry_planar
    end if
    if (mesh_dst%is_topology_periodic()) then
      topology = topology_fully_periodic
    else if (mesh_dst%is_topology_non_periodic()) then
      topology = topology_non_periodic
    end if

    ! Main loop over fields to be processed
    call iter%initialise(source_fields)
    do
      ! Locate the field to be processed in the field collections
      if ( .not.iter%has_next() ) exit
      field => iter%next()
      field_name = field%get_name()

      call source_fields%get_field(field_name, field_src)
      call target_fields%get_field(field_name, field_dst)

      write(log_scratch_space, '(A,A)') "Processing lfric field ", &
                                           trim(field_name)
      call log_event(log_scratch_space, log_level_info)

      ! Convert W2 fields to a set of W3 and Wtheta fields
      fs_id = field_src%which_function_space()
      if (fs_id == W2) then
        call u_in_w3_src%initialise(vector_space=fs_w3_src,   &
                                    name="u_in_w3_src")
        call v_in_w3_src%initialise(vector_space=fs_w3_src,   &
                                    name="v_in_w3_src")
        call w_in_wth_src%initialise(vector_space=fs_wth_src, &
                                    name="w_in_wth_src")
        call u_in_w3_dst%initialise(vector_space=fs_w3_dst,   &
                                    name="u_in_w3_dst")
        call v_in_w3_dst%initialise(vector_space=fs_w3_dst,   &
                                    name="v_in_w3_dst")
        call w_in_wth_dst%initialise(vector_space=fs_wth_dst, &
                                    name="w_in_wth_dst")

        call interp_w2_to_w3wth_alg(field_src, u_in_w3_src,   &
                                    v_in_w3_src, w_in_wth_src)
      end if

      ! Regrid source field depending on regrid method
      select case (regrid_method)
        case (regrid_method_map)
          if (fs_id == W2) then
            call lfric2lfric_map_regrid(u_in_w3_dst, u_in_w3_src)
            call lfric2lfric_map_regrid(v_in_w3_dst, v_in_w3_src)
            call lfric2lfric_map_regrid(w_in_wth_dst, w_in_wth_src)
          else
            call lfric2lfric_map_regrid(field_dst, field_src)
          end if

        case (regrid_method_lfric2lfric)
          write(log_scratch_space, '(A)') &
                              'Regrid method lfric2lfric not implemented yet'
          call log_event(log_scratch_space, log_level_info)

          call lfric2lfric_no_regrid(field_dst)

        case (regrid_method_oasis)
#ifdef MCT
          if (fs_id == W2) then
            call lfric2lfric_oasis_regrid(modeldb, oasis_clock, &
                                  u_in_w3_dst, u_in_w3_src)
            call lfric2lfric_oasis_regrid(modeldb, oasis_clock, &
                                  v_in_w3_dst, v_in_w3_src)
            call lfric2lfric_oasis_regrid(modeldb, oasis_clock, &
                                  w_in_wth_dst, w_in_wth_src)
          else
            call lfric2lfric_oasis_regrid(modeldb, oasis_clock, &
                                          field_dst, field_src)
          end if
#endif
      end select

      ! Rebuild the W2 fields from a set of W3 and Wtheta fields
      if (fs_id == W2) then
        call interp_w3wth_to_w2_alg(field_dst, u_in_w3_dst,    &
                                    v_in_w3_dst, w_in_wth_dst, &
                                    geometry, topology)
      end if
    end do

  end subroutine lfric2lfric_regrid

end module lfric2lfric_regrid_mod
