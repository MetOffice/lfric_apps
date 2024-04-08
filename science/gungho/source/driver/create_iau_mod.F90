!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Create IAU fields.
!> @details Create IAU field collection and add fields.
module create_iau_mod

  use constants_mod,              only : i_def, l_def, str_def
  use log_mod,                    only : log_event,             &
                                         log_scratch_space,     &
                                         LOG_LEVEL_INFO
  use field_mod,                  only : field_type
  use field_parent_mod,           only : read_interface
  use field_collection_mod,       only : field_collection_type
  use fs_continuity_mod,          only : W3, Wtheta
  use mesh_mod,                   only : mesh_type
  use gungho_modeldb_mod,         only : modeldb_type
  use lfric_xios_read_mod,        only : read_field_generic
  use init_time_axis_mod,         only : setup_field
#ifdef UM_PHYSICS
  use iau_config_mod,             only : iau_wet_density
  use nlsizes_namelist_mod,       only : sm_levels
  use jules_control_init_mod,     only : n_land_tile
  use jules_physics_init_mod,     only : snow_lev_tile
#endif

  implicit none

  public  :: create_iau_fields

  contains

  !> @brief   Create and add iau fields.
  !> @details Create IAU  field collection.
  !!          On every timestep these fields will be updated and used by the  model.
  !> @param[in]     mesh              The current 3d mesh identifier.
  !> @param[in]     twod_mesh         The current 2d mesh identifier
  !> @param[in,out] depository        Main collection of all fields in memory.
  !> @param[in,out] prognostic_fields The prognostic variables in the model.
  !> @param[in,out] modeldb           The model database
  subroutine create_iau_fields( mesh, twod_mesh, depository,   &
                                prognostic_fields, modeldb )

    implicit none

    type(mesh_type),             intent(in), pointer :: mesh
    type(mesh_type),             intent(in), pointer :: twod_mesh
    type(field_collection_type), intent(inout)       :: depository
    type(field_collection_type), intent(inout)       :: prognostic_fields
    type(modeldb_type),          intent(inout)       :: modeldb

    type(field_collection_type), pointer             :: iau_fields => null()
    type(field_collection_type), pointer             :: iau_surf_fields => null()
    logical(l_def)                                   :: checkpoint_restart_flag
    logical(l_def),              parameter           :: cyclic=.false.
    logical(l_def),              parameter           :: interp_flag=.true.
    procedure(read_interface),   pointer             :: read_behaviour => null()

    call log_event( 'Create IAU fields', LOG_LEVEL_INFO )

    call modeldb%fields%add_empty_field_collection("iau_fields")
    iau_fields => modeldb%fields%get_field_collection("iau_fields")

    ! Checkpoint all necessary IAU fields.
    ! Only the option to read in iau fields from file for now.

    checkpoint_restart_flag = .false.
    read_behaviour => read_field_generic

    call setup_field( iau_fields, depository, prognostic_fields, &
        "u_in_w3_inc", W3, mesh, checkpoint_restart_flag, &
        read_behaviour=read_behaviour )

    call setup_field( iau_fields, depository, prognostic_fields, &
        "v_in_w3_inc", W3, mesh, checkpoint_restart_flag, &
         read_behaviour=read_behaviour )

    call setup_field( iau_fields, depository, prognostic_fields, &
        "q_inc", Wtheta, mesh, checkpoint_restart_flag, &
        read_behaviour=read_behaviour )

    call setup_field( iau_fields, depository, prognostic_fields, &
        "qcl_inc", Wtheta, mesh, checkpoint_restart_flag, &
        read_behaviour=read_behaviour )

    call setup_field( iau_fields, depository, prognostic_fields, &
        "qcf_inc", Wtheta, mesh, checkpoint_restart_flag, &
        read_behaviour=read_behaviour )

#ifdef UM_PHYSICS
    if ( iau_wet_density ) then
      call setup_field( iau_fields, depository, prognostic_fields, &
          "rho_r2_inc", W3, mesh, checkpoint_restart_flag, &
          read_behaviour=read_behaviour )
    else
      call setup_field( iau_fields, depository, prognostic_fields, &
          "rho_inc", W3, mesh, checkpoint_restart_flag, &
          read_behaviour=read_behaviour )
    end if
#endif

    call setup_field( iau_fields, depository, prognostic_fields, &
        "theta_inc", Wtheta, mesh, checkpoint_restart_flag, &
        read_behaviour=read_behaviour )

    call setup_field( iau_fields, depository, prognostic_fields, &
        "exner_inc", W3, mesh, checkpoint_restart_flag, &
        read_behaviour=read_behaviour )

    call log_event( 'Create IAU land surface fields', LOG_LEVEL_INFO )

    call modeldb%fields%add_empty_field_collection("iau_surf_fields")
    iau_surf_fields => modeldb%fields%get_field_collection("iau_surf_fields")

    checkpoint_restart_flag = .false.
    read_behaviour => read_field_generic

#ifdef UM_PHYSICS
    call setup_field( iau_surf_fields, depository, prognostic_fields, &
         "soil_temperature_inc", W3, mesh, checkpoint_restart_flag, &
         twod_mesh, read_behaviour=read_behaviour, twod=.true., ndata=sm_levels )

    call setup_field( iau_surf_fields, depository, prognostic_fields, &
         "soil_moisture_inc", W3, mesh, checkpoint_restart_flag, &
         twod_mesh, read_behaviour=read_behaviour, twod=.true., ndata=sm_levels )

    call setup_field( iau_surf_fields, depository, prognostic_fields, &
         "snow_layer_temp_inc", W3, mesh, checkpoint_restart_flag, &
         twod_mesh, read_behaviour=read_behaviour, twod=.true., ndata=snow_lev_tile )

    call setup_field( iau_surf_fields, depository, prognostic_fields, &
         "tile_temperature_inc", W3, mesh, checkpoint_restart_flag, &
         twod_mesh, read_behaviour=read_behaviour, twod=.true., ndata=n_land_tile )
#endif

    nullify( iau_fields, iau_surf_fields )

   end subroutine create_iau_fields

end module create_iau_mod
