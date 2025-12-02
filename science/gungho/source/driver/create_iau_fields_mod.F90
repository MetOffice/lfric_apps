!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Create IAU fields.
!> @details Create IAU field collection and add fields.
module create_iau_fields_mod

  use constants_mod,                 only : i_def, l_def, str_def
  use log_mod,                       only : log_event,         &
                                            log_scratch_space, &
                                            LOG_LEVEL_INFO
  use field_mod,                     only : field_type
  use field_parent_mod,              only : read_interface
  use field_collection_mod,          only : field_collection_type
  use finite_element_config_mod,     only : element_order_h, element_order_v
  use function_space_collection_mod, only : function_space_collection
  use function_space_mod,            only : function_space_type
  use fs_continuity_mod,             only : W3, Wtheta
  use mesh_mod,                      only : mesh_type
  use driver_modeldb_mod,            only : modeldb_type
  use lfric_xios_read_mod,           only : read_field_generic
  use init_time_axis_mod,            only : setup_field
  use section_choice_config_mod,     only : iau_surf, &
                                            iau_sst
#ifdef UM_PHYSICS
  use iau_config_mod,                only : iau_use_addinf,  &
                                            iau_use_bcorr,   &
                                            iau_use_pertinc, &
                                            iau_wet_density
  use nlsizes_namelist_mod,          only : sm_levels
  use jules_control_init_mod,        only : n_land_tile
  use jules_physics_init_mod,        only : snow_lev_tile
#endif

  implicit none

  private

  public  :: create_iau_fields
#ifdef UM_PHYSICS
  private :: create_iau_additional_fields
#endif

  contains

  !> @brief   Create and add iau fields.
  !> @details Create IAU field collections to hold increments.
  !> @param[in]     mesh              The current 3d mesh identifier.
  !> @param[in]     twod_mesh         The current 2d mesh identifier
  !> @param[in,out] depository        Main collection of all fields in memory.
  !> @param[in,out] prognostic_fields The prognostic variables in the model.
  !> @param[in,out] modeldb           The model database
  subroutine create_iau_fields( mesh, twod_mesh, depository, &
                                prognostic_fields, modeldb )

    implicit none

    type(mesh_type),             intent(in), pointer :: mesh
    type(mesh_type),             intent(in), pointer :: twod_mesh
    type(field_collection_type), intent(inout)       :: depository
    type(field_collection_type), intent(inout)       :: prognostic_fields
    type(modeldb_type),          intent(inout)       :: modeldb

    type(field_collection_type), pointer             :: iau_fields
    type(field_collection_type), pointer             :: iau_surf_fields
    type(field_collection_type), pointer             :: iau_sst_fields
#ifdef UM_PHYSICS
    type(field_collection_type), pointer             :: iau_pert_fields
    character(str_def)                               :: rho_pert_inc
    character(str_def)                               :: rho_inc_name
#endif
    logical(l_def)                                   :: checkpoint_restart_flag
    logical(l_def),              parameter           :: cyclic=.false.
    logical(l_def),              parameter           :: interp_flag=.true.
    procedure(read_interface),   pointer             :: read_behaviour

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

#ifdef UM_PHYSICS
    ! IAU control-pert pert increment fields
    if ( iau_use_pertinc ) then
      call log_event( 'Create IAU pert fields', LOG_LEVEL_INFO )

      call modeldb%fields%add_empty_field_collection("iau_pert_fields")
      iau_pert_fields => modeldb%fields%get_field_collection("iau_pert_fields")

      checkpoint_restart_flag = .false.
      read_behaviour => read_field_generic

      call setup_field( iau_pert_fields, depository, prognostic_fields, &
          "u_in_w3_pert_inc", W3, mesh, checkpoint_restart_flag, &
          read_behaviour=read_behaviour )

      call setup_field( iau_pert_fields, depository, prognostic_fields, &
          "v_in_w3_pert_inc", W3, mesh, checkpoint_restart_flag, &
           read_behaviour=read_behaviour )

      call setup_field( iau_pert_fields, depository, prognostic_fields, &
          "q_pert_inc", Wtheta, mesh, checkpoint_restart_flag, &
          read_behaviour=read_behaviour )

      call setup_field( iau_pert_fields, depository, prognostic_fields, &
          "qcl_pert_inc", Wtheta, mesh, checkpoint_restart_flag, &
          read_behaviour=read_behaviour )

      call setup_field( iau_pert_fields, depository, prognostic_fields, &
          "qcf_pert_inc", Wtheta, mesh, checkpoint_restart_flag, &
          read_behaviour=read_behaviour )

      if ( iau_wet_density ) then
        rho_pert_inc = "rho_r2_pert_inc"
      else
        rho_pert_inc = "rho_pert_inc"
      end if

      call setup_field( iau_pert_fields, depository, prognostic_fields, &
          trim(rho_pert_inc), W3, mesh, checkpoint_restart_flag, &
          read_behaviour=read_behaviour )

      call setup_field( iau_pert_fields, depository, prognostic_fields, &
          "theta_pert_inc", Wtheta, mesh, checkpoint_restart_flag, &
          read_behaviour=read_behaviour )

      call setup_field( iau_pert_fields, depository, prognostic_fields, &
          "exner_pert_inc", W3, mesh, checkpoint_restart_flag, &
          read_behaviour=read_behaviour )
    end if ! ( iau_use_pertinc )
#endif

    if (iau_surf) then
      ! IAU land surface increment fields
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
    end if

    ! IAU sst increment fields
    if ( iau_sst ) then
      call log_event( 'Create IAU sst fields', LOG_LEVEL_INFO )

      call modeldb%fields%add_empty_field_collection("iau_sst_fields")
      iau_sst_fields => modeldb%fields%get_field_collection("iau_sst_fields")

      checkpoint_restart_flag = .false.
      read_behaviour => read_field_generic

      call setup_field( iau_sst_fields, depository, prognostic_fields, &
         "tstar_sea_inc", W3, mesh, checkpoint_restart_flag, &
         twod_mesh, read_behaviour=read_behaviour, twod=.true. )

    end if ! ( iau_sst )

    nullify( iau_fields, iau_surf_fields )

#ifdef UM_PHYSICS
    ! Set up name of density increments
    if ( iau_wet_density ) then
      rho_inc_name = "rho_r2_inc"
    else
      rho_inc_name = "rho_inc"
    end if

    ! Create field collection and fields to hold additive inflation increments
    if ( iau_use_addinf ) then
      call log_event( 'Create IAU additive inflation fields', LOG_LEVEL_INFO )
      call create_iau_additional_fields ( mesh, modeldb, "iau_addinf_fields", &
                                          "u_in_w3_inc", "v_in_w3_inc",       &
                                          "q_inc", "qcl_inc", "qcf_inc",      &
                                          "theta_inc", "exner_inc", rho_inc_name )
    end if

    ! Create field collection and fields to hold bias correction increments
    if ( iau_use_bcorr ) then
      call log_event( 'Create IAU bias correction fields', LOG_LEVEL_INFO )
      call create_iau_additional_fields ( mesh, modeldb, "iau_bcorr_fields", &
                                          "u_in_w3_inc", "v_in_w3_inc",      &
                                          "q_inc", "qcl_inc", "qcf_inc",     &
                                           "theta_inc", "exner_inc", rho_inc_name )
    end if

    ! Create field collection and fields to hold weighted aggregated increments
    call log_event( 'Create IAU total weighted increment fields', LOG_LEVEL_INFO )
    if ( iau_wet_density ) then
      rho_inc_name = "rho_r2_tot_inc"
    else
      rho_inc_name = "rho_tot_inc"
    end if
    call create_iau_additional_fields ( mesh, modeldb, "iau_tot_inc",         &
                                        "u_in_w3_tot_inc", "v_in_w3_tot_inc", &
                                        "q_tot_inc", "qcl_tot_inc",           &
                                        "qcf_tot_inc", "theta_tot_inc",       &
                                        "exner_tot_inc", rho_inc_name )
#endif

   end subroutine create_iau_fields

#ifdef UM_PHYSICS
  !> @brief   Create additional types of iau increment fields.
  !> @details Create IAU field collection for additional increment
  !>          types associated with additive inflation and the aggregated
  !>          total increment to be added to the prognostic state on any
  !>          given timestep.
  !> @param[in]     mesh              The current 3d mesh identifier
  !> @param[in,out] modeldb           The model database
  !> @param[in]     iau_incs          Type of increment to create fields for
  !> @param[in]     u_in_w3_inc_name  Name of u_in_w3 increment field
  !> @param[in]     v_in_w3_inc_name  Name of v_in_w3 increment field
  !> @param[in]     q_inc_name        Name of q increment field
  !> @param[in]     qcl_inc_name      Name of qcl increment field
  !> @param[in]     qcf_inc_name      Name of qcf increment field
  !> @param[in]     theta_inc_name    Name of theta increment field
  !> @param[in]     exner_inc_name    Name of exner increment field
  !> @param[in]     rho_inc_name      Name of rho increment field
  subroutine create_iau_additional_fields( mesh, modeldb, iau_incs, &
                                           u_in_w3_inc_name, &
                                           v_in_w3_inc_name, &
                                           q_inc_name, &
                                           qcl_inc_name, &
                                           qcf_inc_name, &
                                           theta_inc_name, &
                                           exner_inc_name, &
                                           rho_inc_name )

    implicit none

    type(mesh_type),    intent(in), pointer :: mesh
    type(modeldb_type), intent(inout) :: modeldb
    character(*),       intent(in)    :: iau_incs
    character(*),       intent(in)    :: u_in_w3_inc_name
    character(*),       intent(in)    :: v_in_w3_inc_name
    character(*),       intent(in)    :: q_inc_name
    character(*),       intent(in)    :: qcl_inc_name
    character(*),       intent(in)    :: qcf_inc_name
    character(*),       intent(in)    :: theta_inc_name
    character(*),       intent(in)    :: exner_inc_name
    character(*),       intent(in)    :: rho_inc_name

    type(field_collection_type), pointer :: iau_inc_fields
    type(function_space_type),   pointer :: w3_fs
    type(function_space_type),   pointer :: wtheta_fs

    type(field_type) :: u_in_w3_inc
    type(field_type) :: v_in_w3_inc
    type(field_type) :: q_inc
    type(field_type) :: qcl_inc
    type(field_type) :: qcf_inc
    type(field_type) :: theta_inc
    type(field_type) :: exner_inc
    type(field_type) :: rho_inc

    procedure(read_interface), pointer :: tmp_read_ptr

    call modeldb%fields%add_empty_field_collection(iau_incs)
    iau_inc_fields => modeldb%fields%get_field_collection(iau_incs)

    tmp_read_ptr => read_field_generic
    w3_fs => function_space_collection%get_fs( mesh,            &
                                               element_order_h, &
                                               element_order_v, &
                                               W3 )

    wtheta_fs => function_space_collection%get_fs( mesh,        &
                                               element_order_h, &
                                               element_order_v, &
                                               Wtheta )

    ! u inc
    call u_in_w3_inc%initialise( vector_space = w3_fs, &
                                 name=trim(u_in_w3_inc_name))
    call u_in_w3_inc%set_read_behaviour(tmp_read_ptr)
    call iau_inc_fields%add_field(u_in_w3_inc)

    ! v inc
    call v_in_w3_inc%initialise( vector_space = w3_fs, &
                                 name=trim(v_in_w3_inc_name))
    call v_in_w3_inc%set_read_behaviour(tmp_read_ptr)
    call iau_inc_fields%add_field(v_in_w3_inc)

    ! q inc
    call q_inc%initialise( vector_space = wtheta_fs, name=trim(q_inc_name) )
    call q_inc%set_read_behaviour(tmp_read_ptr)
    call iau_inc_fields%add_field(q_inc)

    ! qcl inc
    call qcl_inc%initialise( vector_space = wtheta_fs, name=trim(qcl_inc_name) )
    call qcl_inc%set_read_behaviour(tmp_read_ptr)
    call iau_inc_fields%add_field(qcl_inc)

    ! qcf inc
    call qcf_inc%initialise( vector_space = wtheta_fs, name=trim(qcf_inc_name) )
    call qcf_inc%set_read_behaviour(tmp_read_ptr)
    call iau_inc_fields%add_field(qcf_inc)

    ! theta inc
    call theta_inc%initialise( vector_space = wtheta_fs, &
                               name=trim(theta_inc_name) )
    call theta_inc%set_read_behaviour(tmp_read_ptr)
    call iau_inc_fields%add_field(theta_inc)

    ! exner inc
    call exner_inc%initialise( vector_space = w3_fs, name=trim(exner_inc_name) )
    call exner_inc%set_read_behaviour(tmp_read_ptr)
    call iau_inc_fields%add_field(exner_inc)

    ! rho inc
    call rho_inc%initialise( vector_space = w3_fs, name=trim(rho_inc_name) )
    call rho_inc%set_read_behaviour(tmp_read_ptr)
    call iau_inc_fields%add_field(rho_inc)

   end subroutine create_iau_additional_fields
#endif

end module create_iau_fields_mod
