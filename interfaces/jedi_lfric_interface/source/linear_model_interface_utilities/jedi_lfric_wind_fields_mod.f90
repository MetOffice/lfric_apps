!------------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!> @brief   Contains routines to deal with the scalar and vector winds.
!>
!> @detail  Contains helper routines for use with the LFRIC-JEDI and
!>          jedi_lfric_tests emulator interfaces. These routines are required
!>          because the aforementioned model interfaces store only cell-centred
!>          Atlas fields.
!>
module jedi_lfric_wind_fields_mod

  use constants_mod,                 only : str_def, i_def
  use field_collection_mod,          only : field_collection_type
  use field_mod,                     only : field_type
  use fs_continuity_mod,             only : W2, W3, Wtheta
  use function_space_collection_mod, only : function_space_collection
  use function_space_mod,            only : function_space_type
  use interpolation_alg_mod,         only : interp_w3wth_to_w2_alg
  use log_mod,                       only : log_event,         &
                                            log_scratch_space, &
                                            LOG_LEVEL_ERROR,   &
                                            LOG_LEVEL_DEBUG
  use mesh_collection_mod,           only : mesh_collection
  use mesh_mod,                      only : mesh_type
  use namelist_mod,                  only : namelist_type
  use pure_abstract_field_mod,       only : pure_abstract_field_type

  implicit none

  private
  public :: create_scalar_winds, setup_vector_wind

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Create scalar wind fields in the modeldb depository and prognostics
  !>        field collections
  !>
  !> @detail Creates and adds scalar winds to the modeldb depository and
  !>         prognostic field collections. This provides cell-centred fields
  !>         that the model interfaces (LFRIC-JEDI, Jedi emulator) can link and
  !>         copy to/from.
  !>
  !> @param[inout] modeldb    The modeldb that contains the model_data
  !> @param[in]    mesh       The current 3d mesh
  subroutine create_scalar_winds( modeldb, mesh )

    use driver_modeldb_mod, only : modeldb_type

    implicit none

    type(modeldb_type), target, intent(inout) :: modeldb

    ! Local
    type( mesh_type ),                 pointer :: mesh
    type( field_type ),                pointer :: field_ptr
    class( pure_abstract_field_type ), pointer :: tmp_ptr
    type( field_collection_type ),     pointer :: depository
    type( field_collection_type ),     pointer :: prognostic_fields
    type( function_space_type ),       pointer :: w3_fs
    type( function_space_type ),       pointer :: wtheta_fs
    type( namelist_type ),             pointer :: base_mesh_nml
    character( len=str_def )                   :: prime_mesh_name
    integer( kind=i_def ),           parameter :: element_order_h = 0_i_def
    integer( kind=i_def ),           parameter :: element_order_v = 0_i_def

    ! Temporary fields to create prognostics
    type( field_type ) :: u_in_w3
    type( field_type ) :: v_in_w3
    type( field_type ) :: w_in_wth

    call log_event( 'jedi-lfric: Creating and adding scalar winds', &
                    LOG_LEVEL_DEBUG )

    nullify( mesh, field_ptr, tmp_ptr, depository, prognostic_fields )
    nullify( wtheta_fs, w3_fs, base_mesh_nml )

    depository => modeldb%fields%get_field_collection("depository")
    prognostic_fields => modeldb%fields%get_field_collection("prognostic_fields")

    ! Get the mesh
    base_mesh_nml => modeldb%configuration%get_namelist('base_mesh')
    call base_mesh_nml%get_value( 'prime_mesh_name', prime_mesh_name )
    mesh => mesh_collection%get_mesh(prime_mesh_name)

    ! Create prognostic fields
    w3_fs     => function_space_collection%get_fs(mesh, element_order_h,       &
                                                  element_order_v, W3)
    wtheta_fs => function_space_collection%get_fs(mesh, element_order_h,       &
                                                  element_order_v, Wtheta)

    ! Populate the depository if fields are not present
    if (.not.depository%field_exists('u_in_w3')) then
      call u_in_w3%initialise( vector_space=w3_fs, name="u_in_w3" )
      call depository%add_field( u_in_w3 )
    endif
    if (.not.depository%field_exists('v_in_w3')) then
      call v_in_w3%initialise( vector_space=w3_fs, name="v_in_w3" )
      call depository%add_field( v_in_w3 )
    endif
    if (.not.depository%field_exists('w_in_wth')) then
      call w_in_wth%initialise( vector_space=wtheta_fs, name="w_in_wth" )
      call depository%add_field( w_in_wth )
    endif

    ! Populate the prognostic field collection
    call depository%get_field('u_in_w3', field_ptr)
    tmp_ptr => field_ptr
    call prognostic_fields%add_reference_to_field(tmp_ptr)

    call depository%get_field('v_in_w3', field_ptr)
    tmp_ptr => field_ptr
    call prognostic_fields%add_reference_to_field(tmp_ptr)

    call depository%get_field('w_in_wth', field_ptr)
    tmp_ptr => field_ptr
    call prognostic_fields%add_reference_to_field(tmp_ptr)

  end subroutine create_scalar_winds

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Setup a W2 vector wind field in the linear state
  !>
  !> @detail Create a W2 vector wind and add it to the linear state. The vector
  !>         wind is populated by interpolating the scalar winds fields. After
  !>         the W2 field is created and populated, the scaler winds are
  !>         removed as they are no longer required.
  !>
  !> @param [in, out] linear_state a field collection the includes scalar winds
  subroutine setup_vector_wind( linear_state )

    implicit none

    type(field_collection_type), intent(inout) :: linear_state

    ! Local
    type(field_type),          pointer :: u_in_w3
    type(field_type),          pointer :: v_in_w3
    type(field_type),          pointer :: w_in_wth
    type(mesh_type),           pointer :: mesh
    type(function_space_type), pointer :: fs
    type(field_type),      allocatable :: vector_wind
    integer( kind=i_def ),   parameter :: element_order_h = 0_i_def
    integer( kind=i_def ),   parameter :: element_order_v = 0_i_def

    nullify( u_in_w3, v_in_w3, w_in_wth, mesh, fs )

    ! Check scaler fields are present
    if ( linear_state%field_exists("u_in_w3") .and. &
         linear_state%field_exists("v_in_w3") .and. &
         linear_state%field_exists("w_in_wth") ) then

      ! Get the scaler winds
      call linear_state%get_field("u_in_w3", u_in_w3)
      call linear_state%get_field("v_in_w3", v_in_w3)
      call linear_state%get_field("w_in_wth", w_in_wth)

      ! Create a vector wind
      mesh => u_in_w3%get_mesh()
      fs => function_space_collection%get_fs(mesh, element_order_h,            &
                                             element_order_v, W2)
      allocate(vector_wind)
      call vector_wind%initialise(fs, "u")

      ! Interpolate cell-centre to edge
      call interp_w3wth_to_w2_alg(vector_wind, u_in_w3, v_in_w3, w_in_wth)

      ! Remove the scalar winds
      call linear_state%remove_field("u_in_w3")
      call linear_state%remove_field("v_in_w3")
      call linear_state%remove_field("w_in_wth")

      ! Add vector wind to collection
      call linear_state%add_field( vector_wind )
      deallocate( vector_wind )

    else
      write(log_scratch_space, '(A)') &
        "The three scalar wind components are required but at least one is missing."
      call log_event( log_scratch_space, LOG_LEVEL_ERROR )

    end if

  end subroutine setup_vector_wind

end module jedi_lfric_wind_fields_mod
