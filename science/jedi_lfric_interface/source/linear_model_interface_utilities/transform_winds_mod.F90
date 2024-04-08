!-----------------------------------------------------------------------------
! (C) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!>@brief  Wind transform code used to interpolate between vector and scalar
!>        winds
!>
!>@detail This code will be replaced by code in #3734 when it is ready.
!>
module transform_winds_mod

  use constants_mod,                  only: i_def
  use field_collection_mod,           only: field_collection_type
  use field_mod,                      only: field_type
  use fs_continuity_mod,              only: W2H, W3
  use function_space_mod,             only: function_space_type
  use function_space_collection_mod,  only: function_space_collection
  use map_fd_to_prognostics_alg_mod,  only: set_wind
  use mesh_mod,                       only: mesh_type
  use physics_mappings_alg_mod,       only: map_physics_winds, &
                                            split_wind_alg

  implicit none

  private
  public :: wind_scalar_to_vector, wind_vector_to_scalar

  contains

  !> @brief   Interpolate the cell-centre scalar winds to the edge based vector wind
  !>
  !> @param[in,out] fields  A field collection that contains the wind fields
  subroutine wind_scalar_to_vector( fields )

    implicit none

    type( field_collection_type ), intent(inout) :: fields

    ! Local
    type( field_type ), pointer :: u_in_w3
    type( field_type ), pointer :: v_in_w3
    type( field_type ), pointer :: w_in_wth
    type( field_type ), pointer :: vector_wind

    nullify( u_in_w3, v_in_w3, w_in_wth, vector_wind )

    ! Get the fields
    call fields%get_field('u_in_w3', u_in_w3)
    call fields%get_field('v_in_w3', v_in_w3)
    call fields%get_field('w_in_wth', w_in_wth)
    call fields%get_field('u', vector_wind)

    ! Interpolate cell centred zonal/meridional winds
    call set_wind( vector_wind, u_in_w3, v_in_w3, w_in_wth )

  end subroutine wind_scalar_to_vector

  !> @brief  Interpolate the edge based vector wind to the cell centre scalar winds
  !>
  !> @param[in,out] fields  A field collection that contains the wind fields
  subroutine wind_vector_to_scalar( fields )

    implicit none

    type( field_collection_type ), intent(inout) :: fields

    ! Local
    type( function_space_type ), pointer :: fs
    type( mesh_type ),           pointer :: mesh
    type( field_type ),          pointer :: u_in_w3
    type( field_type ),          pointer :: v_in_w3
    type( field_type ),          pointer :: w_in_wth
    type( field_type ),          pointer :: vector_wind
    type( field_type )                   :: dummy_w2h
    type( field_type )                   :: dummy_w3
    integer( kind=i_def ),     parameter :: element_order = 0

    nullify( fs, mesh, u_in_w3, v_in_w3, w_in_wth, vector_wind )

    ! Get the fields from the prognostic_fields collection
    call fields%get_field('u_in_w3', u_in_w3)
    call fields%get_field('v_in_w3', v_in_w3)
    call fields%get_field('w_in_wth', w_in_wth)
    call fields%get_field('u', vector_wind)

    ! Setup dummy fields
    mesh => vector_wind%get_mesh()

    ! Get W3 function space from the field
    fs => function_space_collection%get_fs(mesh, element_order, W3)
    call dummy_w3%initialise(fs, "dummy_w3")

    ! Get W2 function space from the field
    fs => function_space_collection%get_fs(mesh, element_order, W2H)
    call dummy_w2h%initialise(fs, "dummy_w2h")

    ! Get horizontal
    call map_physics_winds( u_in_w3, v_in_w3, dummy_w3, vector_wind )

    ! Get vertical
    call split_wind_alg( dummy_w2h, dummy_w2h, w_in_wth, vector_wind, mesh )

  end subroutine wind_vector_to_scalar

end module transform_winds_mod
