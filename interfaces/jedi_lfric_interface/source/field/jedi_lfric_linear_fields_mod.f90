!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!> @brief A module providing the linear variables and a method to create them.
!>
!> @details This module provides a list of the linear model prognostic variable
!>          names and corresponding function space enumerators to be used with
!>          JEDI-LFRIC. In addition, a method is provided that returns a field
!>          collection that contains the set of fields defined in the list. The
!>          list differs from the prognostic set of variables in the linear
!>          model in the following ways:
!>            1. Cell centred winds (W3 of U/V and Wtheta for W) are included
!>               instead for the W2 wind field. This is because the Atlas fields
!>               are not stored on the W2 function space.
!>            2. Some moist fields are omitted from the list because in the
!>               linear model they are computed analytically including:
!>               i)  mixing ratios, namely "m_g" and "m_ci".
!>               ii) all three "moist_dyn" fields.
!>
module jedi_lfric_linear_fields_mod
  use constants_mod,                 only : i_def, str_def, r_def, l_def
  use field_mod,                     only : field_type
  use field_collection_mod,          only : field_collection_type
  use fs_continuity_mod,             only : W3, Wtheta, W2
  use function_space_collection_mod, only : function_space_collection
  use mesh_mod,                      only : mesh_type

  implicit none

  private

  integer( kind=i_def ), parameter :: element_order_h = 0
  integer( kind=i_def ), parameter :: element_order_v = 0
  integer( kind=i_def ), parameter :: nvars = 10
  integer( kind=i_def ), parameter :: ls_nvars = 11
  character( len=str_def ), parameter, public :: &
                                     variable_names(nvars) =  (/'theta   ', &
                                                                'exner   ', &
                                                                'rho     ', &
                                                                'u_in_w3 ', &
                                                                'v_in_w3 ', &
                                                                'w_in_wth', &
                                                                'm_v     ', &
                                                                'm_cl    ', &
                                                                'm_r     ', &
                                                                'm_s     '/)
  character( len=str_def ), parameter, public :: &
                               ls_variable_names(ls_nvars) =  (/'theta        ', &
                                                                'exner        ', &
                                                                'rho          ', &
                                                                'u_in_w3      ', &
                                                                'v_in_w3      ', &
                                                                'w_in_wth     ', &
                                                                'm_v          ', &
                                                                'm_cl         ', &
                                                                'm_r          ', &
                                                                'm_s          ', &
                                                                'land_fraction'/)

  integer( kind=i_def ), parameter, public :: &
                            variable_function_spaces(nvars) = (/Wtheta, &
                                                                W3,     &
                                                                W3,     &
                                                                W3,     &
                                                                W3,     &
                                                                Wtheta, &
                                                                Wtheta, &
                                                                Wtheta, &
                                                                Wtheta, &
                                                                Wtheta/)
  integer( kind=i_def ), parameter, public :: &
                      ls_variable_function_spaces(ls_nvars) = (/Wtheta, &
                                                                W3,     &
                                                                W3,     &
                                                                W3,     &
                                                                W3,     &
                                                                Wtheta, &
                                                                Wtheta, &
                                                                Wtheta, &
                                                                Wtheta, &
                                                                Wtheta, &
                                                                W3/)

  logical( kind=l_def ), parameter, public :: &
                                ls_variable_is_2d(ls_nvars) = (/.false., &
                                                                .false., &
                                                                .false., &
                                                                .false., &
                                                                .false., &
                                                                .false., &
                                                                .false., &
                                                                .false., &
                                                                .false., &
                                                                .false., &
                                                                .true./)

  public :: create_linear_fields

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
contains

!> @brief    Create a field collection that includes the linear model variables
!>
!> @param [in]  mesh           Pointer to a mesh object
!> @param [in]  twod_mesh      Pointer to a 2D mesh object
!> @param [out] linear_fields  A field collection that includes the linear
!>                             fields
subroutine create_linear_fields( mesh, twod_mesh, linear_fields )

  implicit none

  type( mesh_type ), pointer,     intent(in) :: mesh
  type( mesh_type ), pointer,     intent(in) :: twod_mesh
  type( field_collection_type ), intent(out) :: linear_fields

  ! Local
  type( field_type )              :: field
  type( mesh_type ), pointer      :: mesh_for_field
  character( len=str_def )        :: variable_name
  integer                         :: i

  ! Setup the field_collection
  call linear_fields%initialise(name = 'linear_state_trajectory', table_len = ls_nvars)

  ! Create and add the fields defined in the list of variable names
  do i = 1, ls_nvars

    variable_name = trim(ls_variable_names(i))

    if (ls_variable_is_2d(i)) then
      mesh_for_field => twod_mesh
    else
      mesh_for_field => mesh
    end if

    call field%initialise( &
           vector_space = function_space_collection%get_fs(mesh_for_field,     &
                                                           element_order_h,    &
                                                           element_order_v,    &
                                                           ls_variable_function_spaces(i)), &
           name = variable_name )

    call linear_fields%add_field( field )

  end do

end subroutine create_linear_fields

end module jedi_lfric_linear_fields_mod
