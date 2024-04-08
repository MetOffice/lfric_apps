!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief init functionality for the miniapp skeleton

!> @details Handles init of prognostic fields and through the call to
!>          runtime_constants the coordinate fields and fem operators

module init_solver_miniapp_mod

  use constants_mod,                  only : i_def, r_def
  use field_mod,                      only : field_type
  use finite_element_config_mod,      only : element_order
  use function_space_collection_mod,  only : function_space_collection_type, &
                                             function_space_collection
  use fs_continuity_mod,              only : W0, W3
  use init_solver_fields_alg_mod,     only : init_solver_fields_alg
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO, &
                                             LOG_LEVEL_ERROR, &
                                             log_scratch_space
  use mesh_mod,                       only : mesh_type
  use sci_field_vector_mod,           only : field_vector_type

  implicit none


contains
  !> initialise the fields and field vector for the miniapp.
  !> @param[in]    mesh  The primary mesh
  subroutine init_solver_miniapp( mesh, fv )

    implicit none

    type( mesh_type ),  intent(in), pointer :: mesh
    ! prognostic fields
    type( field_vector_type ), intent(inout) :: fv
    type( field_type )                       :: f1, f2


    call log_event( 'solver miniapp: initialisation...', LOG_LEVEL_INFO )


    ! Create prognostic fields
    ! Create a field in the W0 function space (fully continuous field)
    call f1%initialise( vector_space = &
         function_space_collection%get_fs(mesh, element_order, W0) )
    ! Create a field in the W3 function space (fully discontinuous field)
    call f2%initialise( vector_space = &
         function_space_collection%get_fs(mesh, element_order, W3) )

    ! Initialise the fields
    call init_solver_fields_alg(f1, f2)

    ! make the vector
    fv = field_vector_type(2_i_def)
    call fv%import_field(f1,1)
    call fv%import_field(f2,2)

    write(log_scratch_space,'(A,E16.8)') "W0, W3 Fvector initialised to 0.5/1.0 norm=",fv%norm()
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    call log_event( 'solver miniapp initialised', LOG_LEVEL_INFO )

  end subroutine init_solver_miniapp

end module init_solver_miniapp_mod
