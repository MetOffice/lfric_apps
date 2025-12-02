!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Performs the intial conditions output for Gungho
!>
module initial_output_mod

  use constants_mod,              only : i_def, l_def
  use io_config_mod,              only : write_diag, use_xios_io
  use io_context_mod,             only : io_context_type
  use gungho_diagnostics_driver_mod, &
                                  only : gungho_diagnostics_driver
  use driver_modeldb_mod,         only : modeldb_type
  use log_mod,                    only : log_event, log_level_error
  use mesh_mod,                   only : mesh_type
  use io_context_mod,             only : io_context_type
  use limited_area_constants_mod, only : write_masks
  use boundaries_config_mod,      only : limited_area, lbc_method, &
                                         lbc_method_onion_layer
  implicit none

  private
  public :: write_initial_output

contains

  !> @brief Outputs simple diagnostics from Gungho/LFRic
  !>
  !> @param[in,out] modeldb             Working data set of model.
  !> @param[in]     mesh                The primary mesh
  !> @param[in]     twod_mesh           The 2d mesh
  !> @param[in]     io_context_name     The name of the IO context for writing
  !> @param[in]     nodal_output_on_w3  Flag that determines if vector fields
  !>                  should be projected to W3 for nodal output
  !>
  subroutine write_initial_output( modeldb, mesh, twod_mesh, &
                                   io_context_name, nodal_output_on_w3 )

#ifdef USE_XIOS
    !>@todo This will be removed by #3321
    use lfric_xios_context_mod, only : lfric_xios_context_type
    use lfric_xios_action_mod,  only : advance
#endif

    implicit none

    class(modeldb_type),      intent(inout) :: modeldb
    type(mesh_type), pointer, intent(in)    :: mesh
    type(mesh_type), pointer, intent(in)    :: twod_mesh
    character(*),             intent(in)    :: io_context_name
    logical(l_def),           intent(in)    :: nodal_output_on_w3
#ifdef USE_XIOS
    type(lfric_xios_context_type), pointer :: io_context
    ! Call clock initial step before initial conditions output
    if (modeldb%clock%is_initialisation() .and. use_xios_io) then
      call modeldb%io_contexts%get_io_context(io_context_name, io_context)
      call advance(io_context, modeldb%clock)
    end if
#endif

    if (modeldb%clock%is_initialisation() .and. write_diag) then
      ! Calculation and output of initial conditions
      call gungho_diagnostics_driver( modeldb,   &
                                      mesh,      &
                                      twod_mesh, &
                                      nodal_output_on_w3 )

      if ( limited_area )then
        if ( lbc_method == lbc_method_onion_layer ) call write_masks()
      end if
    end if

  end subroutine write_initial_output

end module initial_output_mod
