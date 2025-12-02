!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Outputs diagnostics from linear model

module linear_diagnostics_driver_mod

  use constants_mod,             only : i_def, str_def
  use diagnostics_io_mod,        only : write_scalar_diagnostic, &
                                        write_vector_diagnostic
  use field_collection_mod,      only : field_collection_type
  use driver_modeldb_mod,        only : modeldb_type
  use field_array_mod,           only : field_array_type
  use field_mod,                 only : field_type
  use formulation_config_mod,    only : moisture_formulation,    &
                                        moisture_formulation_dry
  use mesh_mod,                  only : mesh_type
  use mr_indices_mod,            only : nummr, mr_names
  use initialization_config_mod, only : ls_option, &
                                        ls_option_file
  use log_mod,                   only : log_event, &
                                        LOG_LEVEL_INFO
  use linear_config_mod,         only : ls_read_w2h

  implicit none

  private
  public linear_diagnostics_driver

contains

  !> @brief Outputs simple diagnostics from Linear model
  !> @param[in] mesh        The primary mesh
  !> @param[in] modeldb     The working data set for the model run
  !> @param[in] nodal_output_on_w3 Flag that determines if vector fields
  !>                  should be projected to W3 for nodal output
  subroutine linear_diagnostics_driver( mesh,    &
                                        modeldb, &
                                        nodal_output_on_w3 )

    implicit none

    type(mesh_type),      intent(in), pointer :: mesh
    type(modeldb_type),   intent(in), target  :: modeldb
    logical,              intent(in)          :: nodal_output_on_w3
    type(field_collection_type), pointer :: moisture_fields => null()
    type(field_array_type), pointer      :: ls_mr_array => null()

    type( field_collection_type ), pointer :: ls_fields
    type( field_type ),            pointer :: ls_mr(:) => null()

    type( field_type), pointer :: ls_theta => null()
    type( field_type), pointer :: ls_u => null()
    type( field_type), pointer :: ls_rho => null()
    type( field_type), pointer :: ls_exner => null()
    type( field_type), pointer :: ls_v_u => null()
    type( field_type), pointer :: ls_h_u => null()

    integer :: i

    call log_event("Linear: writing diagnostic output", LOG_LEVEL_INFO)

    ls_fields => modeldb%fields%get_field_collection("ls_fields")
    moisture_fields => modeldb%fields%get_field_collection("moisture_fields")
    call moisture_fields%get_field("ls_mr", ls_mr_array)
    ls_mr => ls_mr_array%bundle

    call ls_fields%get_field('ls_theta', ls_theta)
    call ls_fields%get_field('ls_u', ls_u)
    call ls_fields%get_field('ls_rho', ls_rho)
    call ls_fields%get_field('ls_exner', ls_exner)

    ! Scalar fields
    call write_scalar_diagnostic('ls_rho', ls_rho, &
                                 modeldb%clock, mesh, nodal_output_on_w3)
    call write_scalar_diagnostic('ls_theta', ls_theta, &
                                 modeldb%clock, mesh, nodal_output_on_w3)
    call write_scalar_diagnostic('ls_exner', ls_exner, &
                                 modeldb%clock, mesh, nodal_output_on_w3)

    ! Vector fields
    call write_vector_diagnostic('ls_u', ls_u, &
                                 modeldb%clock, mesh, nodal_output_on_w3)


    ! Fluxes - horizontal and vertical (if reading linearisation
    ! state from file)
    if (ls_option == ls_option_file) then
      if (ls_read_w2h) then
        call ls_fields%get_field('ls_v_u', ls_v_u)
        call ls_fields%get_field('ls_h_u', ls_h_u)
        call write_scalar_diagnostic('readls_v_u', ls_v_u, &
                                     modeldb%clock, mesh, nodal_output_on_w3)
        call write_vector_diagnostic('readls_h_u', ls_h_u, &
                                     modeldb%clock, mesh, nodal_output_on_w3)
      end if
    end if

    ! Moisture fields
    if (moisture_formulation /= moisture_formulation_dry) then
      do i=1,nummr
        call write_scalar_diagnostic( 'ls_'//trim(mr_names(i)), ls_mr(i), &
                                      modeldb%clock, mesh, nodal_output_on_w3 )
      end do
    end if

  end subroutine linear_diagnostics_driver

end module linear_diagnostics_driver_mod
