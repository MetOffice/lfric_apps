!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Outputs diagnostics from linear model

module linear_diagnostics_driver_mod

  use constants_mod,             only : i_def, str_def, l_def
  use diagnostic_alg_mod,        only : column_total_diagnostics_alg,          &
                                        calc_wbig_diagnostic_alg,              &
                                        pressure_diag_alg
  use diagnostics_io_mod,        only : write_scalar_diagnostic, &
                                        write_vector_diagnostic
  use field_collection_mod,      only : field_collection_type
  use field_collection_iterator_mod, &
                                 only : field_collection_iterator_type
  use driver_modeldb_mod,        only : modeldb_type
  use field_array_mod,           only : field_array_type
  use field_mod,                 only : field_type
  use field_parent_mod,          only : field_parent_type, write_interface
  use formulation_config_mod,    only : use_physics,             &
                                        moisture_formulation,    &
                                        moisture_formulation_dry
  use fs_continuity_mod,         only : W3, Wtheta
  use mesh_mod,                  only : mesh_type
  use mr_indices_mod,            only : nummr, mr_names
  use initialization_config_mod, only : ls_option, &
                                        ls_option_file
  use initialise_diagnostics_mod, &
                                 only : diagnostic_to_be_sampled
  use io_value_mod,              only : io_value_type, get_io_value
  use io_config_mod,             only : use_xios_io, write_fluxes, &
                                        write_diag
  use log_mod,                   only : log_event, &
                                        LOG_LEVEL_INFO
  use linear_config_mod,         only : ls_read_w2h
  use lfric_xios_write_mod,      only : write_field_generic
  use sci_geometric_constants_mod,                                             &
                                 only : get_panel_id

  implicit none

  private
  public linear_write_initial_output, &
         linear_diagnostics_driver

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
  subroutine linear_write_initial_output( modeldb, mesh, twod_mesh, &
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
      call linear_diagnostics_driver( modeldb,   &
                                      mesh,      &
                                      twod_mesh, &
                                      nodal_output_on_w3 )

    end if

  end subroutine linear_write_initial_output

  !> @brief Outputs simple diagnostics from linear model
  !> @param[in,out] modeldb             Working data set of model run.
  !> @param[in]     mesh                The primary mesh
  !> @param[in]     twod_mesh           The 2d mesh
  !> @param[in]     nodal_output_on_w3  Flag that determines if vector fields
  !>                should be projected to W3 for nodal output
  subroutine linear_diagnostics_driver( modeldb,   &
                                        mesh,      &
                                        twod_mesh, &
                                        nodal_output_on_w3 )

    implicit none

    type(mesh_type),      intent(in), pointer :: mesh
    type(mesh_type),      intent(in), pointer :: twod_mesh
    type(modeldb_type),   intent(in), target  :: modeldb
    logical,              intent(in)          :: nodal_output_on_w3

    type(field_collection_type), pointer :: prognostic_fields
    type(field_collection_type), pointer :: moisture_fields => null()
    type(field_type),            pointer :: mr(:)
    type(field_collection_type), pointer :: derived_fields
    type(field_array_type),      pointer :: ls_mr_array => null()
    type(field_collection_type), pointer :: ls_fields
    type(field_type),            pointer :: ls_mr(:) => null()

    type(field_type), pointer :: theta
    type(field_type), pointer :: u
    type(field_type), pointer :: rho
    type(field_type), pointer :: exner
    type(field_type), pointer :: panel_id
    type(field_type), pointer :: u_in_w2h
    type(field_type), pointer :: v_in_w2h
    type(field_type), pointer :: w_in_wth
    type(field_type), pointer :: exner_in_wth

    type(field_array_type), pointer :: mr_array

    type(field_type), pointer :: ls_theta => null()
    type(field_type), pointer :: ls_u => null()
    type(field_type), pointer :: ls_rho => null()
    type(field_type), pointer :: ls_exner => null()
    type(field_type), pointer :: ls_v_u => null()
    type(field_type), pointer :: ls_h_u => null()

    ! Iterator for field collection
    type(field_collection_iterator_type)  :: iterator

    ! A pointer used for retrieving fields from collections
    ! when iterating over them
    class(field_parent_type),   pointer :: field_ptr
    procedure(write_interface), pointer :: tmp_write_ptr
    type(io_value_type),        pointer :: temp_corr_io_value

    integer :: i, fs
    character(len=str_def) :: name

    call log_event("Linear: writing diagnostic output", LOG_LEVEL_INFO)

    prognostic_fields => modeldb%fields%get_field_collection("prognostic_fields")
    ls_fields => modeldb%fields%get_field_collection("ls_fields")
    moisture_fields => modeldb%fields%get_field_collection("moisture_fields")
    call moisture_fields%get_field("mr", mr_array)
    mr => mr_array%bundle
    call moisture_fields%get_field("ls_mr", ls_mr_array)
    ls_mr => ls_mr_array%bundle
    derived_fields => modeldb%fields%get_field_collection("derived_fields")
    panel_id => get_panel_id(mesh%get_id())

    call prognostic_fields%get_field('theta', theta)
    call prognostic_fields%get_field('u', u)
    call prognostic_fields%get_field('rho', rho)
    call prognostic_fields%get_field('exner', exner)

    call ls_fields%get_field('ls_theta', ls_theta)
    call ls_fields%get_field('ls_u', ls_u)
    call ls_fields%get_field('ls_rho', ls_rho)
    call ls_fields%get_field('ls_exner', ls_exner)

    ! Scalar fields
    call write_scalar_diagnostic('rho', rho, &
                                 modeldb%clock, mesh, nodal_output_on_w3)
    call write_scalar_diagnostic('theta', theta, &
                                 modeldb%clock, mesh, nodal_output_on_w3)
    call write_scalar_diagnostic('exner', exner, &
                                 modeldb%clock, mesh, nodal_output_on_w3)
    call write_scalar_diagnostic('ls_rho', ls_rho, &
                                 modeldb%clock, mesh, nodal_output_on_w3)
    call write_scalar_diagnostic('ls_theta', ls_theta, &
                                 modeldb%clock, mesh, nodal_output_on_w3)
    call write_scalar_diagnostic('ls_exner', ls_exner, &
                                 modeldb%clock, mesh, nodal_output_on_w3)

    ! Vector fields
    if (use_physics .and. use_xios_io .and. .not. write_fluxes) then
      ! These have already been calculated, so no need to recalculate them
      call derived_fields%get_field('u_in_w2h', u_in_w2h)
      call derived_fields%get_field('v_in_w2h', v_in_w2h)
      call derived_fields%get_field('w_in_wth', w_in_wth)
      tmp_write_ptr => write_field_generic
      call u_in_w2h%set_write_behaviour(tmp_write_ptr)
      call v_in_w2h%set_write_behaviour(tmp_write_ptr)
      if (modeldb%clock%is_initialisation()) then
        if (diagnostic_to_be_sampled("init_u_in_w2h")) then
          call u_in_w2h%write_field("init_u_in_w2h")
        end if
        if (diagnostic_to_be_sampled("init_v_in_w2h")) then
          call v_in_w2h%write_field("init_v_in_w2h")
        end if
        if (diagnostic_to_be_sampled("init_w_in_wth")) then
          call w_in_wth%write_field("init_w_in_wth")
        end if
      else
        if (diagnostic_to_be_sampled("u_in_w2h")) then
          call u_in_w2h%write_field("u_in_w2h")
        end if
        if (diagnostic_to_be_sampled("v_in_w2h")) then
          call v_in_w2h%write_field("v_in_w2h")
        end if
        if (diagnostic_to_be_sampled("w_in_wth")) then
          call w_in_wth%write_field("w_in_wth")
        end if
      end if
    else
      call write_vector_diagnostic('u', u, &
                                   modeldb%clock, mesh, nodal_output_on_w3)
    end if

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
        call write_scalar_diagnostic( trim(mr_names(i)), mr(i),           &
                                      modeldb%clock, mesh, nodal_output_on_w3 )
      end do
    end if

    ! Derived physics fields (only those on W3 or Wtheta)
    if (use_physics .and. use_xios_io .and. .not. modeldb%clock%is_initialisation()) then
      field_ptr => null()
      call iterator%initialise(derived_fields)
      do
        if ( .not.iterator%has_next() ) exit
        field_ptr => iterator%next()
        select type(field_ptr)
          type is (field_type)
            fs = field_ptr%which_function_space()
            if ( fs == W3 .or. fs == Wtheta ) then
              name = trim(adjustl( field_ptr%get_name() ))
              if (diagnostic_to_be_sampled(trim(name))) &
                   call field_ptr%write_field(trim(name))
            end if
        end select
      end do
      field_ptr => null()

      ! Get w_in_wth for WBig calculation
      call derived_fields%get_field('w_in_wth', w_in_wth)
      call calc_wbig_diagnostic_alg(w_in_wth, mesh)

      ! Pressure diagnostics
      call prognostic_fields%get_field('exner', exner)
      call pressure_diag_alg(exner)

      call derived_fields%get_field('exner_in_wth', exner_in_wth)
      call pressure_diag_alg(exner_in_wth)

      temp_corr_io_value => get_io_value( modeldb%values, 'temperature_correction_io_value')
      call column_total_diagnostics_alg(rho, mr, derived_fields, exner, &
                                        mesh, twod_mesh,                &
                                        temp_corr_io_value%data(1))

    end if

  end subroutine linear_diagnostics_driver

end module linear_diagnostics_driver_mod
