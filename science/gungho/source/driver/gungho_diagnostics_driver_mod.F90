!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Outputs diagnostics from gungho/lfric_atm

!> @details Calls the routine that generates diagnostic output for
!>          gungho/lfric_atm. This is only a temporary
!>          hard-coded solution in lieu of a proper dianostic system

module gungho_diagnostics_driver_mod

  use constants_mod,             only : i_def, r_def, str_def
  use boundaries_config_mod,     only : limited_area, output_lbcs
  use diagnostic_alg_mod,        only : column_total_diagnostics_alg, &
                                        calc_wbig_diagnostic_alg, &
                                        pressure_diag_alg
  use diagnostics_io_mod,        only : write_scalar_diagnostic, &
                                        write_vector_diagnostic
  use diagnostics_calc_mod,      only : write_divergence_diagnostic, &
                                        write_hydbal_diagnostic,     &
                                        write_vorticity_diagnostic,  &
                                        write_pv_diagnostic
  use initialise_diagnostics_mod, only : diagnostic_to_be_sampled
  use field_array_mod,           only : field_array_type
  use field_collection_iterator_mod, &
                                 only : field_collection_iterator_type
  use field_collection_mod,      only : field_collection_type
  use field_mod,                 only : field_type
  use field_parent_mod,          only : field_parent_type, write_interface
  use io_value_mod,              only : io_value_type, get_io_value
  use lfric_xios_write_mod,      only : write_field_generic
  use formulation_config_mod,    only : use_physics,             &
                                        moisture_formulation,    &
                                        moisture_formulation_dry
  use fs_continuity_mod,         only : W3, Wtheta
  use integer_field_mod,         only : integer_field_type
  use initialization_config_mod, only : ls_option,          &
                                        ls_option_analytic, &
                                        ls_option_file
  use mesh_mod,                  only : mesh_type
  use moist_dyn_mod,             only : num_moist_factors
  use mr_indices_mod,            only : nummr, mr_names
  use log_mod,                   only : log_event, &
                                        LOG_LEVEL_DEBUG
  use sci_geometric_constants_mod,      &
                                 only : get_panel_id, get_height_fe, &
                                        get_height_fv, get_da_msl_proj
  use io_config_mod,             only : use_xios_io, write_fluxes
  use timer_mod,                 only : timer
  use timing_mod,                only : start_timing, stop_timing, tik, LPROF
  use transport_config_mod,      only : transport_ageofair
  use driver_modeldb_mod,        only : modeldb_type

#ifdef UM_PHYSICS
  use pres_lev_diags_alg_mod,    only : pres_lev_diags_alg
  use pmsl_alg_mod,              only : pmsl_alg
  use rh_diag_alg_mod,           only : rh_diag_alg
  use freeze_lev_alg_mod,        only : freeze_lev_alg
#endif

  implicit none

  private
  public gungho_diagnostics_driver

contains

  !> @brief Outputs simple diagnostics from Gungho/LFRic
  !>
  !> @param[in,out] modeldb             Working data set of model run.
  !> @param[in]     mesh                The primary mesh
  !> @param[in]     twod_mesh           The 2d mesh
  !> @param[in]     nodal_output_on_w3  Flag that determines if vector fields
  !>                  should be projected to W3 for nodal output
  subroutine gungho_diagnostics_driver( modeldb,   &
                                        mesh,      &
                                        twod_mesh, &
                                        nodal_output_on_w3 )

    implicit none

    class(modeldb_type), intent(inout), target  :: modeldb
    type(mesh_type),     intent(in),    pointer :: mesh
    type(mesh_type),     intent(in),    pointer :: twod_mesh
    logical,             intent(in)             :: nodal_output_on_w3

    type( field_collection_type ), pointer :: prognostic_fields => null()
    type( field_collection_type ), pointer :: con_tracer_last_outer
    type( field_collection_type ), pointer :: lbc_fields
    type( field_collection_type ), pointer :: moisture_fields => null()
    type( field_type ),            pointer :: mr(:) => null()
    type( field_type ),            pointer :: moist_dyn(:) => null()
    type( field_collection_type ), pointer :: derived_fields

    type( field_type), pointer :: theta => null()
    type( field_type), pointer :: u => null()
    type( field_type), pointer :: h_u => null()
    type( field_type), pointer :: v_u => null()
    type( field_type), pointer :: rho => null()
    type( field_type), pointer :: exner => null()
    type( field_type), pointer :: panel_id => null()
    type( field_type), pointer :: height_w3 => null()
    type( field_type), pointer :: height_wth => null()
    type( field_type), pointer :: lbc_u => null()
    type( field_type), pointer :: lbc_theta => null()
    type( field_type), pointer :: lbc_rho => null()
    type( field_type), pointer :: lbc_exner => null()
    type( field_type), pointer :: lbc_m_v=> null()
    type( field_type), pointer :: lbc_q=> null()
    type( field_type), pointer :: u_in_w2h => null()
    type( field_type), pointer :: v_in_w2h => null()
    type( field_type), pointer :: w_in_wth => null()
    type( field_type), pointer :: ageofair => null()
    type( field_type), pointer :: exner_in_wth => null()
    type( field_type), pointer :: dA => null()

    type(field_array_type), pointer :: mr_array => null()
    type(field_array_type), pointer :: moist_dyn_array => null()

    ! Iterator for field collection
    type(field_collection_iterator_type)  :: iterator

    ! A pointer used for retrieving fields from collections
    ! when iterating over them
    class( field_parent_type ), pointer :: field_ptr  => null()

    type(io_value_type), pointer :: temp_corr_io_value

    character(str_def) :: name

    integer :: fs
    integer :: element_order_h, element_order_v

    procedure(write_interface), pointer  :: tmp_write_ptr => null()

    integer :: i
    integer(tik)  :: id
    if ( LPROF ) call start_timing( id, 'gungho_diagnostics_driver' )
    call log_event("Gungho: writing diagnostic output", LOG_LEVEL_DEBUG)

    ! Get pointers to field collections for use downstream
    prognostic_fields => modeldb%fields%get_field_collection("prognostic_fields")
    lbc_fields => modeldb%fields%get_field_collection("lbc_fields")
    moisture_fields => modeldb%fields%get_field_collection("moisture_fields")
    call moisture_fields%get_field("mr", mr_array)
    mr => mr_array%bundle
    call moisture_fields%get_field("moist_dyn", moist_dyn_array)
    moist_dyn => moist_dyn_array%bundle
    derived_fields => modeldb%fields%get_field_collection("derived_fields")
    panel_id => get_panel_id(mesh%get_id())
    con_tracer_last_outer =>  modeldb%fields%get_field_collection("con_tracer_last_outer")

    ! Can't just iterate through the prognostic/diagnostic collections as
    ! some fields are scalars and some fields are vectors, so explicitly
    ! extract all fields from the collections and output each of them in turn
    call prognostic_fields%get_field('theta', theta)
    call prognostic_fields%get_field('u', u)
    call prognostic_fields%get_field('rho', rho)
    call prognostic_fields%get_field('exner', exner)

    ! Get element orders and get the finite element or finite volume height
    element_order_h = theta%get_element_order_h()
    element_order_v = theta%get_element_order_v()

    if (element_order_h > 0 .or. element_order_v > 0) then
      ! Get the finite element height
      height_w3 => get_height_fe(W3, mesh%get_id())
      height_wth => get_height_fe(Wtheta, mesh%get_id())
    else
      ! Get the finite volume height
      height_w3 => get_height_fv(W3, mesh%get_id())
      height_wth => get_height_fv(Wtheta, mesh%get_id())
    end if

    ! Scalar fields
    call write_scalar_diagnostic('rho', rho, &
                                 modeldb%clock, mesh, nodal_output_on_w3)
    call write_scalar_diagnostic('theta', theta, &
                                 modeldb%clock, mesh, nodal_output_on_w3)
    call write_scalar_diagnostic('exner', exner, &
                                 modeldb%clock, mesh, nodal_output_on_w3)
    call write_scalar_diagnostic('height_w3', height_w3, &
                                 modeldb%clock, mesh, nodal_output_on_w3)
    call write_scalar_diagnostic('height_wth', height_wth, &
                                 modeldb%clock, mesh, nodal_output_on_w3)

    if (transport_ageofair) then
      call con_tracer_last_outer%get_field('ageofair',ageofair)
      call write_scalar_diagnostic('ageofair', ageofair, &
                                   modeldb%clock, mesh, nodal_output_on_w3)
    end if

    ! Write out grid_cell area at initialisation only
    if ( use_physics .and. use_xios_io .and. modeldb%clock%is_initialisation() &
         .and. diagnostic_to_be_sampled("init_area_at_msl") ) then
      dA => get_da_msl_proj(twod_mesh%get_id())
      tmp_write_ptr => write_field_generic
      call dA%set_write_behaviour(tmp_write_ptr)
      call dA%write_field("init_area_at_msl")
    endif

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
        call u_in_w2h%write_field("init_u_in_w2h")
        call v_in_w2h%write_field("init_v_in_w2h")
        call w_in_wth%write_field("init_w_in_wth")
      else
        call u_in_w2h%write_field("u_in_w2h")
        call v_in_w2h%write_field("v_in_w2h")
        call w_in_wth%write_field("w_in_wth")
      end if
    else
      call write_vector_diagnostic('u', u, &
                                   modeldb%clock, mesh, nodal_output_on_w3)
    end if
    call write_vorticity_diagnostic( u, exner, modeldb%clock )
#ifdef UM_PHYSICS
    call write_pv_diagnostic( u, theta, rho, exner, modeldb%clock )
#else
    call write_pv_diagnostic( u, theta, rho, modeldb%clock )
#endif

    ! Moisture fields
    if ( moisture_formulation /= moisture_formulation_dry ) then
      do i=1,nummr
        call write_scalar_diagnostic( trim(mr_names(i)), mr(i), &
                                      modeldb%clock, mesh, nodal_output_on_w3 )
      end do
    end if

    if (limited_area) then
      if (output_lbcs) then
        call lbc_fields%get_field('lbc_theta', lbc_theta)
        call lbc_fields%get_field('lbc_u', lbc_u)
        call lbc_fields%get_field('lbc_rho', lbc_rho)
        call lbc_fields%get_field('lbc_exner', lbc_exner)

        call lbc_fields%get_field('lbc_h_u', h_u)
        call lbc_fields%get_field('lbc_v_u', v_u)

        ! Scalar fields
        call write_scalar_diagnostic('lbc_rho', lbc_rho, &
                                     modeldb%clock, mesh, nodal_output_on_w3)
        call write_scalar_diagnostic('lbc_theta', lbc_theta, &
                                     modeldb%clock, mesh, nodal_output_on_w3)
        call write_scalar_diagnostic('lbc_exner', lbc_exner, &
                                     modeldb%clock, mesh, nodal_output_on_w3)
        call write_scalar_diagnostic('readlbc_v_u', v_u, &
                                     modeldb%clock, mesh, nodal_output_on_w3)

        if ( moisture_formulation /= moisture_formulation_dry ) then
          call lbc_fields%get_field('lbc_m_v', lbc_m_v)
          call write_scalar_diagnostic('lbc_m_v', lbc_m_v, &
                                       modeldb%clock, mesh, nodal_output_on_w3)
          call lbc_fields%get_field('lbc_q', lbc_q)
          call write_scalar_diagnostic('lbc_q', lbc_q, &
                                       modeldb%clock, mesh, nodal_output_on_w3)
          call lbc_fields%get_field('lbc_rho_r2', lbc_rho)
          call write_scalar_diagnostic('lbc_rho_r2', lbc_rho, &
                                       modeldb%clock, mesh, nodal_output_on_w3)
        end if

        ! Vector fields
        call write_vector_diagnostic('lbc_u', lbc_u, &
                                     modeldb%clock, mesh, nodal_output_on_w3)
        call write_vector_diagnostic('readlbc_h_u', h_u, &
                                     modeldb%clock, mesh, nodal_output_on_w3)
      endif
    endif

    ! Derived physics fields (only those on W3 or Wtheta)
    if (use_physics .and. use_xios_io .and. .not. modeldb%clock%is_initialisation()) then

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

#ifdef UM_PHYSICS
      ! RH diagnostics
      call rh_diag_alg(exner_in_wth, theta, mr)
      ! Call PMSL algorithm
      call pmsl_alg(exner, derived_fields, theta, twod_mesh)
      ! Pressure level diagnostics
      call pres_lev_diags_alg(derived_fields, theta, exner, mr, moist_dyn)
      ! Wet bulb freezing level
      call freeze_lev_alg(theta, mr, moist_dyn, exner_in_wth)
#endif

      temp_corr_io_value => get_io_value( modeldb%values, 'temperature_correction_io_value')
      call column_total_diagnostics_alg(rho, mr, derived_fields, exner, &
                                        mesh, twod_mesh,             &
                                        temp_corr_io_value%data(1))

    end if

    if (ls_option /= ls_option_file .and. ls_option /= ls_option_analytic) then
      ! Other derived diagnostics with special pre-processing
      ! Don't output for the tangent linear model
      call write_divergence_diagnostic( u, modeldb%clock, mesh )
      call write_hydbal_diagnostic( theta, moist_dyn, exner, mesh )
    end if
    if ( LPROF ) call stop_timing( id, 'gungho_diagnostics_driver' )
  end subroutine gungho_diagnostics_driver

end module gungho_diagnostics_driver_mod
