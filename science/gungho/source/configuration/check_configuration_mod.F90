!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

module check_configuration_mod

  use constants_mod,        only: i_def, l_def
  use mixing_config_mod,    only: viscosity,                                   &
                                  viscosity_mu
  use transport_config_mod, only: operators,                                   &
                                  operators_fv,                                &
                                  operators_fem,                               &
                                  consistent_metric,                           &
                                  fv_horizontal_order,                         &
                                  fv_vertical_order,                           &
                                  cheap_update,                                &
                                  si_outer_transport,                          &
                                  si_outer_transport_none,                     &
                                  profile_size,                                &
                                  scheme,                                      &
                                  splitting,                                   &
                                  horizontal_method,                           &
                                  vertical_method,                             &
                                  reversible,                                  &
                                  log_space,                                   &
                                  max_vert_cfl_calc,                           &
                                  max_vert_cfl_calc_dep_point,                 &
                                  equation_form,                               &
                                  panel_edge_treatment,                        &
                                  panel_edge_treatment_extended_mesh,          &
                                  panel_edge_treatment_special_edges,          &
                                  panel_edge_treatment_remapping,              &
                                  panel_edge_treatment_none,                   &
                                  panel_edge_high_order,                       &
                                  dry_field_name,                              &
                                  field_names,                                 &
                                  use_density_predictor,                       &
                                  enforce_min_value,                           &
                                  horizontal_monotone,                         &
                                  vertical_monotone,                           &
                                  ffsl_splitting,                              &
                                  ffsl_vertical_order,                         &
                                  ffsl_inner_order,                            &
                                  ffsl_outer_order,                            &
                                  dep_pt_stencil_extent,                       &
                                  substep_transport,                           &
                                  substep_transport_off,                       &
                                  adjust_vhv_wind,                             &
                                  ffsl_unity_3d,                   &
                                  wind_mono_top
  use transport_enumerated_types_mod,                                          &
                            only: scheme_mol_3d,                               &
                                  scheme_ffsl_3d,                              &
                                  scheme_split,                                &
                                  split_method_mol,                            &
                                  split_method_ffsl,                           &
                                  split_method_sl,                             &
                                  split_method_null,                           &
                                  equation_form_advective,                     &
                                  equation_form_conservative,                  &
                                  equation_form_consistent,                    &
                                  splitting_strang_hvh,                        &
                                  splitting_strang_vhv,                        &
                                  splitting_none,                              &
                                  monotone_koren,                              &
                                  monotone_relaxed,                            &
                                  monotone_strict,                             &
                                  monotone_clipping,                           &
                                  monotone_qm_pos,                             &
                                  ffsl_splitting_swift, &
                                  ffsl_splitting_cosmic

  implicit none

  private

  public :: check_configuration
  public :: check_any_scheme_mol
  public :: check_any_horizontal_method_mol
  public :: check_any_vertical_method_mol
  public :: check_any_scheme_split
  public :: check_any_scheme_ffsl
  public :: check_any_scheme_split_ffsl
  public :: check_any_scheme_slice
  public :: check_any_hori_scheme_sl
  public :: check_any_hori_scheme_ffsl
  public :: check_any_vert_scheme_sl
  public :: check_any_reversible_sl
  public :: check_any_splitting_hvh
  public :: check_any_splitting_vhv
  public :: check_any_consistent_swift
  public :: check_any_advective_swift
  public :: check_any_consistent_cosmic
  public :: check_any_shifted
  public :: check_wind_shifted
  public :: check_any_eqn_consistent
  public :: check_any_wt_eqn_conservative
  public :: check_moisture_advective
  public :: check_horz_dep_pts
  public :: check_vert_dep_pts
  public :: check_transport_name
  public :: get_required_stencil_depth

contains

  !> @brief Check the namelist configuration for unsupported combinations
  !>        of options and flag up errors and warnings
  !> @param[in] modeldb  The model database object to check
  subroutine check_configuration(modeldb)

    use driver_modeldb_mod,          only: modeldb_type
    use log_mod,                     only: log_event,                          &
                                           log_scratch_space,                  &
                                           LOG_LEVEL_ERROR,                    &
                                           LOG_LEVEL_WARNING,                  &
                                           LOG_LEVEL_INFO
    use constants_mod,               only: EPS,                                &
                                           r_def,                              &
                                           i_def
    use finite_element_config_mod,   only: cellshape,                          &
                                           cellshape_triangle,                 &
                                           element_order_h,                    &
                                           element_order_v,                    &
                                           rehabilitate,                       &
                                           coord_order,                        &
                                           coord_system,                       &
                                           coord_system_native
    use formulation_config_mod,      only: use_physics,                        &
                                           use_wavedynamics,                   &
                                           dlayer_on
    use io_config_mod,               only: write_diag,                         &
                                           use_xios_io
    use planet_config_mod,           only: gravity,                            &
                                           omega,                              &
                                           rd,                                 &
                                           cp,                                 &
                                           p_zero,                             &
                                           scaling_factor
    use timestepping_config_mod,     only: method,                             &
                                           method_semi_implicit,               &
                                           dt,                                 &
                                           alpha,                              &
                                           outer_iterations,                   &
                                           inner_iterations
    use base_mesh_config_mod,        only: geometry,                           &
                                           geometry_spherical,                 &
                                           geometry_planar,                    &
                                           topology,                           &
                                           topology_fully_periodic,            &
                                           topology_non_periodic,              &
                                           prime_mesh_name
    use departure_points_config_mod, only: horizontal_limit,                   &
                                           horizontal_limit_none,              &
                                           horizontal_limit_cap
    use damping_layer_config_mod,    only: dl_base,                            &
                                           dl_str
    use extrusion_config_mod,        only: domain_height, planet_radius
    use mixed_solver_config_mod,     only: reference_reset_time
    use helmholtz_solver_config_mod, only:                                     &
                            helmholtz_solver_preconditioner => preconditioner, &
                            preconditioner_tridiagonal
    implicit none

      class(modeldb_type),  intent(in) :: modeldb

      logical(kind=l_def) :: any_scheme_mol, any_horz_dep_pts
      integer(kind=i_def) :: i, dry_field_splitting

      call log_event( 'Checking gungho configuration...', LOG_LEVEL_INFO )

      ! Check the options in the finite element namelist
      if ( cellshape == cellshape_triangle ) then
        write( log_scratch_space, '(A)' ) 'Triangular elements are unsupported'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( element_order_h < 0 ) then
        write( log_scratch_space, '(A,I4,A)' ) 'Invalid choice: element_order_h ', &
          element_order_h, ' must be non-negative'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( element_order_v < 0 ) then
        write( log_scratch_space, '(A,I4,A)' ) 'Invalid choice: element_order_v ', &
          element_order_v, ' must be non-negative'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( .not. rehabilitate ) then
        write( log_scratch_space, '(A)' ) 'Only rehabilitated W3 function space is allowed'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( coord_order < 0 ) then
        write( log_scratch_space, '(A,I4,A)' ) 'Invalid choice: coordinate order ', &
          coord_order, ' must be non-negative'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( coord_order == 0 .and. &
           ( topology /= topology_non_periodic .or. geometry == geometry_planar ) ) then
        write( log_scratch_space, '(A)' ) 'For planar geometry or periodic meshes, coordinate order must be positive'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      ! Check the options in the formulation namelist
      if ( .not. use_physics .and. .not. use_wavedynamics ) then
        write( log_scratch_space, '(A)' ) 'Wave dynamics and physics turned off'
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if

      ! Check the io namelist
      if ( .not. write_diag ) then
        write( log_scratch_space, '(A)' ) 'Diagnostic output not enabled'
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if

      ! Check the planet namelist
      if ( gravity < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative gravity: ', &
          gravity
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( planet_radius < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative radius: ', &
          planet_radius
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( omega < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative omega: ', &
          omega
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( rd < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative Rd: ', &
          rd
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( cp < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative cp: ', &
          cp
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( p_zero < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative p_zero: ', &
          p_zero
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( scaling_factor < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Invalid choice: Zero or negative scaling factor: ', &
          scaling_factor
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      ! Check the timestepping namelist
      if ( dt < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Zero or negative dt: ', &
          dt
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( method == method_semi_implicit ) then
        if( alpha < 0.5_r_def ) then
          write( log_scratch_space, '(A,E16.8)' ) 'alpha < 1/2 likely to be unstable: ',&
            alpha
          call log_event( log_scratch_space, LOG_LEVEL_WARNING )
        end if
        if ( outer_iterations < 1 ) then
          write( log_scratch_space, '(A,I4)' ) 'Invalid Choice: outer_iterations must be at least 1:', &
          outer_iterations
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( inner_iterations < 1 ) then
          write( log_scratch_space, '(A,I4)' ) 'Invalid Choice: inner_iterations must be at least 1:', &
          inner_iterations
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if

      ! Check the transport namelist
      if ( geometry == geometry_spherical .and.  consistent_metric) then
        write( log_scratch_space, '(A)' ) 'Consistent metric option only valid for planar geometries'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      any_scheme_mol = check_any_scheme_mol()
      if (any_scheme_mol) then
        ! Check that flux orders are even
        if ( mod(fv_horizontal_order,2_i_def) /= 0_i_def ) then
          write( log_scratch_space, '(A)' ) 'fv_horizontal_order must be even'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( mod(fv_vertical_order,2_i_def) /= 0_i_def ) then
          write( log_scratch_space, '(A)' ) 'fv_vertical_order must be even'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
      ! Check FFSL orders are in allowed range
      if (ffsl_inner_order < 0 .or. ffsl_inner_order > 2) then
        write( log_scratch_space, '(A)' ) 'ffsl_inner_order must be between 0 and 2'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if (ffsl_outer_order < 0 .or. ffsl_outer_order > 2) then
        write( log_scratch_space, '(A)' ) 'ffsl_outer_order must be between 0 and 2'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( panel_edge_treatment == panel_edge_treatment_extended_mesh ) then
        if ( coord_system /=  coord_system_native ) then
          write( log_scratch_space, '(A)' )                                    &
              'extended_mesh panel_edge_treatment only valid for native coordinates'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( coord_order /= 1 ) then
          write( log_scratch_space, '(A)' )                                    &
              'extended_mesh panel_edge_treatment only valid for linear coord_order'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if
      if ( panel_edge_treatment /= panel_edge_treatment_none ) then
        if ( geometry /= geometry_spherical ) then
          write( log_scratch_space, '(A)' ) 'panel_edge_treatment only valid for spherical geometry'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( topology /=  topology_fully_periodic) then
          write( log_scratch_space, '(A)' ) 'panel_edge_treatment only valid for fully periodic topology'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( any_scheme_mol ) then
          if ( coord_system /=  coord_system_native ) then
            write( log_scratch_space, '(A)' )                                    &
                'MoL panel_edge_treatment only valid for native coordinates'
            call log_event( log_scratch_space, LOG_LEVEL_ERROR )
          end if
          if ( coord_order /= 1 ) then
            write( log_scratch_space, '(A)' )                                    &
                'MoL panel_edge_treatment only valid for linear coord_order'
            call log_event( log_scratch_space, LOG_LEVEL_ERROR )
          end if
        end if
      end if
      if ( si_outer_transport /= si_outer_transport_none .AND. cheap_update) then
        write( log_scratch_space, '(A)' ) 'Cheap update cannot be used with si_outer_transport options'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if
      if ( substep_transport /= substep_transport_off .AND. cheap_update) then
        write( log_scratch_space, '(A)' ) 'Cheap update cannot be used with substep_transport'
        call log_event( log_scratch_space, LOG_LEVEL_ERROR )
      end if

      ! Find splitting used by dry field
      dry_field_splitting = -1
      do i = 1, profile_size
        if ( trim(field_names(i)) == trim(dry_field_name) ) then
          dry_field_splitting = splitting(i)
          if ( equation_form(i) /= equation_form_conservative ) then
            write( log_scratch_space, '(A)') trim(field_names(i)) // ' variable ' // &
              'is the specified dry field, but is not using the conservative ' // &
              'form of the transport equation'
            call log_event(log_scratch_space, LOG_LEVEL_ERROR)
          end if
          exit
        end if
      end do

      ! Check some combinations of options, variable-by-variable
      do i = 1, profile_size
        if ( splitting(i) /= splitting_none .AND. scheme(i) /= scheme_split ) then
          write( log_scratch_space, '(A)') trim(field_names(i)) // ' variable ' // &
            'is not being transported with a split transport scheme, so it ' //    &
            'must have its splitting option set to none'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        if ( substep_transport /= substep_transport_off .and. &
             (scheme(i) == scheme_mol_3d .or.                 &
              vertical_method(i) == split_method_mol .or.     &
              horizontal_method(i) == split_method_mol) ) then
          call log_event('Substepping the whole transport is only valid with FFSL ' // &
                         'for all variables. Set substep_transport=off to ' //         &
                         'run with MoL', LOG_LEVEL_ERROR)
        end if

        if ( scheme(i) == scheme_mol_3d .and. vertical_method(i) /= split_method_mol ) then
          write( log_scratch_space, '(A)') trim(field_names(i)) // ' variable ' // &
            'is set to be transported with 3D MoL, so its vertical method must also be MoL'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        if ( scheme(i) == scheme_mol_3d .and. horizontal_method(i) /= split_method_mol ) then
          write( log_scratch_space, '(A)') trim(field_names(i)) // ' variable ' // &
            'is set to be transported with 3D MoL, so its horizontal method must also be MoL'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        if ( scheme(i) == scheme_ffsl_3d .and. vertical_method(i) /= split_method_ffsl ) then
          write( log_scratch_space, '(A)') trim(field_names(i)) // ' variable ' // &
            'is set to be transported with 3D FFSL, so its vertical method must also be FFSL'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        if ( scheme(i) == scheme_ffsl_3d .and. horizontal_method(i) /= split_method_ffsl ) then
          write( log_scratch_space, '(A)') trim(field_names(i)) // ' variable ' // &
            'is set to be transported with 3D FFSL, so its horizontal method must also be FFSL'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        if ( vertical_method(i) == split_method_ffsl .AND. ffsl_vertical_order(i) == 2    &
            .AND. .NOT. reversible(i) ) then
          write( log_scratch_space, '(A)') trim(field_names(i)) // ' variable ' // &
            'is being transported with a reversible form of the FFSL scheme, ' // &
            'so it must also have the "reversible" option set to .true.'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        if ( horizontal_method(i) == split_method_ffsl                         &
             .AND. ffsl_splitting(i) == ffsl_splitting_cosmic                  &
             .AND. panel_edge_treatment == panel_edge_treatment_remapping ) then
          write( log_scratch_space, '(A)') 'Remapping treatment of panel ' //  &
            'edges is not implemented for the COSMIC FFSL splitting'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        if ( operators == operators_fem .and. enforce_min_value(i) ) then
          write( log_scratch_space, '(A)') trim(field_names(i)) // ' variable ' // &
            'is set to use the enforce_min_value option, but this cannot ' // &
            'be used with the FEM operators'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        if ( (vertical_monotone(i) == monotone_strict .or.    &
              vertical_monotone(i) == monotone_relaxed) .and. &
              .not. (vertical_method(i) == split_method_ffsl .or.      &
                     vertical_method(i) == split_method_sl) ) then
          write( log_scratch_space, '(A)') trim(field_names(i)) // ' variable ' // &
            'is set to use strict/relaxed vertical monotonicity, but this is ' // &
            'incompatible with the choice of vertical method'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        if ( (horizontal_monotone(i) == monotone_strict .or.    &
              horizontal_monotone(i) == monotone_relaxed) .and. &
              .not. (horizontal_method(i) == split_method_ffsl .or.        &
                     horizontal_method(i) == split_method_sl) ) then
          write( log_scratch_space, '(A)') trim(field_names(i)) // ' variable ' // &
            'is set to use strict/relaxed horizontal monotonicity, but this is ' // &
            'incompatible with the choice of horizontal method'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        if ( (vertical_monotone(i) == monotone_qm_pos) .and. &
              .not. (vertical_method(i) == split_method_ffsl)   .and. &
              .not. (ffsl_vertical_order(i) == 2) ) then
          write( log_scratch_space, '(A)') trim(field_names(i)) // ' variable ' // &
            'is set to use quasi-monotone positive vertical monotonicity, but this is ' // &
            'incompatible with the choice of vertical method'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        if ( (horizontal_monotone(i) == monotone_qm_pos) .and. &
              .not. (horizontal_method(i) == split_method_ffsl) ) then
          write( log_scratch_space, '(A)') trim(field_names(i)) // ' variable ' // &
            'is set to use quasi-monotone positive horizontal monotonicity, but this is ' // &
            'incompatible with the choice of horizontal method'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        if ( (vertical_monotone(i) == monotone_koren .or.      &
              vertical_monotone(i) == monotone_clipping) .and. &
              .not. (vertical_method(i) == split_method_mol) ) then
          write( log_scratch_space, '(A)') trim(field_names(i)) // ' variable ' // &
            'is set to use Koren/clipping vertical monotonicity, but this is ' // &
            'incompatible with the choice of vertical method'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        if ( (horizontal_monotone(i) == monotone_koren .or.      &
              horizontal_monotone(i) == monotone_clipping) .and. &
              .not. (horizontal_method(i) == split_method_mol) ) then
          write( log_scratch_space, '(A)') trim(field_names(i)) // ' variable ' // &
            'is set to use Koren/clipping horizontal monotonicity, but this is ' // &
            'incompatible with the choice of horizontal method'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        if ( equation_form(i) == equation_form_consistent .and. &
             splitting(i) /= dry_field_splitting ) then
          write( log_scratch_space, '(A)') trim(field_names(i)) // ' variable ' // &
            'is set to use consistent transport, but it is using a different ' // &
            'transport splitting to the specified dry field.'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        if ( equation_form(i) == equation_form_consistent .and. &
             scheme(i) == scheme_ffsl_3d ) then
          write( log_scratch_space, '(A)') trim(field_names(i)) //  &
            'variable is set to use consistent transport, but this is ' // &
            'not yet implemented with 3D FFSL.'
          call log_event(log_scratch_space, LOG_LEVEL_ERROR)
        end if

        ! For 3D unity transport, require all transported fields to use the
        ! same FFSL-FFSL splitting
        if ( ffsl_unity_3d ) then
          if ( splitting(i) /= dry_field_splitting ) then
            call log_event(                                                    &
              '3D unity transport can only be used when all variables '        &
              // 'are transported with the same splitting', LOG_LEVEL_ERROR)
          else if ( vertical_method(i) /= split_method_ffsl                    &
                    .or. horizontal_method(i) /= split_method_ffsl ) then
            call log_event(                                                    &
              '3D unity transport can only be used when all variables '        &
              // 'are using FFSL for vertical and horizontal transport', LOG_LEVEL_ERROR)
          end if
        end if

        ! For applying monotonicity to wind at the model top, require vertical
        ! wind transport to use FFSL
        if ( wind_mono_top ) then
          if ( trim(field_names(i)) == 'wind' .and. .not.                      &
              ( vertical_method(i) == split_method_ffsl .and.                  &
                .not. reversible(i) .and. ffsl_vertical_order(i) == 1 ) ) then
            call log_event(                                                    &
                  'The wind_mono_top option only be used when the wind ' //    &
                  'components are transported vertically with FFSL',           &
                  LOG_LEVEL_ERROR                                              &
            )
          end if
        end if

        ! Can only adjust the V-H-V wind when using FFSL with 3D unity transport
        ! and Strang V-H-V splitting
        if ( adjust_vhv_wind ) then
          if ( .not. ffsl_unity_3d ) then
            call log_event(                                                    &
              'adjust_vhv_wind only implemented with ffsl_unity_3d set to true', LOG_LEVEL_ERROR)
          end if
          if ( splitting(i) /= splitting_strang_vhv ) then
            call log_event(                                                    &
              'adjust_vhv_wind only implemented for Strang VHV splitting', LOG_LEVEL_ERROR)
          end if
        end if

      end do

      ! Check the departure points namelist
      any_horz_dep_pts = check_horz_dep_pts()
      if ( any_horz_dep_pts ) then
        ! Warn that having no limit on horizontal departure points may lead to model
        ! failure, and that capping the departure distance will cap the advecting wind
        if ( horizontal_limit == horizontal_limit_none ) then
          write( log_scratch_space, '(A)' ) &
          'If the maximum horizontal departure distance exceeds dep_pt_stencil_extent this will result in the model failure'
          call log_event( log_scratch_space, LOG_LEVEL_WARNING )
        else if ( horizontal_limit == horizontal_limit_cap ) then
          write( log_scratch_space, '(A)' ) &
          'The maximum horizontal departure distance and thus the advecting wind will be capped'
          call log_event( log_scratch_space, LOG_LEVEL_WARNING )
        end if
      end if

      ! Check the mixing namelist
      if ( viscosity .and. geometry == geometry_spherical ) then
        write( log_scratch_space, '(A)' ) 'Viscosity might not work in spherical domains'
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( viscosity .and. viscosity_mu < EPS ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Negative viscosity coefficient: Anti-diffusion: ', &
          viscosity_mu
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if

      ! Check the damping layer namelist
      if ( dlayer_on .and. (dl_base < 0.0_r_def .or. dl_base > domain_height) ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Damping layer base lies outside of domain: ',&
          dl_base
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if
      if ( dlayer_on .and. dl_str < 0.0_r_def ) then
        write( log_scratch_space, '(A,E16.8)' ) 'Damping layer strength is negative: ',&
          dl_str
        call log_event( log_scratch_space, LOG_LEVEL_WARNING )
      end if

      ! Check for options that are invalid with higher order elements
      if ( element_order_h > 0 .or. element_order_v > 0 ) then
        if ( operators == operators_fv ) then
          write( log_scratch_space, '(A)' ) &
          'FV transport operators only valid for element_order_h = 0 and element_order_v = 0'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( viscosity ) then
          write( log_scratch_space, '(A)' ) &
          'Viscosity only valid for element_order_h = 0 and element_order_v = 0'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( helmholtz_solver_preconditioner == preconditioner_tridiagonal ) then
          write( log_scratch_space, '(A)' ) &
          'Tridiagonal helmholtz preconditioner only valid for  element_order_h = 0 and element_order_v = 0'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
        if ( use_xios_io ) then
          write( log_scratch_space, '(A)' ) 'xios output may not work with element order > 0'
          call log_event( log_scratch_space, LOG_LEVEL_WARNING )
        end if
      end if

      if ( method == method_semi_implicit ) then
        ! Check the mixed solver namelist
        if (reference_reset_time < dt ) then
          write( log_scratch_space, '(A)' ) 'reference_reset_time must be greater than or equal to time step size dt'
          call log_event( log_scratch_space, LOG_LEVEL_ERROR )
        end if
      end if

      call log_event( '...Check gungho config done', LOG_LEVEL_INFO )

  end subroutine check_configuration


  !> @brief   Determine required stencil depth for the current configuration.
  !> @details Depending on the choice of science schemes the required local
  !>          mesh needs to support the anticipated stencils. This function
  !>          returns required stencil depth that needs to be supported.
  !> @return  stencil_depth
  !>
  !===========================================================================
  function get_required_stencil_depth() result(stencil_depth)

    implicit none

    integer(kind=i_def) :: stencil_depth
    integer(kind=i_def) :: sl_depth, special_edge_pts
    logical(kind=l_def) :: any_horz_dep_pts

    stencil_depth = 2

    if (operators == operators_fv) then
      ! Need larger halos for fv operators
      stencil_depth  = max( stencil_depth, fv_horizontal_order/2 )
    end if

    any_horz_dep_pts = check_horz_dep_pts()

    if ( any_horz_dep_pts ) then
      ! When an SL scheme is used, the halo depth should be large enough to
      ! encompass the largest anticipated Courant number (effectively the
      ! departure distance in the SL scheme), plus any extra cells required for
      ! the reconstruction of the field at the departure point.
      ! The maximum anticipated Courant number is set by the user through the
      ! dep_pt_stencil_extent namelist variable. Other considerations are:
      ! - the order of reconstruction
      ! - whether special edge treatment is used (this shifts the stencil by 1)

      if ( panel_edge_treatment == panel_edge_treatment_special_edges          &
           .AND. panel_edge_high_order ) then
        special_edge_pts = 1
      else
        special_edge_pts = 0
      end if

      sl_depth = (                                  &
        dep_pt_stencil_extent                       & ! max anticipated Courant
        + max( ffsl_inner_order, ffsl_outer_order ) & ! max reconstruction order
        + special_edge_pts                          & ! special edge treatment
      )

      if ( panel_edge_treatment == panel_edge_treatment_remapping ) then
        if ( panel_edge_high_order ) then
          sl_depth = max( sl_depth, 3 )
        else
          sl_depth = max( sl_depth, 2 )
        end if
      end if

      stencil_depth = max( stencil_depth, sl_depth )
    end if

  end function get_required_stencil_depth

  !> @brief   Determine whether any of the transport schemes are MoL
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any are using the Method of Lines
  !>          scheme
  !> @return  any_scheme_mol
  function check_any_scheme_mol() result(any_scheme_mol)

    implicit none

    logical(kind=l_def) :: any_scheme_mol

    any_scheme_mol = (check_any_horizontal_method_mol() .or. &
                      check_any_vertical_method_mol())

  end function check_any_scheme_mol

  !> @brief   Determine whether any transport scheme uses horizontal MoL
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any use horizontal MoL
  !> @return  any_horizontal_mol
  function check_any_horizontal_method_mol() result(any_horizontal_mol)

    implicit none

    logical(kind=l_def) :: any_horizontal_mol
    integer(kind=i_def) :: i

    any_horizontal_mol = .false.

    do i = 1, profile_size
      if ( ( scheme(i) == scheme_mol_3d ) .or.                      &
           ( scheme(i) == scheme_split .and.                        &
             horizontal_method(i) == split_method_mol ) ) then
        any_horizontal_mol = .true.
        exit
      end if
    end do

  end function check_any_horizontal_method_mol

  !> @brief   Determine whether any transport scheme uses vertical MoL
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any use vertical MoL
  !> @return  any_vertical_mol
  function check_any_vertical_method_mol() result(any_vertical_mol)

    implicit none

    logical(kind=l_def) :: any_vertical_mol
    integer(kind=i_def) :: i

    any_vertical_mol = .false.

    do i = 1, profile_size
      if ( ( scheme(i) == scheme_mol_3d ) .or.                      &
           ( scheme(i) == scheme_split .and.                        &
             vertical_method(i) == split_method_mol ) ) then
        any_vertical_mol = .true.
        exit
      end if
    end do

  end function check_any_vertical_method_mol

  !> @brief   Determine whether any of the transport schemes are split
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any are using the split vertical-
  !>          horizontal scheme
  !> @return  any_scheme_split
  function check_any_scheme_split() result(any_scheme_split)

    implicit none

    logical(kind=l_def) :: any_scheme_split
    integer(kind=i_def) :: i

    any_scheme_split = .false.

    do i = 1, profile_size
      if ( scheme(i) == scheme_split ) then
        any_scheme_split = .true.
        exit
      end if
    end do

  end function check_any_scheme_split

  !> @brief   Determine whether any of the transport schemes are FFSL
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any are using the Flux-Form
  !>          Semi-Lagrangian scheme
  !> @return  any_scheme_ffsl
  function check_any_scheme_ffsl() result(any_scheme_ffsl)

    implicit none

    logical(kind=l_def) :: any_scheme_ffsl
    integer(kind=i_def) :: i

    any_scheme_ffsl = .false.

    do i = 1, profile_size
      if ( ( scheme(i) == scheme_ffsl_3d ) .or.                      &
           ( scheme(i) == scheme_split .and.                        &
             ( vertical_method(i) == split_method_ffsl .or.          &
               horizontal_method(i) == split_method_ffsl ) ) ) then
        any_scheme_ffsl = .true.
        exit
      end if
    end do

  end function check_any_scheme_ffsl

  !> @brief   Determine whether any of the transport schemes are split 3D FFSL
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any are using the split Flux-Form
  !>          Semi-Lagrangian scheme
  !> @return  any_scheme_ffsl
  function check_any_scheme_split_ffsl() result(any_scheme_ffsl)

    implicit none

    logical(kind=l_def) :: any_scheme_ffsl
    integer(kind=i_def) :: i

    any_scheme_ffsl = .false.

    do i = 1, profile_size
      if ( scheme(i) == scheme_split .and.                        &
           (vertical_method(i) == split_method_ffsl .or.          &
            vertical_method(i) == split_method_null) .and.        &
           (horizontal_method(i) == split_method_ffsl .or.        &
            horizontal_method(i) == split_method_null) ) then
        any_scheme_ffsl = .true.
        exit
      end if
    end do

  end function check_any_scheme_split_ffsl

  !> @brief   Determine whether any of the vertical transport schemes are SLICE
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any are using SLICE
  !> @return  any_scheme_slice
  function check_any_scheme_slice() result(any_scheme_slice)

    implicit none

    logical(kind=l_def) :: any_scheme_slice
    integer(kind=i_def) :: i

    any_scheme_slice = .false.

    do i = 1, profile_size
      if ( scheme(i) == scheme_split .and.                       &
           vertical_method(i) == split_method_sl .and.           &
           (equation_form(i)  == equation_form_conservative .or. &
            equation_form(i)  == equation_form_consistent) ) then
        any_scheme_slice = .true.
        exit
      end if
    end do

  end function check_any_scheme_slice

  !> @brief   Determine whether any of the horizontal transport schemes are SL
  !> @details Loops through the horizontal transport schemes specified for different
  !>          variables and determines whether any are using the semi-Lagrangian scheme.
  !> @return  any_hori_scheme_sl
  function check_any_hori_scheme_sl() result(any_hori_scheme_sl)

    implicit none

    logical(kind=l_def) :: any_hori_scheme_sl
    integer(kind=i_def) :: i

    any_hori_scheme_sl = .false.

    do i = 1, profile_size
      if ( scheme(i) == scheme_split .and. &
           horizontal_method(i) == split_method_sl ) then
        any_hori_scheme_sl = .true.
        exit
      end if
    end do

  end function check_any_hori_scheme_sl

  !> @brief   Determine whether any of the horizontal transport schemes are FFSL
  !> @details Loops through the horizontal transport schemes specified for
  !!          different variables and determines whether any are using the
  !!          Flux-Form Semi-Lagrangian scheme.
  !> @return  Logical for whether a horizontal scheme is FFSL
  function check_any_hori_scheme_ffsl() result(any_hori_scheme_ffsl)

    implicit none

    logical(kind=l_def) :: any_hori_scheme_ffsl
    integer(kind=i_def) :: i

    any_hori_scheme_ffsl = .false.

    do i = 1, profile_size
      if ( scheme(i) == scheme_split .and. &
           horizontal_method(i) == split_method_ffsl ) then
        any_hori_scheme_ffsl = .true.
        exit
      end if
    end do

  end function check_any_hori_scheme_ffsl

  !> @brief   Determine whether any of the vertical transport schemes are SL
  !> @details Loops through the vertical transport schemes specified for different
  !>          variables and determines whether any are using the semi-Lagrangian scheme.
  !> @return  any_vert_scheme_sl
  function check_any_vert_scheme_sl() result(any_vert_scheme_sl)

    implicit none

    logical(kind=l_def) :: any_vert_scheme_sl
    integer(kind=i_def) :: i

    any_vert_scheme_sl = .false.

    do i = 1, profile_size
      if ( scheme(i) == scheme_split .and. &
           vertical_method(i) == split_method_sl ) then
        any_vert_scheme_sl = .true.
        exit
      end if
    end do

  end function check_any_vert_scheme_sl

  !> @brief   Determine whether any of the vertical transport schemes are reversible SL
  !> @details Loops through the vertical transport schemes specified for different
  !>          variables and determines whether any are using the reversible semi-Lagrangian scheme.
  !> @return  any_reversible_sl
  function check_any_reversible_sl() result(any_reversible_sl)

    implicit none

    logical(kind=l_def) :: any_reversible_sl
    integer(kind=i_def) :: i

    any_reversible_sl = .false.

    do i = 1, profile_size
      if ( scheme(i) == scheme_split .and.             &
           vertical_method(i) == split_method_sl .and. &
           reversible(i) ) then
        any_reversible_sl = .true.
        exit
      end if
    end do

  end function check_any_reversible_sl

  !> @brief   Determine whether any of the split transport schemes use
  !!          Strang HVH splitting
  !> @details Loops through the transport splitting specified for different
  !!          variables and determines whether any are using the Strang
  !!          horizontal-vertical-horizontal splitting
  !> @return  any_splitting_hvh
  function check_any_splitting_hvh() result(any_splitting_hvh)

    implicit none

    logical(kind=l_def) :: any_splitting_hvh
    integer(kind=i_def) :: i

    any_splitting_hvh = .false.

    do i = 1, profile_size
      if ( scheme(i) == scheme_split .AND. &
           splitting(i) == splitting_strang_hvh ) then
        any_splitting_hvh = .true.
        exit
      end if
    end do

  end function check_any_splitting_hvh

  !> @brief   Determine whether any of the split transport schemes use
  !!          Strang VHV splitting
  !> @details Loops through the transport splitting specified for different
  !!          variables and determines whether any are using the Strang
  !!          vertical-horizontal-vertical splitting
  !> @return  any_splitting_vhv
  function check_any_splitting_vhv() result(any_splitting_vhv)

    implicit none

    logical(kind=l_def) :: any_splitting_vhv
    integer(kind=i_def) :: i

    any_splitting_vhv = .false.

    do i = 1, profile_size
      if ( scheme(i) == scheme_split .AND. &
           splitting(i) == splitting_strang_vhv ) then
        any_splitting_vhv = .true.
        exit
      end if
    end do

  end function check_any_splitting_vhv

  !> @brief   Determine if SWIFT splitting is used for conservative tracers
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any are using SWIFT
  !> @return  Logical for whether SWIFT splitting is being used
  function check_any_consistent_swift() result(any_swift)

    implicit none

    logical(kind=l_def) :: any_swift
    integer(kind=i_def) :: i

    any_swift = .false.

    do i = 1, profile_size
      if ( equation_form(i) == equation_form_consistent .and. &
           horizontal_method(i) == split_method_ffsl    .and. &
           scheme(i) == scheme_split                    .and. &
           ffsl_splitting(i) == ffsl_splitting_swift ) then
        any_swift = .true.
        exit
      end if
    end do

  end function check_any_consistent_swift

  !> @brief   Determine if SWIFT splitting is being used for advective fields
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any are using SWIFT
  !> @return  Logical for whether SWIFT splitting is being used
  function check_any_advective_swift() result(any_swift)

    implicit none

    logical(kind=l_def) :: any_swift
    integer(kind=i_def) :: i

    any_swift = .false.

    do i = 1, profile_size
      if ( (equation_form(i) == equation_form_advective     .or.  &
            equation_form(i) == equation_form_conservative) .and. &
            horizontal_method(i) == split_method_ffsl       .and. &
            scheme(i) == scheme_split                       .and. &
            ffsl_splitting(i) == ffsl_splitting_swift ) then
        any_swift = .true.
        exit
      end if
    end do

  end function check_any_advective_swift

  !> @brief   Determine whether consistent COSMIC splitting is being used
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any are using COSMIC splitting
  !> @return  Logical for whether consistent COSMIC splitting is being used
  function check_any_consistent_cosmic() result(any_cosmic)

    implicit none

    logical(kind=l_def) :: any_cosmic
    integer(kind=i_def) :: i

    any_cosmic = .false.

    do i = 1, profile_size
      if ( equation_form(i) == equation_form_consistent .and. &
           horizontal_method(i) == split_method_ffsl    .and. &
           scheme(i) == scheme_split                    .and. &
           ffsl_splitting(i) == ffsl_splitting_cosmic ) then
        any_cosmic = .true.
        exit
      end if
    end do

  end function check_any_consistent_cosmic

  !> @brief   Determine whether any of the transport schemes need shifted mesh
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any need the shifted mesh
  !> @return  any_shifted
  function check_any_shifted() result(any_shifted)

    use formulation_config_mod, only: moisture_formulation, &
                                      moisture_formulation_dry
    use io_config_mod,          only: write_conservation_diag

    implicit none

    logical(kind=l_def)      :: any_shifted
    integer(kind=i_def)      :: i

    any_shifted = .false.

    ! Need a shifted mesh if:
    ! (a) a variable uses the conservative or consistent transport equation
    ! (b) a Wtheta variable uses FFSL
    do i = 1, profile_size
      ! Check for a variable using conservative/consistent equation
      ! (but don't include "dry_field" which will never use shifted grid)
      if ( (equation_form(i) == equation_form_conservative .and. &
            field_names(i) /= dry_field_name) .or.               &
            equation_form(i) == equation_form_consistent ) then
          any_shifted = .true.
          return
      end if

      ! Check if there is a transport scheme using FFSL
      select case (scheme(i))
        ! It could be either 3D FFSL or split scheme using FFSL
        case (scheme_ffsl_3d)
          if (equation_form(i) == equation_form_advective) then
            any_shifted = .true.
            return
          end if

        case (scheme_split)
          if ( vertical_method(i) == split_method_ffsl &
                .or. horizontal_method(i) == split_method_ffsl ) then
            if (equation_form(i) == equation_form_advective) then
              any_shifted = .true.
              return
            end if
          end if

      end select
    end do

    if (write_conservation_diag .and. &
        moisture_formulation /= moisture_formulation_dry) any_shifted = .true.

  end function check_any_shifted

  !> @brief   Determine whether the wind is needed on the shifted mesh
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether the transporting wind is needed
  !!          on the shifted mesh
  !> @return  wind_shifted
  function check_wind_shifted() result(wind_shifted)

    implicit none

    logical(kind=l_def)      :: wind_shifted
    integer(kind=i_def)      :: i

    wind_shifted = .false.

    ! Need to compute a transporting wind on the shifted mesh if:
    ! (a) there is a conservative variable that isn't the "dry_field"
    ! (b) an advected Wtheta variable uses FFSL
    do i = 1, profile_size
      ! Check for a variable using conservative/consistent equation
      ! (but don't include "dry_field" which will never use shifted grid)
      if ( equation_form(i) == equation_form_conservative &
            .and. field_names(i) /= dry_field_name ) then
        wind_shifted = .true.
        return
      else if ( equation_form(i) == equation_form_advective           &
                .and. (vertical_method(i) == split_method_ffsl .or.   &
                       horizontal_method(i) == split_method_ffsl .or. &
                       scheme(i) == scheme_ffsl_3d) ) then
        wind_shifted = .true.
        return
      end if
    end do

  end function check_wind_shifted

  !> @brief   Determine whether any of the transport equations are consistent
  !> @details Loops through the transport equations specified for different
  !>          variables and determines whether any are using consistent form
  !> @return  any_eqn_consistent
  function check_any_eqn_consistent() result(any_eqn_consistent)

    implicit none

    logical(kind=l_def) :: any_eqn_consistent
    integer(kind=i_def) :: i

    any_eqn_consistent = .false.

    do i = 1, profile_size
      if ( equation_form(i) == equation_form_consistent ) then
        any_eqn_consistent = .true.
        exit
      end if
    end do

  end function check_any_eqn_consistent

  !> @brief   Determine whether moisture transport is advective
  !> @details Loops through the transport equations specified for different
  !!          variables and determines whether moisture transport is advective,
  !!          as opposed to the two conservative options.
  !> @return  moisture_advective
  function check_moisture_advective() result(moisture_advective)

    implicit none

    logical(kind=l_def) :: moisture_advective
    integer(kind=i_def) :: i

    moisture_advective = .false.

    do i = 1, profile_size
      if ( field_names(i) == 'mr' .and. &
           equation_form(i) == equation_form_advective ) then
        moisture_advective = .true.
        exit
      end if
    end do

  end function check_moisture_advective

  !> @brief   Determine whether any of the Wtheta transport eqns are conservative
  !> @details Loops through the transport equations specified for different
  !!          variables and determines whether any are using conservative form
  !> @return  any_wt_eqn_conservative
  function check_any_wt_eqn_conservative() result(any_wt_eqn_conservative)

    implicit none

    logical(kind=l_def) :: any_wt_eqn_conservative
    integer(kind=i_def) :: i

    any_wt_eqn_conservative = .false.

    do i = 1, profile_size
      if ( equation_form(i) == equation_form_conservative &
           .and. field_names(i) /= dry_field_name ) then
        any_wt_eqn_conservative = .true.
        exit
      end if
    end do

  end function check_any_wt_eqn_conservative

  !> @brief   Determine whether horizontal departure points need computing
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any are using a scheme that
  !>          requires horizontal departure points to be computed
  !> @return  any_horz_dep_pts
  function check_horz_dep_pts() result(any_horz_dep_pts)

    implicit none

    logical(kind=l_def) :: any_horz_dep_pts
    integer(kind=i_def) :: i

    any_horz_dep_pts = .false.

    do i = 1, profile_size
      if ( ( scheme(i) == scheme_ffsl_3d ) .or.                     &
           ( scheme(i) == scheme_split .and.                        &
             horizontal_method(i) == split_method_ffsl ) .or.       &
           ( scheme(i) == scheme_split .and.                        &
             horizontal_method(i) == split_method_sl ) ) then
        any_horz_dep_pts = .true.
        exit
      end if
    end do

  end function check_horz_dep_pts

  !> @brief   Determine whether vertical departure points need computing
  !> @details Loops through the transport schemes specified for different
  !>          variables and determines whether any are using a scheme that
  !>          requires vertical departure points to be computed
  !> @return  any_vert_dep_pts
  function check_vert_dep_pts() result(any_vert_dep_pts)

    implicit none

    logical(kind=l_def) :: any_vert_dep_pts
    integer(kind=i_def) :: i

    any_vert_dep_pts = .false.

    do i = 1, profile_size
      if ( ( scheme(i) == scheme_ffsl_3d ) .or.                     &
           ( scheme(i) == scheme_split .and.                        &
             vertical_method(i) /= split_method_mol ) .or.          &
           ( max_vert_cfl_calc == max_vert_cfl_calc_dep_point ) ) then
        any_vert_dep_pts = .true.
        exit
      end if
    end do

  end function check_vert_dep_pts


  !> @brief     Determine whether any transport of field name is requested
  !> @details   Loops through the transport names and returns true if
  !>            input name is present.
  !>
  !> @param[in] name       The transport_metadata name
  !> @return    name_exists
  function check_transport_name(name) result(name_exists)

    implicit none

    character(*), intent(in) :: name
    logical(kind=l_def) :: name_exists
    integer(kind=i_def) :: i

    name_exists = .false.

    do i = 1, profile_size
      if ( trim(field_names(i)) == trim(name) ) then
        name_exists = .true.
        exit
      end if
    end do

  end function check_transport_name

end module check_configuration_mod
