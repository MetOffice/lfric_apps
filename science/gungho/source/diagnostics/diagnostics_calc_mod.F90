!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!>  @brief Module for computing and outputting derived diagnostics
!!
!!  @details Computes various derived diagnostics that are written out
!!           by the diagnostic system. Also computes norms of some fields
!!           which are written to the output log.
!-------------------------------------------------------------------------------
module diagnostics_calc_mod

  use constants_mod,                 only: i_def, r_def, str_max_filename, l_def
  use diagnostic_alg_mod,            only: divergence_diagnostic_alg,   &
                                           hydbal_diagnostic_alg,       &
                                           vorticity_diagnostic_alg,    &
                                           potential_vorticity_diagnostic_alg
  use io_config_mod,                 only: use_xios_io,          &
                                           nodal_output_on_w3
  use files_config_mod,              only: diag_stem_name
  use sci_project_output_mod,        only: project_output
  use io_mod,                        only: ts_fname, &
                                           nodal_write_field
  use lfric_xios_write_mod,          only: write_field_generic
  use diagnostics_io_mod,            only: write_scalar_diagnostic,     &
                                           write_vector_diagnostic
  use initialise_diagnostics_mod,    only: diagnostic_to_be_sampled, &
                                           init_diag => init_diagnostic_field
  use field_mod,                     only: field_type
  use field_parent_mod,              only: write_interface
  use fs_continuity_mod,             only: W3
  use model_clock_mod,               only: model_clock_type
  use moist_dyn_mod,                 only: num_moist_factors
  use log_mod,                       only: log_event,         &
                                           log_set_level,     &
                                           log_scratch_space, &
                                           LOG_LEVEL_ERROR,   &
                                           LOG_LEVEL_INFO,    &
                                           LOG_LEVEL_DEBUG,   &
                                           LOG_LEVEL_TRACE
  use mesh_mod,                      only: mesh_type
  use sci_geometric_constants_mod,   only: get_coordinates,      &
                                           get_panel_id

  implicit none
  private
  public :: write_divergence_diagnostic, &
            write_pv_diagnostic,         &
            write_hydbal_diagnostic,     &
            write_vorticity_diagnostic

contains

!-------------------------------------------------------------------------------
!>  @brief    Handles divergence diagnostic processing
!!
!!  @details  Handles divergence diagnostic processing
!!
!!> @param[in] u_field     The u field
!!> @param[in] ts          Timestep
!!> @param[in] mesh        Mesh
!-------------------------------------------------------------------------------

subroutine write_divergence_diagnostic(u_field, clock, mesh)
  implicit none

  type(field_type),        intent(in)    :: u_field
  class(model_clock_type), intent(in)    :: clock
  type(mesh_type),         intent(in), pointer :: mesh

  type(field_type)                :: div_field
  real(r_def)                     :: l2_norm


  procedure(write_interface), pointer  :: tmp_write_ptr

  if (diagnostic_to_be_sampled('divergence') .or. &
      diagnostic_to_be_sampled('init_divergence') ) then

    ! Create the divergence diagnostic
    call divergence_diagnostic_alg( div_field, l2_norm, u_field, mesh )

    write( log_scratch_space, '(A,E16.8)' )  &
         'L2 of divergence =',l2_norm
    call log_event( log_scratch_space, LOG_LEVEL_INFO )

    if (use_xios_io) then
      !If using XIOS, we need to set a field I/O method appropriately
      tmp_write_ptr => write_field_generic
      call div_field%set_write_behaviour(tmp_write_ptr)
    end if

    call write_scalar_diagnostic( 'divergence', div_field, &
                                  clock, mesh, .false. )

    nullify(tmp_write_ptr)

  end if

end subroutine write_divergence_diagnostic

!-------------------------------------------------------------------------------
!>  @brief    Handles hydrostatic balance diagnostic processing
!!
!!  @details  Handles hydrostatic balance diagnostic processing
!!
!!> @param[in] theta_field   The theta field
!!> @param[in] exner_field   The exner field
!!> @param[in] mesh          Mesh
!-------------------------------------------------------------------------------

subroutine write_hydbal_diagnostic( theta_field, moist_dyn_field, exner_field,  &
                                    mesh )

  use logging_config_mod, only: run_log_level, run_log_level_debug

  implicit none

  type(field_type), intent(in)    :: theta_field
  type(field_type), intent(in)    :: moist_dyn_field(num_moist_factors)
  type(field_type), intent(in)    :: exner_field
  type(mesh_type),  intent(in), pointer :: mesh

  real(r_def)                     :: l2_norm = 0.0_r_def

  if (run_log_level == run_log_level_debug) then
    call hydbal_diagnostic_alg(l2_norm, theta_field, moist_dyn_field,          &
                               exner_field, mesh)

    write( log_scratch_space, '(A,E16.8)' )  &
        'L2 of hydrostatic imbalance =', l2_norm
    call log_event( log_scratch_space, LOG_LEVEL_DEBUG )
  end if

end subroutine write_hydbal_diagnostic


!-------------------------------------------------------------------------------
!>  @brief    Handles vorticity diagnostic processing.
!>  @details  Calculates the 3D vorticity diagnostic and optionally outputs it.
!!            Then calculates the vorticity on pressure levels diagnostic and
!!            optionally outputs it.
!> @param[in] u_field   The wind field
!> @param[in] exner     Exner pressure
!> @param[in] clock     Model clock object
!-------------------------------------------------------------------------------

subroutine write_vorticity_diagnostic(u_field, exner, clock)
#ifdef UM_PHYSICS
  use pres_lev_diags_alg_mod, only: pres_lev_field_alg
#endif
  implicit none

  type(field_type),        intent(in) :: u_field, exner
  class(model_clock_type), intent(in) :: clock

  type(field_type) :: vorticity, projected_field(3)
  type(field_type), pointer :: chi(:)
  type(field_type), pointer :: panel_id
  type(mesh_type),  pointer :: mesh
#ifdef UM_PHYSICS
  type(field_type) :: plev_xi3
#endif
  procedure(write_interface), pointer  :: tmp_write_ptr

  integer(i_def) :: output_dim, i
  character(len=1)                :: uchar
  character(len=str_max_filename) :: field_name_new

  logical(l_def) :: xi_modlev_flag, plev_xi3_flag
  logical(l_def), parameter :: xi3_axis = .true.

  xi_modlev_flag = diagnostic_to_be_sampled('xi1') .or. &
                   diagnostic_to_be_sampled('xi2') .or. &
                   diagnostic_to_be_sampled('xi3') .or. &
                   diagnostic_to_be_sampled('init_xi1') .or. &
                   diagnostic_to_be_sampled('init_xi2') .or. &
                   diagnostic_to_be_sampled('init_xi3')

#ifdef UM_PHYSICS
  plev_xi3_flag = init_diag(plev_xi3, 'plev__xi3')
#else
  plev_xi3_flag = .false.
#endif

  if (xi_modlev_flag .or. plev_xi3_flag) then

    call vorticity_diagnostic_alg(vorticity, u_field)

    if (use_xios_io) then

      output_dim = 3_i_def

      mesh     => vorticity%get_mesh()
      chi      => get_coordinates( mesh%get_id() )
      panel_id => get_panel_id( mesh%get_id() )

      ! Project the field to the output field
      call project_output( vorticity, projected_field, chi, panel_id, W3 )

#ifdef UM_PHYSICS
      if (plev_xi3_flag) then
        call pres_lev_field_alg(projected_field(3), exner, plev_xi3, xi3_axis)
      end if
#endif

      if (xi_modlev_flag) then

        ! Set up correct I/O handler for Xi projected to W3
        tmp_write_ptr => write_field_generic

        do i = 1, output_dim

          call projected_field(i)%set_write_behaviour(tmp_write_ptr)

          ! Write the component number into a new field name
          write(uchar,'(i1)') i
          field_name_new = trim('xi'//uchar)

          ! Check if we need to write an initial field
          if (clock%is_initialisation()) then
            call projected_field(i)%write_field(trim('init_'//field_name_new))
          else
            call projected_field(i)%write_field(trim(field_name_new))
          end if

        end do

      end if

    else

      call write_vector_diagnostic('xi', vorticity, clock, &
                                 vorticity%get_mesh(), nodal_output_on_w3)

    end if

  end if

end subroutine write_vorticity_diagnostic

#ifdef UM_PHYSICS
!-------------------------------------------------------------------------------
!> @brief     Potential vorticity diagnostic processing and output.
!> @details   Optionally calculate and output both model level and pressure
!>            level diagnostics.
!> @param[in] u_field   The wind field
!> @param[in] theta     The potential temperature field (K)
!> @param[in] rho       The density field (kg/m3)
!> @param[in] exner     The exner pressure (Pa)
!> @param[in] clock     The model clock object
!-------------------------------------------------------------------------------
subroutine write_pv_diagnostic(u_field, theta, rho, exner, clock)

  use pres_lev_diags_alg_mod, only: pres_lev_field_alg

  implicit none

  type(field_type),        intent(in) :: u_field
  type(field_type),        intent(in) :: theta
  type(field_type),        intent(in) :: rho
  type(field_type),        intent(in) :: exner
  class(model_clock_type), intent(in) :: clock

  type(field_type) :: pv
  type(field_type) :: plev_pv
  logical(l_def) :: pv_modlev_flag, plev_pv_flag
  logical(l_def), parameter :: xi3_axis = .false.
  logical(l_def), parameter :: add_W3_version = .false.

  pv_modlev_flag = diagnostic_to_be_sampled('potential_vorticity') .or. &
                   diagnostic_to_be_sampled('init_potential_vorticity')

  plev_pv_flag = init_diag(plev_pv, 'plev__pv')

  if (pv_modlev_flag .or. plev_pv_flag) then

    call potential_vorticity_diagnostic_alg(pv, u_field, theta, rho)

    if (plev_pv_flag) call pres_lev_field_alg(pv, exner, plev_pv, xi3_axis)

    if (pv_modlev_flag) then
      call write_scalar_diagnostic('potential_vorticity', pv, clock, &
                                    pv%get_mesh(), add_W3_version)
    end if

  end if

end subroutine write_pv_diagnostic
#else
!-------------------------------------------------------------------------------
!> @brief     Potential vorticity diagnostic processing and output.
!> @details   Optionally calculate and output model level diagnostic.
!> @param[in] u_field   The wind field
!> @param[in] theta     The potential temperature field (K)
!> @param[in] rho       The density field (kg/m3)
!> @param[in] clock     The model clock object
!-------------------------------------------------------------------------------
subroutine write_pv_diagnostic(u_field, theta, rho, clock)

  implicit none

  type(field_type),        intent(in) :: u_field
  type(field_type),        intent(in) :: theta
  type(field_type),        intent(in) :: rho
  class(model_clock_type), intent(in) :: clock

  type(field_type) :: pv
  logical(l_def) :: pv_modlev_flag
  logical(l_def), parameter :: add_W3_version = .false.

  pv_modlev_flag = diagnostic_to_be_sampled('potential_vorticity') .or. &
                   diagnostic_to_be_sampled('init_potential_vorticity')

  if (pv_modlev_flag) then

    call potential_vorticity_diagnostic_alg(pv, u_field, theta, rho)

    call write_scalar_diagnostic('potential_vorticity', pv, clock, &
                                  pv%get_mesh(), add_W3_version)

  end if

end subroutine write_pv_diagnostic
#endif
end module diagnostics_calc_mod
