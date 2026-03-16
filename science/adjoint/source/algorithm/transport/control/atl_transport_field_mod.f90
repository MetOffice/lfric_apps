!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Contains adjoint of central routine for transporting fields.

module atl_transport_field_mod

  use r_tran_field_mod,                 only: r_tran_field_type
  use constants_mod,                    only: r_def, i_def, r_tran, l_def
  use field_mod,                        only: field_type
  use log_mod,                          only: log_event, LOG_LEVEL_ERROR
  use transport_metadata_mod,           only: transport_metadata_type
  use transport_enumerated_types_mod,   only: scheme_mol_3d,              &
                                              scheme_ffsl_3d,             &
                                              scheme_split,               &
                                              direction_3d,               &
                                              equation_form_conservative, &
                                              equation_form_advective
  use adj_trans_lookup_cache_mod,       only: adj_trans_lookup_cache_type
  use atl_mol_conservative_alg_mod,     only: atl_mol_conservative_alg
  use atl_mol_advective_alg_mod,        only: atl_mol_advective_alg
  use atl_split_transport_mod,          only: atl_split_transport_control
  use transport_controller_mod,         only: transport_controller_type
  use tl_transport_controller_mod,      only: tl_transport_controller_type
  use transport_counter_mod,            only: transport_counter_type

  implicit none

  private

  public :: atl_transport_field

contains

  !=============================================================================
  !> @brief Adjoint of central routine for transporting fields in tangent-linear field.
  !> @param[in,out] field_np1                ACTIVE  Field to return at end of transport step
  !> @param[in,out] field_n                  ACTIVE  Field at the start of the transport step
  !> @param[in]     ls_field_n               PASSIVE Linear field at the start of step
  !> @param[in,out] tl_transport_controller  Object controlling transport by perturbed wind
  !> @param[in]     transport_metadata       Contains the configuration options for
  !!                                         transporting these fields
  !> @param[in]     adj_lookup_table_cache   Lookup table cache
  subroutine atl_transport_field(field_np1, field_n, ls_field_n,              &
                                 tl_transport_controller, transport_metadata, &
                                 adj_lookup_table_cache)

    implicit none

    ! Arguments
    type(field_type),                   intent(inout) :: field_np1
    type(field_type),                   intent(inout) :: field_n
    type(field_type),                   intent(in)    :: ls_field_n
    type(tl_transport_controller_type), intent(inout) :: tl_transport_controller
    type(transport_metadata_type),      intent(inout) :: transport_metadata
    type(adj_trans_lookup_cache_type),  intent(in)    :: adj_lookup_table_cache

    ! Internal variables
    type(transport_counter_type),    pointer :: transport_counter
    type(transport_controller_type), pointer :: transport_controller

    ! Initialise the counter and set metadata in the transport controller
    ! This also updates the transport metadata, depending on outer loop
    call tl_transport_controller%before_transport_field(transport_metadata)

    ! Get the transport counter, and set field_n for each of them
    ! TODO: this may not be strictly necessary, as we may be able to set it for
    ! a single transport counter, but it is the safest thing for now
    transport_controller => tl_transport_controller%get_ls_wind_ls_rho_controller()
    transport_counter => transport_controller%get_transport_counter()
    call transport_counter%set_field_n(field_n)

    transport_controller => tl_transport_controller%get_ls_wind_pert_rho_controller()
    transport_counter => transport_controller%get_transport_counter()
    call transport_counter%set_field_n(field_n)

    transport_controller => tl_transport_controller%get_pert_wind_ls_rho_controller()
    transport_counter => transport_controller%get_transport_counter()
    call transport_counter%set_field_n(field_n)

    ! First choose scheme, and for full 3D schemes then choose equation
    select case ( transport_metadata%get_scheme() )

    ! -------------------------------------------------------------------------!
    ! Full 3D Method of Lines scheme
    ! -------------------------------------------------------------------------!
    case ( scheme_mol_3d )
      ! Choose form of transport equation
      select case ( transport_metadata%get_equation_form() )
      case ( equation_form_conservative )
        call atl_mol_conservative_alg(                                                          &
                field_np1, field_n, ls_field_n, tl_transport_controller, adj_lookup_table_cache &
        )

      case ( equation_form_advective )
        call atl_mol_advective_alg(                                                             &
                field_np1, field_n, ls_field_n, tl_transport_controller, adj_lookup_table_cache &
        )

      case default
        call log_event('Trying to solve unrecognised form of transport equation', &
                        LOG_LEVEL_ERROR)

      end select

    ! -------------------------------------------------------------------------!
    ! Full 3D Flux-Form Semi-Lagrangian scheme
    ! -------------------------------------------------------------------------!
    case ( scheme_ffsl_3d )
      call log_event('FFSL not implemented for tangent-linear transport', &
                      LOG_LEVEL_ERROR)

    ! -------------------------------------------------------------------------!
    ! Some split horizontal/vertical transport scheme
    ! -------------------------------------------------------------------------!
    case ( scheme_split )
      call atl_split_transport_control(field_np1, field_n, ls_field_n, &
                                       tl_transport_controller, adj_lookup_table_cache)

    case default
      call log_event('Trying to transport with unrecognised scheme', &
                      LOG_LEVEL_ERROR)

    end select

    ! Reset any metadata after transporting the field, if necessary
    call tl_transport_controller%after_transport_field()

  end subroutine atl_transport_field

end module atl_transport_field_mod
