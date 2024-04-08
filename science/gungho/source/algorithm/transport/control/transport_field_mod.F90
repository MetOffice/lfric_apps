!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Contains central routine for transporting fields.
!> @details Contains routine to transport a field pointing to
!!          particular routines based on the specified transport options.

module transport_field_mod

  use constants_mod,                    only: r_def, r_tran, i_def
  use field_mod,                        only: field_type
  use r_tran_field_mod,                 only: r_tran_field_type
  use log_mod,                          only: log_event, LOG_LEVEL_ERROR
  use ffsl_control_alg_mod,             only: ffsl_control
  use split_transport_mod,              only: split_transport_control
  use transport_metadata_mod,           only: transport_metadata_type
  use transport_runtime_alg_mod,        only: transport_runtime_type
  use transport_runtime_collection_mod, only: get_transport_runtime
  use transport_enumerated_types_mod,   only: scheme_mol_3d,              &
                                              scheme_ffsl_3d,             &
                                              scheme_split,               &
                                              direction_3d,               &
                                              equation_form_conservative, &
                                              equation_form_advective,    &
                                              equation_form_consistent
  use mol_conservative_alg_mod,         only: mol_conservative_alg
  use mol_advective_alg_mod,            only: mol_advective_alg
  use mol_consistent_alg_mod,           only: mol_consistent_alg
  use psykal_lite_mod,                  only: invoke_copy_rtran_to_rdef, &
                                              invoke_copy_to_rtran

  implicit none

  private

  public :: transport_field, transport_field_r_tran


contains

  !> @brief Central routine for transporting fields.
  !> @details Performs a whole time step, solving the transport equation for
  !!          a field.
  !> @param[in,out] field_np1_rdef     Field to return at end of transport step
  !> @param[in]     field_n_rdef       Field at the start of the transport step
  !> @param[in]     model_dt_rdef      Time difference across time step
  !> @param[in]     transport_metadata Contains the configuration options for
  !!                                   transporting these fields

  subroutine transport_field(field_np1_rdef, field_n_rdef, model_dt_rdef, transport_metadata)

    implicit none

    ! Arguments
    type(field_type),              intent(inout) :: field_np1_rdef
    type(field_type),              intent(in)    :: field_n_rdef
    real(kind=r_def),              intent(in)    :: model_dt_rdef
    type(transport_metadata_type), intent(in)    :: transport_metadata

    ! Local conversions to the r_tran precision
    type(r_tran_field_type) :: field_np1
    type(r_tran_field_type) :: field_n
    real(kind=r_tran)       :: model_dt

    ! transfer r_def input to r_tran fields
    call field_np1%initialise( vector_space = field_np1_rdef%get_function_space()  )
    call field_n%initialise( vector_space = field_n_rdef%get_function_space()  )
    call invoke_copy_to_rtran( field_n, field_n_rdef )
    model_dt = real( model_dt_rdef, r_tran )

    call transport_field_r_tran(field_np1, field_n, model_dt, transport_metadata)

    ! Transfer back to r_def output
    call invoke_copy_rtran_to_rdef( field_np1_rdef, field_np1 )

  end subroutine transport_field


  subroutine transport_field_r_tran(field_np1, field_n, model_dt, transport_metadata)

    implicit none

    ! Arguments
    type(r_tran_field_type),       intent(inout) :: field_np1
    type(r_tran_field_type),       intent(in)    :: field_n
    real(kind=r_tran),             intent(in)    :: model_dt
    type(transport_metadata_type), intent(in)    :: transport_metadata
    type(transport_runtime_type),  pointer       :: transport_runtime => null()

    ! Reset the counter for tracer transport steps and store nth level field
    transport_runtime => get_transport_runtime(field_n%get_mesh())
    call transport_runtime%reset_tracer_step_ctr()
    call transport_runtime%set_field_n(field_n)
    nullify( transport_runtime )

    ! First choose scheme, and for full 3D schemes then choose equation
    select case ( transport_metadata%get_scheme() )

    ! -------------------------------------------------------------------------!
    ! Full 3D Method of Lines scheme
    ! -------------------------------------------------------------------------!
    case ( scheme_mol_3d )
      ! Choose form of transport equation
      select case ( transport_metadata%get_equation_form() )
      case ( equation_form_conservative )
         call mol_conservative_alg(field_np1, field_n, &
                                   direction_3d, transport_metadata)

      case ( equation_form_advective )
         call mol_advective_alg(field_np1, field_n, &
                                direction_3d, transport_metadata)

      case ( equation_form_consistent )
         call mol_consistent_alg(field_np1, field_n, &
                                 direction_3d, transport_metadata)

      case default
        call log_event('Trying to solve unrecognised form of transport equation', &
                        LOG_LEVEL_ERROR)

      end select

    ! -------------------------------------------------------------------------!
    ! Full 3D Flux-Form Semi-Lagrangian scheme
    ! -------------------------------------------------------------------------!
    case ( scheme_ffsl_3d )
      ! Choose form of transport equation
      select case ( transport_metadata%get_equation_form() )
      case ( equation_form_conservative, equation_form_advective )
        call ffsl_control(field_np1, field_n, direction_3d, &
                          model_dt, transport_metadata)

      case ( equation_form_consistent )
        call log_event('Consistent 3D FFSL not implemented', LOG_LEVEL_ERROR)

      case default
        call log_event('Trying to solve unrecognised form of transport equation', &
                        LOG_LEVEL_ERROR)

      end select

    ! -------------------------------------------------------------------------!
    ! Some split horizontal/vertical transport scheme
    ! -------------------------------------------------------------------------!
    case ( scheme_split )
      call split_transport_control(field_np1, field_n, model_dt, &
                                   transport_metadata)

    case default
      call log_event('Trying to transport with unrecognised scheme', &
                      LOG_LEVEL_ERROR)

    end select

  end subroutine transport_field_r_tran



end module transport_field_mod
