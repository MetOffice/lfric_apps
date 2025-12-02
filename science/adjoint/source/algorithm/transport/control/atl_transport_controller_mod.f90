!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Adjoint of the tl_transport_controller initialise routine
module atl_transport_controller_mod

  use adj_transport_controller_mod, only: adj_ls_wind_pert_rho_initialiser, &
                                          adj_pert_wind_ls_rho_initialiser
  use field_mod,                    only: field_type
  use model_clock_mod,              only: model_clock_type
  use tl_transport_controller_mod,  only: tl_transport_controller_type
  use transport_controller_mod,     only: transport_controller_type

  implicit none

  private
  public :: atl_transport_controller_initialiser

  contains
  !=============================================================================
  !> @brief Adjoint of the initialisation of the tl_transport_controller object
  !> @param[in,out] tl_transport_controller  The tangent linear transport controller object
  !> @param[in]     model_clock              Tracks the time within the model
  !> @param[in,out] ref_field_rdef           Perturbed reference density to use for
  !!                                         conservative transport of tracers.
  !> @param[in,out] wind_n_rdef              Perturbed wind at n-th time level.
  !> @param[in,out] wind_np1_rdef            Perturbed wind at (n+1)-th time level.
  subroutine atl_transport_controller_initialiser( tl_transport_controller, &
                                                   model_clock,             &
                                                   ref_field_rdef,          &
                                                   wind_n_rdef,             &
                                                   wind_np1_rdef )
    implicit none

    ! Arguments
    type(tl_transport_controller_type), intent(inout) :: tl_transport_controller
    type(model_clock_type),                intent(in) :: model_clock
    type(field_type),                   intent(inout) :: ref_field_rdef
    type(field_type),                   intent(inout) :: wind_n_rdef
    type(field_type),                   intent(inout) :: wind_np1_rdef

    ! Internal arguments
    type(transport_controller_type),          pointer :: pert_wind_ls_rho_controller
    type(transport_controller_type),          pointer :: ls_wind_pert_rho_controller

    pert_wind_ls_rho_controller => tl_transport_controller%get_pert_wind_ls_rho_controller()
    ls_wind_pert_rho_controller => tl_transport_controller%get_ls_wind_pert_rho_controller()

    call adj_pert_wind_ls_rho_initialiser( pert_wind_ls_rho_controller, &
                                           model_clock,                 &
                                           wind_n_rdef,                 &
                                           wind_np1_rdef )
    call adj_ls_wind_pert_rho_initialiser( ls_wind_pert_rho_controller, &
                                           ref_field_rdef )

  end subroutine atl_transport_controller_initialiser

end module atl_transport_controller_mod
