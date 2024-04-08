!----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!----------------------------------------------------------------------------
!> @brief Controls the setting of UM high level variables which are either
!>         fixed in LFRic or derived from LFRic inputs.

module um_control_init_mod

  ! LFRic namelists which have bee read
  use extrusion_config_mod,        only : number_of_layers
  use timestepping_config_mod,     only : outer_iterations

  implicit none

  private
  public :: um_control_init

contains

  !>@brief Initialise UM high levels variables with are either fixed in LFRic
  !>        or derived from LFRic inputs.
  !>@details Nothing in this file is ever likely to be promoted to the LFRic
  !>          namelist. Everything is either set from an LFRic variable
  !>          already in the namelist, or is "fixed" from the perspective
  !>          of LFRic (but cannot be made a parameter because it is required
  !>          to be variable in the UM and therefore declared as such in the
  !>          UM modules which contains it).
  subroutine um_control_init()

    ! UM modules containing things that need setting
    use dynamics_input_mod, only: numcycles
    use gen_phys_inputs_mod, only: l_mr_physics
    use model_domain_mod, only: model_type, mt_single_column
    use nlsizes_namelist_mod, only: model_levels, cloud_levels, n_cca_lev, &
         tr_vars

    implicit none

    ! ----------------------------------------------------------------
    ! Model dimensions - contained in UM module nlsizes_namelist_mod
    ! ----------------------------------------------------------------
    ! Vertical dimensions set from LFRic number of layers.
    model_levels = number_of_layers
    cloud_levels = number_of_layers
    n_cca_lev    = number_of_layers
    ! Number of tracer variables - not implemented in LFRic so fixed at 0.
    ! This may change if we implement tracers ever.
    tr_vars = 0

    ! ----------------------------------------------------------------
    ! Model type - contained in UM module model_domain_mod
    ! ----------------------------------------------------------------
    ! Currently set to single column, as the physics kernels are called 1
    ! cell at a time. This may change if we pass multiple cells to kernels.
    model_type = mt_single_column

    ! Number of outer iterations of the dynamics - contained in UM
    ! dynamics_input_mod.
    ! Set from LFRic input number of outer iterations.
    numcycles = outer_iterations

    ! Mixing ratio flag - contained in UM gen_phys_inputs.
    ! Set to true as LFRic only supports mixing ratios.
    l_mr_physics = .true.

  end subroutine um_control_init

end module um_control_init_mod
