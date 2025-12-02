!----------------------------------------------------------------------------
! (c) Crown copyright 2020-2024 Met Office. All rights reserved.
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
    use gen_phys_inputs_mod, only: l_mr_physics, l_vol_interp_rho
    use model_domain_mod, only: model_type, mt_lfric
    use nlsizes_namelist_mod, only: model_levels, n_cca_lev

    implicit none

    ! ----------------------------------------------------------------
    ! Model dimensions - contained in UM module nlsizes_namelist_mod
    ! ----------------------------------------------------------------
    ! Vertical dimensions set from LFRic number of layers.
    model_levels = number_of_layers
    n_cca_lev    = number_of_layers

    ! ----------------------------------------------------------------
    ! Model type - contained in UM module model_domain_mod
    ! ----------------------------------------------------------------
    ! Previous usage of mt_single_column has been removed to purge
    ! ambiguous use in UM codebase
    model_type = mt_lfric

    ! Mixing ratio flag - contained in UM gen_phys_inputs.
    ! Set to true as LFRic only supports mixing ratios.
    l_mr_physics = .true.

    ! This is only used in microphysics air density calculations
    ! It should be false for consistency with how the dynamics maps
    ! between meshes
    l_vol_interp_rho = .false.

  end subroutine um_control_init

end module um_control_init_mod
