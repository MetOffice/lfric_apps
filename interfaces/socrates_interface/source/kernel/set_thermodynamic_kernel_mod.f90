!-------------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Setup the thermodynamic fields for the Socrates radiation calls
module set_thermodynamic_kernel_mod

  use argument_mod, only : arg_type,            &
                           GH_FIELD, GH_SCALAR, &
                           GH_REAL, GH_READ,    &
                           GH_WRITE,            &
                           CELL_COLUMN,         &
                           ANY_DISCONTINUOUS_SPACE_1
  use fs_continuity_mod, only : W3, Wtheta
  use constants_mod,     only : r_def, i_def
  use kernel_mod,        only : kernel_type

  implicit none

  private

  type, public, extends(kernel_type) :: set_thermodynamic_kernel_type
    private
    type(arg_type) :: meta_args(15) = (/ &
      arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3                       ), & ! exner
      arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta                   ), & ! exner_in_wth
      arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta                   ), & ! theta_in_wth
      arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3                       ), & ! theta_in_w3
      arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta                   ), & ! rho_in_wth
      arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta                   ), & ! dz_in_wth
      arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta                   ), & ! temperature_wth
      arg_type(GH_FIELD,  GH_REAL, GH_READ,  Wtheta                   ), & ! moist_mass_fac
      arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta                   ), & ! d_mass
      arg_type(GH_FIELD,  GH_REAL, GH_WRITE, Wtheta                   ), & ! layer_heat_capacity
      arg_type(GH_FIELD,  GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! t_layer_boundaries
      arg_type(GH_SCALAR, GH_REAL, GH_READ                            ), & ! p_zero
      arg_type(GH_SCALAR, GH_REAL, GH_READ                            ), & ! kappa
      arg_type(GH_SCALAR, GH_REAL, GH_READ                            ), & ! gravity
      arg_type(GH_SCALAR, GH_REAL, GH_READ                            )  & ! cp
      /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: set_thermodynamic_code
  end type set_thermodynamic_kernel_type

  public :: set_thermodynamic_code

contains
  !> @param[in]    nlayers                    Number of layers
  !> @param[in]    exner_w3                   Exner pressure field in density space
  !> @param[in]    exner_wtheta               Exner pressure field in potential temperature space
  !> @param[in]    theta_wtheta               Potential temperature field
  !> @param[in]    theta_w3                   Potential temperature field in density space
  !> @param[in]    rho_wtheta                 Density field in potential temperature space
  !> @param[in]    dz_wtheta                  Depth of temperature space levels
  !> @param[in]    temperature_wth            Temperature in Wtheta space
  !> @param[in]    moist_mass_fac             1 + total water content
  !> @param[inout] d_mass_wtheta              Mass per square metre of Socrates layers
  !> @param[inout] layer_heat_capacity_wtheta Heat capacity of Socrates layers
  !> @param[inout] t_layer_boundaries         Temperature at layer boundaries
  !> @param[in]    p_zero                     Reference surface pressure
  !> @param[in]    kappa                      Ratio of gas constant and specific heat
  !> @param[in]    gravity                    Gravity
  !> @param[in]    cp                         Specific heat of dry air
  !> @param[in]    ndf_w3                     No. DOFs per cell for W3 space
  !> @param[in]    undf_w3                    No. unique DOFs for W3 space
  !> @param[in]    map_w3                     Dofmap for W3 space column base cell
  !> @param[in]    ndf_wtheta                 No. DOFs per cell for Wtheta space
  !> @param[in]    undf_wtheta                No. unique DOFs for Wtheta space
  !> @param[in]    map_wtheta                 Dofmap for Wtheta space column base cell
  !> @param[in]    ndf_flux                   No. DOFs per cell for flux space
  !> @param[in]    undf_flux                  No. unique DOFs for flux space
  !> @param[in]    map_flux                   Dofmap for flux space column base cell
  subroutine set_thermodynamic_code( nlayers, &
                                     exner_w3, &
                                     exner_wtheta, &
                                     theta_wtheta, &
                                     theta_w3, &
                                     rho_wtheta, &
                                     dz_wtheta, &
                                     temperature_wth, &
                                     moist_mass_fac, &
                                     d_mass_wtheta, &
                                     layer_heat_capacity_wtheta, &
                                     t_layer_boundaries, &
                                     p_zero, &
                                     kappa, &
                                     gravity, &
                                     cp, &
                                     ndf_w3, &
                                     undf_w3, &
                                     map_w3, &
                                     ndf_wtheta, &
                                     undf_wtheta, &
                                     map_wtheta, &
                                     ndf_flux, &
                                     undf_flux, &
                                     map_flux )

    implicit none

    ! Dummy variables
    integer(i_def), intent(in)                            :: nlayers
    integer(i_def), intent(in)                            :: ndf_w3
    integer(i_def), intent(in)                            :: ndf_wtheta
    integer(i_def), intent(in)                            :: ndf_flux
    integer(i_def), intent(in)                            :: undf_w3
    integer(i_def), intent(in)                            :: undf_wtheta
    integer(i_def), intent(in)                            :: undf_flux
    integer(i_def), intent(in),    dimension(ndf_w3)      :: map_w3
    integer(i_def), intent(in),    dimension(ndf_wtheta)  :: map_wtheta
    integer(i_def), intent(in),    dimension(ndf_flux)    :: map_flux
    real(r_def),    intent(in),    dimension(undf_w3)     :: exner_w3
    real(r_def),    intent(in),    dimension(undf_wtheta) :: exner_wtheta
    real(r_def),    intent(in),    dimension(undf_wtheta) :: theta_wtheta
    real(r_def),    intent(in),    dimension(undf_w3)     :: theta_w3
    real(r_def),    intent(in),    dimension(undf_wtheta) :: rho_wtheta
    real(r_def),    intent(in),    dimension(undf_wtheta) :: dz_wtheta
    real(r_def),    intent(in),    dimension(undf_wtheta) :: temperature_wth
    real(r_def),    intent(in),    dimension(undf_wtheta) :: moist_mass_fac
    real(r_def),    intent(inout), dimension(undf_wtheta) :: d_mass_wtheta
    real(r_def),    intent(inout), dimension(undf_wtheta) :: layer_heat_capacity_wtheta
    real(r_def),    intent(inout), dimension(undf_flux)   :: t_layer_boundaries
    real(r_def),    intent(in) :: p_zero, kappa, gravity, cp

    ! Internal vaiables
    integer(i_def) :: k

    ! Set bottom level of theta fields to zero so they aren't get filled with NANs
    ! which can cause problems with later use of builtins on the fields. This issue
    ! could be avoided by having a proper extrusion and functionspace for these fields.
    d_mass_wtheta(map_wtheta(1)) = 0.0_r_def
    layer_heat_capacity_wtheta(map_wtheta(1)) = 0.0_r_def

    ! Calculate the dry mass as required when using mixing ratios.
    ! Mass of bottom layer bounded by the surface:
    d_mass_wtheta(map_wtheta(2)) = rho_wtheta(map_wtheta(2)) * (dz_wtheta(map_wtheta(2)) + dz_wtheta(map_wtheta(1)))
    do k=2,nlayers-1
      d_mass_wtheta(map_wtheta(1) + k) = rho_wtheta(map_wtheta(1)+k) * dz_wtheta(map_wtheta(1)+k)
    end do
    ! Hydrostatic approximation for mass of top layer:
    d_mass_wtheta(map_wtheta(2)+nlayers-1) = p_zero * exner_w3(map_w3(1) + nlayers-1)**(1.0_r_def/kappa) / gravity

    ! Heat capacity of layers
    do k = 1,nlayers
      layer_heat_capacity_wtheta(map_wtheta(1) + k) = &
       d_mass_wtheta(map_wtheta(1) + k) * cp * moist_mass_fac(map_wtheta(1) + k)
    end do

    ! Calculate temperature at layer boundaries
    t_layer_boundaries(map_flux(1)) = theta_wtheta(map_wtheta(1)) * exner_wtheta(map_wtheta(1))
    do k=1,nlayers-1
      t_layer_boundaries(map_flux(1) + k) = theta_w3(map_w3(1) + k) * exner_w3(map_w3(1) + k)
    end do
    t_layer_boundaries(map_flux(1) + nlayers) = temperature_wth(map_wtheta(1) + nlayers )

  end subroutine set_thermodynamic_code

end module set_thermodynamic_kernel_mod
