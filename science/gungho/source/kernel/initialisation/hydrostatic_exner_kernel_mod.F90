!-----------------------------------------------------------------------------
! Copyright (c) 2018,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Computes Exner pressure distribution via hydrostatic balance equation (lowest order only)

module hydrostatic_exner_kernel_mod

use argument_mod,               only : arg_type, func_type,       &
                                       GH_FIELD, GH_REAL,         &
                                       GH_SCALAR,                 &
                                       GH_READ, GH_WRITE,         &
                                       ANY_SPACE_9, GH_BASIS,     &
                                       ANY_DISCONTINUOUS_SPACE_3, &
                                       CELL_COLUMN, GH_EVALUATOR
use constants_mod,              only : r_def, i_def
use idealised_config_mod,       only : test
use kernel_mod,                 only : kernel_type
use fs_continuity_mod,          only : Wtheta, W3
use formulation_config_mod,     only : init_exner_bt, shallow

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: hydrostatic_exner_kernel_type
  private
  type(arg_type) :: meta_args(12) = (/                                     &
       arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),                        &
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta),                    &
       arg_type(GH_FIELD*3, GH_REAL, GH_READ,  Wtheta),                    &
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta),                    &
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),                        &
       arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),               &
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
       arg_type(GH_SCALAR,  GH_REAL, GH_READ),                             &
       arg_type(GH_SCALAR,  GH_REAL, GH_READ),                             &
       arg_type(GH_SCALAR,  GH_REAL, GH_READ),                             &
       arg_type(GH_SCALAR,  GH_REAL, GH_READ),                             &
       arg_type(GH_SCALAR,  GH_REAL, GH_READ)                              &
       /)

  type(func_type) :: meta_funcs(1) = (/                                    &
       func_type(ANY_SPACE_9, GH_BASIS)                                    &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
  integer :: gh_evaluator_targets(1) = (/ Wtheta /)
contains
  procedure, nopass :: hydrostatic_exner_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: hydrostatic_exner_code

contains

!> @brief Computes hydrostatic Exner function
!> @param[in] nlayers Number of layers
!> @param[in,out] exner Exner pressure field
!> @param[in] theta Potential temperature field
!> @param[in] moist_dyn_gas Moist dynamics factor in gas law
!> @param[in] moist_dyn_tot Moist dynamics total mass factor
!> @param[in] moist_dyn_fac Moist dynamics water factor
!> @param[in] height_wt Height coordinate in Wtheta
!> @param[in] height_w3 Height coordinate in W3
!> @param[in] chi_1 First component of the chi coordinate field
!> @param[in] chi_2 Second component of the chi coordinate field
!> @param[in] chi_3 Third component of the chi coordinate field
!> @param[in] panel_id A field giving the ID for mesh panels
!> @param[in] gravity The planet gravity
!> @param[in] p_zero Reference surface pressure
!> @param[in] rd Gas constant for dry air
!> @param[in] cp Specific heat of dry air at constant pressure
!> @param[in] radius Radius of planet at surface
!> @param[in] ndf_w3 Number of degrees of freedom per cell for W3
!> @param[in] undf_w3 Number of unique degrees of freedom for W3
!> @param[in] map_w3 Dofmap for the cell at the base of the column for W3
!> @param[in] ndf_wt Number of degrees of freedom per cell for Wtheta
!> @param[in] undf_wt Number of unique degrees of freedom  for Wtheta
!> @param[in] map_wt Dofmap for the cell at the base of the column for Wtheta
!> @param[in] ndf_chi Number of degrees of freedom per cell for Wchi
!> @param[in] undf_chi Number of unique degrees of freedom  for Wchi
!> @param[in] map_chi Dofmap for the cell at the base of the column for Wchi
!> @param[in] basis_chi_on_wt Basis functions for Wchi evaluated at
!!                            Wtheta nodal points
!> @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!> @param[in] undf_pid Number of unique degrees of freedom for panel_id
!> @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
subroutine hydrostatic_exner_code(nlayers, exner, theta,         &
                                  moist_dyn_gas, moist_dyn_tot,  &
                                  moist_dyn_fac,                 &
                                  height_wt, height_w3,          &
                                  chi_1, chi_2, chi_3, panel_id, &
                                  gravity, p_zero, rd, cp,       &
                                  radius,                        &
                                  ndf_w3, undf_w3, map_w3,       &
                                  ndf_wt, undf_wt, map_wt,       &
                                  ndf_chi, undf_chi, map_chi,    &
                                  basis_chi_on_wt,               &
                                  ndf_pid, undf_pid, map_pid     )

  use analytic_pressure_profiles_mod, only : analytic_pressure
  use sci_chi_transform_mod,          only : chi2xyz

  implicit none

  ! Arguments
  integer(kind=i_def),                              intent(in)    :: nlayers
  integer(kind=i_def),                              intent(in)    :: ndf_w3, undf_w3
  integer(kind=i_def),                              intent(in)    :: ndf_wt, undf_wt
  integer(kind=i_def),                              intent(in)    :: ndf_chi, undf_chi
  integer(kind=i_def),                              intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), dimension(ndf_w3),           intent(in)    :: map_w3
  integer(kind=i_def), dimension(ndf_wt),           intent(in)    :: map_wt
  integer(kind=i_def), dimension(ndf_chi),          intent(in)    :: map_chi
  integer(kind=i_def), dimension(ndf_pid),          intent(in)    :: map_pid
  real(kind=r_def),    dimension(undf_w3),          intent(inout) :: exner
  real(kind=r_def),    dimension(undf_wt),          intent(in)    :: theta
  real(kind=r_def),    dimension(undf_wt),          intent(in)    :: moist_dyn_gas, &
                                                                     moist_dyn_tot, &
                                                                     moist_dyn_fac
  real(kind=r_def),    dimension(undf_wt),          intent(in)    :: height_wt
  real(kind=r_def),    dimension(undf_w3),          intent(in)    :: height_w3
  real(kind=r_def),    dimension(undf_chi),         intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def),    dimension(undf_pid),         intent(in)    :: panel_id
  real(kind=r_def),    dimension(1,ndf_chi,ndf_wt), intent(in)    :: basis_chi_on_wt
  real(kind=r_def),                                 intent(in)    :: gravity
  real(kind=r_def),                                 intent(in)    :: p_zero
  real(kind=r_def),                                 intent(in)    :: rd
  real(kind=r_def),                                 intent(in)    :: cp
  real(kind=r_def),                                 intent(in)    :: radius

  ! Internal variables
  integer(kind=i_def)                  :: k, dfc, layers_offset, wt_dof, ipanel
  real(kind=r_def)                     :: dz
  real(kind=r_def)                     :: theta_moist
  real(kind=r_def)                     :: g_local(0:nlayers-1)
  real(kind=r_def)                     :: coords(3), xyz(3)
  real(kind=r_def)                     :: exner_start
  real(kind=r_def), dimension(ndf_chi) :: chi_1_e, chi_2_e, chi_3_e

  ipanel = int(panel_id(map_pid(1)), i_def)

  if (init_exner_bt) then
    layers_offset = 0
    wt_dof = 1
  else
    layers_offset = nlayers - 1
    wt_dof = 2
  end if

  do dfc = 1, ndf_chi
    chi_1_e(dfc) = chi_1( map_chi(dfc) + layers_offset )
    chi_2_e(dfc) = chi_2( map_chi(dfc) + layers_offset )
    chi_3_e(dfc) = chi_3( map_chi(dfc) + layers_offset )
  end do

  ! Horizontal coordinates of cell bottom or top
  coords(:) = 0.0_r_def
  do dfc = 1, ndf_chi
    coords(1) = coords(1) + chi_1_e(dfc)*basis_chi_on_wt(1,dfc,wt_dof)
    coords(2) = coords(2) + chi_2_e(dfc)*basis_chi_on_wt(1,dfc,wt_dof)
    coords(3) = coords(3) + chi_3_e(dfc)*basis_chi_on_wt(1,dfc,wt_dof)
  end do

  call chi2xyz(coords(1), coords(2), coords(3), ipanel, xyz(1), xyz(2), xyz(3))

  ! Exner at the model surface or top
  exner_start = analytic_pressure( xyz, test, 0.0_r_def)

  ! Set local value of gravity at each vertical level
  if (shallow) then
    do k = 0, nlayers-1
      g_local(k) = gravity
    end do
  else
    do k = 0, nlayers-1
      g_local(k) = gravity * (radius / (radius + height_wt(map_wt(1)+k))) ** 2
    end do
  end if

  if (init_exner_bt) then

    ! Bottom-up initialization
    ! Exner at the bottom level
    dz = height_w3(map_w3(1))-height_wt(map_wt(1))
    theta_moist = moist_dyn_gas(map_wt(1)) * theta(map_wt(1)) /   &
                  moist_dyn_tot(map_wt(1))
    exner(map_w3(1)) = exner_start - g_local(0) * dz / (cp * theta_moist)

    ! Exner on other levels
    do k = 1, nlayers-1
      dz = height_w3(map_w3(1)+k)-height_w3(map_w3(1)+k-1)
      theta_moist = moist_dyn_gas(map_wt(1)+k) * theta(map_wt(1)+k) /   &
                    moist_dyn_tot(map_wt(1)+k)
      exner(map_w3(1)+k) = exner(map_w3(1)+k-1) - g_local(k) * dz / (cp * theta_moist)
    end do

  else

    ! Top-down initialization
    ! Exner at the top level
    dz = height_wt(map_wt(1)+nlayers) - height_w3(map_w3(1)+nlayers-1)
    theta_moist = moist_dyn_gas(map_wt(1)+nlayers) * theta(map_wt(1)+nlayers) / &
                  moist_dyn_tot(map_wt(1)+nlayers)
    exner(map_w3(1)+nlayers-1) = exner_start + g_local(nlayers-1) * dz / (cp * theta_moist)

    ! Exner on other levels
    do k = nlayers-1, 1, -1
      dz = height_w3(map_w3(1)+k)-height_w3(map_w3(1)+k-1)
      theta_moist = moist_dyn_gas(map_wt(1)+k) * theta(map_wt(1)+k) /   &
                    moist_dyn_tot(map_wt(1)+k)
      exner(map_w3(1)+k-1) = exner(map_w3(1)+k) + g_local(k-1) * dz / (cp * theta_moist)
    end do

  end if

end subroutine hydrostatic_exner_code

end module hydrostatic_exner_kernel_mod
