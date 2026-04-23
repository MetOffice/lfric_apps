!-----------------------------------------------------------------------------
! Copyright (c) 2026,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Computes vertical balance for simultaneous initialisation of
!>        temperature and vapour profiles

module vert_balance_kernel_mod

use argument_mod,           only: arg_type, func_type,       &
                                  GH_FIELD, GH_REAL,         &
                                  GH_SCALAR,                 &
                                  GH_READ, GH_WRITE,         &
                                  ANY_SPACE_9, GH_BASIS,     &
                                  ANY_DISCONTINUOUS_SPACE_3, &
                                  CELL_COLUMN, GH_EVALUATOR
use constants_mod,          only: r_def, i_def
use idealised_config_mod,   only: test
use initial_temperature_config_mod, &
                            only: temp_variable => profile_variable, &
                                  profile_variable_absolute,         &
                                  profile_variable_potential
use initial_vapour_config_mod, &
                            only: vapour_variable => profile_variable, &
                                  profile_variable_mr,                 &
                                  profile_variable_rh
use kernel_mod,             only: kernel_type
use extrusion_config_mod,   only: planet_radius
use fs_continuity_mod,      only: Wtheta, W3
use formulation_config_mod, only: shallow
use planet_config_mod,      only: gravity, Rd, cp, p_zero, kappa, &
                                  recip_epsilon
use physics_common_mod,     only: qsaturation

use log_mod,                only: log_event, log_scratch_space, &
                                  LOG_LEVEL_INFO, LOG_LEVEL_ERROR

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: vert_balance_kernel_type
  private
  type(arg_type) :: meta_args(9) = (/                                     &
       arg_type(GH_FIELD,   GH_REAL, GH_WRITE, Wtheta),                   &
       arg_type(GH_FIELD,   GH_REAL, GH_WRITE, Wtheta),                   &
       arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),                       &
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta),                   &
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta),                   &
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta),                   &
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),                       &
       arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),              &
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3) &
       /)

  type(func_type) :: meta_funcs(1) = (/                                   &
       func_type(ANY_SPACE_9, GH_BASIS)                                   &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
  integer :: gh_evaluator_targets(1) = (/ Wtheta /)
contains
  procedure, nopass :: vert_balance_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: vert_balance_code

contains

!> @brief Computes vertical balance
!> @param[in] nlayers Number of layers
!> @param[inout] theta Potential temperature field
!> @param[inout] mr_v  Water vapour mixing ratio
!> @param[inout] exner Exner pressure field
!> @param[in] height_wt Height coordinate in Wtheta
!> @param[in] height_w3 Height coordinate in W3
!> @param[in] chi_1 First component of the chi coordinate field
!> @param[in] chi_2 Second component of the chi coordinate field
!> @param[in] chi_3 Third component of the chi coordinate field
!> @param[in] panel_id A field giving the ID for mesh panels
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
subroutine vert_balance_code( nlayers, theta, mr_v, exner,      &
                              temp_specified, vapour_specified, &
                              height_wt, height_w3,             &
                              chi_1, chi_2, chi_3, panel_id,    &
                              ndf_wt, undf_wt, map_wt,          &
                              ndf_w3, undf_w3, map_w3,          &
                              ndf_chi, undf_chi, map_chi,       &
                              basis_chi_on_wt,                  &
                              ndf_pid, undf_pid, map_pid )

  use analytic_pressure_profiles_mod, only : analytic_pressure
  use sci_chi_transform_mod,          only : chi2xyz

  implicit none

  ! Arguments
  integer(kind=i_def),                              intent(in)    :: nlayers
  integer(kind=i_def),                              intent(in)    :: ndf_wt, undf_wt
  integer(kind=i_def),                              intent(in)    :: ndf_w3, undf_w3
  integer(kind=i_def),                              intent(in)    :: ndf_chi, undf_chi
  integer(kind=i_def),                              intent(in)    :: ndf_pid, undf_pid
  integer(kind=i_def), dimension(ndf_wt),           intent(in)    :: map_wt
  integer(kind=i_def), dimension(ndf_w3),           intent(in)    :: map_w3
  integer(kind=i_def), dimension(ndf_chi),          intent(in)    :: map_chi
  integer(kind=i_def), dimension(ndf_pid),          intent(in)    :: map_pid
  real(kind=r_def),    dimension(undf_wt),          intent(inout) :: theta
  real(kind=r_def),    dimension(undf_wt),          intent(inout) :: mr_v
  real(kind=r_def),    dimension(undf_w3),          intent(inout) :: exner
  real(kind=r_def),    dimension(undf_wt),          intent(in)    :: temp_specified
  real(kind=r_def),    dimension(undf_wt),          intent(in)    :: vapour_specified
  real(kind=r_def),    dimension(undf_wt),          intent(in)    :: height_wt
  real(kind=r_def),    dimension(undf_w3),          intent(in)    :: height_w3
  real(kind=r_def),    dimension(undf_chi),         intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def),    dimension(undf_pid),         intent(in)    :: panel_id
  real(kind=r_def),    dimension(1,ndf_chi,ndf_wt), intent(in)    :: basis_chi_on_wt

  ! Internal variables
  integer(kind=i_def)                  :: k, dfc, wt_dof, ipanel
  real(kind=r_def)                     :: dz(0:nlayers)
  real(kind=r_def)                     :: g_local(0:nlayers-1)
  real(kind=r_def)                     :: wt(0:nlayers)
  real(kind=r_def)                     :: coords(3), xyz(3)
  real(kind=r_def)                     :: exner_itn(0:nlayers)
  real(kind=r_def), dimension(ndf_chi) :: chi_1_e, chi_2_e, chi_3_e

  real(kind=r_def),    parameter :: eps = 1.0e-5_r_def
  integer(kind=i_def), parameter :: nitns = 100_i_def
  real(kind=r_def)    :: tol
  real(kind=r_def)    :: balance_plus, balance_minus, balance_zero
  real(kind=r_def)    :: delta_exner, dbalance_dexner
  real(kind=r_def)    :: exner_top, p_top, t_abs_top
  integer(kind=i_def) :: itn

  tol = 8.0_r_def * spacing(1.0_r_def)

  ipanel = int(panel_id(map_pid(1)), i_def)

  wt_dof = 1

  do dfc = 1, ndf_chi
    chi_1_e(dfc) = chi_1( map_chi(dfc) )
    chi_2_e(dfc) = chi_2( map_chi(dfc) )
    chi_3_e(dfc) = chi_3( map_chi(dfc) )
  end do

  ! Horizontal coordinates of cell bottom or top
  coords(:) = 0.0_r_def
  do dfc = 1, ndf_chi
    coords(1) = coords(1) + chi_1_e(dfc)*basis_chi_on_wt(1,dfc,wt_dof)
    coords(2) = coords(2) + chi_2_e(dfc)*basis_chi_on_wt(1,dfc,wt_dof)
    coords(3) = coords(3) + chi_3_e(dfc)*basis_chi_on_wt(1,dfc,wt_dof)
  end do

  call chi2xyz(coords(1), coords(2), coords(3), ipanel, xyz(1), xyz(2), xyz(3))

  ! Exner at the model surface
  exner_itn(0) = analytic_pressure( xyz, test, 0.0_r_def)

  ! Set local value of gravity at each vertical level
  if (shallow) then
    do k = 0, nlayers-1
      g_local(k) = gravity
    end do
  else
    do k = 0, nlayers-1
      g_local(k) = gravity * (planet_radius / (planet_radius + height_wt(map_wt(1)+k))) ** 2
    end do
  end if

  wt(0) = 0.0_r_def
  dz(0) = height_w3(map_w3(1))-height_wt(map_wt(1))
  do k = 1, nlayers - 1
    dz(k) = height_w3(map_w3(1)+k)-height_w3(map_w3(1)+k-1)
    wt(k) = ( height_wt(map_wt(1)+k) - height_w3(map_w3(1)+k-1) ) / dz(k)
  end do
  k = nlayers - 1
  dz(k+1) = height_wt(map_wt(2)+k) - height_w3(map_w3(1)+k)
  wt(k+1) = ( height_wt(map_wt(2)+k) - height_w3(map_w3(1)+k-1) ) / dz(k)


  ! Initial guess for Exner
  select case( temp_variable )
  case( profile_variable_potential )
    theta(:) = temp_specified(:)
    do k = 0, nlayers-1
      exner_itn(k+1) = exner_itn(k) - g_local(k) * dz(k) / (cp * theta(map_wt(1)+k))
    end do
  case( profile_variable_absolute )
    exner_itn(0) = log(exner_itn(0))
    do k = 0, nlayers-1
      exner_itn(k+1) = exner_itn(k) - g_local(k) * dz(k) / &
                                      (cp * temp_specified(map_wt(1)+k))
    end do
    exner_itn(:) = exp(exner_itn(:))
  case default
    call log_event( 'Unknown profile variable in initial_temperature namelist.', &
                    LOG_LEVEL_ERROR)
  end select

  do k = 0, nlayers - 1

    do itn = 1, nitns

      balance_zero = calc_balance( exner_itn(k+1), exner_itn(k), &
                                   theta(map_wt(1)+k), mr_v(map_wt(1)+k), &
                                   temp_specified(map_wt(1)+k), &
                                   vapour_specified(map_wt(1)+k), &
                                   dz(k), wt(k), g_local(k) )
      balance_minus = calc_balance( (1.0_r_def-eps)*exner_itn(k+1), &
                                    exner_itn(k), &
                                    theta(map_wt(1)+k), mr_v(map_wt(1)+k), &
                                    temp_specified(map_wt(1)+k), &
                                    vapour_specified(map_wt(1)+k), &
                                    dz(k), wt(k), g_local(k) )
      balance_plus = calc_balance( (1.0_r_def+eps)*exner_itn(k+1), &
                                   exner_itn(k), &
                                   theta(map_wt(1)+k), mr_v(map_wt(1)+k), &
                                   temp_specified(map_wt(1)+k), &
                                   vapour_specified(map_wt(1)+k), &
                                   dz(k), wt(k), g_local(k) )

      dbalance_dexner = ( balance_plus - balance_minus ) / &
                        ( 2.0_r_def * exner_itn(k+1) * eps )

      delta_exner = -balance_zero / dbalance_dexner

      exner_itn(k+1) = exner_itn(k+1) + delta_exner

      if ( abs( delta_exner ) <= tol ) goto 10
    end do

    write(log_scratch_space,'(''Itn fail in vetr_balance at k = '', i0)') k
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)

10  continue

    exner(map_w3(1)+k) = exner_itn(k+1)

  end do

  ! Ensure that theta and mr_v are consistent with final iterate of exner
  exner_top = exp( (1.0_r_def - wt(nlayers)) * log(exner_itn(nlayers-1)) + &
                   wt(nlayers) * log(exner_itn(nlayers)) )

  select case( temp_variable )
    case( profile_variable_absolute )
      t_abs_top = temp_specified(map_wt(2)+nlayers-1)
      theta(map_wt(2)+nlayers-1) =  t_abs_top / exner_top
    case( profile_variable_potential )
      t_abs_top = exner_top * theta(map_wt(2)+nlayers-1)
  end select

  if ( vapour_variable == profile_variable_rh ) then
    p_top = p_zero * exner_top ** (1.0_r_def / kappa)
    mr_v(map_wt(2)+nlayers-1) = vapour_specified(map_wt(2)+nlayers-1) * &
                                qsaturation( t_abs_top, 0.01_r_def * p_top )
  end if

  do k = 0, nlayers - 1
    balance_zero = calc_balance( exner_itn(k+1), exner_itn(k),          &
                                 theta(map_wt(1)+k), mr_v(map_wt(1)+k), &
                                 temp_specified(map_wt(1)+k),           &
                                 vapour_specified(map_wt(1)+k),         &
                                 dz(k), wt(k), g_local(k) )
  end do

end subroutine vert_balance_code

function calc_balance( exner_above, exner_below, theta, mr_v, &
                       temp_specified, vapour_specified, dz, wt, g ) result( balance )

implicit none

real( kind=r_def ), intent(in)    :: exner_above
real( kind=r_def ), intent(in)    :: exner_below
real( kind=r_def ), intent(inout) :: theta
real( kind=r_def ), intent(inout) :: mr_v
real( kind=r_def ), intent(in)    :: temp_specified
real( kind=r_def ), intent(in)    :: vapour_specified
real( kind=r_def ), intent(in)    :: dz, wt, g
real( kind=r_def )                :: balance

real( kind=r_def ) :: theta_virtual, exner_theta, temp, p

exner_theta = (1.0_r_def - wt) * exner_below + wt * exner_above

select case( temp_variable )
  case( profile_variable_potential )
    theta = temp_specified
    temp = theta * exner_theta
  case( profile_variable_absolute )
    temp = temp_specified
    theta = temp / exner_theta
end select

select case( vapour_variable )
  case( profile_variable_mr )
    mr_v = vapour_specified
  case( profile_variable_rh )
    p = p_zero * exner_theta ** (1.0_r_def / kappa)
    mr_v = vapour_specified * qsaturation(temp, 0.01_r_def*p)
  case default
    call log_event( 'Unknown profile_variable in initial_vapour namelist.', &
                    LOG_LEVEL_ERROR )
end select

theta_virtual = ( 1.0_r_def + recip_epsilon * mr_v ) * theta / (1.0_r_def + mr_v )

balance = g + cp * theta_virtual * ( exner_above - exner_below ) / dz

end function calc_balance

end module vert_balance_kernel_mod
