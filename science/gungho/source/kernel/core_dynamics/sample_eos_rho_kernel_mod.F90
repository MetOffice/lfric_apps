!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Computes density from equation of state

module sample_eos_rho_kernel_mod

use argument_mod,               only : arg_type, func_type,   &
                                       GH_FIELD, GH_REAL,     &
                                       GH_SCALAR,             &
                                       GH_READ, GH_WRITE,     &
                                       ANY_SPACE_1, GH_BASIS, &
                                       CELL_COLUMN, GH_EVALUATOR
use constants_mod,              only : r_def, i_def
use idealised_config_mod,       only : test
use fs_continuity_mod,          only : Wtheta, W3
use kernel_mod,                 only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: sample_eos_rho_kernel_type
  private
  type(arg_type) :: meta_args(7) = (/                       &
       arg_type(GH_FIELD,  GH_REAL, GH_WRITE, W3),          &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  W3),          &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_SPACE_1), &
       arg_type(GH_FIELD,  GH_REAL, GH_READ,  ANY_SPACE_1), &
       arg_type(GH_SCALAR, GH_REAL, GH_READ),               &
       arg_type(GH_SCALAR, GH_REAL, GH_READ),               &
       arg_type(GH_SCALAR, GH_REAL, GH_READ)                &
       /)
  type(func_type) :: meta_funcs(2) = (/                    &
       func_type(W3,          GH_BASIS),                   &
       func_type(ANY_SPACE_1, GH_BASIS)                    &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: sample_eos_rho_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: sample_eos_rho_code
contains

!> @brief Computes density from equation of state
!! @param[in] nlayers Number of layers
!! @param[in,out] rho Density field
!! @param[in] exner Exner pressure field
!! @param[in] theta Potential temperature field
!! @param[in] moist_dyn_gas Moist dynamics factor
!! @param[in] kappa Ratio of rd and cp
!! @param[in] rd Specific heat of dry air at constant density
!! @param[in] p_zero Reference surface pressure
!! @param[in] ndf_w3 Number of degrees of freedom per cell for W3
!! @param[in] undf_w3 Number of unique degrees of freedom for W3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for W3
!! @param[in] basis_3 Basis functions evaluated at degrees of freedom for W3
!! @param[in] ndf_wt Number of degrees of freedom per cell for Wtheta
!! @param[in] undf_wt Number of unique degrees of freedom for Wtheta
!! @param[in] map_wt Dofmap for the cell at the base of the column for Wtheta
!! @param[in] basis_t Basis functions evaluated at degrees of freedom for W3
subroutine sample_eos_rho_code(nlayers, rho, exner,              &
                               theta, moist_dyn_gas,             &
                               kappa, rd, p_zero,                &
                               ndf_w3, undf_w3, map_w3, basis_3, &
                               ndf_wt, undf_wt, map_wt, basis_t)

  use analytic_temperature_profiles_mod, only : analytic_temperature

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf_w3, undf_w3, ndf_wt, undf_wt
  integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt

  real(kind=r_def), dimension(undf_w3),  intent(inout)       :: rho
  real(kind=r_def), dimension(undf_w3),  intent(in)          :: exner
  real(kind=r_def), dimension(undf_wt),  intent(in)          :: theta
  real(kind=r_def), dimension(undf_wt),  intent(in)          :: moist_dyn_gas
  real(kind=r_def), dimension(1,ndf_w3,ndf_w3),  intent(in)  :: basis_3
  real(kind=r_def), dimension(1,ndf_wt,ndf_w3),  intent(in)  :: basis_t

  real(kind=r_def), intent(in) :: kappa
  real(kind=r_def), intent(in) :: rd
  real(kind=r_def), intent(in) :: p_zero

  ! Internal variables
  integer(kind=i_def)                  :: k, df, dft, df3
  real(kind=r_def), dimension(ndf_w3)  :: exner_e
  real(kind=r_def), dimension(ndf_wt)  :: theta_vd_e
  real(kind=r_def)                     :: exner_cell, theta_vd_cell

  ! Compute density from eqn of state
  do k = 0, nlayers-1

    do df3 = 1, ndf_w3
      exner_e(df3) = exner( map_w3(df3) + k)
    end do

    do dft = 1, ndf_wt
      theta_vd_e(dft) = theta( map_wt(dft) + k) * moist_dyn_gas( map_wt(dft) + k)
    end do

    do df = 1, ndf_w3

      exner_cell = 0.0_r_def
      do df3 = 1, ndf_w3
        exner_cell = exner_cell + exner_e(df3)*basis_3(1,df3,df)
      end do

      theta_vd_cell = 0.0_r_def
      do dft = 1, ndf_wt
        theta_vd_cell = theta_vd_cell + theta_vd_e(dft)*basis_t(1,dft,df)
      end do

      rho(map_w3(df)+k) = (p_zero*exner_cell**((1.0_r_def - kappa)/kappa))/(rd*theta_vd_cell)
    end do

  end do

end subroutine sample_eos_rho_code

end module sample_eos_rho_kernel_mod
