!-----------------------------------------------------------------------------
! Copyright (c) 2021,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Computes RHS dry density by galerkin projection from equation of state

module project_eos_rho_kernel_mod

use argument_mod,               only : arg_type, func_type,                    &
                                       GH_FIELD, GH_READ, GH_WRITE, GH_REAL,   &
                                       ANY_SPACE_2, ANY_DISCONTINUOUS_SPACE_3, &
                                       GH_BASIS, GH_DIFF_BASIS, GH_SCALAR,     &
                                       CELL_COLUMN, GH_QUADRATURE_XYoZ
use constants_mod,              only : r_def, i_def
use idealised_config_mod,       only : test
use fs_continuity_mod,          only : WTHETA, W3
use kernel_mod,                 only : kernel_type

implicit none

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: project_eos_rho_kernel_type
  private
  type(arg_type) :: meta_args(9) = (/                                   &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3),                       &
       arg_type(GH_FIELD, GH_REAL, GH_READ, W3),                        &
       arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                    &
       arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                    &
       arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_2),             &
       arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3), &
       arg_type(GH_SCALAR, GH_REAL, GH_READ),                           &
       arg_type(GH_SCALAR, GH_REAL, GH_READ),                           &
       arg_type(GH_SCALAR, GH_REAL, GH_READ)                            &
       /)
  type(func_type) :: meta_funcs(3) = (/                     &
       func_type(W3, GH_BASIS),                             &
       func_type(WTHETA, GH_BASIS),                         &
       func_type(ANY_SPACE_2, GH_BASIS, GH_DIFF_BASIS)      &
       /)
       integer :: operates_on = CELL_COLUMN
       integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass ::project_eos_rho_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public project_eos_rho_code
contains

!> @brief Computes RHS density from equation of state
!! @param[in] nlayers Number of layers
!! @param[in,out] rhs_rho RHS of density field
!! @param[in] exner Exner pressure field
!! @param[in] theta Potential temperature field
!! @param[in] moist_dyn_gas Moist dynamics factor
!! @param[in] chi1 1st coordinate field in Wchi
!! @param[in] chi2 2nd coordinate field in Wchi
!! @param[in] chi3 3rd coordinate field in Wchi
!! @param[in] panel_id Field giving the ID for mesh panels.
!! @param[in] kappa Ratio of rd and cp
!! @param[in] rd Specific heat of dry air at constant density
!! @param[in] p_zero Reference surface pressure
!! @param[in] ndf_w3 Number of degrees of freedom per cell for w3
!! @param[in] undf_w3 Number of unique degrees of freedom  for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_wt Number of degrees of freedom per cell for wtheta
!! @param[in] undf_wt Number of unique degrees of freedom  for wtheta
!! @param[in] map_wt Dofmap for the cell at the base of the column for wt
!! @param[in] wt_basis Basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_chi Number of degrees of freedom per cell for chi space
!! @param[in] undf_chi Number of unique degrees of freedom  for chi space
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi space
!! @param[in] chi_basis Wchi basis functions evaluated at gaussian quadrature points.
!! @param[in] chi_diff_basis Derivatives of Wchi basis functions
!!                           evaluated at gaussian quadrature points
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h horizontal quadrature weights
!! @param[in] wqp_v vertical quadrature weights
subroutine project_eos_rho_code(nlayers,                           &
                                rhs_rho, exner, theta,             &
                                moist_dyn_gas,                     &
                                chi1, chi2, chi3,                  &
                                panel_id,                          &
                                kappa, rd, p_zero,                 &
                                ndf_w3, undf_w3, map_w3, w3_basis, &
                                ndf_wt, undf_wt, map_wt, wt_basis, &
                                ndf_chi, undf_chi, map_chi,        &
                                chi_basis, chi_diff_basis,         &
                                ndf_pid, undf_pid, map_pid,        &
                                nqp_h, nqp_v, wqp_h, wqp_v         &
                                )

  use sci_coordinate_jacobian_mod, only: coordinate_jacobian

  implicit none

  !Arguments
  integer(kind=i_def), intent(in) :: nlayers, nqp_h, nqp_v
  integer(kind=i_def), intent(in) :: ndf_w3, undf_w3,  ndf_wt, undf_wt
  integer(kind=i_def), intent(in) :: ndf_chi, undf_chi, ndf_pid, undf_pid
  integer(kind=i_def), dimension(ndf_w3), intent(in)  :: map_w3
  integer(kind=i_def), dimension(ndf_wt), intent(in)  :: map_wt
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid

  real(kind=r_def), dimension(undf_w3),  intent(inout)       :: rhs_rho
  real(kind=r_def), dimension(undf_w3),  intent(in)          :: exner
  real(kind=r_def), dimension(undf_wt),  intent(in)          :: theta
  real(kind=r_def), dimension(undf_wt),  intent(in)          :: moist_dyn_gas
  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v),  intent(in)  :: w3_basis
  real(kind=r_def), dimension(1,ndf_wt,nqp_h,nqp_v),  intent(in)  :: wt_basis

  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_diff_basis
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v), intent(in) :: chi_basis

  real(kind=r_def), dimension(undf_chi), intent(in) :: chi1, chi2, chi3
  real(kind=r_def), dimension(undf_pid), intent(in) :: panel_id
  real(kind=r_def), dimension(nqp_h),    intent(in) :: wqp_h
  real(kind=r_def), dimension(nqp_v),    intent(in) :: wqp_v
  real(kind=r_def),                      intent(in) :: kappa
  real(kind=r_def),                      intent(in) :: rd
  real(kind=r_def),                      intent(in) :: p_zero

  !Internal variables
  integer(kind=i_def) :: k, df, dft, df3, ipanel
  integer(kind=i_def) :: qp1, qp2

  real(kind=r_def), dimension(ndf_w3)          :: exner_e
  real(kind=r_def), dimension(ndf_wt)          :: theta_vd_e
  real(kind=r_def)                             :: exner_at_quad, theta_vd_at_quad
  real(kind=r_def), dimension(ndf_w3)          :: rhs_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(ndf_chi)         :: chi1_e, chi2_e, chi3_e
  real(kind=r_def)                             :: integrand

  ipanel = int(panel_id(map_pid(1)), i_def)

  ! Compute density from eqn of state
  do k = 0, nlayers-1

    do df = 1, ndf_chi
      chi1_e(df) = chi1( map_chi(df) + k )
      chi2_e(df) = chi2( map_chi(df) + k )
      chi3_e(df) = chi3( map_chi(df) + k )
    end do

    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v,             &
                             chi1_e, chi2_e, chi3_e,            &
                             ipanel, chi_basis, chi_diff_basis, &
                             jac, dj )

    do df3 = 1, ndf_w3
      exner_e(df3) = exner( map_w3(df3) + k)
    end do

    do dft = 1, ndf_wt
      theta_vd_e(dft) = theta( map_wt(dft) + k) * moist_dyn_gas( map_wt(dft) + k)
    end do

    ! Compute RHS
    do df = 1, ndf_w3
      rhs_e(df) = 0.0_r_def
      do qp2 = 1, nqp_v
        do qp1 = 1, nqp_h

          exner_at_quad = 0.0_r_def
          do df3 = 1, ndf_w3
            exner_at_quad  = exner_at_quad + exner_e(df3)*w3_basis(1,df3,qp1,qp2)
          end do

          theta_vd_at_quad = 0.0_r_def
          do dft = 1, ndf_wt
            theta_vd_at_quad = theta_vd_at_quad + theta_vd_e(dft)*wt_basis(1,dft,qp1,qp2)
          end do

          integrand =  w3_basis(1,df,qp1,qp2) * (p_zero*exner_at_quad**((1.0_r_def - kappa)/kappa)) &
                                                 /(rd*theta_vd_at_quad) * dj(qp1,qp2)
          rhs_e(df) = rhs_e(df) + wqp_h(qp1)*wqp_v(qp2)*integrand
        end do
      end do
    end do
    do df = 1,ndf_w3
      rhs_rho(map_w3(df)+k) = rhs_e(df)
    end do

  end do
end subroutine project_eos_rho_code

end module project_eos_rho_kernel_mod
