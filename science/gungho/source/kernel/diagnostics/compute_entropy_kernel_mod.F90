!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes cell integrated dry entropy (physical and dynamic).
!>
!> @details The kernel computes the cell integrated entropy (physical and
!>          dynamic as:
!> physical: \f[ \int( \rho*\log(\theta) )dV \f]
!> dynamic: \f[ \int( 1/2*\rho*\theta*\theta )dV \f]
!>
module compute_entropy_kernel_mod

  use argument_mod,      only : arg_type, func_type,         &
                                GH_FIELD, GH_WRITE, GH_READ, &
                                GH_REAL, ANY_SPACE_9,        &
                                ANY_DISCONTINUOUS_SPACE_3,   &
                                GH_BASIS, GH_DIFF_BASIS,     &
                                CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W3, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none

  private

!---------------------------------------------------------------------------
! Public types
!---------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the
!> Psy layer.
!>
type, public, extends(kernel_type) :: compute_entropy_kernel_type
  private
  type(arg_type) :: meta_args(6) = (/                                      &
       arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),                        &
       arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),                        &
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),                        &
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta),                    &
       arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),               &
       arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3)  &
       /)
  type(func_type) :: meta_funcs(3) = (/                                   &
       func_type(W3,          GH_BASIS),                                  &
       func_type(Wtheta,      GH_BASIS),                                  &
       func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                    &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_QUADRATURE_XYoZ
contains
  procedure, nopass :: compute_entropy_code
end type

!---------------------------------------------------------------------------
! Contained functions/subroutines
!---------------------------------------------------------------------------
public :: compute_entropy_code
contains

!> @brief Compute the cell integrated entropy
!! @param[in] nlayers The number of layers
!! @param[in,out] The physical entropy, \rho\log(\theta)
!! @param[in,out] The dynamic entropy,  0.5\rho\theta^2
!! @param[in] rho The density
!! @param[in] theta Potential temperature
!! @param[in] chi_1 1st coordinate field in Wchi
!! @param[in] chi_2 2nd coordinate field in Wchi
!! @param[in] chi_3 3rd coordinate field in Wchi
!! @param[in] panel_id A field giving the ID for mesh panels
!! @param[in] ndf_w3 The number of degrees of freedom per cell for w3
!! @param[in] undf_w3 The number of unique degrees of freedom  for w3
!! @param[in] map_w3 Dofmap for the cell at the base of the column for w3
!! @param[in] w3_basis 4-dim array holding basis functions evaluated at gaussian quadrature points
!! @param[in] ndf_w2 The number of degrees of freedom per cell for w2
!! @param[in] undf_w2 The number of unique degrees of freedom  for w2
!! @param[in] map_w2 Dofmap for the cell at the base of the column for w2
!! @param[in] w2_basis 4-dim array holding basis functions evaluated at quadrature points
!! @param[in] ndf_wtheta The number of degrees of freedom per cell for wtheta
!! @param[in] undf_wtheta The number of unique degrees of freedom  for wtheta
!! @param[in] map_wtheta Dofmap for the cell at the base of the column for wtheta
!! @param[in] wtheta_basis 4-dim array holding basis functions evaluated at quadrature points
!! @param[in] ndf_chi The number of degrees of freedom per cell for chi
!! @param[in] undf_chi The number of unique degrees of freedom  for chi
!! @param[in] map_chi Dofmap for the cell at the base of the column for chi
!! @param[in] chi_basis Wchi basis functions evaluated at gaussian quadrature points
!! @param[in] chi_diff_basis Differential Wchi basis functions evaluated at gaussian quadrature point
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of quadrature points in the horizontal
!! @param[in] nqp_v Number of quadrature points in the vertical
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine compute_entropy_code( nlayers,                             &
                                 e_phys_dry, e_dyn_dry,               &
                                 rho, theta,                          &
                                 chi_1, chi_2, chi_3, panel_id,       &
                                 ndf_w3, undf_w3, map_w3, w3_basis,   &
                                 ndf_wtheta, undf_wtheta, map_wtheta, &
                                 wtheta_basis,                        &
                                 ndf_chi, undf_chi, map_chi,          &
                                 chi_basis, chi_diff_basis,           &
                                 ndf_pid, undf_pid, map_pid,          &
                                 nqp_h, nqp_v, wqp_h, wqp_v           )

  use sci_coordinate_jacobian_mod, only: coordinate_jacobian

  implicit none

  ! Arguments
  integer(kind=i_def),                        intent(in) :: nlayers, nqp_h, nqp_v
  integer(kind=i_def),                        intent(in) :: ndf_w3, ndf_wtheta
  integer(kind=i_def),                        intent(in) :: ndf_pid, ndf_chi
  integer(kind=i_def),                        intent(in) :: undf_w3, undf_wtheta
  integer(kind=i_def),                        intent(in) :: undf_pid, undf_chi
  integer(kind=i_def), dimension(ndf_w3),     intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_wtheta), intent(in) :: map_wtheta
  integer(kind=i_def), dimension(ndf_chi),    intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid),    intent(in) :: map_pid

  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v),     intent(in) :: w3_basis
  real(kind=r_def), dimension(1,ndf_wtheta,nqp_h,nqp_v), intent(in) :: wtheta_basis
  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v),    intent(in) :: chi_basis
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v),    intent(in) :: chi_diff_basis

  real(kind=r_def), dimension(undf_w3),     intent(inout) :: e_phys_dry, &
                                                             e_dyn_dry
  real(kind=r_def), dimension(undf_w3),     intent(in)    :: rho
  real(kind=r_def), dimension(undf_wtheta), intent(in)    :: theta
  real(kind=r_def), dimension(undf_chi),    intent(in)    :: chi_1, chi_2, chi_3
  real(kind=r_def), dimension(undf_pid),    intent(in)    :: panel_id

  real(kind=r_def), dimension(nqp_h),       intent(in)  ::  wqp_h
  real(kind=r_def), dimension(nqp_v),       intent(in)  ::  wqp_v

  ! Internal variables
  integer(kind=i_def)   :: df, k, loc
  integer(kind=i_def)   :: qp1, qp2
  integer(kind=i_def)   :: ipanel

  real(kind=r_def), dimension(ndf_chi)         :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac
  real(kind=r_def), dimension(ndf_w3)          :: e_phys_dry_e, e_dyn_dry_e
  real(kind=r_def), dimension(ndf_wtheta)      :: theta_e
  real(kind=r_def), dimension(ndf_w3)          :: rho_e

  real(kind=r_def) :: rho_at_quad, theta_at_quad
  real(kind=r_def) :: weight, e_phys_dry_term, e_dyn_dry_term

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
    ! Extract element arrays of chi
    do df = 1, ndf_chi
      loc = map_chi(df) + k
      chi_1_e(df) = chi_1( loc )
      chi_2_e(df) = chi_2( loc )
      chi_3_e(df) = chi_3( loc )
    end do
    call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi_1_e, chi_2_e, chi_3_e,  &
                             ipanel, chi_basis, chi_diff_basis, jac, dj)

    do df = 1, ndf_wtheta
      theta_e(df) = theta( map_wtheta(df) + k )
    end do
    do df = 1, ndf_w3
      rho_e(df)      = rho( map_w3(df) + k )
      e_phys_dry_e(df) = 0.0_r_def
      e_dyn_dry_e(df)  = 0.0_r_def
    end do
    ! Compute the entropy integrated over one cell
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        rho_at_quad   = 0.0_r_def
        do df = 1, ndf_w3
          rho_at_quad  = rho_at_quad + rho_e(df)*w3_basis(1,df,qp1,qp2)
        end do
        theta_at_quad = 0.0_r_def
        do df = 1, ndf_wtheta
          theta_at_quad = theta_at_quad                         &
                        + theta_e(df)*wtheta_basis(1,df,qp1,qp2)
        end do

        e_phys_dry_term = rho_at_quad*log(theta_at_quad)
        e_dyn_dry_term  = 0.5_r_def*rho_at_quad*theta_at_quad*theta_at_quad

        do df = 1, ndf_w3
          weight = wqp_h(qp1)*wqp_v(qp2)*rho_at_quad*dj(qp1,qp2)
          e_phys_dry_e(df) = e_phys_dry_e(df) + weight*e_phys_dry_term
          e_dyn_dry_e(df)  = e_dyn_dry_e(df)  + weight*e_dyn_dry_term
        end do
      end do
    end do
    do df = 1, ndf_w3
      e_phys_dry( map_w3(df) + k ) = e_phys_dry_e(df)
      e_dyn_dry( map_w3(df) + k )  = e_dyn_dry_e(df)
    end do
  end do

end subroutine compute_entropy_code

end module compute_entropy_kernel_mod
