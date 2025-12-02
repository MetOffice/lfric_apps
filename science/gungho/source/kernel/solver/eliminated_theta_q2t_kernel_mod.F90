!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Compute the q2t matrix for analytic elimination of theta by Galerkin Projection.
!!        The family of qXY matrices come from elimination of theta from the mixed solver.
!> @details Operator to map the residual from the thermodynamic equation
!!          to the momentum equation: q2t = - const * norm_u * <v, k dexner/dchi3 * w>
!!          where v is a test function in the W2 space,
!!          w is a basis function in the Wtheta space and k is unit vector
!!          in the vertical direction of the reference cell.
!!          For more details, see the solver section of
!!          https://code.metoffice.gov.uk/trac/lfric/wiki/GhaspSupport/Documentation
module eliminated_theta_q2t_kernel_mod

  use argument_mod,      only: arg_type, func_type,     &
                               GH_OPERATOR, GH_FIELD,   &
                               GH_REAL, GH_SCALAR,      &
                               GH_READ, GH_WRITE,       &
                               GH_BASIS, GH_DIFF_BASIS, &
                               CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,     only: i_def, r_def, r_solver
  use fs_continuity_mod, only: W2, Wtheta
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------

  type, public, extends(kernel_type) :: eliminated_theta_q2t_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                       &
        arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W2, Wtheta), &
        arg_type(GH_FIELD,    GH_REAL, GH_READ,  W2),         &
        arg_type(GH_FIELD,    GH_REAL, GH_READ,  Wtheta),     &
        arg_type(GH_SCALAR,   GH_REAL, GH_READ)               &
        /)
    type(func_type) :: meta_funcs(2) = (/                     &
        func_type(W2,     GH_BASIS),                          &
        func_type(Wtheta, GH_BASIS, GH_DIFF_BASIS)            &
        /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: eliminated_theta_q2t_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public eliminated_theta_q2t_code

contains

!> @brief Compute the q2t matrix that arises from analytic elimination of theta:
!!        q2t = - const * norm_u * <v, k dexner/dchi3 * w>.
!> @param[in]     cell          Horizontal cell index.
!> @param[in]     nlayers       Number of layers.
!> @param[in]     ncell_3d      Number of cells in the 3D mesh
!> @param[in,out] q2t_op        Projection matrix
!> @param[in]     norm_u        Normalisation for the momentum equation
!> @param[in]     exner         Exner pressure in Wtheta space
!> @param[in]     const         Constant scalar to multiply operator by
!> @param[in]     ndf_w2        Degrees of freedom per cell for the velocity space
!> @param[in]     undf_w2       Total degrees of freedom for the velocity space
!> @param[in]     map_w2        Cell dofmap for the velocity space
!> @param[in]     basis_w2      Vector basis function for the velocity space
!!                              evaluated on quadrature points
!> @param[in]     ndf_wt        Degrees of freedom per cell for the theta space
!> @param[in]     undf_wt       Total degrees of freedom for the theta space
!> @param[in]     map_wt        Cell dofmap for the theta space
!> @param[in]     basis_wt      Basis function for the theta space
!!                              evaluated on quadrature points
!> @param[in]     diff_basis_wt Differential basis function for the theta space
!!                              evaluated on quadrature points
!> @param[in]     nqp_h         Number of horizontal quadrature points
!> @param[in]     nqp_v         Number of vertical quadrature points
!> @param[in]     wqp_h         Horizontal quadrature weights
!> @param[in]     wqp_v         Vertical quadrature weights
subroutine eliminated_theta_q2t_code(cell, nlayers, ncell_3d, &
                                     q2t_op,                  &
                                     norm_u, exner,           &
                                     const,                   &
                                     ndf_w2, undf_w2, map_w2, &
                                     basis_w2,                &
                                     ndf_wt, undf_wt, map_wt, &
                                     basis_wt, diff_basis_wt, &
                                     nqp_h, nqp_v, wqp_h, wqp_v)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ncell_3d, cell
  integer(kind=i_def), intent(in) :: nqp_h, nqp_v
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt

  integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2

  real(kind=r_def), dimension(3, ndf_w2, nqp_h, nqp_v), intent(in) :: basis_w2
  real(kind=r_def), dimension(1, ndf_wt, nqp_h, nqp_v), intent(in) :: basis_wt
  real(kind=r_def), dimension(3, ndf_wt, nqp_h, nqp_v), intent(in) :: diff_basis_wt

  real(kind=r_solver), dimension(ncell_3d, ndf_w2, ndf_wt), intent(inout) :: q2t_op

  real(kind=r_solver), dimension(undf_wt), intent(in) :: exner
  real(kind=r_solver), dimension(undf_w2), intent(in) :: norm_u
  real(kind=r_solver),                     intent(in) :: const

  real(kind=r_def), dimension(nqp_h),   intent(in) :: wqp_h
  real(kind=r_def), dimension(nqp_v),   intent(in) :: wqp_v

  ! Internal variables
  integer(kind=i_def) :: dft, df2, k, ik
  integer(kind=i_def) :: qp1, qp2

  real(kind=r_solver) :: dexnerdz_q, integrand, wt
  real(kind=r_solver), dimension(3, ndf_w2, nqp_h, nqp_v) :: rsol_basis_w2
  real(kind=r_solver), dimension(1, ndf_wt, nqp_h, nqp_v) :: rsol_basis_wt
  real(kind=r_solver), dimension(3, ndf_wt, nqp_h, nqp_v) :: rsol_diff_basis_wt

  rsol_basis_w2 = real(basis_w2, r_solver)
  rsol_basis_wt = real(basis_wt, r_solver)
  rsol_diff_basis_wt = real(diff_basis_wt, r_solver)


  do k = 0, nlayers-1
     ik = 1 + k + (cell-1)*nlayers

    q2t_op(ik, :, :) = 0.0_r_solver
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h

        dexnerdz_q = 0.0_r_solver
        do dft = 1, ndf_wt
          dexnerdz_q = dexnerdz_q + exner(map_wt(dft)+k)*rsol_diff_basis_wt(3, dft, qp1, qp2)
        end do

        wt = real(wqp_h(qp1) * wqp_v(qp2), r_solver)
        do dft = 1, ndf_wt
          integrand = wt * const  &
                    * dexnerdz_q * rsol_basis_wt(1,dft,qp1,qp2)
          do df2 = 1, ndf_w2
            q2t_op(ik,df2,dft) = q2t_op(ik,df2,dft)    &
                               - norm_u(map_w2(df2)+k) &
                                *rsol_basis_w2(3,df2,qp1,qp2)*integrand
          end do
        end do
      end do
    end do
  end do

end subroutine eliminated_theta_q2t_code

end module eliminated_theta_q2t_kernel_mod
