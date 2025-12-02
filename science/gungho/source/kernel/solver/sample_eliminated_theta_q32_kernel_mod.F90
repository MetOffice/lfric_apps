!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Compute the q32 matrix for analytic elimination of theta. The family
!!        of qXY matrices come from elimination of theta from the mixed solver.
!> @details Operator to map the velocity into the left hand side of the
!!          equation of state: q32 = q32 +  const*<sigma, (dtheta/dchi3)/(theta*detJ) * k.v>
!!          where v is a basis function in the W2 space,
!!          sigma is a test function in the W3 space and k is unit vector
!!          in the vertical direction of the reference cell.
!!          For more details, see the solver section of
!!          https://code.metoffice.gov.uk/trac/lfric/wiki/GhaspSupport/Documentation
module sample_eliminated_theta_q32_kernel_mod

  use argument_mod,      only: arg_type, func_type,       &
                               GH_OPERATOR, GH_FIELD,     &
                               GH_REAL, GH_SCALAR,        &
                               GH_READ, GH_WRITE,         &
                               GH_BASIS, GH_DIFF_BASIS,   &
                               CELL_COLUMN, GH_EVALUATOR
  use constants_mod,     only: i_def, r_def, r_solver
  use fs_continuity_mod, only: W3, W2, Wtheta
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  type, public, extends(kernel_type) :: sample_eliminated_theta_q32_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/                                      &
        arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W3, W2),                    &
        arg_type(GH_FIELD,    GH_REAL, GH_READ,  Wtheta),                    &
        arg_type(GH_FIELD, GH_REAL, GH_READ, W3),                            &
        arg_type(GH_SCALAR,   GH_REAL, GH_READ)                              &
        /)
    type(func_type) :: meta_funcs(3) = (/                 &
        func_type(W3,     GH_BASIS),                      &
        func_type(W2,     GH_BASIS),                      &
        func_type(Wtheta, GH_BASIS, GH_DIFF_BASIS)        &
        /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: sample_eliminated_theta_q32_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public sample_eliminated_theta_q32_code

contains

!> @brief Compute the q32 matrix that arises from analytic elimination of theta
!!        in the equation of state:
!!        q32 = q32 +  const*<sigma, (dtheta/dchi3)/(theta*detJ) * k.v>.
!> @param[in]     cell           Horizontal cell index.
!> @param[in]     nlayers        Number of layers
!> @param[in]     ncell_3d       Number of cells in the 3D mesh
!> @param[in,out] q32_op         sampection matrix
!> @param[in]     theta          Potential temperature field
!> @param[in]     detj_at_w3     Det J evaluated at W3 DoFs
!> @param[in]     const          Constant scalar to multiply operator by
!> @param[in]     ndf_w3         Degrees of freedom per cell for the pressure space
!> @param[in]     undf_w3        Total degrees of freedom for the pressure space
!> @param[in]     map_w3         Cell dofmap for the pressure space
!> @param[in]     basis_w3       Basis function for the pressure space
!!                               evaluated at W3 DoFs
!> @param[in]     ndf_w2         Degrees of freedom per cell for the velocity space
!> @param[in]     basis_w2       Vector basis function for the velocity space
!!                               evaluated at W3 DoFs
!> @param[in]     ndf_wt         Degrees of freedom per cell for the theta space
!> @param[in]     undf_wt        Total degrees of freedom for the theta space
!> @param[in]     map_wt         Cell dofmap for the theta space
!> @param[in]     basis_wt       Basis function for the theta space
!!                               evaluated at W3 DoFs
!> @param[in]     diff_basis_wt  Differential basis function for the theta space
!!                               evaluated at W3 DoFs
subroutine sample_eliminated_theta_q32_code(cell, nlayers, ncell_3d,   &
                                            q32_op,                    &
                                            theta,                     &
                                            detj_at_w3,                &
                                            const,                     &
                                            ndf_w3, undf_w3,           &
                                            map_w3, basis_w3,          &
                                            ndf_w2, basis_w2,          &
                                            ndf_wt, undf_wt, map_wt,   &
                                            basis_wt, diff_basis_wt)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ncell_3d, cell
  integer(kind=i_def), intent(in) :: ndf_w2, ndf_w3, undf_w3
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt

  integer(kind=i_def), dimension(ndf_wt), intent(in)  :: map_wt
  integer(kind=i_def), dimension(ndf_w3), intent(in)  :: map_w3

  real(kind=r_def), dimension(3, ndf_w2, ndf_w3), intent(in) :: basis_w2
  real(kind=r_def), dimension(1, ndf_w3, ndf_w3), intent(in) :: basis_w3
  real(kind=r_def), dimension(1, ndf_wt, ndf_w3), intent(in) :: basis_wt
  real(kind=r_def), dimension(3, ndf_wt, ndf_w3), intent(in) :: diff_basis_wt

  real(kind=r_solver), dimension(ncell_3d, ndf_w3, ndf_w2),  intent(inout) :: q32_op
  real(kind=r_solver), dimension(undf_wt), intent(in)                      :: theta
  real(kind=r_solver), dimension(undf_w3), intent(in)                      :: detj_at_w3
  real(kind=r_solver),                     intent(in)                      :: const

  ! Internal variables
  integer(kind=i_def) :: df, df3, df2, dft, k, ik

  real(kind=r_solver)    :: dthetadz_q, theta_q, detj_e, prod

  real(kind=r_solver), dimension(ndf_w3, ndf_w2) :: samp

  real(kind=r_solver), dimension(3, ndf_w2, ndf_w3) :: rsol_basis_w2
  real(kind=r_solver), dimension(1, ndf_w3, ndf_w3) :: rsol_basis_w3
  real(kind=r_solver), dimension(1, ndf_wt, ndf_w3) :: rsol_basis_wt
  real(kind=r_solver), dimension(3, ndf_wt, ndf_w3) :: rsol_diff_basis_wt

  rsol_basis_w2      = real(basis_w2, r_solver)
  rsol_basis_w3      = real(basis_w3, r_solver)
  rsol_basis_wt      = real(basis_wt, r_solver)
  rsol_diff_basis_wt = real(diff_basis_wt, r_solver)


  do k = 0, nlayers-1
    ik = 1 + k + (cell-1)*nlayers

    samp(:,:) = 0.0_r_solver
    do df = 1, ndf_w3
      theta_q = 0.0_r_solver
      dthetadz_q = 0.0_r_solver
      detj_e = 0.0_r_solver
      do dft = 1, ndf_wt
        dthetadz_q = dthetadz_q + theta(map_wt(dft)+k)*rsol_diff_basis_wt(3, dft, df)
        theta_q    = theta_q    + theta(map_wt(dft)+k)*rsol_basis_wt(1, dft, df)
      end do
      do df3 = 1, ndf_w3
        detj_e    = detj_e   + detj_at_w3(map_w3(df3)+k)*rsol_basis_w3(1, df3, df)
      end do
      ! Ensure that dtheta/dz (and hence the static stability, N^2 =
      ! g/theta*dtheta/dz) is positive
      dthetadz_q = max(1.0_r_solver, dthetadz_q)
      do df2 = 1, ndf_w2
        prod = const* dthetadz_q/theta_q * rsol_basis_w2(3,df2,df)
        do df3 = 1, ndf_w3
          samp(df3,df2) = samp(df3,df2) + rsol_basis_w3(1,df3,df)* prod/detj_e
        end do
      end do
    end do
    ! Incrementally add sampled operator
    q32_op(ik,:,:) = q32_op(ik,:,:) + samp(:,:)
  end do

end subroutine sample_eliminated_theta_q32_code

end module sample_eliminated_theta_q32_kernel_mod
