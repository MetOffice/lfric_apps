!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Compute the normalised operators for the left hand side of the
!>        equation of state.
!>
!> @details Compute the normalised operators for the semi-implicit left hand
!>          side of the equation of state. These are:
!>          m3exner = M3^{-1}*(1-kappa)/kappa*E*<sigma,sigma/exner*det(J)>
!>          m3rho   = M3^{-1}*<sigma,sigma/rho*det(J)>
!>          p3theta = M3^{-1}*<sigma,gamma/theta*det(J)>
!>          for functions sigma in W3 and gamma in the theta space
!>
module project_eos_operators_kernel_mod

  use argument_mod,            only: arg_type, func_type,         &
                                     GH_OPERATOR, GH_FIELD,       &
                                     GH_SCALAR, GH_REAL, GH_READ, &
                                     GH_WRITE, ANY_SPACE_1,       &
                                     ANY_DISCONTINUOUS_SPACE_3,   &
                                     GH_BASIS, GH_DIFF_BASIS,     &
                                     CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,           only: r_solver, r_def, i_def
  use sci_coordinate_jacobian_mod, only: pointwise_coordinate_jacobian
  use fs_continuity_mod,       only: W3, Wtheta
  use kernel_mod,              only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  type, public, extends(kernel_type) :: project_eos_operators_kernel_type
    private
    type(arg_type) :: meta_args(12) = (/                                      &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W3, W3),                    &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W3, W3),                    &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W3, Wtheta),                &
         arg_type(GH_OPERATOR, GH_REAL, GH_READ, W3, W3),                     &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W3),                        &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W3),                        &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  Wtheta),                    &
         arg_type(GH_FIELD*3,  GH_REAL, GH_READ,  ANY_SPACE_1),               &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_SCALAR,   GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR,   GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR,   GH_REAL, GH_READ)                              &
         /)
    type(func_type) :: meta_funcs(3) = (/                                     &
         func_type(W3,          GH_BASIS),                                    &
         func_type(Wtheta,      GH_BASIS),                                    &
         func_type(ANY_SPACE_1, GH_BASIS, GH_DIFF_BASIS)                      &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: project_eos_operators_code
  end type project_eos_operators_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: project_eos_operators_code

contains

!> @brief Computes the equation of state operators
!! @param[in] cell Cell number
!! @param[in] nlayers Number of layers
!! @param[in] ncell_3d1 ncell*nlayers
!! @param[in,out] m3exner W3 mass matrix weighted by reference pressure
!! @param[in] ncell_3d2 ncell*nlayers
!! @param[in,out] m3rho W3 mass matrix weighted by reference density
!! @param[in] ncell_3d3 ncell*nlayers
!! @param[in,out] p3theta Projection matrix weighted by reference potential temperature
!! @param[in] ncell_3d4 ncell*nlayers
!! @param[in] m3_inv Inverse W3 mass matrix
!! @param[in] exner Reference pressure
!! @param[in] rho Reference density
!! @param[in] theta Reference potential temperature
!! @param[in] chi1 1st coordinate field in Wchi
!! @param[in] chi2 2nd coordinate field in Wchi
!! @param[in] chi3 3rd coordinate field in Wchi
!! @param[in] panel_id Field giving the ID for mesh panels
!! @param[in] kappa Ratio of rd and cp
!! @param[in] rd Specific heat of dry air at constant density
!! @param[in] p_zero Reference surface pressure
!! @param[in] ndf_w3 Number of degrees of freedom per cell for the operator space
!! @param[in] undf_w3 Total number of degrees of freedom for the W3 space
!! @param[in] map_w3 Dofmap for the bottom layer in the W3 space
!! @param[in] basis_w3 Basis functions evaluated at quadrature points
!! @param[in] ndf_wt Number of degrees of freedom per cell for the theta space
!! @param[in] undf_wt Total number of degrees of freedom for the theta space
!! @param[in] map_wt Dofmap for the bottom layer in the theta space
!! @param[in] basis_wt Basis functions evaluated at quadrature points
!! @param[in] ndf_chi Number of degrees of freedom per cell for the coordinate field
!! @param[in] undf_chi Number of unique degrees of freedom for chi field
!! @param[in] map_chi Dofmap for the cell at the base of the column
!! @param[in] basis_chi Wchi basis functions evaluated at quadrature points
!! @param[in] diff_basis_chi Wchi differential basis functions evaluated at quadrature points
!! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
!! @param[in] undf_pid Number of unique degrees of freedom for panel_id
!! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
!! @param[in] nqp_h Number of horizontal quadrature points
!! @param[in] nqp_v Number of vertical quadrature points
!! @param[in] wqp_h Horizontal quadrature weights
!! @param[in] wqp_v Vertical quadrature weights
subroutine project_eos_operators_code(cell, nlayers,                      &
                                      ncell_3d1, m3exner,                 &
                                      ncell_3d2, m3rho,                   &
                                      ncell_3d3, p3theta,                 &
                                      ncell_3d4, m3_inv,                  &
                                      exner, rho, theta,                  &
                                      chi1, chi2, chi3,                   &
                                      panel_id,                           &
                                      kappa, rd, p_zero,                  &
                                      ndf_w3, undf_w3, map_w3, basis_w3,  &
                                      ndf_wt, undf_wt, map_wt, basis_wt,  &
                                      ndf_chi, undf_chi,                  &
                                      map_chi, basis_chi, diff_basis_chi, &
                                      ndf_pid, undf_pid, map_pid,         &
                                      nqp_h, nqp_v, wqp_h, wqp_v)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)     :: cell, nqp_h, nqp_v
  integer(kind=i_def), intent(in)     :: nlayers
  integer(kind=i_def), intent(in)     :: ndf_w3, ndf_chi, ndf_wt, ndf_pid
  integer(kind=i_def), intent(in)     :: undf_chi, undf_w3, undf_wt, undf_pid
  integer(kind=i_def), intent(in)     :: ncell_3d1,  ncell_3d2, ncell_3d3, ncell_3d4

  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_pid), intent(in) :: map_pid
  integer(kind=i_def), dimension(ndf_wt),  intent(in) :: map_wt

  real(kind=r_solver), dimension(ncell_3d1,ndf_w3,ndf_w3),  intent(inout)  :: m3exner
  real(kind=r_solver), dimension(ncell_3d2,ndf_w3,ndf_w3),  intent(inout)  :: m3rho
  real(kind=r_solver), dimension(ncell_3d3,ndf_w3,ndf_wt),  intent(inout)  :: p3theta

  real(kind=r_solver), dimension(ncell_3d4,ndf_w3,ndf_w3),  intent(in)  :: m3_inv

  real(kind=r_def), dimension(1,ndf_chi,nqp_h,nqp_v), intent(in) :: basis_chi
  real(kind=r_def), dimension(3,ndf_chi,nqp_h,nqp_v), intent(in) :: diff_basis_chi
  real(kind=r_def), dimension(1,ndf_w3, nqp_h,nqp_v), intent(in) :: basis_w3
  real(kind=r_def), dimension(1,ndf_wt, nqp_h,nqp_v), intent(in) :: basis_wt

  real(kind=r_solver), dimension(undf_w3),  intent(in) :: rho, exner
  real(kind=r_solver), dimension(undf_wt),  intent(in) :: theta
  real(kind=r_def),    dimension(undf_chi), intent(in) :: chi1
  real(kind=r_def),    dimension(undf_chi), intent(in) :: chi2
  real(kind=r_def),    dimension(undf_chi), intent(in) :: chi3
  real(kind=r_def),    dimension(undf_pid), intent(in) :: panel_id
  real(kind=r_def),                         intent(in) :: kappa, rd, p_zero

  real(kind=r_def), dimension(nqp_h), intent(in) :: wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in) :: wqp_v

  ! Internal variables
  integer(kind=i_def)                          :: df, df1, df2, k, ik
  integer(kind=i_def)                          :: qp1, qp2
  real(kind=r_solver), dimension(ndf_chi)      :: chi1_e, chi2_e, chi3_e
  real(kind=r_solver)                          :: rho_quad, exner_quad, theta_quad
  real(kind=r_solver)                          :: integrand1, integrand2, integrand3
  real(kind=r_solver), dimension(3,3)          :: jac
  real(kind=r_solver)                          :: dj
  real(kind=r_solver)                          :: p0_over_rd, onemk_over_k, wt

  real(kind=r_solver), dimension(1,ndf_chi,nqp_h,nqp_v) :: rsol_basis_chi
  real(kind=r_solver), dimension(3,ndf_chi,nqp_h,nqp_v) :: rsol_diff_basis_chi
  real(kind=r_solver), dimension(1,ndf_w3, nqp_h,nqp_v)  :: rsol_basis_w3
  real(kind=r_solver), dimension(1,ndf_wt, nqp_h,nqp_v)  :: rsol_basis_wt

  integer(kind=i_def) :: ipanel

  rsol_basis_chi      = real(basis_chi, r_solver)
  rsol_diff_basis_chi = real(diff_basis_chi, r_solver)
  rsol_basis_w3       = real(basis_w3, r_solver)
  rsol_basis_wt       = real(basis_wt, r_solver)

  p0_over_rd = real(p_zero/Rd, r_solver)
  onemk_over_k = real((1.0_r_def - kappa)/kappa, r_solver)

  ipanel = int(panel_id(map_pid(1)), i_def)

  do k = 0, nlayers-1
    do df = 1, ndf_chi
      chi1_e(df) = real(chi1(map_chi(df) + k), r_solver)
      chi2_e(df) = real(chi2(map_chi(df) + k), r_solver)
      chi3_e(df) = real(chi3(map_chi(df) + k), r_solver)
    end do
    ik = 1 + k + (cell-1)*nlayers
    m3exner(ik,:,:) = 0.0_r_solver
    m3rho(ik,:,:)   = 0.0_r_solver
    p3theta(ik,:,:) = 0.0_r_solver
    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        call pointwise_coordinate_jacobian(ndf_chi, chi1_e, chi2_e, chi3_e,     &
                                           ipanel, rsol_basis_chi(:,:,qp1,qp2), &
                                           rsol_diff_basis_chi(:,:,qp1,qp2),    &
                                           jac, dj                          )

        exner_quad = 0.0_r_solver
        rho_quad = 0.0_r_solver
        do df = 1,ndf_w3
          exner_quad = exner_quad + exner(map_w3(df)+k)*rsol_basis_w3(1,df,qp1,qp2)
          rho_quad   = rho_quad   + rho(map_w3(df)+k)  *rsol_basis_w3(1,df,qp1,qp2)
        end do
        theta_quad = 0.0_r_solver
        do df = 1,ndf_wt
          theta_quad = theta_quad + theta(map_wt(df)+k)*rsol_basis_wt(1,df,qp1,qp2)
        end do
        wt = real(wqp_h(qp1)*wqp_v(qp2), r_solver)*dj
        integrand1 = wt*onemk_over_k*(p0_over_rd*exner_quad**onemk_over_k &
            /(rho_quad*theta_quad))/exner_quad
        integrand2 = wt/rho_quad
        integrand3 = wt/theta_quad
        do df2 = 1, ndf_w3
          do df1 = 1, ndf_w3
            m3exner(ik,df1,df2) = m3exner(ik,df1,df2)              &
                              + integrand1*rsol_basis_w3(1,df1,qp1,qp2) &
                                          *rsol_basis_w3(1,df2,qp1,qp2)
            m3rho(ik,df1,df2) = m3rho(ik,df1,df2)                  &
                              + integrand2*rsol_basis_w3(1,df1,qp1,qp2) &
                                          *rsol_basis_w3(1,df2,qp1,qp2)
          end do
        end do
        do df2 = 1, ndf_wt
          do df1 = 1, ndf_w3
            p3theta(ik,df1,df2) = p3theta(ik,df1,df2)                &
                                + integrand3*rsol_basis_w3(1,df1,qp1,qp2) &
                                            *rsol_basis_wt(1,df2,qp1,qp2)
          end do
        end do
      end do
    end do

    ! Normalise by inverse W3 mass matrix
    p3theta(ik,:,:) = matmul(m3_inv(ik,:,:), p3theta(ik,:,:))
    m3exner(ik,:,:) = matmul(m3_inv(ik,:,:), m3exner(ik,:,:))
    m3rho(ik,:,:)   = matmul(m3_inv(ik,:,:), m3rho(ik,:,:))
  end do

end subroutine project_eos_operators_code

end module project_eos_operators_kernel_mod
