!-----------------------------------------------------------------------------
! (C) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Compute the operators for the left hand side of the equation of
!>        state.
!>
!> @details Compute the operators for the semi-implicit left hand side of the
!>          equation of state. These are:
!>          m3exner = (1-kappa)/kappa*E*sigma*1/exner
!>          m3rho   = sigma*1/rho
!>          p3theta = gamma*1/theta
!>          for functions sigma in W3 and gamma in the theta space
!>
module sample_eos_operators_kernel_mod

  use argument_mod,            only: arg_type, func_type,         &
                                     GH_OPERATOR, GH_FIELD,       &
                                     GH_SCALAR, GH_REAL, GH_READ, &
                                     GH_WRITE, ANY_SPACE_1, GH_BASIS, &
                                     CELL_COLUMN, GH_EVALUATOR
  use constants_mod,           only: r_def, r_solver, i_def
  use fs_continuity_mod,       only: W3, Wtheta
  use kernel_mod,              only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  type, public, extends(kernel_type) :: sample_eos_operators_kernel_type
    private
    type(arg_type) :: meta_args(9) = (/                                       &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W3, W3),                    &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W3, W3),                    &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W3, Wtheta),                &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W3),                        &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  W3),                        &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  Wtheta),                    &
         arg_type(GH_SCALAR,   GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR,   GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR,   GH_REAL, GH_READ)                              &
         /)
    type(func_type) :: meta_funcs(2) = (/                                     &
         func_type(W3,          GH_BASIS),                                    &
         func_type(Wtheta,      GH_BASIS)                                     &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: sample_eos_operators_code
  end type sample_eos_operators_kernel_type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: sample_eos_operators_code

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
!! @param[in] exner Reference pressure
!! @param[in] rho Reference density
!! @param[in] theta Reference potential temperature
!! @param[in] kappa Ratio of rd and cp
!! @param[in] rd Specific heat of dry air at constant density
!! @param[in] p_zero Reference surface pressure
!! @param[in] ndf_w3 Number of degrees of freedom per cell for the operator space
!! @param[in] undf_w3 Total number of degrees of freedom for the W3 space
!! @param[in] map_w3 Dofmap for the bottom layer in the W3 space
!! @param[in] basis_w3 Basis functions evaluated at W3 DoFs
!! @param[in] ndf_wt Number of degrees of freedom per cell for the theta space
!! @param[in] undf_wt Total number of degrees of freedom for the theta space
!! @param[in] map_wt Dofmap for the bottom layer in the theta space
!! @param[in] basis_wt Basis functions evaluated W3 DoFs
subroutine sample_eos_operators_code(cell, nlayers,                      &
                                     ncell_3d1, m3exner,                 &
                                     ncell_3d2, m3rho,                   &
                                     ncell_3d3, p3theta,                 &
                                     exner, rho, theta,                  &
                                     kappa, rd, p_zero,                  &
                                     ndf_w3, undf_w3, map_w3, basis_w3,  &
                                     ndf_wt, undf_wt, map_wt, basis_wt)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in)     :: cell
  integer(kind=i_def), intent(in)     :: nlayers
  integer(kind=i_def), intent(in)     :: ndf_w3, ndf_wt
  integer(kind=i_def), intent(in)     :: undf_w3, undf_wt
  integer(kind=i_def), intent(in)     :: ncell_3d1,  ncell_3d2, ncell_3d3

  integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_wt),  intent(in) :: map_wt

  real(kind=r_solver), dimension(ncell_3d1,ndf_w3,ndf_w3),  intent(inout)  :: m3exner
  real(kind=r_solver), dimension(ncell_3d2,ndf_w3,ndf_w3),  intent(inout)  :: m3rho
  real(kind=r_solver), dimension(ncell_3d3,ndf_w3,ndf_wt),  intent(inout)  :: p3theta

  real(kind=r_def), dimension(1,ndf_w3,ndf_w3),  intent(in) :: basis_w3
  real(kind=r_def), dimension(1,ndf_wt,ndf_w3),  intent(in) :: basis_wt

  real(kind=r_solver), dimension(undf_w3),  intent(in)           :: rho, exner
  real(kind=r_solver), dimension(undf_wt),  intent(in)           :: theta
  real(kind=r_def),                         intent(in)           :: kappa, rd, p_zero

  ! Internal variables
  integer(kind=i_def) :: df, df3, dft, k, ik
  real(kind=r_solver) :: rho_cell, exner_cell, theta_cell
  real(kind=r_solver) :: p0_over_rd, onemk_over_k

  real(kind=r_solver), dimension(1,ndf_w3,ndf_w3) :: rsol_basis_w3
  real(kind=r_solver), dimension(1,ndf_wt,ndf_w3) :: rsol_basis_wt

  rsol_basis_w3 = real(basis_w3, r_solver)
  rsol_basis_wt = real(basis_wt, r_solver)

  p0_over_rd = real(p_zero/Rd, r_solver)
  onemk_over_k = real((1.0_r_def - kappa)/kappa, r_solver)

  do k = 0, nlayers-1
    ik = 1 + k + (cell-1)*nlayers
    m3exner(ik,:,:) = 0.0_r_solver
    m3rho(ik,:,:)   = 0.0_r_solver
    p3theta(ik,:,:) = 0.0_r_solver
    do df = 1, ndf_w3
      exner_cell = 0.0_r_solver
      rho_cell = 0.0_r_solver
      do df3 = 1, ndf_w3
        exner_cell = exner_cell + exner(map_w3(df3)+k)*rsol_basis_w3(1,df3,df)
        rho_cell   = rho_cell   + rho(map_w3(df3)+k)  *rsol_basis_w3(1,df3,df)
      end do

      theta_cell = 0.0_r_solver
      do dft = 1,ndf_wt
        theta_cell = theta_cell + theta(map_wt(dft)+k)*rsol_basis_wt(1,dft,df)
      end do

      do df3 = 1, ndf_w3
        m3exner(ik,df,df3) = m3exner(ik,df,df3)              &
                             + onemk_over_k*          &
                             (p0_over_rd * exner_cell**onemk_over_k &
                              /(rho_cell*theta_cell))*rsol_basis_w3(1,df3,df)/exner_cell
        m3rho(ik,df,df3) = m3rho(ik,df,df3)                          &
                           + rsol_basis_w3(1,df3,df)/rho_cell
      end do

      do dft = 1, ndf_wt
        p3theta(ik,df,dft) = p3theta(ik,df,dft)                &
                             + rsol_basis_wt(1,dft,df)/theta_cell
      end do

    end do
  end do

end subroutine sample_eos_operators_code

end module sample_eos_operators_kernel_mod
