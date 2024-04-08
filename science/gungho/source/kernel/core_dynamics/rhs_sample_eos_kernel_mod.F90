!-----------------------------------------------------------------------------
! Copyright (c) 2021,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes rhs of the equation of state for the nonlinear equations via
!>          sampling
!>
!> The kernel computes the rhs of the equation of state for the nonlinear
!> equations,
!> That is: \f[ rhs_{\Pi} = 1 - p0/Rd * exner ^ (1-kappa)/kappa /(rho*theta_vd) \f]
!>
module rhs_sample_eos_kernel_mod

  use argument_mod,               only : arg_type, func_type,   &
                                         GH_FIELD, GH_REAL,     &
                                         GH_READ, GH_WRITE,     &
                                         GH_BASIS, GH_SCALAR,   &
                                         CELL_COLUMN, GH_EVALUATOR
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W3, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: rhs_sample_eos_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                                     &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),                       &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),                       &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  W3),                       &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta),                   &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta),                   &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                              &
         /)
    type(func_type) :: meta_funcs(2) = (/                                   &
         func_type(W3,          GH_BASIS),                                  &
         func_type(Wtheta,      GH_BASIS)                                   &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: rhs_sample_eos_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: rhs_sample_eos_code

contains

!> @brief Computes rhs of the equation of state for the nonlinear equations
!! @param[in] nlayers Number of layers
!! @param[in,out] rhs_eos RHS array for the equation of state
!! @param[in] exner Pressure
!! @param[in] rho Density
!! @param[in] theta Potential temperature
!! @param[in] moist_dyn_gas Moist dynamics factor in gas law
!! @param[in] kappa         Ratio of rd and cp
!! @param[in] rd            Specific heat of dry air at constant density
!! @param[in] p_zero        Reference surface pressure
!! @param[in] ndf_w3 Number of degrees of freedom per cell for W3
!! @param[in] undf_w3 Number of (local) unique degrees of freedom
!! @param[in] map_w3 Dofmap for the cell at the base of the column for W3
!! @param[in] w3_basis Basis functions evaluated at the W3 DoFs
!! @param[in] ndf_wt Number of degrees of freedom per cell for wt
!! @param[in] undf_wt Number of (local) unique degrees of freedom
!! @param[in] map_wt Dofmap for the cell at the base of the column for wt
!! @param[in] wt_basis Basis functions evaluated at the W3 DoFs
subroutine rhs_sample_eos_code(nlayers,                                       &
                               rhs_eos, exner, rho, theta, moist_dyn_gas,     &
                               kappa, rd, p_zero,                             &
                               ndf_w3, undf_w3, map_w3, w3_basis,             &
                               ndf_wt, undf_wt, map_wt, wt_basis)

  implicit none
  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wt, ndf_w3
  integer(kind=i_def), intent(in) :: undf_wt, undf_w3
  integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_wt),  intent(in) :: map_wt

  real(kind=r_def), dimension(1,ndf_w3,ndf_w3),  intent(in) :: w3_basis
  real(kind=r_def), dimension(1,ndf_wt,ndf_w3),  intent(in) :: wt_basis

  real(kind=r_def), dimension(undf_w3), intent(inout) :: rhs_eos
  real(kind=r_def), dimension(undf_wt), intent(in)    :: theta
  real(kind=r_def), dimension(undf_wt), intent(in)    :: moist_dyn_gas
  real(kind=r_def), dimension(undf_w3), intent(in)    :: rho, exner

  real(kind=r_def), intent(in) :: kappa
  real(kind=r_def), intent(in) :: rd
  real(kind=r_def), intent(in) :: p_zero

  ! Internal variables
  integer(kind=i_def) :: df, df3, dft, k

  real(kind=r_def), dimension(ndf_wt)  :: theta_vd_e
  real(kind=r_def), dimension(ndf_w3)  :: rho_e, exner_e
  real(kind=r_def)                     :: rho_cell, theta_vd_cell, exner_cell
  real(kind=r_def)                     :: p0_over_rd, onemk_over_k

  p0_over_rd = p_zero/Rd
  onemk_over_k = (1.0_r_def - kappa)/kappa


  do k = 0, nlayers-1

    do dft = 1, ndf_wt
      theta_vd_e(dft) = theta(map_wt(dft) + k) * moist_dyn_gas(map_wt(dft) + k)
    end do

    do df3 = 1, ndf_w3
      exner_e(df3)  = exner(map_w3(df3) + k)
      rho_e(df3)    = rho(map_w3(df3) + k)
      rhs_eos(map_w3(df3)+k) = 0.0_r_def
    end do

    do df = 1, ndf_w3

      theta_vd_cell = 0.0_r_def
      do dft = 1, ndf_wt
        theta_vd_cell = theta_vd_cell + theta_vd_e(dft)*wt_basis(1,dft,df)
      end do

      exner_cell = 0.0_r_def
      rho_cell = 0.0_r_def
      do df3 = 1, ndf_w3
        exner_cell = exner_cell + exner_e(df3)*w3_basis(1,df3,df)
        rho_cell  = rho_cell   + rho_e(df3)*w3_basis(1,df3,df)
      end do

      rhs_eos(map_w3(df)+k) = rhs_eos(map_w3(df)+k) - (1.0_r_def - (p0_over_rd * exner_cell**onemk_over_k) &
          /(rho_cell*theta_vd_cell))

    end do
  end do

end subroutine rhs_sample_eos_code

end module rhs_sample_eos_kernel_mod
