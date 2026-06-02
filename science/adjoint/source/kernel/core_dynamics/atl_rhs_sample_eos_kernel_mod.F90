!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the adjoint of the equation of state by sampling.
module atl_rhs_sample_eos_kernel_mod

  use argument_mod,      only : arg_type, func_type,   &
                                GH_FIELD, GH_REAL,     &
                                GH_READ, GH_READWRITE, &
                                GH_BASIS,              &
                                GH_SCALAR,             &
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
  type, public, extends(kernel_type) :: atl_rhs_sample_eos_kernel_type
    private
    type(arg_type) :: meta_args(12) = (/                                    &
         arg_type(GH_FIELD,   GH_REAL, GH_READWRITE, W3),                   &
         arg_type(GH_FIELD,   GH_REAL, GH_READWRITE, W3),                   &
         arg_type(GH_FIELD,   GH_REAL, GH_READWRITE, W3),                   &
         arg_type(GH_FIELD,   GH_REAL, GH_READWRITE, Wtheta),               &
         arg_type(GH_FIELD,   GH_REAL, GH_READWRITE, Wtheta),               &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,      W3),                   &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,      W3),                   &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,      Wtheta),               &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,      Wtheta),               &
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
    procedure, nopass :: atl_rhs_sample_eos_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: atl_rhs_sample_eos_code

  contains

  !> @brief Computes rhs of the equation of state for the nonlinear equations
  !! @param[in]     nlayers            Number of layers
  !! @param[in,out] rhs_eos            RHS array for the tangent linear equation of state
  !! @param[in,out] exner              Change in Pressure
  !! @param[in,out] rho                Change in Density
  !! @param[in,out] theta              Change in Potential temperature
  !! @param[in,out] moist_dyn_gas      Change in Moist dynamics factor in gas law
  !! @param[in]     ls_exner           Lin state for Pressure
  !! @param[in]     ls_rho             Lin state for Density
  !! @param[in]     ls_theta           Lin state for Potential temperature
  !! @param[in]     ls_moist_dyn_gas   Lin state for Moist dynamics factor in gas law
  !! @param[in]     kappa              Ratio of rd and cp
  !! @param[in]     rd                 Specific heat of dry air at constant density
  !! @param[in]     p_zero             Reference surface pressure
  !! @param[in]     ndf_w3             Number of degrees of freedom per cell for W3
  !! @param[in]     undf_w3            Number of (local) unique degrees of freedom
  !! @param[in]     map_w3             Dofmap for the cell at the base of the column for W3
  !! @param[in]     w3_basis           W3 to W3 basis functions evaluated at quadrature points
  !! @param[in]     w3_wt_basis        W3 to Wtheta basis functions evaluated at quadrature points
  !! @param[in]     ndf_wt             Number of degrees of freedom per cell for wt
  !! @param[in]     undf_wt            Number of (local) unique degrees of freedom
  !! @param[in]     map_wt             Dofmap for the cell at the base of the column for wt
  !! @param[in]     wt_basis           Wtheta to W3 basis functions evaluated at quadrature points
  !! @param[in]     wt_wt_basis        Wtheta to Wtheta basis functions evaluated at quadrature points
  subroutine atl_rhs_sample_eos_code( nlayers,                                        &
                                      rhs_eos, exner, rho, theta, moist_dyn_gas,      &
                                      ls_exner, ls_rho, ls_theta, ls_moist_dyn_gas,   &
                                      kappa, rd, p_zero,                              &
                                      ndf_w3, undf_w3, map_w3, w3_basis, w3_wt_basis, &
                                      ndf_wt, undf_wt, map_wt, wt_basis, wt_wt_basis )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wt, ndf_w3
    integer(kind=i_def), intent(in) :: undf_wt, undf_w3
    integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3
    integer(kind=i_def), dimension(ndf_wt),  intent(in) :: map_wt

    real(kind=r_def), dimension(1,ndf_w3,ndf_w3), intent(in) :: w3_basis
    real(kind=r_def), dimension(1,ndf_w3,ndf_wt), intent(in) :: w3_wt_basis
    real(kind=r_def), dimension(1,ndf_wt,ndf_w3), intent(in) :: wt_basis
    real(kind=r_def), dimension(1,ndf_wt,ndf_wt), intent(in) :: wt_wt_basis

    real(kind=r_def), dimension(undf_w3), intent(inout) :: rhs_eos
    real(kind=r_def), dimension(undf_wt), intent(inout) :: theta
    real(kind=r_def), dimension(undf_wt), intent(inout) :: moist_dyn_gas
    real(kind=r_def), dimension(undf_w3), intent(inout) :: rho, exner
    real(kind=r_def), dimension(undf_wt), intent(in)    :: ls_theta
    real(kind=r_def), dimension(undf_wt), intent(in)    :: ls_moist_dyn_gas
    real(kind=r_def), dimension(undf_w3), intent(in)    :: ls_rho, ls_exner

    real(kind=r_def), intent(in) :: kappa
    real(kind=r_def), intent(in) :: rd
    real(kind=r_def), intent(in) :: p_zero

    ! Internal variables
    integer(kind=i_def) :: df, df3, dft, k

    real(kind=r_def), dimension(ndf_wt)  :: theta_vd_e
    real(kind=r_def), dimension(ndf_wt)  :: ls_theta_vd_e
    real(kind=r_def), dimension(ndf_w3)  :: rho_e, exner_e
    real(kind=r_def)                     :: rho_cell, theta_vd_cell, exner_cell
    real(kind=r_def), dimension(ndf_w3)  :: ls_rho_e, ls_exner_e
    real(kind=r_def)                     :: ls_rho_cell, ls_theta_vd_cell
    real(kind=r_def)                     :: ls_exner_cell
    real(kind=r_def)                     :: p0_over_rd, onemk_over_k
    real(kind=r_def)                     :: ls_eos

    p0_over_rd = p_zero / rd
    onemk_over_k = -1 + 1.0_r_def / kappa
    do k = nlayers - 1, 0, -1
      do df = 1, ndf_wt, 1
        ls_theta_vd_e(df) = ls_moist_dyn_gas(k + map_wt(df)) * ls_theta(k + map_wt(df))
      end do

      do df = 1, ndf_w3, 1
        ls_exner_e(df) = ls_exner(map_w3(df) + k)
        ls_rho_e(df) = ls_rho(map_w3(df) + k)
      end do

      exner_e = 0.0_r_def
      rho_e = 0.0_r_def
      theta_vd_e = 0.0_r_def
      do df = ndf_w3, 1, -1
        ls_theta_vd_cell = 0.0_r_def
        do dft = 1, ndf_wt, 1
          ls_theta_vd_cell = ls_theta_vd_cell + ls_theta_vd_e(dft) * wt_basis(1,dft,df)
        end do

        ls_exner_cell = 0.0_r_def
        ls_rho_cell = 0.0_r_def
        do df3 = 1, ndf_w3, 1
          ls_exner_cell = ls_exner_cell + ls_exner_e(df3) * w3_basis(1,df3,df)
          ls_rho_cell = ls_rho_cell + ls_rho_e(df3) * w3_basis(1,df3,df)
        end do

        ls_eos = -ls_exner_cell ** onemk_over_k * p0_over_rd / (ls_rho_cell * ls_theta_vd_cell)
        exner_cell = (-rhs_eos(map_w3(df) + k) * ls_eos * onemk_over_k / ls_exner_cell)
        theta_vd_cell = ls_eos * rhs_eos(map_w3(df) + k) / ls_theta_vd_cell
        rho_cell = ls_eos * rhs_eos(map_w3(df) + k) / ls_rho_cell
        rhs_eos(map_w3(df) + k) = 0.0_r_def

        do df3 = ndf_w3, 1, -1
          rho_e(df3) = rho_e(df3) + rho_cell * w3_basis(1,df3,df)
          exner_e(df3) = exner_e(df3) + exner_cell * w3_basis(1,df3,df)
        end do

        do dft = ndf_wt, 1, -1
          theta_vd_e(dft) = theta_vd_e(dft) + theta_vd_cell * wt_basis(1,dft,df)
        end do
      end do ! df

      do df = ndf_w3, 1, -1
        rho(map_w3(df) + k) = rho(map_w3(df) + k) + rho_e(df)
        exner(map_w3(df) + k) = exner(map_w3(df) + k) + exner_e(df)
      end do

      do df = ndf_wt, 1, -1
        theta(k + map_wt(df)) = theta(k + map_wt(df)) + ls_moist_dyn_gas(k + map_wt(df)) * theta_vd_e(df)
        moist_dyn_gas(k + map_wt(df)) = moist_dyn_gas(k + map_wt(df)) + ls_theta(k + map_wt(df)) * theta_vd_e(df)
      end do
    end do ! k

  end subroutine atl_rhs_sample_eos_code

end module atl_rhs_sample_eos_kernel_mod
