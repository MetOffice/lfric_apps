!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Adjoint for computing the sampling of the pressure field
!!        into the same space as density.
module atl_sample_eos_pressure_kernel_mod

  use argument_mod,      only : arg_type, func_type,       &
                                GH_FIELD,                  &
                                GH_READ, GH_READWRITE,     &
                                GH_REAL, GH_BASIS,         &
                                CELL_COLUMN, GH_EVALUATOR
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : W3, Wtheta
  use kernel_mod,        only : kernel_type
  use planet_config_mod, only : kappa, rd, p_zero

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: atl_sample_eos_pressure_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                      &
         arg_type(GH_FIELD,    GH_REAL, GH_READWRITE, W3),                   &
         arg_type(GH_FIELD,    GH_REAL, GH_READWRITE, W3),                   &
         arg_type(GH_FIELD,    GH_REAL, GH_READWRITE, Wtheta),               &
         arg_type(GH_FIELD,    GH_REAL, GH_READWRITE, Wtheta),               &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,      W3),                   &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,      Wtheta),               &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,      Wtheta)                &
         /)
    type(func_type) :: meta_funcs(2) = (/                                    &
         func_type(W3,          GH_BASIS),                                   &
         func_type(Wtheta,      GH_BASIS)                                    &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: atl_sample_eos_pressure_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: atl_sample_eos_pressure_code

  contains

  !> @brief Compute the tangent linear for sample pressure
  !! @param[in]      nlayers           Number of layers
  !! @param[in,out]  exner             Change in pressure field
  !! @param[in,out]  rho               Change in density
  !! @param[in,out]  theta             Change in potential temperature
  !! @param[in,out]  moist_dyn_gas     Change in moist dynamics factor
  !! @param[in]      ls_rho            Lin state for density
  !! @param[in]      ls_theta          Lin state for potential temperature
  !! @param[in]      ls_moist_dyn_gas  Lin state for moist dynamics factor
  !! @param[in]      ndf_w3            Number of degrees of freedom per cell for w3
  !! @param[in]      undf_w3           Number of unique degrees of freedom for w3
  !! @param[in]      map_w3            Dofmap for the cell at the base of the column for w3
  !! @param[in]      w3_basis          W3 to W3 basis functions evaluated at Gaussian quadrature points
  !! @param[in]      w3_wt_basis       W3 to Wtheta basis functions evaluated at Gaussian quadrature points
  !! @param[in]      ndf_wt            Number of degrees of freedom per cell for theta space
  !! @param[in]      undf_wt           Number of unique degrees of freedom for theta space
  !! @param[in]      map_wt            Dofmap for the cell at the base of the column for theta space
  !! @param[in]      wt_basis          Wtheta to W3 basis functions evaluated at Gaussian quadrature points
  !! @param[in]      wt_wt_basis       Wtheta to Wtheta basis functions evaluated at Gaussian quadrature points
  subroutine atl_sample_eos_pressure_code( nlayers,                                        &
                                           exner, rho, theta, moist_dyn_gas,               &
                                           ls_rho, ls_theta, ls_moist_dyn_gas,             &
                                           ndf_w3, undf_w3, map_w3, w3_basis, w3_wt_basis, &
                                           ndf_wt, undf_wt, map_wt, wt_basis, wt_wt_basis )

    use sci_coordinate_jacobian_mod, only : coordinate_jacobian

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wt, ndf_w3
    integer(kind=i_def), intent(in) :: undf_wt, undf_w3
    integer(kind=i_def), dimension(ndf_wt),  intent(in) :: map_wt
    integer(kind=i_def), dimension(ndf_w3),  intent(in) :: map_w3

    real(kind=r_def), dimension(1,ndf_w3, ndf_w3), intent(in) :: w3_basis
    real(kind=r_def), dimension(1,ndf_w3, ndf_wt), intent(in) :: w3_wt_basis
    real(kind=r_def), dimension(1,ndf_wt, ndf_w3), intent(in) :: wt_basis
    real(kind=r_def), dimension(1,ndf_wt, ndf_wt), intent(in) :: wt_wt_basis

    real(kind=r_def), dimension(undf_w3),  intent(inout) :: exner
    real(kind=r_def), dimension(undf_w3),  intent(inout) :: rho
    real(kind=r_def), dimension(undf_wt),  intent(inout) :: theta
    real(kind=r_def), dimension(undf_wt),  intent(inout) :: moist_dyn_gas
    real(kind=r_def), dimension(undf_w3),  intent(in)    :: ls_rho
    real(kind=r_def), dimension(undf_wt),  intent(in)    :: ls_theta
    real(kind=r_def), dimension(undf_wt),  intent(in)    :: ls_moist_dyn_gas

    ! Internal variables
    integer(kind=i_def) :: df, df3, dft, k

    real(kind=r_def), dimension(ndf_w3)          :: rho_e
    real(kind=r_def), dimension(ndf_wt)          :: theta_vd_e
    real(kind=r_def), dimension(ndf_w3)          :: ls_rho_e
    real(kind=r_def), dimension(ndf_wt)          :: ls_theta_vd_e
    real(kind=r_def)                             :: ls_exner

    real(kind=r_def) :: rho_cell, theta_vd_cell
    real(kind=r_def) :: ls_rho_cell, ls_theta_vd_cell

    do k = nlayers - 1, 0, -1
      do df = 1, ndf_w3, 1
        ls_rho_e(df) = ls_rho(map_w3(df) + k)
      end do

      do df = 1, ndf_wt, 1
        ls_theta_vd_e(df) = ls_moist_dyn_gas(k + map_wt(df)) * ls_theta(k + map_wt(df))
      end do

      theta_vd_e = 0.0_r_def
      rho_e = 0.0_r_def
      do df = ndf_w3, 1, -1
        ls_rho_cell = 0.0_r_def
        do df3 = 1, ndf_w3, 1
          ls_rho_cell = ls_rho_cell + ls_rho_e(df3) * w3_basis(1,df3,df)
        end do

        ls_theta_vd_cell = 0.0_r_def
        do dft = 1, ndf_wt, 1
          ls_theta_vd_cell = ls_theta_vd_cell + ls_theta_vd_e(dft) * wt_basis(1,dft,df)
        end do

        ls_exner = (ls_rho_cell * ls_theta_vd_cell * rd / p_zero) ** (kappa / (1.0_r_def - kappa))
        rho_cell = kappa * exner(map_w3(df) + k) * ls_exner / (-kappa * ls_rho_cell + 1.0_r_def * ls_rho_cell)
        theta_vd_cell = kappa * exner(map_w3(df) + k) * ls_exner / (-kappa * ls_theta_vd_cell + 1.0_r_def * ls_theta_vd_cell)
        exner(map_w3(df) + k) = 0.0_r_def

        do dft = ndf_wt, 1, -1
          theta_vd_e(dft) = theta_vd_e(dft) + theta_vd_cell * wt_basis(1,dft,df)
        end do

        do df3 = ndf_w3, 1, -1
          rho_e(df3) = rho_e(df3) + rho_cell * w3_basis(1,df3,df)
        end do
      end do ! df

      do df = ndf_wt, 1, -1
        theta(k + map_wt(df)) = theta(k + map_wt(df)) + ls_moist_dyn_gas(k + map_wt(df)) * theta_vd_e(df)
        moist_dyn_gas(k + map_wt(df)) = moist_dyn_gas(k + map_wt(df)) + ls_theta(k + map_wt(df)) * theta_vd_e(df)
      end do

      do df = ndf_w3, 1, -1
        rho(map_w3(df) + k) = rho(map_w3(df) + k) + rho_e(df)
      end do
    end do ! k

  end subroutine atl_sample_eos_pressure_code

end module atl_sample_eos_pressure_kernel_mod
