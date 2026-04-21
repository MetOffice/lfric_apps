module atl_sample_eos_pressure_kernel_mod
  use argument_mod, only : arg_type, cell_column, func_type, gh_basis, gh_evaluator, gh_field, gh_read, gh_readwrite, gh_real, &
&gh_write
  use constants_mod, only : i_def, r_def
  use fs_continuity_mod, only : w3, wtheta
  use kernel_mod, only : kernel_type
  use planet_config_mod, only : kappa, rd, p_zero
  implicit none
  type, public, extends(kernel_type) :: atl_sample_eos_pressure_kernel_type
  type(ARG_TYPE) :: META_ARGS(7) = (/ &
    arg_type(gh_field, gh_real, gh_readwrite, w3), &
    arg_type(gh_field, gh_real, gh_readwrite, w3), &
    arg_type(gh_field, gh_real, gh_readwrite, wtheta), &
    arg_type(gh_field, gh_real, gh_readwrite, wtheta), &
    arg_type(gh_field, gh_real, gh_read, w3), &
    arg_type(gh_field, gh_real, gh_read, wtheta), &
    arg_type(gh_field, gh_real, gh_read, wtheta)/)
  type(FUNC_TYPE) :: META_FUNCS(2) = (/ &
    func_type(w3, gh_basis), &
    func_type(wtheta, gh_basis)/)
  INTEGER :: GH_SHAPE = gh_evaluator
  INTEGER :: OPERATES_ON = cell_column
  CONTAINS
    PROCEDURE, NOPASS :: atl_sample_eos_pressure_code
END TYPE atl_sample_eos_pressure_kernel_type

  private

  public :: atl_sample_eos_pressure_code

  contains
  subroutine atl_sample_eos_pressure_code(nlayers, exner, rho, theta, moist_dyn_gas, ls_rho, ls_theta, ls_moist_dyn_gas, ndf_w3, &
&undf_w3, map_w3, w3_basis, basis_w3_on_wt, ndf_wt, undf_wt, map_wt, wt_basis, basis_wt_on_wt)
    use sci_coordinate_jacobian_mod, only : coordinate_jacobian
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wt
    integer(kind=i_def), intent(in) :: ndf_w3
    integer(kind=i_def), intent(in) :: undf_wt
    integer(kind=i_def), intent(in) :: undf_w3
    integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt
    integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
    real(kind=r_def), dimension(1,ndf_w3,ndf_w3), intent(in) :: w3_basis
    real(kind=r_def), dimension(1,ndf_w3,ndf_wt), intent(in) :: basis_w3_on_wt
    real(kind=r_def), dimension(1,ndf_wt,ndf_w3), intent(in) :: wt_basis
    real(kind=r_def), dimension(1,ndf_wt,ndf_wt), intent(in) :: basis_wt_on_wt
    real(kind=r_def), dimension(undf_w3), intent(inout) :: exner
    real(kind=r_def), dimension(undf_w3), intent(inout) :: rho
    real(kind=r_def), dimension(undf_wt), intent(inout) :: theta
    real(kind=r_def), dimension(undf_wt), intent(inout) :: moist_dyn_gas
    real(kind=r_def), dimension(undf_w3), intent(in) :: ls_rho
    real(kind=r_def), dimension(undf_wt), intent(in) :: ls_theta
    real(kind=r_def), dimension(undf_wt), intent(in) :: ls_moist_dyn_gas
    integer(kind=i_def) :: df
    integer(kind=i_def) :: df3
    integer(kind=i_def) :: dft
    integer(kind=i_def) :: k
    real(kind=r_def), dimension(ndf_w3) :: rho_e
    real(kind=r_def), dimension(ndf_wt) :: theta_vd_e
    real(kind=r_def), dimension(ndf_w3) :: ls_rho_e
    real(kind=r_def), dimension(ndf_wt) :: ls_theta_vd_e
    real(kind=r_def) :: rho_cell
    real(kind=r_def) :: theta_vd_cell
    real(kind=r_def) :: ls_rho_cell
    real(kind=r_def) :: ls_theta_vd_cell
    real(kind=r_def) :: tmp_ls_exner
    real(kind=r_def) :: tmp_exner

    tmp_exner = 0.0_r_def
    theta_vd_cell = 0.0_r_def
    theta_vd_e = 0.0_r_def
    rho_cell = 0.0_r_def
    rho_e = 0.0_r_def
    do k = nlayers - 1, 0, -1
      do df = 1, ndf_w3, 1
        ls_rho_e(df) = ls_rho(map_w3(df) + k)
      enddo
      do df = 1, ndf_wt, 1
        ls_theta_vd_e(df) = ls_moist_dyn_gas(k + map_wt(df)) * ls_theta(k + map_wt(df))
      enddo
      do df = ndf_w3, 1, -1
        ls_rho_cell = 0.0_r_def
        do df3 = 1, ndf_w3, 1
          ls_rho_cell = ls_rho_cell + ls_rho_e(df3) * w3_basis(1,df3,df)
        enddo
        ls_theta_vd_cell = 0.0_r_def
        do dft = 1, ndf_wt, 1
          ls_theta_vd_cell = ls_theta_vd_cell + ls_theta_vd_e(dft) * wt_basis(1,dft,df)
        enddo
        tmp_ls_exner = (ls_rho_cell * ls_theta_vd_cell * rd / p_zero) ** (kappa / (1.0 - kappa))
        tmp_exner = tmp_exner + exner(map_w3(df) + k)
        exner(map_w3(df) + k) = 0.0
        rho_cell = rho_cell + kappa * tmp_exner * tmp_ls_exner / (-kappa * ls_rho_cell + 1.0 * ls_rho_cell)
        theta_vd_cell = theta_vd_cell + kappa * tmp_exner * tmp_ls_exner / (-kappa * ls_theta_vd_cell + 1.0 * ls_theta_vd_cell)
        tmp_exner = 0.0
        do dft = ndf_wt, 1, -1
          theta_vd_e(dft) = theta_vd_e(dft) + theta_vd_cell * wt_basis(1,dft,df)
        enddo
        theta_vd_cell = 0.0
        do df3 = ndf_w3, 1, -1
          rho_e(df3) = rho_e(df3) + rho_cell * w3_basis(1,df3,df)
        enddo
        rho_cell = 0.0
      enddo
      do df = ndf_wt, 1, -1
        theta(k + map_wt(df)) = theta(k + map_wt(df)) + ls_moist_dyn_gas(k + map_wt(df)) * theta_vd_e(df)
        moist_dyn_gas(k + map_wt(df)) = moist_dyn_gas(k + map_wt(df)) + ls_theta(k + map_wt(df)) * theta_vd_e(df)
        theta_vd_e(df) = 0.0
      enddo
      do df = ndf_w3, 1, -1
        rho(map_w3(df) + k) = rho(map_w3(df) + k) + rho_e(df)
        rho_e(df) = 0.0
      enddo
    enddo

  end subroutine atl_sample_eos_pressure_code

end module atl_sample_eos_pressure_kernel_mod
