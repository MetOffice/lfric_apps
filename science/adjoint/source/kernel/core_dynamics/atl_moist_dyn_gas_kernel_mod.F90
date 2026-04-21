module atl_moist_dyn_gas_kernel_mod
  use argument_mod, only : arg_type, cell_column, gh_field, gh_read, gh_readwrite, gh_real, gh_scalar, gh_write
  use constants_mod, only : i_def, r_def
  use fs_continuity_mod, only : wtheta
  use kernel_mod, only : kernel_type
  use planet_config_mod, only : recip_epsilon
  implicit none
  type, public, extends(kernel_type) :: atl_moist_dyn_gas_kernel_type
  type(ARG_TYPE) :: META_ARGS(2) = (/ &
    arg_type(gh_field, gh_real, gh_readwrite, wtheta), &
    arg_type(gh_field, gh_real, gh_readwrite, wtheta)/)
  INTEGER :: OPERATES_ON = cell_column
  CONTAINS
    PROCEDURE, NOPASS :: atl_moist_dyn_gas_code
END TYPE atl_moist_dyn_gas_kernel_type

  private

  public :: atl_moist_dyn_gas_code

  contains
  subroutine atl_moist_dyn_gas_code(nlayers, moist_dyn_gas, mr_v, ndf_wtheta, undf_wtheta, map_wtheta)
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wtheta
    integer(kind=i_def), intent(in) :: undf_wtheta
    real(kind=r_def), dimension(undf_wtheta), intent(inout) :: moist_dyn_gas
    real(kind=r_def), dimension(undf_wtheta), intent(inout) :: mr_v
    integer(kind=i_def), dimension(ndf_wtheta), intent(in) :: map_wtheta
    integer(kind=i_def) :: k
    integer(kind=i_def) :: df
    real(kind=r_def) :: mr_v_at_dof

    mr_v_at_dof = 0.0_r_def
    do k = nlayers - 1, 0, -1
      do df = ndf_wtheta, 1, -1
        mr_v_at_dof = mr_v_at_dof + moist_dyn_gas(map_wtheta(df) + k) * recip_epsilon
        moist_dyn_gas(map_wtheta(df) + k) = 0.0
        mr_v(map_wtheta(df) + k) = mr_v(map_wtheta(df) + k) + mr_v_at_dof
        mr_v_at_dof = 0.0
      enddo
    enddo

  end subroutine atl_moist_dyn_gas_code

end module atl_moist_dyn_gas_kernel_mod
