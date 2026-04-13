module atl_moist_dyn_mass_kernel_mod
  use argument_mod, only : arg_type, cell_column, gh_field, gh_read, gh_readwrite, gh_real, gh_scalar, gh_write
  use constants_mod, only : i_def, r_def
  use fs_continuity_mod, only : wtheta
  use kernel_mod, only : kernel_type
  implicit none
  type, public, extends(kernel_type) :: atl_moist_dyn_mass_kernel_type
  type(ARG_TYPE) :: META_ARGS(2) = (/ &
    arg_type(gh_field, gh_real, gh_readwrite, wtheta), &
    arg_type(gh_field*6, gh_real, gh_readwrite, wtheta)/)
  INTEGER :: OPERATES_ON = cell_column
  CONTAINS
    PROCEDURE, NOPASS :: atl_moist_dyn_mass_code
END TYPE atl_moist_dyn_mass_kernel_type

  private

  public :: atl_moist_dyn_mass_code

  contains
  subroutine atl_moist_dyn_mass_code(nlayers, moist_dyn_tot, mr_v, mr_cl, mr_r, mr_s, mr_g, mr_ci, ndf_wtheta, undf_wtheta, &
&map_wtheta)
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wtheta
    integer(kind=i_def), intent(in) :: undf_wtheta
    real(kind=r_def), dimension(undf_wtheta), intent(inout) :: moist_dyn_tot
    real(kind=r_def), dimension(undf_wtheta), intent(inout) :: mr_v
    real(kind=r_def), dimension(undf_wtheta), intent(inout) :: mr_cl
    real(kind=r_def), dimension(undf_wtheta), intent(inout) :: mr_r
    real(kind=r_def), dimension(undf_wtheta), intent(inout) :: mr_ci
    real(kind=r_def), dimension(undf_wtheta), intent(inout) :: mr_s
    real(kind=r_def), dimension(undf_wtheta), intent(inout) :: mr_g
    integer(kind=i_def), dimension(ndf_wtheta), intent(in) :: map_wtheta
    integer(kind=i_def) :: k
    integer(kind=i_def) :: df
    real(kind=r_def) :: mr_v_at_dof
    real(kind=r_def) :: mr_cl_at_dof
    real(kind=r_def) :: mr_r_at_dof
    real(kind=r_def) :: mr_ci_at_dof
    real(kind=r_def) :: mr_s_at_dof
    real(kind=r_def) :: mr_g_at_dof

    mr_v_at_dof = 0.0_r_def
    mr_cl_at_dof = 0.0_r_def
    mr_r_at_dof = 0.0_r_def
    mr_s_at_dof = 0.0_r_def
    mr_g_at_dof = 0.0_r_def
    mr_ci_at_dof = 0.0_r_def
    do k = nlayers - 1, 0, -1
      do df = ndf_wtheta, 1, -1
        mr_ci_at_dof = mr_ci_at_dof + moist_dyn_tot(map_wtheta(df) + k)
        mr_cl_at_dof = mr_cl_at_dof + moist_dyn_tot(map_wtheta(df) + k)
        mr_g_at_dof = mr_g_at_dof + moist_dyn_tot(map_wtheta(df) + k)
        mr_r_at_dof = mr_r_at_dof + moist_dyn_tot(map_wtheta(df) + k)
        mr_s_at_dof = mr_s_at_dof + moist_dyn_tot(map_wtheta(df) + k)
        mr_v_at_dof = mr_v_at_dof + moist_dyn_tot(map_wtheta(df) + k)
        moist_dyn_tot(map_wtheta(df) + k) = 0.0
        mr_g(map_wtheta(df) + k) = mr_g(map_wtheta(df) + k) + mr_g_at_dof
        mr_g_at_dof = 0.0
        mr_s(map_wtheta(df) + k) = mr_s(map_wtheta(df) + k) + mr_s_at_dof
        mr_s_at_dof = 0.0
        mr_ci(map_wtheta(df) + k) = mr_ci(map_wtheta(df) + k) + mr_ci_at_dof
        mr_ci_at_dof = 0.0
        mr_r(map_wtheta(df) + k) = mr_r(map_wtheta(df) + k) + mr_r_at_dof
        mr_r_at_dof = 0.0
        mr_cl(map_wtheta(df) + k) = mr_cl(map_wtheta(df) + k) + mr_cl_at_dof
        mr_cl_at_dof = 0.0
        mr_v(map_wtheta(df) + k) = mr_v(map_wtheta(df) + k) + mr_v_at_dof
        mr_v_at_dof = 0.0
      enddo
    enddo

  end subroutine atl_moist_dyn_mass_code

end module atl_moist_dyn_mass_kernel_mod
