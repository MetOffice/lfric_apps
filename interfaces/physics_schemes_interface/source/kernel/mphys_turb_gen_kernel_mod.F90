!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to UM code for mixedphase turbulence cloud generation

module mphys_turb_gen_kernel_mod

use argument_mod,         only: arg_type, GH_REAL,   &
                                GH_FIELD, GH_READ,   &
                                GH_READWRITE, DOMAIN
use fs_continuity_mod,    only: WTHETA, W3
use kernel_mod,           only: kernel_type
use constants_mod,        only: r_def, i_def, r_um, i_um

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
type, public, extends(kernel_type) :: mphys_turb_gen_kernel_type
  private
  type(arg_type) :: meta_args(22) = (/           &
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! theta
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! mr_v
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! mr_cl
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! mr_ci
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! mr_s
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! ns_mphys
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! ni_mphys
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! cfl
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! cff
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! bcf
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! exner
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! rho
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! wvar
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! dz_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! dtheta
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! dmr_v
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! dmr_cl
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! dmr_ci
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! dmr_s
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! dcfl
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! dcff
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA)  & ! dbcf
       /)
  integer :: operates_on = DOMAIN
contains
  procedure, nopass :: mphys_turb_gen_code
end type mphys_turb_gen_kernel_type

public mphys_turb_gen_code
contains

!> @brief Call the UM code for mixedphase turbulence cloud generation
!> @param[in]     nlayers    Number of layers
!> @param[in]     theta      Potential temperature (K)
!> @param[in]     mr_v       Vapour mixing ratio (kg/kg)
!> @param[in]     mr_cl      Liquid mixing ratio (kg/kg)
!> @param[in]     mr_ci      Ice mixing ratio (kg/kg)
!> @param[in]     mr_s       Snow mixing ratio (kg/kg)
!> @param[in]     ns_mphys   Cloud ice number mixing ratio in wth
!> @param[in]     ni_mphys   Snow number mixing ratio in wth
!> @param[in]     cfl        Liquid cloud fraction
!> @param[in]     cff        Frozen cloud fraction
!> @param[in]     bcf        Total cloud fraction
!> @param[in]     exner      Exner pressure (Pa)
!> @param[in]     rho        Dry density (kg/m3)
!> @param[in]     wvar       Vertical velocity variance
!> @param[in]     dz_in_wth  Depth of wth space
!> @param[in,out] dtheta     Potential temperature increment from microphysics
!> @param[in,out] dmr_v      Vapour mixing ratio increment from microphysics
!> @param[in,out] dmr_cl     Liquid mixing ratio increment from microphysics
!> @param[in,out] dmr_ci     Ice mixing ratio increment from microphysics
!> @param[in,out] dmr_s      Snow mixing ratio increment from microphysics
!> @param[in,out] dcfl       Liquid cloud fraction inc from microphysics
!> @param[in,out] dcff       Ice cloud fraction increment from microphysics
!> @param[in,out] dbcf       Total cloud fraction inc from microphysics
!> @param[in]     ndf_wth    Number of dofs per cell for theta space
!> @param[in]     undf_wth   Number of unique degrees of freedom in segment for theta space
!> @param[in]     map_wth    Dofmap for a segment of the theta space

  subroutine mphys_turb_gen_code( nlayers,   &
                                  seg_len,   &
                                  theta,     &
                                  mr_v,      &
                                  mr_cl,     &
                                  mr_ci,     &
                                  mr_s,      &
                                  ns_mphys,  &
                                  ni_mphys,  &
                                  cfl,       &
                                  cff,       &
                                  bcf,       &
                                  exner,     &
                                  rho,       &
                                  wvar,      &
                                  dz_in_wth, &
                                  dtheta,    &
                                  dmr_v,     &
                                  dmr_cl,    &
                                  dmr_ci,    &
                                  dmr_s,     &
                                  dcfl,      &
                                  dcff,      &
                                  dbcf,      &
                                  ndf_wth, undf_wth, map_wth)

    !---------------------------------------
    ! UM modules
    !---------------------------------------
    use mphys_turb_gen_mixed_phase_mod, only: mphys_turb_gen_mixed_phase
    use nlsizes_namelist_mod, only: bl_levels
    use planet_constants_mod, only: p_zero, kappa
    use microphysics_config_mod, only: microphysics_casim

    implicit none

    integer(i_def), intent(in) :: nlayers, seg_len
    integer(i_def), intent(in) :: ndf_wth,  undf_wth
    integer(i_def), intent(in), dimension(ndf_wth, seg_len) :: map_wth

    real(r_def), intent(in), dimension(undf_wth) :: theta
    real(r_def), intent(in), dimension(undf_wth) :: mr_v
    real(r_def), intent(in), dimension(undf_wth) :: mr_cl
    real(r_def), intent(in), dimension(undf_wth) :: mr_ci
    real(r_def), intent(in), dimension(undf_wth) :: mr_s
    real(r_def), intent(in), dimension(undf_wth) :: ns_mphys
    real(r_def), intent(in), dimension(undf_wth) :: ni_mphys
    real(r_def), intent(in), dimension(undf_wth) :: cfl
    real(r_def), intent(in), dimension(undf_wth) :: cff
    real(r_def), intent(in), dimension(undf_wth) :: bcf
    real(r_def), intent(in), dimension(undf_wth) :: exner
    real(r_def), intent(in), dimension(undf_wth) :: rho
    real(r_def), intent(in), dimension(undf_wth) :: wvar
    real(r_def), intent(in), dimension(undf_wth) :: dz_in_wth

    real(r_def), intent(inout), dimension(undf_wth) :: dtheta
    real(r_def), intent(inout), dimension(undf_wth) :: dmr_v
    real(r_def), intent(inout), dimension(undf_wth) :: dmr_cl
    real(r_def), intent(inout), dimension(undf_wth) :: dmr_ci
    real(r_def), intent(inout), dimension(undf_wth) :: dmr_s
    real(r_def), intent(inout), dimension(undf_wth) :: dcfl
    real(r_def), intent(inout), dimension(undf_wth) :: dcff
    real(r_def), intent(inout), dimension(undf_wth) :: dbcf

    ! local variables
    real(r_um), dimension(seg_len,1,nlayers) :: &
                                      q_work, qcl_work, qcf_work, t_work,      &
                                      qcf2_work, cff_work, cfl_work, cf_work,  &
                                      bl_w_var, rhodz_dry, rhodz_moist, deltaz,&
                                      t_inc, dqcl_mp, qcl_mpt, tau_d, inv_prt, &
                                      disprate, inv_mt, si_avg, dcfl_mp,       &
                                      sigma2_s

    real(r_um), dimension(seg_len,1,0:nlayers) :: &
                                        q_n, cfl_n, cf_n, press_wth, &
                                        q_inc, qcl_inc, cfl_inc, cf_inc, &
                                        icenumber, snownumber

    integer(i_um) :: i, j, k

    do i = 1, seg_len
      do k = 1, nlayers
        ! Only single ice, no numbers required
        qcf_work(i,1,k) = mr_s(map_wth(1,i) + k) + dmr_s(map_wth(1,i) + k)
      end do
    end do
    if (microphysics_casim) then
      do i = 1, seg_len
        do k = 1, nlayers
          ! Set ice and snow number, and qcf2
          snownumber(i,1,k) = ns_mphys(map_wth(1,i) + k)
          icenumber(i,1,k) = ni_mphys(map_wth(1,i) + k)
          qcf2_work(i,1,k) = mr_ci(map_wth(1,i) + k) + dmr_ci(map_wth(1,i) + k)
        end do
      end do
    end if

    j = 1
    do i = 1, seg_len
      do k = 1, nlayers
        ! Time level n quantities
        q_n(i,j,k) = mr_v(map_wth(1,i) + k)
        cfl_n(i,j,k) = cfl(map_wth(1,i) + k)
        cf_n(i,j,k) = bcf(map_wth(1,i) + k)

        press_wth(i,j,k) = p_zero * (exner(map_wth(1,i) + k)) ** (1.0_r_um/kappa)
        rhodz_dry(i,j,k) = rho(map_wth(1,i) + k)

        ! Updated values after microphysics
        t_work(i,j,k) = (theta(map_wth(1,i) + k) + dtheta(map_wth(1,i) + k)) * &
                    exner(map_wth(1,i) + k)
        q_work(i,j,k) = mr_v(map_wth(1,i) + k) + dmr_v(map_wth(1,i) + k)
        qcl_work(i,j,k) = mr_cl(map_wth(1,i) + k) + dmr_cl(map_wth(1,i) + k)
        cff_work(i,j,k) = cff(map_wth(1,i) + k) + dcff(map_wth(1,i) + k)
        cfl_work(i,j,k) = cfl(map_wth(1,i) + k) + dcfl(map_wth(1,i) + k)
        cf_work(i,j,k) = bcf(map_wth(1,i) + k) + dbcf(map_wth(1,i) + k)

        !Increments from microphysics
        t_inc(i,j,k) = dtheta(map_wth(1,i) + k) * exner(map_wth(1,i) + k)
        q_inc(i,j,k) = dmr_v(map_wth(1,i) + k)
        qcl_inc(i,j,k) = dmr_cl(map_wth(1,i) + k)
        cfl_inc(i,j,k) = dcfl(map_wth(1,i) + k)
        cf_inc(i,j,k) = dbcf(map_wth(1,i) + k)

        ! Vertical velocity variance
        bl_w_var(i,j,k) = wvar(map_wth(1,i) + k)

        ! Layer depth
        deltaz(i,j,k) = dz_in_wth(map_wth(1,i) + k)
      end do
    end do

    j = 1
    do i = 1, seg_len
      ! Lowest layer needs level 0 adding
      deltaz(i,j,1) = deltaz(i,j,1) + dz_in_wth(map_wth(1,i))
    end do

    rhodz_dry = rhodz_dry * deltaz

    call mphys_turb_gen_mixed_phase( q_work, t_work, qcl_work, qcf_work,     &
                                     q_inc, qcl_inc, cfl_inc,  cf_inc,       &
                                     t_inc,  dqcl_mp, bl_levels,             &
                                     bl_w_var, cff_work, cfl_work, cf_work,  &
                                     q_n, cfl_n, cf_n, press_wth,            &
                                     rhodz_dry, rhodz_moist, deltaz,         &
                                     qcl_mpt, tau_d, inv_prt, disprate,      &
                                     inv_mt, si_avg, dcfl_mp, sigma2_s,      &
                                     qcf2_work, icenumber, snownumber )

    j = 1
    do i = 1, seg_len
      do k = 1, nlayers
        dtheta(map_wth(1,i) + k) = t_inc(i,j,k) / exner(map_wth(1,i) + k)
        dmr_v(map_wth(1,i) + k) = q_inc(i,j,k)
        dmr_cl(map_wth(1,i) + k) = qcl_inc(i,j,k)
        dcfl(map_wth(1,i) + k) = cfl_inc(i,j,k)
        dbcf(map_wth(1,i) + k) = cf_inc(i,j,k)
      end do
      dtheta(map_wth(1,i) + 0) = dtheta(map_wth(1,i) + 1)
      dmr_v(map_wth(1,i) + 0)   = dmr_v(map_wth(1,i) + 1)
      dmr_cl(map_wth(1,i) + 0)   = dmr_cl(map_wth(1,i) + 1)
      dbcf(map_wth(1,i) + 0) = dbcf(map_wth(1,i) + 1)
      dcfl(map_wth(1,i) + 0) = dcfl(map_wth(1,i) + 1)
    end do

  end subroutine mphys_turb_gen_code

end module mphys_turb_gen_kernel_mod
