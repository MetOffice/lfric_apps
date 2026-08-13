!-----------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some of the content of this file has been produced with the assistance of
! Anthropic Claude Opus 5 (Claude Code).
!-----------------------------------------------------------------------------
!> @brief Mechanistic activation of the CASIM cloud droplet number.

module casim_arg_act_kernel_mod

use argument_mod,      only: arg_type,          &
                             GH_FIELD, GH_REAL, &
                             GH_READ,           &
                             GH_READWRITE,      &
                             CELL_COLUMN
use fs_continuity_mod, only: WTHETA
use kernel_mod,        only: kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: casim_arg_act_kernel_type
  private
  type(arg_type) :: meta_args(31) = (/                                      &
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! m_cl
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! m_cl_pre_fast
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! rho_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! theta
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! exner_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! cf_liq
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! cf_liq_pre_fast
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! w_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! wvar
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! nl_mphys
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! n_ait_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! ait_sol_su
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! ait_sol_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! ait_sol_om
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! n_acc_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! acc_sol_su
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! acc_sol_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! acc_sol_om
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! acc_sol_ss
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! n_cor_sol
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! cor_sol_su
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! cor_sol_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! cor_sol_om
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! cor_sol_ss
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! n_ait_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! ait_ins_bc
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! ait_ins_om
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! n_acc_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! acc_ins_du
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                   & ! n_cor_ins
       arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA)                    & ! cor_ins_du
       /)
   integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: casim_arg_act_code
end type

public :: casim_arg_act_code

contains

!> @brief Mechanistic activation of the CASIM cloud droplet number
!> @details Liquid cloud can be created or removed by processes outside the
!!          CASIM microphysics. For the mechanistic activation options, the
!!          cloud droplet number for any newly created cloud is worked out
!!          from the GLOMAP modal aerosol using the CASIM activation scheme,
!!          and the number is removed again where the cloud has gone.
!> @param[in]     nlayers          Number of layers
!> @param[in]     m_cl             Liquid cloud mass mixing ratio
!> @param[in]     m_cl_pre_fast    Liquid cloud mass mixing ratio as it was on
!!                                  entry to the fast physics
!> @param[in]     rho_in_wth       Dry density in potential temperature space
!> @param[in]     theta            Potential temperature field
!> @param[in]     exner_in_wth     Exner pressure in potential temperature space
!> @param[in]     cf_liq           Liquid cloud fraction
!> @param[in]     cf_liq_pre_fast  Liquid cloud fraction as it was on entry to
!!                                  the fast physics
!> @param[in]     w_in_wth         'Vertical' wind in theta space
!> @param[in]     wvar             Vertical velocity variance
!> @param[in,out] nl_mphys         CASIM cloud-droplet number concentration
!> @param[in]     n_ait_sol        Soluble Aitken mode number mixing ratio
!> @param[in]     ait_sol_su       Soluble Aitken mode H2SO4 mass mixing ratio
!> @param[in]     ait_sol_bc       Soluble Aitken mode black carbon m.m.r.
!> @param[in]     ait_sol_om       Soluble Aitken mode organic m.m.r.
!> @param[in]     n_acc_sol        Soluble accumulation mode number m.r.
!> @param[in]     acc_sol_su       Soluble accumulation mode H2SO4 m.m.r.
!> @param[in]     acc_sol_bc       Soluble accumulation mode black carbon m.m.r.
!> @param[in]     acc_sol_om       Soluble accumulation mode organic m.m.r.
!> @param[in]     acc_sol_ss       Soluble accumulation mode sea salt m.m.r.
!> @param[in]     n_cor_sol        Soluble coarse mode number mixing ratio
!> @param[in]     cor_sol_su       Soluble coarse mode H2SO4 mass mixing ratio
!> @param[in]     cor_sol_bc       Soluble coarse mode black carbon m.m.r.
!> @param[in]     cor_sol_om       Soluble coarse mode organic m.m.r.
!> @param[in]     cor_sol_ss       Soluble coarse mode sea salt m.m.r.
!> @param[in]     n_ait_ins        Insoluble Aitken mode number mixing ratio
!> @param[in]     ait_ins_bc       Insoluble Aitken mode black carbon m.m.r.
!> @param[in]     ait_ins_om       Insoluble Aitken mode organic m.m.r.
!> @param[in]     n_acc_ins        Insoluble accumulation mode number m.r.
!> @param[in]     acc_ins_du       Insoluble accumulation mode dust m.m.r.
!> @param[in]     n_cor_ins        Insoluble coarse mode number mixing ratio
!> @param[in]     cor_ins_du       Insoluble coarse mode dust m.m.r.
!> @param[in]     ndf_wth          Number of degrees of freedom per cell for
!!                                  potential temperature space
!> @param[in]     undf_wth         Number unique of degrees of freedom for
!!                                  potential temperature space
!> @param[in]     map_wth          Dofmap for the cell at the base of the
!!                                  column for potential temperature space

subroutine casim_arg_act_code( nlayers,                            &
                               m_cl, m_cl_pre_fast,                &
                               rho_in_wth, theta, exner_in_wth,    &
                               cf_liq, cf_liq_pre_fast,            &
                               w_in_wth, wvar,                     &
                               nl_mphys,                           &
                               n_ait_sol, ait_sol_su,              &
                               ait_sol_bc, ait_sol_om,             &
                               n_acc_sol, acc_sol_su,              &
                               acc_sol_bc, acc_sol_om,             &
                               acc_sol_ss,                         &
                               n_cor_sol, cor_sol_su,              &
                               cor_sol_bc, cor_sol_om,             &
                               cor_sol_ss,                         &
                               n_ait_ins, ait_ins_bc,              &
                               ait_ins_om,                         &
                               n_acc_ins, acc_ins_du,              &
                               n_cor_ins, cor_ins_du,              &
                               ndf_wth, undf_wth, map_wth          )

    use constants_mod,              only: r_def, i_def, r_um, i_um

    !---------------------------------------
    ! UM modules
    !---------------------------------------

    use planet_constants_mod,       only: p_zero, kappa
    use nlsizes_namelist_mod,       only: tr_ukca
    use mphys_inputs_mod,           only: wvarfac
    use casim_activation_in_um_mod, only: activate_column_ukca
    use casim_ukca_tracer_mod,      only: casim_ukca_tracer_column

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth
    integer(kind=i_def), intent(in) :: undf_wth

    real(kind=r_def), intent(in),  dimension(undf_wth) :: m_cl
    real(kind=r_def), intent(in),  dimension(undf_wth) :: m_cl_pre_fast
    real(kind=r_def), intent(in),  dimension(undf_wth) :: rho_in_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: theta
    real(kind=r_def), intent(in),  dimension(undf_wth) :: exner_in_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cf_liq
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cf_liq_pre_fast
    real(kind=r_def), intent(in),  dimension(undf_wth) :: w_in_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: wvar

    real(kind=r_def), intent(inout), dimension(undf_wth) :: nl_mphys

    real(kind=r_def), intent(in),  dimension(undf_wth) :: n_ait_sol
    real(kind=r_def), intent(in),  dimension(undf_wth) :: ait_sol_su
    real(kind=r_def), intent(in),  dimension(undf_wth) :: ait_sol_bc
    real(kind=r_def), intent(in),  dimension(undf_wth) :: ait_sol_om
    real(kind=r_def), intent(in),  dimension(undf_wth) :: n_acc_sol
    real(kind=r_def), intent(in),  dimension(undf_wth) :: acc_sol_su
    real(kind=r_def), intent(in),  dimension(undf_wth) :: acc_sol_bc
    real(kind=r_def), intent(in),  dimension(undf_wth) :: acc_sol_om
    real(kind=r_def), intent(in),  dimension(undf_wth) :: acc_sol_ss
    real(kind=r_def), intent(in),  dimension(undf_wth) :: n_cor_sol
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cor_sol_su
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cor_sol_bc
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cor_sol_om
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cor_sol_ss
    real(kind=r_def), intent(in),  dimension(undf_wth) :: n_ait_ins
    real(kind=r_def), intent(in),  dimension(undf_wth) :: ait_ins_bc
    real(kind=r_def), intent(in),  dimension(undf_wth) :: ait_ins_om
    real(kind=r_def), intent(in),  dimension(undf_wth) :: n_acc_ins
    real(kind=r_def), intent(in),  dimension(undf_wth) :: acc_ins_du
    real(kind=r_def), intent(in),  dimension(undf_wth) :: n_cor_ins
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cor_ins_du

    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth

    ! Local variables for the kernel
    ! Bounds on the activating vertical velocity, as used by the UM
    real(r_um), parameter :: min_velocity = 0.01_r_um ! ms-1
    real(r_um), parameter :: max_velocity = 4.00_r_um ! ms-1
    real(r_um), parameter :: eps = epsilon(1.0_r_um)

    real(r_um), dimension(nlayers) ::                                          &
         cloud_mass_post_qtbal, cloud_mass_pre_ap2, rho_col, t_col, p_col,     &
         cf_liquid, cf_liquid_pre_ap2, w_tke, cloudnumber_col

    ! UKCA tracer array holding the GLOMAP modes for this column
    real(r_um), allocatable :: tracer_ukca(:,:,:,:)

    integer(i_um) :: k

    !-------------------------------------------------------------------------
    ! End of Declarations
    !-------------------------------------------------------------------------

    allocate( tracer_ukca(1,1,0:nlayers,tr_ukca) )

    call casim_ukca_tracer_column( nlayers, undf_wth, map_wth(1),              &
                                   n_ait_sol, ait_sol_su, ait_sol_bc,          &
                                   ait_sol_om,                                 &
                                   n_acc_sol, acc_sol_su, acc_sol_bc,          &
                                   acc_sol_om, acc_sol_ss,                     &
                                   n_cor_sol, cor_sol_su, cor_sol_bc,          &
                                   cor_sol_om, cor_sol_ss,                     &
                                   n_ait_ins, ait_ins_bc, ait_ins_om,          &
                                   n_acc_ins, acc_ins_du,                      &
                                   n_cor_ins, cor_ins_du,                      &
                                   tracer_ukca )

    do k = 1, nlayers

      cloud_mass_post_qtbal(k) = m_cl(map_wth(1) + k)
      cloud_mass_pre_ap2(k)    = m_cl_pre_fast(map_wth(1) + k)
      rho_col(k)               = rho_in_wth(map_wth(1) + k)
      t_col(k)                 = exner_in_wth(map_wth(1) + k) *                &
                                 theta(map_wth(1) + k)
      p_col(k)                 = p_zero *                                      &
                          ( exner_in_wth(map_wth(1) + k) )**(1.0_r_um / kappa)
      cf_liquid(k)             = cf_liq(map_wth(1) + k)
      cf_liquid_pre_ap2(k)     = cf_liq_pre_fast(map_wth(1) + k)
      cloudnumber_col(k)       = nl_mphys(map_wth(1) + k)

      ! Activating vertical velocity, following the UM
      if ( wvar(map_wth(1) + k) > eps ) then
        w_tke(k) = w_in_wth(map_wth(1) + k) +                                  &
                   wvarfac * sqrt( wvar(map_wth(1) + k) )
      else
        w_tke(k) = w_in_wth(map_wth(1) + k)
      end if

      w_tke(k) = min( max_velocity, max( min_velocity, w_tke(k) ) )

    end do

    call activate_column_ukca( cloud_mass_post_qtbal,                          &
                               cloud_mass_pre_ap2,                             &
                               rho_col, t_col, p_col,                          &
                               cf_liquid, cf_liquid_pre_ap2, w_tke,            &
                               cloudnumber_col,                                &
                               tracer_ukca )

    do k = 1, nlayers
      nl_mphys(map_wth(1) + k) = cloudnumber_col(k)
    end do

    ! Set level 0 the same as level 1 (as done in the UM)
    nl_mphys(map_wth(1)) = nl_mphys(map_wth(1) + 1)

    deallocate( tracer_ukca )

end subroutine casim_arg_act_code

end module casim_arg_act_kernel_mod
