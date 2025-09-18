!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to cloud scheme.

module pc2_initiation_kernel_mod

use argument_mod,      only: arg_type,          &
                             GH_FIELD, GH_REAL, &
                             GH_INTEGER,        &
                             GH_READ, GH_WRITE, &
                             DOMAIN,            &
                             ANY_DISCONTINUOUS_SPACE_9, &
                             ANY_DISCONTINUOUS_SPACE_1
use fs_continuity_mod, only: WTHETA, W3
use kernel_mod,        only: kernel_type
use empty_data_mod,    only: empty_real_data

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: pc2_initiation_kernel_type
  private
  type(arg_type) :: meta_args(40) = (/                                   &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! mv_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! ml_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! mi_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! ms_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! cfl_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! cff_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! bcf_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! theta_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! exner_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        & ! exner_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! dsldzm
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! mix_len_bm
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! wvar
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! gradrinr
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! zh
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! zhsc
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! inv_depth
       arg_type(GH_FIELD, GH_INTEGER,GH_READ,ANY_DISCONTINUOUS_SPACE_9), & ! bl_type_ind
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! tau_dec_bm
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! tau_hom_bm
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! tau_mph_bm
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! mv_n_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! ml_n_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! theta_n_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! exner_n_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! zlcl_mixed
       arg_type(GH_FIELD, GH_INTEGER,GH_READ,ANY_DISCONTINUOUS_SPACE_1), & ! cumulus
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                    & ! height_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                    & ! dtheta_inc
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                    & ! dmv_inc_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                    & ! dmcl_inc_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                    & ! dmci_inc_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                    & ! dms_inc_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                    & ! dcfl_inc_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                    & ! dcff_inc_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                    & ! dbcf_inc_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                    & ! sskew_bm
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                    & ! svar_bm
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                    & ! svar_tb
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA)                     & ! rh_crit_wth
       /)
   integer :: operates_on = DOMAIN
contains
  procedure, nopass :: pc2_initiation_code
end type

public :: pc2_initiation_code

contains

!> @brief Interface to pc2 initiation
!> @details Calculates whether cloud should be created from clear-sky conditions
!>          or whether overcast conditions should be broken up.
!>          More info is in UMDP 30.
!> @param[in]     nlayers        Number of layers
!> @param[in]     mv_wth         Vapour mass mixing ratio
!> @param[in]     ml_wth         Liquid cloud mass mixing ratio
!> @param[in]     mi_wth         Ice cloud mass mixing ratio
!> @param[in]     ms_wth         Snow cloud mass mixing ratio
!> @param[in]     cfl_wth        Liquid cloud fraction
!> @param[in]     cff_wth        Ice cloud fraction
!> @param[in]     bcf_wth        Bulk cloud fraction
!> @param[in]     theta_wth      Potential temperature field
!> @param[in]     exner_wth      Exner pressure in theta space
!> @param[in]     exner_w3       Exner pressure in w3 space
!> @param[in]     dsldzm         Liquid potential temperature gradient in wth
!> @param[in]     mix_len_bm     Turb length-scale for bimodal in wth
!> @param[in]     wvar           Vertical velocity variance in wth
!> @param[in]     gradrinr       Gradient Richardson Number in wth
!> @param[in]     zh             Mixed-layer height
!> @param[in]     zhsc           Decoupled layer height
!> @param[in]     inv_depth      Depth of BL top inversion layer
!> @param[in]     bl_type_ind    Diagnosed BL types
!> @param[in]     tau_dec_bm     Decorrelation time scale in wth
!> @param[in]     tau_hom_bm     Homogenisation time scale in wth
!> @param[in]     tau_mph_bm     Phase-relaxation time scale in wth
!> @param[in]     mv_n_wth       Start of timestep vapour mass mixing ratio
!> @param[in]     ml_n_wth       Start of timestep liquid cloud mass mixing ratio
!> @param[in]     theta_n_wth    Start of timestep theta in theta space
!> @param[in]     exner_n_wth    Start of timestep exner in theta space
!> @param[in]     zlcl_mixed     The height of the lifting condensation level in the mixed layer
!> @param[in]     cumulus        The logical cumulus flag
!> @param[in]     height_wth     Height of wth levels above mean sea level
!> @param[in,out] dtheta_inc_wth Increment to theta in theta space
!> @param[in,out] dmv_inc_wth    Increment to water vapour in theta space
!> @param[in,out] dmcl_inc_wth   Increment to liquid water content in theta space
!> @param[in,out] dmci_inc_wth   Increment to ice water content in theta space
!> @param[in,out] dms_inc_wth    Increment to snow content in theta space
!> @param[in,out] dcfl_inc_wth   Increment to liquid cloud fraction in theta space
!> @param[in,out] dcff_inc_wth   Increment to ice cloud fraction in theta space
!> @param[in,out] dbcf_inc_wth   Increment to bulk cloud fraction in theta space
!> @param[in,out] sskew_bm       Bimodal skewness of SD PDF
!> @param[in,out] svar_bm        Bimodal variance of SD PDF
!> @param[in,out] svar_tb        Unimodal variance of SD PDF
!> @param[in]     rh_crit_wth    Critical relative humidity in theta space
!> @param[in]     ndf_wth        Number of degrees of freedom per cell for theta space
!> @param[in]     undf_wth       Number unique of degrees of freedom  for theta space
!> @param[in]     map_wth        Dofmap for the cell at the base of the column for theta space
!> @param[in]     ndf_w3         Number of degrees of freedom per cell for density space
!> @param[in]     undf_w3        Number unique of degrees of freedom  for density space
!> @param[in]     map_w3         Dofmap for the cell at the base of the column for density space
!> @param[in]     ndf_2d         Number of degrees of freedom per cell for density space
!> @param[in]     undf_2d        Number of unique degrees of freedom for density space
!> @param[in]     map_2d         Dofmap for the cell at the base of the column for density space
!> @param[in]     ndf_bl         Number of DOFs per cell for BL types
!> @param[in]     undf_bl        Number of total DOFs for BL types
!> @param[in]     map_bl         Dofmap for cell for BL types

subroutine pc2_initiation_code( nlayers, seg_len,                  &
                                mv_wth,                            &
                                ml_wth,                            &
                                mi_wth,                            &
                                ms_wth,                            &
                                cfl_wth,                           &
                                cff_wth,                           &
                                bcf_wth,                           &
                                theta_wth,                         &
                                exner_wth,                         &
                                exner_w3,                          &
                                dsldzm,                            &
                                mix_len_bm,                        &
                                wvar,                              &
                                gradrinr,                          &
                                zh,                                &
                                zhsc,                              &
                                inv_depth,                         &
                                bl_type_ind,                       &
                                tau_dec_bm,                        &
                                tau_hom_bm,                        &
                                tau_mph_bm,                        &
                                ! Start of timestep values for RHt
                                mv_n_wth,                          &
                                ml_n_wth,                          &
                                theta_n_wth,                       &
                                exner_n_wth,                       &
                                zlcl_mixed,                        &
                                cumulus,                           &
                                height_wth,                        &
                                ! Responses
                                dtheta_inc_wth,                    &
                                dmv_inc_wth,                       &
                                dmcl_inc_wth,                      &
                                dmci_inc_wth,                      &
                                dms_inc_wth,                       &
                                dcfl_inc_wth,                      &
                                dcff_inc_wth,                      &
                                dbcf_inc_wth,                      &
                                sskew_bm,                          &
                                svar_bm,                           &
                                svar_tb,                           &
                                rh_crit_wth,                       &
                                ! Other
                                ndf_wth, undf_wth, map_wth,        &
                                ndf_w3,  undf_w3,  map_w3,         &
                                ndf_2d,  undf_2d,  map_2d,         &
                                ndf_bl,  undf_bl,  map_bl)

    use constants_mod,    only: r_def, i_def, r_um, i_um
    use cloud_config_mod, only: i_bm_ez_opt, i_bm_ez_opt_entpar

    !---------------------------------------
    ! UM modules
    !---------------------------------------

    use pc2_initiation_ctl_mod,     only: pc2_initiation_ctl
    use planet_constants_mod,       only: p_zero, kappa, lcrcp, planet_radius
    use gen_phys_inputs_mod,        only: l_mr_physics

    use free_tracers_inputs_mod,    only: l_wtrac, n_wtrac
    use water_tracers_mod,          only: wtrac_type

    ! Redirect routine names to avoid clash with existing qsat routines
    use qsat_mod, only: qsat_wat_mix

    implicit none

    ! Arguments

    integer(kind=i_def), intent(in) :: nlayers, seg_len
    integer(kind=i_def), intent(in) :: ndf_wth , ndf_w3
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: ndf_bl, undf_bl

    real(kind=r_def), intent(in), dimension(undf_wth) :: mv_wth
    real(kind=r_def), intent(in), dimension(undf_wth) :: ml_wth
    real(kind=r_def), intent(in), dimension(undf_wth) :: mi_wth
    real(kind=r_def), intent(in), dimension(undf_wth) :: ms_wth
    real(kind=r_def), intent(in), dimension(undf_wth) :: bcf_wth
    real(kind=r_def), intent(in), dimension(undf_wth) :: cfl_wth
    real(kind=r_def), intent(in), dimension(undf_wth) :: cff_wth
    real(kind=r_def), intent(in), dimension(undf_wth) :: theta_wth
    real(kind=r_def), intent(in), dimension(undf_wth) :: exner_wth
    real(kind=r_def), intent(in), dimension(undf_w3)  :: exner_w3

    real(kind=r_def), intent(in), dimension(undf_wth) :: dsldzm
    real(kind=r_def), intent(in), dimension(undf_wth) :: mix_len_bm
    real(kind=r_def), intent(in), dimension(undf_wth) :: wvar
    real(kind=r_def), intent(in), dimension(undf_wth) :: gradrinr
    real(kind=r_def), intent(in), dimension(undf_wth) :: tau_dec_bm
    real(kind=r_def), intent(in), dimension(undf_wth) :: tau_hom_bm
    real(kind=r_def), intent(in), dimension(undf_wth) :: tau_mph_bm

    real(kind=r_def), intent(inout), dimension(:), pointer :: sskew_bm
    real(kind=r_def), intent(inout), dimension(:), pointer :: svar_bm
    real(kind=r_def), intent(inout), dimension(:), pointer :: svar_tb

    real(kind=r_def), intent(in), dimension(undf_2d)  :: zh
    real(kind=r_def), intent(in), dimension(undf_2d)  :: zhsc
    real(kind=r_def), intent(in), dimension(undf_2d)  :: inv_depth
    integer(kind=i_def), intent(in), dimension(undf_bl) :: bl_type_ind

    real(kind=r_def), intent(in), dimension(undf_2d) :: zlcl_mixed
    integer(kind=i_def), intent(in), dimension(undf_2d) :: cumulus
    real(kind=r_def), intent(in), dimension(undf_wth) :: height_wth

    ! Start of timestep values
    real(kind=r_def), intent(in),     dimension(undf_wth) :: mv_n_wth
    real(kind=r_def), intent(in),     dimension(undf_wth) :: ml_n_wth
    real(kind=r_def), intent(in),     dimension(undf_wth) :: theta_n_wth
    real(kind=r_def), intent(in),     dimension(undf_wth) :: exner_n_wth

    integer(kind=i_def), intent(in), dimension(ndf_wth,seg_len) :: map_wth
    integer(kind=i_def), intent(in), dimension(ndf_w3,seg_len)  :: map_w3
    integer(kind=i_def), intent(in), dimension(ndf_2d,seg_len)  :: map_2d
    integer(kind=i_def), intent(in), dimension(ndf_bl,seg_len)  :: map_bl

    ! The changes to the fields as a result
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dtheta_inc_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmv_inc_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmcl_inc_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmci_inc_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dms_inc_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dcfl_inc_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dcff_inc_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dbcf_inc_wth

    real(kind=r_def), intent(in), dimension(undf_wth) :: rh_crit_wth

    logical, dimension(seg_len,1) :: l_cumulus

    real(r_um), dimension(seg_len,1,nlayers) :: qv_work, qcl_work, qcf_work,   &
         cfl_work, cff_work, bcf_work, t_work, theta_work, rhts, t_incr,       &
         qv_incr, qcl_incr, qcf_incr, cfl_incr, cff_incr, bcf_incr, rhcpt,     &
         zeros, tgrad_in, mix_len_in, tau_dec_in, tau_hom_in, tau_mph_in,      &
         z_theta, wvar_in,                                                     &
         gradrinr_in, tlts, qtts, ptts, qsl_tl, p_theta_levels, sskew_out,     &
         svar_turb_out, svar_bm_out, qcf2_work, qcf2_incr

    real(r_um), dimension(seg_len,1) :: zh_in, zhsc_in, dzh_in, bl_type_7_in,  &
         p_star, zlcl_mix

    real(r_um), dimension(seg_len,1,nlayers+1) :: p_rho_levels

    real(r_um), dimension(seg_len,1,0:nlayers) :: r_theta_levels

    real(r_um) :: t_n

    integer(i_um) :: k, i

    ! Water tracer field which is not currently used but is required by
    ! UM routine
    type(wtrac_type), dimension(n_wtrac) :: wtrac

    ! Hardwired things for PC2
    !
    integer(i_um), parameter :: nSCMDpkgs=15
    logical,       parameter :: l_scmdiags(nscmdpkgs) = .false.
    logical,       parameter :: calculate_increments  = .false.

    ! Convert cumulus flag from integer to logical
    do i = 1, seg_len
      l_cumulus(i,1) = (cumulus(map_2d(1,i)) == 1_i_def)
      zlcl_mix(i,1) = zlcl_mixed(map_2d(1,i))
      do k = 0, nlayers
        r_theta_levels(i,1,k) = height_wth(map_wth(1,i)+k) + planet_radius
      end do
    end do

    do k = 1, nlayers
      do i = 1, seg_len
        zeros(i,1,k)=0.0_r_um
      end do
    end do
    !-----------------------------------------------------------------------
    ! Initialisation of prognostic variables and arrays
    !-----------------------------------------------------------------------

    do i = 1, seg_len
      p_star(i,1) = p_zero*(exner_wth(map_wth(1,i) + 0))     &
                                   **(1.0_r_def/kappa)
      do k = 1, nlayers
        ! Calculate temperature, this array will be updated.
        t_work(i,1,k)  = theta_wth(map_wth(1,i) + k) *                  &
                         exner_wth(map_wth(1,i) + k)

        ! Pressure at centre of theta levels
        p_theta_levels(i,1,k)    = p_zero*(exner_wth(map_wth(1,i) + k)) &
                                           **(1.0_r_def/kappa)

        ! pressure at layer boundaries
        p_rho_levels(i,1,k) = p_zero*( exner_w3(map_w3(1,i) + k-1))     &
                                           **(1.0_r_def/kappa)

        ! Bimodal cloud scheme inputs
        wvar_in(i,1,k)    = wvar(map_wth(1,i) + k)
        gradrinr_in(i,1,k)= gradrinr(map_wth(1,i) + k)
        tau_dec_in(i,1,k) = tau_dec_bm(map_wth(1,i) + k)
        tau_hom_in(i,1,k) = tau_hom_bm(map_wth(1,i) + k)
        tau_mph_in(i,1,k) = tau_mph_bm(map_wth(1,i) + k)
        z_theta(i,1,k) = height_wth(map_wth(1,i) + k) - height_wth(map_wth(1,i) + 0)

        ! Moist prognostics
        qv_work(i,1,k)   = mv_wth(map_wth(1,i) + k)
        qcl_work(i,1,k)  = ml_wth(map_wth(1,i) + k)
        qcf_work(i,1,k)  = ms_wth(map_wth(1,i) + k)
        qcf2_work(i,1,k)  = mi_wth(map_wth(1,i) + k)

        ! Critical relative humidity
        rhcpt(i,1,k)     = rh_crit_wth(map_wth(1,i) + k)

        ! Some start of timestep fields:
        ! Pressure on theta levels at start of time-step
        ptts(i,1,k) = p_zero*(exner_n_wth(map_wth(1,i) + k))     &
                                           **(1.0_r_def/kappa)

        ! Temperature at start of timestep
        t_n    = theta_n_wth(map_wth(1,i) + k) *                &
                      exner_n_wth(map_wth(1,i) + k)

        ! Q total at start of time-step
        qtts(i,1,k) = mv_n_wth(map_wth(1,i) + k) + ml_n_wth(map_wth(1,i) + k)

        ! Liquid temperature
        tlts(i,1,k) = t_n - ( lcrcp * ml_n_wth(map_wth(1,i) + k) )

      end do     ! k
    end do

    if (i_bm_ez_opt == i_bm_ez_opt_entpar) then
      ! Length-scale used for entraining parcel mode construction method
      do i = 1, seg_len
        do k = 1, nlayers
          mix_len_in(i,1,k)  = mix_len_bm(map_wth(1,i) + k)
        end do
      end do
    else
      ! SL-gradient used for stable-layer mode construction method
      do i = 1, seg_len
        do k = 1, nlayers
          tgrad_in(i,1,k)    = dsldzm(map_wth(1,i) + k)
        end do
      end do
    end if

    do i = 1, seg_len
      ! 2d fields for gradrinr-based entrainment zones
      zh_in(i,1)        = zh(map_2d(1,i))
      zhsc_in(i,1)      = zhsc(map_2d(1,i))
      dzh_in(i,1)       = real(inv_depth(map_2d(1,i)), r_um)
      bl_type_7_in(i,1) = bl_type_ind(map_bl(1,i)+6)
    end do

    ! Calculate qsat(TL) with respect to liquid water, operate on whole column.
    call qsat_wat_mix(qsl_tl, tlts, ptts, seg_len, 1, nlayers )
    do k = 1, nlayers
      do i = 1, seg_len
        ! Total relative humidity (using start of time-step liquid temperature).
        rhts(i,1,k) = qtts(i,1,k) / qsl_tl(i,1,k)
      end do
    end do

    do i = 1, seg_len
      do k = 1, nlayers
        ! Recast LFRic cloud fractions onto cloud fraction work arrays.
        bcf_work(i,1,k) = bcf_wth(map_wth(1,i) + k)
        cfl_work(i,1,k) = cfl_wth(map_wth(1,i) + k)
        cff_work(i,1,k) = cff_wth(map_wth(1,i) + k)
      end do
    end do

    ! Initialize
    do k = 1, nlayers
      do i = 1, seg_len
        t_incr(i,1,k)   = 0.0_r_um
        qv_incr(i,1,k)  = 0.0_r_um
        qcl_incr(i,1,k) = 0.0_r_um
        qcf_incr(i,1,k) = 0.0_r_um
        qcf2_incr(i,1,k) = 0.0_r_um
        bcf_incr(i,1,k) = 0.0_r_um
        cfl_incr(i,1,k) = 0.0_r_um
        cff_incr(i,1,k) = 0.0_r_um
      end do
    end do

    call pc2_initiation_ctl(                               &
                            ! Dimensions of Rh crit array
                            seg_len,                       &
                            1,                             &
                            ! Pass in zlcl_mix
                            zlcl_mix,                      &
                            r_theta_levels,                &
                            ! Model switches
                            l_mr_physics, l_wtrac,         &
                            ! SCM diagnostics switches
                            nSCMDpkgs,                     &
                            L_SCMDiags,                    &
                            ! Primary fields passed IN/OUT
                            t_work,                        &
                            qv_work,                       &
                            qcl_work,                      &
                            qcf_work,                      &
                            qcf2_work,                     &
                            bcf_work,                      &
                            cfl_work,                      &
                            cff_work,                      &
                            rhts,                          &
                            tlts,                          &
                            qtts,                          &
                            ptts,                          &
                            zeros,                         &
                            ! Primary fields passed IN
                            p_rho_levels,                  &
                            p_star,                        &
                            p_theta_levels,                &
                            l_cumulus,                     &
                            rhcpt,                         &
                            ! Bimodal inputs
                            tgrad_in,                      &
                            wvar_in,                       &
                            tau_dec_in,                    &
                            tau_hom_in,                    &
                            tau_mph_in,                    &
                            z_theta,                       &
                            gradrinr_in,                   &
                            mix_len_in,                    &
                            zh_in,                         &
                            zhsc_in,                       &
                            dzh_in,                        &
                            bl_type_7_in,                  &
                            calculate_increments,          &
                            ! Output increments
                            t_incr,                        &
                            qv_incr,                       &
                            qcl_incr,                      &
                            qcf_incr,                      &
                            qcf2_incr,                     &
                            bcf_incr,                      &
                            cfl_incr,                      &
                            cff_incr,                      &
                            sskew_out,                     &
                            svar_turb_out,                 &
                            svar_bm_out,                   &
                            wtrac)

    ! Recast back to LFRic space
    do k = 1, nlayers
      do i = 1, seg_len
        ! *_work arrays have been updated

        ! New theta found from new temperature.
        theta_work(i,1,k) = t_work(i,1,k) /                  &
                          exner_wth(map_wth(1,i) + k)
        ! All increments found from difference between the updated *_work values
        ! and the values that were intent in.
        dtheta_inc_wth(map_wth(1,i) + k) = theta_work(i,1,k)   &
                                       - theta_wth(map_wth(1,i) + k)
        !
        dmv_inc_wth (map_wth(1,i)+k) = qv_work(i,1,k)  - mv_wth(map_wth(1,i) + k)
        dmcl_inc_wth(map_wth(1,i)+k) = qcl_work(i,1,k) - ml_wth(map_wth(1,i) + k)
        dmci_inc_wth(map_wth(1,i)+k) = qcf2_work(i,1,k) - mi_wth(map_wth(1,i) + k)
        dms_inc_wth(map_wth(1,i)+k) = qcf_work(i,1,k) - ms_wth(map_wth(1,i) + k)
        dcfl_inc_wth(map_wth(1,i)+k) = cfl_work(i,1,k) - cfl_wth(map_wth(1,i) + k)
        dcff_inc_wth(map_wth(1,i)+k) = cff_work(i,1,k) - cff_wth(map_wth(1,i) + k)
        dbcf_inc_wth(map_wth(1,i)+k) = bcf_work(i,1,k) - bcf_wth(map_wth(1,i) + k)
      end do
    end do
    if (.not. associated(sskew_bm, empty_real_data)) then
      do k = 1, nlayers
        do i = 1, seg_len
          sskew_bm(map_wth(1,i)+k)     = sskew_out(i,1,k)
        end do
      end do
    end if
    if (.not. associated(svar_bm, empty_real_data)) then
      do k = 1, nlayers
        do i = 1, seg_len
          svar_bm(map_wth(1,i)+k)      = svar_bm_out(i,1,k)
        end do
      end do
    end if
    if (.not. associated(svar_tb, empty_real_data)) then
      do k = 1, nlayers
        do i = 1, seg_len
          svar_tb(map_wth(1,i)+k)      = svar_turb_out(i,1,k)
        end do
      end do
    end if
    do i = 1, seg_len
      dtheta_inc_wth(map_wth(1,i)+0) = dtheta_inc_wth(map_wth(1,i)+1)
      dmv_inc_wth   (map_wth(1,i)+0) = dmv_inc_wth   (map_wth(1,i)+1)
      dmcl_inc_wth  (map_wth(1,i)+0) = dmcl_inc_wth  (map_wth(1,i)+1)
      dmci_inc_wth  (map_wth(1,i)+0) = dmci_inc_wth  (map_wth(1,i)+1)
      dms_inc_wth   (map_wth(1,i)+0) = dms_inc_wth   (map_wth(1,i)+1)
      dcfl_inc_wth  (map_wth(1,i)+0) = dcfl_inc_wth  (map_wth(1,i)+1)
      dcff_inc_wth  (map_wth(1,i)+0) = dcff_inc_wth  (map_wth(1,i)+1)
      dbcf_inc_wth  (map_wth(1,i)+0) = dbcf_inc_wth  (map_wth(1,i)+1)
    end do

end subroutine pc2_initiation_code

end module pc2_initiation_kernel_mod
