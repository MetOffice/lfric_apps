!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to microphysics scheme.

module mphys_kernel_mod

use argument_mod,            only: arg_type,                  &
                                   GH_FIELD, GH_REAL,         &
                                   GH_READ, GH_WRITE,         &
                                   GH_READWRITE,              &
                                   ANY_DISCONTINUOUS_SPACE_1, &
                                   ANY_DISCONTINUOUS_SPACE_2, &
                                   DOMAIN
use fs_continuity_mod,       only: WTHETA, W3
use kernel_mod,              only: kernel_type
use empty_data_mod,          only: empty_real_data
use microphysics_config_mod, only: prog_tnuc
use aerosol_config_mod,      only: murk_prognostic

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel.
!> Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: mphys_kernel_type
  private
  type(arg_type) :: meta_args(48) = (/                                      &
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! mv_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! ml_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! ms_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! mr_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! mg_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! cf_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! cfl_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! cff_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                           & ! u_in_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                           & ! v_in_w3,
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! theta_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! exner_in_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                           & ! wetrho_in_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                           & ! dry_rho_in_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                           & ! height_w3
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! height_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                       & ! cloud_drop_no_conc
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE,  WTHETA),                  & ! murk
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1),    & ! sd_orog
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! dmv_wth
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! dml_wth
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! dms_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! dmr_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! dmg_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! ls_rain_2d
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! ls_snow_2d
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! ls_graup_2d
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),    & ! lsca_2d
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! precfrac
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! ls_rain_3d
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! ls_snow_3d
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! autoconv
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! accretion
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! rim_cry
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! rim_agg
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! dtheta
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! dcfl_wth
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! dcff_wth
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                   & ! dbcf_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2),    & ! f_arr_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ, WTHETA),                        & ! tnuc
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! superc_liq_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! superc_rain_wth
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! sfwater
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! sfrain
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! sfsnow
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                       & ! refl_tot
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1)     & ! refl_1km
       /)
   integer :: operates_on = DOMAIN
contains
  procedure, nopass :: mphys_code
end type

public :: mphys_code

contains

!> @brief Interface to the microphysics scheme
!> @param[in]     nlayers             Number of layers
!> @param[in]     seg_len             Number of horizontal cells in segment
!> @param[in]     mv_wth              Vapour mass mixing ratio
!> @param[in]     ml_wth              Liquid cloud mass mixing ratio
!> @param[in]     ms_wth              Snow cloud mass mixing ratio
!> @param[in]     mr_wth              Rain mass mixing ratio
!> @param[in]     mg_wth              Graupel mass mixing ratio
!> @param[in]     cf_wth              Bulk cloud fraction
!> @param[in]     cfl_wth             Liquid cloud fraction
!> @param[in]     cff_wth             Ice cloud fraction
!> @param[in]     u_in_w3             'Zonal' wind in density space
!> @param[in]     v_in_w3             'Meridional' wind in density space
!> @param[in]     theta_in_wth        Potential temperature field
!> @param[in]     exner_in_wth        Exner pressure in potential temperature space
!> @param[in]     wetrho_in_w3        Wet density in density space
!> @param[in]     dry_rho_in_w3       Dry density in density space
!> @param[in]     height_w3           Height of density space levels above surface
!> @param[in]     height_wth          Height of potential temperature space levels
!!                                     above surface
!> @param[in]     cloud_drop_no_conc  Cloud Droplet Number Concentration
!> @param[in]     sd_orog             Standard deviation of sub-grid orog
!> @param[in,out] dmv_wth             Increment to vapour mass mixing ratio
!> @param[in,out] dml_wth             Increment to liquid cloud mass mixing ratio
!> @param[in,out] dms_wth             Increment to snow cloud mass mixing ratio
!> @param[in,out] dmr_wth             Increment to rain mass mixing ratio
!> @param[in,out] dmg_wth             Increment to graupel mass mixing ratio
!> @param[in,out] ls_rain_2d          Large scale rain from twod_fields
!> @param[in,out] ls_snow_2d          Large scale snow from twod_fields
!> @param[in,out] ls_graup_2d         Large scale graupel from twod_fields
!> @param[in,out] lsca_2d             Large scale cloud amount (2d)
!> @param[in,out] precfrac            3D precipitation fraction
!> @param[in,out] ls_rain_3d          Large scale rain from 3d fields (kg m-2 s-1)
!> @param[in,out] ls_snow_3d          Large scale snow from 3d fields (kg m-2 s-1)
!> @param[in,out] autoconv            Rain autoconversion rate (kg kg-1 s-1)
!> @param[in,out] accretion           Rain accretion rate (kg kg-1 s-1)
!> @param[in,out] rim_cry             Riming rate for ice crystals (kg kg-1 s-1)
!> @param[in,out] rim_agg             Riming rate for ice aggregates (kg kg-1 s-1)
!> @param[in,out] dtheta              Increment to theta
!> @param[in,out] dcfl_wth            Increment to liquid cloud fraction
!> @param[in,out] dcff_wth            Increment to ice cloud fraction
!> @param[in,out] dbcf_wth            Increment to bulk cloud fraction
!> @param[in]     f_arr_wth           Parameters related to fractional standard
!!                                     deviation of condensate
!> @param[in]     tnuc                Temperature of nucleation (K)
!> @param[in,out] superc_liq_wth      Supercooled cloud liquid water content
!> @param[in,out] superc_rain_wth     Supercooled rain water content
!> @param[in,out] sfwater             Sub-grid orographic water produced by Seeder Feeder scheme
!> @param[in,out] sfrain              Extra rain produced by Seeder Feeder scheme
!> @param[in,out] sfsnow              Extra snow produced by Seeder Feeder scheme
!> @param[in,out] refl_tot            Total (3d) radar reflectivity (dBZ)
!> @param[in,out] refl_1km            Radar reflectivity at 1km above the surface
!> @param[in]     ndf_wth             Number of degrees of freedom per cell for
!!                                     potential temperature space
!> @param[in]     undf_wth            Number unique of degrees of freedom in segment for
!!                                     potential temperature space
!> @param[in]     map_wth             Dofmap for a segment of any potential temperature space field
!> @param[in]     ndf_w3              Number of degrees of freedom per cell for
!!                                     density space
!> @param[in]     undf_w3             Number unique of degrees of freedom in segment for
!!                                     density space
!> @param[in]     map_w3              Dofmap for a segment of any density space field
!> @param[in]     ndf_2d              Number of degrees of freedom per cell for
!!                                     2D fields
!> @param[in]     undf_2d             Number unique of degrees of freedom in segment for
!!                                     2D fields
!> @param[in]     map_2d              Dofmap for a segment of any 2D field
!> @param[in]     ndf_farr            Number of degrees of freedom per cell for fsd array
!> @param[in]     undf_farr           Number unique of degrees of freedom for fsd array
!> @param[in]     map_farr            Dofmap for a segment of the fsd array field

subroutine mphys_code( nlayers, seg_len,            &
                       mv_wth,   ml_wth,   ms_wth,  &
                       mr_wth,   mg_wth,            &
                       cf_wth,   cfl_wth,  cff_wth, &
                       u_in_w3, v_in_w3,            &
                       theta_in_wth,                &
                       exner_in_wth, wetrho_in_w3,  &
                       dry_rho_in_w3,               &
                       height_w3, height_wth,       &
                       cloud_drop_no_conc, murk,    &
                       sd_orog,                     &
                       dmv_wth,  dml_wth,  dms_wth, &
                       dmr_wth,  dmg_wth,           &
                       ls_rain_2d, ls_snow_2d,      &
                       ls_graup_2d, lsca_2d,        &
                       precfrac,                    &
                       ls_rain_3d, ls_snow_3d,      &
                       autoconv, accretion,         &
                       rim_cry, rim_agg,            &
                       dtheta,                      &
                       dcfl_wth, dcff_wth, dbcf_wth,&
                       f_arr_wth, tnuc,             &
                       superc_liq_wth,              &
                       superc_rain_wth,             &
                       sfwater, sfrain, sfsnow,     &
                       refl_tot, refl_1km,          &
                       ndf_wth, undf_wth, map_wth,  &
                       ndf_w3,  undf_w3,  map_w3,   &
                       ndf_2d,  undf_2d,  map_2d,   &
                       ndf_farr,undf_farr,map_farr  )

    use constants_mod,              only: r_def, i_def, r_um, i_um
    use cloud_config_mod,           only: cld_fsd_hill

    !---------------------------------------
    ! UM modules
    !---------------------------------------

    use mphys_inputs_mod,           only: l_mcr_qrain,                         &
                                          l_mcr_qgraup, l_mcr_precfrac,        &
                                          l_subgrid_graupel_frac
    use mphys_constants_mod,        only: mprog_min
    use mphys_bypass_mod,           only: l_ref_diag
    use mphys_radar_mod,            only: ref_lim

    use cloud_inputs_mod,           only: i_cld_vn, rhcrit
    use pc2_constants_mod,          only: i_cld_off, i_cld_pc2

    use ls_ppn_mod,                 only: ls_ppn

    use planet_constants_mod,       only: p_zero, kappa, planet_radius
    use water_constants_mod,        only: tm
    use arcl_mod,                   only: npd_arcl_compnts
    use def_easyaerosol,            only: t_easyaerosol_cdnc


    use mphys_diags_mod,            only: praut, pracw, piacw, psacw,          &
                                          dbz_tot, dbz_g, dbz_i, dbz_i2,       &
                                          dbz_r, dbz_l,                        &
                                          sfwater_um => sfwater,               &
                                          sfrain_um => sfrain,                 &
                                          sfsnow_um => sfsnow,                 &
                                          l_praut_diag, l_pracw_diag,          &
                                          l_piacw_diag, l_psacw_diag,          &
                                          l_refl_tot, l_refl_1km,              &
                                          l_sfwater_diag, l_sfrain_diag,       &
                                          l_sfsnow_diag

    use microphysics_config_mod,    only: orog_rain, orog_block, nsigmasf

    use lsp_froude_moist_mod,       only: lsp_froude_moist

    use free_tracers_inputs_mod,    only: n_wtrac
    use wtrac_atm_step_mod,         only: atm_step_wtrac_type
    use wtrac_mphys_mod,            only: mp_wtrac_type

    implicit none

    ! Arguments

    integer(kind=i_def), intent(in) :: nlayers, seg_len
    integer(kind=i_def), intent(in) :: ndf_wth,  ndf_w3,  ndf_2d,  ndf_farr
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3, undf_2d, undf_farr

    real(kind=r_def), intent(in),  dimension(undf_wth) :: mv_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: ml_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: ms_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: mr_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: mg_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cf_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cfl_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cff_wth
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: u_in_w3
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: v_in_w3
    real(kind=r_def), intent(in),  dimension(undf_wth) :: theta_in_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: exner_in_wth
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: wetrho_in_w3
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: dry_rho_in_w3
    real(kind=r_def), intent(in),  dimension(undf_w3)  :: height_w3
    real(kind=r_def), intent(in),  dimension(undf_wth) :: height_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: cloud_drop_no_conc
    real(kind=r_def), intent(inout), dimension(undf_wth) :: murk
    real(kind=r_def), intent(in),  dimension(undf_2d)  :: sd_orog
    real(kind=r_def), intent(in),  dimension(undf_farr):: f_arr_wth
    real(kind=r_def), intent(in),  dimension(undf_wth) :: tnuc

    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmv_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dml_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dms_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmr_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dmg_wth
    real(kind=r_def), intent(inout), dimension(undf_2d)  :: ls_rain_2d
    real(kind=r_def), intent(inout), dimension(undf_2d)  :: ls_snow_2d
    real(kind=r_def), intent(inout), dimension(undf_2d)  :: ls_graup_2d
    real(kind=r_def), intent(inout), dimension(undf_2d)  :: lsca_2d
    real(kind=r_def), intent(inout), dimension(undf_wth) :: precfrac
    real(kind=r_def), intent(inout), dimension(undf_wth) :: ls_rain_3d
    real(kind=r_def), intent(inout), dimension(undf_wth) :: ls_snow_3d
    real(kind=r_def), intent(inout), dimension(undf_wth) :: autoconv
    real(kind=r_def), intent(inout), dimension(undf_wth) :: accretion
    real(kind=r_def), intent(inout), dimension(undf_wth) :: rim_cry
    real(kind=r_def), intent(inout), dimension(undf_wth) :: rim_agg
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dtheta
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dcfl_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dcff_wth
    real(kind=r_def), intent(inout), dimension(undf_wth) :: dbcf_wth

    real(kind=r_def), pointer, intent(inout) :: superc_liq_wth(:)
    real(kind=r_def), pointer, intent(inout) :: superc_rain_wth(:)

    real(kind=r_def), pointer, intent(inout) :: sfwater(:)
    real(kind=r_def), pointer, intent(inout) :: sfrain(:)
    real(kind=r_def), pointer, intent(inout) :: sfsnow(:)
    real(kind=r_def), pointer, intent(inout) :: refl_tot(:)
    real(kind=r_def), pointer, intent(inout) :: refl_1km(:)

    integer(kind=i_def), intent(in), dimension(ndf_wth, seg_len) :: map_wth
    integer(kind=i_def), intent(in), dimension(ndf_w3, seg_len)  :: map_w3
    integer(kind=i_def), intent(in), dimension(ndf_2d, seg_len)  :: map_2d
    integer(kind=i_def), intent(in), dimension(ndf_farr, seg_len):: map_farr

    ! Local variables for the kernel

    real(r_um), dimension(seg_len,1,nlayers) ::                                &
         u_on_p, v_on_p, q_work, qcl_work, qcf_work, deltaz, cfl_work,         &
         cff_work, cf_work, rhodz_dry, rhodz_moist, t_n, t_work,               &
         p_theta_levels, ls_rain3d, ls_snow3d, ls_graup3d, rainfrac3d,         &
         n_drop_pot, n_drop_3d, so4_accu_work, so4_diss_work,                  &
         aged_bmass_work, cloud_bmass_work, aged_ocff_work, cloud_ocff_work,   &
         nitr_acc_work, nitr_diss_work, aerosol_work, biogenic, rho_r2,        &
         dry_rho, ukca_cdnc_array, tnuc_new, theta, z_rho, z_theta,            &
         r_rho_levels

    real(r_um), dimension(seg_len, 1, nlayers, 1) :: arcl

    real(r_um), dimension(seg_len, 1, 0:nlayers) :: flash_pot, r_theta_levels

    real(r_um), dimension(seg_len, 1) :: ls_rain, ls_snow, ls_graup,           &
                                         snow_depth, land_frac, hmteff, zb,    &
                                         cos_theta_latitude

    real(r_um), dimension(:,:,:), allocatable :: qrain_work, qcf2_work,        &
                                                 qgraup_work, precfrac_work

    real(r_um), dimension(nlayers) :: rhcpt
    real(r_um), dimension(3, seg_len, 1, nlayers) :: f_arr

    real(r_um), dimension(1,1,1) :: sea_salt_film, sea_salt_jet

    real(r_um), parameter :: alt_1km = 1000.0_r_um ! metres

    logical, dimension(seg_len,1) :: land_sea_mask

    integer(i_um) :: i,k,n
    integer(i_um), parameter :: j=1

    integer(i_um), dimension(npd_arcl_compnts) :: i_arcl_compnts

    integer(i_um), dimension(seg_len) :: land_index

    real(r_um), dimension(seg_len) :: ls_rainfrac

    integer(i_um) :: lspice_dim1, lspice_dim2, lspice_dim3,                    &
                     salt_dim1, salt_dim2, salt_dim3, cdnc_dim1, cdnc_dim2,    &
                     cdnc_dim3, rhc_row_length, rhc_rows,                      &
                     n_arcl_compnts, land_points

    logical :: l_cosp_lsp

    type (t_easyaerosol_cdnc) :: easyaerosol_cdnc

    ! Water tracer fields which are not currently used but are required by
    ! UM routine
    type (atm_step_wtrac_type), dimension(n_wtrac) :: wtrac_as
    type (mp_wtrac_type), dimension(n_wtrac) :: wtrac_mp

    !-----------------------------------------------------------------------
    ! Initialisation of non-prognostic variables and arrays
    !-----------------------------------------------------------------------

    ! These must be set as below to match the declarations above
    lspice_dim1 = int(seg_len, i_um)
    lspice_dim2 = 1_i_um
    lspice_dim3 = nlayers

    allocate ( easyaerosol_cdnc % cdnc(1,1,1) )
    easyaerosol_cdnc % dim1 = 1_i_um
    easyaerosol_cdnc % dim2 = 1_i_um
    easyaerosol_cdnc % dim3 = 1_i_um

    cdnc_dim1   = int(seg_len, i_um)
    cdnc_dim2   = 1_i_um
    cdnc_dim3   = nlayers
    do i = 1, seg_len
      do k = 1, nlayers
        ukca_cdnc_array(i,j,k) = cloud_drop_no_conc(map_wth(1,i) + k)
      end do
    end do

    salt_dim1   = 1_i_um
    salt_dim2   = 1_i_um
    salt_dim3   = 1_i_um

    rhc_row_length = 1_i_um
    rhc_rows       = 1_i_um
    n_arcl_compnts = 1_i_um

    land_points = int(seg_len, i_um)
    do i = 1, seg_len
      land_index(i)  = int(i, i_um)
    end do

    if (murk_prognostic) then
      do i = 1, seg_len
        do k = 1, nlayers
          aerosol_work(i,j,k) = murk(map_wth(1,i) + k)
        end do
      end do
    end if

    land_sea_mask = .false.

    l_cosp_lsp = .false.

    if (prog_tnuc) then
      do i = 1, seg_len
        do k = 1, nlayers
          tnuc_new(i,j,k) = real(tnuc(map_wth(1,i) + k),kind=r_um)
        end do ! k
      end do ! i
    end if

    !-----------------------------------------------------------------------
    ! Initialisation of prognostic variables and arrays
    !-----------------------------------------------------------------------

    ! This assumes that map_wth(1) points to level 0
    ! and map_w3(1) points to level 1
    do i = 1, seg_len
      do k = 1, nlayers
        ! height of levels from centre of planet
        r_rho_levels(i,j,k)   = height_w3(map_w3(1,i) + k-1) + planet_radius

        rho_r2(i,j,k) = wetrho_in_w3(map_w3(1,i) + k-1) *                      &
                        ( r_rho_levels(i,j,k)**2 )
        dry_rho(i,j,k) = dry_rho_in_w3(map_w3(1,i) + k-1)

        u_on_p(i,j,k) = u_in_w3(map_w3(1,i) + k-1)
        v_on_p(i,j,k) = v_in_w3(map_w3(1,i) + k-1)
        t_n(i,j,k)    = theta_in_wth(map_wth(1,i) + k) *                       &
                        exner_in_wth(map_wth(1,i) + k)
        ! N.B. dtheta is actually a temperature increment when passed in
        t_work(i,j,k) = t_n(i,j,k) + dtheta(map_wth(1,i) + k)

        ! pressure on theta levels
        p_theta_levels(i,j,k)    = p_zero*(exner_in_wth(map_wth(1,i) + k))     &
                                        **(1.0_r_um/kappa)
        ! Compulsory moist prognostics
        q_work(i,j,k)    = mv_wth(map_wth(1,i) + k) + dmv_wth(map_wth(1,i) + k )
        qcl_work(i,j,k)  = ml_wth(map_wth(1,i) + k) + dml_wth(map_wth(1,i) + k )
        qcf_work(i,j,k)  = ms_wth(map_wth(1,i) + k) + dms_wth(map_wth(1,i) + k )

      end do ! k
    end do ! i

    do i = 1, seg_len
      do k = 0, nlayers
        r_theta_levels(i,j,k) = height_wth(map_wth(1,i) + k) + planet_radius
      end do
    end do ! i

    ! Optional moist prognostics

    ! Perform allocation of the qcf2 variable as it is required in the UM
    ! microphysics, even if it is not actually used.
    allocate(qcf2_work(1,1,1))

    if (l_mcr_qrain) then
      allocate (qrain_work (seg_len, 1, nlayers) )
      do i = 1, seg_len
        do k = 1, nlayers
          qrain_work(i,j,k) = mr_wth(map_wth(1,i) + k)
        end do ! k
      end do ! i
    else
      allocate (qrain_work(1,1,1))
    end if

    if (l_mcr_qgraup) then
      allocate (qgraup_work (seg_len, 1, nlayers) )
      do i = 1, seg_len
        do k = 1, nlayers
          qgraup_work(i,j,k) = mg_wth(map_wth(1,i) + k)
        end do ! k
      end do ! i
    else
      allocate(qgraup_work(1,1,1))
    end if

    if ( l_mcr_precfrac ) then
      ! Prognostic precipitation fraction...
      ! set from LFRic input prognostic
      allocate (precfrac_work (seg_len, 1, nlayers) )
      do i = 1, seg_len
        do k = 1, nlayers
          precfrac_work(i,j,k) = precfrac(map_wth(1,i)+k)
        end do ! k
      end do ! i

    else  ! ( l_mcr_precfrac )
      ! Prognostic precipitation fraction switched off; minimal allocation

      allocate(precfrac_work(1,1,1))

    end if  ! ( l_mcr_precfrac )

    if ( cld_fsd_hill ) then
      ! Parameters from fractional standard deviation (FSD) parametrization
      ! There are 3 parameters used in the empirical fit, each stored as a different
      ! element in the f_arr array.
      do i = 1, seg_len
        do n = 1, 3
          do k = 1, nlayers
            f_arr(n,i,j,k) = f_arr_wth(map_farr(1,i) + (n-1)*(nlayers+1) + k)
          end do
        end do
      end do
    end if

    ! Set this to 1 to account for quasi-uniform grid
    do i = 1, seg_len
      cos_theta_latitude(i,1) = 1.0_r_um
    end do

    ! Note: need other options once Smith scheme is in use.
    if ( i_cld_vn == i_cld_off ) then
      do k = 1, nlayers
        do i = 1, seg_len
          if (qcl_work(i,j,k) >= mprog_min) then
            cfl_work(i,j,k) = 1.0_r_um
          else
            cfl_work(i,j,k) = 0.0_r_um
          end if

          if (qcf_work(i,j,k) >= mprog_min) then
            cff_work(i,j,k) = 1.0_r_um
          else
            cff_work(i,j,k) = 0.0_r_um
          end if

          cf_work(i,j,k) = max( cff_work(i,j,k), cfl_work(i,j,k) )

        end do

        rhcpt(k) = 1.0_r_um

      end do

    else ! i_cld_vn > 0
      do k = 1, nlayers
        do i = 1, seg_len
          cf_work(i,j,k)  = cf_wth( map_wth(1,i) + k) + dbcf_wth( map_wth(1,i) + k)
          cfl_work(i,j,k) = cfl_wth(map_wth(1,i) + k) + dcfl_wth( map_wth(1,i) + k)
          cff_work(i,j,k) = cff_wth(map_wth(1,i) + k) + dcff_wth( map_wth(1,i) + k)
        end do

        rhcpt(k) = rhcrit(k)

      end do
    end if ! i_cld_vn

    l_sfwater_diag = .not. associated(sfwater, empty_real_data)
    l_sfrain_diag = .not. associated(sfrain, empty_real_data)
    l_sfsnow_diag = .not. associated(sfsnow, empty_real_data)

    l_refl_tot = .not. associated(refl_tot, empty_real_data)
    l_refl_1km = .not. associated(refl_1km, empty_real_data)

    ! Calculate low-level blocking for the  Seeder Feeder scheme
    if (orog_rain) then

      do i = 1, seg_len
        zb(i,j) = 0.0_r_um
        hmteff(i,j) = sd_orog(map_2d(1,i)) * nsigmasf
      end do

      if (orog_block) then

        ! Theta and level height above the surface
        do i = 1, seg_len
          do k = 1, nlayers
            theta(i,j,k)    = theta_in_wth(map_wth(1,i) + k)
            z_rho(i,j,k)    = r_rho_levels(i,j,k) - r_theta_levels(i,j,0)
            z_theta(i,j,k)  = r_theta_levels(i,j,k) - r_theta_levels(i,j,0)
          end do ! i
        end do ! k

        CALL lsp_froude_moist( nlayers, seg_len,                              &
                          u_on_p, v_on_p, theta,                              &
                          z_rho, z_theta,                                     &
                          hmteff, zb )

      endif

    else
      ! Needed so conversion to single precision works
      zb = 0.0_r_um
      hmteff = 0.0_r_um
    endif

    ! Allocate arrays for diagnostics
    if (l_praut_diag) then
      allocate(praut( seg_len, 1, nlayers ))
      do k = 1, nlayers
        do i = 1, seg_len
          praut(i,j,k) = 0.0_r_um
        end do
      end do
    endif
    if (l_pracw_diag) then
      allocate(pracw( seg_len, 1, nlayers ))
      do k = 1, nlayers
        do i = 1, seg_len
          pracw(i,j,k) = 0.0_r_um
        end do
      end do
    endif
    if (l_piacw_diag) then
      allocate(piacw( seg_len, 1, nlayers ))
      do k = 1, nlayers
        do i = 1, seg_len
          piacw(i,j,k) = 0.0_r_um
        end do
      end do
    endif
    if (l_psacw_diag) then
      allocate(psacw( seg_len, 1, nlayers ))
      do k = 1, nlayers
        do i = 1, seg_len
          psacw(i,j,k) = 0.0_r_um
        end do
      end do
    endif
    if (l_sfwater_diag) then
      allocate(sfwater_um( seg_len, 1, nlayers ))
    endif
    if (l_sfrain_diag) then
      allocate(sfrain_um( seg_len, 1, nlayers ))
    endif
    if (l_sfsnow_diag) then
      allocate(sfsnow_um( seg_len, 1, nlayers ))
    endif

    if (l_refl_tot .or. l_refl_1km) then

      ! Allocate space for radar reflectivity fields. All radar fields
      ! currently need to be allocated in order for the UM routines to not
      ! produce segmentation faults.
      allocate(dbz_tot( seg_len, 1, nlayers ))
      allocate(dbz_g(   seg_len, 1, nlayers ))
      allocate(dbz_r(   seg_len, 1, nlayers ))
      allocate(dbz_i(   seg_len, 1, nlayers ))
      allocate(dbz_i2(  seg_len, 1, nlayers ))
      allocate(dbz_l(   seg_len, 1, nlayers ))

      ! Initialise fields to radar reflectivity limit (not zero as these
      ! fields are in log space and can go negative).
      dbz_tot = ref_lim
      dbz_g = ref_lim
      dbz_r = ref_lim
      dbz_i = ref_lim
      dbz_i2 = ref_lim
      dbz_l = ref_lim

      ! The following flag tells the UM microphysics to call the radar
      ! reflectivity generating routine.
      l_ref_diag = .true.

    end if ! l_refl_tot or l_refl_1km

    ! CALL to ls_ppn
    call ls_ppn(                                                               &
                p_theta_levels,                                                &
                land_sea_mask, deltaz, r_theta_levels, r_rho_levels,           &
                cf_work, cfl_work, cff_work, precfrac_work,                    &
                rhcpt, f_arr, cos_theta_latitude,                              &
                lspice_dim1,lspice_dim2,lspice_dim3,                           &
                rho_r2, dry_rho, q_work, qcf_work, qcl_work, t_work,           &
                qcf2_work, qrain_work, qgraup_work,                            &
                u_on_p, v_on_p,                                                &
                sea_salt_film, sea_salt_jet,                                   &
                salt_dim1, salt_dim2, salt_dim3,                               &
                ukca_cdnc_array,                                               &
                cdnc_dim1, cdnc_dim2, cdnc_dim3,                               &
                easyaerosol_cdnc,                                              &
                biogenic,                                                      &
                snow_depth, land_frac,                                         &
                so4_accu_work,                                                 &
                so4_diss_work, aged_bmass_work, cloud_bmass_work,              &
                aged_ocff_work, cloud_ocff_work, nitr_acc_work,                &
                nitr_diss_work, aerosol_work,                                  &
                n_arcl_compnts, i_arcl_compnts, arcl,                          &
                ls_rain, ls_snow, ls_graup,                                    &
                ls_rain3d, ls_snow3d, ls_graup3d, rainfrac3d,                  &
                n_drop_pot, n_drop_3d,                                         &
                rhc_row_length, rhc_rows,                                      &
                rhodz_dry, rhodz_moist,                                        &
                ls_rainfrac, land_points, land_index,                          &
                l_cosp_lsp,                                                    &
                hmteff, zb, tnuc_new,                                          &
                wtrac_as, wtrac_mp)

  ! Update theta and compulsory prognostic variables
  do i = 1, seg_len
    do k = 1, nlayers
      dtheta(map_wth(1,i) + k) = ( t_work(i,j,k) - t_n(i,j,k) )             &
                                / exner_in_wth(map_wth(1,i) + k)

      dmv_wth(map_wth(1,i) + k ) = q_work(i,j,k)   - mv_wth( map_wth(1,i) + k )
      dml_wth(map_wth(1,i) + k ) = qcl_work(i,j,k) - ml_wth( map_wth(1,i) + k )
      dms_wth(map_wth(1,i) + k ) = qcf_work(i,j,k) - ms_wth( map_wth(1,i) + k )
    end do
  end do ! k (nlayers)

  ! Increment level 0 the same as level 1
  !  (as done in the UM)
  do i = 1, seg_len
    dtheta(map_wth(1,i) + 0)    = dtheta(map_wth(1,i) + 1)
    dmv_wth(map_wth(1,i) + 0)   = dmv_wth(map_wth(1,i) + 1)
    dml_wth(map_wth(1,i) + 0)   = dml_wth(map_wth(1,i) + 1)
    dms_wth(map_wth(1,i) + 0)   = dms_wth(map_wth(1,i) + 1)
  end do

  ! Update optional additional prognostic variables
  ! No need for else statements here as dms_wth and associated variables
  ! should have already been initialised to zero.

  if (l_mcr_qrain) then
    do i = 1, seg_len
      ! Update level 0 to be the same as level 1 (as per UM)
      dmr_wth(map_wth(1,i) + 0) = qrain_work(i,j,1) - mr_wth(map_wth(1,i) + 0)
      do k = 1, nlayers
        dmr_wth( map_wth(1,i) + k) = qrain_work(i,j,k) - mr_wth( map_wth(1,i) + k )
      end do
    end do
  end if

  if (l_mcr_qgraup) then
    do i = 1, seg_len
      ! Update level 0 to be the same as level 1 (as per UM)
      dmg_wth(map_wth(1,i) + 0) = qgraup_work(i,j,1) - mg_wth(map_wth(1,i) + 0)
      do k = 1, nlayers
        dmg_wth( map_wth(1,i) + k) = qgraup_work(i,j,k) - mg_wth( map_wth(1,i) + k )
      end do
    end do
  end if

  if (murk_prognostic) then
    do k = 1, nlayers
      do i = 1, seg_len
        murk(map_wth(1,i) + k) = aerosol_work(i,j,k)
      end do
    end do
    do i = 1, seg_len
      murk(map_wth(1,i) + 0) = murk(map_wth(1,i) + 1)
    end do
  end if

  ! Cloud fraction increments
  ! Always calculate them, but only add them on in slow_physics if using PC2.
  do i = 1, seg_len
    do k = 1, nlayers
      dbcf_wth( map_wth(1,i) + k) = cf_work(i,j,k)  - cf_wth(  map_wth(1,i) + k )
      dcfl_wth( map_wth(1,i) + k) = cfl_work(i,j,k) - cfl_wth( map_wth(1,i) + k )
      dcff_wth( map_wth(1,i) + k) = cff_work(i,j,k) - cff_wth( map_wth(1,i) + k )
    end do
    ! Set level 0 as the same as level 1
    dbcf_wth(map_wth(1,i) + 0) = dbcf_wth(map_wth(1,i) + 1)
    dcfl_wth(map_wth(1,i) + 0) = dcfl_wth(map_wth(1,i) + 1)
    dcff_wth(map_wth(1,i) + 0) = dcff_wth(map_wth(1,i) + 1)

    ! Copy ls_rain and ls_snow
    ls_rain_2d(map_2d(1,i))  = ls_rain(i,j)
    ls_snow_2d(map_2d(1,i))  = ls_snow(i,j)
    ls_graup_2d(map_2d(1,i)) = ls_graup(i,j)

    lsca_2d(map_2d(1,i))     = ls_rainfrac(i)
    do k = 1, nlayers
      ls_rain_3d(map_wth(1,i) + k) = ls_rain3d(i,j,k)
      ls_snow_3d(map_wth(1,i) + k) = ls_snow3d(i,j,k)
    end do ! model levels
  end do

  if ( l_mcr_precfrac ) then
    do k = 1, nlayers
      do i = 1, seg_len
        precfrac(map_wth(1,i)+k) = precfrac_work(i,j,k)
      end do
    end do ! k
    do i = 1, seg_len
      precfrac(map_wth(1,i)+0) = precfrac(map_wth(1,i)+1)
    end do
  end if

  ! Copy diagnostics if selected: autoconversion, accretion, riming rates &
  ! radar reflectivity.
  if (l_praut_diag) then
    do k = 1, nlayers
      do i = 1, seg_len
        autoconv(map_wth(1,i) + k) = praut(i,j,k)
      end do
    end do
  end if
  if (l_pracw_diag) then
    do k = 1, nlayers
      do i = 1, seg_len
        accretion(map_wth(1,i) + k) = pracw(i,j,k)
      end do
    end do
  end if
  if (l_piacw_diag) then
    do k = 1, nlayers
      do i = 1, seg_len
        rim_cry(map_wth(1,i) + k) = piacw(i,j,k)
      end do
    end do
  end if
  if (l_psacw_diag) then
    do k = 1, nlayers
      do i = 1, seg_len
        rim_agg(map_wth(1,i) + k) = psacw(i,j,k)
      end do
    end do
  end if
  if (l_refl_tot) then
    do i = 1, seg_len
      refl_tot(map_wth(1,i)) = ref_lim ! 0 level
    end do
    do k = 1, nlayers
      do i = 1, seg_len
        refl_tot(map_wth(1,i) + k) = dbz_tot(i,j,k)
      end do
    end do
  end if

  if (l_refl_1km) then
    do i = 1, seg_len
      do k = 1, nlayers
        ! Select the first altitude above 1km (following what the UM does).
        if (height_wth(map_wth(1,i) + k) >= alt_1km ) then
          refl_1km(map_2d(1,i)) = dbz_tot(i,j,k)
          exit
        end if
      end do
    end do
  end if

  if (.not. associated(superc_liq_wth, empty_real_data) ) then
    do i = 1, seg_len
      do k = 1, nlayers
        if (t_n(i,j,k) < tm) then
          ! Following the UM, have taken start of timestep quantity for the
          ! supercooled liquid cloud. This is where the model cloud should
          ! be in a steady-state.
          superc_liq_wth( map_wth(1,i) + k) = ml_wth(map_wth(1,i) + k)
        else
          superc_liq_wth( map_wth(1,i) + k) = 0.0_r_um
        end if
      end do ! nlayers
    end do ! seg_len
  end if ! not assoc. superc_liq_wth

  if (.not. associated(superc_rain_wth, empty_real_data) ) then
    do i = 1, seg_len
      do k = 1, nlayers
        if (t_n(i,j,k) < tm) then
          ! Following the UM, have taken start of timestep quantity for the
          ! supercooled rain. This is where the model cloud should
          ! be in a steady-state.
          superc_rain_wth( map_wth(1,i) + k) = mr_wth(map_wth(1,i) + k)
        else
          superc_rain_wth( map_wth(1,i) + k) = 0.0_r_um
        end if
      end do ! nlayers
    end do ! seg_len
  end if ! not assoc. superc_rain_wth

  if (l_sfwater_diag) then
    do k=1,nlayers
      do i=1,seg_len
        sfwater(map_wth(1,i)+k) = sfwater_um(i,1,k)
      end do
    end do
  end if
  if (l_sfrain_diag) then
    do k=1,nlayers
      do i=1,seg_len
        sfrain(map_wth(1,i)+k) = sfrain_um(i,1,k)
      end do
    end do
  end if
  if (l_sfsnow_diag) then
    do k=1,nlayers
      do i=1,seg_len
        sfsnow(map_wth(1,i)+k) = sfsnow_um(i,1,k)
      end do
    end do
  end if

  if (allocated(psacw)) deallocate (psacw)
  if (allocated(piacw)) deallocate (piacw)
  if (allocated(pracw)) deallocate (pracw)
  if (allocated(praut)) deallocate (praut)
  if (allocated(sfwater_um)) deallocate (sfwater_um)
  if (allocated(sfrain_um)) deallocate (sfrain_um)
  if (allocated(sfsnow_um)) deallocate (sfsnow_um)
  if (allocated(dbz_tot)) deallocate(dbz_tot)
  if (allocated(dbz_g)) deallocate(dbz_g)
  if (allocated(dbz_r)) deallocate(dbz_r)
  if (allocated(dbz_i)) deallocate(dbz_i)
  if (allocated(dbz_i2)) deallocate(dbz_i2)
  if (allocated(dbz_l)) deallocate(dbz_l)

  deallocate( precfrac_work )
  deallocate( qgraup_work )
  deallocate( qrain_work  )
  deallocate( qcf2_work   )
  deallocate( easyaerosol_cdnc % cdnc )

  ! Reset the radar reflectivity call flag so the fields are
  ! only generated when needed.
  l_ref_diag = .false.

  ! N.B. Calls to aerosol code (rainout, mass_calc) etc have been omitted
  ! as it is expected these will be retired for LFRic/GHASP.
  ! Call to diagnostics also omitted here, as it will probably have a different
  ! structure under LFRic.

end subroutine mphys_code

end module mphys_kernel_mod
