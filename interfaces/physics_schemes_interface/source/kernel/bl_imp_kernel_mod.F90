!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the implicit UM boundary layer scheme.
!>
module bl_imp_kernel_mod

  use argument_mod,           only : arg_type,                  &
                                     GH_FIELD, GH_SCALAR,       &
                                     GH_INTEGER, GH_REAL,       &
                                     GH_READ, GH_READWRITE,     &
                                     GH_WRITE, DOMAIN,          &
                                     ANY_DISCONTINUOUS_SPACE_1, &
                                     ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,             only : i_def, i_um, r_def, r_um, r_bl
  use fs_continuity_mod,         only : W3, Wtheta
  use kernel_mod,                only : kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: bl_imp_kernel_type
    private
    type(arg_type) :: meta_args(30) = (/                                          &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                &! loop
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! theta_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! exner_in_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_v_n
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_cl_n
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_cf_n
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! theta_star
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! height_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! height_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! dt_conv
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_v
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_cl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! m_cf
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   &! dtrdz_tq_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! rdz_tq_bl
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),&! blend_height_tq
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2),&! bl_type_ind
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       &! rhokh_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! moist_flux_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, W3),                       &! heat_flux_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTheta),                   &! dqw_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, WTheta),                   &! dtl_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     WTheta),                   &! dqw_nt_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     WTheta),                   &! dtl_nt_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     WTheta),                   &! qw_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     WTheta),                   &! tl_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     WTheta),                   &! ct_ctq_wth
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! dqw1_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),&! dtl1_2d
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1) &! ct_ctq1_2d
         /)

    integer :: operates_on = DOMAIN
  contains
    procedure, nopass :: bl_imp_code
  end type

  public :: bl_imp_code

contains

  !> @brief Interface to the implicit UM BL scheme
  !> @details The UM Boundary Layer scheme does:
  !>             vertical mixing of heat, momentum and moisture,
  !>             as documented in UMDP24
  !> @param[in]     nlayers              Number of layers
  !> @param[in]     loop                 BL predictor-corrector counter
  !> @param[in]     theta_in_wth         Potential temperature field
  !> @param[in]     exner_in_wth         Exner pressure field in wth space
  !> @param[in]     m_v_n                Vapour mixing ratio at time level n
  !> @param[in]     m_cl_n               Cloud liq mixing ratio at time level n
  !> @param[in]     m_cf_n               Cloud fro mixing ratio at time level n
  !> @param[in]     theta_star           Potential temperature after advection
  !> @param[in]     height_w3            Height of density space above surface
  !> @param[in]     height_wth           Height of theta space above surface
  !> @param[in]     dt_conv              Convection temperature increment
  !> @param[in]     m_v                  Vapour mixing ration after advection
  !> @param[in]     m_cl                 Cloud liq mixing ratio after advection
  !> @param[in]     m_cf                 Cloud fro mixing ratio after advection
  !> @param[in]     dtrdz_tq_bl          dt/(rho*r*r*dz) in wth
  !> @param[in]     rdz_tq_bl            1/dz in w3
  !> @param[in]     blend_height_tq      Blending height for wth levels
  !> @param[in]     bl_type_ind          Boundary layer type indicator
  !> @param[in]     rhokh_bl             Heat eddy diffusivity on BL levels
  !> @param[in,out] moist_flux_bl        Vertical moisture flux on BL levels
  !> @param[in,out] heat_flux_bl         Vertical heat flux on BL levels
  !> @param[in,out] dqw_wth              Total increment to total water
  !> @param[in,out] dtl_wth              Total increment to liquid temperature
  !> @param[in,out] dqw_nt_wth           Non-turbulent increment to total water
  !> @param[in,out] dtl_nt_wth           Non-turbulent increment to liquid temperature
  !> @param[in,out] qw_wth               Total water
  !> @param[in,out] tl_wth               Liquid temperature
  !> @param[in,out] ct_ctq_wth           Coefficient in predictor-corrector
  !> @param[in,out] dqw1_2d              Total water increment at blending height
  !> @param[in,out] dtl1_2d              Liquid temperature increment at blending height
  !> @param[in,out] ct_ctq1_2d           Coefficient at blending height
  !> @param[in]     ndf_wth              Number of DOFs per cell for potential temperature space
  !> @param[in]     undf_wth             Number of unique DOFs for potential temperature space
  !> @param[in]     map_wth              Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_w3               Number of DOFs per cell for density space
  !> @param[in]     undf_w3              Number of unique DOFs for density space
  !> @param[in]     map_w3               Dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_2d               Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d              Number of unique DOFs for 2D fields
  !> @param[in]     map_2d               Dofmap for the cell at the base of the column for 2D fields
  !> @param[in]     ndf_bl               Number of DOFs per cell for BL types
  !> @param[in]     undf_bl              Number of total DOFs for BL types
  !> @param[in]     map_bl               Dofmap for cell for BL types
  subroutine bl_imp_code(nlayers, seg_len,                   &
                         loop,                               &
                         theta_in_wth,                       &
                         exner_in_wth,                       &
                         m_v_n,                              &
                         m_cl_n,                             &
                         m_cf_n,                             &
                         theta_star,                         &
                         height_w3,                          &
                         height_wth,                         &
                         dt_conv,                            &
                         m_v,                                &
                         m_cl,                               &
                         m_cf,                               &
                         dtrdz_tq_bl,                        &
                         rdz_tq_bl,                          &
                         blend_height_tq,                    &
                         bl_type_ind,                        &
                         rhokh_bl,                           &
                         moist_flux_bl,                      &
                         heat_flux_bl,                       &
                         dqw_wth, dtl_wth, dqw_nt_wth,       &
                         dtl_nt_wth, qw_wth, tl_wth,         &
                         ct_ctq_wth, dqw1_2d, dtl1_2d,       &
                         ct_ctq1_2d,                         &
                         ndf_wth,                            &
                         undf_wth,                           &
                         map_wth,                            &
                         ndf_w3,                             &
                         undf_w3,                            &
                         map_w3,                             &
                         ndf_2d,                             &
                         undf_2d,                            &
                         map_2d,                             &
                         ndf_bl, undf_bl, map_bl)

    !---------------------------------------
    ! UM modules containing switches or global constants
    !---------------------------------------
    use bl_option_mod, only: l_noice_in_turb, alpha_cd, puns, pstb, &
         flux_bc_opt, specified_fluxes_only
    use nlsizes_namelist_mod, only: bl_levels
    use planet_constants_mod, only: planet_radius

    ! subroutines used
    use bdy_impl3_mod, only: bdy_impl3

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, seg_len
    integer(kind=i_def), intent(in) :: loop

    integer(kind=i_def), intent(in) :: ndf_wth, ndf_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: undf_wth, undf_w3
    integer(kind=i_def), intent(in) :: map_wth(ndf_wth,seg_len)
    integer(kind=i_def), intent(in) :: map_w3(ndf_w3,seg_len)
    integer(kind=i_def), intent(in) :: map_2d(ndf_2d,seg_len)
    integer(kind=i_def), intent(in) :: ndf_bl, undf_bl
    integer(kind=i_def), intent(in) :: map_bl(ndf_bl,seg_len)

    real(kind=r_def), dimension(undf_w3),  intent(inout) :: moist_flux_bl,     &
                                                            heat_flux_bl
    real(kind=r_def), dimension(undf_wth), intent(inout) :: dqw_wth, dtl_wth,  &
                                                            dqw_nt_wth,        &
                                                            dtl_nt_wth, qw_wth,&
                                                            tl_wth, ct_ctq_wth
    real(kind=r_def), dimension(undf_w3),  intent(in)   :: height_w3,          &
                                                           rdz_tq_bl,          &
                                                           rhokh_bl
    real(kind=r_def), dimension(undf_wth), intent(in)   :: theta_in_wth,       &
                                                           exner_in_wth,       &
                                                           m_v_n, m_cl_n,      &
                                                           m_cf_n,             &
                                                           theta_star,         &
                                                           height_wth,         &
                                                           dt_conv,            &
                                                           dtrdz_tq_bl,        &
                                                           m_v, m_cl, m_cf
    integer(kind=i_def), dimension(undf_2d), intent(in) :: blend_height_tq
    integer(kind=i_def), dimension(undf_bl), intent(in) :: bl_type_ind

    real(kind=r_def), intent(inout) :: dqw1_2d(undf_2d)
    real(kind=r_def), intent(inout) :: dtl1_2d(undf_2d)
    real(kind=r_def), intent(inout) :: ct_ctq1_2d(undf_2d)

    !-----------------------------------------------------------------------
    ! Local variables for the kernel
    !-----------------------------------------------------------------------
    ! loop counters etc
    integer(i_def) :: k, i

    ! local switches and scalars
    logical :: l_correct

    ! profile fields from level 1 upwards
    real(r_bl), dimension(seg_len,1,nlayers) ::  t_latest, q_latest, &
         qcl_latest, qcf_latest, t, r_rho_levels

    ! profile field on boundary layer levels
    real(r_bl), dimension(seg_len,1,bl_levels) :: fqw, ftl, rhokh,       &
         dtrdz_charney_grid, rdz_charney_grid, qw, tl, dqw, dtl, ct_ctq, &
         dqw_nt, dtl_nt

    ! profile fields on u/v points and all levels
    real(r_bl), dimension(seg_len,1,nlayers) :: r_u, r_v

    ! profile fields on u/v points and BL levels
    real(r_bl), dimension(seg_len,1,bl_levels) :: taux, tauy, &
         dtrdz_u, dtrdz_v, rhokm_u, rhokm_v, cq_cm_u, cq_cm_v

    ! profile fields from level 2 upwards
    real(r_bl), dimension(seg_len,1,2:bl_levels) :: rdz_u, rdz_v

    ! profile fields from level 0 upwards
    real(r_bl), dimension(seg_len,1,0:nlayers) :: q, qcl, qcf, r_theta_levels

    ! single level real fields
    real(r_bl), dimension(seg_len,1) :: gamma1, gamma2, ctctq1_1, &
         dqw1_1, dtl1_1, cq_cm_u_1, du_1, cq_cm_v_1, dv_1

    ! single level integer fields
    integer(i_um), dimension(seg_len,1) :: k_blend_tq, k_blend_uv

    ! parameters for new BL solver
    real(r_bl) :: pnonl,p1,p2
    real(r_bl), dimension(seg_len) :: i1, e1, e2
    real(r_bl), parameter :: sqrt2 = sqrt(2.0_r_bl)

    !-----------------------------------------------------------------------
    ! Mapping of LFRic fields into UM variables
    !-----------------------------------------------------------------------
    ! assuming map_wth(1) points to level 0
    ! and map_w3(1) points to level 1
    !-----------------------------------------------------------------------
    do i = 1, seg_len
      do k = 1, nlayers
        ! Temperature
        t(i,1,k) = theta_in_wth(map_wth(1,i) + k) * &
                   exner_in_wth(map_wth(1,i) + k)
        ! height of rho levels from centre of planet
        r_rho_levels(i,1,k) = height_w3(map_w3(1,i) + k-1) + planet_radius
        ! water vapour mixing ratio
        q(i,1,k) = m_v_n(map_wth(1,i) + k)
        ! cloud liquid mixing ratio
        qcl(i,1,k) = m_cl_n(map_wth(1,i) + k)
      end do
    end do
    if (l_noice_in_turb) then
      do k = 1, nlayers
        do i = 1, seg_len
          qcf(i,1,k) = 0.0_r_bl
        end do
      end do
    else
      do i = 1, seg_len
        do k = 1, nlayers
          ! cloud ice mixing ratio
          qcf(i,1,k) = m_cf_n(map_wth(1,i) + k)
        end do
      end do
    end if

    ! surface height
    do i = 1, seg_len
      r_theta_levels(i,1,0) = height_wth(map_wth(1,i) + 0) + planet_radius
    end do

    !-----------------------------------------------------------------------
    ! Things passed from other parametrization schemes on this timestep
    !-----------------------------------------------------------------------
    do i = 1, seg_len
      do k = 1, bl_levels
        rhokh(i,1,k) = rhokh_bl(map_w3(1,i) + k-1)
        fqw(i,1,k) = moist_flux_bl(map_w3(1,i) + k-1)
        ftl(i,1,k) = heat_flux_bl(map_w3(1,i) + k-1)
        dtrdz_charney_grid(i,1,k) = dtrdz_tq_bl(map_wth(1,i) + k)
        rdz_charney_grid(i,1,k) = rdz_tq_bl(map_w3(1,i) + k-1)
      end do
      k_blend_tq(i,1) = blend_height_tq(map_2d(1,i))
    end do
    if (loop == 2) then
      do i = 1, seg_len
        do k = 1, bl_levels
          dqw(i,1,k) = dqw_wth(map_wth(1,i) + k)
          dtl(i,1,k) = dtl_wth(map_wth(1,i) + k)
        end do
      end do
    else
      do i = 1, seg_len
        dqw(i,1,1) = 0.0_r_bl
        dtl(i,1,1) = 0.0_r_bl
        ct_ctq(i,1,1) = 0.0_r_bl
      end do
    end if

    !-----------------------------------------------------------------------
    ! Increments / fields with increments added
    !-----------------------------------------------------------------------
    do i = 1, seg_len
      do k = 1, nlayers
        t_latest(i,1,k) = theta_star(map_wth(1,i) + k) * &
                          exner_in_wth(map_wth(1,i) + k) + &
                          dt_conv(map_wth(1,i) + k)
        q_latest(i,1,k)   = m_v(map_wth(1,i) + k)
        qcl_latest(i,1,k) = m_cl(map_wth(1,i) + k)
      end do
    end do
    if (l_noice_in_turb) then
      do k = 1, nlayers
        do i = 1, seg_len
          qcf_latest(i,1,k) = 0.0_r_bl
        end do
      end do
    else
      do i = 1, seg_len
        do k = 1, nlayers
          qcf_latest(i,1,k) = m_cf(map_wth(1,i) + k)
        end do
      end do
    end if

    do i = 1, seg_len
      p1=bl_type_ind(map_bl(1,i)+0)*pstb + &
           (1.0_r_bl-bl_type_ind(map_bl(1,i)+0))*puns
      p2=bl_type_ind(map_bl(1,i)+1)*pstb + &
           (1.0_r_bl-bl_type_ind(map_bl(1,i)+1))*puns
      pnonl=max(p1,p2)
      i1(i) = (1.0_r_bl+1.0_r_bl/sqrt2)*(1.0_r_bl+pnonl)
      e1(i) = (1.0_r_bl+1.0_r_bl/sqrt2)*( pnonl + (1.0_r_bl/sqrt2) + &
              sqrt(pnonl*(sqrt2-1.0_r_bl)+0.5_r_bl) )
      e2(i) = (1.0_r_bl+1.0_r_bl/sqrt2)*( pnonl+(1.0_r_bl/sqrt2) - &
              sqrt(pnonl*(sqrt2-1.0_r_bl)+0.5_r_bl))
      gamma1(i,1) = i1(i)
    end do

    if (loop == 1) then
      do i = 1, seg_len
        gamma2(i,1) = i1(i) - e1(i)
      end do
      l_correct = .false.
    else
      do i = 1, seg_len
        gamma2(i,1) = i1(i) - e2(i)
      end do
      l_correct = .true.
    end if

    call bdy_impl3 (                                                         &
         ! IN levels/switches
         bl_levels, l_correct,                                               &
         ! IN fields
         q, qcl, qcf, q_latest, qcl_latest, qcf_latest, t, t_latest,         &
         dtrdz_charney_grid, dtrdz_u, dtrdz_v, rhokh, rhokm_u, rhokm_v,      &
         rdz_charney_grid, rdz_u, rdz_v, gamma1, gamma2, real(alpha_cd,r_bl),&
         r_u, r_v, r_theta_levels, r_rho_levels, k_blend_tq,                 &
         k_blend_uv, k_blend_uv,                                             &
         ! INOUT fields
         fqw, ftl, taux, tauy, r_u, r_v, dqw, dtl,                           &
         ! OUT fields
         dqw_nt, dtl_nt, qw, tl, ct_ctq, cq_cm_u, cq_cm_v, cq_cm_u_1,        &
         cq_cm_v_1, du_1, dv_1, dqw1_1,dtl1_1,ctctq1_1                       &
         )

    do k = 1, bl_levels
      do i = 1, seg_len
        dqw_wth(map_wth(1,i) + k) = dqw(i,1,k)
        dtl_wth(map_wth(1,i) + k) = dtl(i,1,k)
        dqw_nt_wth(map_wth(1,i) + k) = dqw_nt(i,1,k)
        dtl_nt_wth(map_wth(1,i) + k) = dtl_nt(i,1,k)
        ct_ctq_wth(map_wth(1,i) + k) = ct_ctq(i,1,k)
      end do
    end do
    if (loop == 1) then
      do k = 1, bl_levels
        do i = 1, seg_len
          qw_wth(map_wth(1,i) + k) = qw(i,1,k)
          tl_wth(map_wth(1,i) + k) = tl(i,1,k)
        end do
      end do
      do i = 1, seg_len
        dqw1_2d(map_2d(1,i)) = dqw1_1(i,1)
        dtl1_2d(map_2d(1,i)) = dtl1_1(i,1)
        ct_ctq1_2d(map_2d(1,i)) = ctctq1_1(i,1)
      end do
    else !loop = 2
      do k = 0, bl_levels-1
        do i = 1, seg_len
          heat_flux_bl(map_w3(1,i)+k) = ftl(i,1,k+1)
          moist_flux_bl(map_w3(1,i)+k) = fqw(i,1,k+1)
        end do
      end do
      ! With specified fluxes only, we need to update the surface temperature
      ! using the final dtl value at level 1. We use the dtl1_2d field to
      ! carry this to jules_imp as it's normally only used on the 1st loop
      if (flux_bc_opt == specified_fluxes_only) then
        do i = 1, seg_len
          dtl1_2d(map_2d(1,i)) = dtl(i,1,1)
        end do
      end if
    end if

  end subroutine bl_imp_code

end module bl_imp_kernel_mod
