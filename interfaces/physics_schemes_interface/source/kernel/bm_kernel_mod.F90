!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the bimodal cloud scheme
!>
module bm_kernel_mod

  use argument_mod,       only : arg_type,              &
                                 GH_FIELD, GH_REAL,     &
                                 GH_READ, GH_READWRITE, &
                                 GH_INTEGER,            &
                                 ANY_DISCONTINUOUS_SPACE_1, &
                                 ANY_DISCONTINUOUS_SPACE_9, &
                                 GH_WRITE, DOMAIN
  use constants_mod,      only : r_def, i_def, i_um, r_um
  use fs_continuity_mod,  only : W3, Wtheta
  use kernel_mod,         only : kernel_type
  use empty_data_mod,     only : empty_real_data

  implicit none

  private

  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: bm_kernel_type
    private
    type(arg_type) :: meta_args(26) = (/                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! theta_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,      W3),     & ! exner_in_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! exner_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! dsldzm
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! mix_len_bm
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! wvar
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! tau_dec_bm
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! tau_hom_bm
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! tau_mph_bm
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! height_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! gradrinr
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! zh
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! zhsc
         arg_type(GH_FIELD, GH_REAL, GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! inv_depth
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,   ANY_DISCONTINUOUS_SPACE_9), & ! bl_type_ind
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! m_v
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! m_cl
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA), & ! m_cf
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! cf_area
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! cf_ice
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! cf_liq
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA), & ! cf_bulk
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     WTHETA), & ! theta_inc
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     WTHETA), & ! sskew_bm
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     WTHETA), & ! svar_bm
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     WTHETA)  & ! svar_tb
         /)
    integer :: operates_on = DOMAIN
  contains
    procedure, nopass :: bm_code
  end type

  public :: bm_code

contains

  !> @brief Interface to the bimodal cloud scheme
  !> @details The UM large-scale cloud scheme:
  !>          determines the fraction of the grid-box that is covered by ice and
  !>          liquid cloud and the amount of liquid condensate present in those clouds.
  !>          Here there is an interface to the bimodal cloud scheme. Which is the
  !>          scheme described in UMDP 39.
  !> @param[in]     nlayers       Number of layers
  !> @param[in]     theta_in_wth  Predicted theta in its native space
  !> @param[in]     exner_in_w3   Pressure in the w3 space
  !> @param[in]     exner_in_wth  Exner Pressure in the theta space
  !> @param[in]     dsldzm        Liquid potential temperature gradient in wth
  !> @param[in]     mix_len_bm    Turb length-scale for bimodal in wth
  !> @param[in]     wvar          Vertical velocity variance in wth
  !> @param[in]     tau_dec_bm    Decorrelation time scale in wth
  !> @param[in]     tau_hom_bm    Homogenisation time scale in wth
  !> @param[in]     tau_mph_bm    Phase-relaxation time scale in wth
  !> @param[in]     height_wth    Theta-level height in wth
  !> @param[in]     gradrinr      Gradient Richardson Number in wth
  !> @param[in]     zh            Mixed-layer height
  !> @param[in]     zhsc          Decoupled layer height
  !> @param[in]     inv_depth     Depth of BL top inversion layer
  !> @param[in]     bl_type_ind   Diagnosed BL types
  !> @param[in,out] m_v           Vapour mixing ratio in wth
  !> @param[in,out] m_cl          Cloud liquid mixing ratio in wth
  !> @param[in]     m_cf          Frozen liquid mixing ratio in wth
  !> @param[in,out] cf_area       Area cloud fraction
  !> @param[in,out] cf_ice        Ice cloud fraction
  !> @param[in,out] cf_liq        Liquid cloud fraction
  !> @param[in,out] cf_bulk       Bulk cloud fraction
  !> @param[in,out] theta_inc     Increment to theta
  !> @param[in,out] sskew_bm      Bimodal skewness of SD PDF
  !> @param[in,out] svar_bm       Bimodal variance of SD PDF
  !> @param[in,out] svar_tb       Unimodal variance of SD PDF
  !> @param[in]     ndf_wth       Number of degrees of freedom per cell for potential temperature space
  !> @param[in]     undf_wth      Number unique of degrees of freedom  for potential temperature space
  !> @param[in]     map_wth       Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_w3        Number of degrees of freedom per cell for density space
  !> @param[in]     undf_w3       Number unique of degrees of freedom  for density space
  !> @param[in]     map_w3        Dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_2d        Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d       Number of unique DOFs for 2D fields
  !> @param[in]     map_2d        Dofmap for the cell at the base of the column for 2D fields
  !> @param[in]     ndf_bl        Number of DOFs per cell for BL types
  !> @param[in]     undf_bl       Number of total DOFs for BL types
  !> @param[in]     map_bl        Dofmap for cell for BL types

  subroutine bm_code(nlayers,      &
                     seg_len,      &
                     theta_in_wth, &
                     exner_in_w3,  &
                     exner_in_wth, &
                     dsldzm,       &
                     mix_len_bm,   &
                     wvar,         &
                     tau_dec_bm,   &
                     tau_hom_bm,   &
                     tau_mph_bm,   &
                     height_wth,   &
                     gradrinr,     &
                     zh,           &
                     zhsc,         &
                     inv_depth,    &
                     bl_type_ind,  &
                     m_v,          &
                     m_cl,         &
                     m_cf,         &
                     cf_area,      &
                     cf_ice,       &
                     cf_liq,       &
                     cf_bulk,      &
                     theta_inc,    &
                     sskew_bm,     &
                     svar_bm,      &
                     svar_tb,      &
                     ndf_wth,      &
                     undf_wth,     &
                     map_wth,      &
                     ndf_w3,       &
                     undf_w3,      &
                     map_w3,       &
                     ndf_2d,       &
                     undf_2d,      &
                     map_2d,       &
                     ndf_bl,       &
                     undf_bl,      &
                     map_bl)

    !---------------------------------------
    ! UM modules
    !---------------------------------------
    ! Structures holding diagnostic arrays - not used

    ! Other modules containing stuff passed to CLD
    use cloud_config_mod,     only: i_bm_ez_opt, i_bm_ez_opt_entpar
    use planet_constants_mod, only: p_zero, kappa, cp
    use water_constants_mod,  only: lc
    use bm_ctl_mod,           only: bm_ctl
    use gen_phys_inputs_mod,  only: l_mr_physics

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)     :: nlayers, seg_len
    integer(kind=i_def), intent(in)     :: ndf_wth, ndf_w3, ndf_2d
    integer(kind=i_def), intent(in)     :: undf_wth, undf_w3, undf_2d
    integer(kind=i_def), intent(in)     :: ndf_bl, undf_bl

    integer(kind=i_def), intent(in),    dimension(ndf_wth,seg_len)  :: map_wth
    integer(kind=i_def), intent(in),    dimension(ndf_w3,seg_len)   :: map_w3
    integer(kind=i_def), intent(in),    dimension(ndf_2d,seg_len)   :: map_2d
    integer(kind=i_def), intent(in),    dimension(ndf_bl,seg_len)   :: map_bl

    real(kind=r_def),    intent(in),    dimension(undf_w3)  :: exner_in_w3
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: exner_in_wth
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: theta_in_wth
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: dsldzm
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: mix_len_bm
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: wvar
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: gradrinr
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: tau_dec_bm
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: tau_hom_bm
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: tau_mph_bm
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: height_wth
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: m_v
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: m_cl
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: m_cf
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: cf_area
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: cf_ice
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: cf_liq
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: cf_bulk
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: theta_inc
    real(kind=r_def),    intent(inout), dimension(:), pointer :: sskew_bm
    real(kind=r_def),    intent(inout), dimension(:), pointer :: svar_bm
    real(kind=r_def),    intent(inout), dimension(:), pointer :: svar_tb

    real(kind=r_def),    intent(in),    dimension(undf_2d)  :: zh
    real(kind=r_def),    intent(in),    dimension(undf_2d)  :: zhsc
    real(kind=r_def),    intent(in),    dimension(undf_2d)  :: inv_depth
    integer(kind=i_def), intent(in),    dimension(undf_bl)  :: bl_type_ind

    ! Local variables for the kernel
    integer(i_um) :: k, i

    real(r_def) :: dmv1(seg_len)

    ! profile fields from level 1 upwards
    real(r_um), dimension(seg_len,1,nlayers) :: cf_inout, cfl_inout, cff_inout,&
         area_cloud_fraction, qt, qcl_out, qcf_in, tl,                         &
         tgrad_in, mix_len_in, tau_dec_in, tau_hom_in, tau_mph_in,             &
         z_theta, wvar_in, gradrinr_in, sskew_out,                             &
         svar_turb_out, svar_bm_out

    real(r_um), dimension(seg_len,1) :: zh_in, zhsc_in, dzh_in, bl_type_7_in

    ! profile fields from level 1 upwards
    real(r_um), dimension(seg_len,1,nlayers) :: p_theta_levels

    ! Extra bimodal scheme diagnostics not yet implemented in LFRic
    real(r_um), dimension(1,1,1,1) :: entzone,                                 &
                                      sl_modes, qw_modes, rh_modes, sd_modes
    ! Switch set to false to indicate not to calculate the above
    logical, parameter :: l_calc_diag = .FALSE.

    ! error status
    integer(i_um) :: errorstatus

    errorstatus = 0

    do i = 1, seg_len
      do k = 1, nlayers
        ! liquid temperature on theta levels
        tl(i,1,k) = (theta_in_wth(map_wth(1,i) + k) &
                    * exner_in_wth(map_wth(1,i)+ k)) - &
                    (lc * m_cl(map_wth(1,i) + k)) / cp
        ! total water and ice water on theta levels
        qt(i,1,k) =  m_v(map_wth(1,i) + k) + m_cl(map_wth(1,i) + k)
        qcf_in(i,1,k) = m_cf(map_wth(1,i) + k)
        ! cloud fields
        cf_inout(i,1,k) = cf_bulk(map_wth(1,i) + k)
        cff_inout(i,1,k) = cf_ice(map_wth(1,i) + k)
        cfl_inout(i,1,k) = cf_liq(map_wth(1,i) + k)
        area_cloud_fraction(i,1,k) = cf_area(map_wth(1,i) + k)
        tau_dec_in(i,1,k) = tau_dec_bm(map_wth(1,i) + k)
        tau_hom_in(i,1,k) = tau_hom_bm(map_wth(1,i) + k)
        tau_mph_in(i,1,k) = tau_mph_bm(map_wth(1,i) + k)
        z_theta(i,1,k) = height_wth(map_wth(1,i) + k) &
                       - height_wth(map_wth(1,i) + 0)
        wvar_in(i,1,k)     = wvar(map_wth(1,i) + k)
        gradrinr_in(i,1,k) = gradrinr(map_wth(1,i) + k)
      end do
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

    ! 2d fields for gradrinr-based entrainment zones
    do i = 1, seg_len
      zh_in(i,1)        = zh(map_2d(1,i))
      zhsc_in(i,1)       = zhsc(map_2d(1,i))
      dzh_in(i,1)       = real(inv_depth(map_2d(1,i)), r_um)
      bl_type_7_in(i,1) = bl_type_ind(map_bl(1,i)+6)
    end do

    do i = 1, seg_len
      do k = 1, nlayers
        ! pressure on theta levels
        p_theta_levels(i,1,k) = p_zero*(exner_in_wth(map_wth(1,i) + k))**(1.0_r_def/kappa)
      end do
    end do

    call bm_ctl( p_theta_levels, tgrad_in, wvar_in,                   &
                 tau_dec_in, tau_hom_in, tau_mph_in, z_theta,         &
                 gradrinr_in, mix_len_in, zh_in, zhsc_in, dzh_in,     &
                 bl_type_7_in,                                        &
                 nlayers, l_mr_physics,                               &
                 tl, cf_inout, qt, qcf_in, qcl_out,                   &
                 cfl_inout, cff_inout,                                &
                 sskew_out, svar_turb_out, svar_bm_out, entzone,      &
                 sl_modes, qw_modes, rh_modes, sd_modes,              &
                 l_calc_diag, errorstatus )

    ! update main model prognostics
    !-----------------------------------------------------------------------
    ! Save old value of m_v at level 1 for level 0 increment
    do i = 1, seg_len
      dmv1(i) = m_v(map_wth(1,i) + 1)
    end do

    do k = 1, nlayers
      do i = 1, seg_len
        ! potential temperature increment on theta levels
        theta_inc(map_wth(1,i) + k) = tl(i,1,k)/exner_in_wth(map_wth(1,i) + k) &
                                    - theta_in_wth(map_wth(1,i) + k)
        ! water vapour on theta levels
        m_v(map_wth(1,i) + k)       = qt(i,1,k)
        ! cloud liquid water on theta levels
        m_cl(map_wth(1,i) + k)      = qcl_out(i,1,k)
        ! cloud fractions on theta levels
        cf_bulk(map_wth(1,i) + k)   = cf_inout(i,1,k)
        cf_liq(map_wth(1,i) + k)    = cfl_inout(i,1,k)
        cf_ice(map_wth(1,i) + k)    = cff_inout(i,1,k)
        cf_area(map_wth(1,i) + k)   = cf_inout(i,1,k)
      end do
    end do
    if (.not. associated(sskew_bm, empty_real_data)) then
      do k = 1, nlayers
        do i = 1, seg_len
          sskew_bm(map_wth(1,i) + k)  = sskew_out(i,1,k)
        end do
      end do
    end if
    if (.not. associated(svar_bm, empty_real_data)) then
      do k = 1, nlayers
        do i = 1, seg_len
          svar_bm(map_wth(1,i) + k)   = svar_bm_out(i,1,k)
        end do
      end do
    end if
    if (.not. associated(svar_tb, empty_real_data)) then
      do k = 1, nlayers
        do i = 1, seg_len
          svar_tb(map_wth(1,i) + k)   = svar_turb_out(i,1,k)
        end do
      end do
    end if
    do i = 1, seg_len
      theta_inc(map_wth(1,i) + 0) = theta_inc(map_wth(1,i) + 1)
      m_v(map_wth(1,i) + 0)       = m_v(map_wth(1,i) + 0) &
                                  + m_v(map_wth(1,i) + 1) - dmv1(i)
      m_cl(map_wth(1,i) + 0)      = m_cl(map_wth(1,i) + 1)
      cf_bulk(map_wth(1,i) + 0)   = cf_bulk(map_wth(1,i) + 1)
      cf_liq(map_wth(1,i) + 0)    = cf_liq(map_wth(1,i) + 1)
      cf_ice(map_wth(1,i) + 0)    = cf_ice(map_wth(1,i) + 1)
      cf_area(map_wth(1,i) + 0)   = cf_area(map_wth(1,i) + 1)
    end do

  end subroutine bm_code

end module bm_kernel_mod
