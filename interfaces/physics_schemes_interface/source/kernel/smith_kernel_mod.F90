!-----------------------------------------------------------------------------
! (c) Crown copyright 2017 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to large scale cloud scheme.
!>
module smith_kernel_mod

  use argument_mod,       only : arg_type,              &
                                 GH_FIELD, GH_REAL,     &
                                 GH_READ, GH_READWRITE, &
                                 GH_WRITE, DOMAIN, &
                                 ANY_DISCONTINUOUS_SPACE_1, &
                                 GH_INTEGER
  use constants_mod,      only : r_def, r_double, i_def, i_um, r_um
  use fs_continuity_mod,  only : W3, Wtheta
  use kernel_mod,         only : kernel_type

  implicit none

  private

  !> Kernel metadata type.
  !>
  type, public, extends(kernel_type) :: smith_kernel_type
    private
    type(arg_type) :: meta_args(14) = (/                                       &
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                    & ! theta_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,      W3),                        & ! exner_in_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                    & ! exner_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                    & ! rh_crit_wth
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,   ANY_DISCONTINUOUS_SPACE_1), & ! ntml_2d
         arg_type(GH_FIELD, GH_INTEGER, GH_READ,   ANY_DISCONTINUOUS_SPACE_1), & ! cumulus_2d
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                    & ! m_v
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                    & ! m_cl
         arg_type(GH_FIELD, GH_REAL, GH_READ,      WTHETA),                    & ! m_cf
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                    & ! cf_area
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                    & ! cf_ice
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                    & ! cf_liq
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, WTHETA),                    & ! cf_bulk
         arg_type(GH_FIELD, GH_REAL, GH_WRITE,     WTHETA)                     & ! theta_inc
        /)
    integer :: operates_on = DOMAIN
  contains
    procedure, nopass :: smith_code
  end type

  public :: smith_code

contains

  !> @brief Interface to the cloud scheme
  !> @details The UM large-scale cloud scheme:
  !>          determines the fraction of the grid-box that is covered by cloud
  !>          and the amount of liquid and ice condensate present in those clouds.
  !>          Here there is an interface to the Smith cloud scheme. Which is the
  !>          scheme described in UMDP 29.
  !> @param[in]     nlayers       Number of layers
  !> @param[in]     theta_in_wth  Predicted theta in its native space
  !> @param[in]     exner_in_w3   Pressure in the w3 space
  !> @param[in]     exner_in_wth  Exner Pressure in the theta space
  !> @param[in]     rh_crit_wth   Critical Relative Humidity
  !> @param[in]     ntml_2d       Boundary layer top level
  !> @param[in]     cumulus_2d    Cumulus flag
  !> @param[in,out] m_v           Vapour mixing ratio in wth
  !> @param[in,out] m_cl          Cloud liquid mixing ratio in wth
  !> @param[in]     m_cf          Total frozen mixing ratio incl snow in wth
  !> @param[in,out] cf_area       Area cloud fraction
  !> @param[in,out] cf_ice        Ice cloud fraction
  !> @param[in,out] cf_liq        Liquid cloud fraction
  !> @param[in,out] cf_bulk       Bulk cloud fraction
  !> @param[in,out] theta_inc     Increment to theta
  !> @param[in]     ndf_wth       Number of degrees of freedom per cell for potential temperature space
  !> @param[in]     undf_wth      Number unique of degrees of freedom  for potential temperature space
  !> @param[in]     map_wth       Dofmap for the cell at the base of the column for potential temperature space
  !> @param[in]     ndf_w3        Number of degrees of freedom per cell for density space
  !> @param[in]     undf_w3       Number unique of degrees of freedom  for density space
  !> @param[in]     map_w3        Dofmap for the cell at the base of the column for density space
  !> @param[in]     ndf_2d        Number of degrees of freedom per cell for 2D fields
  !> @param[in]     undf_2d       Number unique of degrees of freedom  for 2D fields
  !> @param[in]     map_2d        Dofmap for the cell at the base of the column for 2D fields

  subroutine smith_code(nlayers,      &
                        seg_len,      &
                        theta_in_wth, &
                        exner_in_w3,  &
                        exner_in_wth, &
                        rh_crit_wth,  &
                        ntml_2d,      &
                        cumulus_2d,   &
                        m_v,          &
                        m_cl,         &
                        m_cf,         &
                        cf_area,      &
                        cf_ice,       &
                        cf_liq,       &
                        cf_bulk,      &
                        theta_inc,    &
                        ndf_wth,      &
                        undf_wth,     &
                        map_wth,      &
                        ndf_w3,       &
                        undf_w3,      &
                        map_w3,       &
                        ndf_2d,       &
                        undf_2d,      &
                        map_2d)

    !---------------------------------------
    ! UM modules
    !---------------------------------------
    ! Structures holding diagnostic arrays - not used

    ! Other modules containing stuff passed to CLD
    use nlsizes_namelist_mod, only: bl_levels
    use planet_constants_mod, only: p_zero, kappa, cp
    use water_constants_mod,  only: lc
    use ls_arcld_mod,         only: ls_arcld
    use gen_phys_inputs_mod,  only: l_mr_physics

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)     :: nlayers, seg_len
    integer(kind=i_def), intent(in)     :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in)     :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in)     :: ndf_2d, undf_2d

    integer(kind=i_def), intent(in),    dimension(ndf_wth,seg_len)  :: map_wth
    integer(kind=i_def), intent(in),    dimension(ndf_w3,seg_len)   :: map_w3
    integer(kind=i_def), intent(in),    dimension(ndf_2d,seg_len)   :: map_2d

    real(kind=r_def),    intent(in),    dimension(undf_w3)  :: exner_in_w3
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: exner_in_wth
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: theta_in_wth
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: rh_crit_wth
    integer(kind=i_def), intent(in),    dimension(undf_2d)  :: ntml_2d, cumulus_2d
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: m_v
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: m_cl
    real(kind=r_def),    intent(in),    dimension(undf_wth) :: m_cf
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: cf_area
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: cf_ice
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: cf_liq
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: cf_bulk
    real(kind=r_def),    intent(inout), dimension(undf_wth) :: theta_inc

    ! Local variables for the kernel
    integer(i_um) :: k, i

    real(r_def) :: dmv1(seg_len)

    ! profile fields from level 1 upwards
    real(r_um), dimension(seg_len,1,nlayers) ::                      &
         cf_inout, cfl_inout, cff_inout, area_cloud_fraction, rhcpt, &
         qt, qcl_out, qcf_in,tl

    ! profile fields from level 0 upwards
    real(r_um), dimension(seg_len,1,0:nlayers) :: &
         p_rho_minus_one, p_theta_levels

    ! single level fields
    real(r_um), dimension(seg_len,1) :: fv_cos_theta_latitude

    logical, dimension(seg_len,1) :: cumulus
    ! single level fields
    integer(i_um), dimension(seg_len,1) :: ntml
    integer(i_um) :: errorstatus, large_levels, levels_per_level

    ! Determine number of sublevels for vertical gradient area cloud
    ! Want an odd number of sublevels per level: 3 is hardwired in do loops
    levels_per_level = 3
    large_levels = ((nlayers - 2)* levels_per_level) + 2

    do i = 1, seg_len
      cumulus(i,1) = (cumulus_2d(map_2d(1,i)) == 1_i_def)
      ntml(i,1) = ntml_2d(map_2d(1,i))
    end do

    errorstatus=0

    do i = 1, seg_len
      do k = 1, nlayers
        ! liquid temperature on theta levels
        tl(i,1,k) = ( theta_in_wth(map_wth(1,i) + k)   &
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
        ! 3D RH_crit field
        rhcpt(i,1,k) = rh_crit_wth(map_wth(1,i) + k)
      end do
    end do

    do i = 1, seg_len
      ! pressure on theta levels
      p_theta_levels(i,1,0) = p_zero*(exner_in_wth(map_wth(1,i) + 0))**(1.0_r_def/kappa)
      ! pressure on rho levels
      p_rho_minus_one(i,1,0) = p_zero*(exner_in_wth(map_wth(1,i) + 0))**(1.0_r_def/kappa)
      do k = 1, nlayers-1
        ! pressure on theta levels
        p_theta_levels(i,1,k) = p_zero*(exner_in_wth(map_wth(1,i) + k))**(1.0_r_def/kappa)
        ! pressure on rho levels without level 1
        p_rho_minus_one(i,1,k) = p_zero*(exner_in_w3(map_w3(1,i) + k))**(1.0_r_def/kappa)
      end do
      ! pressure on theta levels
      p_theta_levels(i,1,nlayers) = p_zero*(exner_in_wth(map_wth(1,i) + nlayers))** &
                    (1.0_r_def/kappa)
      ! pressure on rho levels
      p_rho_minus_one(i,1,nlayers) = 0.0_r_def
    end do

    call ls_arcld( p_theta_levels, rhcpt, p_rho_minus_one,            &
                 seg_len, 1, bl_levels,                               &
                 levels_per_level, large_levels,                      &
                 fv_cos_theta_latitude,                               &
                 ntml, cumulus, l_mr_physics, qcf_in,                 &
                 tl, qt, qcl_out,                                     &
                 area_cloud_fraction,  cf_inout,                      &
                 cfl_inout, cff_inout ,                               &
                 errorstatus)

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
        m_v(map_wth(1,i) + k)  = qt(i,1,k)
        ! cloud liquid water on theta levels
        m_cl(map_wth(1,i) + k) = qcl_out(i,1,k)
        ! cloud fractions on theta levels
        cf_bulk(map_wth(1,i) + k) = cf_inout(i,1,k)
        cf_liq(map_wth(1,i) + k) = cfl_inout(i,1,k)
        cf_ice(map_wth(1,i) + k) = cff_inout(i,1,k)
        cf_area(map_wth(1,i) + k) = area_cloud_fraction(i,1,k)
      end do
    end do
    do i = 1, seg_len
      theta_inc(map_wth(1,i) + 0) = theta_inc(map_wth(1,i) + 1)
      m_v(map_wth(1,i) + 0)  = m_v(map_wth(1,i) + 0) + m_v(map_wth(1,i) + 1) &
                             - dmv1(i)
      m_cl(map_wth(1,i) + 0) = m_cl(map_wth(1,i) + 1)
      cf_bulk(map_wth(1,i) + 0) = cf_bulk(map_wth(1,i) + 1)
      cf_liq(map_wth(1,i) + 0) = cf_liq(map_wth(1,i) + 1)
      cf_ice(map_wth(1,i) + 0) = cf_ice(map_wth(1,i) + 1)
      cf_area(map_wth(1,i) + 0) = cf_area(map_wth(1,i) + 1)
    end do

  end subroutine smith_code

end module smith_kernel_mod
