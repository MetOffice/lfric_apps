!-----------------------------------------------------------------------------
! (c) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Interface to the UM orographic gravity wave and blocking drag scheme
!>
!>
module orographic_drag_kernel_mod

  use argument_mod,               only: arg_type,          &
                                        GH_FIELD, GH_REAL, &
                                        GH_READ, GH_WRITE, &
                                        CELL_COLUMN,       &
                                        ANY_DISCONTINUOUS_SPACE_1

  use planet_constants_mod,       only: p_zero, kappa
  use constants_mod,              only: r_def, r_um, i_def, i_um, pi
  use empty_data_mod,             only: empty_real_data
  use fs_continuity_mod,          only: W3, Wtheta
  use kernel_mod,                 only: kernel_type
  use orographic_drag_config_mod, only: cd_flow_blocking,            &
                                        gwd_scaling,                 &
                                        fr_crit_gwd,                 &
                                        fr_sat_gwd,                  &
                                        mountain_height_scaling,     &
                                        orographic_gwd_heating,      &
                                        orographic_blocking_heating, &
                                        vertical_smoothing

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> Metadata describing the kernel to PSyclone
  !>
  type, public, extends(kernel_type) :: orographic_drag_kernel_type
    private
    type(arg_type) :: meta_args(24) = (/                                   &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3),                        & ! du_orog_blk, u wind increment blocking
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3),                        & ! dv_orog_blk, v wind increment blocking
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3),                        & ! du_orog_gwd, u wind increment gwd
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, W3),                        & ! dv_orog_gwd, v wind increment gwd
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta),                    & ! dtemp_orog_blk, T increment blocking
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta),                    & ! dtemp_orog_gwd, T increment gwd
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        & ! u_in_w3, zonal wind
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        & ! v_in_w3, meridional wind
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        & ! wetrho_in_w3, wet density in w3
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! theta, theta in wtheta
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! exner_in_wth
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! sd_orog
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! grad_xx_orog
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! grad_xy_orog
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), & ! grad_yy_orog
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! mr_v
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! mr_cl
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! mr_cf
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        & ! height_w3
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    & ! height_wth
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta),                    & ! taux_orog_blk
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta),                    & ! tauy_orog_blk
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta),                    & ! taux_orog_gwd
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta)                     & ! tauy_orog_gwd
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: orographic_drag_kernel_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: orographic_drag_kernel_code

contains

  !> @brief   Call the UM orographic gravity wave and blocking drag scheme
  !> @details This code calls the UM orographic gravity wave and blocking
  !>          drag scheme, which calculates the zonal and meridional winds and
  !>          temperature increments from parametrized orographic drag.
  !! @param[in]     nlayers        Integer the number of layers
  !! @param[in,out] du_orog_blk    u increment from blocking
  !! @param[in,out] dv_orog_blk    v increment from blocking
  !! @param[in,out] du_orog_gwd    u increment from gravity wave drag
  !! @param[in,out] dv_orog_gwd    v increment from gravity wave drag
  !! @param[in,out] dtemp_orog_blk T increment from blocking
  !! @param[in,out] dtemp_orog_gwd T increment from gravity wave drag
  !! @param[in]     u_in_w3        Zonal wind
  !! @param[in]     v_in_w3        Meridional wind
  !! @param[in]     wetrho_in_w3   Moist density
  !! @param[in]     theta_in_wth   Potential temperature
  !! @param[in]     exner_in_wth   Exner pressure
  !! @param[in]     sd_orog        Standard deviation of sub-grid orog
  !! @param[in]     grad_xx_orog   (dh/dx)**2
  !! @param[in]     grad_xy_orog   (dh/dx)*(dh/dy)
  !! @param[in]     grad_yy_orog   (dh/dy)**2
  !! @param[in]     mr_v           Water vapour mixing ratio
  !! @param[in]     mr_cl          Cloud liquid mixing ratio
  !! @param[in]     mr_cf          Cloud frozen mixing ratio
  !! @param[in]     height_w3      Height at rho levels
  !! @param[in]     height_wth     Height at theta levels
  !! @param[in,out] taux_orog_blk  x-stress from blocking
  !! @param[in,out] tauy_orog_blk  y-stress from blocking
  !! @param[in,out] taux_orog_gwd  x-stress from orographic gwd
  !! @param[in,out] tauy_orog_gwd  y-stress from orographic gwd
  !! @param[in]     ndf_w3         Number of degrees of freedom per cell for wth
  !! @param[in]     undf_w3        Number of unique degrees of freedom for wth
  !! @param[in]     map_w3         Dofmap for the cell at w3
  !! @param[in]     ndf_wth        Number of degrees of freedom per cell for wth
  !! @param[in]     undf_wth       Number of unique degrees of freedom for wth
  !! @param[in]     map_wth        Dofmap for the cell at wth
  !! @param[in]     ndf_2d         Number of degrees of freedom per cell for 2d
  !! @param[in]     undf_2d        Number of unique degrees of freedom for 2d
  !! @param[in]     map_2d         Dofmap for the 2d cell
  !>
  subroutine orographic_drag_kernel_code(                                  &
                        max_pos_cells, nlayers, du_orog_blk, dv_orog_blk,  &
                        du_orog_gwd, dv_orog_gwd,                          &
                        dtemp_orog_blk, dtemp_orog_gwd, u_in_w3, v_in_w3,  &
                        wetrho_in_w3, theta_in_wth, exner_in_wth, sd_orog, &
                        grad_xx_orog, grad_xy_orog, grad_yy_orog,          &
                        mr_v, mr_cl, mr_cf,                                &
                        height_w3, height_wth,                             &
                        ! Diagnostics
                        taux_orog_blk, tauy_orog_blk,                      &
                        taux_orog_gwd, tauy_orog_gwd,                      &
                        ! Spatial information
                        ndf_w3, undf_w3, map_w3, ndf_wth, undf_wth,        &
                        map_wth, ndf_2d, undf_2d, map_2d, seg_len, cell_index)

    !---------------------------------------
    ! UM modules
    !---------------------------------------
    use timestep_mod, only: timestep
    use gw_block_mod, only: gw_block
    use gw_setup_mod, only: gw_setup
    use gw_wave_mod,  only: gw_wave

    implicit none

    !----------------------------------------------------------------------
    ! Arguments
    !----------------------------------------------------------------------
    integer(i_def), intent(in) :: nlayers, max_pos_cells, seg_len !nlayers, total cells (whether orogrogaphically relevant or not), segment length
    integer(i_def), intent(in), dimension(seg_len) :: cell_index 
    integer(i_def), intent(in) :: ndf_w3, ndf_wth, ndf_2d
    integer(i_def), intent(in) :: undf_w3, undf_wth, undf_2d
    integer(i_def), intent(in), dimension(ndf_w3, max_pos_cells)  :: map_w3 
    integer(i_def), intent(in), dimension(ndf_wth, max_pos_cells) :: map_wth
    integer(i_def), intent(in), dimension(ndf_2d, max_pos_cells)  :: map_2d

    real(r_def), intent(inout), dimension(undf_w3)  :: du_orog_blk, du_orog_gwd, &
                                                       dv_orog_blk, dv_orog_gwd
    real(r_def), intent(inout), dimension(undf_wth) :: dtemp_orog_blk, dtemp_orog_gwd
    real(r_def), intent(in), dimension(undf_w3)     :: u_in_w3, v_in_w3, &
                                                       wetrho_in_w3
    real(r_def), intent(in), dimension(undf_wth)  :: theta_in_wth, exner_in_wth
    real(r_def), intent(in), dimension(undf_wth)  :: mr_v, mr_cl, mr_cf
    real(r_def), intent(in), dimension(undf_2d)   :: sd_orog,      &
                                                     grad_xx_orog, &
                                                     grad_xy_orog, &
                                                     grad_yy_orog

    real(r_def), intent(in), dimension(undf_w3)   :: height_w3
    real(r_def), intent(in), dimension(undf_wth)  :: height_wth

    real(r_def), pointer, intent(inout)::      &
            taux_orog_blk(:), tauy_orog_blk(:),&
            taux_orog_gwd(:), tauy_orog_gwd(:)

    !----------------------------------------------------------------------
    ! Local variables for input to the kernel
    !----------------------------------------------------------------------

    ! At present, only seg_len = 1 is permitted. This kernel will require work in
    ! order to make this generalisable to seg_len > 1.

    real(r_um), dimension(seg_len, nlayers) :: &
                    u_on_p, v_on_p,            & ! u and v at p points
                    theta,                     & ! potential temperature
                    temp,                      & ! temperature
                    press,                     & ! pressure
                    q,                         & ! water vapour mixing ratio
                    qcl,                       & ! cloud liquid mixing ratio
                    qcf,                       & ! cloud ice mixing ratio
                    dtemp_dt_blk,              & ! Temperature tendency from blocking
                    dtemp_dt_orog_gwd,         & ! Temperature tendency from gwd
                    wetrho,                    & ! wetrho
                    du_dt_blk, du_dt_orog_gwd, & ! zonal wind tendencies
                    dv_dt_blk, dv_dt_orog_gwd, & ! meridional wind tendencies
                    z_rho_levels,              & ! height above the surface
                    z_theta_levels

    real(r_um), dimension(seg_len, 0:nlayers) :: &
                       tau_x_blk, tau_y_blk,     &
                       tau_x_orog_gwd, tau_y_orog_gwd

    ! The following variables are used in the nonhydrostatic option of the scheme
    ! which is hardcoded to .false. in the UM code.
    ! We set them here just so that they can be passed in to the subroutines.
    real(r_um), parameter :: delta_lambda=1.0_r_um, delta_phi=1.0_r_um
    real(r_um), dimension(seg_len) :: latitude
    logical, parameter :: nonhydro=.false., dynbeta=.false.

    ! Output from gw_setup
    real(r_um), dimension(seg_len, nlayers) :: &
                    nsq,       & ! moist Brunt-Vaisala frequency
                    nsq_dry,   & ! dry Brunt-Vaisala frequency
                    nsq_unsat, & ! unsaturated moist Brunt-Vaisala frequency
                    nsq_sat,   & ! saturated moist Brunt-Vaisala frequency
                    dzcond       ! Ascent to Lifting Condensation Level

    logical, dimension(seg_len, nlayers) :: l_lapse ! logical array

    real(r_um), dimension(seg_len) :: & ! Low-level averaged ...
                ulow, vlow,           & ! ... Zonal and meridional wind
                rholow,               & ! ... Density
                nlow,                 & ! ... Brunt-Vaisala frequency
                psilow,               & ! ... Angle between wind and major axis of topography
                psi1,                 & ! ... Wind angle relative to x-axis
                modu                    ! ... Wind speed

    ! Input subgrid orographic characteristics
    real(r_um), dimension(seg_len) ::               &
                sd, grad_xx, grad_xy, grad_yy,      &
                orog_f1, orog_f2, orog_f3, orog_amp

    ! Computed subgrid orographic characteristics
    real(r_um), dimension(seg_len) :: mt_high, slope, anis, banis,  &
                                      canis, mtdir

    ! Computed orographic variables
    real(r_um), dimension(seg_len) :: zb

    ! local namelist inputs
    real(r_um) :: fbcd, gsharp, gwd_frc, gwd_fsat, nsigma
    logical :: l_fb_heating, l_gw_heating, l_smooth

    integer(i_um) :: k, i

    integer(i_um), dimension(seg_len) :: ktop, kbot

    ! Flags for diagnostics that are not used in LFRic
    logical, parameter ::                                                     &
               u_s_d_on = .false.,  v_s_d_on= .false., nsq_s_d_on= .false.,   &
               du_dt_diag_on = .false., dv_dt_diag_on = .false.,              &
               stress_u_on = .false., stress_v_on = .false.,                  &
               fr_d_on = .false., bld_d_on= .false., tausx_d_on = .false.,    &
               tausy_d_on = .false., bldt_d_on= .false.

    ! Flags for diagnostics that are used in LFRic
    logical :: tau_x_blk_flag, tau_y_blk_flag,     &
               tau_x_orog_gwd_flag, tau_y_orog_gwd_flag

    ! Flag for determining if scheme needs to be computed
    logical, dimension(seg_len) :: drag

    ! Diagnostics not used
    real(r_um), dimension(seg_len) :: &
               u_s_d, v_s_d, nsq_s_d, &
               fr_d, bld_d, bldt_d,   &
               tausx_d, tausy_d


    ! Diagnostics not used
    real(r_um), dimension(seg_len, nlayers) :: &
                       du_dt_diag, dv_dt_diag, &
                       stress_u, stress_v


    !dhc - cant be parameter
    latitude=1.0_r_um
    !-----------------------------------------------------------------------
    ! Initialise arrays and call UM code
    !-----------------------------------------------------------------------
    drag = .true.

    ! Increments need to be intialised to zero because they are added onto
    ! previous increments in UM code (not overwritten).

    do i= 1, seg_len
      do k = 1, nlayers
        du_dt_orog_gwd(i,k) = 0.0_r_um
        dv_dt_orog_gwd(i,k) = 0.0_r_um
        du_dt_blk(i,k) = 0.0_r_um
        dv_dt_blk(i,k) = 0.0_r_um
        dtemp_dt_blk(i,k) = 0.0_r_um
        dtemp_dt_orog_gwd(i,k) = 0.0_r_um
      end do
    end do

    !need loop i (points)
    ! Recasting fields to UM precision
    do i= 1, seg_len
      do k = 1, nlayers
        u_on_p(i,k) = real(u_in_w3(map_w3(1,cell_index(i)) + k-1), r_um)
        v_on_p(i,k) = real(v_in_w3(map_w3(1,cell_index(i)) + k-1), r_um)

        wetrho(i,k) = real(wetrho_in_w3(map_w3(1,cell_index(i)) + k-1), r_um)

        theta(i,k)  = real(theta_in_wth(map_wth(1,cell_index(i)) + k), r_um)
        temp(i,k)  = real(exner_in_wth(map_wth(1,cell_index(i)) + k)*theta_in_wth(map_wth(1,cell_index(i)) + k), r_um)

        ! Pressure on layer boundaries (note, top layer is set to zero below)
        press(i,k) = real(p_zero*(exner_in_wth(map_wth(1,cell_index(i)) + k)) &
                                      **(1.0_r_um/kappa), r_um)

        ! water vapour mixing ratio
        q(i,k) = mr_v(map_wth(1,cell_index(i)) + k)
        ! cloud liquid mixing ratio
        qcl(i,k) = mr_cl(map_wth(1,cell_index(i)) + k)
        ! cloud ice mixing ratio
        qcf(i,k) = mr_cf(map_wth(1,cell_index(i)) + k)

        l_lapse(i,k) = .false.

        z_rho_levels(i,k)   = real(height_w3(map_w3(1,cell_index(i)) + k-1) - height_wth(map_wth(1,cell_index(i)) + 0), r_um)
        z_theta_levels(i,k) = real(height_wth(map_wth(1,cell_index(i)) + k) - height_wth(map_wth(1,cell_index(i)) + 0), r_um)
      end do   ! k
    end do !i

    do i=1, seg_len
      press(i,nlayers) = 0.0_r_um

      sd(i)      = real(sd_orog(map_2d(1,cell_index(i))), r_um)
      grad_xx(i) = real(grad_xx_orog(map_2d(1,cell_index(i))), r_um)
      grad_xy(i) = real(grad_xy_orog(map_2d(1,cell_index(i))), r_um)
      grad_yy(i) = real(grad_yy_orog(map_2d(1,cell_index(i))), r_um)

      ! Scale aware inputs (not currently used in LFRic)
      orog_f1(i)  = 0.0_r_um
      orog_f2(i)  = 0.0_r_um
      orog_f3(i)  = 0.0_r_um
      orog_amp(i) = 0.0_r_um
    end do !i
    
    ! Recasting of LFRic to UM namelist inputs
    fbcd         = real(cd_flow_blocking, r_um)
    gsharp       = real(gwd_scaling, r_um)
    gwd_frc      = real(fr_crit_gwd, r_um)
    gwd_fsat     = real(fr_sat_gwd, r_um)
    nsigma       = real(mountain_height_scaling, r_um)
    l_fb_heating = orographic_blocking_heating
    l_gw_heating = orographic_gwd_heating
    l_smooth     = vertical_smoothing

    ! Set stash flags and arrays
    if (.not. associated(taux_orog_blk, empty_real_data) ) then
      tau_x_blk_flag = .true.
    else
      tau_x_blk_flag = .false.
    end if
    if (.not. associated(tauy_orog_blk, empty_real_data) ) then
      tau_y_blk_flag = .true.
    else
      tau_y_blk_flag = .false.
    end if
    if (.not. associated(taux_orog_gwd, empty_real_data) ) then
      tau_x_orog_gwd_flag = .true.
    else
      tau_x_orog_gwd_flag = .false.
    end if
    if (.not. associated(tauy_orog_gwd, empty_real_data) ) then
      tau_y_orog_gwd_flag = .true.
    else
      tau_y_orog_gwd_flag = .false.
    end if

    ! Call routine to setup orographic drag fields
    call gw_setup(nlayers, seg_len, seg_len,                     &
                  ! Inputs
                  u_on_p, v_on_p,  wetrho, theta,                    &
                  ! Inputs to calculate moist buoyancy frequency
                  temp, q, qcl, qcf, press,                          &
                  ! Outputs from moist buoyancy frequency calculation
                  nsq, nsq_dry, nsq_unsat, nsq_sat,                  &
                  dzcond, l_lapse, kbot,                             &
                  ulow, vlow, rholow, nlow, psilow, psi1, modu,      &
                  ! Time-independent input
                  z_theta_levels, sd,                                &
                  grad_xx, grad_xy, grad_yy,                         &
                  ! Time-independent output
                  mt_high, slope, anis, banis, canis, mtdir, ktop,   &
                  drag, nsigma,                                      &
                  ! diagnostics
                  u_s_d, u_s_d_on, seg_len,                          &
                  v_s_d, v_s_d_on, seg_len)

    ! Call routine to compute orographic blocking depth and drag
    call gw_block(nlayers,seg_len,seg_len,timestep,u_on_p,v_on_p,    &
                  wetrho,nsq,ulow, vlow, modu,                           &
                  z_rho_levels,z_theta_levels,mt_high,                   &
                  sd,slope,anis,mtdir,zb,banis,canis,                    &
                  du_dt_blk,dv_dt_blk,                                   &
                  dtemp_dt_blk,fbcd,gwd_frc,drag,                        &
                  l_fb_heating,                                          &
                  ! diagnostics (not used)
                  du_dt_diag, seg_len,du_dt_diag_on,                     &
                  dv_dt_diag, seg_len,dv_dt_diag_on,                     &
                  stress_u, stress_u_on,seg_len,stress_u_on,             &
                  stress_v, stress_v_on,seg_len,                         &
                  tau_x_blk, seg_len, tau_x_blk_flag,                    &
                  tau_y_blk, seg_len, tau_y_blk_flag,                    &
                  fr_d,fr_d_on, seg_len,                                 &
                  bld_d,bld_d_on, seg_len,                               &
                  bldt_d,bldt_d_on, seg_len,                             &
                  tausx_d,tausx_d_on, seg_len,                           &
                  tausy_d,tausy_d_on, seg_len)

    ! Call routine to compute orographic gravity wave drag
    call gw_wave(nlayers,seg_len,seg_len,u_on_p,v_on_p,wetrho,   &
                 nsq_dry, nsq_unsat, nsq_sat, dzcond, l_lapse, kbot, &
                 ulow,vlow,rholow,psi1,psilow,nlow,modu,ktop,        &
                 z_rho_levels,z_theta_levels,delta_lambda,delta_phi, &
                 latitude,mt_high,sd,slope,zb,banis,canis,           &
                 orog_f1,orog_f2,orog_f3,orog_amp,                   &
                 du_dt_orog_gwd,dv_dt_orog_gwd,dtemp_dt_orog_gwd,    &
                 dynbeta,nonhydro,l_smooth,                          &
                 gwd_fsat,gsharp,drag,l_gw_heating,                  &
                 ! diagnostics (not used)
                 du_dt_diag,seg_len,du_dt_diag_on,du_dt_diag_on,     &
                 dv_dt_diag, seg_len,dv_dt_diag_on,                  &
                 stress_u, seg_len, stress_u_on, stress_u_on,        &
                 stress_v, seg_len, stress_v_on,                     &
                 ! diagnostics (used)
                 tau_x_orog_gwd, seg_len, tau_x_orog_gwd_flag,       &
                 tau_y_orog_gwd, seg_len, tau_y_orog_gwd_flag,       &
                 ! diagnostics (not used)
                 tausx_d, tausx_d_on, seg_len,                       &
                 tausy_d, tausy_d_on, seg_len,                       &
                 nsq_s_d, nsq_s_d_on, seg_len)

    ! Map variables back
    do i = 1, seg_len
      do k = 1, nlayers
        du_orog_blk(map_w3(1,cell_index(i)) + k-1) = real(du_dt_blk(i,k)*timestep, r_def)
        dv_orog_blk(map_w3(1,cell_index(i)) + k-1) = real(dv_dt_blk(i,k)*timestep, r_def)

        du_orog_gwd(map_w3(1,cell_index(i)) + k-1) = real(du_dt_orog_gwd(i,k)*timestep, r_def)
        dv_orog_gwd(map_w3(1,cell_index(i)) + k-1) = real(dv_dt_orog_gwd(i,k)*timestep, r_def)

        dtemp_orog_blk(map_wth(1,cell_index(i)) + k)      = real(dtemp_dt_blk(i,k)*timestep, r_def)
        dtemp_orog_gwd(map_wth(1,cell_index(i)) + k) = real(dtemp_dt_orog_gwd(i,k)*timestep, r_def)
      end do ! k
    end do !i
    
    ! Set level 0 increment such that theta increment will equal level 1
    do i = 1, seg_len
      dtemp_orog_blk(map_wth(1,cell_index(i)) + 0) = real(dtemp_dt_blk(i,1)*timestep, r_def)   &
                                   * exner_in_wth(map_wth(1,cell_index(i)) + 0)              &
                                   / exner_in_wth(map_wth(1,cell_index(i)) + 1)
      dtemp_orog_gwd(map_wth(1,cell_index(i)) + 0) = real(dtemp_dt_orog_gwd(i,1)*timestep, r_def)&
                                   * exner_in_wth(map_wth(1,cell_index(i)) + 0)                &
                                   / exner_in_wth(map_wth(1,cell_index(i)) + 1)

      ! Map diagnostics back
      if (.not. associated(taux_orog_blk, empty_real_data) ) then
        do k = 0, nlayers
          taux_orog_blk(map_wth(1,cell_index(i)) + k) = real(tau_x_blk(i,k), r_def)
        end do
      end if
      if (.not. associated(tauy_orog_blk, empty_real_data) ) then
        do k = 0, nlayers
          tauy_orog_blk(map_wth(1,cell_index(i)) + k) = real(tau_y_blk(i,k), r_def)
        end do
      end if
      if (.not. associated(taux_orog_gwd, empty_real_data) ) then
        do k = 0, nlayers
          taux_orog_gwd(map_wth(1,cell_index(i)) + k) = real(tau_x_orog_gwd(i,k), r_def)
        end do
      end if
      if (.not. associated(tauy_orog_gwd, empty_real_data) ) then
        do k = 0, nlayers
          tauy_orog_gwd(map_wth(1,cell_index(i)) + k) = real(tau_y_orog_gwd(i,k), r_def)
        end do
      end if
    end do !i
    
  end subroutine orographic_drag_kernel_code

end module orographic_drag_kernel_mod
