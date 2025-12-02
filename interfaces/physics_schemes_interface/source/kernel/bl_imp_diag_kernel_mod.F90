!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculate boundary layer diagnostics.
!>
module bl_imp_diag_kernel_mod

  use argument_mod,           only: arg_type,                  &
                                    GH_FIELD,                  &
                                    ANY_DISCONTINUOUS_SPACE_1, &
                                    GH_REAL,                   &
                                    GH_READ, GH_WRITE,         &
                                    CELL_COLUMN
  use fs_continuity_mod,      only: W3, WTHETA
  use constants_mod,          only: i_def, r_def
  use kernel_mod,             only: kernel_type
  use planet_config_mod,      only: gravity, cp

  implicit none
  private

  type, public, extends(kernel_type) :: bl_imp_diag_kernel_type
    private
    type(arg_type) :: meta_args(11) = (/                                          &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   & ! bt_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   & ! bq_bl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       & ! ftl
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       & ! fqw
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       & ! taux
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       & ! tauy
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W3),                       & ! wetrho_in_w3
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),& ! theta_star_surf
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),& ! qv_star_surf
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),& ! ustar_implicit
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1) & ! zh
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: bl_imp_diag_code
  end type

  public :: bl_imp_diag_code

contains

  !> @brief Calculate boundary layer diagnostics.
  !> @param[in]    nlayers          Number of layers
  !> @param[in]    bt_bl            Buoyancy parameter for heat
  !> @param[in]    bq_bl            Buoyancy parameter for moisture
  !> @param[in]    ftl              Vertical heat flux on BL levels
  !> @param[in]    fqw              Vertical moisture flux on BL levels
  !> @param[in]    taux             'Zonal' momentum stress
  !> @param[in]    tauy             'Meridional' momentum stress
  !> @param[in]    wetrho_in_w3     Wet density field in w3 space
  !> @param[out]   theta_star_surf  Atmospheric stability via surface heat flux
  !> @param[out]   qv_star_surf     Atmospheric stability via surface moist flux
  !> @param[out]   ustar_implicit   Implicit surface friction velocity
  !> @param[in]    zh               Boundary layer depth
  !> @param[in]    ndf_wth          Number of DOFs per cell for Wtheta fields
  !> @param[in]    undf_wth         Number of unique DOFs for Wtheta fields
  !> @param[in]    map_wth          Dofmap for Wtheta fields
  !> @param[in]    ndf_w3           Number of DOFs per cell for W3 fields
  !> @param[in]    undf_w3          Number of unique DOFs for W3 fields
  !> @param[in]    map_w3           Dofmap for Wtheta fields
  !> @param[in]    ndf_2d           Number of DOFs per cell for 2D fields
  !> @param[in]    undf_2d          Number of unique DOFs for 2D fields
  !> @param[in]    map_2d           Dofmap for 2D fields
  subroutine bl_imp_diag_code(nlayers,                              &
                              bt_bl,                                &
                              bq_bl,                                &
                              ftl,                                  &
                              fqw,                                  &
                              taux,                                 &
                              tauy,                                 &
                              wetrho_in_w3,                         &
                              theta_star_surf,                      &
                              qv_star_surf,                         &
                              ustar_implicit,                       &
                              zh,                                   &
                              ndf_wth,                              &
                              undf_wth,                             &
                              map_wth,                              &
                              ndf_w3,                               &
                              undf_w3,                              &
                              map_w3,                               &
                              ndf_2d,                               &
                              undf_2d,                              &
                              map_2d)
    implicit none

    integer(i_def), intent(in) :: nlayers
    integer(i_def), intent(in) :: ndf_wth, undf_wth
    integer(i_def), intent(in) :: ndf_w3, undf_w3
    integer(i_def), intent(in) :: ndf_2d, undf_2d
    integer(i_def), intent(in) :: map_wth(ndf_wth)
    integer(i_def), intent(in) :: map_w3(ndf_w3)
    integer(i_def), intent(in) :: map_2d(ndf_2d)
    real(r_def), intent(in) :: bt_bl(ndf_wth)
    real(r_def), intent(in) :: bq_bl(ndf_wth)
    real(r_def), intent(in) :: ftl(ndf_w3)
    real(r_def), intent(in) :: fqw(ndf_w3)
    real(r_def), intent(in) :: taux(ndf_w3)
    real(r_def), intent(in) :: tauy(ndf_w3)
    real(r_def), intent(in) :: wetrho_in_w3(ndf_w3)
    real(r_def), intent(out) :: theta_star_surf(ndf_2d)
    real(r_def), intent(out) :: qv_star_surf(ndf_2d)
    real(r_def), intent(out) :: ustar_implicit(ndf_2d)
    real(r_def), intent(in) :: zh(ndf_2d)
    real(r_def) :: taux_surf, tauy_surf, wm
    real(r_def) :: ftl_surf, fqw_surf, fb_surf
    real(r_def), parameter :: one_quarter = 1.0_r_def/4.0_r_def
    real(r_def), parameter :: one_third = 1.0_r_def/3.0_r_def
    real(r_def), parameter :: c_ws = 0.25_r_def

    taux_surf = taux(map_w3(1)) / wetrho_in_w3(map_w3(1))
    tauy_surf = tauy(map_w3(1)) / wetrho_in_w3(map_w3(1))
    ustar_implicit(map_2d(1)) = ( taux_surf*taux_surf + &
                                  tauy_surf*tauy_surf )**one_quarter

    ftl_surf = ftl(map_w3(1)) / cp
    fqw_surf = fqw(map_w3(1))
    fb_surf = gravity * ( bt_bl(map_wth(1)) * ftl_surf +   &
                          bq_bl(map_wth(1)) * fqw_surf ) / &
                          wetrho_in_w3(map_w3(1))
    if ( ftl_surf > 0.0_r_def .and. fb_surf > 0.0_r_def ) then
      wm = ( ustar_implicit(map_2d(1))**3 + &
             c_ws * zh(map_2d(1)) * fb_surf )**one_third
      ! Eg: ftl=15 W/m2 and zh=500m gives wm~0.4 and
      !     not interested in very weakly forced convection
      wm = max( 0.4_r_def, wm )
      theta_star_surf(map_2d(1)) = ftl_surf / (wm * wetrho_in_w3(map_w3(1)))
      qv_star_surf(map_2d(1)) = max( 0.0_r_def, fqw_surf ) / &
                                   (wm * wetrho_in_w3(map_w3(1)))
    else
      theta_star_surf(map_2d(1)) = 0.0_r_def
      qv_star_surf(map_2d(1)) = 0.0_r_def
    end if

  end subroutine bl_imp_diag_code
end module bl_imp_diag_kernel_mod
