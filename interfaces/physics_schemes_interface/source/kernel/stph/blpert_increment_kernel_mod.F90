!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Compute increment of boundary layer perturbation.
!>
module blpert_increment_kernel_mod

  use argument_mod,           only: arg_type,                  &
                                    GH_FIELD, GH_SCALAR,       &
                                    ANY_DISCONTINUOUS_SPACE_1, &
                                    GH_INTEGER, GH_REAL,       &
                                    GH_READ, GH_WRITE,         &
                                    GH_READWRITE,              &
                                    CELL_COLUMN
  use fs_continuity_mod,      only: Wtheta
  use constants_mod,          only: i_def, r_def
  use kernel_mod,             only: kernel_type
  use stochastic_physics_config_mod, only: &
                                           ! Switches to use different
                                           ! parametrizations
                                           blpert_type, &
                                           blpert_type_theta_mag, &
                                           blpert_type_theta_star, &
                                           blpert_type_theta_and_moist, &
                                           ! peturbation detail
                                           blpert_add_vertical_shape, &
                                           blpert_max_magnitude

  implicit none
  private

  type, public, extends(kernel_type) :: blpert_increment_kernel_type
    private
    type(arg_type) :: meta_args(14) = (/                                          &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     WTHETA),                   & ! dtheta_blpert
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     WTHETA),                   & ! dmv_blpert
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   & ! mv
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      WTHETA),                   & ! height_wth
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),& ! blpert_sw
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),& ! ntml
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),& ! dtheta_base
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),& ! dmv_base
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),& ! theta_star_surf
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),& ! qv_star_surf
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),& ! blpert_rand_fld
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ                                ),& ! pert_lev_bot
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ                                ),& ! pert_lev_top
         arg_type(GH_SCALAR, GH_REAL,    GH_READ                                ) & ! tfac
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: blpert_increment_code
  end type

  public :: blpert_increment_code

contains

  !> @brief Compute increment of boundary layer perturbation.
  !> @param[in]    nlayers          Number of layers
  !> @param[out]   dtheta_blpert    Increments for theta
  !> @param[out]   dmv_blpert       Increments for mv
  !> @param[in]    mv               Moisture field
  !> @param[in]    height_wth       Height of wth levels above surface
  !> @param[in]    blpert_sw        Swich to boot perturbation (0:off, 1:on)
  !> @param[in]    ntml             Number of turbulently mixed levels
  !> @param[out]   dtheta_base      2d base field of theta perturbation
  !> @param[out]   dmv_base         2d base field of mv perturbation
  !> @param[in]    theta_star_surf  Atmospheric stability via surface heat flux
  !> @param[in]    qv_star_surf     Atmospheric stability via surface moist flux
  !> @param[in]    blpert_rand_fld  Latest random number field (-1.0 ~ 1.0)
  !> @param[in]    pert_lev_bot     Number of lowest layer to be perturbed
  !> @param[in]    pert_lev_top     Number of highest layer to be perturbed
  !> @param[in]    tfac             Adjustment factor depending on timestep
  !> @param[in]    ndf_wth          Number of DOFs per cell for Wtheta fields
  !> @param[in]    undf_wth         Number of unique DOFs for Wtheta fields
  !> @param[in]    map_wth          Dofmap for Wtheta fields
  !> @param[in]    ndf_2d           Number of DOFs per cell for 2D fields
  !> @param[in]    undf_2d          Number of unique DOFs for 2D fields
  !> @param[in]    map_2d           Dofmap for 2D fields
  subroutine blpert_increment_code(nlayers,                           &
                                   dtheta_blpert,                     &
                                   dmv_blpert,                        &
                                   mv,                                &
                                   height_wth,                        &
                                   blpert_sw,                         &
                                   ntml,                              &
                                   dtheta_base,                       &
                                   dmv_base,                          &
                                   theta_star_surf,                   &
                                   qv_star_surf,                      &
                                   blpert_rand_fld,                   &
                                   pert_lev_bot,                      &
                                   pert_lev_top,                      &
                                   tfac,                              &
                                   ndf_wth,                           &
                                   undf_wth,                          &
                                   map_wth,                           &
                                   ndf_2d,                            &
                                   undf_2d,                           &
                                   map_2d)
    implicit none

    integer(i_def), intent(in) :: nlayers
    integer(i_def), intent(in) :: ndf_wth, undf_wth
    integer(i_def), intent(in) :: ndf_2d, undf_2d
    integer(i_def), intent(in) :: map_wth(ndf_wth)
    integer(i_def), intent(in) :: map_2d(ndf_2d)
    real(r_def), intent(out) :: dtheta_blpert(ndf_wth)
    real(r_def), intent(out) :: dmv_blpert(ndf_wth)
    real(r_def), intent(in) :: mv(ndf_wth)
    real(r_def), intent(in) :: height_wth(ndf_wth)
    integer(i_def), intent(in) :: blpert_sw(ndf_2d)
    integer(i_def), intent(in) :: ntml(ndf_2d)
    real(r_def), intent(out) :: dtheta_base(ndf_2d)
    real(r_def), intent(out) :: dmv_base(ndf_2d)
    real(r_def), intent(in) :: theta_star_surf(ndf_2d)
    real(r_def), intent(in) :: qv_star_surf(ndf_2d)
    real(r_def), intent(in) :: blpert_rand_fld(ndf_2d)
    integer(i_def), intent(in) :: pert_lev_bot
    integer(i_def), intent(in) :: pert_lev_top
    real(r_def), intent(in) :: tfac
    integer(i_def) :: k, kmax_pert
    real(r_def) :: blpert_max_mv_mag
    real(r_def) :: zfac, zmix

    if (blpert_sw(map_2d(1)) /= 1_i_def) then
      dtheta_base(map_2d(1)) = 0.0_r_def
      dmv_base(map_2d(1)) = 0.0_r_def
      do k = 0, nlayers
        dtheta_blpert(map_wth(1)+k) = 0.0_r_def
        dmv_blpert(map_wth(1)+k) = 0.0_r_def
      end do
    else
      !----------------------------------------------------------------!
      ! Compute perturbation base
      !----------------------------------------------------------------!
      select case(blpert_type)
        case (blpert_type_theta_mag)
          dtheta_base(map_2d(1)) = blpert_rand_fld(map_2d(1)) * &
              tfac * blpert_max_magnitude
          dmv_base(map_2d(1)) = 0.0_r_def
        case(blpert_type_theta_star)
          dtheta_base(map_2d(1)) = blpert_rand_fld(map_2d(1)) * &
              min( tfac * theta_star_surf(map_2d(1)), blpert_max_magnitude )
          dmv_base(map_2d(1)) = 0.0_r_def
        case (blpert_type_theta_and_moist)
          dtheta_base(map_2d(1)) = blpert_rand_fld(map_2d(1)) * &
              min( tfac * theta_star_surf(map_2d(1)), blpert_max_magnitude )
          blpert_max_mv_mag = 0.1_r_def * mv(map_wth(1))
          dmv_base(map_2d(1)) = blpert_rand_fld(map_2d(1)) * &
              min( tfac * qv_star_surf(map_2d(1)), blpert_max_mv_mag )
      end select

      !----------------------------------------------------------------!
      ! Compute perturbation profile
      !----------------------------------------------------------------!
      do k = 0, pert_lev_bot - 1
        dtheta_blpert(map_wth(1)+k) = 0.0_r_def
        dmv_blpert(map_wth(1)+k) = 0.0_r_def
      end do

      select case(blpert_type)
        case (blpert_type_theta_mag)
          kmax_pert = min( pert_lev_top, 2*ntml(map_2d(1))/3 )
          do k = pert_lev_bot, kmax_pert
            dtheta_blpert(map_wth(1)+k) = dtheta_base(map_2d(1))
            dmv_blpert(map_wth(1)+k) = 0.0_r_def
          end do
        case(blpert_type_theta_star)
          kmax_pert = min( pert_lev_top, 2*ntml(map_2d(1))/3 )
          do k = pert_lev_bot, kmax_pert
            dtheta_blpert(map_wth(1)+k) = dtheta_base(map_2d(1))
            dmv_blpert(map_wth(1)+k) = 0.0_r_def
          end do
        case (blpert_type_theta_and_moist)
          if (blpert_add_vertical_shape) then
            kmax_pert = min( pert_lev_top, ntml(map_2d(1)) )
            zmix = height_wth(map_wth(1)+ntml(map_2d(1))+1)
            do k = pert_lev_bot, kmax_pert
              zfac = min( 2.0_r_def * height_wth(map_wth(1)+k)/zmix, &
                  2.0_r_def * (1.0_r_def - height_wth(map_wth(1)+k)/zmix) )
              dtheta_blpert(map_wth(1)+k) = dtheta_base(map_2d(1)) * zfac
              dmv_blpert(map_wth(1)+k) = dmv_base(map_2d(1)) * zfac
            end do
          else
            kmax_pert = min( pert_lev_top, 2*ntml(map_2d(1))/3 )
            do k = pert_lev_bot, kmax_pert
              dtheta_blpert(map_wth(1)+k) = dtheta_base(map_2d(1))
              dmv_blpert(map_wth(1)+k) = dmv_base(map_2d(1))
            end do
          end if
      end select ! blpert_type

      do k = kmax_pert + 1, nlayers
        dtheta_blpert(map_wth(1)+k) = 0.0_r_def
        dmv_blpert(map_wth(1)+k) = 0.0_r_def
      end do

    end if ! blpert_sw

  end subroutine blpert_increment_code
end module blpert_increment_kernel_mod
