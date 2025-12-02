!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Update prognostic field.
!>
module blpert_update_rand_kernel_mod

  use argument_mod,           only: arg_type,                  &
                                    GH_FIELD, GH_SCALAR,       &
                                    ANY_DISCONTINUOUS_SPACE_1, &
                                    GH_INTEGER, GH_REAL,       &
                                    GH_READ, GH_WRITE,         &
                                    GH_READWRITE,              &
                                    CELL_COLUMN
  use constants_mod,          only: i_def, r_def, l_def
  use kernel_mod,             only: kernel_type
  use stochastic_physics_config_mod, only: blpert_type, &
                                           blpert_type_theta_star, &
                                           blpert_type_theta_and_moist, &
                                           blpert_noncumulus_points, &
                                           blpert_time_correlation

  implicit none
  private

  type, public, extends(kernel_type) :: blpert_update_rand_kernel_type
    private
    type(arg_type) :: meta_args(9) = (/                                           &
         arg_type(GH_FIELD,  GH_INTEGER, GH_WRITE,     ANY_DISCONTINUOUS_SPACE_1),& ! blpert_sw
         arg_type(GH_FIELD,  GH_INTEGER, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),& ! blpert_flag
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),& ! blpert_area
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_1),& ! cumulus
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),& ! blpert_rand_fld
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),& ! theta_star_surf
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1),& ! rand_numb
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ                                ),& ! outer
         arg_type(GH_SCALAR, GH_REAL,    GH_READ                                ) & ! auto_corr_coeff
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: blpert_update_rand_code
  end type

  public :: blpert_update_rand_code

contains

  !> @brief Update prognostic field.
  !> @param[in]     nlayers          Number of layers
  !> @param[out]    blpert_sw        Swich to perturbate field (0:off, 1:on)
  !> @param[in,out] blpert_flag      Flag to show whether this point has ever been perturbed more than once
  !> @param[in]     blpert_area      Flag of Perturbation target area
  !> @param[in]     cumulus          Cumulus flag
  !> @param[in,out] blpert_rand_fld  Latest random number field (-1.0 ~ 1.0)
  !> @param[in]     theta_star_surf  Atmospheric stability via surface heat flux
  !> @param[in]     rand_numb        Newly generated Random number (-1.0 ~ 1.0)
  !> @param[in]     outer            Outer loop counter
  !> @param[in]     auto_corr_coeff  Time decorrelation coefficient
  !> @param[in]     ndf_2d           Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d          Number of unique DOFs for 2D fields
  !> @param[in]     map_2d           Dofmap for 2D fields
  subroutine blpert_update_rand_code(nlayers,                           &
                                     blpert_sw,                         &
                                     blpert_flag,                       &
                                     blpert_area,                       &
                                     cumulus,                           &
                                     blpert_rand_fld,                   &
                                     theta_star_surf,                   &
                                     rand_numb,                         &
                                     outer,                             &
                                     auto_corr_coeff,                   &
                                     ndf_2d,                            &
                                     undf_2d,                           &
                                     map_2d)
    implicit none

    integer(i_def), intent(in) :: nlayers
    integer(i_def), intent(in) :: ndf_2d, undf_2d
    integer(i_def), intent(in) :: map_2d(ndf_2d)
    integer(i_def), intent(out) :: blpert_sw(ndf_2d)
    integer(i_def), intent(inout) :: blpert_flag(ndf_2d)
    integer(i_def), intent(in) :: blpert_area(ndf_2d)
    integer(i_def), intent(in) :: cumulus(ndf_2d)
    real(r_def), intent(inout) :: blpert_rand_fld(ndf_2d)
    real(r_def), intent(in) :: theta_star_surf(ndf_2d)
    real(r_def), intent(in) :: rand_numb(ndf_2d)
    integer(i_def), intent(in) :: outer
    real(r_def), intent(in) :: auto_corr_coeff
    logical(l_def) :: area_flag, cumulus_flag, sflux_flag
    real(r_def), parameter :: rand_min = -1.0_r_def
    real(r_def), parameter :: rand_max = 1.0_r_def
    real(r_def), parameter :: shock_amp = 1.0_r_def

    !----------------------------------------------------------------!
    ! Compute perturbaton flag
    !----------------------------------------------------------------!
    area_flag = blpert_area(map_2d(1)) == 1_i_def

    cumulus_flag = blpert_noncumulus_points .or. cumulus(map_2d(1)) == 1_i_def

    if (blpert_type == blpert_type_theta_star .or. &
        blpert_type == blpert_type_theta_and_moist) then
      sflux_flag = theta_star_surf(map_2d(1)) > 0.0_r_def
    else
      sflux_flag = .true.
    end if

    if (area_flag .and. cumulus_flag .and. sflux_flag) then
      blpert_sw(map_2d(1)) = 1_i_def
    else
      blpert_sw(map_2d(1)) = 0_i_def
    end if

    !----------------------------------------------------------------!
    ! Update random number field
    !----------------------------------------------------------------!
    if (blpert_sw(map_2d(1)) == 1_i_def) then
      if ( blpert_time_correlation .and. &
           blpert_flag(map_2d(1)) == 1_i_def ) then
        if ( outer == 1_i_def ) then
          ! Evolve perturbation with an AR1 process
          blpert_rand_fld(map_2d(1)) = &
              auto_corr_coeff * blpert_rand_fld(map_2d(1)) + &
              sqrt(1.0_r_def - auto_corr_coeff**2) * &
              rand_numb(map_2d(1)) * shock_amp
          ! Reflect perturbation at specified maximum and minimum
          if (blpert_rand_fld(map_2d(1)) > rand_max) then
            blpert_rand_fld(map_2d(1)) = 2.0_r_def * rand_max - &
                                         blpert_rand_fld(map_2d(1))
          end if
          if (blpert_rand_fld(map_2d(1)) < rand_min) then
            blpert_rand_fld(map_2d(1)) = 2.0_r_def * rand_min - &
                                         blpert_rand_fld(map_2d(1))
          end if
        end if
      else
        ! Random correlation or newly perturbed case
        blpert_rand_fld(map_2d(1)) = rand_numb(map_2d(1))
      end if
    end if

    !----------------------------------------------------------------!
    ! Update blpert_flag
    !----------------------------------------------------------------!
    ! Note that this flag means whether this point has ever been
    ! perturbed more than once so far.
    if (blpert_sw(map_2d(1)) == 1_i_def) then
      blpert_flag(map_2d(1)) = 1_i_def
    end if

  end subroutine blpert_update_rand_code
end module blpert_update_rand_kernel_mod
