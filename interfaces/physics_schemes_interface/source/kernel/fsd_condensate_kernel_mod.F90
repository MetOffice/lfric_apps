!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates fractional standard deviation (FSD) of condensate

module fsd_condensate_kernel_mod

use argument_mod,         only: arg_type,              &
                                GH_FIELD, GH_REAL,     &
                                GH_READ, GH_READWRITE, &
                                GH_WRITE, CELL_COLUMN, &
                                ANY_DISCONTINUOUS_SPACE_1
use extrusion_config_mod, only: planet_radius
use fs_continuity_mod,    only: WTHETA
use kernel_mod,           only: kernel_type
use cloud_config_mod,     only: use_fsd_eff_res, fsd_min_conv_frac, &
     fsd_conv_const, fsd_nonconv_liq_const, fsd_nonconv_ice_const
use fsd_parameters_mod,   only: fsd_eff_lam, f_cons

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer

type, public, extends(kernel_type) :: fsd_condensate_kernel_type
  private
  type(arg_type) :: meta_args(7) = (/                                   &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, WTHETA),                   & ! sigma_mc
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),& ! f_arr_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                   & ! acf_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                   & ! cff_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                   & ! cca_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA),                   & ! dz_wth
       arg_type(GH_FIELD, GH_REAL, GH_READ,  WTHETA)                    & ! delta
       /)
   integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: fsd_condensate_code
end type

public :: fsd_condensate_code

contains

!> @brief Calculates fractional standard deviation (FSD) of condensate
!> @details Uses cloud fraction, grid box size and convective cloudiness
!>          to determine variability in condensate for use by
!>          radiation and microphysics scheme.
!> @param[in]     nlayers    Number of layers
!> @param[in,out] sigma_mc   Fractional standard deviation of condensate
!> @param[in,out] f_arr_wth  Parameters related to fractional standard deviation of condensate
!> @param[in]     acf_wth    Area cloud fraction
!> @param[in]     cff_wth    Frozen cloud fraction
!> @param[in]     cca_wth    Convective cloud amount
!> @param[in]     dz_wth     Delta z at wtheta levels
!> @param[in]     delta      Edge length on wtheta points
!> @param[in]     ndf_wth    Number of degrees of freedom per cell for potential temperature space
!> @param[in]     undf_wth   Number unique of degrees of freedom  for potential temperature space
!> @param[in]     map_wth    Dofmap for the cell at the base of the column for potential temperature space
!> @param[in]     ndf_farr   Number of degrees of freedom per cell for fsd array
!> @param[in]     undf_farr  Number unique of degrees of freedom for fsd array
!> @param[in]     map_farr   Dofmap for the cell at the base of the column for fsd array

subroutine fsd_condensate_code( nlayers,                    &
                                           ! InOut
                                sigma_mc_wth,               &
                                f_arr_wth,                  &
                                           ! In
                                acf_wth,                    &
                                cff_wth,                    &
                                cca_wth,                    &
                                dz_wth,                     &
                                delta,                      &
                                           ! DoF and Map info
                                ndf_wth, undf_wth, map_wth, &
                                ndf_farr, undf_farr, map_farr)

    use constants_mod, only: r_def, i_def, r_um, i_um

    !---------------------------------------
    ! UM modules
    !---------------------------------------

    use rad_input_mod, only: two_d_fsd_factor ! Condensate variability is
    ! parametrized using an empirical function fitting it to one-dimensional
    ! data. The two_d_fsd_factor rescales this value to make it representative
    ! of two-dimensional varibility. See Hill et al (2015) DOI: 10.1002/qj.2506

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth,  undf_wth
    integer(kind=i_def), intent(in) :: ndf_farr, undf_farr

    integer(kind=i_def), intent(in), dimension(ndf_wth)   :: map_wth
    integer(kind=i_def), intent(in), dimension(ndf_farr)  :: map_farr

    real(kind=r_def), intent(in),    dimension(undf_wth)  :: acf_wth
    real(kind=r_def), intent(in),    dimension(undf_wth)  :: cff_wth
    real(kind=r_def), intent(in),    dimension(undf_wth)  :: cca_wth
    real(kind=r_def), intent(in),    dimension(undf_wth)  :: dz_wth
    real(kind=r_def), intent(in),    dimension(undf_wth)  :: delta
    real(kind=r_def), intent(inout), dimension(undf_wth)  :: sigma_mc_wth
    real(kind=r_def), intent(inout), dimension(undf_farr) :: f_arr_wth

    real(kind=r_def), dimension(  nlayers) :: x_in_km ! Grid box size
    real(kind=r_def), dimension(3,nlayers) :: f_arr   ! FSD parameters
    real(kind=r_def)                       :: conv_thick_part ! Convective contrib
    real(kind=r_def)                       :: cloud_scale     ! Size of cld
    real(kind=r_def)                       :: col_max_cca     ! Maximum CCA
    real(kind=r_def), parameter            :: one_third       = 1.0_r_def/3.0_r_def
    real(kind=r_def), parameter            :: m_to_km         = 0.001_r_def

    integer(i_um) :: k, n

    if (use_fsd_eff_res) then
      ! Use a fixed effective resolution in FSD parametrization.
      do k = 1, nlayers
        x_in_km(k) = fsd_eff_lam * planet_radius * m_to_km
      end do
    else
      ! Use edge length as measure of horizontal distance in FSD parametrization
      ! (and convert from m to km).
      do k = 1, nlayers
        x_in_km(k) = delta(map_wth(1) + k) * m_to_km
      end do
    end if

    ! Compute column max convective cloud amount
    col_max_cca = cca_wth(map_wth(1) + 1)
    do k = 2, nlayers
      col_max_cca = max(col_max_cca, cca_wth(map_wth(1) + k))
    end do

    ! Equations taken from Hill et al (2015, DOI: 10.1002/qj.2506)
    ! modified to include different FSD in non-convective liquid and ice cloud
    do k = 1, nlayers

      if (cca_wth(map_wth(1) + k) > fsd_min_conv_frac) then
        conv_thick_part = fsd_conv_const * (x_in_km(k)**(-0.12_r_def))         &
                        * ( dz_wth(map_wth(1)+k) * m_to_km ) ** 0.07_r_def
      else if (cff_wth(map_wth(1)+k) > 0.0_r_def .and.                         &
               col_max_cca < 0.02_r_def) then
        conv_thick_part = fsd_nonconv_ice_const * (x_in_km(k)**0.002_r_def)    &
                        * ( dz_wth(map_wth(1)+k) * m_to_km ) ** 0.12_r_def
      else
        conv_thick_part = fsd_nonconv_liq_const * (x_in_km(k)**0.002_r_def)    &
                        * ( dz_wth(map_wth(1)+k) * m_to_km ) ** 0.12_r_def
      end if

      ! There are 3 parameters in the FSD array, so define all 3
      f_arr(1,k) = 0.12_r_def * conv_thick_part
      f_arr(2,k) = 0.23_r_def * conv_thick_part
      f_arr(3,k) = 0.05_r_def * conv_thick_part

      if (acf_wth(map_wth(1)+k) < 1.0_r_def) then
        cloud_scale = acf_wth(map_wth(1)+k) * x_in_km(k)
        sigma_mc_wth(map_wth(1)+k) = two_d_fsd_factor *                        &
                     (f_arr(2,k)-(f_arr(3,k)*acf_wth(map_wth(1)+k)))           &
                     * (cloud_scale ** one_third)                              &
                     * ((((f_cons(1) * cloud_scale) ** f_cons(2)) + 1.0_r_def) &
                     ** (f_cons(3)))
      else
        sigma_mc_wth(map_wth(1)+k) = two_d_fsd_factor * f_arr(1,k)             &
                     * (x_in_km(k) ** one_third)                               &
                     * ((((f_cons(1) * x_in_km(k)) ** f_cons(2)) + 1.0_r_def)  &
                     ** (f_cons(3)))
      end if

    end do ! k

    ! Recast into LFRic space
    do n = 1, 3
      do k = 1, nlayers
        ! There are 3 parameters in the FSD array, so loop over all 3
        f_arr_wth(map_farr(1) + (n-1)*(nlayers+1) + k) = f_arr(n,k)
      end do
    end do

    ! Deal with bottom level
    sigma_mc_wth(map_wth(1)+0) = sigma_mc_wth(map_wth(1)+1)
    do n = 1, 3
      ! There are 3 parameters in the FSD array, so loop over all 3
      f_arr_wth(map_farr(1) + (n-1)*(nlayers+1) + 0) = &
      f_arr_wth(map_farr(1) + (n-1)*(nlayers+1) + 1)
    end do

end subroutine fsd_condensate_code

end module fsd_condensate_kernel_mod
