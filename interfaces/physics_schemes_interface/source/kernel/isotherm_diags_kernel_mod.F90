!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some content in this file was generated or refactored with assistance from
! - Claude Code (Claude Opus 5.5), 2026-10-16.
!> @brief Height and pressure of an isotherm (dry-bulb temperature).

module isotherm_diags_kernel_mod

  use argument_mod,      only: arg_type,                  &
                               GH_FIELD, GH_SCALAR,       &
                               GH_REAL, GH_LOGICAL,       &
                               GH_READ, GH_WRITE,         &
                               CELL_COLUMN,               &
                               ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,     only: r_def, i_def, l_def, rmdi
  use empty_data_mod,    only: empty_real_data
  use fs_continuity_mod, only: Wtheta
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: isotherm_diags_kernel_type
    private
    ! theta comes first so that nlayers is that of Wtheta
    type(arg_type) :: meta_args(10) = (/                                  &
         arg_type(GH_FIELD, GH_REAL, GH_READ, Wtheta),                    &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ, Wtheta),                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ, Wtheta),                    &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                           &
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                        &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                           &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                           &
         arg_type(GH_SCALAR, GH_REAL, GH_READ)                            &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: isotherm_diags_code
  end type isotherm_diags_kernel_type

  public :: isotherm_diags_code

contains

  !> @brief Finds the height and pressure of the highest level below 23 km
  !>        at which the temperature falls through a reference value.
  !> @details Each output is only written if it was requested, i.e. it is
  !>          not empty_real_data. Points where no level is found are set
  !>          to rmdi. If the lowest level is already colder than the
  !>          reference, the isotherm is placed at the surface, or, when
  !>          extrapolate is true, below it using a fixed lapse rate.
  !> @param[in]     nlayers            Number of layers
  !> @param[in]     theta              Potential temperature (K)
  !> @param[in,out] isotherm_height    Height of the isotherm above mean sea
  !>                                   level (m)
  !> @param[in,out] isotherm_pressure  Pressure at the isotherm (Pa)
  !> @param[in]     exner              Exner pressure in Wtheta
  !> @param[in]     height_wth         Height of Wtheta levels above mean sea
  !>                                   level (m)
  !> @param[in]     t_ref              Temperature of the isotherm (K)
  !> @param[in]     extrapolate        Extrapolate below the surface if the
  !>                                   lowest level is colder than t_ref
  !> @param[in]     p_zero             Reference pressure for Exner (Pa)
  !> @param[in]     recip_kappa        Ratio of cp to R
  !> @param[in]     g_over_r           Gravity divided by the gas constant
  !> @param[in]     ndf_wth            Number of DOFs per cell for Wtheta
  !> @param[in]     undf_wth           Number of unique DOFs for Wtheta
  !> @param[in]     map_wth            Dofmap for the Wtheta column base
  !> @param[in]     ndf_2d             Number of DOFs per cell for 2D fields
  !> @param[in]     undf_2d            Number of unique DOFs for 2D fields
  !> @param[in]     map_2d             Dofmap for the 2D fields
  subroutine isotherm_diags_code( nlayers,              &
                                  theta,                &
                                  isotherm_height,      &
                                  isotherm_pressure,    &
                                  exner,                &
                                  height_wth,           &
                                  t_ref,                &
                                  extrapolate,          &
                                  p_zero,               &
                                  recip_kappa,          &
                                  g_over_r,             &
                                  ndf_wth, undf_wth, map_wth, &
                                  ndf_2d, undf_2d, map_2d )

    implicit none

    ! Arguments added automatically in call to kernel
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d) :: map_2d

    ! Arguments passed explicitly from algorithm
    real(kind=r_def), intent(in), dimension(undf_wth) :: theta
    real(kind=r_def), intent(inout), pointer :: isotherm_height(:)
    real(kind=r_def), intent(inout), pointer :: isotherm_pressure(:)
    real(kind=r_def), intent(in), dimension(undf_wth) :: exner
    real(kind=r_def), intent(in), dimension(undf_wth) :: height_wth
    real(kind=r_def),    intent(in) :: t_ref
    logical(kind=l_def), intent(in) :: extrapolate
    real(kind=r_def),    intent(in) :: p_zero
    real(kind=r_def),    intent(in) :: recip_kappa
    real(kind=r_def),    intent(in) :: g_over_r

    ! Lapse rate (K/m) and pressure gradient (Pa/m) assumed below the surface
    real(kind=r_def), parameter :: lapse_below_surface = 10.0e-3_r_def
    real(kind=r_def), parameter :: dpdz_below_surface = -10.0_r_def
    ! Only levels below this height (m) are searched. Above the stratosphere
    ! the temperature falls again, and the search keeps the highest crossing.
    real(kind=r_def), parameter :: height_cutoff = 23000.0_r_def

    ! Internal variables
    integer(kind=i_def) :: k, k_cross
    logical(kind=l_def) :: l_height, l_pressure
    real(kind=r_def)    :: t_lower, t_upper, dz, lapse
    real(kind=r_def)    :: z_iso, p_iso

    l_height   = .not. associated( isotherm_height, empty_real_data )
    l_pressure = .not. associated( isotherm_pressure, empty_real_data )

    z_iso = rmdi
    p_iso = rmdi

    ! Lowest level already colder than the isotherm
    t_lower = theta(map_wth(1) + 1) * exner(map_wth(1) + 1)
    if ( t_lower < t_ref ) then
      if ( extrapolate ) then
        dz = ( t_lower - t_ref ) / lapse_below_surface
        z_iso = max( height_wth(map_wth(1)) + dz, 0.0_r_def )
        p_iso = p_zero * exner(map_wth(1))**recip_kappa                    &
              + dpdz_below_surface * dz
      else
        z_iso = height_wth(map_wth(1))
        p_iso = p_zero * exner(map_wth(1))**recip_kappa
      end if
    end if

    ! Highest level below the cut-off where the temperature falls through
    ! the isotherm. Heights increase up the column, so the search can stop
    ! at the first level at or above the cut-off.
    k_cross = 0
    do k = 2, nlayers
      if ( height_wth(map_wth(1) + k) >= height_cutoff ) exit
      t_upper = theta(map_wth(1) + k) * exner(map_wth(1) + k)
      if ( t_upper < t_ref .and. t_lower > t_ref ) k_cross = k
      t_lower = t_upper
    end do

    ! Interpolate within the crossing layer assuming a constant lapse rate,
    ! which is positive because of the strict tests above
    if ( k_cross > 0 ) then
      k = k_cross
      t_lower = theta(map_wth(1) + k - 1) * exner(map_wth(1) + k - 1)
      t_upper = theta(map_wth(1) + k) * exner(map_wth(1) + k)
      lapse = ( t_lower - t_upper )                                        &
            / ( height_wth(map_wth(1) + k) - height_wth(map_wth(1) + k - 1) )
      z_iso = height_wth(map_wth(1) + k - 1) + ( t_lower - t_ref ) / lapse
      p_iso = p_zero * exner(map_wth(1) + k - 1)**recip_kappa              &
            * ( t_ref / t_lower )**( g_over_r / lapse )
    end if

    if ( l_height ) isotherm_height(map_2d(1)) = z_iso
    if ( l_pressure ) isotherm_pressure(map_2d(1)) = p_iso

  end subroutine isotherm_diags_code

end module isotherm_diags_kernel_mod
