!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! Some content in this file was generated or refactored with assistance from
! - Met Office Claude Code via AI Foundry (Claude Opus 5.5).

!> @brief Column maximum of the orographic drag stress in the height band of
!>        the mountain wave turbulence predictor.

module mtn_wave_max_stress_kernel_mod

  use argument_mod,      only: arg_type,                  &
                               GH_FIELD, GH_REAL,         &
                               GH_READ, GH_WRITE,         &
                               CELL_COLUMN,               &
                               ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,     only: r_def, i_def
  use fs_continuity_mod, only: W3, Wtheta
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: mtn_wave_max_stress_kernel_type
    private
    type(arg_type) :: meta_args(7) = (/                                    &
         ! height_w3 comes first so that nlayers is that of the 3D mesh
         arg_type(GH_FIELD, GH_REAL, GH_READ,  W3),                        &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta),                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  Wtheta)                     &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: mtn_wave_max_stress_code
  end type mtn_wave_max_stress_kernel_type

  public :: mtn_wave_max_stress_code

  ! Empirical height band approximating levels 23 to 47 of a 70-level set
  ! with an 80 km lid. It is raised over high ground by fractions of the
  ! surface height, so that the highest ranges still give a signal.
  !> Height of the bottom of the band over a sea-level surface (m)
  real(kind=r_def), parameter :: band_bottom = 3767.0_r_def
  !> Height of the top of the band over a sea-level surface (m)
  real(kind=r_def), parameter :: band_top = 15248.0_r_def
  !> Fraction of the surface height added to the bottom of the band
  real(kind=r_def), parameter :: orog_bottom = 0.5_r_def
  !> Fraction of the surface height added to the top of the band
  real(kind=r_def), parameter :: orog_top = 1.5_r_def

contains

  !> @brief   Column maximum of the orographic drag stress in the height band
  !>          of the mountain wave turbulence predictor.
  !> @details The stress is the magnitude of the total, blocking plus
  !>          gravity wave, orographic drag stress. Level k is in the band
  !>          when the height of layer k's centre, less a fraction of the
  !>          surface height, lies strictly inside the band. That height is
  !>          half a layer below the stress on level k; the original method
  !>          pairs them this way and this kernel keeps it. Heights increase
  !>          with k, so the search stops at the first level above the band.
  !> @param[in]     nlayers        Number of layers
  !> @param[in]     height_w3      Height of the layer centres above mean
  !>                               sea level (m)
  !> @param[in,out] max_stress     Column maximum of the stress in the band
  !>                               (Pa)
  !> @param[in]     taux_orog_blk  Zonal stress from orographic blocking (Pa)
  !> @param[in]     tauy_orog_blk  Meridional stress from orographic blocking
  !>                               (Pa)
  !> @param[in]     taux_orog_gwd  Zonal stress from orographic gravity wave
  !>                               drag (Pa)
  !> @param[in]     tauy_orog_gwd  Meridional stress from orographic gravity
  !>                               wave drag (Pa)
  !> @param[in]     height_wth     Height of the theta levels above mean sea
  !>                               level (m)
  !> @param[in]     ndf_w3         Number of degrees of freedom per cell for W3
  !> @param[in]     undf_w3        Number of unique degrees of freedom for W3
  !> @param[in]     map_w3         Dofmap for the cell for W3
  !> @param[in]     ndf_2d         Number of degrees of freedom per cell for
  !>                               the 2D output
  !> @param[in]     undf_2d        Number of unique degrees of freedom for the
  !>                               2D output
  !> @param[in]     map_2d         Dofmap for the cell for the 2D output
  !> @param[in]     ndf_wth        Number of degrees of freedom per cell for
  !>                               Wtheta
  !> @param[in]     undf_wth       Number of unique degrees of freedom for
  !>                               Wtheta
  !> @param[in]     map_wth        Dofmap for the cell for Wtheta
  subroutine mtn_wave_max_stress_code(nlayers,                        &
                                      height_w3,                      &
                                      max_stress,                     &
                                      taux_orog_blk, tauy_orog_blk,   &
                                      taux_orog_gwd, tauy_orog_gwd,   &
                                      height_wth,                     &
                                      ndf_w3, undf_w3, map_w3,        &
                                      ndf_2d, undf_2d, map_2d,        &
                                      ndf_wth, undf_wth, map_wth)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w3, undf_w3
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
    integer(kind=i_def), intent(in), dimension(ndf_w3)  :: map_w3
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d
    integer(kind=i_def), intent(in), dimension(ndf_wth) :: map_wth

    real(kind=r_def), intent(in),    dimension(undf_w3)  :: height_w3
    real(kind=r_def), intent(inout), dimension(undf_2d)  :: max_stress
    real(kind=r_def), intent(in),    dimension(undf_wth) :: taux_orog_blk
    real(kind=r_def), intent(in),    dimension(undf_wth) :: tauy_orog_blk
    real(kind=r_def), intent(in),    dimension(undf_wth) :: taux_orog_gwd
    real(kind=r_def), intent(in),    dimension(undf_wth) :: tauy_orog_gwd
    real(kind=r_def), intent(in),    dimension(undf_wth) :: height_wth

    ! Internal variables
    integer(kind=i_def) :: k
    real(kind=r_def)    :: surface, height, stress, column_max

    surface = height_wth(map_wth(1))
    column_max = 0.0_r_def

    do k = 1, nlayers
      height = height_w3(map_w3(1) + k - 1)
      if ( height - orog_top * surface >= band_top ) exit
      if ( height - orog_bottom * surface <= band_bottom ) cycle
      stress = sqrt( ( taux_orog_blk(map_wth(1) + k)                         &
                     + taux_orog_gwd(map_wth(1) + k) )**2                    &
                   + ( tauy_orog_blk(map_wth(1) + k)                         &
                     + tauy_orog_gwd(map_wth(1) + k) )**2 )
      column_max = max(column_max, stress)
    end do

    max_stress(map_2d(1)) = column_max

  end subroutine mtn_wave_max_stress_code

end module mtn_wave_max_stress_kernel_mod
