!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Stratospheric aerosol climatology

module strat_aerosol_kernel_mod

use argument_mod,      only : arg_type,          &
                              GH_FIELD, GH_REAL, &
                              GH_READ, GH_WRITE, &
                              CELL_COLUMN,       &
                              GH_INTEGER,        &
                              ANY_DISCONTINUOUS_SPACE_1
use fs_continuity_mod, only:  Wtheta, W3
use constants_mod,     only : r_def, i_def
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
! The type declaration for the kernel.
! Contains the metadata needed by the PSy layer.
type, public, extends(kernel_type) :: strat_aerosol_kernel_type
  private
  type(arg_type) :: meta_args(3) = (/                                    &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, Wtheta),                    & ! sulphuric
       arg_type(GH_FIELD, GH_INTEGER,GH_READ,ANY_DISCONTINUOUS_SPACE_1), & ! trop_level
       arg_type(GH_FIELD, GH_REAL, GH_READ,  W3)                         & ! exner
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: strat_aerosol_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: strat_aerosol_code

contains

!> @param[in]     nlayers               Number of layers
!> @param[in,out] sulphuric             Sulphuric acid aerosol MMR
!> @param[in]     trop_level            Level of tropopause
!> @param[in]     exner                 Exner pressure on w3 space
!> @param[in]     ndf_wth               No. DOFs per cell for wth space
!> @param[in]     undf_wth              No. unique DOFs for wth space
!> @param[in]     map_wth               Dofmap for wth space column base cell
!> @param[in]     ndf_2d                No. DOFs per cell for 2D space
!> @param[in]     undf_2d               No. unique DOFs for 2D space
!> @param[in]     map_2d                Dofmap for 2D space column base cell
!> @param[in]     ndf_w3                No. DOFs per cell for w3 space
!> @param[in]     undf_w3               No. unique DOFs for w3 space
!> @param[in]     map_w3                Dofmap for w3 space column base cell
subroutine strat_aerosol_code(nlayers,                    &
                              sulphuric,                  &
                              trop_level,                 &
                              exner,                      &
                              ndf_wth, undf_wth, map_wth, &
                              ndf_2d, undf_2d, map_2d,    &
                              ndf_w3, undf_w3, map_w3)

  use planet_config_mod, only: p_zero, kappa, gravity
  use aerosol_config_mod, only: sulphuric_strat_column

  implicit none

  ! Arguments
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: ndf_wth, ndf_2d, ndf_w3
  integer(i_def), intent(in) :: undf_wth, undf_2d, undf_w3

  integer(i_def), dimension(ndf_wth), intent(in) :: map_wth
  integer(i_def), dimension(ndf_2d),  intent(in) :: map_2d
  integer(i_def), dimension(ndf_w3),  intent(in) :: map_w3

  real(r_def), dimension(undf_wth), intent(inout) :: sulphuric
  integer(i_def), dimension(undf_2d),  intent(in) :: trop_level
  real(r_def), dimension(undf_w3),  intent(in) :: exner

  ! Local variables
  integer(i_def) :: i_trop
  real(r_def) :: sulphuric_mmr


  ! Level of tropopause
  i_trop = trop_level(map_2d(1))

  ! Sulphuric mass mixing ratio (kg/kg) above tropopause using the ratio of
  ! the stratospheric column amount of sulphuric acid aerosol (kg m-2) to
  ! the hydrostatic mass of the atmosphere above the tropopause:
  sulphuric_mmr = sulphuric_strat_column * gravity &
                / ( p_zero * exner(map_w3(1)+i_trop-1)**(1.0_r_def/kappa) )

  ! MMR is zero in the troposphere
  sulphuric(map_wth(1):map_wth(1)+i_trop-1) = 0.0_r_def
  ! Constant MMR above the tropopause
  sulphuric(map_wth(1)+i_trop:map_wth(1)+nlayers) = sulphuric_mmr

end subroutine strat_aerosol_code

end module strat_aerosol_kernel_mod
