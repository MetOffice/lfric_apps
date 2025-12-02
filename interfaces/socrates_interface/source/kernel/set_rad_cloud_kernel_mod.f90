!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Set up cloud fields to be used in radiative transfer
module set_rad_cloud_kernel_mod

use argument_mod,      only : arg_type, &
                              GH_FIELD, GH_SCALAR, GH_REAL, &
                              GH_READ, GH_WRITE, &
                              CELL_COLUMN
use fs_continuity_mod, only : Wtheta
use constants_mod,     only : r_def, i_def
use kernel_mod,        only : kernel_type

implicit none

private

type, public, extends(kernel_type) :: set_rad_cloud_kernel_type
  private
  type(arg_type) :: meta_args(13) = (/ &
    arg_type(GH_FIELD,  GH_REAL,   GH_WRITE, Wtheta),                    & ! radiative_cloud_fraction
    arg_type(GH_FIELD,  GH_REAL,   GH_WRITE, Wtheta),                    & ! radiative_conv_fraction
    arg_type(GH_FIELD,  GH_REAL,   GH_WRITE, Wtheta),                    & ! conv_liquid_fraction
    arg_type(GH_FIELD,  GH_REAL,   GH_WRITE, Wtheta),                    & ! conv_frozen_fraction
    arg_type(GH_FIELD,  GH_REAL,   GH_WRITE, Wtheta),                    & ! conv_liquid_mmr
    arg_type(GH_FIELD,  GH_REAL,   GH_WRITE, Wtheta),                    & ! conv_frozen_mmr
    arg_type(GH_FIELD,  GH_REAL,   GH_READ,  Wtheta),                    & ! area_fraction
    arg_type(GH_FIELD,  GH_REAL,   GH_READ,  Wtheta),                    & ! cca
    arg_type(GH_FIELD,  GH_REAL,   GH_READ,  Wtheta),                    & ! ccw
    arg_type(GH_FIELD,  GH_REAL,   GH_READ,  Wtheta),                    & ! temperature_in_wth
    arg_type(GH_SCALAR, GH_REAL,   GH_READ         ),                    & ! allicetdegc
    arg_type(GH_SCALAR, GH_REAL,   GH_READ         ),                    & ! starticetkelvin
    arg_type(GH_SCALAR, GH_REAL,   GH_READ         )                     & ! zerodegc
    /)

  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: set_rad_cloud_code

end type set_rad_cloud_kernel_type

public :: set_rad_cloud_code

contains
!> @param[in]    nlayers                  Number of layers in Wtheta field
!> @param[inout] radiative_cloud_fraction Large scale cloud fraction
!> @param[inout] radiative_conv_fraction  Convective cloud fraction
!> @param[inout] conv_liquid_fraction     Convective liquid cloud fraction
!> @param[inout] conv_frozen_fraction     Convective frozen cloud fraction
!> @param[inout] conv_liquid_mmr          Convective liquid gridbox MMR
!> @param[inout] conv_frozen_mmr          Convective frozen gridbox MMR
!> @param[in]    area_fraction            Total cloud area fraction field
!> @param[in]    cca                      Convective cloud amount (fraction)
!> @param[in]    ccw                      Convective cloud water (kg/kg) (can be ice or liquid)
!> @param[in]    temperature_in_wth       Temperature in Wtheta space
!> @param[in]    allicetdegc
!> @param[in]    starticetkelvin
!> @param[in]    zerodegc
!> @param[in]    ndf_wtheta               No. DOFs per cell for Wtheta space
!> @param[in]    undf_wtheta              No. unique DOFs for Wtheta space
!> @param[in]    map_wtheta               Dofmap for Wtheta space column base cell
subroutine set_rad_cloud_code( nlayers, &
                               radiative_cloud_fraction, &
                               radiative_conv_fraction, &
                               conv_liquid_fraction, &
                               conv_frozen_fraction, &
                               conv_liquid_mmr, &
                               conv_frozen_mmr, &
                               area_fraction, &
                               cca, &
                               ccw, &
                               temperature_in_wth, &
                               allicetdegc, &
                               starticetkelvin, &
                               zerodegc, &
                               ndf_wtheta, &
                               undf_wtheta, &
                               map_wtheta )

  implicit none

  ! Dummy arguments
  integer(i_def), intent(in) :: nlayers
  integer(i_def), intent(in) :: ndf_wtheta
  integer(i_def), intent(in), dimension(ndf_wtheta) :: map_wtheta
  integer(i_def), intent(in) :: undf_wtheta
  real(r_def), intent(inout), dimension(undf_wtheta) :: &
    radiative_cloud_fraction, radiative_conv_fraction, &
    conv_liquid_fraction, conv_frozen_fraction, &
    conv_liquid_mmr, conv_frozen_mmr
  real(r_def), intent(in), dimension(undf_wtheta) :: &
    area_fraction, cca, ccw, temperature_in_wth
  real(r_def), intent(in) :: allicetdegc, starticetkelvin, zerodegc

  ! Internal constants and variables
  integer(i_def) :: j, k

  k = map_wtheta(1)
  radiative_cloud_fraction(k) = 0.0_r_def
  radiative_conv_fraction(k) = 0.0_r_def
  conv_liquid_fraction(k) = 0.0_r_def
  conv_frozen_fraction(k) = 0.0_r_def
  conv_liquid_mmr(k) = 0.0_r_def
  conv_frozen_mmr(k) = 0.0_r_def

  do j=1, nlayers
    k = map_wtheta(1) + j
    ! Squeeze large scale cloud into area not filled by convective cloud
    radiative_cloud_fraction(k) = area_fraction(k) * (1.0_r_def - cca(k))
    ! Convective cloud amount unaltered
    radiative_conv_fraction(k) = cca(k)
    ! Set liquid fraction of convective cloud using layer temperature
    conv_liquid_fraction(k) &
      = 1.0_r_def - ( (temperature_in_wth(k) - starticeTKelvin) &
                    / (alliceTdegC - (starticeTKelvin - zerodegc)) )
    conv_liquid_fraction(k) &
      = MIN(MAX(0.0_r_def, conv_liquid_fraction(k)), 1.0_r_def)
    conv_frozen_fraction(k) = 1.0_r_def - conv_liquid_fraction(k)
    ! Convert liquid and ice fractions to gridbox fractions
    conv_liquid_fraction(k) &
      = conv_liquid_fraction(k) * radiative_conv_fraction(k)
    conv_frozen_fraction(k) &
      = conv_frozen_fraction(k) * radiative_conv_fraction(k)
    ! Calculate conserved gridbox mean liquid and ice MMRs
    conv_liquid_mmr(k) = ccw(k) * conv_liquid_fraction(k)
    conv_frozen_mmr(k) = ccw(k) * conv_frozen_fraction(k)

  end do

end subroutine set_rad_cloud_code

end module set_rad_cloud_kernel_mod
