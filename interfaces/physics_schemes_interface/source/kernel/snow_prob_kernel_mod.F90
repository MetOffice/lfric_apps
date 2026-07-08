!-------------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Interface to snow probability (Boyden method) calculation.

module snow_prob_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD, GH_READ, GH_WRITE,      &
                               GH_REAL, CELL_COLUMN,                       &
                               ANY_DISCONTINUOUS_SPACE_1,                  &
                               ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,     only: r_def, i_def
  use kernel_mod,        only: kernel_type

  implicit none
  private

  !> Kernel metadata for PSyclone
  type, public, extends(kernel_type) :: snow_prob_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: snow_prob_code
  end type snow_prob_kernel_type

  public :: snow_prob_code

contains

  !> @details Snow probability (%) by the Boyden method:
  !!            P = clamp( 5220 + 3.86666*Z1000 - 4*Z850 , 0, 100 )
  !!          with Z1000/Z850 the geopotential height (m) of the 1000/850 hPa
  !!          surfaces. Constants are dimensional in {%, m} and calibrated
  !!          against percent output; do not rescale (Boyden 1964,
  !!          Met.Mag. 93, 353-365; UM lineage: STASH m01s20i028,
  !!          pws_snow_prob_diag.F90). The UM's legacy positive-missing-data
  !!          guard (pws_geopht_1000 >= 1.0e8 -> 0) is deliberately dropped:
  !!          Z1000/Z850 are always computed by the calling algorithm, so the
  !!          failure mode it (ineffectively) targeted cannot occur.
  !> @param[in]     nlayers    Number of layers (1: 2d fields)
  !> @param[in]     geopot     Geopotential height (m) at 1000 hPa (slot 1)
  !!                           and 850 hPa (slot 2)
  !> @param[in,out] snow_prob  Snow probability (%), clamped to [0,100]
  !> @param[in]     ndf_geo    Number of degrees of freedom per cell, geopot
  !> @param[in]     undf_geo   Number of total degrees of freedom, geopot
  !> @param[in]     map_geo    Dofmap for the cell at the base of the column
  !> @param[in]     ndf_2d     Number of degrees of freedom per cell, output
  !> @param[in]     undf_2d    Number of total degrees of freedom, output
  !> @param[in]     map_2d     Dofmap for the cell at the base of the column
  subroutine snow_prob_code(nlayers,                    &
                            geopot,                     &
                            snow_prob,                  &
                            ndf_geo, undf_geo, map_geo, &
                            ndf_2d, undf_2d, map_2d)

    implicit none

    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_geo, undf_geo
    integer(kind=i_def), intent(in), dimension(ndf_geo) :: map_geo
    integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
    integer(kind=i_def), intent(in), dimension(ndf_2d)  :: map_2d

    real(kind=r_def), intent(in),    dimension(undf_geo) :: geopot
    real(kind=r_def), intent(inout), dimension(undf_2d)  :: snow_prob

    ! Boyden (1964) regression constants
    real(kind=r_def), parameter :: boyden_intercept  = 5220.0_r_def   ! %
    real(kind=r_def), parameter :: boyden_coeff_1000 = 3.86666_r_def  ! % m-1
    real(kind=r_def), parameter :: boyden_coeff_850  = 4.0_r_def      ! % m-1

    ! Multidata slots fixed by plevs_snow in snow_prob_alg_mod:
    ! map_geo(1) -> Z1000, map_geo(1)+1 -> Z850
    snow_prob(map_2d(1)) = boyden_intercept                        &
                         + boyden_coeff_1000 * geopot(map_geo(1))  &
                         - boyden_coeff_850  * geopot(map_geo(1)+1)

    snow_prob(map_2d(1)) = min( max(snow_prob(map_2d(1)), 0.0_r_def), &
                                100.0_r_def )

  end subroutine snow_prob_code

end module snow_prob_kernel_mod
