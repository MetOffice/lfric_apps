!-----------------------------------------------------------------------------
! (c) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates the reciprocal of ocean fraction from sea ice fractions

module set_r_ocean_fraction_mod

use kernel_mod,              only : kernel_type
use argument_mod,            only : arg_type, GH_REAL,         &
                                    GH_INTEGER, GH_FIELD,      &
                                    GH_SCALAR,                 &
                                    GH_WRITE, GH_READ,         &
                                    ANY_DISCONTINUOUS_SPACE_1, &
                                    ANY_DISCONTINUOUS_SPACE_2, &
                                    CELL_COLUMN
use constants_mod,           only : r_def, i_def

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: set_r_ocean_fraction_type
  private
  type(arg_type) :: meta_args(3) = (/                                    &
       arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! r_ocean_fraction
       arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), & ! sea_ice_fraction
       arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                          & ! n_sea_ice_tile,
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: set_r_ocean_fraction_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: set_r_ocean_fraction_code
contains

!> @brief Kernel which calculates the reciprocal of the ocean fraction from
!> @brief the sea ice fractions
!! @param[in] nlayers              Number of layers
!! @param[in,out] r_ocean_fraction Output recipricol of the ocean fraction
!! @param[in] sea_ice_fraction     Input sea ice fractions
!!                                 (as fraction of marine portion of grid box)
!! @param[in] n_sea_ice_tile       Number of sea ice tiles
!! @param[in] ndf_sea              Number of degrees of freedom for ocean
!! @param[in] undf_sea             Number of unique degrees of freedom for ocean
!! @param[in] map_sea              Dofmap for ocean
!! @param[in] ndf_sea_ice          Number of degrees of freedom for sea ice fraction
!! @param[in] undf_sea_ice         Number of unique degrees of freedom for sea ice fraction
!! @param[in] map_sea_ice          Dofmap for sea ice fraction


subroutine set_r_ocean_fraction_code(nlayers,                        &
                      r_ocean_fraction, sea_ice_fraction,            &
                      n_sea_ice_tile,                                &
                      ndf_sea, undf_sea, map_sea,                    &
                      ndf_sea_ice, undf_sea_ice, map_sea_ice)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_sea
  integer(kind=i_def), intent(in) :: undf_sea
  integer(kind=i_def), intent(in) :: map_sea(ndf_sea)
  integer(kind=i_def), intent(in) :: ndf_sea_ice
  integer(kind=i_def), intent(in) :: undf_sea_ice
  integer(kind=i_def), intent(in) :: map_sea_ice(ndf_sea_ice)
  real(kind=r_def), intent(inout) :: r_ocean_fraction(undf_sea)
  real(kind=r_def), intent(in)    :: sea_ice_fraction(undf_sea_ice)
  integer(kind=i_def), intent(in) :: n_sea_ice_tile

  ! Local variables
  integer(kind=i_def) :: i
  real(kind=r_def)    :: total_sea_ice_fraction
  real(kind=r_def)    :: ocean_fraction

  ! Calculate the total sea ice fraction.
  total_sea_ice_fraction = 0.0_r_def
  do i = 1, n_sea_ice_tile
    total_sea_ice_fraction = total_sea_ice_fraction +                & 
                             sea_ice_fraction(map_sea_ice(1) + i-1)
  end do

  ! Calculate the ocean fraction
  ocean_fraction = 1.0_r_def - total_sea_ice_fraction

  ! Calculate the recipricol of the ocean fraction
  if (ocean_fraction > tiny(0.0_r_def)) then
    r_ocean_fraction(map_sea(1)) = 1.0_r_def / ocean_fraction
  else
    r_ocean_fraction(map_sea(1)) = 0.0_r_def
  end if

end subroutine set_r_ocean_fraction_code

end module set_r_ocean_fraction_mod
