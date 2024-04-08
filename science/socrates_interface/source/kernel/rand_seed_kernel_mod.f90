!-------------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Set up random seed field for Socrates
!> @details Uses date and time of the model run and a value for the number of
!> unique seeds per degree of latitude / longitude (for MCICA). 100 per degree
!> equates to approximately a 1km square area of the globe (for the Earth).
module rand_seed_kernel_mod

use argument_mod,  only : arg_type,                  &
                          GH_FIELD, GH_SCALAR,       &
                          GH_REAL, GH_INTEGER,       &
                          GH_READ, GH_WRITE,         &
                          CELL_COLUMN,               &
                          ANY_DISCONTINUOUS_SPACE_1
use constants_mod, only : r_def, i_def
use kernel_mod,    only : kernel_type

implicit none

private

type, public, extends(kernel_type) :: rand_seed_kernel_type
  private
  type(arg_type) :: meta_args(7) = (/ &
    arg_type(GH_FIELD,    GH_INTEGER, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1 ), & ! rand_seed
    arg_type(GH_FIELD,    GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1 ), & ! latitude
    arg_type(GH_FIELD,    GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1 ), & ! longitude
    arg_type(GH_SCALAR,   GH_INTEGER, GH_READ                             ), & ! year
    arg_type(GH_SCALAR,   GH_INTEGER, GH_READ                             ), & ! day_of_year
    arg_type(GH_SCALAR,   GH_INTEGER, GH_READ                             ), & ! second_of_day
    arg_type(GH_SCALAR,   GH_INTEGER, GH_READ                             )  & ! seeds_per_degree
    /)

  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: rand_seed_code
end type rand_seed_kernel_type

public :: rand_seed_code

contains

  !> @param[in]    nlayers                    Number of layers
  !> @param[inout] rand_seed                  Field of random seeds
  !> @param[in]    latitude                   Latitude field
  !> @param[in]    longitude                  Longitude field
  !> @param[in]    year                       Year
  !> @param[in]    day_of_year                Day of the Year
  !> @param[in]    second_of_day              Second of the day
  !> @param[in]    seeds_per_degree           Number of unique seeds per degree
  !> @param[in]    ndf                        No. DOFs per cell for W3 2D space
  !> @param[in]    undf                       No. unique DOFs for W3 2D space
  !> @param[in]    map                        Dofmap for W3 space column base cell
  subroutine rand_seed_code( nlayers,          &
                             rand_seed,        &
                             latitude,         &
                             longitude,        &
                             year,             &
                             day_of_year,      &
                             second_of_day,    &
                             seeds_per_degree, &
                             ndf,              &
                             undf,             &
                             map )
    use constants_mod, only : r_def, i_def, radians_to_degrees

    implicit none

    integer(kind = i_def), intent(in) :: nlayers
    integer(kind = i_def), intent(in) :: ndf
    integer(kind = i_def), intent(in), dimension(ndf) :: map
    integer(kind = i_def), intent(in) :: undf
    integer(kind = i_def), intent(in) :: day_of_year
    integer(kind = i_def), intent(in) :: seeds_per_degree, year, second_of_day
    integer(kind = i_def), intent(inout), dimension(undf) :: rand_seed
    real(kind = r_def),    intent(in), dimension(undf) :: latitude
    real(kind = r_def),    intent(in), dimension(undf) :: longitude

    rand_seed(map(1)) &
      = int((latitude(map(1)) * radians_to_degrees + 90.0_r_def) * seeds_per_degree, i_def) &
      + int(longitude(map(1)) * radians_to_degrees * seeds_per_degree, i_def) * 180_i_def * seeds_per_degree &
      + ((abs(year - 2000_i_def) * 366_i_def) + day_of_year) * 86400_i_def + second_of_day

  end subroutine rand_seed_code

end module rand_seed_kernel_mod