!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculates the location of the minimim and maximum values of a field
!>
module field_stats_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD, GH_REAL, GH_READ, &
                               CELL_COLUMN, GH_SCALAR, GH_WRITE,     &
                               ANY_DISCONTINUOUS_SPACE_1,            &
                               ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,     only: r_def, i_def
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: field_stats_kernel_type
    private
    type(arg_type) :: meta_args(16) = (/                                   &
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1),  &
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1),  &
         arg_type(GH_SCALAR,GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR,GH_REAL, GH_READ),                             &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: field_stats_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: field_stats_code

contains

!> @details       For the input field, tests if the current column contains the
!>                minimum or maximum values of the field. If it does, write the
!>                location to the output arrays. If it does not, exit without
!>                doing anything.
!> @param[in]     nlayers    the number of layers
!> @param[in]     field      field to calculate stats of
!> @param[in]     height     height of field above sphere
!> @param[in]     fmax       maximum value of the field
!> @param[in]     fmin       minimum value of the field
!> @param[in]     latitude   latitude of field
!> @param[in]     longitude  longitude of field
!> @param[in,out] max_lev    level of maximum value
!> @param[in,out] min_lev    level of minimum value
!> @param[in,out] max_count  number of times the maximum occurs
!> @param[in,out] min_count  number of times the minimum occurs
!> @param[in,out] max_lat    latitude of maximum value
!> @param[in,out] min_lat    latitude of minimum value
!> @param[in,out] max_lon    longitude of maximum value
!> @param[in,out] min_lon    longitude of minimum value
!> @param[in,out] max_height height of maximum value
!> @param[in,out] min_height height of minimum value
!> @param[in]     ndf_3d     The number of dofs per cell for 3d field
!> @param[in]     undf_3d    The number of unique dofs for 3d field
!> @param[in]     map_3d     array holding the dofmap for 3d field
!> @param[in]     ndf_2d     The number of dofs per cell for 2d field
!> @param[in]     undf_2d    The number of unique dofs for 2d field
!> @param[in]     map_2d     array holding the dofmap for 2d field
subroutine field_stats_code(nlayers,                    &
                            field, height,              &
                            fmax, fmin,                 &
                            latitude, longitude,        &
                            max_lev, min_lev,           &
                            max_count, min_count,       &
                            max_lat, min_lat,           &
                            max_lon, min_lon,           &
                            max_height, min_height,     &
                            ndf_3d, undf_3d, map_3d,    &
                            ndf_2d, undf_2d, map_2d     &
                            )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_3d, undf_3d
  integer(kind=i_def), intent(in) :: ndf_2d, undf_2d

  real(r_def), intent(in) :: fmax, fmin
  real(kind=r_def), dimension(undf_3d), intent(in) :: field, height
  real(kind=r_def), dimension(undf_2d), intent(in) :: latitude, longitude
  real(kind=r_def), dimension(undf_2d), intent(inout) :: max_lev, min_lev
  real(kind=r_def), dimension(undf_2d), intent(inout) :: max_count, min_count
  real(kind=r_def), dimension(undf_2d), intent(inout) :: max_lat, min_lat
  real(kind=r_def), dimension(undf_2d), intent(inout) :: max_lon, min_lon
  real(kind=r_def), dimension(undf_2d), intent(inout) :: max_height, min_height
  integer(kind=i_def), dimension(ndf_3d),  intent(in) :: map_3d
  integer(kind=i_def), dimension(ndf_2d),  intent(in) :: map_2d

  ! Internal variables
  integer(kind=i_def) :: k

  ! We loop to nlayers+ndf_3d-2 to ensure that the kernel works correctly for
  ! both w3 and wtheta fields, i.e. we need to loop to nlayers for wtheta
  ! fields (ndf_3d=2), and nlayers-1 for w3 fields (ndf_3d=1)
  do k = 0, nlayers+ndf_3d-2

    if (field(map_3d(1) + k) >= fmax) then
      ! If the maximum is at this location, write its information
      max_lev(map_2d(1))    = real(k, r_def)
      max_count(map_2d(1))  = 1.0_r_def
      max_lat(map_2d(1))    = latitude(map_2d(1))
      max_lon(map_2d(1))    = longitude(map_2d(1))
      max_height(map_2d(1)) = height(map_3d(1) + k)
      exit
    end if

  end do

  ! We loop to nlayers+ndf_3d-2 to ensure that the kernel works correctly for
  ! both w3 and wtheta fields, i.e. we need to loop to nlayers for wtheta
  ! fields (ndf_3d=2), and nlayers-1 for w3 fields (ndf_3d=1)
  do k = 0, nlayers+ndf_3d-2

    if (field(map_3d(1) + k) <= fmin) then
      ! If the minimum is at this location, write its information
      min_lev(map_2d(1))    = real(k, r_def)
      min_count(map_2d(1))  = 1.0_r_def
      min_lat(map_2d(1))    = latitude(map_2d(1))
      min_lon(map_2d(1))    = longitude(map_2d(1))
      min_height(map_2d(1)) = height(map_3d(1) + k)
      exit
    end if

  end do

end subroutine field_stats_code

end module field_stats_kernel_mod
