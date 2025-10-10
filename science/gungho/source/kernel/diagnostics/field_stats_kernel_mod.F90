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
  use constants_mod,     only: r_def, i_def, r_single, r_double
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
    type(arg_type) :: meta_args(20) = (/                                   &
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1),  &
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1),  &
         arg_type(GH_SCALAR,GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR,GH_REAL, GH_READ),                             &
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
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_2)  &
         /)
    integer :: operates_on = CELL_COLUMN
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: field_stats_code

  interface field_stats_code
    module procedure &
      field_stats_code_r_single, &
      field_stats_code_r_double
    end interface

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
!> @param[in]     fmax_excess threshold above which to count values of the field
!> @param[in]     fmin_excess threshold above which to count values of the field
!> @param[in]     latitude   latitude of field
!> @param[in]     longitude  longitude of field
!> @param[in,out] max_lev    level of maximum value
!> @param[in,out] min_lev    level of minimum value
!> @param[in,out] max_excess_count  num of times the upper threshold is exceeded
!> @param[in,out] min_excess_count  num of times the lower threshold is exceeded
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
subroutine field_stats_code_r_single(nlayers,                            &
                                     field, height,                      &
                                     fmax, fmin,                         &
                                     fmax_excess, fmin_excess,           &
                                     latitude, longitude,                &
                                     max_lev, min_lev,                   &
                                     max_excess_count, min_excess_count, &
                                     max_count, min_count,               &
                                     max_lat, min_lat,                   &
                                     max_lon, min_lon,                   &
                                     max_height, min_height,             &
                                     ndf_3d, undf_3d, map_3d,            &
                                     ndf_2d, undf_2d, map_2d             &
                                    )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_3d, undf_3d
  integer(kind=i_def), intent(in) :: ndf_2d, undf_2d

  real(kind=r_single), intent(in) :: fmax, fmin, fmax_excess, fmin_excess
  real(kind=r_single), dimension(undf_3d), intent(in)    :: field
  real(kind=r_def),    dimension(undf_3d), intent(in)    :: height
  real(kind=r_def),    dimension(undf_2d), intent(in)    :: latitude, longitude
  real(kind=r_def),    dimension(undf_2d), intent(inout) :: max_lev, min_lev
  real(kind=r_def),    dimension(undf_2d), intent(inout) :: max_count, min_count
  real(kind=r_def),    dimension(undf_2d), intent(inout) :: max_excess_count
  real(kind=r_def),    dimension(undf_2d), intent(inout) :: min_excess_count
  real(kind=r_def),    dimension(undf_2d), intent(inout) :: max_lat, min_lat
  real(kind=r_def),    dimension(undf_2d), intent(inout) :: max_lon, min_lon
  real(kind=r_def),    dimension(undf_2d), intent(inout) :: max_height, min_height
  integer(kind=i_def), dimension(ndf_3d),  intent(in)    :: map_3d
  integer(kind=i_def), dimension(ndf_2d),  intent(in)    :: map_2d

  ! Internal variables
  integer(kind=i_def) :: df, k, k_max, df_max

  if (ndf_3d == 4) then
    df_max = 4
  else
    df_max = 1
  end if

  if (ndf_3d == 2) then
    k_max = nlayers
  else
    k_max = nlayers - 1
  end if

  ! Count number of points which exceed thresholds -----------------------------
  do df = 1, df_max
    max_excess_count(map_2d(df)) = 0.0_r_def
    min_excess_count(map_2d(df)) = 0.0_r_def
    do k = 0, k_max
      if (field(map_3d(df) + k) >= fmax_excess) then
        ! If the field exceeds the specified threshold, increment the counter
        max_excess_count(map_2d(df)) = max_excess_count(map_2d(df)) + 1.0_r_def
      end if

      if (field(map_3d(df) + k) <= fmin_excess) then
        ! If the field exceeds the specified threshold, increment the counter
        min_excess_count(map_2d(df)) = min_excess_count(map_2d(df)) + 1.0_r_def
      end if
    end do
  end do

  ! Search for location of global maximum --------------------------------------
  do df = 1, df_max
    do k = 0, k_max

      if (field(map_3d(df) + k) >= fmax) then
        ! If the maximum is at this location, write its information
        max_lev(map_2d(df))    = real(k, r_def)
        max_lat(map_2d(df))    = latitude(map_2d(df))
        max_lon(map_2d(df))    = longitude(map_2d(df))
        max_height(map_2d(df)) = height(map_3d(df) + k)
        max_count(map_2d(df)) = 1.0_r_def
        exit
      end if
    end do
  end do

  ! Search for location of global minimum --------------------------------------
  do df = 1, df_max
    do k = 0, k_max
      if (field(map_3d(df) + k) <= fmin) then
        ! If the minimum is at this location, write its information
        min_lev(map_2d(df))    = real(k, r_def)
        min_lat(map_2d(df))    = latitude(map_2d(df))
        min_lon(map_2d(df))    = longitude(map_2d(df))
        min_height(map_2d(df)) = height(map_3d(df) + k)
        min_count(map_2d(df)) = 1.0_r_def
        exit
      end if
    end do
  end do

end subroutine field_stats_code_r_single

subroutine field_stats_code_r_double(nlayers,                            &
                                     field, height,                      &
                                     fmax, fmin,                         &
                                     fmax_excess, fmin_excess,           &
                                     latitude, longitude,                &
                                     max_lev, min_lev,                   &
                                     max_excess_count, min_excess_count, &
                                     max_count, min_count,               &
                                     max_lat, min_lat,                   &
                                     max_lon, min_lon,                   &
                                     max_height, min_height,             &
                                     ndf_3d, undf_3d, map_3d,            &
                                     ndf_2d, undf_2d, map_2d             &
                                    )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_3d, undf_3d
  integer(kind=i_def), intent(in) :: ndf_2d, undf_2d

  real(kind=r_double), intent(in) :: fmax, fmin, fmax_excess, fmin_excess
  real(kind=r_double), dimension(undf_3d), intent(in)    :: field
  real(kind=r_def),    dimension(undf_3d), intent(in)    :: height
  real(kind=r_def),    dimension(undf_2d), intent(in)    :: latitude, longitude
  real(kind=r_def),    dimension(undf_2d), intent(inout) :: max_lev, min_lev
  real(kind=r_def),    dimension(undf_2d), intent(inout) :: max_count, min_count
  real(kind=r_def),    dimension(undf_2d), intent(inout) :: max_excess_count
  real(kind=r_def),    dimension(undf_2d), intent(inout) :: min_excess_count
  real(kind=r_def),    dimension(undf_2d), intent(inout) :: max_lat, min_lat
  real(kind=r_def),    dimension(undf_2d), intent(inout) :: max_lon, min_lon
  real(kind=r_def),    dimension(undf_2d), intent(inout) :: max_height, min_height
  integer(kind=i_def), dimension(ndf_3d),  intent(in)    :: map_3d
  integer(kind=i_def), dimension(ndf_2d),  intent(in)    :: map_2d

  ! Internal variables
  integer(kind=i_def) :: df, k, k_max, df_max

  if (ndf_3d == 4) then
    df_max = 4
  else
    df_max = 1
  end if

  if (ndf_3d == 2) then
    k_max = nlayers
  else
    k_max = nlayers - 1
  end if

  ! Count number of points which exceed thresholds -----------------------------
  do df = 1, df_max
    max_excess_count(map_2d(df)) = 0.0_r_def
    min_excess_count(map_2d(df)) = 0.0_r_def
    do k = 0, k_max
      if (field(map_3d(df) + k) >= fmax_excess) then
        ! If the field exceeds the specified threshold, increment the counter
        max_excess_count(map_2d(df)) = max_excess_count(map_2d(df)) + 1.0_r_def
      end if

      if (field(map_3d(df) + k) <= fmin_excess) then
        ! If the field exceeds the specified threshold, increment the counter
        min_excess_count(map_2d(df)) = min_excess_count(map_2d(df)) + 1.0_r_def
      end if
    end do
  end do

  ! Search for location of global maximum --------------------------------------
  do df = 1, df_max
    do k = 0, k_max

      if (field(map_3d(df) + k) >= fmax) then
        ! If the maximum is at this location, write its information
        max_lev(map_2d(df))    = real(k, r_def)
        max_lat(map_2d(df))    = latitude(map_2d(df))
        max_lon(map_2d(df))    = longitude(map_2d(df))
        max_height(map_2d(df)) = height(map_3d(df) + k)
        max_count(map_2d(df)) = 1.0_r_def
        exit
      end if
    end do
  end do

  ! Search for location of global minimum --------------------------------------
  do df = 1, df_max
    do k = 0, k_max
      if (field(map_3d(df) + k) <= fmin) then
        ! If the minimum is at this location, write its information
        min_lev(map_2d(df))    = real(k, r_def)
        min_lat(map_2d(df))    = latitude(map_2d(df))
        min_lon(map_2d(df))    = longitude(map_2d(df))
        min_height(map_2d(df)) = height(map_3d(df) + k)
        min_count(map_2d(df)) = 1.0_r_def
        exit
      end if
    end do
  end do

end subroutine field_stats_code_r_double

end module field_stats_kernel_mod
