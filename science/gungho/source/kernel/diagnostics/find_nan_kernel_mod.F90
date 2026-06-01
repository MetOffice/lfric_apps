!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Finds the locations of NaN values of a field
!>
module find_nan_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD, GH_REAL, GH_READ, &
                               CELL_COLUMN, GH_SCALAR, GH_WRITE,     &
                               ANY_DISCONTINUOUS_SPACE_1,            &
                               ANY_DISCONTINUOUS_SPACE_2,            &
                               ANY_DISCONTINUOUS_SPACE_3,            &
                               ANY_DISCONTINUOUS_SPACE_4
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
  type, public, extends(kernel_type) :: find_nan_kernel_type
    private
    type(arg_type) :: meta_args(9) = (/                                    &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_2), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_4), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_4), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_4), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_4), &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_4)  &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: find_nan_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: find_nan_code

contains

!> @details       For the input field, tests if the current column contains the
!>                minimum or maximum values of the field. If it does, write the
!>                location to the output arrays. If it does not, exit without
!>                doing anything.
!> @param[in]     nlayers    the number of layers
!> @param[in]     dummy      dummy field to allow psyclone to compile
!> @param[in]     field      field to calculate stats of
!> @param[in]     height     height of field above sphere
!> @param[in]     latitude   latitude of field
!> @param[in]     longitude  longitude of field
!> @param[in,out] nan_lev    level of NaN
!> @param[in,out] nan_count  number of times the NaN occurs
!> @param[in,out] nan_lat    latitude of NaN
!> @param[in,out] nan_lon    longitude of NaN
!> @param[in,out] nan_height height of NaN
!> @param[in]     ndf_3d     The number of dofs per cell for 3d field
!> @param[in]     undf_3d    The number of unique dofs for 3d field
!> @param[in]     map_3d     array holding the dofmap for 3d field
!> @param[in]     ndf_2d     The number of dofs per cell for 2d field
!> @param[in]     undf_2d    The number of unique dofs for 2d field
!> @param[in]     map_2d     array holding the dofmap for 2d field
subroutine find_nan_code(nlayers,                    &
                         field, height,              &
                         latitude, longitude,        &
                         nan_lev,                    &
                         nan_count,                  &
                         nan_lat,                    &
                         nan_lon,                    &
                         nan_height,                 &
                         ndf_3d, undf_3d, map_3d,    &
                         ndf_h, undf_h, map_h,       &
                         ndf_2d, undf_2d, map_2d,    &
                         ndf_n, undf_n, map_n )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_3d, undf_3d
  integer(kind=i_def), intent(in) :: ndf_2d, undf_2d
  integer(kind=i_def), intent(in) :: ndf_h, undf_h
  integer(kind=i_def), intent(in) :: ndf_n, undf_n

  real(kind=r_def),    dimension(undf_3d), intent(in) :: field, height
  real(kind=r_def),    dimension(undf_2d), intent(in) :: latitude, longitude
  real(kind=r_def),    dimension(undf_n),  intent(inout) :: nan_lev
  real(kind=r_def),    dimension(undf_n),  intent(inout) :: nan_count
  real(kind=r_def),    dimension(undf_n),  intent(inout) :: nan_lat
  real(kind=r_def),    dimension(undf_n),  intent(inout) :: nan_lon
  real(kind=r_def),    dimension(undf_n),  intent(inout) :: nan_height
  integer(kind=i_def), dimension(ndf_3d),  intent(in) :: map_3d
  integer(kind=i_def), dimension(ndf_2d),  intent(in) :: map_2d
  integer(kind=i_def), dimension(ndf_h),   intent(in) :: map_h
  integer(kind=i_def), dimension(ndf_n),   intent(in) :: map_n

  ! Internal variables
  integer(kind=i_def) :: k, df, dofs_to_loop, top

  ! We loop to nlayers+ndf_3d-2 to ensure that the kernel works correctly for
  ! both w3 and wtheta fields, i.e. we need to loop to nlayers for wtheta
  ! fields (ndf_3d=2), and nlayers-1 for w3 fields (ndf_3d=1)
  select case (ndf_3d)
  case (1, 2)
    dofs_to_loop = 1
    top = nlayers+ndf_3d-2
  case (4, 6)
    dofs_to_loop = 4
    top = nlayers-1
  end select

  do df = 1, dofs_to_loop
    do k = 0, top
      if (field(map_3d(df) + k) /= field(map_3d(df) + k)) then
        ! If the maximum is at this location, write its information
        nan_lev(map_n(1))    = real(k, r_def)
        nan_count(map_n(1))  = 1.0_r_def
        nan_lat(map_n(1))    = latitude(map_2d(df))
        nan_lon(map_n(1))    = longitude(map_2d(df))
        nan_height(map_n(1)) = height(map_h(df) + k)
        exit
      end if
    end do
  end do

  ! If this is a W2 field, loop over vertical DoFs
  if (ndf_3d == 6) then
    df = 5
    do k = 0, nlayers
      if (field(map_3d(df) + k) /= field(map_3d(df) + k)) then
        ! If the maximum is at this location, write its information
        nan_lev(map_n(1))    = real(k, r_def)
        nan_count(map_n(1))  = 1.0_r_def
        nan_lat(map_n(1))    = latitude(map_2d(df))
        nan_lon(map_n(1))    = longitude(map_2d(df))
        nan_height(map_n(1)) = height(map_h(df) + k)
        exit
      end if
    end do
  end if

end subroutine find_nan_code

end module find_nan_kernel_mod
