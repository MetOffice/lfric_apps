!-----------------------------------------------------------------------------
! (c) Crown copyright 2021 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Module for interpolating profile data
module profile_interp_kernel_mod

use argument_mod,      only: arg_type,                  &
                             GH_FIELD, GH_SCALAR,       &
                             GH_REAL,  GH_LOGICAL,      &
                             GH_WRITE, GH_READ,         &
                             ANY_DISCONTINUOUS_SPACE_1, &
                             CELL_COLUMN

use constants_mod,     only: r_def, i_def, l_def
use kernel_mod,        only: kernel_type

implicit none

private

real(r_def), public    :: profile_data(100)
real(r_def), public    :: profile_heights(100)
integer(i_def), public :: profile_size

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: profile_interp_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                        &
         arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1), &
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                              &
         /)
    integer :: operates_on = CELL_COLUMN
contains
    procedure, nopass :: profile_interp_code
end type profile_interp_kernel_type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: profile_interp_code

contains
!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
!> @brief Interpolate profile data onto a field
!> @param[in] nlayers Number of layers
!> @param[in,out] field Field to be initialised
!> @param[in] height Height coordinate
!> @param[in] ndf Number of degrees of freedom per cell
!> @param[in] undf Number of unique degrees of freedom
!> @param[in] map Dofmap for the cell at the base of the column
subroutine profile_interp_code( nlayers, field, height, &
                                profile_extrap,         &
                                ndf, undf, map )

implicit none

! Arguments
integer(kind=i_def), intent(in)                  :: nlayers
integer(kind=i_def), intent(in)                  :: ndf, undf
integer(kind=i_def), dimension(ndf), intent(in)  :: map

real(kind=r_def), dimension(undf), intent(inout) :: field
real(kind=r_def), dimension(undf), intent(in)    :: height

logical(kind=l_def), intent(in)                  :: profile_extrap

! Loop counter over output profile
integer(kind=i_def) :: k

! map(1) + k
integer(kind=i_def) :: kp

! Index of uppermost dof
integer(kind=i_def) :: top_df

! For each output point with index k, a search is conducted to
! find the index, input_pt, of the input profile coordinate, such
! that output point k lies between input points input_pt and
! input_pt+1
integer(kind=i_def) :: input_pt

! Lowest and highest values of the input profile coordinate, which
! may be strictly monotonically increasing or decreasing
real(kind=r_def)    :: coord_lowest, coord_highest

! Data values of input profile at lowest and highest coordinate positions
real(kind=r_def)    :: data_at_lowest, data_at_highest

! This value is non-negative only when point k of the output profile
! lies between points input_pt and input_pt+1 of the input profile
real(kind=r_def)    :: cell_check

! Linear interpolation weight
real(kind=r_def)    :: interp_weight

! Check whether coordinate is increasing (eg. height), or
! decreasing (eg. pressure)
if ( profile_heights(profile_size) < profile_heights(1) ) then
  ! Decreasing coordinate - lowest value is at upper boundary;
  !                         hightest value is at lower boundary
  coord_lowest    = profile_heights(profile_size)
  data_at_lowest  = profile_data(profile_size)
  coord_highest   = profile_heights(1)
  data_at_highest = profile_data(1)
else
  ! Increasing coordinate - lowest value is at lower boundary
  !                         highest value is at upper boundary
  coord_lowest    = profile_heights(1)
  data_at_lowest  = profile_data(1)
  coord_highest   = profile_heights(profile_size)
  data_at_highest = profile_data(profile_size)
end if

! For Wtheta ndf = 2, loop k = 0, nalyers
! For W3 ndf = 1, loop k = 0, nlayers - 1
top_df = nlayers + ndf - 2

! If the profile size specified every level, just copy the data across
if (top_df+1 == profile_size) then

  do k = 0, top_df
    kp = map(1) + k
    field(kp) = profile_data(k+1)
  end do

else

  ! Loop through output points
  do k = 0, top_df

    kp = map(1) + k

    ! If output point is outside the range of the input coordinate, use
    ! constant extrapolation
    if ( height(kp) <= coord_lowest ) then
      if ( profile_extrap ) field(kp) = data_at_lowest
    else if ( height(kp) >= coord_highest ) then
      if ( profile_extrap ) field(kp) = data_at_highest
    else
      ! Locate the position of the output point on the input coordinate
      do input_pt = 1, profile_size - 1
        cell_check = (profile_heights(input_pt + 1) - height(kp)) * &
                     (height(kp) - profile_heights(input_pt))
        ! Exit search if output point k lies between input_pt and input_pt+1
        if (cell_check >= 0.0_r_def) exit
      end do
      ! Linearly interpolate the input data to the output point
      interp_weight = ( height(kp) - profile_heights(input_pt) ) / &
                      ( profile_heights(input_pt+1) - profile_heights(input_pt) )
      field(kp) = (1.0_r_def - interp_weight) * profile_data(input_pt) + &
                                        interp_weight * profile_data(input_pt+1)
    end if

  end do

end if

end subroutine profile_interp_code

end module profile_interp_kernel_mod
