! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_nearest_neighbour_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int32, int64, real64

use log_mod, only: log_scratch_space, log_event, LOG_LEVEL_ERROR, LOG_LEVEL_INFO

implicit none

private

public :: find_nn_on_um_grid

contains

subroutine find_nn_on_um_grid(um_mask, um_mask_grid_type, lon_ref, lat_ref,    &
                              idx_lon_nn, idx_lat_nn)
!
! This routine finds the indices of the nearest neighbour on a UM
! mask to a specified point with a given the lon and lat using an ever
! increasing search radius around the specified point.
!
use coord_transform_mod,     only: central_angle
use lfricinp_get_latlon_mod, only: get_um_grid_coords
use lfricinp_um_grid_mod,    only: um_grid
use constants_mod,           only: PI, degrees_to_radians
!
! Argument(s)
!
logical,             intent(in)  :: um_mask(:,:)
character(len=*),    intent(in)  :: um_mask_grid_type
real(kind=real64),   intent(in)  :: lat_ref, lon_ref
integer(kind=int32), intent(out) :: idx_lon_nn, idx_lat_nn

!
! Local variables
!
real(kind=real64)   :: phi_min, phi, lat, lon
real(kind=real64)   :: circle_ang_rad_o, circle_ang_rad, delta_ang_rad
integer(kind=int64) :: ix, iy, ix_w, i
integer(kind=int64) :: ix_start, iy_start, ix_end, iy_end
integer(kind=int64) :: ix_start_o, iy_start_o, ix_end_o, iy_end_o
integer(kind=int64) :: search_point_num, num_on_mask_points
integer(kind=int64), allocatable :: search_points(:,:)

! Allocate array that will contain UM grid (mask) indices that needs to be
! checked during each iteration of the spiral search.
num_on_mask_points = size((um_mask .eqv. .true.))
allocate(search_points(num_on_mask_points,2))

! Make and initial guess to the angular search radius from the reference point
circle_ang_rad_o = sqrt((um_grid%spacing_x)**2 + (um_grid%spacing_y)**2)
circle_ang_rad_o = circle_ang_rad_o * degrees_to_radians
circle_ang_rad   = circle_ang_rad_o

! Set initial value of minimum angular distance of a on mask point from the
! reference point. This variable will be used in a comparison test and updated
! appropriately during the search
phi_min = 2.0_real64 * PI

! Initial limits of the "previous" search circle UM grid bounding box
ix_start_o = 0
iy_start_o = 0
ix_end_o   = 0
iy_end_o   = 0

do

  ! Find indices limits of the UM grid bounding box that contains the current
  ! search circle. Note the indices in the longitude direction can be out of
  ! range (e.g. negative). These must be "wrapped" back onto the proper UM grid
  ! index range. An out of range (e.g. negative) index simply means the search
  ! circle crosses the prime meridian.
  call get_lat_lon_index_limits(um_mask_grid_type, lat_ref, lon_ref,           &
                                circle_ang_rad, iy_start, iy_end,              &
                                ix_start, ix_end)

  ! Find any NEW on mask points that needs checking inside the current search
  ! circle bounding box, whilst ignoring previously searched bounding box points
  !
  ! Re-initialise points to be searched for in this iteration
  search_point_num = 0
  search_points = 0
  !
  ! Loops over all points in the bounding box
  do iy = iy_start, iy_end
    do ix = ix_start, ix_end
      ! Only consider points not in previous searched bounding boxes
      if ( (ix < ix_start_o) .or. (iy < iy_start_o) .or.                       &
           (ix > ix_end_o) .or. (iy > iy_end_o) ) then
        ! Wrap longitude indices to account for case where the search circle
        ! crosses the prime meridian
        ix_w = wrap_lon_indices(um_mask_grid_type, ix)
        ! Add point to current search/check list if its on mask
        if (um_mask(ix_w,iy)) then
          search_point_num = search_point_num + 1
          search_points(search_point_num,1) = ix_w
          search_points(search_point_num,2) = iy
        end if
      end if
    end do
  end do

  ! Loop over all points that needs their distances checked during this
  ! iteration
  do i = 1, search_point_num
    ix = search_points(i,1)
    iy = search_points(i,2)
    call get_um_grid_coords(um_mask_grid_type, ix, iy, lon, lat)
    call central_angle(lon_ref, lat_ref, lon, lat, phi)
    ! If the current point distance is less than the current on mask nearest
    ! neighbour found, update the nearest neighbour with current point
    if (phi < phi_min) then
      idx_lon_nn = int(ix, kind=int32)
      idx_lat_nn = int(iy, kind=int32)
      phi_min = phi
    end if
  end do

  write(log_scratch_space,*) '*! ', lon_ref, lat_ref, search_point_num,        &
                                    num_on_mask_points, phi_min, circle_ang_rad
  call log_event(log_scratch_space, LOG_LEVEL_INFO)

  ! If current nearest neighbour distance estimate is GREATER than the current
  ! search circle radius, then adjust the search circle radius appropriately.
  ! Note this usually occur during the first number of iterations where no on
  ! mask points are found within the search radius. This can also happen when
  ! an on mask point is found on the bounding box, but it lies outside the
  ! search circle itself, e.g. if the on mask point lies on a corner of the
  ! bounding box.
  !
  ! If current nearest neighbour distance estimate is LESS than the current
  ! search circle radius, that means the NN on mask point has successfully been
  ! found!
  if ( phi_min > circle_ang_rad ) then
    ! Adjust the search radius. It is done so that the search radius initially
    ! increases rapidly, from a presumably small value - circle_ang_rad_o - and
    ! plateaus to not more than 20% of its current value. However if the current
    ! best nearest NN candidate is closer, then the search radius is increased
    ! to only slightly above the current best NN distance, as there is no need
    ! to look at larger distances.
    delta_ang_rad = max(0.20_real64*circle_ang_rad, circle_ang_rad_o)
    if ( (phi_min-circle_ang_rad) > delta_ang_rad) then
      circle_ang_rad = circle_ang_rad + delta_ang_rad
    else
      circle_ang_rad = 1.05_real64 * phi_min
    end if
    ! Store current UM grid bounding box boundaries for use in next iteration
    ix_start_o = ix_start
    iy_start_o = iy_start
    ix_end_o   = ix_end
    iy_end_o   = iy_end
  else
    ! If this section is reached it means the on mask nearest neighbour has been
    ! found, so terminate the search by exiting the search/iteration loop
    exit
  end if

end do

deallocate(search_points)

end subroutine find_nn_on_um_grid


function wrap_lon_indices(grid_type, idx) result(idx_w)
!
! Wrap out of bounds (e.g. negative) longitude indices back on to
! the proper UM grid index range
!
use lfricinp_um_grid_mod, only: um_grid

implicit none

character(len=*),    intent(in)  :: grid_type
integer(kind=int64), intent(in)  :: idx
integer(kind=int64)              :: idx_w, nx

if (grid_type == 'p' ) then
  nx  = um_grid%num_p_points_x
else if (grid_type == 'u') then
  nx  = um_grid%num_u_points_x
else if (grid_type == 'v') then
  nx  = um_grid%num_v_points_x
else
  write(log_scratch_space,'(A)') 'UM grid code ' // grid_type // ' was not '// &
                                'recognised in routine WRAP_LON_INDICES'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

if (idx > nx) then
  idx_w = mod(idx, nx)
else if (idx < 1) then
  idx_w = mod(idx, nx) + nx
else
  idx_w = idx
end if

end function wrap_lon_indices


subroutine get_lat_lon_index_limits(grid_type, lat, lon, circle_ang_rad,       &
                                    iy_start, iy_end, ix_start, ix_end)

!
! Given a circle with angular radius CIRCLE_AND_RAD centered on
! the point with coordinates LAT & LON find the limits of the smallest
! UM grid bounding box that contains the circle.
!
use constants_mod,        only: radians_to_degrees
use lfricinp_um_grid_mod, only: um_grid

implicit none

!
! Argument(s)
!
character(len=*),    intent(in)  :: grid_type
real(kind=real64),   intent(in)  :: lat, lon
real(kind=real64),   intent(in)  :: circle_ang_rad
integer(kind=int64), intent(out) :: iy_start, iy_end, ix_start, ix_end

!
! Local variables
real(kind=real64)   :: lat_deg, lon_deg
real(kind=real64)   :: origin_x, origin_y, phi_x, phi_y
integer(kind=int64) :: nx, ny


!
! Convert longitude and latitude to degrees
lon_deg = lon * radians_to_degrees
lat_deg = lat * radians_to_degrees

! Readjust longitude range from [-180,180] to [0,360]
if (lon_deg < 0.0_real64) lon_deg = lon_deg + 360.0_real64

! Set the number points int he latitude and longitude direction on the grid and
! the lat lon origin values
if (grid_type == 'p' ) then

  nx  = um_grid%num_p_points_x
  ny  = um_grid%num_p_points_y
  origin_x = um_grid%p_origin_x
  origin_y = um_grid%p_origin_y

else if (grid_type == 'u') then

  nx  = um_grid%num_u_points_x
  ny  = um_grid%num_u_points_y
  origin_x = um_grid%u_origin_x
  origin_y = um_grid%u_origin_y

else if (grid_type == 'v') then

  nx  = um_grid%num_v_points_x
  ny  = um_grid%num_v_points_y
  origin_x = um_grid%v_origin_x
  origin_y = um_grid%v_origin_y

else

  write(log_scratch_space,'(A)') 'UM grid code ' // grid_type // ' was not '// &
                                'recognised in routine GET_LAT_LON_INDEX_LIMITS'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)

end if

! Note the maximum change in longitude away from the centre of the search
! circle, here denoted by phi_x has not been strictly verified. This needs to
! be independently confirmed in future. For now it appears to work. The
! corresponding variable for latitude, phi_y, is self-evident however as
! latitudes are strictly measured along a meridian which is a great circle.
phi_x = (asin(circle_ang_rad) / cos(lat)) * radians_to_degrees
phi_y = circle_ang_rad * radians_to_degrees

! Find longitude index limits of the UM grid box that contains the search
! circle. Note if the search cicle includes a pole the WHOLE longitude range
! is used. Also these "relative" longitude indices can be negative. An out of
! range (e.g negative) index simply indicates the circle has crossed the prime
! meridian. This is not "corrected", to preserve the implied ORDERING when
! looping over the indices. The calling routine will have to "wrap" the out of
! range indices back on to the proper UM grid index range, e.g. when they are
! used inside the body of said loop for instance
if ( ((lat_deg+phi_y) > 90.0_real64) .or. ((lat_deg-phi_y) < -90.0_real64) ) then
  ix_start = 1
  ix_end   = nx
else
  ix_start = ceiling(((lon_deg - phi_x) - origin_x) / um_grid%spacing_x) + 1
  ix_end   = ceiling(((lon_deg + phi_x) - origin_x) / um_grid%spacing_x)
ENDIF

! Find latitude index limits of the UM grid bounding box that contains the
! search circle. Note the latitude indices are limited between maximum and
! minimum values on the grid, unlike the case for the longitude indices above
! which are allowed the be out of range.
iy_start = ceiling(((lat_deg - phi_y) - origin_y) / um_grid%spacing_y) + 1
if(iy_start < 1) iy_start = 1
iy_end   = ceiling(((lat_deg + phi_y) - origin_y) / um_grid%spacing_y)
if(iy_end > ny)  iy_end  = ny

end subroutine get_lat_lon_index_limits

end module lfricinp_nearest_neighbour_mod
