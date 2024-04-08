! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_nearest_neighbour_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int32, int64, real64

USE log_mod, ONLY: log_scratch_space, log_event, LOG_LEVEL_ERROR, LOG_LEVEL_INFO

IMPLICIT NONE

PRIVATE

PUBLIC :: find_nn_on_um_grid

CONTAINS

SUBROUTINE find_nn_on_um_grid(um_mask, um_mask_grid_type, lon_ref, lat_ref,    &
                              idx_lon_nn, idx_lat_nn)
!
! This routine finds the indices of the nearest neighbour on a UM
! mask to a specified point with a given the lon and lat using an ever
! increasing search radius around the specified point.
!
USE coord_transform_mod,     ONLY: central_angle
USE lfricinp_get_latlon_mod, ONLY: get_um_grid_coords
USE lfricinp_um_grid_mod,    ONLY: um_grid
USE constants_mod,           ONLY: PI, degrees_to_radians
!
! Argument(s)
!
LOGICAL,             INTENT(IN)  :: um_mask(:,:)
CHARACTER(LEN=*),    INTENT(IN)  :: um_mask_grid_type
REAL(KIND=real64),   INTENT(IN)  :: lat_ref, lon_ref
INTEGER(KIND=int32), INTENT(OUT) :: idx_lon_nn, idx_lat_nn

!
! Local variables
!
REAL(KIND=real64)   :: phi_min, phi, lat, lon
REAL(KIND=real64)   :: circle_ang_rad_o, circle_ang_rad, delta_ang_rad
INTEGER(KIND=int64) :: ix, iy, ix_w, i
INTEGER(KIND=int64) :: ix_start, iy_start, ix_end, iy_end
INTEGER(KIND=int64) :: ix_start_o, iy_start_o, ix_end_o, iy_end_o
INTEGER(KIND=int64) :: search_point_num, num_on_mask_points
INTEGER(KIND=int64), ALLOCATABLE :: search_points(:,:)

! Allocate array that will contain UM grid (mask) indices that needs to be
! checked during each iteration of the spiral search.
num_on_mask_points = SIZE((um_mask .eqv. .TRUE.))
ALLOCATE(search_points(num_on_mask_points,2))

! Make and initial guess to the angular search radius from the reference point
circle_ang_rad_o = SQRT((um_grid%spacing_x)**2 + (um_grid%spacing_y)**2)
circle_ang_rad_o = circle_ang_rad_o * degrees_to_radians
circle_ang_rad   = circle_ang_rad_o

! Set initial value of minimum angular distance of a on mask point from the
! reference point. This variable will be used in a comparison test and updated
! appropriately during the search
phi_min = 2.0 * PI

! Initial limits of the "previous" search circle UM grid bounding box
ix_start_o = 0
iy_start_o = 0
ix_end_o   = 0
iy_end_o   = 0

DO

  ! Find indices limits of the UM grid bounding box that contains the current
  ! search circle. Note the indices in the longitude direction can be out of
  ! range (e.g. negative). These must be "wrapped" back onto the proper UM grid
  ! index range. An out of range (e.g. negative) index simply means the search
  ! circle crosses the prime meridian.
  CALL get_lat_lon_index_limits(um_mask_grid_type, lat_ref, lon_ref,           &
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
  DO iy = iy_start, iy_end
    DO ix = ix_start, ix_end
      ! Only consider points not in previous searched bounding boxes
      IF ( (ix < ix_start_o) .OR. (iy < iy_start_o) .OR.                       &
           (ix > ix_end_o) .OR. (iy > iy_end_o) ) THEN
        ! Wrap longitude indices to account for case where the search circle
        ! crosses the prime meridian
        ix_w = wrap_lon_indices(um_mask_grid_type, ix)
        ! Add point to current search/check list if its on mask
        IF (um_mask(ix_w,iy)) THEN
          search_point_num = search_point_num + 1
          search_points(search_point_num,1) = ix_w
          search_points(search_point_num,2) = iy
        END IF
      END IF
    END DO
  END DO

  ! Loop over all points that needs their distances checked during this
  ! iteration
  DO i = 1, search_point_num
    ix = search_points(i,1)
    iy = search_points(i,2)
    CALL get_um_grid_coords(um_mask_grid_type, ix, iy, lon, lat)
    CALL central_angle(lon_ref, lat_ref, lon, lat, phi)
    ! If the current point distance is less than the current on mask nearest
    ! neighbour found, update the nearest neighbour with current point
    IF (phi < phi_min) THEN
      idx_lon_nn = INT(ix, KIND=int32)
      idx_lat_nn = INT(iy, KIND=int32)
      phi_min = phi
    END IF
  END DO

  WRITE(log_scratch_space,*) '*! ', lon_ref, lat_ref, search_point_num,        &
                                    num_on_mask_points, phi_min, circle_ang_rad
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

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
  IF ( phi_min > circle_ang_rad ) THEN
    ! Adjust the search radius. It is done so that the search radius initially
    ! increases rapidly, from a presumably small value - circle_ang_rad_o - and
    ! plateaus to not more than 20% of its current value. However if the current
    ! best nearest NN candidate is closer, then the search radius is increased
    ! to only slightly above the current best NN distance, as there is no need
    ! to look at larger distances.
    delta_ang_rad = MAX(0.20*circle_ang_rad, circle_ang_rad_o)
    IF ( (phi_min-circle_ang_rad) > delta_ang_rad) THEN
      circle_ang_rad = circle_ang_rad + delta_ang_rad
    ELSE
      circle_ang_rad = 1.05 * phi_min
    END IF
    ! Store current UM grid bounding box boundaries for use in next iteration
    ix_start_o = ix_start
    iy_start_o = iy_start
    ix_end_o   = ix_end
    iy_end_o   = iy_end
  ELSE
    ! If this section is reached it means the on mask nearest neighbour has been
    ! found, so terminate the search by exiting the search/iteration loop
    EXIT
  END IF

END DO

DEALLOCATE(search_points)

END SUBROUTINE find_nn_on_um_grid


FUNCTION wrap_lon_indices(grid_type, idx) RESULT(idx_w)
!
! Wrap out of bounds (e.g. negative) longitude indices back on to
! the proper UM grid index range
!
USE lfricinp_um_grid_mod, ONLY: um_grid

IMPLICIT NONE

CHARACTER(LEN=*),    INTENT(IN)  :: grid_type
INTEGER(KIND=int64), INTENT(IN)  :: idx
INTEGER(KIND=int64)              :: idx_w, nx

IF (grid_type == 'p' ) THEN
  nx  = um_grid%num_p_points_x
ELSE IF (grid_type == 'u') THEN
  nx  = um_grid%num_u_points_x
ELSE IF (grid_type == 'v') THEN
  nx  = um_grid%num_v_points_x
ELSE
  WRITE(log_scratch_space,'(A)') 'UM grid code ' // grid_type // ' was not '// &
                                'recognised in routine WRAP_LON_INDICES'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

IF (idx > nx) THEN
  idx_w = MOD(idx, nx)
ELSE IF (idx < 1) THEN
  idx_w = MOD(idx, nx) + nx
ELSE
  idx_w = idx
END IF

END FUNCTION wrap_lon_indices


SUBROUTINE get_lat_lon_index_limits(grid_type, lat, lon, circle_ang_rad,       &
                                    iy_start, iy_end, ix_start, ix_end)

!
! Given a circle with angular radius CIRCLE_AND_RAD centered on
! the point with coordinates LAT & LON find the limits of the smallest
! UM grid bounding box that contains the circle.
!
USE constants_mod,        ONLY: radians_to_degrees
USE lfricinp_um_grid_mod, ONLY: um_grid

IMPLICIT NONE

!
! Argument(s)
!
CHARACTER(LEN=*),    INTENT(IN)  :: grid_type
REAL(KIND=real64),   INTENT(IN)  :: lat, lon
REAL(KIND=real64),   INTENT(IN)  :: circle_ang_rad
INTEGER(KIND=int64), INTENT(OUT) :: iy_start, iy_end, ix_start, ix_end

!
! Local variables
REAL(KIND=real64)   :: lat_deg, lon_deg
REAL(KIND=real64)   :: origin_x, origin_y, phi_x, phi_y
INTEGER(KIND=int64) :: nx, ny


!
! Convert longitude and latitude to degrees
lon_deg = lon * radians_to_degrees
lat_deg = lat * radians_to_degrees

! Readjust longitude range from [-180,180] to [0,360]
IF (lon_deg < 0.0) lon_deg = lon_deg + 360.0

! Set the number points int he latitude and longitude direction on the grid and
! the lat lon origin values
IF (grid_type == 'p' ) THEN

  nx  = um_grid%num_p_points_x
  ny  = um_grid%num_p_points_y
  origin_x = um_grid%p_origin_x
  origin_y = um_grid%p_origin_y

ELSE IF (grid_type == 'u') THEN

  nx  = um_grid%num_u_points_x
  ny  = um_grid%num_u_points_y
  origin_x = um_grid%u_origin_x
  origin_y = um_grid%u_origin_y

ELSE IF (grid_type == 'v') THEN

  nx  = um_grid%num_v_points_x
  ny  = um_grid%num_v_points_y
  origin_x = um_grid%v_origin_x
  origin_y = um_grid%v_origin_y

ELSE

  WRITE(log_scratch_space,'(A)') 'UM grid code ' // grid_type // ' was not '// &
                                'recognised in routine GET_LAT_LON_INDEX_LIMITS'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)

END IF

! Note the maximum change in longitude away from the centre of the search
! circle, here denoted by phi_x has not been strictly verified. This needs to
! be independently confirmed in future. For now it appears to work. The
! corresponding variable for latitude, phi_y, is self-evident however as
! latitudes are strictly measured along a meridian which is a great circle.
phi_x = (ASIN(circle_ang_rad) / COS(lat)) * radians_to_degrees
phi_y = circle_ang_rad * radians_to_degrees

! Find longitude index limits of the UM grid box that contains the search
! circle. Note if the search cicle includes a pole the WHOLE longitude range
! is used. Also these "relative" longitude indices can be negative. An out of
! range (e.g negative) index simply indicates the circle has crossed the prime
! meridian. This is not "corrected", to preserve the implied ORDERING when
! looping over the indices. The calling routine will have to "wrap" the out of
! range indices back on to the proper UM grid index range, e.g. when they are
! used inside the body of said loop for instance
IF ( ((lat_deg+phi_y) > 90.0) .OR. ((lat_deg-phi_y) < -90.0) ) THEN
  ix_start = 1
  ix_end   = nx
ELSE
  ix_start = CEILING(((lon_deg - phi_x) - origin_x) / um_grid%spacing_x) + 1
  ix_end   = CEILING(((lon_deg + phi_x) - origin_x) / um_grid%spacing_x)
ENDIF

! Find latitude index limits of the UM grid bounding box that contains the
! search circle. Note the latitude indices are limited between maximum and
! minimum values on the grid, unlike the case for the longitude indices above
! which are allowed the be out of range.
iy_start = CEILING(((lat_deg - phi_y) - origin_y) / um_grid%spacing_y) + 1
IF(iy_start < 1) iy_start = 1
iy_end   = CEILING(((lat_deg + phi_y) - origin_y) / um_grid%spacing_y)
IF(iy_end > ny)  iy_end  = ny

END SUBROUTINE get_lat_lon_index_limits

END MODULE lfricinp_nearest_neighbour_mod
