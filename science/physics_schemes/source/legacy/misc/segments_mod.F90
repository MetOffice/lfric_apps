! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Defines derived types and subroutines for segmentation control.
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Misc

module segments_mod

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none
private ! Everything is private, unless explicitly exposed later

!------------------------------------------------------------------------------
! Public derived types
!------------------------------------------------------------------------------

! Derived type to contain information about individual segments

!------------------------------------------------------------------------------
! name
!   segment_type
!
! COMPONENTS
!   fp          -- Index of the first point of the segment, where the index runs
!                  from 1 to the number of points in the collapsed array.
!   first_x     -- East-West index of the first point in the segment, where the
!                  index runs from 1 to the number of columns.
!   first_y     -- North-South index of the first point in the segment, where
!                  the index runs from 1 to the number of rows.
!   start_index -- The index in list_points on which this segment starts. Where
!                  the list_points array is not used, start_index and fp should
!                  be numerically equal.
!   end_index   -- The index in list_points on which this segment ends.
!   use_points  -- The number of points in this segment. The last segment size
!                  is calculated from the number of grid points to handle for
!                  this MPI task.
!   seg_points  -- The number of points in this segment. The last segment size
!                  is calculated from the total number of horizontal data points
!                  on the MPI domain on this MPI task.
!   start       -- Index of the first point in the segment of collapsed array,
!                  where the index runs from 1 to the number of points in the
!                  collapsed array.
!   finish      -- Index of the last point in the segment of collapsed array.
!
! NOMENCLATURE
!   collapsed array -- Refers to arrays of rank 2 (or rank 2 slices of
!                      higher-ranked arrays) representing horizontal data
!                      points as if they were collapsed into a single rank 1
!                      array with the standard Fortran data ordering pattern
!                      (column-major).
!
! NOTES
!   The components fp and start are numerically equivalent. They have been kept
!   as two separate names because they were named differently and calculated
!   separately in segmentation control code when it resided in individual model
!   code sections themselves.  The different names help the developer to see
!   that the original segmentation functionality and the functionality in this
!   module are equivalent. It may be wise to remove fp in the future as this
!   module becomes more established.
!
!   There may be a case for reviewing the need for both use_points and
!   seg_points.
!------------------------------------------------------------------------------
type :: segment_type
  integer             :: fp
  integer             :: first_x
  integer             :: first_y
  integer             :: start_index
  integer             :: end_index
  integer             :: use_points
  integer             :: seg_points
  integer             :: start
  integer             :: finish
end type segment_type

! Derived type to contain information about segments as a whole, not individual
! segments.
type :: meta_segment_type
  integer :: num_parallel !The number of members in the parallel team.
  integer :: npnts !The number of grid points to be handled by
                   !each member of the parallel team.
  integer :: fp    !The first grid point to be handled by each
                   !member of the parallel team.
  integer :: lp    !The last grid point to be handled by each
                   !member of the parallel team.
  integer :: num_segments !The total number of segments to be
                          !handled by each member in the parallel team.
  integer :: step   !Step from a point in one
                             !segment to the equivalent point in the next.
end type meta_segment_type

! routines, particularly where arrays with levels are passed through as
! subroutine arguments.

! Make types publicly accessible
public :: segment_type
public :: meta_segment_type

!------------------------------------------------------------------------------
! Public routines
!------------------------------------------------------------------------------

public :: segments_mod_seg_meta        ! Initialises meta_segments variables
public :: segments_mod_segments        ! Initialises segments variables

!------------------------------------------------------------------------------
! Private variables
!------------------------------------------------------------------------------

!DrHook-related parameters
integer(kind=jpim), parameter :: zhook_in     = 0
integer(kind=jpim), parameter :: zhook_out    = 1

character(len=*), parameter, private :: ModuleName='SEGMENTS_MOD'

contains

!-----------------------------------------------------------------------------
! Calculate segment properties
!-----------------------------------------------------------------------------

!------------------------------------------------------------------------------
! SYNOPSIS
!   segments_mod_seg_meta(meta_segments, ipar, num_parallel,  &
!                         grid_points, segment_size, num_segments)
!
! DESCRIPTION
!   Sets up meta-information about the segmentation. The information is placed
!   into a variable of type meta_segment_type.
!
! ARGUMENTS
!   meta_segments -- Output high-level information about segmentation, such as
!                    the number of segments.
!   ipar          -- The ID of the calling parallel team member, counting
!                    from 1.
!   num_parallel  -- The number of members in the parallel team.
!   grid_points   -- The number of grid points to be handled by this MPI task.
!   segment_size  -- The target segment size, if greater than 0.
!   num_segments  -- The number of segments into which to divide the data.
!                    This argument has no effect if segment_size has been set
!                    to a value greater than zero.
!------------------------------------------------------------------------------

subroutine segments_mod_seg_meta(meta_segments, ipar, num_parallel,            &
                                grid_points, segment_size, num_segments)
implicit none

!Arguments
integer,                  intent(in) :: ipar
integer,                  intent(in) :: num_parallel
integer,                  intent(in) :: grid_points
integer,                  intent(in) :: segment_size
integer,                  intent(in) :: num_segments
type(meta_segment_type), intent(out) :: meta_segments

!Internal variables
integer :: num_points
integer :: remainder_points
integer :: i
integer :: running_fp
integer :: running_lp

!DrHook
real(kind=jprb) :: zhook_handle

character(len=*), parameter :: RoutineName='SEGMENTS_MOD_SEG_META'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!Set the number of parallel team members in the derived type.
meta_segments%num_parallel = num_parallel

!Divide the total number of grid points between the number of members of the
!parallel team.
num_points = grid_points/num_parallel
meta_segments%npnts = num_points

!Find the remainder points
remainder_points = mod(grid_points, num_parallel)

!Find the start point and the end point for this parallel team member. The
!Set up the first one first.
running_fp = 0
running_lp = 0

!Loop over the remaining team members, accumulating the first and last points
!as we go. If there is only one team member, this loop will not execute at
!all.
do i = 1,ipar
  running_fp = running_lp + 1
  running_lp = running_lp + num_points
  if (i <= remainder_points) then
    running_lp = running_lp + 1
  end if
end do

!Set up first points and last points in meta_segments type.
meta_segments%fp = running_fp
meta_segments%lp = running_lp

!If this parallel team member is numbered up-to and including the number of
!remainder points, add another point to this member. Just an exercise in
!spreading out the remainder points over different team members.
meta_segments%npnts = num_points
if (ipar <= remainder_points) then
  meta_segments%npnts = num_points + 1
end if

!Different behaviour depending on the segment size
if (segment_size > 0) then
  meta_segments%step = segment_size
  meta_segments%num_segments =                                                 &
    ceiling( real(meta_segments%npnts)                                         &
    / real(meta_segments%step ))
else
  meta_segments%num_segments =                                                 &
    min(meta_segments%npnts, num_segments)

  !Catch the case where the number of points is less than the number of
  !parallel team members. The highest-numbered team members will otherwise
  !attempt to divide by zero, having not acquired any remainder points and
  !hence having no segments. The value of meta_segments%step is otherwise
  !unaffected by the max.
  meta_segments%step = meta_segments%npnts                                     &
    / max(meta_segments%num_segments,1)

end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine segments_mod_seg_meta

!------------------------------------------------------------------------------
! SYNOPSIS
!   segments_mod_segments(segments,meta_segments,row_length,list_points)
!
! DESCRIPTION
!   Calculates the number of points and the start and end points for each
!   segment, and stores them in a derived type array. Uses some pre-calculated
!   meta-information. The information about each segment is stored in the
!   segments array.
!
! ARGUMENTS
!   segments      -- The output derived-type array containing information about
!                    each individual segment.
!   meta_segments -- Derived-type variable containing high-level information
!                    about segmentation, such as the number of segments.
!   row_length    -- The number row length in the whole domain for this
!                    MPI task.
!   list_points (optional) -- Array of points on which to work.
!
! NOTES
!   If a list_points array (optional) is supplied to the routine, then the
!   total number of grid points may not be the same as the sum of the segment
!   sizes.
!------------------------------------------------------------------------------

subroutine segments_mod_segments (segments, meta_segments,                     &
                            row_length, list_points)
implicit none

!Arguments
type(meta_segment_type), intent(in)    :: meta_segments
integer                , intent(in)    :: row_length
type(segment_type) , intent(in out) :: segments(1:meta_segments%num_segments)
integer, optional                      :: list_points(:)

!Local variables
logical :: use_list_points
integer :: first_point
integer :: seg_start
integer :: start_index
integer :: end_index
integer :: seg_points
integer :: use_points
integer :: i

!DrHook
real(kind=jprb) :: zhook_handle

character(len=*), parameter :: RoutineName='SEGMENTS_MOD_SEGMENTS'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!Determine whether a list_points array was provided to the subroutine,
!and set a logical flag accordingly.
use_list_points = .false.
if (present(list_points)) then
  use_list_points = .true.
  first_point = list_points(meta_segments%fp)
else
  first_point = meta_segments%fp
end if

seg_start   = first_point

!Handle segmentation for the calling parallel team member only.
do i = 1, meta_segments%num_segments

  !Note that end_index may be different for the last segment, and is modified
  !later.
  start_index = meta_segments%fp + (i-1) * meta_segments%step
  end_index   = start_index + meta_segments%step -1

  !Find the number of points in this segment, and also the number of points to
  !use in this segment. This if block may reduce further, but kept as it is
  !for the sake of simplicity.
  if (use_list_points) then
    if (i == meta_segments%num_segments) then
      use_points = meta_segments%npnts -                                       &
                meta_segments%step*(meta_segments%num_segments-1)
      end_index  = start_index + use_points - 1
    else
      use_points = meta_segments%step
    end if

    seg_points = list_points(end_index) - seg_start + 1

  else !not using list_points array
    if (i == meta_segments%num_segments) then
      use_points = meta_segments%npnts -                                       &
        meta_segments%step * (meta_segments%num_segments-1)
      end_index  = start_index + use_points - 1
    else
      use_points = meta_segments%step
    end if

    seg_points = use_points

  end if

  !Set first points
  segments(i)%fp         = first_point
  segments(i)%first_y    = (first_point-1)/row_length + 1
  segments(i)%first_x    = first_point-(segments(i)%first_y-1)*row_length

  !Set number of points in segment
  segments(i)%use_points = use_points
  segments(i)%seg_points = seg_points

  !Set start and end point in segments array
  segments(i)%start  = seg_start
  segments(i)%finish = seg_start + seg_points - 1

  !Set start and end indices
  segments(i)%start_index = start_index
  segments(i)%end_index   = end_index

  !Increments
  seg_start   = seg_start   + seg_points
  first_point = first_point + seg_points

end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine segments_mod_segments

end module segments_mod
