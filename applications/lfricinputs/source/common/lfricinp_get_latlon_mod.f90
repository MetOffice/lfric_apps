! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_get_latlon_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int32, int64, real64

use log_mod, only: log_scratch_space, log_event, LOG_LEVEL_ERROR

implicit none

private

public :: get_um_grid_coords, get_lfric_mesh_coords

contains

subroutine get_um_grid_coords(grid_type, idx, idy, lon, lat)
!
! This routine returns the latitude and longitude of a UM grid point with
! indices IDX and IDY. Note this routine is not compatible with variable
! resolution LAMs!
!
use constants_mod,        only: degrees_to_radians
use lfricinp_um_grid_mod, only: um_grid

implicit none

!
! Argument(s)
!
character(len=*),    intent(in)  :: grid_type
integer(kind=int64), intent(in)  :: idx, idy
real(kind=real64),   intent(out) :: lat, lon

if (grid_type == 'p' ) then

  lon = um_grid%p_origin_x + real(idx-1) * um_grid%spacing_x
  lat = um_grid%p_origin_y + real(idy-1) * um_grid%spacing_y

else if (grid_type == 'u') then

  lon = um_grid%u_origin_x + real(idx-1) * um_grid%spacing_x
  lat = um_grid%u_origin_y + real(idy-1) * um_grid%spacing_y

else if (grid_type == 'v') then

  lon = um_grid%v_origin_x + real(idx-1) * um_grid%spacing_x
  lat = um_grid%v_origin_y + real(idy-1) * um_grid%spacing_y

else

  write(log_scratch_space,'(A)') 'UM grid code ' // grid_type // ' was not ' //&
                                 'recognised in routine GET_UM_GRID_COORDS'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)

end if

! Readjust longitude range from [0,360] to [-180,180]
if (lon > 180.0_real64) lon = lon - 360.0_real64

! Convert longitude and latitude to radians
lon = lon * degrees_to_radians
lat = lat * degrees_to_radians

end subroutine get_um_grid_coords


subroutine get_lfric_mesh_coords(cell_lid, lon, lat)
!
! This routine returns the latitude and longitude at the centre of a LFRic mesh
! cell with local cell id CELL_LID
!
use local_mesh_mod,             only: local_mesh_type
use lfricinp_lfric_driver_mod,  only: mesh

implicit none

!
! Argument(s)
!
integer,           intent(in)  :: cell_lid
real(kind=real64), intent(out) :: lon, lat

! Local variables
type(local_mesh_type),  pointer  :: local_mesh => null()
integer(kind=int32)              :: nverts
integer(kind=int32)              :: i
real(kind=real64)                :: vert_coords(2)

local_mesh => mesh%get_local_mesh()

!
! NOTE: Below we take the cell centre coordinates to be the average of the
! vertices' coordinates making up the cell. This is due to the fact the LFRic
! infrastucture does not provide a means of accessing the cell centre
! coordinates from the global mesh object.
!
lon = 0.0_real64
lat = 0.0_real64
nverts = local_mesh%get_nverts_per_cell()
do i = 1, nverts
    call local_mesh%get_vert_coords( local_mesh%get_vert_on_cell(i, cell_lid), &
                                     vert_coords )
    lat = lat + vert_coords(2)
    lon = lon + vert_coords(1)
end do
lat = lat / real(nverts)
lon = lon / real(nverts)

nullify(local_mesh)

end subroutine get_lfric_mesh_coords

end module lfricinp_get_latlon_mod
