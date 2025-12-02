! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_reorder_snow_field_mod

!use lfricinp_um_parameters_mod, only: fnamelen, um_rmdi

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : real64, int64


implicit none

private

public :: lfricinp_reorder_snow_field

contains

subroutine lfricinp_reorder_snow_field(field, um_grid)
! Take a snow layer field and converts from UM to
! LFRic ordering
use lfricinp_grid_type_mod, only: lfricinp_grid_type

use log_mod, only: log_event, LOG_LEVEL_ERROR, log_scratch_space, &
                   LOG_LEVEL_INFO
! Arguments
real(kind=real64), allocatable, intent(in out) :: field(:,:)
! Note 1st index is 2D field, 2nd index is pseudo level number
type(lfricinp_grid_type), intent(in):: um_grid

! Local variables
real(kind=real64), allocatable :: field_temp(:,:)
integer(kind=int64) :: total_lev
integer(kind=int64) :: field_size_2d
integer(kind=int64) :: tile_num
integer(kind=int64) :: layer_num, i
integer(kind=int64) :: pseud_lev_count

total_lev = size(field, 2)
field_size_2d = size(field, 1)

if ( total_lev /= um_grid%num_snow_layers * um_grid%num_surface_types ) then
  write(log_scratch_space, '(2(A,I0))') "Mismatch between total_lev = ", &
       total_lev, " and snow_layers * surface types = ", &
       um_grid%num_snow_layers * um_grid%num_surface_types
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

allocate(field_temp(field_size_2d, total_lev))

! UM ordering:
!   All tiles for snow layer 1, followed by all tiles for snow
!   layer 2, then all tiles for snow layer 3. Ie the snow layers
!   for a given tile type are not contiguous
! LFRic ordering:
!   Tile 1, snow layer 1,2,3. Tile 2, snow layer 1,2,3. Ie the
!   snow layers are grouped together
pseud_lev_count = 1
do tile_num = 1, um_grid%num_surface_types
  do layer_num = 1, um_grid%num_snow_layers
    do i = 1, field_size_2d
      field_temp(i, pseud_lev_count) =  field(i, tile_num + &
           ((layer_num - 1) * um_grid%num_surface_types))
    end do
    pseud_lev_count = pseud_lev_count + 1
  end do
end do

call move_alloc(field_temp, field)

end subroutine lfricinp_reorder_snow_field

end module lfricinp_reorder_snow_field_mod
