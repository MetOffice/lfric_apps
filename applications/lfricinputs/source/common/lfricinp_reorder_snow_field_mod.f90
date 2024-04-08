! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_reorder_snow_field_mod

!USE lfricinp_um_parameters_mod, ONLY: fnamelen, um_rmdi

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : real64, int64


IMPLICIT NONE

PRIVATE

PUBLIC :: lfricinp_reorder_snow_field

CONTAINS

SUBROUTINE lfricinp_reorder_snow_field(field, um_grid)
! Take a snow layer field and converts from UM to
! LFRic ordering
USE lfricinp_grid_type_mod, ONLY: lfricinp_grid_type

USE log_mod, ONLY: log_event, LOG_LEVEL_ERROR, log_scratch_space, &
    LOG_LEVEL_INFO
! Arguments
REAL(KIND=real64), ALLOCATABLE, INTENT(IN OUT) :: field(:,:)
! Note 1st index is 2D field, 2nd index is pseudo level number
TYPE(lfricinp_grid_type), INTENT(IN):: um_grid

! Local variables
REAL(KIND=real64), ALLOCATABLE :: field_temp(:,:)
INTEGER(KIND=int64) :: total_lev
INTEGER(KIND=int64) :: field_size_2d
INTEGER(KIND=int64) :: tile_num
INTEGER(KIND=int64) :: layer_num, i
INTEGER(KIND=int64) :: pseud_lev_count

total_lev = SIZE(field, 2)
field_size_2d = SIZE(field, 1)

IF ( total_lev /= um_grid%num_snow_layers * um_grid%num_surface_types ) THEN
  WRITE(log_scratch_space, '(2(A,I0))') "Mismatch between total_lev = ", &
       total_lev, " and snow_layers * surface types = ", &
       um_grid%num_snow_layers * um_grid%num_surface_types
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

ALLOCATE(field_temp(field_size_2d, total_lev))

! UM ordering:
!   All tiles for snow layer 1, followed by all tiles for snow
!   layer 2, then all tiles for snow layer 3. Ie the snow layers
!   for a given tile type are not contiguous
! LFRic ordering:
!   Tile 1, snow layer 1,2,3. Tile 2, snow layer 1,2,3. Ie the
!   snow layers are grouped together
pseud_lev_count = 1
DO tile_num = 1, um_grid%num_surface_types
  DO layer_num = 1, um_grid%num_snow_layers
    DO i = 1, field_size_2d
      field_temp(i, pseud_lev_count) =  field(i, tile_num + &
           ((layer_num - 1) * um_grid%num_surface_types))
    END DO
    pseud_lev_count = pseud_lev_count + 1
  END DO
END DO

CALL MOVE_ALLOC(field_temp, field)

END SUBROUTINE lfricinp_reorder_snow_field

END MODULE lfricinp_reorder_snow_field_mod
