! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_um_grid_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int64

! LFRic modules
USE log_mod,  ONLY: log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR

! UM2LFRic modules
USE lfricinp_grid_type_mod, ONLY: lfricinp_grid_type

IMPLICIT NONE

! Declare basic structured grid object
TYPE(lfricinp_grid_type), PUBLIC, SAVE :: um_grid

PRIVATE

PUBLIC :: lfricinp_set_grid_from_file

CONTAINS

!-----------------------------------------------------------

SUBROUTINE lfricinp_set_grid_from_file(um_input_file, num_snow_layers, &
                                        num_surface_types)
! Description:
!  Extracts grid information from UM input file to populate grid_info object

! Shumlib modules
USE f_shum_file_mod,   ONLY: shum_file_type

IMPLICIT NONE

! Arguments
TYPE(shum_file_type), INTENT(INOUT) :: um_input_file
INTEGER(KIND=int64), INTENT(IN) :: num_snow_layers
INTEGER(KIND=int64), INTENT(IN) :: num_surface_types

! Create UM grid info object, pass um_input_file to constructor
um_grid = lfricinp_grid_type(um_input_file)
CALL log_event("Printing grid information from UM input dump", LOG_LEVEL_INFO)
CALL um_grid%print_grid_coords()

! Set pseudo level information
um_grid%num_snow_layers = num_snow_layers
um_grid%num_surface_types = num_surface_types

END SUBROUTINE lfricinp_set_grid_from_file

END MODULE lfricinp_um_grid_mod
