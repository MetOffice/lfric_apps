! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_um_grid_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int64

! LFRic modules
use log_mod,  only: log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR

! UM2LFRic modules
use lfricinp_grid_type_mod, only: lfricinp_grid_type

implicit none

! Declare basic structured grid object
type(lfricinp_grid_type), public, save :: um_grid

private

public :: lfricinp_set_grid_from_file

contains

!-----------------------------------------------------------

subroutine lfricinp_set_grid_from_file(um_input_file, num_snow_layers, &
                                        num_surface_types, num_ice_cats)
! Description:
!  Extracts grid information from UM input file to populate grid_info object

! Shumlib modules
use f_shum_file_mod,   only: shum_file_type

implicit none

! Arguments
type(shum_file_type), intent(INOUT) :: um_input_file
integer(kind=int64), intent(in) :: num_snow_layers
integer(kind=int64), intent(in) :: num_surface_types
integer(kind=int64), intent(in) :: num_ice_cats

! Create UM grid info object, pass um_input_file to constructor
um_grid = lfricinp_grid_type(um_input_file)
call log_event("Printing grid information from UM input dump", LOG_LEVEL_INFO)
call um_grid%print_grid_coords()

! Set pseudo level information
um_grid%num_snow_layers = num_snow_layers
um_grid%num_surface_types = num_surface_types
um_grid%num_ice_cats = num_ice_cats

end subroutine lfricinp_set_grid_from_file

end module lfricinp_um_grid_mod
