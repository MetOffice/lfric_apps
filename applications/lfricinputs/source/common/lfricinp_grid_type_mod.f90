! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_grid_type_mod
! NEED to REMEMBER to rename THIS file AS WELL AS THE module!!! - do AFTER
! CODE REVIEW
! Intrinsic modules
use, intrinsic :: iso_fortran_env, only: real64, int64

! LFRic modules
use log_mod,  only: log_event, LOG_LEVEL_ERROR, log_scratch_space
use constants_mod,    only: imdi, rmdi

implicit none

private

! Type to contain information relating to the structured grid
type, public :: lfricinp_grid_type
  ! Horizontal grid
  integer(kind=int64) :: num_arakawa_cells_x = imdi
  integer(kind=int64) :: num_arakawa_cells_y = imdi
  integer(kind=int64) :: num_p_points_x = imdi
  integer(kind=int64) :: num_p_points_y = imdi
  integer(kind=int64) :: num_u_points_x = imdi
  integer(kind=int64) :: num_u_points_y = imdi
  integer(kind=int64) :: num_v_points_x = imdi
  integer(kind=int64) :: num_v_points_y = imdi
  integer(kind=int64) :: num_snow_layers = imdi
  integer(kind=int64) :: num_surface_types = imdi
  integer(kind=int64) :: num_ice_cats = imdi
  real(kind=real64)   :: spacing_x = rmdi
  real(kind=real64)   :: spacing_y = rmdi
  real(kind=real64)   :: grid_origin_x = rmdi
  real(kind=real64)   :: grid_origin_y = rmdi
  real(kind=real64)   :: p_origin_x = rmdi
  real(kind=real64)   :: p_origin_y = rmdi
  real(kind=real64)   :: u_origin_x = rmdi
  real(kind=real64)   :: u_origin_y = rmdi
  real(kind=real64)   :: v_origin_x = rmdi
  real(kind=real64)   :: v_origin_y = rmdi
  ! Hardwired parameters for now
  logical             :: rotated_pole = .false.
  real(kind=real64)   :: pole_lat = 90.0_real64
  real(kind=real64)   :: pole_long = 0.0_real64
  integer(kind=int64) :: horiz_grid_type = 0 ! Global
  contains
  procedure :: print_grid_coords
  procedure :: set_grid_coords
end type lfricinp_grid_type

interface lfricinp_grid_type
  module procedure lfricinp_grid_info_constructor_shumlib
end interface

contains

!-----------------------------------------------------------------

function lfricinp_grid_info_constructor_shumlib(um_input_file) result (self)
! Description:
!  Creates the lfricinp_grid_type object based on the information
!  passed in from a UM file.

! Shumlib modules
use f_shum_file_mod,       only: shum_file_type
use f_shum_fixed_length_header_indices_mod,       &
                           only: horiz_grid_type, &
                                 grid_staggering
use f_shum_fieldsfile_mod, only: f_shum_fixed_length_header_len

! lfricinputs modules
use lfricinp_check_shumlib_status_mod, &
                           only: shumlib

! DEPENDS ON: c_shum_byteswap.o
! This is required to force fcm-make to compile the C code; whilst the built-in
! dependency analyser successfully works out that it needs to compile the
! Fortran side of the byte-swapping code, it requires an explicit statement
! to force it to compile the C part of the byte-swapping code. This is
! currently the approved way of linking Fortran and C in fcm-make.

implicit none

! Arguments
type(shum_file_type), intent(INOUT) :: um_input_file
!
type(lfricinp_grid_type) :: self

! Local variables
! Grid sizes - integer header
integer(kind=int64), parameter :: ih_num_p_points_x = 6
integer(kind=int64), parameter :: ih_num_p_points_y = 7
! Grid sizes - real header
integer(kind=int64), parameter :: rh_grid_spacing_x = 1
integer(kind=int64), parameter :: rh_grid_spacing_y = 2
integer(kind=int64), parameter :: rh_grid_origin_y_coord = 3
integer(kind=int64), parameter :: rh_grid_origin_x_coord = 4

! Fixed length header
integer(kind=int64) :: um_file_fixed_length_header(      &
                            f_shum_fixed_length_header_len)
! Integer constants
integer(kind=int64), allocatable  :: um_file_integer_constants(:)
! Real constants
real(kind=real64), allocatable :: um_file_real_constants(:)

character(len=*), parameter :: routinename = &
     'lfricinp_grid_info_constructor_shumlib'

! Get fixed length header
call shumlib(routinename//'::get_fixed_length_header', &
     um_input_file % get_fixed_length_header(um_file_fixed_length_header))

! Get integer constants
call shumlib(routinename//'::get_integer_constants', &
     um_input_file % get_integer_constants(um_file_integer_constants))

! Get real constants
call shumlib(routinename//'::get_real_constants', &
      um_input_file % get_real_constants(um_file_real_constants))

call self%set_grid_coords(                                              &
       grid_staggering = um_file_fixed_length_header(grid_staggering),&
       num_p_points_x = um_file_integer_constants(ih_num_p_points_x), &
       num_p_points_y = um_file_integer_constants(ih_num_p_points_y), &
       grid_spacing_x = um_file_real_constants(rh_grid_spacing_x),    &
       grid_spacing_y = um_file_real_constants(rh_grid_spacing_y),    &
       grid_origin_x =  um_file_real_constants(rh_grid_origin_x_coord),&
       grid_origin_y =  um_file_real_constants(rh_grid_origin_y_coord) )

end function lfricinp_grid_info_constructor_shumlib

!------------------------------------------------------------------

subroutine set_grid_coords(self,                                 &
                         grid_staggering,  num_p_points_x,     &
                         num_p_points_y,   grid_spacing_x,     &
                         grid_spacing_y,   grid_origin_x,      &
                         grid_origin_y )

! Description: Set grid info from argument list

implicit none

class(lfricinp_grid_type), intent(INOUT) :: self

! Arguments
integer(kind=int64), intent(in) :: grid_staggering
integer(kind=int64), intent(in) :: num_p_points_x
integer(kind=int64), intent(in) :: num_p_points_y
real(kind=real64), intent(in)   :: grid_spacing_x
real(kind=real64), intent(in)   :: grid_spacing_y
real(kind=real64), intent(in)   :: grid_origin_x
real(kind=real64), intent(in)   :: grid_origin_y

! Grid staggering indicator values
integer(kind=int64), parameter :: arakawa_C_endgame = 6
integer(kind=int64), parameter :: arakawa_C_nd = 3

! Calculate the number of Arakawa cells that make up the grid.
if (grid_staggering == arakawa_C_endgame) then
  !      EG 2x2 grid example
  !
  !      .___V___.___V___.
  !      |       |       |    P points are on stagger location CENTER
  !      U___P___U___P___|
  !      |       |       |    U points are on stagger location EDGE1
  !      .___V___.___V___.
  !      |       |       |    V points are on stagger location EDGE2
  !      U___P___U___P___|
  !      |       |       |    o - Grid origin
  !      o___V___.___V___.
  !
  ! If the stagger adjusted Arakawa C with v at poles is used then the
  ! number of P points is equal to to number of cells
  self%num_arakawa_cells_x = num_p_points_x
  self%num_arakawa_cells_y = num_p_points_y
  ! Set number of p points
  self%num_p_points_x = num_p_points_x
  self%num_p_points_y = num_p_points_y
  ! Set number of u points - same as p points in both directions
  self%num_u_points_x = self%num_p_points_x
  self%num_u_points_y = self%num_p_points_y
  ! Set number of v points - same as p points in x, one extra in y
  self%num_v_points_x = self%num_p_points_x
  self%num_v_points_y = self%num_p_points_y + 1
  ! Set grid spacing
  self%spacing_x = grid_spacing_x
  self%spacing_y = grid_spacing_y
  ! Set base grid origin
  self%grid_origin_x = grid_origin_x
  self%grid_origin_y = grid_origin_y
  ! Set origin points for p, u and v points
  self%p_origin_x = grid_origin_x + (0.5_real64 * self%spacing_x)
  self%p_origin_y = grid_origin_y + (0.5_real64 * self%spacing_y)
  self%u_origin_x = grid_origin_x
  self%u_origin_y = grid_origin_y + (0.5_real64 * self%spacing_y)
  self%v_origin_x = grid_origin_x + (0.5_real64 * self%spacing_x)
  self%v_origin_y = grid_origin_y
else if(grid_staggering == arakawa_C_nd) then
! If a standard Arakawa C grid is used then the x direction the number of P
! points is equal to the number of cells but the y direction the number of
! P points is one greater than the number of cells.
  call log_event("New Dynamics Arakawa C grid not supported", LOG_LEVEL_ERROR)
else
  write(log_scratch_space, '(A,I0)') "Unsupported value, grid staggering = ", &
       grid_staggering
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

end subroutine set_grid_coords

!------------------------------------------------------------------

subroutine print_grid_coords(self)
! Description:
!  Prints out contents of lfricinp_grid_type
use log_mod,       only: log_event, LOG_LEVEL_INFO
implicit none
class(lfricinp_grid_type), intent(in) :: self
write(log_scratch_space, '(A)') "|============ Grid Diagnostics ============|"
call log_event(log_scratch_space, LOG_LEVEL_INFO)
write(log_scratch_space, '(A,I0)') &
     "|  num_arakawa_cells_x: ", self%num_arakawa_cells_x
call log_event(log_scratch_space, LOG_LEVEL_INFO)
write(log_scratch_space, '(A,I0)') &
     "|  num_arakawa_cells_y: ", self%num_arakawa_cells_y
call log_event(log_scratch_space, LOG_LEVEL_INFO)
write(log_scratch_space, '(A,I0)') "|  num_p_points_x: ", self%num_p_points_x
call log_event(log_scratch_space, LOG_LEVEL_INFO)
write(log_scratch_space, '(A,I0)') "|  num_p_points_y: ", self%num_p_points_y
call log_event(log_scratch_space, LOG_LEVEL_INFO)
write(log_scratch_space, '(A,I0)') "|  num_u_points_x: ", self%num_u_points_x
call log_event(log_scratch_space, LOG_LEVEL_INFO)
write(log_scratch_space, '(A,I0)') "|  num_u_points_y: ", self%num_u_points_y
call log_event(log_scratch_space, LOG_LEVEL_INFO)
write(log_scratch_space, '(A,I0)') "|  num_v_points_x: ", self%num_v_points_x
call log_event(log_scratch_space, LOG_LEVEL_INFO)
write(log_scratch_space, '(A,I0)') "|  num_v_points_y: ", self%num_v_points_y
call log_event(log_scratch_space, LOG_LEVEL_INFO)
write(log_scratch_space, '(A,F0.6)') "|  spacing_x: ", self%spacing_x
call log_event(log_scratch_space, LOG_LEVEL_INFO)
write(log_scratch_space, '(A,F0.6)') "|  spacing_y: ", self%spacing_y
call log_event(log_scratch_space, LOG_LEVEL_INFO)
write(log_scratch_space, '(A,F0.6)') "|  p_origin_x: ", self%p_origin_x
call log_event(log_scratch_space, LOG_LEVEL_INFO)
write(log_scratch_space, '(A,F0.6)') "|  p_origin_y: ", self%p_origin_y
call log_event(log_scratch_space, LOG_LEVEL_INFO)
write(log_scratch_space, '(A,F0.6)') "|  u_origin_x: ", self%u_origin_x
call log_event(log_scratch_space, LOG_LEVEL_INFO)
write(log_scratch_space, '(A,F0.6)') "|  u_origin_y: ", self%u_origin_y
call log_event(log_scratch_space, LOG_LEVEL_INFO)
write(log_scratch_space, '(A,F0.6)') "|  v_origin_x: ", self%v_origin_x
call log_event(log_scratch_space, LOG_LEVEL_INFO)
write(log_scratch_space, '(A,F0.6)') "|  v_origin_y: ", self%v_origin_y
call log_event(log_scratch_space, LOG_LEVEL_INFO)
write(log_scratch_space, '(A)')   "|==========================================|"
call log_event(log_scratch_space, LOG_LEVEL_INFO)

end subroutine print_grid_coords

!------------------------------------------------------------------

end module lfricinp_grid_type_mod
