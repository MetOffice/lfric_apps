! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_grid_type_mod
! NEED TO REMEMBER TO RENAME THIS FILE AS WELL AS THE MODULE!!! - DO AFTER
! CODE REVIEW
! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY: real64, int64

! LFRic modules
USE log_mod,  ONLY: log_event, LOG_LEVEL_ERROR, log_scratch_space
USE constants_mod,    ONLY: imdi, rmdi

IMPLICIT NONE

PRIVATE

! Type to contain information relating to the structured grid
TYPE, PUBLIC :: lfricinp_grid_type
  ! Horizontal grid
  INTEGER(KIND=int64) :: num_arakawa_cells_x = imdi
  INTEGER(KIND=int64) :: num_arakawa_cells_y = imdi
  INTEGER(KIND=int64) :: num_p_points_x = imdi
  INTEGER(KIND=int64) :: num_p_points_y = imdi
  INTEGER(KIND=int64) :: num_u_points_x = imdi
  INTEGER(KIND=int64) :: num_u_points_y = imdi
  INTEGER(KIND=int64) :: num_v_points_x = imdi
  INTEGER(KIND=int64) :: num_v_points_y = imdi
  INTEGER(KIND=int64) :: num_snow_layers = imdi
  INTEGER(KIND=int64) :: num_surface_types = imdi
  REAL(KIND=real64)   :: spacing_x = rmdi
  REAL(KIND=real64)   :: spacing_y = rmdi
  REAL(KIND=real64)   :: grid_origin_x = rmdi
  REAL(KIND=real64)   :: grid_origin_y = rmdi
  REAL(KIND=real64)   :: p_origin_x = rmdi
  REAL(KIND=real64)   :: p_origin_y = rmdi
  REAL(KIND=real64)   :: u_origin_x = rmdi
  REAL(KIND=real64)   :: u_origin_y = rmdi
  REAL(KIND=real64)   :: v_origin_x = rmdi
  REAL(KIND=real64)   :: v_origin_y = rmdi
  ! Hardwired parameters for now
  LOGICAL             :: rotated_pole = .FALSE.
  REAL(KIND=real64)   :: pole_lat = 90.0
  REAL(KIND=real64)   :: pole_long = 0.0
  INTEGER(KIND=int64) :: horiz_grid_type = 0 ! Global
  CONTAINS
  PROCEDURE :: print_grid_coords
  PROCEDURE :: set_grid_coords
END TYPE lfricinp_grid_type

INTERFACE lfricinp_grid_type
  MODULE PROCEDURE lfricinp_grid_info_constructor_shumlib
END INTERFACE

CONTAINS

!-----------------------------------------------------------------

FUNCTION lfricinp_grid_info_constructor_shumlib(um_input_file) RESULT (self)
! Description:
!  Creates the lfricinp_grid_type object based on the information
!  passed in from a UM file.

! Shumlib modules
USE f_shum_file_mod,   ONLY: shum_file_type
USE f_shum_fixed_length_header_indices_mod, ONLY: &
    horiz_grid_type, grid_staggering
USE f_shum_fieldsfile_mod, ONLY:f_shum_fixed_length_header_len

! lfricinputs modules
USE lfricinp_check_shumlib_status_mod, ONLY: shumlib

! DEPENDS ON: c_shum_byteswap.o
! This is required to force fcm-make to compile the C code; whilst the built-in
! dependency analyser successfully works out that it needs to compile the
! Fortran side of the byte-swapping code, it requires an explicit statement
! to force it to compile the C part of the byte-swapping code. This is
! currently the approved way of linking Fortran and C in fcm-make.

IMPLICIT NONE

! Arguments
TYPE(shum_file_type), INTENT(INOUT) :: um_input_file
!
TYPE(lfricinp_grid_type) :: self

! Local variables
! Grid sizes - integer header
INTEGER(KIND=int64), PARAMETER :: ih_num_p_points_x = 6
INTEGER(KIND=int64), PARAMETER :: ih_num_p_points_y = 7
! Grid sizes - real header
INTEGER(KIND=int64), PARAMETER :: rh_grid_spacing_x = 1
INTEGER(KIND=int64), PARAMETER :: rh_grid_spacing_y = 2
INTEGER(KIND=int64), PARAMETER :: rh_grid_origin_y_coord = 3
INTEGER(KIND=int64), PARAMETER :: rh_grid_origin_x_coord = 4

! Fixed length header
INTEGER(KIND=int64) :: um_file_fixed_length_header(      &
                            f_shum_fixed_length_header_len)
! Integer constants
INTEGER(KIND=int64), ALLOCATABLE  :: um_file_integer_constants(:)
! Real constants
REAL(KIND=real64), ALLOCATABLE :: um_file_real_constants(:)

CHARACTER(LEN=*), PARAMETER :: routinename = &
     'lfricinp_grid_info_constructor_shumlib'

! Get fixed length header
CALL shumlib(routinename//'::get_fixed_length_header', &
     um_input_file % get_fixed_length_header(um_file_fixed_length_header))

! Get integer constants
CALL shumlib(routinename//'::get_integer_constants', &
     um_input_file % get_integer_constants(um_file_integer_constants))

! Get real constants
CALL shumlib(routinename//'::get_real_constants', &
      um_input_file % get_real_constants(um_file_real_constants))

CALL self%set_grid_coords(                                              &
       grid_staggering = um_file_fixed_length_header(grid_staggering),&
       num_p_points_x = um_file_integer_constants(ih_num_p_points_x), &
       num_p_points_y = um_file_integer_constants(ih_num_p_points_y), &
       grid_spacing_x = um_file_real_constants(rh_grid_spacing_x),    &
       grid_spacing_y = um_file_real_constants(rh_grid_spacing_y),    &
       grid_origin_x =  um_file_real_constants(rh_grid_origin_x_coord),&
       grid_origin_y =  um_file_real_constants(rh_grid_origin_y_coord) )

END FUNCTION lfricinp_grid_info_constructor_shumlib

!------------------------------------------------------------------

SUBROUTINE set_grid_coords(self,                                 &
                         grid_staggering,  num_p_points_x,     &
                         num_p_points_y,   grid_spacing_x,     &
                         grid_spacing_y,   grid_origin_x,      &
                         grid_origin_y )

! Description: Set grid info from argument list

IMPLICIT NONE

CLASS(lfricinp_grid_type), INTENT(INOUT) :: self

! Arguments
INTEGER(KIND=int64), INTENT(IN) :: grid_staggering
INTEGER(KIND=int64), INTENT(IN) :: num_p_points_x
INTEGER(KIND=int64), INTENT(IN) :: num_p_points_y
REAL(KIND=real64), INTENT(IN)   :: grid_spacing_x
REAL(KIND=real64), INTENT(IN)   :: grid_spacing_y
REAL(KIND=real64), INTENT(IN)   :: grid_origin_x
REAL(KIND=real64), INTENT(IN)   :: grid_origin_y

! Grid staggering indicator values
INTEGER(KIND=int64), PARAMETER :: arakawa_C_endgame = 6
INTEGER(KIND=int64), PARAMETER :: arakawa_C_nd = 3

! Calculate the number of Arakawa cells that make up the grid.
IF (grid_staggering == arakawa_C_endgame) THEN
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
  self%p_origin_x = grid_origin_x + (0.5 * self%spacing_x)
  self%p_origin_y = grid_origin_y + (0.5 * self%spacing_y)
  self%u_origin_x = grid_origin_x
  self%u_origin_y = grid_origin_y + (0.5 * self%spacing_y)
  self%v_origin_x = grid_origin_x + (0.5 * self%spacing_x)
  self%v_origin_y = grid_origin_y
ELSE IF(grid_staggering == arakawa_C_nd) THEN
! If a standard Arakawa C grid is used then the x direction the number of P
! points is equal to the number of cells but the y direction the number of
! P points is one greater than the number of cells.
  CALL log_event("New Dynamics Arakawa C grid not supported", LOG_LEVEL_ERROR)
ELSE
  WRITE(log_scratch_space, '(A,I0)') "Unsupported value, grid staggering = ", &
       grid_staggering
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

END SUBROUTINE set_grid_coords

!------------------------------------------------------------------

SUBROUTINE print_grid_coords(self)
! Description:
!  Prints out contents of lfricinp_grid_type
USE log_mod,       ONLY: log_event, LOG_LEVEL_INFO
IMPLICIT NONE
CLASS(lfricinp_grid_type), INTENT(IN) :: self
WRITE(log_scratch_space, '(A)') "|============ Grid Diagnostics ============|"
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
WRITE(log_scratch_space, '(A,I0)') &
     "|  num_arakawa_cells_x: ", self%num_arakawa_cells_x
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
WRITE(log_scratch_space, '(A,I0)') &
     "|  num_arakawa_cells_y: ", self%num_arakawa_cells_y
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
WRITE(log_scratch_space, '(A,I0)') "|  num_p_points_x: ", self%num_p_points_x
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
WRITE(log_scratch_space, '(A,I0)') "|  num_p_points_y: ", self%num_p_points_y
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
WRITE(log_scratch_space, '(A,I0)') "|  num_u_points_x: ", self%num_u_points_x
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
WRITE(log_scratch_space, '(A,I0)') "|  num_u_points_y: ", self%num_u_points_y
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
WRITE(log_scratch_space, '(A,I0)') "|  num_v_points_x: ", self%num_v_points_x
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
WRITE(log_scratch_space, '(A,I0)') "|  num_v_points_y: ", self%num_v_points_y
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
WRITE(log_scratch_space, '(A,F0.6)') "|  spacing_x: ", self%spacing_x
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
WRITE(log_scratch_space, '(A,F0.6)') "|  spacing_y: ", self%spacing_y
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
WRITE(log_scratch_space, '(A,F0.6)') "|  p_origin_x: ", self%p_origin_x
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
WRITE(log_scratch_space, '(A,F0.6)') "|  p_origin_y: ", self%p_origin_y
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
WRITE(log_scratch_space, '(A,F0.6)') "|  u_origin_x: ", self%u_origin_x
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
WRITE(log_scratch_space, '(A,F0.6)') "|  u_origin_y: ", self%u_origin_y
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
WRITE(log_scratch_space, '(A,F0.6)') "|  v_origin_x: ", self%v_origin_x
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
WRITE(log_scratch_space, '(A,F0.6)') "|  v_origin_y: ", self%v_origin_y
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
WRITE(log_scratch_space, '(A)')   "|==========================================|"
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

END SUBROUTINE print_grid_coords

!------------------------------------------------------------------

END MODULE lfricinp_grid_type_mod
