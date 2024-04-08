! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_regrid_weights_type_mod

USE lfricinp_um_parameters_mod, ONLY: fnamelen, um_rmdi

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int32, int64, real64


IMPLICIT NONE

PRIVATE

PUBLIC :: lfricinp_regrid_weights_type

! Type to contain weights information read from weights file
TYPE :: lfricinp_regrid_weights_type
  CHARACTER(LEN=fnamelen) :: filename

  ! Weights order, ie num_wgts = 1, first order, =2 second order
  INTEGER(KIND=int32) :: num_wgts
  ! Number of links is the total number of connections between source
  ! and destination points
  INTEGER(KIND=int32) :: num_links
  ! Number of points on each grid/mesh
  INTEGER(KIND=int32) :: num_points_src
  INTEGER(KIND=int32) :: num_points_dst

  ! Remap matrix to contain weights indices
  REAL(KIND=real64), ALLOCATABLE :: remap_matrix(:,:)
  ! Source and destination addresses, used to link each remap matrix
  ! elements to a particular source and destination point
  INTEGER(KIND=int32), ALLOCATABLE :: src_address(:)
  INTEGER(KIND=int32), ALLOCATABLE :: src_address_2d(:,:)
  INTEGER(KIND=int32), ALLOCATABLE :: dst_address(:)
  INTEGER(KIND=int32), ALLOCATABLE :: dst_address_2d(:,:)
CONTAINS
  PROCEDURE :: load
  PROCEDURE :: populate_src_address_2d
  PROCEDURE :: populate_dst_address_2d
  PROCEDURE :: regrid_src_2d_dst_1d
  PROCEDURE :: regrid_src_1d_dst_2d
  PROCEDURE :: validate_src
  PROCEDURE :: validate_dst
END TYPE lfricinp_regrid_weights_type

CONTAINS
!---------------------------------------------------------
! Start of type bound procedures
!---------------------------------------------------------

SUBROUTINE regrid_src_2d_dst_1d(self, src, dst)

! Apply weights to perform regrid between source and destination
! Uses real arrays, could be overloaded for different types,
! precision and shapes
REAL(KIND=real64), INTENT(IN) :: src(:,:)
REAL(KIND=real64), INTENT(IN OUT) :: dst(:)

CLASS(lfricinp_regrid_weights_type) :: self

INTEGER(KIND=int32) :: l, w

dst(:) = 0.0
! For now assume that any normalisation of weights has already been performed
! offline so there is no need to divide the destination by the fractional area
DO w = 1, self%num_wgts
  DO l = 1, self%num_links
    dst(self%dst_address(l)) = dst(self%dst_address(l)) + &
         (self%remap_matrix(w,l) * src(self%src_address_2d(l,1), &
                                       self%src_address_2d(l,2)))
  END DO
END DO

END SUBROUTINE regrid_src_2d_dst_1d

!---------------------------------------------------------

SUBROUTINE regrid_src_1d_dst_2d(self, src, dst)

! Apply weights to perform copy between source and destination
! Uses real arrays, could be overloaded for different types,
! precision and shapes
REAL(KIND=real64), INTENT(IN) :: src(:)
REAL(KIND=real64), INTENT(IN OUT) :: dst(:,:)

CLASS(lfricinp_regrid_weights_type) :: self

INTEGER(KIND=int32) :: l, w

dst(:,:) = 0.0
! For now assume that any normalisation of weights has already been performed
! offline so there is no need to divide the destination by the fractional area
DO w = 1, self%num_wgts
  DO l = 1, self%num_links
    dst(self%dst_address_2d(l,1), self%dst_address_2d(l,2)) = &
         dst(self%dst_address_2d(l,1), self%dst_address_2d(l,2)) + &
         (self%remap_matrix(w,l) * src(self%src_address(l)))
  END DO
END DO

END SUBROUTINE regrid_src_1d_dst_2d

!---------------------------------------------------------

SUBROUTINE load(self, filename)
! Description:
!  Loads a SCRIP style weights file into a weights derived type
! LFRic modules
USE log_mod, ONLY: log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR, &
    log_scratch_space
! lfricinputs modules
USE lfricinp_check_stat_ncdf_mod, ONLY: check_stat_ncdf
! External libraries
USE netcdf, ONLY: NF90_OPEN, NF90_NOWRITE, NF90_GET_VAR, NF90_INQ_DIMID, &
     NF90_INQUIRE_DIMENSION, NF90_INQ_VARID, NF90_CLOSE

IMPLICIT NONE

CLASS(lfricinp_regrid_weights_type) :: self
CHARACTER(LEN=fnamelen) :: filename

INTEGER(KIND=int32) :: dim_id
INTEGER(KIND=int32) :: var_id
INTEGER(KIND=int32) :: file_id

dim_id = -1
var_id = -1
file_id = -1

self % filename = TRIM(filename)
CALL log_event("Loading weights file: " // TRIM(self%filename), &
     LOG_LEVEL_INFO)

! Open weights file
CALL check_stat_ncdf( NF90_OPEN (self%filename, NF90_NOWRITE, &
     file_id))

! Get dimensions
CALL check_stat_ncdf( NF90_INQ_DIMID (file_id, "num_wgts", &
     dim_id))
CALL check_stat_ncdf( NF90_INQUIRE_DIMENSION (file_id, dim_id, &
     len=self%num_wgts))

! Number of links between source and destination grids/mesh
CALL check_stat_ncdf( NF90_INQ_DIMID (file_id, "num_links", &
     dim_id))
CALL check_stat_ncdf( NF90_INQUIRE_DIMENSION (file_id, dim_id, &
     len=self%num_links))

! Number of points in source grid/mesh
CALL check_stat_ncdf( NF90_INQ_DIMID (file_id, "src_grid_size", &
     dim_id))
CALL check_stat_ncdf( NF90_INQUIRE_DIMENSION (file_id, dim_id, &
     len=self%num_points_src))

! Number of points in destination grids/mesh
CALL check_stat_ncdf( NF90_INQ_DIMID (file_id, "dst_grid_size", &
     dim_id))
CALL check_stat_ncdf( NF90_INQUIRE_DIMENSION (file_id, dim_id, &
     len=self%num_points_dst))

! Allocate size of remap matrix
IF (ALLOCATED(self%remap_matrix)) DEALLOCATE(self%remap_matrix)
ALLOCATE(self%remap_matrix( self%num_wgts, self%num_links))

WRITE(log_scratch_space, '(2(A,I0),A)' ) "Allocated remap_matrix(num_wgts = ", &
     self%num_wgts, " , num_links = ", self%num_links, " )"
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

! Read remap matrix from file
CALL check_stat_ncdf( NF90_INQ_VARID (file_id, "remap_matrix", var_id))
CALL check_stat_ncdf( NF90_GET_VAR (file_id, var_id, self%remap_matrix(:,:)))

! Allocate source and destination address arrays
IF (ALLOCATED(self%src_address)) DEALLOCATE(self%src_address)
ALLOCATE( self%src_address( self%num_links ))
IF (ALLOCATED(self%dst_address)) DEALLOCATE(self%dst_address)
ALLOCATE( self%dst_address( self%num_links ))

! Read address arrays
CALL check_stat_ncdf( NF90_INQ_VARID (file_id, "src_address", var_id))
CALL check_stat_ncdf( NF90_GET_VAR (file_id, var_id, self%src_address))
CALL check_stat_ncdf( NF90_INQ_VARID (file_id, "dst_address", var_id))
CALL check_stat_ncdf( NF90_GET_VAR (file_id, var_id, self%dst_address))

! Close the file and release associated resources
CALL check_stat_ncdf( NF90_CLOSE (file_id))

END SUBROUTINE load

!---------------------------------------------------------

SUBROUTINE populate_src_address_2d (self, nx)
IMPLICIT NONE
! Convert the 1D addresses into 2d space
INTEGER(KIND=int32), INTENT(IN) :: nx
CLASS(lfricinp_regrid_weights_type) :: self

INTEGER(KIND=int32) :: l

IF(ALLOCATED(self%src_address_2d)) DEALLOCATE(self%src_address_2d)
ALLOCATE(self%src_address_2d(self%num_links, 2))
DO l = 1, self%num_links
  ! y coordinate
  self%src_address_2d(l,2) = &
       ((self%src_address(l) -1) / nx) + 1
  ! x coordinate
  self%src_address_2d(l,1) = self%src_address(l) - &
       nx*(self%src_address_2d(l,2) - 1)
END DO

END SUBROUTINE populate_src_address_2d

!---------------------------------------------------------

SUBROUTINE populate_dst_address_2d (self, nx)
IMPLICIT NONE
! Convert the 1D addresses into 2d space
INTEGER(KIND=int32), INTENT(IN) :: nx
CLASS(lfricinp_regrid_weights_type) :: self

INTEGER(KIND=int32) :: l

IF(ALLOCATED(self%dst_address_2d)) DEALLOCATE(self%dst_address_2d)
ALLOCATE(self%dst_address_2d(self%num_links, 2))
DO l = 1, self%num_links
  ! y coordinate
  self%dst_address_2d(l,2) = &
       ((self%dst_address(l) -1) / nx) + 1
  ! x coordinate
  self%dst_address_2d(l,1) = self%dst_address(l) - &
       nx*(self%dst_address_2d(l,2) - 1)
END DO

END SUBROUTINE populate_dst_address_2d

!---------------------------------------------------------

SUBROUTINE validate_src(self, num_points_src)
USE log_mod, ONLY: log_scratch_space, LOG_LEVEL_ERROR, log_event

INTEGER(KIND=int32), INTENT(IN) :: num_points_src

CLASS(lfricinp_regrid_weights_type) :: self

IF( num_points_src /= self%num_points_src) THEN
  WRITE(log_scratch_space, '(2(A,I0))') &
       "Size mismatch between weights file and " // &
       "code grid definition source grid sizes. Num points weights file: ", &
       self%num_points_src, " Num points grid definition: ", num_points_src
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

END SUBROUTINE validate_src

!---------------------------------------------------------

SUBROUTINE validate_dst(self, num_points_dst)
USE log_mod, ONLY: log_scratch_space, LOG_LEVEL_ERROR, log_event

INTEGER(KIND=int32), INTENT(IN) :: num_points_dst

CLASS(lfricinp_regrid_weights_type) :: self

IF( num_points_dst /= self%num_points_dst) THEN
  WRITE(log_scratch_space, '(2(A,I0))') &
       "Size mismatch between weights file and " // &
       "code grid definition source grid sizes. Num points weights file: ", &
       self%num_points_dst, " Num points grid definition: ", num_points_dst
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

END SUBROUTINE validate_dst
!---------------------------------------------------------
! End of type bound procedures
!---------------------------------------------------------

END MODULE lfricinp_regrid_weights_type_mod
