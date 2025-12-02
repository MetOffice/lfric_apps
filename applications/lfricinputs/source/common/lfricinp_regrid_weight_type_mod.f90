! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!> @brief   Module for regridding weights TYPE
!> @details Holds data from a SCRIP style weights file, and proceedures for
!>          reading the file and regridding
module lfricinp_regrid_weights_type_mod

use lfricinp_um_parameters_mod, only: fnamelen, um_rmdi

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int32, int64, real64


implicit none

private

public :: lfricinp_regrid_weights_type

! Type to contain weights information read from weights file
type :: lfricinp_regrid_weights_type
  character(len=fnamelen) :: filename

  ! Weights order, ie num_wgts = 1, first order, =2 second order
  integer(kind=int32) :: num_wgts
  ! Number of links is the total number of connections between source
  ! and destination points
  integer(kind=int32) :: num_links
  ! Number of points on each grid/mesh
  integer(kind=int32) :: num_points_src
  integer(kind=int32) :: num_points_dst

  ! Remap matrix to contain weights indices
  real(kind=real64), allocatable :: remap_matrix(:,:)
  ! Source and destination addresses, used to link each remap matrix
  ! elements to a particular source and destination point
  integer(kind=int32), allocatable :: src_address(:)
  integer(kind=int32), allocatable :: src_address_2d(:,:)
  integer(kind=int32), allocatable :: dst_address(:)
  integer(kind=int32), allocatable :: dst_address_2d(:,:)
contains
  procedure :: load
  procedure :: populate_src_address_2d
  procedure :: populate_dst_address_2d
  procedure :: regrid_src_2d_dst_1d
  procedure :: regrid_src_1d_dst_2d
  procedure :: validate_src
  procedure :: validate_dst
end type lfricinp_regrid_weights_type

contains
!---------------------------------------------------------
! Start of type bound procedures
!---------------------------------------------------------
!> @brief Applies weights to perform 2d to 1d regridding from source to
!!        destination
!> @param[in]     src     2 dimensional source array
!> @param[inout]  dst     1 dimensional destination array
subroutine regrid_src_2d_dst_1d(self, src, dst)

! Uses real arrays, could be overloaded for different types,
! precision and shapes
real(kind=real64), intent(in) :: src(:,:)
real(kind=real64), intent(in out) :: dst(:)

class(lfricinp_regrid_weights_type) :: self

integer(kind=int32) :: l, w

dst(:) = 0.0_real64
! For now assume that any normalisation of weights has already been performed
! offline so there is no need to divide the destination by the fractional area
do w = 1, self%num_wgts
  do l = 1, self%num_links
    dst(self%dst_address(l)) = dst(self%dst_address(l)) + &
         (self%remap_matrix(w,l) * src(self%src_address_2d(l,1), &
                                       self%src_address_2d(l,2)))
  end do
end do

end subroutine regrid_src_2d_dst_1d

!---------------------------------------------------------
!> @brief Applies weights to perform 2d to 1d regridding from source to
!!        destination
!> @param[in]     src     1 dimensional source array
!> @param[inout]  dst     2 dimensional destination array
subroutine regrid_src_1d_dst_2d(self, src, dst)

! Uses real arrays, could be overloaded for different types,
! precision and shapes
real(kind=real64), intent(in) :: src(:)
real(kind=real64), intent(in out) :: dst(:,:)

class(lfricinp_regrid_weights_type) :: self

integer(kind=int32) :: l, w

dst(:,:) = 0.0_real64
! For now assume that any normalisation of weights has already been performed
! offline so there is no need to divide the destination by the fractional area
do w = 1, self%num_wgts
  do l = 1, self%num_links
    dst(self%dst_address_2d(l,1), self%dst_address_2d(l,2)) = &
         dst(self%dst_address_2d(l,1), self%dst_address_2d(l,2)) + &
         (self%remap_matrix(w,l) * src(self%src_address(l)))
  end do
end do

end subroutine regrid_src_1d_dst_2d

!---------------------------------------------------------
!> @brief    Loads a SCRIP style weights file into the weights derived type
!> @param[in]  filename  File to be loaded
subroutine load(self, filename)
! Description:
!  Loads a SCRIP style weights file into a weights derived type
! LFRic modules
use log_mod, only: log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR, &
                   log_scratch_space
! lfricinputs modules
use lfricinp_check_stat_ncdf_mod, only: check_stat_ncdf
! External libraries
use netcdf, only: NF90_OPEN, NF90_NOWRITE, NF90_GET_VAR, NF90_INQ_DIMID, &
                  NF90_INQUIRE_DIMENSION, NF90_INQ_VARID, NF90_CLOSE

implicit none

class(lfricinp_regrid_weights_type) :: self
character(len=fnamelen) :: filename

integer(kind=int32) :: dim_id
integer(kind=int32) :: var_id
integer(kind=int32) :: file_id

dim_id = -1
var_id = -1
file_id = -1

self % filename = trim(filename)
call log_event("Loading weights file: " // trim(self%filename), &
     LOG_LEVEL_INFO)

! Open weights file
call check_stat_ncdf( NF90_OPEN (self%filename, NF90_NOWRITE, &
     file_id))

! Get dimensions
call check_stat_ncdf( NF90_INQ_DIMID (file_id, "num_wgts", &
     dim_id))
call check_stat_ncdf( NF90_INQUIRE_DIMENSION (file_id, dim_id, &
     len=self%num_wgts))

! Number of links between source and destination grids/mesh
call check_stat_ncdf( NF90_INQ_DIMID (file_id, "num_links", &
     dim_id))
call check_stat_ncdf( NF90_INQUIRE_DIMENSION (file_id, dim_id, &
     len=self%num_links))

! Number of points in source grid/mesh
call check_stat_ncdf( NF90_INQ_DIMID (file_id, "src_grid_size", &
     dim_id))
call check_stat_ncdf( NF90_INQUIRE_DIMENSION (file_id, dim_id, &
     len=self%num_points_src))

! Number of points in destination grids/mesh
call check_stat_ncdf( NF90_INQ_DIMID (file_id, "dst_grid_size", &
     dim_id))
call check_stat_ncdf( NF90_INQUIRE_DIMENSION (file_id, dim_id, &
     len=self%num_points_dst))

! Allocate size of remap matrix
if (allocated(self%remap_matrix)) deallocate(self%remap_matrix)
allocate(self%remap_matrix( self%num_wgts, self%num_links))

write(log_scratch_space, '(2(A,I0),A)' ) "Allocated remap_matrix(num_wgts = ", &
     self%num_wgts, " , num_links = ", self%num_links, " )"
call log_event(log_scratch_space, LOG_LEVEL_INFO)

! Read remap matrix from file
call check_stat_ncdf( NF90_INQ_VARID (file_id, "remap_matrix", var_id))
call check_stat_ncdf( NF90_GET_VAR (file_id, var_id, self%remap_matrix(:,:)))

! Allocate source and destination address arrays
if (allocated(self%src_address)) deallocate(self%src_address)
allocate( self%src_address( self%num_links ))
if (allocated(self%dst_address)) deallocate(self%dst_address)
allocate( self%dst_address( self%num_links ))

! Read address arrays
call check_stat_ncdf( NF90_INQ_VARID (file_id, "src_address", var_id))
call check_stat_ncdf( NF90_GET_VAR (file_id, var_id, self%src_address))
call check_stat_ncdf( NF90_INQ_VARID (file_id, "dst_address", var_id))
call check_stat_ncdf( NF90_GET_VAR (file_id, var_id, self%dst_address))

! Close the file and release associated resources
call check_stat_ncdf( NF90_CLOSE (file_id))

end subroutine load

!---------------------------------------------------------
!> @brief     Convert 1d src_address into 2d src_address_2d
!> @param[in] nx Number of x coordinates
subroutine populate_src_address_2d (self, nx)
implicit none
! Convert the 1D addresses into 2d space
integer(kind=int32), intent(in) :: nx
class(lfricinp_regrid_weights_type) :: self

integer(kind=int32) :: l

if(allocated(self%src_address_2d)) deallocate(self%src_address_2d)
allocate(self%src_address_2d(self%num_links, 2))
do l = 1, self%num_links
  ! y coordinate
  self%src_address_2d(l,2) = &
       ((self%src_address(l) -1) / nx) + 1
  ! x coordinate
  self%src_address_2d(l,1) = self%src_address(l) - &
       nx*(self%src_address_2d(l,2) - 1)
end do

end subroutine populate_src_address_2d

!---------------------------------------------------------
!> @brief     Convert 1d dst_address into 2d dst_address_2d
!> @param[in] nx Number of x coordinates
subroutine populate_dst_address_2d (self, nx)
implicit none
! Convert the 1D addresses into 2d space
integer(kind=int32), intent(in) :: nx
class(lfricinp_regrid_weights_type) :: self

integer(kind=int32) :: l

if(allocated(self%dst_address_2d)) deallocate(self%dst_address_2d)
allocate(self%dst_address_2d(self%num_links, 2))
do l = 1, self%num_links
  ! y coordinate
  self%dst_address_2d(l,2) = &
       ((self%dst_address(l) -1) / nx) + 1
  ! x coordinate
  self%dst_address_2d(l,1) = self%dst_address(l) - &
       nx*(self%dst_address_2d(l,2) - 1)
end do

end subroutine populate_dst_address_2d

!---------------------------------------------------------
!> @brief Checks the number of points in the source grid against suppied value
!> @param[in] num_points_src The number of grid points expected in the source
subroutine validate_src(self, num_points_src)
use log_mod, only: log_scratch_space, LOG_LEVEL_ERROR, log_event

integer(kind=int32), intent(in) :: num_points_src

class(lfricinp_regrid_weights_type) :: self

if( num_points_src /= self%num_points_src) then
  write(log_scratch_space, '(2(A,I0))') &
       "Size mismatch between weights file and " // &
       "code grid definition source grid sizes. Num points weights file: ", &
       self%num_points_src, " Num points grid definition: ", num_points_src
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

end subroutine validate_src

!---------------------------------------------------------
!> @brief Checks the number of points in the destination grid against suppied value
!> @param[in] num_points_src The number of grid points expected in the destination
subroutine validate_dst(self, num_points_dst)
use log_mod, only: log_scratch_space, LOG_LEVEL_ERROR, log_event

integer(kind=int32), intent(in) :: num_points_dst

class(lfricinp_regrid_weights_type) :: self

if( num_points_dst /= self%num_points_dst) then
  write(log_scratch_space, '(2(A,I0))') &
       "Size mismatch between weights file and " // &
       "code grid definition source grid sizes. Num points weights file: ", &
       self%num_points_dst, " Num points grid definition: ", num_points_dst
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

end subroutine validate_dst
!---------------------------------------------------------
! End of type bound procedures
!---------------------------------------------------------

end module lfricinp_regrid_weights_type_mod
