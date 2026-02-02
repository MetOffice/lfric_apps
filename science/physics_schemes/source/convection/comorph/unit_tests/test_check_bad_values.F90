! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

#if !defined(LFRIC)
! Unit test which checks whether the bad-value checker works on
! the current platform
program test_check_bad_values

use comorph_constants_mod, only: real_hmprec, real_cvprec,  name_length,       &
                     nx_full, ny_full, k_bot_conv, k_top_conv, k_top_init
use set_dependent_constants_mod, only: set_dependent_constants
use cmpr_type_mod, only: cmpr_type, cmpr_alloc

use check_bad_values_mod, only: check_bad_values_3d,                           &
                                check_bad_values_cmpr

implicit none

! Dimensions for the test array
integer, parameter :: nx = 10
integer, parameter :: ny = 20
integer, parameter :: nz = 5

! Character string input to check_bad_values
! (used to make error messages more informative)
character(len=name_length) :: where_string
character(len=name_length) :: field_name
! Flag to test for negative values.
logical :: l_positive

! Full 3-D field to test
real(kind=real_hmprec) :: field(nx,ny,nz)
! Compression list to test
real(kind=real_cvprec) :: field_cmpr(nx*ny)

! Compression indices for the compression list
type(cmpr_type) :: cmpr

! Lower and upper bounds of the test array
integer :: lb(3), ub(3)

real(kind=real_hmprec), parameter :: zero = 0.0_real_hmprec
real(kind=real_hmprec), parameter :: one  = 1.0_real_hmprec

integer :: i, j, k, ic

! Set array sizes in comorph_constants_mod to the hardwired values used here.
nx_full = nx
ny_full = ny
k_bot_conv = 1
k_top_conv = nz
k_top_init = nz-1
! Need to call this to set newline character for error messages
call set_dependent_constants

! Array lower-bounds are just 1 (no haloes or anything).
lb = [1,1,1]
ub = [nx,ny,nz]

! Initialise array to zero
do k = 1, nz
  do j = 1, ny
    do i = 1, nx
      field(i,j,k) = zero
    end do
  end do
end do

! Setup compression indices for test on compressed data
call cmpr_alloc( cmpr, nx*ny )
cmpr%n_points = nx*ny
do j = 1, ny
  do i = 1, nx
    ic = (j-1)*nx+i
    cmpr%index_i(ic) = i
    cmpr%index_j(ic) = j
  end do
end do

! Set other inputs to check_bad_values
where_string = "Unit test for check_bad_values"
l_positive = .true.


! Test whether routine can detect negative values
do k = 1, nz
  field(k,k,k) = -one
end do
field_name = "field_with_negative"
call check_bad_values_3d( lb, ub, field, where_string,                         &
                          field_name, l_positive )
k = 1
do ic = 1, cmpr%n_points
  field_cmpr(ic) = real(                                                       &
        field( cmpr%index_i(ic), cmpr%index_j(ic), k ), real_cvprec )
end do
call check_bad_values_cmpr( cmpr, k, field_cmpr, where_string,                 &
                            field_name, l_positive )

! Test whether routine can detect div-by-zero
do k = 1, nz
  field(k,k,k) = one / zero
end do
field_name = "field_with_div_by_zero"
call check_bad_values_3d( lb, ub, field, where_string,                         &
                          field_name, l_positive )
k = 1
do ic = 1, cmpr%n_points
  field_cmpr(ic) = real(                                                       &
        field( cmpr%index_i(ic), cmpr%index_j(ic), k ), real_cvprec )
end do
call check_bad_values_cmpr( cmpr, k, field_cmpr, where_string,                 &
                            field_name, l_positive )


! Test whether routine can detect SQRT(negative)
do k = 1, nz
  field(k,k,k) = sqrt(-one)
end do
field_name = "field_with_sqrt_negative"
call check_bad_values_3d( lb, ub, field, where_string,                         &
                          field_name, l_positive )
k = 1
do ic = 1, cmpr%n_points
  field_cmpr(ic) = real(                                                       &
        field( cmpr%index_i(ic), cmpr%index_j(ic), k ), real_cvprec )
end do
call check_bad_values_cmpr( cmpr, k, field_cmpr, where_string,                 &
                            field_name, l_positive )


end program test_check_bad_values
#endif
