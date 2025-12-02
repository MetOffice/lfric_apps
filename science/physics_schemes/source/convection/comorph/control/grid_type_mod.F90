! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module grid_type_mod

use comorph_constants_mod, only: real_cvprec, real_hmprec

implicit none

!----------------------------------------------------------------
! Type definition for a structure to group together the model
! grid heights and related fields
!----------------------------------------------------------------
type :: grid_type

  ! Heights of model-levels where thermodynamic fields are
  ! defined, relative to the surface
  real(kind=real_hmprec), pointer :: height_full(:,:,:) => null()
  ! These are taken to be the centre of the layer, on which
  ! increments are calculated and output

  ! Heights of half-levels, relative to the surface
  real(kind=real_hmprec), pointer :: height_half(:,:,:) => null()
  ! These are taken to be the layer interfaces, at the top
  ! and bottom of the layers:
  ! height_half(k) and height_half(k+1) are the lower and upper
  ! bounds of the layer centred at height_full(k).

  ! Height of surface above the centre of the Earth
  real(kind=real_hmprec), pointer :: r_surf(:,:) => null()
  ! Needed in the case where we have spherical coordinates and
  ! are not making the shallow atmosphere approximation;
  ! the layer masses are scaled by the changing grid area
  ! with height in order to yield conservation.

  ! Pressure on the same levels as thermodynamic fields
  ! (and height_full)
  real(kind=real_hmprec), pointer :: pressure_full(:,:,:)                      &
                                                       => null()

  ! Pressure at half-levels (same heights as height_half)
  real(kind=real_hmprec), pointer :: pressure_half(:,:,:)                      &
                                                       => null()

  ! Dry-density on the same levels as thermodynamic fields
  ! (and height_full)
  real(kind=real_hmprec), pointer :: rho_dry(:,:,:) => null()
  ! This is the dry-density used for the purposes of
  ! conservation.

end type grid_type


!----------------------------------------------------------------
! Addresses of compressed grid fields in a super-array used
! to pass them into the convection calculations
!----------------------------------------------------------------

! Number of fields in the super array
integer, parameter :: n_grid = 2

! Indices of each grid field in the super-array
integer, parameter :: i_height = 1
integer, parameter :: i_pressure = 2


contains


!----------------------------------------------------------------
! Subroutine to nullify the pointers in a grid_type structure
!----------------------------------------------------------------
subroutine grid_nullify( grid )
implicit none
type(grid_type), intent(in out) :: grid

grid % height_full   => null()
grid % height_half   => null()
grid % r_surf        => null()
grid % pressure_full => null()
grid % pressure_half => null()
grid % rho_dry       => null()

return
end subroutine grid_nullify


!----------------------------------------------------------------
! Subroutine to compress grid fields into a super-array to
! pass into the convection calculations
!----------------------------------------------------------------
subroutine grid_compress( grid, cmpr,                                          &
                          k, grid_k_super,                                     &
                          k_half, grid_half_super )

use cmpr_type_mod, only: cmpr_type
use compress_mod, only: compress

implicit none

! Structure containing pointers to the full 3-D grid arrays
type(grid_type), intent(in) :: grid

! Structure containing compression list info
type(cmpr_type), intent(in) :: cmpr

! Model-level k-index to compress from for full-level
integer, optional, intent(in) :: k
! Structure for output compressed copies at full-level
real(kind=real_cvprec), optional, intent(out) :: grid_k_super                  &
                                                 (:,:)

! Model-level k-index to compress from for half-level
integer, optional, intent(in) :: k_half
! Structure for output compressed copies at half-level
real(kind=real_cvprec), optional, intent(out) :: grid_half_super               &
                                                 (:,:)

! Bounds of the full arrays; needed for passing into
! compress to get indexing right when they have halos
integer :: lb(3), ub(3)

! Compress height and pressure at the current full level k
if ( present(k) .and. present(grid_k_super) ) then
  lb = lbound(grid % height_full)
  ub = ubound(grid % height_full)
  call compress( cmpr, lb(1:2), ub(1:2), grid % height_full(:,:,k),            &
                 grid_k_super(:,i_height) )
  lb = lbound(grid % pressure_full)
  ub = ubound(grid % pressure_full)
  call compress( cmpr, lb(1:2), ub(1:2), grid % pressure_full(:,:,k),          &
                 grid_k_super(:,i_pressure) )
end if

! Compress height and pressure at the next layer interface k_half
if ( present(k_half) .and. present(grid_half_super) ) then
  lb = lbound(grid % height_half)
  ub = ubound(grid % height_half)
  call compress( cmpr, lb(1:2), ub(1:2), grid % height_half(:,:,k_half),       &
                 grid_half_super(:,i_height) )
  lb = lbound(grid % pressure_half)
  ub = ubound(grid % pressure_half)
  call compress( cmpr, lb(1:2), ub(1:2), grid % pressure_half(:,:,k_half),     &
                 grid_half_super(:,i_pressure) )
end if

return
end subroutine grid_compress


!----------------------------------------------------------------
! Subroutine to check the input grid fields for bad values
!----------------------------------------------------------------
subroutine grid_check_bad_values( grid, where_string )

use comorph_constants_mod, only: name_length
use check_bad_values_mod, only: check_bad_values_3d

implicit none

type(grid_type), intent(in) :: grid
character(len=name_length), intent(in) :: where_string

character(len=name_length) :: field_name
integer :: lb(3), ub(3)
logical, parameter :: l_positive = .true.

field_name = "height_full"
lb = lbound( grid % height_full )
ub = ubound( grid % height_full )
call check_bad_values_3d( lb, ub, grid % height_full,                          &
                          where_string, field_name,                            &
                          l_positive )

field_name = "height_half"
lb = lbound( grid % height_half )
ub = ubound( grid % height_half )
call check_bad_values_3d( lb, ub, grid % height_half,                          &
                          where_string, field_name,                            &
                          l_positive, l_half=.true. )

field_name = "pressure_full"
lb = lbound( grid % pressure_full )
ub = ubound( grid % pressure_full )
call check_bad_values_3d( lb, ub, grid % pressure_full,                        &
                          where_string, field_name,                            &
                          l_positive )

field_name = "pressure_half"
lb = lbound( grid % pressure_half )
ub = ubound( grid % pressure_half )
call check_bad_values_3d( lb, ub, grid % pressure_half,                        &
                          where_string, field_name,                            &
                          l_positive, l_half=.true. )

field_name = "rho_dry"
lb = lbound( grid % rho_dry )
ub = ubound( grid % rho_dry )
call check_bad_values_3d( lb, ub, grid % rho_dry,                              &
                          where_string, field_name,                            &
                          l_positive )

! Note: no check on grid % r_surf yet; need to make a 2-D
! version of check_bad_values.

return
end subroutine grid_check_bad_values


end module grid_type_mod
