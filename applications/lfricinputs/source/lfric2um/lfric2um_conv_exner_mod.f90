! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfric2um_conv_exner_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only: int64, real64

! lfric2um modules
use lfricinp_lfric_driver_mod, only: local_rank
use lfricinp_um_grid_mod, only: um_grid

implicit none

private
public :: lfric2um_conv_p_exner, lfric2um_conv_exner_p

contains

subroutine lfric2um_conv_p_exner(p_to_exner)
! Description:
!  Convert pressure to Exner pressure
use planet_config_mod, only: p_zero, kappa

implicit none

real(kind=real64), dimension(:,:), intent(inout)  :: p_to_exner

integer(kind=int64) :: ix, iy, nx, ny
character(len=*), parameter :: routinename='lfric2um_conv_p_exner'

nx = um_grid%num_p_points_x
ny = um_grid%num_p_points_y

! Conversion is all done on the root process
if (local_rank == 0) then
  do ix = 1, nx
    do iy = 1, ny
      p_to_exner(ix, iy) = (p_to_exner(ix, iy) / p_zero) ** kappa
    end do
  end do
end if

end subroutine lfric2um_conv_p_exner

subroutine lfric2um_conv_exner_p(exner_to_p)
! Description:
!  Convert Exner pressure to pressure
use planet_config_mod, only: p_zero, one_over_kappa

implicit none

real(kind=real64), dimension(:,:), intent(inout)  :: exner_to_p

integer(kind=int64) :: ix, iy, nx, ny
character(len=*), parameter :: routinename='lfric2um_conv_exner_p'

nx = um_grid%num_p_points_x
ny = um_grid%num_p_points_y

! Conversion is all done on the root process
if (local_rank == 0) then
  do ix = 1, nx
    do iy = 1, ny
      exner_to_p(ix, iy) = p_zero * exner_to_p(ix, iy) ** one_over_kappa
    end do
  end do
end if

end subroutine lfric2um_conv_exner_p

end module lfric2um_conv_exner_mod
