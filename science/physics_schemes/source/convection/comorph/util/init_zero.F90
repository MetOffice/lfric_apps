! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module init_zero_mod

implicit none

contains


! Subroutine to just initialise a 2-D field to zero.
! Used when initialising array pointers
! (which apparently doesn't vectorise very efficiently);
! the fields can be passed into this routine, in which
! they are declared as normal arrays, allowing faster
! vectorisation.
subroutine init_zero_2d( lb, ub, field )

use comorph_constants_mod, only: real_hmprec, nx_full, ny_full

implicit none

! Lower and upper bounds of the array (in case it has halos)
integer, intent(in) :: lb(2), ub(2)
! Array to be initialised to zero
real(kind=real_hmprec), intent(in out) :: field( lb(1):ub(1), lb(2):ub(2) )

real(kind=real_hmprec), parameter :: zero = 0.0_real_hmprec

! Loop counters
integer :: i, j

do j = 1, ny_full
  do i = 1, nx_full
    field(i,j) = zero
  end do
end do

return
end subroutine init_zero_2d


! As above, but for 3-D fields
subroutine init_zero_3d( lb, ub, field )

use comorph_constants_mod, only: real_hmprec,                                  &
                     nx_full, ny_full, k_bot_conv, k_top_conv

implicit none

! Lower and upper bounds of the array (in case it has halos)
integer, intent(in) :: lb(3), ub(3)
! Array to be initialised to zero
real(kind=real_hmprec), intent(in out) :: field                                &
                                      ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )

real(kind=real_hmprec), parameter :: zero = 0.0_real_hmprec

! Loop counters
integer :: i, j, k

!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC) private( i, j, k )            &
!$OMP SHARED( nx_full, ny_full, k_bot_conv, k_top_conv, field )
do k = k_bot_conv, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full
      field(i,j,k) = zero
    end do
  end do
end do
!$OMP end PARALLEL do

return
end subroutine init_zero_3d


! As above, but only at points on a compression list passed in
subroutine init_zero_cmpr( cmpr, lb, ub, field )

use comorph_constants_mod, only: real_hmprec
use cmpr_type_mod, only: cmpr_type

implicit none

! Compression indices of points where initialisation is needed
type(cmpr_type), intent(in) :: cmpr

! Lower and upper bounds of the array (in case it has halos)
integer, intent(in) :: lb(2), ub(2)
! Array to be initialised to zero
real(kind=real_hmprec), intent(in out) :: field( lb(1):ub(1), lb(2):ub(2) )

real(kind=real_hmprec), parameter :: zero = 0.0_real_hmprec

integer :: ic

do ic = 1, cmpr % n_points
  field( cmpr%index_i(ic), cmpr%index_j(ic) ) = zero
end do

return
end subroutine init_zero_cmpr


end module init_zero_mod
