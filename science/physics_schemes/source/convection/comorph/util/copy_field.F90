! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module copy_field_mod

implicit none

contains


! Subroutines to copy data between array pointers
! (which doesn't vectorise very efficiently);
! the fields can be passed into these routines, in which
! they are declared as normal arrays, allowing faster
! vectorisation.

!----------------------------------------------------------------
! Subroutine to copy a full 3-D field
!----------------------------------------------------------------
subroutine copy_field_3d( lb_in, ub_in, field_in, lb_out, ub_out, field_out )

use comorph_constants_mod, only: real_hmprec,                                  &
                     nx_full, ny_full, k_bot_conv, k_top_conv

implicit none

! Lower and upper bounds of the input array (in case it has halos)
integer, intent(in) :: lb_in(3), ub_in(3)
! Input array to be copied from
real(kind=real_hmprec), intent(in) :: field_in                                 &
              ( lb_in(1):ub_in(1),   lb_in(2):ub_in(2),   lb_in(3):ub_in(3)   )

! Lower and upper bounds of the output array (in case it has halos)
integer, intent(in) :: lb_out(3), ub_out(3)
! Output array to be copied into
real(kind=real_hmprec), intent(in out) :: field_out                            &
              ( lb_out(1):ub_out(1), lb_out(2):ub_out(2), lb_out(3):ub_out(3) )

! Loop counters
integer :: i, j, k

!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC) private( i, j, k )            &
!$OMP SHARED( nx_full, ny_full, k_bot_conv, k_top_conv, field_in, field_out )
do k = k_bot_conv, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full
      field_out(i,j,k) = field_in(i,j,k)
    end do
  end do
end do
!$OMP end PARALLEL do

return
end subroutine copy_field_3d


!----------------------------------------------------------------
! Subroutine to copy a full 2-D field only at the points
! specified in an input compression list.
!----------------------------------------------------------------
subroutine copy_field_cmpr( cmpr,                                              &
                            lb_in, ub_in, field_in, lb_out, ub_out, field_out )

use comorph_constants_mod, only: real_hmprec
use cmpr_type_mod, only: cmpr_type

implicit none

! Structure containing compression indices
type(cmpr_type), intent(in) :: cmpr

! Lower and upper bounds of the input array (in case it has halos)
integer, intent(in) :: lb_in(2), ub_in(2)
! Input array to be copied from
real(kind=real_hmprec), intent(in) :: field_in                                 &
                                   ( lb_in(1):ub_in(1),   lb_in(2):ub_in(2)   )

! Lower and upper bounds of the output array (in case it has halos)
integer, intent(in) :: lb_out(2), ub_out(2)
! Output array to be copied into
real(kind=real_hmprec), intent(in out) :: field_out                            &
                                   ( lb_out(1):ub_out(1), lb_out(2):ub_out(2) )
! Needs intent inout to preserve input values at points that
! aren't in the compression list.

! Loop counter
integer :: ic, i, j

do ic = 1, cmpr % n_points
  i = cmpr % index_i(ic)
  j = cmpr % index_j(ic)
  field_out(i,j) = field_in(i,j)
end do

return
end subroutine copy_field_cmpr


end module copy_field_mod
