! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module diff_field_mod

implicit none

contains

! Subroutine to subtract one full 2-D field from another,
! only at points specified in an input compression list.
! Used for calculating increment diagnostics by differencing
! saved values from earlier with the latest values.
! The field and diagnostic arrays are passed around as pointers,
! (for which calculations don't vectorise very efficienty);
! the fields can be passed into this routine, in which
! they are declared as normal arrays, allowing faster
! vectorisation.
subroutine diff_field_cmpr( cmpr,                                              &
                            lb_in, ub_in, field_in, lb_out, ub_out, field_out )

use comorph_constants_mod, only: real_hmprec
use cmpr_type_mod, only: cmpr_type

implicit none

! Structure containing compression indices
type(cmpr_type), intent(in) :: cmpr

! Lower and upper bounds of the input array (in case it has halos)
integer, intent(in) :: lb_in(2), ub_in(2)
! In: latest value of the field
real(kind=real_hmprec), intent(in) :: field_in                                 &
                                   ( lb_in(1):ub_in(1),   lb_in(2):ub_in(2)   )

! Lower and upper bounds of the array to be decremented
integer, intent(in) :: lb_out(2), ub_out(2)
! In: saved previous value of the field
! Out: increment = latest value minus previous value
real(kind=real_hmprec), intent(in out) :: field_out                            &
                                   ( lb_out(1):ub_out(1), lb_out(2):ub_out(2) )

! Loop counter
integer :: ic, i, j

do ic = 1, cmpr % n_points
  i = cmpr % index_i(ic)
  j = cmpr % index_j(ic)
  field_out(i,j) = field_in(i,j) - field_out(i,j)
end do

return
end subroutine diff_field_cmpr

end module diff_field_mod
