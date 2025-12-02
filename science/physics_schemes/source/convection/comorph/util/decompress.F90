! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module decompress_mod

implicit none

contains


!----------------------------------------------------------------
! Subroutine to expand compressed array data into a
! full 2-D array
!----------------------------------------------------------------
! Doing this very trivial operation in a separate subroutine
! serves 2 purposes:
! a) Make code shorter and easier to follow.
! b) Often the full fields are declared as array pointers,
!    for which operations don't vectorise efficiently.  We
!    can pass them into this routine, where they are declared
!    as normal arrays, making the operation faster.
subroutine decompress( cmpr, field_cmpr, lb, ub, field_full )

use comorph_constants_mod, only: real_cvprec, real_hmprec
use cmpr_type_mod, only: cmpr_type

implicit none

! Structure containing compression indices
type(cmpr_type), intent(in) :: cmpr

! Input data to be expanded (has convection scheme precision)
real(kind=real_cvprec), intent(in) :: field_cmpr(cmpr%n_points)

! Lower and upper bounds of of the full 2-D array (in case it has halos)
integer, intent(in) :: lb(2), ub(2)
! Expanded output array (has host-model precision)
real(kind=real_hmprec), intent(in out) :: field_full                           &
                                          ( lb(1):ub(1), lb(2):ub(2) )
! Note: if not all the points in field_full are in the list of
! compressed data in field-cmpr, the input data in field_full
! at points not in the compression list will be retained on output.
! This is intended, so that e.g. when writing updated values of
! fields with convection increments added on at convecting points,
! the values at non-convecting points are left unaltered.

! Loop counter
integer :: ic

! Expand compressed data into the correct indices of the full
! field, converting to the host-model's precision.

! Directives tell the compiler that the stored i,j indices
! used in the loop are in memory-order (be it an uneven stride)
! so that some vectorisation is possible, which speeds things
! up slightly.
! For testing later !DIR$ IVDEP
! For testing later !DIR$ VECTOR ALWAYS
do ic = 1, cmpr%n_points
  field_full( cmpr%index_i(ic), cmpr%index_j(ic) )                             &
    = real( field_cmpr(ic), real_hmprec )
end do

return
end subroutine decompress


end module decompress_mod
