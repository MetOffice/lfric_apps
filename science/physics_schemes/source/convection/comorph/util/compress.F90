! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module compress_mod

implicit none

contains

!----------------------------------------------------------------
! Subroutine to compress a full 2-D field on to convecting points
!----------------------------------------------------------------
! Doing this very trivial operation in a separate subroutine
! serves 2 purposes:
! a) Make code shorter and easier to follow.
! b) Often the full fields are declared as array pointers,
!    for which operations don't vectorise efficiently.  We
!    can pass them into this routine, where they are declared
!    as normal arrays, making the operation faster.
subroutine compress( cmpr, lb, ub, field_full, field_cmpr )

use comorph_constants_mod, only: real_cvprec, real_hmprec
use cmpr_type_mod, only: cmpr_type

implicit none

! Structure containing compression indices
type(cmpr_type), intent(in) :: cmpr

! Lower and upper bounds of the full 2-D array (in case it has halos)
integer, intent(in) :: lb(2), ub(2)

! Full 2-D array to be compressed
real(kind=real_hmprec), intent(in) :: field_full                               &
                                      ( lb(1):ub(1), lb(2):ub(2) )

! Array to store compressed data
real(kind=real_cvprec), intent(out) :: field_cmpr(cmpr%n_points)

! Loop counter
integer :: ic

! Extract data at the required points,
! converting to the convection scheme's precision.

! Directives tell the compiler that the stored i,j indices
! used in the loop are in memory-order (be it an uneven stride)
! so that some vectorisation is possible, which speed things
! up slightly.
! for testing later !DIR$ IVDEP
! for testing later !DIR$ VECTOR ALWAYS
do ic = 1, cmpr%n_points
  field_cmpr(ic)                                                               &
    = real( field_full( cmpr%index_i(ic), cmpr%index_j(ic) ),                  &
            real_cvprec )
end do

return
end subroutine compress

end module compress_mod
