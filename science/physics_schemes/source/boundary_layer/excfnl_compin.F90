! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  subroutine EXCFNL_COMPIN------------------------------------------

!  Purpose: Compute condition for inner interation and
!           compress index array

!  Code Description:
!    Language: FORTRAN 95
!    This code is written to UMDP3 programming standards.

!  Documentation: UMDP No.24

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!--------------------------------------------------------------------
module excfnl_compin_mod

use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName = 'EXCFNL_COMPIN_MOD'
contains

subroutine excfnl_compin (                                                     &
! in fields
 up, wb_ratio, dec_thres, switch,                                              &
! INOUT fields
 c_len_i, ind_todo_i, todo_inner                                               &
 )

use atm_fields_bounds_mod, only: pdims
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

! Intent in:

integer, intent(in) :: up(pdims%i_end*pdims%j_end)
real(kind=r_bl),    intent(in) ::                                              &
 wb_ratio(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end),                &
                                  ! WBN_INT/WBP_INT
 dec_thres(pdims%i_start:pdims%i_end,pdims%j_start:pdims%j_end)
                                  ! Local decoupling threshold

integer, intent(in)           :: switch ! =1 for KSURF, 2 for KTOP

! Intent INOUT:

integer, intent(in out)        :: c_len_i
integer, intent(in out) :: ind_todo_i(pdims%i_end*pdims%j_end)
logical, intent(in out) :: todo_inner(pdims%i_end*pdims%j_end)

! local variables
integer                       :: n,m
integer                       :: l, i1, j1

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='EXCFNL_COMPIN'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!     Check for active elements
select case (switch)
case (1)
          ! For top of KSURF
  do n = 1, c_len_i
    l = ind_todo_i(n)
    j1=(l-1)/pdims%i_end+1
    i1=l-(j1-1)*pdims%i_end

    todo_inner(n) = (                                                          &
   (up(l) == 1 .and. wb_ratio(i1,j1) <  dec_thres(i1,j1)) .or.                 &
              ! keep working up while wb_ratio lt thres
   (up(l) == 0 .and. wb_ratio(i1,j1) >  dec_thres(i1,j1)) )
              ! keep working down while wb_ratio gt thres

  end do

case (2)
          ! For base of KTOP
  do n = 1, c_len_i
    l = ind_todo_i(n)
    j1=(l-1)/pdims%i_end+1
    i1=l-(j1-1)*pdims%i_end

    todo_inner(n) = (                                                          &
   (up(l) == 1 .and. wb_ratio(i1,j1) >  dec_thres(i1,j1)) .or.                 &
              ! keep working up while wb_ratio gt thres
   (up(l) == 0 .and. wb_ratio(i1,j1) <  dec_thres(i1,j1)) )
              ! keep working down while wb_ratio lt thres

  end do

end select

!     Compress

m = 0
!CDIR Nodep
do n = 1, c_len_i
  if (todo_inner(n)) then
    m = m+1
    todo_inner(m) = todo_inner(n)
    ind_todo_i(m) = ind_todo_i(n)
  end if
end do
c_len_i = m

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine excfnl_compin
end module excfnl_compin_mod
