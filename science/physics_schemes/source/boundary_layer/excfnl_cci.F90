! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  subroutine EXCFNL_CCI---------------------------------------------

!  Purpose: Compute Compressed Index

!  Code Description:
!    Language: FORTRAN 95
!    This code is written to UMDP3 programming standards.

!  Documentation: UMDP No.24

!  Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: boundary_layer
!--------------------------------------------------------------------
module excfnl_cci_mod

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName = 'EXCFNL_CCI_MOD'
contains

subroutine excfnl_cci (                                                        &
 c_len, to_do, ind_todo                                                        &
 )

use atm_fields_bounds_mod, only: pdims
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

integer, intent(in out) :: c_len
logical, intent(in out) :: to_do(pdims%i_end*pdims%j_end)
integer, intent(in out) :: ind_todo(pdims%i_end*pdims%j_end)

! local variables
integer                      :: n,m

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='EXCFNL_CCI'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
m = 0
!          Compress index for main loops
!CDIR Nodep
do n = 1, c_len
  if (to_do(n)) then
    m=m+1
    to_do(m)    = to_do(n)
    ind_todo(m) = ind_todo(n)
  end if
end do
c_len = m

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine excfnl_cci
end module excfnl_cci_mod