#if ! defined(DRHOOK)
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module yomhook

!
! Description:
!   Dummy module to replace the DrHook library.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: dummy_libraries
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.1.

use parkind1, only: jpim, jprb
implicit none

logical, parameter :: lhook = .false.

interface dr_hook
  module procedure dr_hook
  module procedure dr_hook_comm
end interface dr_hook

contains

subroutine dr_hook(id,code,handle)
implicit none

!Arguments
character(len=*),   intent(in)    :: id
integer(kind=jpim), intent(in)    :: code
real(kind=jprb),    intent(in out) :: handle

!Nothing to do

return
end subroutine dr_hook

subroutine dr_hook_comm(id,code,handle,lcomm,comm)
implicit none

!Arguments
character(len=*),   intent(in)    :: id
integer(kind=jpim), intent(in)    :: code
logical(kind=jpim), intent(in)    :: lcomm
integer(kind=jpim), intent(in)    :: comm
real(kind=jprb),    intent(in out) :: handle

!Nothing to do

return
end subroutine dr_hook_comm

end module yomhook
#endif
