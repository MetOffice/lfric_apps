! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Purpose : To copy a single diagnostic field from secondary space to
!             the main data array for stash processing
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH
module copydiag_mod

! Module-wide use statements
use um_types, only: real_64, real_32

implicit none

character(len=*), parameter, private :: ModuleName='COPYDIAG_MOD'

! Set everything as private, and only expose the interfaces we want to
private
public :: copydiag

interface copydiag
  module procedure copydiag_dp_real,                                           &
                   copydiag_dp_real_vector,                                    &
                   !Add a single precision real version when required
                   copydiag_integer,                                           &
                   copydiag_logical
end interface

contains

!===============================================================================
subroutine copydiag_dp_real(diagout,diagin,row_length,rows)

use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook

implicit none

integer, intent(in) :: row_length, rows

real(kind=real_64), intent(in)  :: diagin(row_length,rows)

real(kind=real_64), intent(out) :: diagout(row_length * rows)

!Local variables
integer :: i,j,jpart

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='COPYDIAG_DP_REAL'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL do                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(none)                                                            &
!$OMP private(i,j,jpart)                                                       &
!$OMP SHARED(rows, row_length, diagout, diagin)
do j = 1, rows
  jpart = (j-1) * row_length
  do i = 1, row_length
    diagout(i+jpart) = diagin(i,j)
  end do
end do
!$OMP end PARALLEL do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine copydiag_dp_real

!===============================================================================
subroutine copydiag_dp_real_vector(diagout,diagin,length)

use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook

implicit none

integer, intent(in) :: length

real(kind=real_64), intent(in)  :: diagin(length)

real(kind=real_64), intent(out) :: diagout(length)

!Local variables
integer :: l

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='COPYDIAG_DP_REAL_VECTOR'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL do                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(none)                                                            &
!$OMP private(l)                                                               &
!$OMP SHARED(length, diagout, diagin)
do l = 1, length
  diagout(l) = diagin(l)
end do
!$OMP end PARALLEL do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine copydiag_dp_real_vector

!===============================================================================
subroutine copydiag_integer(diagout,diagin,row_length,rows)

use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook

implicit none

integer, intent(in) :: row_length, rows

integer, intent(in)  :: diagin(row_length,rows)

real(kind=real_64), intent(out) :: diagout(row_length * rows)

!Local variables
integer :: i,j,jpart

real, parameter :: real_mould = 0.0_real_64

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='COPYDIAG_INTEGER'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL do                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(none)                                                            &
!$OMP private(i,j,jpart)                                                       &
!$OMP SHARED(rows, row_length, diagout, diagin)
do j = 1, rows
  jpart = (j-1) * row_length
  do i = 1, row_length
    diagout(i+jpart) = transfer(diagin(i,j),real_mould)
  end do
end do
!$OMP end PARALLEL do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine copydiag_integer


!===============================================================================
subroutine copydiag_logical(diagout,diagin,row_length,rows)

use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook

implicit none

integer, intent(in) :: row_length, rows

logical, intent(in)  :: diagin(row_length,rows)

real(kind=real_64), intent(out) :: diagout(row_length * rows)

!Local variables
integer :: i,j,jpart

real, parameter :: real_mould = 0.0_real_64

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='COPYDIAG_LOGICAL'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL do                                                              &
!$OMP SCHEDULE(STATIC)                                                         &
!$OMP DEFAULT(none)                                                            &
!$OMP private(i,j,jpart)                                                       &
!$OMP SHARED(rows, row_length, diagout, diagin)
do j = 1, rows
  jpart = (j-1) * row_length
  do i = 1, row_length
    diagout(i+jpart) = transfer(diagin(i,j),real_mould)
  end do
end do
!$OMP end PARALLEL do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine copydiag_logical
end module copydiag_mod