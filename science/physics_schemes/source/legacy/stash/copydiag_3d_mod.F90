! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   Subroutine COPYDIAG_3D -------------------------------------------
!
!   Purpose : To copy a diagnostic field from secondary space to the
!             main data array for stash processing, and to extend the
!             data to a full horizontal field. Input data of multilevel
!             fields is assumed to be on all model levels. Output data
!             is on the levels required.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH

module copydiag_3d_mod

! Module-wide use statements
use um_types, only: real_64, real_32

implicit none

character(len=*), parameter, private :: ModuleName='COPYDIAG_3D_MOD'

! Set everything as private, and only expose the interfaces we want to
private
public :: copydiag_3d

interface copydiag_3d
  module procedure copydiag_3d_dp_real,                                        &
                   copydiag_3d_sp_real,                                        &
                   copydiag_3d_logical,                                        &
                   copydiag_3d_integer
end interface

contains

subroutine copydiag_3d_dp_real(diagout,diagin,row_length,rows,levels,          &
     stlist,len_stlist,stash_levels,len_stashlevels)

use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook

use set_levels_list_mod, only: set_levels_list

implicit none

integer, intent(in) ::                                                         &
  row_length,                       & ! Number of points in a row
  rows,                             & ! Number of rows
  levels,                           & ! Number of levels in input data
  len_stlist,                                                                  &
  stlist(len_stlist),               & ! Stash list
  len_stashlevels,                                                             &
  stash_levels(len_stashlevels,levels) ! Stash levels list.

real(kind=real_64), intent(in) :: diagin(row_length, rows, levels)

real(kind=real_64), intent(out) :: diagout(row_length * rows * levels)

! Local variables
integer ::                                                                     &
  i,j,k,kout,nout,                                                             &
  kpart,jpart,                                                                 &
  levels_list(levels)         ! list of levels to be copied

logical ::                                                                     &
   list(levels)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='COPYDIAG_3D_DP_REAL'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call set_levels_list(levels,len_stlist,stlist,list,stash_levels,               &
      len_stashlevels)

!  Move data from DIAGIN to DIAGOUT at levels requested

nout = 0
do k = 1, levels
  if (list(k)) then
    nout = nout + 1
    levels_list(nout) = k
  end if
end do

!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED( nout, levels_list, rows, row_length, diagout, diagin)            &
!$OMP private ( i, j, k, kout, kpart, jpart )
!DIR$ IVDEP
do kout = 1, nout
  k     = levels_list(kout)
  kpart = (kout-1) * row_length * rows

  do j = 1, rows
    jpart = (j-1) * row_length
    do i = 1, row_length
      diagout(i+jpart+kpart) = diagin(i,j,k)
    end do ! i
  end do ! j
end do ! k
!$OMP end PARALLEL do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine copydiag_3d_dp_real

!===============================================================================

subroutine copydiag_3d_sp_real(diagout,diagin,row_length,rows,levels,          &
     stlist,len_stlist,stash_levels,len_stashlevels)

use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook

use set_levels_list_mod, only: set_levels_list

implicit none

integer, intent(in) ::                                                         &
  row_length,                       & ! Number of points in a row
  rows,                             & ! Number of rows
  levels,                           & ! Number of levels in input data
  len_stlist,                                                                  &
  stlist(len_stlist),               & ! Stash list
  len_stashlevels,                                                             &
  stash_levels(len_stashlevels,levels) ! Stash levels list.

real(kind=real_32), intent(in) :: diagin(row_length, rows, levels)

real(kind=real_64), intent(out) :: diagout(row_length * rows * levels)

! Local variables
integer ::                                                                     &
  i,j,k,kout,nout,                                                             &
  kpart,jpart,                                                                 &
  levels_list(levels)         ! list of levels to be copied

logical ::                                                                     &
   list(levels)

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='COPYDIAG_3D_SP_REAL'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call set_levels_list(levels,len_stlist,stlist,list,stash_levels,               &
      len_stashlevels)

!  Move data from DIAGIN to DIAGOUT at levels requested

nout = 0
do k = 1, levels
  if (list(k)) then
    nout = nout + 1
    levels_list(nout) = k
  end if
end do

!$OMP PARALLEL do DEFAULT(none)                                                &
!$OMP SHARED( nout, levels_list, rows, row_length, diagout, diagin)            &
!$OMP private ( i, j, k, kout, kpart, jpart )
!DIR$ IVDEP
do kout = 1, nout
  k     = levels_list(kout)
  kpart = (kout-1) * row_length * rows

  do j = 1, rows
    jpart = (j-1) * row_length
    do i = 1, row_length
      diagout(i+jpart+kpart) = real(diagin(i,j,k),kind=real_64)
    end do ! i
  end do ! j
end do ! k
!$OMP end PARALLEL do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine copydiag_3d_sp_real

!===============================================================================

subroutine copydiag_3d_logical(diagout,diagin,row_length,rows,levels,          &
     stlist,len_stlist,stash_levels,len_stashlevels)

use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook

use set_levels_list_mod, only: set_levels_list

implicit none

integer, intent(in) ::                                                         &
  row_length,                       & ! Number of points in a row
  rows,                             & ! Number of rows
  levels,                           & ! Number of levels in input data
  len_stlist,                                                                  &
  stlist(len_stlist),               & ! Stash list
  len_stashlevels,                                                             &
  stash_levels(len_stashlevels,levels) ! Stash levels list.

logical, intent(in) :: diagin(row_length,rows,levels)

real,intent(out) ::   diagout(row_length *  rows * levels)

! Local variables
integer ::                                                                     &
  i,j,k,kout,nout,                                                             &
  kpart,jpart,                                                                 &
  levels_list(levels)         ! list of levels to be copied

logical ::                                                                     &
   list(levels)

real, parameter :: real_mould = 0.0_real_64

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='COPYDIAG_3D_LOGICAL'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call set_levels_list(levels,len_stlist,stlist,list,stash_levels,               &
      len_stashlevels)

!  Move data from DIAGIN to DIAGOUT at levels requested

nout = 0
do k = 1, levels
  if (list(k)) then
    nout = nout + 1
    levels_list(nout) = k
  end if
end do

!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED( nout, levels_list, rows, row_length, diagout, diagin)            &
!$OMP private ( i, j, k, kout, kpart, jpart )
!DIR$ IVDEP
do kout = 1, nout
  k     = levels_list(kout)
  kpart = (kout-1) * row_length * rows

  do j = 1, rows
    jpart = (j-1) * row_length
    do i = 1, row_length
      diagout(i+jpart+kpart) = transfer(diagin(i,j,k),real_mould)
    end do ! i
  end do ! j
end do ! k
!$OMP end PARALLEL do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine copydiag_3d_logical

!===============================================================================

subroutine copydiag_3d_integer(diagout,diagin,row_length,rows,levels,          &
     stlist,len_stlist,stash_levels,len_stashlevels)

use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook

use set_levels_list_mod, only: set_levels_list

implicit none

integer, intent(in) ::                                                         &
  row_length,                       & ! Number of points in a row
  rows,                             & ! Number of rows
  levels,                           & ! Number of levels in input data
  len_stlist,                                                                  &
  stlist(len_stlist),               & ! Stash list
  len_stashlevels,                                                             &
  stash_levels(len_stashlevels,levels) ! Stash levels list.

integer,intent(in) :: diagin(row_length,rows,levels)

real,intent(out) :: diagout(row_length * rows * levels)

! Local variables
integer ::                                                                     &
  i,j,k,kout,nout,                                                             &
  kpart,jpart,                                                                 &
  levels_list(levels)         ! list of levels to be copied

logical ::                                                                     &
   list(levels)

real, parameter :: real_mould = 0.0_real_64

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='COPYDIAG_3D_INTEGER'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call set_levels_list(levels,len_stlist,stlist,list,stash_levels,               &
      len_stashlevels)

!  Move data from DIAGIN to DIAGOUT at levels requested

nout = 0
do k = 1, levels
  if (list(k)) then
    nout = nout + 1
    levels_list(nout) = k
  end if
end do

!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC)                               &
!$OMP SHARED( nout, levels_list, rows, row_length, diagout, diagin)            &
!$OMP private ( i, j, k, kout, kpart, jpart )
!DIR$ IVDEP
do kout = 1, nout
  k     = levels_list(kout)
  kpart = (kout-1) * row_length * rows

  do j = 1, rows
    jpart = (j-1) * row_length
    do i = 1, row_length
      diagout(i+jpart+kpart) = transfer(diagin(i,j,k),real_mould)
    end do ! i
  end do ! j
end do ! k
!$OMP end PARALLEL do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine copydiag_3d_integer

end module copydiag_3d_mod
