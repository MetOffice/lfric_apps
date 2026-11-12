! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SET_LEVELS_LIST
!
! Purpose : To set up a list of levels at which a diagnostic is
!           required, using information in the STASH list.
! Service routine  version for Cray YMP
!
! Programming Standard : Unified Model Documentation paper number 4
!                      : Version no 2, dated 18/01/90
!
! System components covered : D3
!
! System task : P0
!
! Documentation: U.M. Documentation paper number P0,C4
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH

module set_levels_list_mod

implicit none

character(len=*), parameter, private :: ModuleName='SET_LEVELS_LIST_MOD'

contains
subroutine set_levels_list(levels,len_stlist,stlist,list,stash_levels,         &
                           len_stashlevels)

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use ereport_mod, only: ereport
use errormessagelength_mod, only: errormessagelength

implicit none

integer, intent(in) ::                                                         &
  levels,                                                                      &
    ! Number of levels in input data
  len_stlist,                                                                  &
  stlist(len_stlist),                                                          &
    ! STASH list
  len_stashlevels,                                                             &
  stash_levels(len_stashlevels,*)
    ! list of levels required

logical, intent(out) ::                                                        &
  list(levels) ! List of levels required.

!  Local variables
integer :: k,kout

integer :: errorcode
character(len=errormessagelength) :: cmessage

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='SET_LEVELS_LIST'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!  Initialise levels list to false
do k=1,levels
  list(k)= .false.
end do

!  Check for method of levels selection
!  Levels list must be present.

if (stlist(10) <  0) then

  ! Set logical array list to identify levels required.
  do kout=2,stash_levels(1,-stlist(10))+1

    if ((stash_levels(kout,-stlist(10)) >= 1) .and.                            &
        (stash_levels(kout,-stlist(10)) <= levels)) then
      ! LEVEL is in THE range OF LIST.
      list(stash_levels(kout,-stlist(10))) =.true.
    else ! LEVEL is out OF THE range OF LIST.
      errorcode=2
      write(cmessage,'(A,3I0)') 'Level out of range: section/item/level ',     &
        stlist(2), stlist(1), stash_levels(kout,-stlist(10))
      call ereport(RoutineName, errorcode, cmessage)
    end if
  end do

else if (stlist(10) /= 100) then

  !  Set list of levels according to its definition in stlist :
  !l If stlist(10) positive, input on range of model levels starting
  !l at level stlist(10), finishing at level stlist(11)
  do k=stlist(10),stlist(11)
    list(k)= .true.
  end do

else ! Illegal control data
  errorcode=1
  write(cmessage,'(A,4I0)') 'Bad control data: stlist(10/11)/section/item ',   &
    stlist(10) ,stlist(11), stlist(2), stlist(1)
  call ereport(RoutineName, errorcode, cmessage)
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine set_levels_list
end module set_levels_list_mod