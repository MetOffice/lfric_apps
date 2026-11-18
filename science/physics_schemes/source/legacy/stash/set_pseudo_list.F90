! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine SET_PSEUDO_LIST -----------------------------------------
!
! Purpose : To set up a list of pseudo levels at which a diagnostic
!           is required, using information in the STASH list.
!
! Copy of Subroutine SET_LEVELS_LIST (Deck SETLST1) taken and
! adapted for pseudo levels.
!
! Programming Standard : Unified Model Documentation paper number 3
!
! Documentation: U.M. Documentation paper number C4
!
! -----------------------------------------------------------------
!
!   Arguments
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: STASH
module set_pseudo_list_mod

implicit none

character(len=*), parameter, private :: ModuleName='SET_PSEUDO_LIST_MOD'

contains
subroutine set_pseudo_list                                                     &
      (n_levels,len_stlist,stlist,pseudo_list,                                 &
      stash_pseudo_levels,num_stash_pseudo)

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use stparam_mod, only: st_pseudo_in

use ereport_mod, only: ereport
use errormessagelength_mod, only: errormessagelength

implicit none

integer, intent(in) ::                                                         &
n_levels,                                                                      &
  ! Number of possible pseudo levels
len_stlist,                                                                    &
  ! Dimension of STLIST
stlist(len_stlist),                                                            &
  ! STASH list
num_stash_pseudo,                                                              &
  ! Dimension for STASH_PSEUDO_LEVELS
stash_pseudo_levels(num_stash_pseudo+1,*)
  ! Pseudo levels

logical, intent(out) ::                                                        &
  pseudo_list(n_levels) ! ist of pseudo levels required.

!  ---------------------------------------------------------------------

!  Local variables

integer ::                                                                     &
  jlev,                                                                        &
  ! Loop counter over levels
  level_no,                                                                    &
  ! Level no in pseudo list
  list_no
  ! Pseudo level list number

integer :: errorcode
character(len=errormessagelength) :: cmessage

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='SET_PSEUDO_LIST'

!  ---------------------------------------------------------------------

!  Initialise pseudo levels list to false

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
do jlev=1,n_levels
  pseudo_list(jlev)= .false.
end do

!  Get pseudo list number
list_no = stlist(st_pseudo_in)

!  Check that Pseudo list number is valid (should be ge 0)

if (list_no <  0) then
  errorcode=1
  write(cmessage,'(A,4I0)') 'Bad control data: '//                             &
    'st_pseudo_in/stlist(st_pseudo_in)/section/item ',                         &
    st_pseudo_in,stlist(st_pseudo_in),stlist(2), stlist(1)
  call ereport(RoutineName, errorcode, cmessage)
end if

!  Set logical array list to identify pseudo levels required.
if (list_no >  0) then

  do jlev=2,stash_pseudo_levels(1,list_no)+1
    level_no = stash_pseudo_levels(jlev,list_no)
    if (level_no >= 1 .and. level_no <= n_levels) then

      ! Level is within range
      pseudo_list(level_no) =.true.

    else
      ! Level is out of range
      errorcode=2
      write(cmessage,'(A,4I0)') 'Level out of range: ' //                      &
        'list_no/section/item/level ',                                         &
        list_no, stlist(2), stlist(1), level_no
      call ereport(RoutineName, errorcode, cmessage)
    end if
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine set_pseudo_list
end module set_pseudo_list_mod