! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  Subroutine V_TO_P  for calculating variables held at v points
!  at p points.
!
!  This routine does interior points of array not halos,
!  but requires halo information to be set.
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: Grids
module v_to_p_mod

implicit none

character(len=*), parameter, private :: ModuleName = 'V_TO_P_MOD'

contains

subroutine v_to_p                                                              &
  ( array_on_v_points, ini_start, ini_end, inj_start, inj_end                  &
  , outi_start, outi_end, outj_start, outj_end, levels, array_on_p_points )

use model_domain_mod,   only: model_type, mt_single_column, mt_lfric
!$ use compute_chunk_size_mod, only: compute_chunk_size  ! Note OpenMP sentinel

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer, intent(in) :: ini_start, ini_end
integer, intent(in) :: inj_start, inj_end
integer, intent(in) :: outi_start, outi_end
integer, intent(in) :: outj_start, outj_end
integer, intent(in) :: levels

real, intent(in)  :: array_on_v_points( ini_start:ini_end                      &
                                      , inj_start:inj_end                      &
                                      , levels )

real, intent(out) :: array_on_p_points( outi_start:outi_end                    &
                                      , outj_start:outj_end                    &
                                      , levels )

! local variables
integer :: i, j, k
integer :: j0, j1
integer :: ompt_start   ! start of thread's work
integer :: ompt_end     ! end of thread's work

! Parameters
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName = 'V_TO_P'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! this procedure is now thread safe. i.e. it can be called from an
! OMP parallel region. Compute_chunk_size will determine the start and end
! iterations for each thread as if it were a static schedule.
! If called from a serial region, start and end are 1 and levels

ompt_start = 1
ompt_end = levels
! only call compute_chunk_size if compiling with OMP
! The procedure call is protected by the optional compile
! sentinel
!$ call compute_chunk_size(1,levels,ompt_start,ompt_end)

j0 = outj_start
j1 = outj_end

if ((model_type == mt_single_column) .or. (model_type == mt_lfric)) then
  do k=ompt_start, ompt_end
    do j=j0, j1
      do i=outi_start, outi_end
        array_on_p_points(i,j,k) = array_on_v_points(i,j,k)
      end do
    end do
  end do

else !not scm or lfric
  do k=ompt_start, ompt_end
    do j=j0, j1
      do i=outi_start, outi_end
        array_on_p_points(i,j,k) = 0.5 * ( array_on_v_points(i,j,k)            &
                                         + array_on_v_points(i,j-1,k) )
      end do
    end do
  end do

end if !model_type

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine v_to_p

end module v_to_p_mod
