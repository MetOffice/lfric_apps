! *****************************COPYRIGHT********************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT********************************

!  Data module for setting segments sizes used in physics code
! Description:

! Method:
!   segment sizes used by the physics code sections
!   are defined here and assigned default values. These may
!   be overridden by namelist input.

!   Any routine wishing to use these options may do so with the 'use'
!   statement.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: top_level

! Code Description:
!   Language: FORTRAN 95

module tuning_segments_mod

use missing_data_mod, only: rmdi, imdi
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim
use errormessagelength_mod, only: errormessagelength
use chk_opts_mod, only: chk_var, def_src
use bl_option_mod, only: i_bl_vn
use mphys_inputs_mod, only: l_rain, l_casim
use cv_run_mod, only: i_convection_vn, i_cv_comorph, l_param_conv

implicit none

!=======================================================================
! Segments namelist options - ordered as in meta-data sort-key, with namelist
! option and possible values
!=======================================================================

integer :: bl_segment_size = imdi

integer :: gw_seg_size = 32   ! Size of segments for optimising for cache
                              ! use and OpenMP
integer :: ussp_seg_size = imdi ! Segment size for USSP gravity wave drag
integer :: precip_segment_size = imdi

integer :: conv_gr_segment_size = imdi ! Segment size for GR convection scheme

integer :: a_convect_segments = imdi ! No of batches used in convection
                                     ! Original default 1

integer :: a_convect_seg_size = imdi ! Size of convection segments. Can be
                                     ! specified as an alternative to no. of
                                     ! segments. Original default -99
integer :: a_convect_trac_seg_size = imdi ! Size of convection segments for
                                          ! convect. cycles using tracers.
logical :: l_conv_1_seg_per_thread = .false.  ! If true, ignore the above and
                                              ! set number of segments equal
                                              ! to number of OMP threads.

integer :: a_sw_segments = imdi  ! No of batches used in shortwave code
integer :: a_sw_seg_size = -99   ! Size of sw batches, -ve disabled
integer :: a_lw_segments = imdi  ! No of batches used in longwave code
integer :: a_lw_seg_size = -99   ! Size of lw batches, -ve disabled

integer :: ukca_mode_seg_size = imdi  ! No. of columns per chunk
integer :: ukca_chem_seg_size = imdi  ! Segment size for asad chem solver

integer :: gc_radaer_seg_size = imdi  ! No. of columns per chunk

logical :: l_autotune_segments = .false. ! Automatic segment tuning on/off.
logical :: l_autotune_verbose  = .false. ! Prints additional information from
                                         ! the tuning algorithm.
logical :: l_autotune_dry_run  = .false. ! Performs timings etc., but does not
                                         ! change the segment size.
logical :: l_autotune_sync     = .false. ! Perform a barrier during start
                                         ! calliper.
integer :: autotune_trial_radius = imdi  ! Maximum change in segment size
                                         ! (in either direction).
integer :: autotune_min_seg_size = imdi  ! The smallest possible segment size,
                                         ! unless smaller than granularity.
integer :: autotune_max_seg_size = imdi  ! The maximum possible segment size.
integer :: autotune_granularity  = imdi  ! Segment sizes will be a multiple of
                                         ! this value.
integer :: autotune_num_leaders = imdi   ! Num. entries in the leader board.
integer :: autotune_cooling_steps  = imdi ! Num. steps over which the
                                          ! ficitious temperature falls to
                                          ! the specified fraction of its
                                          ! initial value.
real    :: autotune_cooling_fraction = rmdi ! Fraction of the initial
                                            ! ficitious temperature to be
                                            ! reached after the specified
                                            ! number of steps.
real    :: autotune_init_time_coeff = rmdi  ! Proportional increase in initial
                                            ! elapsed time that is to have an
                                            ! acceptance probability of ~50%.

!=======================================================================
!tuning_segments namelist
!=======================================================================

namelist/tuning_segments/ bl_segment_size,gw_seg_size,ussp_seg_size,           &
    precip_segment_size, a_convect_seg_size,                                   &
    a_convect_trac_seg_size, a_convect_segments, l_conv_1_seg_per_thread,      &
    a_lw_seg_size,a_lw_segments, a_sw_seg_size,a_sw_segments,                  &
    ukca_mode_seg_size, ukca_chem_seg_size, gc_radaer_seg_size,                &
    l_autotune_segments, l_autotune_verbose,                                   &
    l_autotune_sync,                                                           &
    l_autotune_dry_run, autotune_trial_radius, autotune_cooling_fraction,      &
    autotune_cooling_steps, autotune_min_seg_size, autotune_max_seg_size,      &
    autotune_granularity, autotune_init_time_coeff, autotune_num_leaders

!DrHook-related parameters
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

character(len=*), parameter, private :: ModuleName='TUNING_SEGMENTS_MOD'

contains

subroutine check_tuning_segments()

! Description:
!   Subroutine to apply logic checks and set control variables based on the
!   options selected in the tuning_segments namelist.

use ereport_mod,  only: ereport
use umPrintMgr,   only: newline
use errormessagelength_mod, only: errormessagelength

implicit none

! Internal error message
character(len=errormessagelength) :: cmessage
! Error code for ereport
integer :: icode

character (len=*), parameter       :: RoutineName = 'CHECK_TUNING_SEGMENTS'
real(kind=jprb)                    :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

if (i_bl_vn > 1)                                                               &
  call chk_var(bl_segment_size, 'bl_segment_size', '[1:999]' )
if (l_rain .and. .not. l_casim)                                                &
  call chk_var(precip_segment_size, 'precip_segment_size', '[1:999]' )
call chk_var(gw_seg_size, 'gw_seg_size', '[1:999]' )
call chk_var(ussp_seg_size, 'ussp_seg_size', '[1:999]' )

if (i_convection_vn /= 0 .and. l_param_conv) then
  if ( l_conv_1_seg_per_thread ) then
    ! Ignore convection segment size / number of segments as they are
    ! overridden with number of segments = number of OMP threads.
    if ( l_autotune_segments ) then
      cmessage =                                                               &
        "Segment size autotuning is switched on, but this cannot " //newline// &
        "affect the segment size used by convection because the "  //newline// &
        "number of segments is overridden to be equal to the "     //newline// &
        "number of OMP threads (l_conv_1_seg_per_thread=.true.)"
      icode = -10
      call ereport( routinename, icode, cmessage )
    end if
    if ( .not. i_convection_vn==i_cv_comorph ) then
      cmessage =                                                               &
        "Option selected to override segment size with the "       //newline// &
        "number of OMP threads (l_conv_1_seg_per_thread=.true.), " //newline// &
        "But this option has only been implemented in the "        //newline// &
        "comorph convection scheme, and a different convection "   //newline// &
        "scheme is in use."
      icode = 10
      call ereport( routinename, icode, cmessage )
    end if
  else  ! ( l_conv_1_seg_per_thread )
    call chk_var(a_convect_seg_size, 'a_convect_seg_size', '[<0,1:999]' )
    call chk_var(a_convect_trac_seg_size, 'a_convect_trac_seg_size',           &
                                                           '[<0,1:999]' )
    if (a_convect_seg_size < 0 .or. a_convect_trac_seg_size < 0) then
      call chk_var(a_convect_segments, 'a_convect_segments', '[<0,1:999]' )
    end if
  end if  ! ( .not. l_conv_1_seg_per_thread )
end if  ! (i_convection_vn /= 0 .and. l_param_conv)

call chk_var(a_lw_seg_size, 'a_lw_seg_size', '[<0,1:999]' )
if ( a_lw_seg_size < 0)                                                        &
  call chk_var(a_lw_segments, 'a_lw_segments', '[<0,1:999]' )
call chk_var(a_sw_seg_size, 'a_sw_seg_size', '[<0,1:999]' )
if ( a_sw_seg_size < 0)                                                        &
  call chk_var(a_sw_segments, 'a_sw_segments', '[<0,1:999]' )
call chk_var(ukca_mode_seg_size, 'ukca_mode_seg_size', '[<0,1:216]' )
call chk_var(ukca_chem_seg_size, 'ukca_chem_seg_size', '[<0,1:999]' )
call chk_var(gc_radaer_seg_size, 'gc_radaer_seg_size', '[<0,1:216]' )

if ( l_autotune_segments ) then
  call chk_var(autotune_trial_radius,                                          &
                    'autotune_trial_radius'    , '[1:999]')
  call chk_var(autotune_granularity,                                           &
                    'autotune_granularity'     , '[1:999]')
  call chk_var(autotune_min_seg_size,                                          &
                    'autotune_min_seg_size'    , '[1:999]')
  call chk_var(autotune_max_seg_size,                                          &
                    'autotune_max_seg_size'    , '[1:999]')
  call chk_var(autotune_num_leaders,                                           &
                    'autotune_num_leaders'     , '[1:999]')
  call chk_var(autotune_cooling_steps,                                         &
                    'autotune_cooling_steps'   , '[1:999]')
  call chk_var(autotune_cooling_fraction,                                      &
                    'autotune_cooling_fraction', '[0.0:1.0]')
  call chk_var(autotune_init_time_coeff,                                       &
                    'autotune_init_time_coeff' , '[1.0:2.0]')
end if

def_src = ''
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine check_tuning_segments

subroutine print_nlist_tuning_segments()
use umPrintMgr, only: umPrint
implicit none
character(len=50000) :: lineBuffer
real(kind=jprb)      :: zhook_handle

character(len=*), parameter :: RoutineName='PRINT_NLIST_TUNING_SEGMENTS'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call umPrint('Contents of namelist tuning_segments',                           &
    src='tuning_segments_mod')

!Segment size / number of segments
write(lineBuffer,'(A,I0)') 'bl_segment_size = ',bl_segment_size
call umPrint(lineBuffer,src='tuning_segments_mod')
write(lineBuffer,'(A,I0)') 'gw_seg_size = ',gw_seg_size
call umPrint(lineBuffer,src='tuning_segments_mod')
write(lineBuffer,'(A,I0)') 'ussp_seg_size = ',ussp_seg_size
call umPrint(lineBuffer,src='tuning_segments_mod')
write(lineBuffer,'(A,I0)') 'precip_segment_size = ',precip_segment_size
call umPrint(lineBuffer,src='tuning_segments_mod')
write(lineBuffer,'(A,I0)') 'a_convect_seg_size = ',a_convect_seg_size
call umPrint(lineBuffer,src='tuning_segments_mod')
write(lineBuffer,'(A,I0)') 'a_convect_trac_seg_size = ',a_convect_trac_seg_size
call umPrint(lineBuffer,src='tuning_segments_mod')
write(lineBuffer,'(A,I0)') 'a_convect_segments = ',a_convect_segments
call umPrint(lineBuffer,src='tuning_segments_mod')
write(lineBuffer,'(A,L1)') 'l_conv_1_seg_per_thread = ',l_conv_1_seg_per_thread
call umPrint(lineBuffer,src='tuning_segments_mod')
write(lineBuffer,'(A,I0)') 'a_lw_seg_size = ',a_lw_seg_size
call umPrint(lineBuffer,src='tuning_segments_mod')
write(lineBuffer,'(A,I0)') 'a_lw_segments = ',a_lw_segments
call umPrint(lineBuffer,src='tuning_segments_mod')
write(lineBuffer,'(A,I0)') 'a_sw_seg_size = ',a_sw_seg_size
call umPrint(lineBuffer,src='tuning_segments_mod')
write(lineBuffer,'(A,I0)') 'a_sw_segments = ',a_sw_segments
call umPrint(lineBuffer,src='tuning_segments_mod')
write(lineBuffer,'(A,I0)') 'ukca_mode_seg_size = ',ukca_mode_seg_size
call umPrint(lineBuffer,src='tuning_segments_mod')
write(lineBuffer,'(A,I0)') 'ukca_chem_seg_size = ',ukca_chem_seg_size
call umPrint(lineBuffer,src='tuning_segments_mod')
write(lineBuffer,'(A,I0)') 'gc_radaer_seg_size = ',gc_radaer_seg_size
call umPrint(lineBuffer,src='tuning_segments_mod')


!Automatic segment size tuning
write(lineBuffer,'(A,L1)') 'l_autotune_segments = ',l_autotune_segments
call umPrint(lineBuffer,src='tuning_segments_mod')

if (l_autotune_segments) then
  write(lineBuffer,'(A,L1)') 'l_autotune_verbose = ',l_autotune_verbose
  call umPrint(lineBuffer,src='tuning_segments_mod')
  write(lineBuffer,'(A,L1)') 'l_autotune_dry_run = ',l_autotune_dry_run
  call umPrint(lineBuffer,src='tuning_segments_mod')
  write(lineBuffer,'(A,L1)') 'l_autotune_sync = ',l_autotune_sync
  call umPrint(lineBuffer,src='tuning_segments_mod')
  write(lineBuffer,'(A,I0)') 'autotune_trial_radius = ',autotune_trial_radius
  call umPrint(lineBuffer,src='tuning_segments_mod')
  write(lineBuffer,'(A,I0)') 'autotune_min_seg_size = ',autotune_min_seg_size
  call umPrint(lineBuffer,src='tuning_segments_mod')
  write(lineBuffer,'(A,I0)') 'autotune_max_seg_size = ',autotune_max_seg_size
  call umPrint(lineBuffer,src='tuning_segments_mod')
  write(lineBuffer,'(A,I0)') 'autotune_granularity = ',autotune_granularity
  call umPrint(lineBuffer,src='tuning_segments_mod')
  write(lineBuffer,'(A,I0)') 'autotune_num_leaders = ',autotune_num_leaders
  call umPrint(lineBuffer,src='tuning_segments_mod')
  write(lineBuffer,'(A,I0)') 'autotune_cooling_steps = ',autotune_cooling_steps
  call umPrint(lineBuffer,src='tuning_segments_mod')
  write(lineBuffer,'(A,ES14.5)') 'autotune_cooling_fraction = ',               &
        autotune_cooling_fraction
  call umPrint(lineBuffer,src='tuning_segments_mod')
  write(lineBuffer,'(A,ES14.5)') 'autotune_init_time_coeff = ',                &
        autotune_init_time_coeff
  call umPrint(lineBuffer,src='tuning_segments_mod')
end if

call umPrint('- - - - - - end of namelist - - - - - -',                        &
    src='tuning_segments_mod')

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine print_nlist_tuning_segments


end module tuning_segments_mod
