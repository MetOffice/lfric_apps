! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************

! Purpose:
!   Module containing runtime options/data used by the electrification scheme

! Method:
!   Switches and associated data values used by the electrification scheme
!   are defined here and assigned default values. These may be overridden
!   by namelist input.

!   A description of what each switch or number refers to is provided
!   with the namelist

!   Any routine wishing to use these options may do so with the 'use'
!   statement.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Electric

module electric_inputs_mod

use missing_data_mod,       only: imdi, rmdi
use yomhook,                only: lhook, dr_hook
use parkind1,               only: jprb, jpim
use errormessagelength_mod, only: errormessagelength
use umPrintMgr,             only: umPrint, newline
use um_types,               only: real_umphys

implicit none

!===========================================================================
! Logical options for the electric namelist (default is .false.)
!===========================================================================

! Use interactive lighting in INFERNO fire model
logical :: l_cg_inferno = .false.

! Use updated parameterisation from Luhar et al. at CSIRO
! in Price & Rind scheme
logical :: l_pr_update = .false.

!===========================================================================
! Integer options for the electric namelist (default is imdi)
!===========================================================================
! Choice of method to run with the electric scheme
integer :: electric_method = imdi

! Parameters for electric_method, used to define which scheme to use
! in flash_rate_mod but not used in namelist
integer, parameter :: no_lightning  = 0
integer, parameter :: em_gwp        = 1
integer, parameter :: em_mccaul     = 2
integer, parameter :: em_price_rind = 3

! Storm definition option
integer :: storm_definition = imdi

! Choices for storm definition: graupel only or graupel/ice
integer, parameter :: graupel_only    = 1
integer, parameter :: graupel_and_ice = 2

!===========================================================================
! Real options for the electric namelist (default is rmdi)
!===========================================================================

! Settings for Mccaul et al (2009) scheme
real(kind=real_umphys) :: k1 = rmdi

real(kind=real_umphys) :: k2 = rmdi

! Settings for Graupel water path based scheme

real(kind=real_umphys) :: g1 = rmdi

real(kind=real_umphys) :: g2 = rmdi

! Total Ice and Graupel Water Path thresholds [kg m-2] used to define
! location of a storm.
real(kind=real_umphys) :: tiwp_thresh = rmdi

real(kind=real_umphys) :: gwp_thresh = rmdi

! Define the RUN_ELECTRIC namelist

namelist/run_electric/                                                         &
electric_method, storm_definition, k1, k2, g1, g2, tiwp_thresh, gwp_thresh,    &
l_cg_inferno, l_pr_update

!DrHook-related parameters
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

character(len=*), parameter, private :: ModuleName='ELECTRIC_INPUTS_MOD'

contains

subroutine print_nlist_run_electric()

implicit none

! Local variables

character(len=50000) :: lineBuffer
real(kind=jprb) :: zhook_handle

character(len=*), parameter :: RoutineName='PRINT_NLIST_RUN_ELECTRIC'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Start of subroutine printing

call umPrint('Contents of namelist run_electric',                              &
              src='electric_inputs_mod')

write(lineBuffer,'(A,I0)')'electric_method =', electric_method
call umPrint(lineBuffer,src='electric_inputs_mod')

write(lineBuffer,'(A,I0)')'storm_definition =', storm_definition
call umPrint(lineBuffer,src='electric_inputs_mod')

write(lineBuffer,'(A,F0.4)')'k1 =', k1
call umPrint(lineBuffer,src='electric_inputs_mod')

write(lineBuffer,'(A,F0.4)')'k2 =', k2
call umPrint(lineBuffer,src='electric_inputs_mod')

write(lineBuffer,'(A,F0.4)')'g1 =', g1
call umPrint(lineBuffer,src='electric_inputs_mod')

write(lineBuffer,'(A,F0.4)')'g2 =', g2
call umPrint(lineBuffer,src='electric_inputs_mod')

write(lineBuffer,'(A,F0.4)')'tiwp_thresh =', tiwp_thresh
call umPrint(lineBuffer,src='electric_inputs_mod')

write(lineBuffer,'(A,F0.4)')'gwp_thresh =', gwp_thresh
call umPrint(lineBuffer,src='electric_inputs_mod')

write(lineBuffer,'(A,L1)')'l_cg_inferno =', l_cg_inferno
call umPrint(lineBuffer,src='electric_inputs_mod')

write(lineBuffer,'(A,L1)')'l_pr_update =', l_pr_update
call umPrint(lineBuffer,src='electric_inputs_mod')

call umPrint('- - - - - - end of namelist - - - - - -',                        &
              src='electric_inputs_mod')

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine print_nlist_run_electric


subroutine check_run_electric

use mphys_inputs_mod,  only: graupel_option, no_graupel
use cv_run_mod,        only: l_param_conv
use chk_opts_mod,      only: chk_var
use ereport_mod,       only: ereport
use jules_vegetation_mod,  only: l_inferno
use science_fixes_mod, only: l_fix_gr_autoc

implicit none

character(len=*), parameter :: RoutineName='CHECK_RUN_ELECTRIC'

character(len=errormessagelength) :: comments
character(len=100) :: ChkStr

real(kind=jprb) :: zhook_handle

integer :: ErrorStatus ! Error status 0 = 0K, < 0 = warning, > 0 = error

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

ErrorStatus = 0 ! No error at start of checking

! Check electric_method is sensible. This must fail if it is not
! a sensible number as everything else in the scheme is dependent
! on it.

comments = 'electric_method should be selected in run_electric'

write(ChkStr,'(4(A,I1),A)') '[',no_lightning,',',em_gwp,',',em_mccaul,','      &
                               ,em_price_rind,']'

call chk_var( electric_method, 'electric_method', ChkStr, cmessage=comments)

! If we get this far, electric_method is OK, so can check other variables
! can now be checked based on the value of electric_method

if ( electric_method == em_gwp .or. electric_method == em_mccaul ) then

  ! Fail here if the user has switched off prognostic graupel
  if (graupel_option == no_graupel) then

    ErrorStatus = 100

    comments = 'Prognostic graupel must be included to run with  '//newline//  &
               'microphysics-based lightning scheme. To continue '//newline//  &
               'please switch this option on in the microphysics '//newline//  &
               '(large-scale precipitation). The variable named  '//newline//  &
               'graupel_option must be set to a non-zero value.  '

    call ereport(RoutineName, ErrorStatus, comments)

  end if ! graupel_option == no_graupel

  ! Check against storm definition
  write(ChkStr,'(2(A,I1),A)') '[',graupel_only,',',graupel_and_ice,']'
  call chk_var( storm_definition, 'storm_definition', ChkStr, cmessage=comments)

  if ( electric_method == em_gwp ) then
    ! Specific checks for the GWP scheme

    ! Need to check g1 and g2
    ! These should mirror exactly what is in the Rose metadata
    call chk_var( g1, 'g1', '[1.0e-9:1.0e-7]' )
    call chk_var( g2, 'g2', '[0.0:0.01]' )

  else if ( electric_method == em_mccaul ) then
    ! Specific checks for the McCaul scheme

    ! Need to check k1 and k2
    ! These should mirror exactly what is in the Rose metadata
    call chk_var( k1, 'k1', '[0.0:50.0]' )
    call chk_var( k2, 'k2', '[0.0:50.0]' )
  end if ! electric_method

  ! Specific checks on the storm thresholds used

  if ( storm_definition == graupel_and_ice ) then
    call chk_var( tiwp_thresh, 'tiwp_thresh', '[>=0.0]' )
  end if

  call chk_var( gwp_thresh,  'gwp_thresh', '[>=0.0]' )

else if ( electric_method == em_price_rind ) then

  ! Fail here if the user has switched off the convection scheme
  if ( .not. l_param_conv) then

    ErrorStatus = 100

    comments = 'Model must have parametrized convection included   '//newline//&
               'to run with the convection-based lightning scheme. '//newline//&
               'To continue, please switch l_param_conv to .true.  '//newline//&
               'in the convection scheme panel.'

    call ereport(RoutineName, ErrorStatus, comments)

  end if ! not l_param_conv

end if ! electric_method

if ( l_cg_inferno .and. .not. l_inferno ) then
  ErrorStatus = 100

  comments = 'l_cg_inferno is .true but l_inferno is .false.       '//newline//&
             ' l_inferno needs to be .true. to run interactive lightning.'

  call ereport(RoutineName,ErrorStatus,comments)

  if (electric_method /= em_price_rind ) then
    ErrorStatus = 100

    comments = 'l_cg_inferno is .true. and electric_method is not  '//newline//&
               'em_price_rind (3).'

    call ereport(RoutineName,ErrorStatus,comments)
  end if
end if  ! ( l_cg_inferno .and. .not. l_inferno )

! Produce a warning due to low amounts of lightning in cases where the temporary
! graupel fix is active, the storm definition is set to graupel only and one of
! the two microphysics-based schemes are in use.
if ( l_fix_gr_autoc .and. storm_definition == graupel_only ) then
  if ( electric_method == em_mccaul .or. electric_method == em_gwp ) then
    ErrorStatus = -10

    comments =                                                                 &
          'Bug fix l_fix_gr_autoc is active, however in namelist  '//newline// &
          'run_electric, the variable storm_definition=1          '//newline// &
          '(based on graupel only). While the model will still run'//newline// &
          'the lightning amounts will be very low. You are advised'//newline// &
          'to set storm_definition=2 to get realistic lightning totals'

    call ereport(RoutineName,ErrorStatus,comments)

  end if
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine check_run_electric

end module electric_inputs_mod
