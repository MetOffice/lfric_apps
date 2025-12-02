! *****************************COPYRIGHT**************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT**************************************
!
! Description:
!   Module containing runtime options/data used by the murk scheme
!
! Method:
!   Switches and associated data values used by the murk scheme
!   are defined here and assigned default values. These may be overridden
!   by namelist input.
!
!   Any routine wishing to use these options may do so with the 'use'
!   statement.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: atmos_service_visibility

module murk_inputs_mod

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim
use missing_data_mod, only: rmdi
use errormessagelength_mod, only: errormessagelength

use um_types, only: real_umphys

implicit none

!===========================================================================
! logical options set from run_murk namelist
!===========================================================================

logical :: l_murk        = .false.  ! Use total aerosol field

logical :: l_murk_advect = .false.  ! Aerosol advection

logical :: l_murk_source = .false.  ! Aerosol source & sink terms

logical :: l_murk_lbc    = .false.  ! Murk aerosol lbcs active

logical :: l_murk_vis    = .false.  ! Murk aerosol used for visibility

!===========================================================================
! real options set from run_murk namelist
!===========================================================================

real(kind=real_umphys) :: murk_source_scale = rmdi
                                   ! Aerosol emissions scaling factor
real(kind=real_umphys) :: m0_murk_scale = rmdi
                                   ! Scaling factor for murk as used in CDNC
!----------------------------------------------------------------------
!----------------------------------------------------------------------
!----------------------------------------------------------------------

!murk washout variables for casim cloud microphysics
real(kind=real_umphys) :: krain_murk = rmdi
real(kind=real_umphys) :: ksnow_murk = rmdi

! Define the run_murk namelist

namelist /run_murk/ l_murk, l_murk_advect, l_murk_source, l_murk_lbc,          &
                    l_murk_vis, murk_source_scale, m0_murk_scale,              &
                    krain_murk, ksnow_murk

!===========================================================================
! logical options not set in namelist
!===========================================================================

logical :: l_murk_bdry   = .false.  ! UK Mes boundary model

logical :: l_murk_rad    = .false.  ! Include radiative effects of aerosol

!DrHook-related parameters
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

character(len=*), parameter, private :: ModuleName='MURK_INPUTS_MOD'

contains

subroutine print_nlist_run_murk()
use umPrintMgr, only: umPrint
implicit none
character(len=50000) :: lineBuffer
real(kind=jprb) :: zhook_handle

character(len=*), parameter :: RoutineName='PRINT_NLIST_RUN_MURK'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call umPrint('Contents of namelist run_murk',                                  &
    src='murk_inputs_mod')

write(lineBuffer,'(A,L1)')' l_murk = ', l_murk
call umPrint(lineBuffer,src='murk_inputs_mod')
write(lineBuffer,'(A,L1)')' l_murk_advect  = ', l_murk_advect
call umPrint(lineBuffer,src='murk_inputs_mod')
write(lineBuffer,'(A,L1)')' l_murk_source = ', l_murk_source
call umPrint(lineBuffer,src='murk_inputs_mod')
write(lineBuffer,'(A,F16.4)')' murk_source_scale = ',murk_source_scale
call umPrint(lineBuffer,src='murk_inputs_mod')
write(lineBuffer,'(A,F16.4)')' m0_murk_scale = ',m0_murk_scale
call umPrint(lineBuffer,src='mphys_inputs_mod')
write(lineBuffer,'(A,L1)')' l_murk_lbc = ', l_murk_lbc
call umPrint(lineBuffer,src='murk_inputs_mod')
write(lineBuffer,'(A,L1)')' l_murk_vis = ', l_murk_vis
call umPrint(lineBuffer,src='murk_inputs_mod')
write(lineBuffer,'(A,F20.6)')' krain_murk = ', krain_murk
call umPrint(lineBuffer,src='murk_inputs_mod')
write(lineBuffer,'(A,F20.6)')' ksnow_murk = ', ksnow_murk
call umPrint(lineBuffer,src='murk_inputs_mod')

call umPrint('- - - - - - end of namelist - - - - - -',                        &
    src='murk_inputs_mod')

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine print_nlist_run_murk


end module murk_inputs_mod
