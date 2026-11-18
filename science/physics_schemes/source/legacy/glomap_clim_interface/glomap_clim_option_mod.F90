! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Module to hold GLOMAP-mode climatology variables in
!          RUN_GLOMAP_AEROCLIM namelist
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: GLOMAP_CLIM
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
module glomap_clim_option_mod

use errormessagelength_mod, only:                                              &
    errormessagelength

use filenamelength_mod,     only:                                              &
    filenamelength

use missing_data_mod,       only:                                              &
    imdi

use parkind1,               only:                                              &
    jpim,                                                                      &
    jprb

use yomhook,                only:                                              &
    dr_hook,                                                                   &
    lhook

implicit none

integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

! Declarations for GLOMAP_CLIM sub-model
! -----------------------------------------------------------------------------

! Namelist items:

! Namelist logicals:
logical :: l_glomap_mode_clim    = .false. ! True when GLOMAP climatology is on

logical :: l_glomap_clim_aie1    = .false. ! True when first aerosol indirect
                                           !  effect required (on radiation)

logical :: l_glomap_clim_aie2    = .false. ! True when second aerosol indirect
                                           !  effect required (on precip.)

logical :: l_glomap_clim_radaer  = .false. ! True when GLOMAP climatology is
                                           ! used in radaer


logical :: l_glomap_clim_radaer_sustrat=.false. ! Use H2SO4 for stratospheric
                                                ! sulphate

! Namelist integers:
integer :: i_glomap_clim_setup = imdi   ! Controls aerosol scheme
                                        ! e.g. SUSSOCBC_5MODE when ==2

                                        ! note SUSSBCOCDU_7mode (8) not
                                        ! available from UM

integer :: i_glomap_clim_nwbins = imdi  ! Controls value of nwbins in Activate
                                        !  See Rosalind West paper for details
                                        !  doi:10.5194/acp-14-6369-2014

integer :: i_glomap_clim_activation_scheme = imdi ! 0 - OFF
                                                  ! 1 - Use Abdul-Razzak and
                                                  !     Ghan Activation method
                                                  !     to calculate CDNC.
                                                  ! 2 - Use Jones method
                                                  !     to calculate CDNC.

! Tuning options for BC absorption
integer :: i_glomap_clim_tune_bc = imdi
! 0 = No tuning. BC density at default value, standard volume-mixing method
!     used for incorporating BC in the refractive index calculation.
! 1 = BC density at tuned value, standard volume-mixing method still used
!     for incorporating BC in the refractive index calculation.
! 2 = BC density at a different tuned value, Maxwell-Garnet method used for
!     incorporating BC in the refractive index calculation.

! Namelist characters:

! RADAER lookup tables and optical properties namelists.
character(len=filenamelength) :: gclmaclw = 'gclmaclw is unset'
                               !  Aitken + Insoluble acc mode (LW)
character(len=filenamelength) :: gclmacsw = 'gclmacsw is unset'
                               !  Aitken + Insoluble acc mode (SW)
character(len=filenamelength) :: gclmanlw = 'gclmanlw is unset'
                               !  Soluble accum mode (LW)
character(len=filenamelength) :: gclmansw = 'gclmansw is unset'
                               !  Soluble accum mode (SW)
character(len=filenamelength) :: gclmcrlw = 'gclmcrlw is unset'
                               !  Coarse mode (LW)
character(len=filenamelength) :: gclmcrsw = 'gclmcrsw is unset'
                               !  Coarse mode (SW)
character(len=filenamelength) :: gclmcnlw = 'gclmcnlw is unset'
                               !  Coarse narrow mode (LW)
character(len=filenamelength) :: gclmcnsw = 'gclmcnsw is unset'
                               !  Coarse narrow mode (SW)
character(len=filenamelength) :: gclmsulw = 'gclmsulw is unset'
                               !  Super-coarse mode (LW)
character(len=filenamelength) :: gclmsusw = 'gclmsusw is unset'
                               !  Super-coarse mode (SW)
character(len=filenamelength) :: gclmprec = 'gclmprec is unset'
                               !  Precomputed values

!Define the RUN_GLOMAP_AEROCLIM namelist
namelist/run_glomap_aeroclim/ l_glomap_mode_clim,                              &
                              l_glomap_clim_aie1,                              &
                              l_glomap_clim_aie2,                              &
                              l_glomap_clim_radaer,                            &
                              l_glomap_clim_radaer_sustrat,                    &
                              i_glomap_clim_setup,                             &
                              i_glomap_clim_nwbins,                            &
                              i_glomap_clim_activation_scheme,                 &
                              i_glomap_clim_tune_bc,                           &
                              gclmaclw, gclmacsw, gclmanlw, gclmansw,          &
                              gclmcrlw, gclmcrsw, gclmcnlw, gclmcnsw,          &
                              gclmsulw, gclmsusw, gclmprec

!===========================================================================
! start of items not set in namelist
!===========================================================================
! Glomap-mode integers:

! case( i_glomap_clim_activation_scheme == 0 )
integer, parameter :: i_gc_activation_off = 0

! case( i_glomap_clim_activation_scheme == 1 )
integer, parameter :: i_gc_activation_arg = 1

! case( i_glomap_clim_activation_scheme == 2 )
integer, parameter :: i_gc_activation_jones = 2

!===========================================================================
! end of items not set in namelist
!===========================================================================

character(len=*), parameter, private :: ModuleName='GLOMAP_CLIM_OPTION_MOD'

contains

!---------------------------------------------------------------------------
subroutine print_nlist_run_glomap_clim ()

use umPrintMgr, only:                                                          &
    umPrint,                                                                   &
    umMessage

implicit none

! Local variables
real(kind=jprb)             :: zhook_handle
character(len=*), parameter :: RoutineName='PRINT_NLIST_RUN_GLOMAP_CLIM'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

write(umMessage,'(A)')    'Contents of namelist run_glomap_aeroclim'
call umPrint(umMessage, src=ModuleName)

write(umMessage,'(A,L1)') ' l_glomap_mode_clim = ', l_glomap_mode_clim
call umPrint(umMessage, src=ModuleName)

write(umMessage,'(A,L1)') ' l_glomap_clim_aie1 = ', l_glomap_clim_aie1
call umPrint(umMessage, src=ModuleName)

write(umMessage,'(A,L1)') ' l_glomap_clim_aie2 = ', l_glomap_clim_aie2
call umPrint(umMessage, src=ModuleName)

write(umMessage,'(A,L1)') ' l_glomap_clim_radaer = ', l_glomap_clim_radaer
call umPrint(umMessage, src=ModuleName)

write(umMessage,'(A,I6)') ' i_glomap_clim_tune_bc = ', i_glomap_clim_tune_bc
call umPrint(umMessage, src=ModuleName)

write(umMessage,'(A,L1)') ' l_glomap_clim_radaer_sustrat = ',                  &
                            l_glomap_clim_radaer_sustrat
call umPrint(umMessage, src=ModuleName)

write(umMessage,'(A,I0)') ' i_glomap_clim_setup = ', i_glomap_clim_setup
call umPrint(umMessage, src=ModuleName)

write(umMessage,'(A,I0)') ' i_glomap_clim_nwbins = ', i_glomap_clim_nwbins
call umPrint(umMessage, src=ModuleName)

write(umMessage,'(A,I0)') ' i_glomap_clim_activation_scheme = ',               &
                            i_glomap_clim_activation_scheme
call umPrint(umMessage, src=ModuleName)

write(umMessage,"(A,A)")  ' gclmaclw = ', trim(gclmaclw)
call umprint(umMessage,src=ModuleName)

write(umMessage,"(A,A)")  ' gclmacsw = ', trim(gclmacsw)
call umprint(umMessage,src=ModuleName)

write(umMessage,"(A,A)")  ' gclmanlw = ', trim(gclmanlw)
call umprint(umMessage,src=ModuleName)

write(umMessage,"(A,A)")  ' gclmansw = ', trim(gclmansw)
call umprint(umMessage,src=ModuleName)

write(umMessage,"(A,A)")  ' gclmcrlw = ', trim(gclmcrlw)
call umprint(umMessage,src=ModuleName)

write(umMessage,"(A,A)")  ' gclmcrsw = ', trim(gclmcrsw)
call umprint(umMessage,src=ModuleName)

write(umMessage,"(A,A)")  ' gclmcnlw = ', trim(gclmcnlw)
call umprint(umMessage,src=ModuleName)

write(umMessage,"(A,A)")  ' gclmcnsw = ', trim(gclmcnsw)
call umprint(umMessage,src=ModuleName)

write(umMessage,"(A,A)")  ' gclmsulw = ', trim(gclmsulw)
call umprint(umMessage,src=ModuleName)

write(umMessage,"(A,A)")  ' gclmsusw = ', trim(gclmsusw)
call umprint(umMessage,src=ModuleName)

write(umMessage,"(A,A)")  ' gclmprec = ', trim(gclmprec)
call umprint(umMessage,src=ModuleName)

write(umMessage,'(A)')    '- - - - - - end of namelist - - - - - -'
call umPrint(umMessage, src = ModuleName)

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine print_nlist_run_glomap_clim

!---------------------------------------------------------------------------
! Description:
!   Subroutine to apply logic checks based on the
!   options selected in the run_glomap_aeroclim namelist.

subroutine check_glomap_clim_options()

use chk_opts_mod, only:                                                        &
    chk_var,                                                                   &
    def_src

use ukca_config_specification_mod,     only:                                   &
    i_sussbcoc_5mode

use umPrintMgr,   only:                                                        &
    newline

implicit none

! Local variables

character(len=*), parameter       :: RoutineName='CHECK_GLOMAP_CLIM_OPTIONS'
real(kind=jprb)                   :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

def_src = RoutineName

if (l_glomap_mode_clim) then

  ! Check if valid value of i_glomap_clim_activation_scheme.
  ! Currently i_glomap_clim_activation_scheme must equal
  ! i_gc_activation_off (OFF) or i_gc_activation_arg (Abdul-Razzak Ghan) or
  ! i_gc_activation_jones (Jones).
  call chk_var ( i_glomap_clim_activation_scheme,                              &
                'i_glomap_clim_activation_scheme',                             &
            [i_gc_activation_off, i_gc_activation_arg, i_gc_activation_jones], &
               cmessage = 'Invalid value of i_glomap_clim_activation_scheme.'  &
              //newline// 'This must == i_gc_activation_off or'                &
              //newline// 'this must == i_gc_activation_arg or'                &
              //newline// 'this must == i_gc_activation_jones.' )

  ! Check if valid value of i_glomap_clim_setup.
  ! Currently i_glomap_clim_setup must equal i_sussbcoc_5mode.
  call chk_var ( i_glomap_clim_setup, 'i_glomap_clim_setup',                   &
                 [i_sussbcoc_5mode],                                           &
                 cmessage = 'Invalid value of i_glomap_clim_setup.'            &
                //newline// 'This must == i_sussbcoc_5mode.' )

end if

def_src = ''

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine check_glomap_clim_options

!---------------------------------------------------------------------------
! Description:
!   Subroutine to apply logic checks based on the
!   options selected in the run_glomap_aeroclim namelist.

subroutine check_run_glomap_clim()

use ereport_mod,      only:                                                    &
    ereport

use mphys_inputs_mod, only:                                                    &
    l_mcr_arcl,                                                                &
    l_autoconv_murk

use ukca_config_specification_mod,     only:                                   &
    i_ukca_activation_arg,                                                     &
    i_ukca_activation_jones
!!!! Should use 'ukca_api_mod' instead after its indirect dependencies on
!!!! 'glomap_clim_option_mod' have been removed

use ukca_option_mod,  only:                                                    &
    i_ukca_activation_scheme,                                                  &
    l_ukca,                                                                    &
    l_ukca_radaer,                                                             &
    l_ukca_aie1,                                                               &
    l_ukca_aie2

use umPrintMgr,       only:                                                    &
    newline

implicit none

! Local variables
character(len=*), parameter       :: RoutineName='CHECK_RUN_GLOMAP_CLIM'
character(len=errormessagelength) :: cmessage
integer                           :: errcode
real(kind=jprb)                   :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Use of ukca and glomap_mode climatology
! This is not currently tested - it may be desirable to use both in future
if (l_ukca .and. l_glomap_mode_clim) then
  cmessage = 'Cannot set both l_ukca & l_glomap_mode_clim to .true.'
  errcode  = 1
  call ereport(RoutineName,errcode,cmessage)
end if

! Direct aerosol effects
if (l_ukca_radaer .and. l_glomap_clim_radaer) then
  cmessage = 'Cannot supply RADAER with aerosol fields from both online'  //   &
             newline// 'GLOMAP-mode (UKCA) and GLOMAP-mode climatology aerosols'
  errcode  = 1
  call ereport(RoutineName,errcode,cmessage)
end if

! Indirect aerosol effects
if ( (l_glomap_clim_aie1 .or. l_glomap_clim_aie2) .and. .not.                  &
     l_glomap_mode_clim ) then
  cmessage = 'Cannot use AIE without GLOMAP-mode aerosols'
  errcode  = 1
  call ereport(RoutineName,errcode,cmessage)
end if

! Activate
if ( ( ( i_ukca_activation_scheme == i_ukca_activation_arg ) .or.              &
       ( i_ukca_activation_scheme == i_ukca_activation_jones ) ) .and.         &
     ( ( i_glomap_clim_activation_scheme == i_gc_activation_arg ) .or.         &
       ( i_glomap_clim_activation_scheme == i_gc_activation_jones ) ) ) then
  cmessage = 'Cannot feed both ukca GLOMAP-mode and climatology aerosols'
  errcode  = 1
  call ereport(RoutineName,errcode,cmessage)
end if

! aerosol indirect effects without CDNC
if ( ( l_glomap_clim_aie1 .or. l_glomap_clim_aie2 ) .and. .not.                &
     ( ( i_glomap_clim_activation_scheme == i_gc_activation_jones ) .or.       &
       ( i_glomap_clim_activation_scheme == i_gc_activation_arg ) ) ) then
  cmessage = 'CDNC required for aerosol indirect effects.' //newline//         &
             'Change i_glomap_clim_activation_scheme so that it is not off.'
  errcode  = 1
  call ereport(RoutineName,errcode,cmessage)
end if

! Warning if Activate or Jones without aerosol indirect effects
if ( ( ( i_glomap_clim_activation_scheme == i_gc_activation_jones ) .or.       &
       ( i_glomap_clim_activation_scheme == i_gc_activation_arg ) ) .and. .not.&
     ( l_glomap_clim_aie1 .or. l_glomap_clim_aie2 ) ) then
  cmessage =  'Both l_glomap_clim_aie1 and l_glomap_clim_aie2 are false;'      &
  //newline// 'you may be calculating a CDNC field'                            &
  //newline// 'which is not required for model evolution.'
  errcode  = -1
  call ereport(RoutineName,errcode,cmessage)
end if

! Aerosol indirect effect
if (l_ukca_aie1 .and. l_glomap_clim_aie1) then
  cmessage = 'Cannot set both l_ukca_aie1 and l_glomap_clim_aie1 to .true.'
  errcode  = 1
  call ereport(RoutineName,errcode,cmessage)
end if

! Aerosol indirect effect
if (l_ukca_aie2 .and. l_glomap_clim_aie2) then
  cmessage = 'Cannot set both l_ukca_aie2 and l_glomap_clim_aie2 to .true.'
  errcode  = 1
  call ereport(RoutineName,errcode,cmessage)
end if

! Aerosol indirect effect
if (l_autoconv_murk .and. l_glomap_clim_aie2) then
  cmessage = 'Cannot set both l_autoconv_murk and l_glomap_clim_aie2 to .true.'
  errcode  = 1
  call ereport(RoutineName,errcode,cmessage)
end if

! Aerosol indirect effect
if (l_mcr_arcl .and. l_glomap_clim_aie2) then
  cmessage = 'Cannot set both l_mcr_arcl and l_glomap_clim_aie2 to .true.'
  errcode  = 1
  call ereport(RoutineName,errcode,cmessage)
end if

if ( ( i_glomap_clim_activation_scheme == i_gc_activation_arg ) .and.          &
     ( i_glomap_clim_nwbins < 1 ) ) then
  cmessage = 'Cannot set i_glomap_clim_nwbins less than one.'
  errcode  = 1
  call ereport(RoutineName,errcode,cmessage)
end if

if ( ( i_glomap_clim_activation_scheme == i_gc_activation_arg ) .and.          &
     ( i_glomap_clim_nwbins > 20 ) ) then
  cmessage = 'Cannot set i_glomap_clim_nwbins greater than twenty.'
  errcode  = 1
  call ereport(RoutineName,errcode,cmessage)
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
end subroutine check_run_glomap_clim

!---------------------------------------------------------------------------

end module glomap_clim_option_mod
