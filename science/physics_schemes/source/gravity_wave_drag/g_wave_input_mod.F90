! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module g_wave_input_mod

! Description:
!       Input/namelist control of GWD scheme.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: gravity_wave_drag

! Code Description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 v8 programming standards.

use missing_data_mod, only: rmdi,imdi
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim
use errormessagelength_mod, only: errormessagelength
use ereport_mod, only: ereport

use um_types, only: real_umphys

implicit none

! Control logicals
logical :: l_gwd      = .false.  ! Use orographic drag scheme
logical :: l_use_ussp = .false.  ! Use spectral GWD scheme, with opaque lid

logical :: l_gwd_40km = .true.   ! Turn off orographic GWD above 40km

logical :: l_nonhydro = .false.  ! Turn on nonhydro scheme
logical :: l_dynbeta  = .false.  ! Turn on dynamically adjusting angle
                                 ! of group vel in nonhydro scheme

logical :: l_taus_scale = .false.! Froude No. Dependent surface stress

logical :: l_fix_gwsatn = .true. ! Bug fixes to gwsatn
integer :: sat_scheme = 0        ! Switch to determine whether to use
                                 ! stress based    (sat_scheme=0) or
                                 ! amplitude based (sat_scheme=1)
                                 ! saturation test

! GWD variable inputs

integer :: i_gwd_vn = imdi    ! Switch to determine version of GWD scheme
                              ! 4 => 4A scheme
                              ! 5 => 5A scheme
integer, parameter :: i_gwd_vn_4a = 4
integer, parameter :: i_gwd_vn_5a = 5

real(kind=real_umphys) :: kay_gwave   =  rmdi
                              ! Surface stress constant for GW drag
real(kind=real_umphys) :: gwd_frc     =  rmdi
                              ! Critical Froude Number for 4A Scheme
real(kind=real_umphys) :: nsigma      =  rmdi
                              ! Number of standard deviations above the
                              ! mean orography of top of sub-grid mountains

! Maximum and minimum values for the STPH_RP scheme
! Gravity Wave drag parameters moved to stochastic_physics_run_mod.F90
! Max and min values for the critical froud number
! Max and min values for the gravity wave constant

real(kind=real_umphys) :: gwd_fsat  = rmdi
                              ! Critical Froude number for wave breaking
                              ! used in the amplitude based saturation test


!New parameters for 5A version of orographic drag scheme
real(kind=real_umphys)    :: gsharp  = rmdi
                              ! Function of mountain sharpness (used to tune
                              ! amplitude of orographic gwd)
real(kind=real_umphys)    :: fbcd = rmdi        ! Flow blocking drag coefficient

logical :: l_smooth = .false. ! Turn on smoothing of acceleration over
                              ! a vertical wavelength

logical :: scale_aware = .false. ! Use scale_aware scheme
real(kind=real_umphys) :: middle =  rmdi
real(kind=real_umphys) :: var    =  rmdi

logical :: l_gw_heating(3) = .false.
                              ! Turn on for heating due to
                              ! (1) flow blocking drag
                              ! (2) mountain-wave dissipation
                              ! (3) non-orographic gravity-wave dissipation

! Switch to determine method for calculating N
integer :: i_moist      = imdi
                        ! 0 for standard dry buoyancy frequency N
                        ! 1 for moist N at low-levels (blocking and total wave
                        ! stress) but dry N for distributing wave stress in the
                        ! vertical
                        ! 2 for moist buoyancy frequency throughout the
                        ! entire Orographic Drag (OD) scheme

! Non-orographic gravity wave (USSP scheme) variable inputs
!
! Switch to determine version of USSP scheme (presently always '1a')
integer :: i_ussp_vn = imdi
!
! 1 => 1A scheme (Standard global invariant source)
integer, parameter :: i_ussp_vn_1a = 1
!
! Factor enhancement for invariant global wave launch amplitude
real(kind=real_umphys)    :: ussp_launch_factor = rmdi
!
! Factor enhancement for conversion from convective rain to GW source flux
real(kind=real_umphys)    :: cgw_scale_factor   = rmdi
!
! Characteristic (spectrum peak) wavelength (m)
real(kind=real_umphys)    :: wavelstar          = rmdi
!
! Switch for precipitation dependent variable source flux calculation
! F : calculate standard USSP isotropic GW launch flux
! T : calculate variable CGW launch flux
logical :: l_add_cgw = .false.

namelist/run_gwd/i_gwd_vn                                                      &
                ,kay_gwave, gwd_frc, nsigma, gwd_fsat                          &
                , gsharp, fbcd, l_smooth, l_gw_heating                         &
                ,scale_aware, middle, var                                      &
                ,ussp_launch_factor, wavelstar, l_gwd, l_use_ussp,             &
                 i_ussp_vn, l_add_cgw, cgw_scale_factor,                       &
                 i_moist

!DrHook-related parameters
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

character(len=*), parameter, private :: ModuleName='G_WAVE_INPUT_MOD'

contains
subroutine print_nlist_run_gwd()
use umPrintMgr, only: umPrint,ummessage,PrNorm
implicit none
real(kind=jprb) :: zhook_handle

character(len=*), parameter :: RoutineName='PRINT_NLIST_RUN_GWD'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call umPrint('Contents of namelist run_gwd',                                   &
    src='g_wave_input_mod')

write(ummessage,'(A,I0)')' i_gwd_vn = ',i_gwd_vn
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,F14.2)')' kay_gwave = ',kay_gwave
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,F14.2)')' gwd_frc = ',gwd_frc
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,F14.2)')' nsigma = ',nsigma
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,F14.2)')' gwd_fsat = ',gwd_fsat
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,F14.2)')' gsharp = ',gsharp
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,F14.2)')' fbcd = ',fbcd
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,L1)')' l_smooth = ',l_smooth
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,L1)')' scale_aware = ',scale_aware
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,F14.2)')' middle = ',middle
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,F14.2)')' var = ',var
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,3L1)')' l_gw_heating = ',l_gw_heating
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,F16.4)')' ussp_launch_factor = ',ussp_launch_factor
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,F16.4)')' wavelstar = ',wavelstar
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,L1)')' l_gwd = ',l_gwd
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,L1)')' l_use_ussp = ',l_use_ussp
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,I0)')' i_ussp_vn = ',i_ussp_vn
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,L1)')' l_add_cgw = ',l_add_cgw
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,F16.4)')' cgw_scale_factor = ',cgw_scale_factor
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')
write(ummessage,'(A,I0)')' i_moist = ',i_moist
call umPrint(ummessage,level=PrNorm,src='g_wave_input_mod')

call umPrint('- - - - - - end of namelist - - - - - -',                        &
    src='g_wave_input_mod')

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine print_nlist_run_gwd


end module g_wave_input_mod
