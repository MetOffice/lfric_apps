! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Purpose: Module to hold all EasyAerosol variables in easyaerosol
!          namelist
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: radiation_control
!
! Code description:
!   Language: FORTRAN 95
!   This code is written to UMDP3 programming standards.
!
! ---------------------------------------------------------------------
module easyaerosol_option_mod

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim
use filenamelength_mod, only: filenamelength
use errormessagelength_mod, only: errormessagelength


implicit none

! Declaration for EasyAerosol
!----------------------------
! Namelist items

! Use EasyAerosol in the SW and/or LW spectrum
logical :: l_easyaerosol_sw = .false.
logical :: l_easyaerosol_lw = .false.

! Use EasyAerosol Cloud Droplet Number Concentrations
! in the radiation scheme (for cloud albedo)
logical :: l_easyaerosol_cdnc = .false.

! Use EasyAerosol Cloud Droplet Number Concentrations
! in the calculation of autoconversion rates.
logical :: l_easyaerosol_autoconv = .false.

! EasyAerosol climatology can be read in as zonal means
logical :: l_easyaerosol_zonal     = .false.

! Number of distributions in each spectrum: extinction, absorption,
! and asymmetry
integer, parameter :: n_easy = 3

! Information on netCDF files for EasyAerosol distributions
integer, parameter :: n_easyaerosol_files = n_easy*2+1 ! = 7 netCDF files
character (len=filenamelength) :: easyaerosol_dir = 'unset'
                                                    ! Directory
character (len=filenamelength) :: easyaerosol_files(n_easyaerosol_files) = [   &
  'unset', 'unset',                                                            &
  'unset', 'unset',                                                            &
  'unset', 'unset',                                                            &
  'unset' ]
                                                    ! Name of EasyAerosol files

! Define the easyaerosol namelist

namelist /easyaerosol/ l_easyaerosol_sw, l_easyaerosol_lw,                     &
                       l_easyaerosol_cdnc, l_easyaerosol_autoconv,             &
                       l_easyaerosol_zonal,                                    &
                       easyaerosol_dir, easyaerosol_files

!DrHook-related parameters
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

character(len=*), parameter, private :: ModuleName='EASYAEROSOL_OPTION_MOD'

contains

subroutine print_nlist_easyaerosol
use umPrintMgr, only: umPrint
implicit none
character(len=50000) :: lineBuffer
integer :: i ! loop counter
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='PRINT_NLIST_EASYAEROSOL'
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call umPrint('Contents of namelist easyaerosol',                               &
    src='easyaerosol_option_mod')

write(lineBuffer,'(A,L1)')' l_easyaerosol_sw = ',l_easyaerosol_sw
call umPrint(lineBuffer,src='easyaerosol_option_mod')
write(lineBuffer,'(A,L1)')' l_easyaerosol_lw = ',l_easyaerosol_lw
call umPrint(lineBuffer,src='easyaerosol_option_mod')
write(lineBuffer,'(A,L1)')' l_easyaerosol_cdnc = ',l_easyaerosol_cdnc
call umPrint(lineBuffer,src='easyaerosol_option_mod')
write(lineBuffer,'(A,L1)')' l_easyaerosol_autoconv = ',l_easyaerosol_autoconv
call umPrint(lineBuffer,src='easyaerosol_option_mod')
write(lineBuffer,'(A,L1)')' l_easyaerosol_zonal = ',l_easyaerosol_zonal
call umPrint(lineBuffer,src='easyaerosol_option_mod')
write(lineBuffer,'(A,A)')' easyaerosol_dir = ',trim(easyaerosol_dir)
call umPrint(lineBuffer,src='easyaerosol_option_mod')
write(lineBuffer,'(A,I0,A,A)')' easyaerosol_files(',1,') = ',                  &
  easyaerosol_files(1)
call umPrint(lineBuffer,src='easyaerosol_option_mod')
do i = 2, n_easyaerosol_files
  if (easyaerosol_files(i) /= 'unset') then
    write(lineBuffer,'(A,I0,A,A)')' easyaerosol_files(',i,') = ',              &
      trim(easyaerosol_files(i))
    call umPrint(lineBuffer,src='easyaerosol_option_mod')
  end if
end do

call umPrint('- - - - - - end of namelist - - - - - -',                        &
             src='easyaerosol_option_mod')

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine print_nlist_easyaerosol


end module easyaerosol_option_mod
