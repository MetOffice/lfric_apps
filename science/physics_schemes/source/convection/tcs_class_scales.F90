! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module defining the tcs warm rain "scales_conv) derived type and
! subroutines for allocating and deallocating an instance of this class.
!
module tcs_class_scales


use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use ereport_mod, only: ereport
use um_types, only: real_umphys

implicit none
!
! Description:
! This module defines the tcs warm rain "scales_conv" derived type and
! subroutines for allocating and deallocating an instance of this class.
!
!
! Method:
!   Variables of type "scales_conv" store arrays of pointers to
! derived scales (e.g. velocity scales, cloud depth, lcl...)
! Each array is defined over convective points.
! Note that assignment of these variables is dealt with in the
! tcs_calc_scales
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.2 programming standards.

type, public :: scales_conv
  integer, pointer :: n_xx
  !
  ! These will have shape (n_xx)
  !
  real(kind=real_umphys), pointer :: wstar_up(:)
                     ! cumulus layer convective velocity scale(m/s)
  real(kind=real_umphys), pointer :: wstar_dn(:)
                     ! subcloud layer velocity scale(m/s)
  real(kind=real_umphys), pointer :: mb(:)
                     ! Cloud base mass flux (m/s)
  real(kind=real_umphys), pointer :: mb_new(:)
                     ! revised mb for incloud calculations (m/s)
  real(kind=real_umphys), pointer :: zcld(:)
                     ! Depth of cloud layer (m)
  real(kind=real_umphys), pointer :: zcld_uv(:)
                     ! Depth of cloud layer (m) CMT cal
  real(kind=real_umphys), pointer :: mb_o_wsc(:)
                     ! mb/wstar_up
  real(kind=real_umphys), pointer :: root_mb_o_wsc(:)
                     ! sqrt of above
  real(kind=real_umphys), pointer :: mb_new_o_wsc(:)
                     ! mb_new/wstar_up
  real(kind=real_umphys), pointer :: root_mb_new_o_wsc(:)
                     ! sqrt of above
  real(kind=real_umphys), pointer :: wstar_up3(:)
                     ! wstar_up**3 * root_mb_new_o_wsc
  real(kind=real_umphys), pointer :: wsc_o_mb(:)
                     ! Convective velocity scale /mb
  real(kind=real_umphys), pointer :: zlcl(:)
                     ! height of the lifting condensation level
  real(kind=real_umphys), pointer :: wup2_cb(:)
                     ! Cloud base updraught velocity squared
                     ! = wup_a1 * wstar_dn(i) * wstar_dn(i)
end type scales_conv

character(len=*), parameter, private :: ModuleName='TCS_CLASS_SCALES'

contains

subroutine allocate_scales(var, n_xx)
  !
  ! Allocates memory to a scales_conv instance "var".
  ! Inputs "n_xx" defines the sizes of the arrays
  !

use errormessagelength_mod, only: errormessagelength

implicit none

type(scales_conv) :: var
integer, optional, target :: n_xx

character(len=*), parameter ::  RoutineName = 'ALLOCATE_SCALES'
character(len=errormessagelength) :: Message
integer :: ErrorStatus ! Return code:
                       !   0 = Normal exit
                       ! +ve = Fatal Error
                       ! -ve = Warning

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
if (present(n_xx)) var%n_xx => n_xx
if (.not. associated(var%n_xx)) then
  ErrorStatus=1
  write(Message, '(A50)')                                                      &
     ' Error allocating scales: Dimensions not defined '

  call Ereport(RoutineName, ErrorStatus, Message)
end if

allocate(var%wstar_up(n_xx))
allocate(var%wstar_dn(n_xx))
allocate(var%mb(n_xx))
allocate(var%mb_new(n_xx))
allocate(var%zcld(n_xx))
allocate(var%zcld_uv(n_xx))
allocate(var%mb_o_wsc(n_xx))
allocate(var%root_mb_o_wsc(n_xx))
allocate(var%mb_new_o_wsc(n_xx))
allocate(var%root_mb_new_o_wsc(n_xx))
allocate(var%wstar_up3(n_xx))
allocate(var%wsc_o_mb(n_xx))
allocate(var%zlcl(n_xx))
allocate(var%wup2_cb(n_xx))
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine allocate_scales

subroutine deallocate_scales(var)
  !
  ! Deallocates memory assocciated with scales instance "var".
  !

implicit none
type(scales_conv) :: var

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DEALLOCATE_SCALES'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
nullify(var%n_xx)

deallocate(var%wup2_cb)
deallocate(var%zlcl)
deallocate(var%wsc_o_mb)
deallocate(var%wstar_up3)
deallocate(var%root_mb_new_o_wsc)
deallocate(var%mb_new_o_wsc)
deallocate(var%root_mb_o_wsc)
deallocate(var%mb_o_wsc)
deallocate(var%zcld_uv)
deallocate(var%zcld)
deallocate(var%mb_new)
deallocate(var%mb)
deallocate(var%wstar_dn)
deallocate(var%wstar_up)
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine deallocate_scales

end module tcs_class_scales
