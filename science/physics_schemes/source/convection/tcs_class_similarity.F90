! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module defining the tcs warm rain "similarity" derived type and
! subroutines for allocating and deallocating an instance of this class.
!
module tcs_class_similarity


use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use ereport_mod, only: ereport
use um_types, only: real_umphys

implicit none
!
! Description:
! This module defines the tcs warm rain "similarity" derived type and
! subroutines for allocating and deallocating an instance of this class.
!
!
! Method:
!   Variables of type "similarity" store arrays of in-cloud similarity
! functions as described in
!   <reference to documentation to go here, once available>
! Each array is defined over convective points and in-cloud levels.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.2 programming standards.


type, public :: similarity
  integer, pointer :: n_xx
  integer, pointer :: nlev
  !
  ! These will have shape (n_xx,nlev)
  !
  !--------------------------------------------------------------
  ! These are similarity functions
  !--------------------------------------------------------------
  real(kind=real_umphys), pointer :: g_func(:,:)
                             ! g function
  real(kind=real_umphys), pointer :: k_func(:,:)
                             ! k function
  real(kind=real_umphys), pointer :: f0_func(:,:)
                             ! f0 function
  real(kind=real_umphys), pointer :: f1_func(:,:)
                             ! f1 function
  real(kind=real_umphys), pointer :: fw_func(:,:)
                             ! fw function
  real(kind=real_umphys), pointer :: ftheta_func(:,:)
                             ! ftheta function
  real(kind=real_umphys), pointer :: cmask(:,:)
                             ! cloud mask
  real(kind=real_umphys), pointer :: pc2_detr(:,:)
                             ! detrainment rate for pc2
  !
  ! Arrays for functions of height - uv levels
  !
  real(kind=real_umphys), pointer :: fql_func_rho(:,:)
                             ! fql function
  real(kind=real_umphys), pointer :: g_func_rho(:,:)
                             ! g function
  real(kind=real_umphys), pointer :: fw_func_rho(:,:)
                             ! fw function
  real(kind=real_umphys), pointer :: fng_func_rho(:,:)
                             ! fng function
  real(kind=real_umphys), pointer :: k_func_rho(:,:)
                             ! k function
  real(kind=real_umphys), pointer :: b_func_rho(:,:)
                             ! b function
  real(kind=real_umphys), pointer :: gql_func_rho(:,:)
                             ! gql function
  real(kind=real_umphys), pointer :: cmask_rho(:,:)
                             ! cloud mask
end type similarity

character(len=*), parameter, private :: ModuleName='TCS_CLASS_SIMILARITY'

contains

subroutine allocate_similarity(var, n_xx, nlev, initval)
!
! Allocates memory to a similarity instance "var".
! Inputs "n_xx" and "nlev" define the sizes of the arrays
! and if present all fields are initialized to the value "initval"
!

use errormessagelength_mod, only: errormessagelength
implicit none

type(similarity), intent(in out) :: var
integer, intent(in), optional, target :: n_xx, nlev
real(kind=real_umphys), intent(in), optional :: initval
logical :: l_init

character(len=*), parameter ::  RoutineName = 'ALLOCATE_SIMILARITY'
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
if (present(nlev)) var%nlev => nlev
if (.not. (associated(var%n_xx)                                                &
           .and. associated(var%nlev))) then
  ErrorStatus=1
  write(Message, '(A53)')                                                      &
     ' Error allocating similarity: Dimensions not defined '

  call Ereport(RoutineName, ErrorStatus, Message)
end if

allocate(var%g_func(var%n_xx,var%nlev))
allocate(var%k_func(var%n_xx,var%nlev))
allocate(var%f0_func(var%n_xx,var%nlev))
allocate(var%f1_func(var%n_xx,var%nlev))
allocate(var%fw_func(var%n_xx,var%nlev))
allocate(var%ftheta_func(var%n_xx,var%nlev))
allocate(var%cmask(var%n_xx,var%nlev))
allocate(var%pc2_detr(var%n_xx,var%nlev))
allocate(var%fql_func_rho(var%n_xx,var%nlev))
allocate(var%g_func_rho(var%n_xx,var%nlev))
allocate(var%fw_func_rho(var%n_xx,var%nlev))
allocate(var%fng_func_rho(var%n_xx,var%nlev))
allocate(var%k_func_rho(var%n_xx,var%nlev))
allocate(var%b_func_rho(var%n_xx,var%nlev))
allocate(var%gql_func_rho(var%n_xx,var%nlev))
allocate(var%cmask_rho(var%n_xx,var%nlev))

! initialise if required.
if (present(initval)) then
  l_init=.true.
else
  l_init=.false.
end if
if (l_init) then
  var%g_func       = real(initval)
  var%k_func       = real(initval)
  var%f0_func      = real(initval)
  var%f1_func      = real(initval)
  var%fw_func      = real(initval)
  var%ftheta_func  = real(initval)
  var%cmask        = real(initval)
  var%pc2_detr     = real(initval)
  var%fql_func_rho = real(initval)
  var%g_func_rho   = real(initval)
  var%fw_func_rho  = real(initval)
  var%fng_func_rho = real(initval)
  var%k_func_rho   = real(initval)
  var%b_func_rho   = real(initval)
  var%gql_func_rho = real(initval)
  var%cmask_rho    = real(initval)
end if
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine allocate_similarity

subroutine deallocate_similarity(var)
  !
  ! Deallocates memory assocciated with similarity instance "var".
  !
implicit none
type(similarity) :: var

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DEALLOCATE_SIMILARITY'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
nullify(var%n_xx)
nullify(var%nlev)

deallocate(var%cmask_rho)
deallocate(var%gql_func_rho)
deallocate(var%b_func_rho)
deallocate(var%k_func_rho)
deallocate(var%fng_func_rho)
deallocate(var%fw_func_rho)
deallocate(var%g_func_rho)
deallocate(var%fql_func_rho)
deallocate(var%pc2_detr)
deallocate(var%cmask)
deallocate(var%ftheta_func)
deallocate(var%fw_func)
deallocate(var%f1_func)
deallocate(var%f0_func)
deallocate(var%k_func)
deallocate(var%g_func)
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine deallocate_similarity

end module tcs_class_similarity
