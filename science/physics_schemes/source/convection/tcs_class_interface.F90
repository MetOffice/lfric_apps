! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module defining the tcs warm rain "interface_input" derived type and
! subroutines for allocating and deallocating an instance of this class
! and also for assigning values.
!
module tcs_class_interface


use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use ereport_mod, only: ereport
use um_types, only: real_umphys

implicit none
!
! Description:
! This module defines the tcs warm rain "cloud_input" derived type and
! subroutines for allocating and deallocating an instance of this class.
!
!
! Method:
!   Variables of type "interface_input" store arrays of pointers to
! arrays sampled on interface levels (e.g.  cloud base, freezing level,
! inversion base...)
! Each array is defined over convective points.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.2 programming standards.

type, public :: interface_input
   !
   ! These are fields evaluated on interface levels,
   ! e.g. cloudbase (ntml), inversion (ntpar) and levels around these
   !
  integer, pointer :: n_xx
  ! Number of convectively active points
  integer, pointer :: ntra
  ! Number of tracers
  ! These fields will all have shape (n_xx)
  integer, pointer :: levels(:) ! The levels on which the field is evaluated
  real(kind=real_umphys), pointer :: theta(:)
  real(kind=real_umphys), pointer :: q_mix(:)
  real(kind=real_umphys), pointer :: qse(:)
  real(kind=real_umphys), pointer :: p_theta(:)  ! = p_layer_centres
  real(kind=real_umphys), pointer :: p_rho(:)    ! = p_layer_boundaries
  real(kind=real_umphys), pointer :: exner_theta(:)   ! = exner_layer_centres
  real(kind=real_umphys), pointer :: exner_rho(:)     ! = exner_layer_boundaries
  real(kind=real_umphys), pointer :: z_theta(:)
  real(kind=real_umphys), pointer :: z_rho(:)
  real(kind=real_umphys), pointer :: rho(:)
  !These will have shape (n_xx,ntra)
  real(kind=real_umphys), pointer :: tracer(:,:)
end type interface_input

character(len=*), parameter, private :: ModuleName='TCS_CLASS_INTERFACE'

contains

subroutine allocate_interface_input(var, n_xx, ntra, initval)
  !
  ! Allocates memory to a interface_input instance "var".
  ! Inputs "n_xx" and "ntra" define the sizes of the arrays
  ! and if present all fields are initialized to the value "initval"
  !
use errormessagelength_mod, only: errormessagelength

implicit none
type(interface_input) :: var
integer, optional, target :: n_xx, ntra
real(kind=real_umphys), intent(in), optional :: initval
logical :: l_init

character(len=*), parameter ::  RoutineName = 'ALLOCATE_INTERFACE_INPUT'
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
if (present(ntra)) var%ntra => ntra
if (.not. (associated(var%n_xx) .and. associated(var%ntra))) then
  ErrorStatus=1
  write(Message, '(A54)')                                                      &
     ' Error allocating interface_input: Dimensions not defined '

  call Ereport(RoutineName, ErrorStatus, Message)
end if

! variables with shape (n_xx)
allocate(var%levels(var%n_xx))
allocate(var%theta(var%n_xx))
allocate(var%q_mix(var%n_xx))
allocate(var%qse(var%n_xx))
allocate(var%p_theta(var%n_xx))
allocate(var%p_rho(var%n_xx))
allocate(var%exner_theta(var%n_xx))
allocate(var%exner_rho(var%n_xx))
allocate(var%z_theta(var%n_xx))
allocate(var%z_rho(var%n_xx))
allocate(var%rho(var%n_xx))

! variables with shape (n_xx, ntra)
allocate(var%tracer(var%n_xx,var%ntra))

! initialise if required.
if (present(initval)) then
  l_init=.true.
else
  l_init=.false.
end if
if (l_init) then
  var%theta       =real(initval)
  var%q_mix       =real(initval)
  var%qse         =real(initval)
  var%p_theta     =real(initval)
  var%p_rho       =real(initval)
  var%exner_theta =real(initval)
  var%exner_rho   =real(initval)
  var%z_theta     =real(initval)
  var%z_rho       =real(initval)
  var%rho         =real(initval)
  var%tracer      =real(initval)
end if
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine allocate_interface_input

subroutine deallocate_interface_input(var)
  !
  ! Deallocates memory assocciated with cloud_input instance "var".
  !

implicit none
type(interface_input) :: var

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DEALLOCATE_INTERFACE_INPUT'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
nullify(var%n_xx)
nullify(var%ntra)
!    nullify(var%levels)

deallocate(var%tracer)

deallocate(var%rho)
deallocate(var%z_rho)
deallocate(var%z_theta)
deallocate(var%exner_rho)
deallocate(var%exner_theta)
deallocate(var%p_rho)
deallocate(var%p_theta)
deallocate(var%qse)
deallocate(var%q_mix)
deallocate(var%theta)
deallocate(var%levels)
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine deallocate_interface_input

subroutine assign_interface_input(var, levels, theta,                          &
   q_mix, qse, p_theta, p_rho, exner_theta, exner_rho,                         &
   z_theta, z_rho, rho, tracer )
  !
  ! Assigns values to a interface_input instance "var".
  !
implicit none
type(interface_input), intent(in out) :: var
integer, intent(in), target :: levels(:)
real(kind=real_umphys), intent(in) :: theta(:,:)
real(kind=real_umphys), intent(in) :: q_mix(:,:)
real(kind=real_umphys), intent(in) :: qse(:,:)
real(kind=real_umphys), intent(in) :: p_theta(:,0:)      ! input starts 0
real(kind=real_umphys), intent(in) :: p_rho(:,0:)        ! input starts 0
real(kind=real_umphys), intent(in) :: exner_theta(:,0:)  ! input starts 0
real(kind=real_umphys), intent(in) :: exner_rho(:,0:)    ! input starts 0
real(kind=real_umphys), intent(in) :: z_theta(:,:)
real(kind=real_umphys), intent(in) :: z_rho(:,:)
real(kind=real_umphys), intent(in) :: rho(:,:)
real(kind=real_umphys), intent(in) :: tracer(:,:,:)

integer :: i

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='ASSIGN_INTERFACE_INPUT'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
do i=1,var%n_xx
  var%levels(i)      = levels(i)
  var%theta(i)       = theta(i,levels(i))
  var%q_mix(i)       = q_mix(i,levels(i))
  var%qse(i)         = qse(i,levels(i))
  var%p_theta(i)     = p_theta(i,levels(i))
  var%p_rho(i)       = p_rho(i,levels(i))
  var%exner_theta(i) = exner_theta(i,levels(i))
  var%exner_rho(i)   = exner_rho(i,levels(i))
  var%z_theta(i)     = z_theta(i,levels(i))
  var%z_rho(i)       = z_rho(i,levels(i))
  var%rho(i)         = rho(i,levels(i))
  ! Working note: need to put on a check that tracer
  ! levels go up to levels(i) - c.f. also assign_cloud_input
  var%tracer(i,:)    = tracer(i,levels(i),:)
end do
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine assign_interface_input

end module tcs_class_interface
