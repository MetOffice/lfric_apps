! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Module defining the tcs warm rain cloud_input derived type and
! subroutines for allocating and deallocating an instance of this class
! and also for assigning values.
!
module tcs_class_cloud


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
!   Variables of type "cloud_input" store arrays of pointers to
! those sections of the existing full model fields, or derivations of
! them, within the cloud layer. Each array is defined over convective
! points and in-cloud levels.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.2 programming standards.

type :: cloud_input
   !--------------------------------------------------------------
   ! These are arrays which are sampled on cloud levels
   !--------------------------------------------------------------
  integer, pointer :: n_xx
                            ! Number of convectively active points
  integer, pointer :: nlev
                            ! Number of levels
  integer, pointer :: ntra
                            ! Number of tracers
  ! These will have shape (n_xx,nlev)
  real(kind=real_umphys), pointer :: eta_theta(:,:)
                            ! Non-dimensional height (theta levels)
  real(kind=real_umphys), pointer :: eta_rho(:,:)
                            ! Non-dimensional height (rho levels)
  real(kind=real_umphys), pointer :: theta(:,:)
                            ! theta
  real(kind=real_umphys), pointer :: q_mix(:,:)
                            ! q - mixing ratio
  real(kind=real_umphys), pointer :: qse(:,:)
                            ! saturation mixing ratio
  real(kind=real_umphys), pointer :: exner_theta(:,:)
                            ! exner on theta levels in cloud
  real(kind=real_umphys), pointer :: exner_rho(:,:)
                            ! exner pressure on rho in cloud
  real(kind=real_umphys), pointer :: p_theta(:,:)
                            ! pressure on theta levels in cloud
  real(kind=real_umphys), pointer :: p_rho(:,:)
                            ! pressure on rho levels in cloud
  real(kind=real_umphys), pointer :: z_theta(:,:)
                            ! height of theta levels
  real(kind=real_umphys), pointer :: z_rho(:,:)
                            ! height of rho levels
  real(kind=real_umphys), pointer :: r2rho(:,:)
                            ! r2*rho rho levels (kg/m)
  real(kind=real_umphys), pointer :: r2rho_theta(:,:)
                            ! r2*rho theta levels (kg/m)
  real(kind=real_umphys), pointer :: rho(:,:)
                            ! density on  rho levels
  ! These will have shape (n_xx,nlev,ntra)
  real(kind=real_umphys), pointer :: tracer(:,:,:)
                            ! tracer in cloud
end type cloud_input

character(len=*), parameter, private :: ModuleName='TCS_CLASS_CLOUD'

contains

subroutine allocate_cloud_input(var, n_xx, nlev, ntra, initval)
  !
  ! Allocates memory to a cloud_input instance "var".
  ! Inputs "n_xx", "nlev" and "ntra" define the sizes of the arrays
  ! and if present all fields are initialized to the value "initval"
  !

use errormessagelength_mod, only: errormessagelength
implicit none

type(cloud_input) :: var
integer, optional, target :: n_xx, nlev, ntra
real(kind=real_umphys), intent(in), optional :: initval
logical :: l_init

character(len=*), parameter ::  RoutineName = 'ALLOCATE_CLOUD_INPUT'
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
if (present(ntra)) var%ntra => ntra
if (.not. (associated(var%n_xx)                                                &
           .and. associated(var%nlev)                                          &
           .and. associated(var%ntra))) then
  ErrorStatus=1
  write(Message, '(A54)')                                                      &
     ' Error allocating cloud_input: Dimensions not defined '

  call Ereport(RoutineName, ErrorStatus, Message)
end if

allocate(var%eta_theta(var%n_xx,var%nlev))
allocate(var%eta_rho(var%n_xx,var%nlev))
allocate(var%theta(var%n_xx,var%nlev))
allocate(var%q_mix(var%n_xx,var%nlev))
allocate(var%qse(var%n_xx,var%nlev))
allocate(var%exner_theta(var%n_xx,var%nlev))
allocate(var%exner_rho(var%n_xx,var%nlev))
allocate(var%p_theta(var%n_xx,var%nlev))
allocate(var%p_rho(var%n_xx,var%nlev))
allocate(var%z_theta(var%n_xx,var%nlev))
allocate(var%z_rho(var%n_xx,var%nlev))
allocate(var%r2rho(var%n_xx,var%nlev))
allocate(var%r2rho_theta(var%n_xx,var%nlev))
allocate(var%rho(var%n_xx,var%nlev))
allocate(var%tracer(var%n_xx,var%nlev,var%ntra))

! initialise if required.
if (present(initval)) then
  l_init=.true.
else
  l_init=.false.
end if
if (l_init) then
  var%eta_theta   =real(initval)
  var%eta_rho     =real(initval)
  var%theta       =real(initval)
  var%q_mix       =real(initval)
  var%qse         =real(initval)
  var%exner_theta =real(initval)
  var%exner_rho   =real(initval)
  var%p_theta     =real(initval)
  var%p_rho       =real(initval)
  var%z_theta     =real(initval)
  var%z_rho       =real(initval)
  var%r2rho       =real(initval)
  var%r2rho_theta =real(initval)
  var%rho         =real(initval)
  var%tracer      =real(initval)
end if
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine allocate_cloud_input

subroutine deallocate_cloud_input(var)
  !
  ! Deallocates memory assocciated with cloud_input instance "var".
  !

implicit none
type(cloud_input) :: var

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='DEALLOCATE_CLOUD_INPUT'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
nullify(var%n_xx)
nullify(var%nlev)
nullify(var%ntra)

deallocate(var%tracer)
deallocate(var%rho)
deallocate(var%r2rho_theta)
deallocate(var%r2rho)
deallocate(var%z_rho)
deallocate(var%z_theta)
deallocate(var%p_rho)
deallocate(var%p_theta)
deallocate(var%exner_rho)
deallocate(var%exner_theta)
deallocate(var%qse)
deallocate(var%q_mix)
deallocate(var%theta)
deallocate(var%eta_theta)
deallocate(var%eta_rho)
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine deallocate_cloud_input

subroutine assign_cloud_input(                                                 &
     var, base_levels,  maxlevs , maxtrlevs,                                   &
     theta, q_mix, qse, p_theta, p_rho, exner_theta, exner_rho,                &
     z_theta, z_rho, r2rho, r2rho_theta, rho, tracer )
  !
  ! Assigns values to a cloud_input instance "var".
  !

use tcs_common_warm , only: scales, cb
use errormessagelength_mod, only: errormessagelength
implicit none
type(cloud_input), intent(in out) :: var
integer, intent(in) :: base_levels(:)  ! Levels at base of region
! Working note: maxlevs would normally be
! maxval(levels%ntpar(:) - levels%ntml(:) + 1), but for now
! we allow this to be input since it differs in certain
! places in the old code.
integer, intent(in) :: maxlevs         ! Maximum number of levels
                                       ! within region
integer, intent(in) :: maxtrlevs       ! Maximum number of levels
                                       ! for tracers within region
real(kind=real_umphys), intent(in) :: theta(:,:)
real(kind=real_umphys), intent(in) :: q_mix(:,:)
real(kind=real_umphys), intent(in) :: qse(:,:)
real(kind=real_umphys), intent(in) :: p_theta(:,0:)      ! input starts 0
real(kind=real_umphys), intent(in) :: p_rho(:,0:)        ! input starts 0
real(kind=real_umphys), intent(in) :: exner_theta(:,0:)  ! input starts 0
real(kind=real_umphys), intent(in) :: exner_rho(:,0:)    ! input starts 0
real(kind=real_umphys), intent(in) :: z_theta(:,:)
real(kind=real_umphys), intent(in) :: z_rho(:,:)
real(kind=real_umphys), intent(in) :: r2rho(:,:)
real(kind=real_umphys), intent(in) :: r2rho_theta(:,:)
real(kind=real_umphys), intent(in) :: rho(:,:)
real(kind=real_umphys), intent(in) :: tracer(:,:,:)

integer :: i, lbase, ltop, maxsize

character(len=*), parameter ::  RoutineName = 'ASSIGN_CLOUD_INPUT'
character(len=errormessagelength) :: Message
integer :: ErrorStatus ! Return code:
!   0 = Normal exit
! +ve = Fatal Error
! -ve = Warning

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

maxsize=size(theta(1,:))
do i=1,var%n_xx
  lbase = base_levels(i)
  ltop  = lbase + maxlevs - 1
  if (ltop >= maxsize) then
    ErrorStatus=1
    write(Message, '(A, 5I3)')                                                 &
       'Too many levels in warm convective layer: '//                          &
       'i, lbase, ltop, maxlevs, maxsize = '                                   &
       , i, lbase, ltop, maxlevs, maxsize

    call Ereport(RoutineName, ErrorStatus, Message)
  end if
  var%eta_theta(i,1:maxlevs) = (z_theta(i,lbase:ltop) - cb%z_rho(i))           &
                       /scales%zcld(i)
  var%eta_rho(i,1:maxlevs) = (z_rho(i,lbase:ltop) - cb%z_rho(i))               &
                     /scales%zcld(i)
  var%theta(i,1:maxlevs)       = theta(i,lbase:ltop)
  var%q_mix(i,1:maxlevs)       = q_mix(i,lbase:ltop)
  var%qse(i,1:maxlevs)         = qse(i,lbase:ltop)
  var%p_theta(i,1:maxlevs)     = p_theta(i,lbase:ltop)
  var%p_rho(i,1:maxlevs)       = p_rho(i,lbase:ltop)
  var%exner_theta(i,1:maxlevs) = exner_theta(i,lbase:ltop)
  var%exner_rho(i,1:maxlevs)   = exner_rho(i,lbase:ltop)
  var%z_theta(i,1:maxlevs)     = z_theta(i,lbase:ltop)
  var%z_rho(i,1:maxlevs)       = z_rho(i,lbase:ltop)
  var%r2rho(i,1:maxlevs)       = r2rho(i,lbase:ltop)
  var%r2rho_theta(i,1:maxlevs) = r2rho_theta(i,lbase:ltop)
  var%rho(i,1:maxlevs)         = rho(i,lbase:ltop)
  ltop  = lbase + maxtrlevs - 1
  var%tracer(i,1:maxtrlevs,:)    = tracer(i,lbase:ltop,:)
end do
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine assign_cloud_input

end module tcs_class_cloud
