! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
! Module to calculate pc2 increments
!
module tcs_pc2


use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
use um_types, only: real_umphys

implicit none
!
! Description:
! Module to calculate pc2 increments
!
! Method:
!   Currenly uses a fairly ad hoc treatment based on in-cloud
!   liquid water flux and detrainment rates defined in
!   tcs_parameters_warm
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: Fortran 90.
!   This code is written to UMDP3 version 8.1 programming standards.
!

character(len=*), parameter, private :: ModuleName='TCS_PC2'

contains

subroutine calc_pc2( n_xx, max_cldlev, timestep, ntml, ntpar,                  &
   mf_h_cld, ql_up_cld, wql_cld, qcl, cf_liquid,                               &
   cld_in, sim, dqclbydt, dqcfbydt, dbcfbydt, dcflbydt, dcffbydt)

  !-------------------------------------------------------------------
  !
  ! Description:
  !   Calculate the increments to PC2 variables.
  !
  !-------------------------------------------------------------------

use tcs_parameters_warm,    only: Aepsilon
use tcs_class_similarity,   only: similarity
use tcs_class_cloud,        only: cloud_input
use tcs_common_warm,        only: scales

implicit none
!----------------------------------------------------------------
! Subroutine Arguments
!----------------------------------------------------------------
integer, intent(in) ::                                                         &
   n_xx                                                                        &
                      ! No. of congestus convection points
   , max_cldlev
                      ! Maximum number of convective cloud levels

real(kind=real_umphys), intent(in) ::                                          &
   timestep
                            ! Model timestep (s)

integer, intent(in) ::                                                         &
   ntml(:)                                                                     &
                            ! Top level of surface mixed layer defined
                            ! relative to theta,q grid
   ,ntpar(:)
! Top level of surface mixed layer defined
! relative to theta,q grid

real(kind=real_umphys), intent(in) ::                                          &
   mf_h_cld(:,:)                                                               &
                            ! mass flux in cld (rho levels)
   ,ql_up_cld(:,:)                                                             &
                            ! liquid water content in the updraught
   ,wql_cld(:,:)
                            ! liquid water flux

real(kind=real_umphys), intent(in out) ::                                      &
   qcl(:,:)                                                                    &
                            ! Liq condensate mix ratio (kg/kg)
   , cf_liquid(:,:)
                            ! Frozen water cloud volume
! Liq water cloud volume

type(cloud_input), intent(in) :: cld_in ! fields on cloud levels

type(similarity), intent(in) :: sim !Similarity functions

real(kind=real_umphys), intent(in out) ::                                      &
   dqclbydt(:,:)                                                               &
   ,dqcfbydt(:,:)                                                              &
   ,dbcfbydt(:,:)                                                              &
   ,dcflbydt(:,:)                                                              &
   ,dcffbydt(:,:)

!-------------------------------------------------------------------
! Variables defined locally
!-------------------------------------------------------------------

real(kind=real_umphys), allocatable :: dmf(:,:)   ! mass flux divergence
real(kind=real_umphys), allocatable :: entr(:,:)  ! entrainment rate
real(kind=real_umphys), allocatable :: detr(:,:)  ! detrainment rate

!-------------------------
! Loop counters
!-------------------------
integer :: i,k,lbase

!=================
! Select method
!=================
integer, parameter :: ipc2_method=2

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CALC_PC2'

!-------------------------------------------------------------------
! 1.0 Allocate and initialise arrays
!-------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

allocate(dmf(n_xx,max_cldlev))
allocate(detr(n_xx,max_cldlev))
allocate(entr(n_xx,max_cldlev))

dqclbydt(:,:) = 0.0
dqcfbydt(:,:) = 0.0
dbcfbydt(:,:) = 0.0
dcflbydt(:,:) = 0.0
dcffbydt(:,:) = 0.0

if (ipc2_method==1) then

  !-----------------------------------------------------------------
  ! 2.0 calculate mass flux divergence
  !-----------------------------------------------------------------

  do k=1,max_cldlev-1
    do i = 1,n_xx
      dmf(i,k) = (mf_h_cld(i,k+1) - mf_h_cld(i,k))/                            &
         (cld_in%z_rho(i,k+1)-cld_in%z_rho(i,k))
    end do
  end do

  k=max_cldlev
  do i = 1,n_xx
    dmf(i,k) = mf_h_cld(i,k)/cld_in%z_rho(i,k)
  end do


  !-----------------------------------------------------------------
  ! 3.0 calculate entrainment and detrainment rates
  !-----------------------------------------------------------------

  do k=1,max_cldlev-1
    do i = 1,n_xx
      entr(i,k) = Aepsilon*scales%wstar_up(i)/scales%mb(i)                     &
         /(scales%zcld(i) * cld_in%eta_theta(i,k))
      detr(i,k) = entr(i,k)-dmf(i,k)
    end do
  end do

  !-----------------------------------------------------------------
  ! 3.1 Put sensible limits on detr
  !-----------------------------------------------------------------

  do k=1,max_cldlev-1
    do i = 1,n_xx
      if (k > ntpar(i) - ntml(i) + 1                                           &
         .or. k < ntml(i) + 1                                                  &
         .or. detr(i,k) < 0.0) detr(i,k) = 0.0
      if (detr(i,k) > .5) detr(i,k) = .5
    end do
  end do


  !-----------------------------------------------------------------
  ! 4.0 calculate detrained water
  !-----------------------------------------------------------------

  ! working note: for now we just calculate detrained liquid - if the
  ! scheme is being run in a cold environment this may not be
  ! appropriate.
  ! working note: we don't detrain anything in the lowest cloud level
  do i = 1,n_xx
    do k=2,ntpar(i)-ntml(i)
      lbase = ntml(i)+k
      dqclbydt(i,lbase) = detr(i,k)*(ql_up_cld(i,k-1)-qcl(i,lbase))            &
         /cld_in%rho(i,k)
      dcflbydt(i,lbase) = detr(i,k)*(1.0-cf_liquid(i,lbase))/cld_in%rho(i,k)
    end do
    ! Enhanced detrainment at cloud top
    k=ntpar(i)-ntml(i)+1
    lbase = ntml(i)+k
    dqclbydt(i,lbase) = detr(i,k)*(ql_up_cld(i,k-1))/cld_in%rho(i,k)
    dcflbydt(i,lbase) = detr(i,k)*(1.0-cf_liquid(i,k))/cld_in%rho(i,k)
  end do


else if (ipc2_method==2) then
  do i = 1,n_xx
    do k=1,ntpar(i)-ntml(i)+1
      lbase = ntml(i)+k
      detr(i,k) = sim%pc2_detr(i,k)
      dqclbydt(i,lbase) = detr(i,k)*wql_cld(i,k)
    end do
  end do
end if

!-------------------------------------------------------------------
! 5.0 Put sensible limits on increments and work out change in cloud
!     fraction
!-------------------------------------------------------------------

do i = 1,n_xx
  do k=1,ntpar(i)-ntml(i)
    lbase = ntml(i)+k
    dqclbydt(i,lbase) = max(dqclbydt(i,lbase), -qcl(i,lbase)/timestep)
    ! dqclbydt(i,lbase) = min(dqclbydt(i,lbase), ql_up_cld(i,k)/timestep)
    dcflbydt(i,lbase) = dqclbydt(i,lbase)*                                     &
       (1.0-cf_liquid(i,k))/(ql_up_cld(i,k)-qcl(i,lbase))
    dcflbydt(i,lbase) = max(min(dcflbydt(i,lbase),1.0),0.0)
  end do
  k=ntpar(i)-ntml(i)+1
  lbase = ntml(i)+k
  dqclbydt(i,lbase) = max(dqclbydt(i,lbase), -qcl(i,lbase)/timestep)
  ! dqclbydt(i,lbase) = min(dqclbydt(i,lbase), ql_up_cld(i,k-1)/timestep)
  dcflbydt(i,lbase) = dqclbydt(i,lbase)*                                       &
     (1.0-cf_liquid(i,k))/(ql_up_cld(i,k)-qcl(i,lbase) + 1e-10)
  dcflbydt(i,lbase) = max(min(dcflbydt(i,lbase),1.0),0.0)
end do

!-------------------------------------------------------------------
! 6.0 Deallocation
!-------------------------------------------------------------------

deallocate(entr)
deallocate(detr)
deallocate(dmf)
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return

end subroutine calc_pc2

end module tcs_pc2
