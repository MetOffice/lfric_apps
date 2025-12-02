! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!    subroutine LSP_SCAV-----------------------------------------------
!    Purpose: Scavenge aerosol by large scale precipitation.
!
!    Programming standard: Unified Model Documentation Paper No 3,
!                          Version 7, dated 11/3/93.
!
!    Logical component covered: Part of P26.
!
!    Documentation: Unified Model Documentation Paper No 26.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation
module lsp_scav_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_SCAV_MOD'

contains

subroutine lsp_scav( points, rain, snow, droplet_flux, aerosol )

use lsprec_mod, only: timestep, one

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,              only: real_lsprec

! Dr Hook Modules
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none
integer ::                                                                     &
                    ! Input integer scalar :-
 points         ! in Number of points to be processed.

real (kind=real_lsprec) ::                                                     &
                    ! Input real arrays :-
 rain(points),                                                                 &
                    ! in Rate of rainfall in this layer from
!                     !       above
!                     !       (kg per sq m per s).
   snow(points),                                                               &
                      ! in Rate of snowfall in this layer from
!                     !       above
!                     !       (kg per sq m per s).
   droplet_flux(points)
                      ! In Rate of droplet settling in this layer
                      !       from above
                      !       (kg per sq m per s).
real (kind=real_lsprec) ::                                                     &
                    ! Updated real arrays :-
 aerosol(points) ! INOUT Aerosol mixing ratio


!-----------------------------------------------------------------------
!  Define local scalars.
!-----------------------------------------------------------------------
!  (a) Reals effectively expanded to workspace by the Cray (using
!      vector registers).
real (kind=real_lsprec) ::                                                     &
                    ! Real workspace.
 krain,ksnow
parameter(krain=1.0e-4_real_lsprec,ksnow=1.0e-4_real_lsprec)
real (kind=real_lsprec) ::                                                     &
                    ! Real workspace.
 rrain,rsnow
!  (b) Others.
integer :: i       ! Loop counter (horizontal field index).

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_SCAV'


! Overall rate = KRAIN*(R) where R is in mm/hr=kg/m2/s*3600.0
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

rrain=krain*timestep*3600.0_real_lsprec
rsnow=ksnow*timestep*3600.0_real_lsprec
do i = 1, points
  aerosol(i)=aerosol(i)/                                                       &
        (one+rrain*rain(i)+rrain*droplet_flux(i)+rsnow*snow(i))
end do
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_scav
end module lsp_scav_mod
