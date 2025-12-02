! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module lsp_init_wtrac_mod

use um_types, only: real_lsprec

implicit none
! Description:
!   Set up required water tracer fields at the start of lsp_ice
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private :: ModuleName = 'LSP_INIT_WTRAC_MOD'
contains

! Subroutine Interface:
subroutine lsp_init_wtrac(ip, timestep, rhodz_dry, rhodz_moist,                &
                           q, qcl, qcf, qrain, t, droplet_flux,                &
                           wtrac_mp_cpr_old, wtrac_mp_cpr)

use mphys_inputs_mod,        only: l_mcr_qrain
use lsprec_mod,              only: zero
use gen_phys_inputs_mod,     only: l_mr_physics
use free_tracers_inputs_mod, only: n_wtrac
use wtrac_mphys_mod,         only: mp_cpr_old_wtrac_type, mp_cpr_wtrac_type

use yomhook,                 only: lhook, dr_hook
use parkind1,                only: jprb, jpim

implicit none

integer, intent(in) :: ip         ! No. of points

real (kind=real_lsprec), intent(in) :: timestep

real (kind=real_lsprec), intent(in) :: rhodz_dry(ip)
                         ! Dry air density*deltaz for dry air (kg m-2)
real (kind=real_lsprec), intent(in) ::  rhodz_moist(ip)
                         ! Moist air density*deltaz for moist air (kg m-2)

! 'Normal' water fields
real (kind=real_lsprec), intent(in) :: q(ip)
real (kind=real_lsprec), intent(in) :: qcl(ip)
real (kind=real_lsprec), intent(in) :: qcf(ip)
real (kind=real_lsprec), intent(in) :: qrain(ip)
real (kind=real_lsprec), intent(in) :: t(ip)
real (kind=real_lsprec), intent(in) :: droplet_flux(ip)

! Structure to store old normal water fields - used for water tracer
! calculations
type(mp_cpr_old_wtrac_type), intent(in out) :: wtrac_mp_cpr_old

! Water tracer fields
type(mp_cpr_wtrac_type), intent(in out) :: wtrac_mp_cpr(n_wtrac)

! Local variables

integer :: i, i_wt    ! Loop counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_INIT_WTRAC'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Allocate arrays in structure to store old values of 'normal' water and
! sublimation.  (These are deallocated in lsp_tidy.)
allocate(wtrac_mp_cpr_old%q(ip))
allocate(wtrac_mp_cpr_old%qcl(ip))
allocate(wtrac_mp_cpr_old%qcf(ip))
allocate(wtrac_mp_cpr_old%qrain(ip))
allocate(wtrac_mp_cpr_old%t(ip))
allocate(wtrac_mp_cpr_old%droplet_flux(ip))
allocate(wtrac_mp_cpr_old%qchange(ip))
allocate(wtrac_mp_cpr_old%het_q(ip))
allocate(wtrac_mp_cpr_old%hom_qr(ip))
allocate(wtrac_mp_cpr_old%depos_v(ip))
allocate(wtrac_mp_cpr_old%depos_l(ip))


do i = 1, ip
  ! Set to current values
  wtrac_mp_cpr_old%q(i)            = q(i)
  wtrac_mp_cpr_old%qcl(i)          = qcl(i)
  wtrac_mp_cpr_old%qcf(i)          = qcf(i)
  wtrac_mp_cpr_old%qrain(i)        = qrain(i)
  wtrac_mp_cpr_old%t(i)            = t(i)
  wtrac_mp_cpr_old%droplet_flux(i) = droplet_flux(i)

  ! Set store of water change amount to zero
  wtrac_mp_cpr_old%qchange(i) = zero
  wtrac_mp_cpr_old%het_q(i)   = zero
  wtrac_mp_cpr_old%hom_qr(i)  = zero
  wtrac_mp_cpr_old%depos_v(i) = zero
  wtrac_mp_cpr_old%depos_l(i) = zero
end do

! Set fluxes to zero
do i_wt = 1, n_wtrac
  do i = 1, ip
    wtrac_mp_cpr(i_wt)%snowratet(i) = zero
    wtrac_mp_cpr(i_wt)%rainratet(i) = zero
  end do
end do

if (.not. l_mcr_qrain) then
  ! Rain is a diagnostic quantity
  do i_wt = 1, n_wtrac
    do i = 1, ip
      if (wtrac_mp_cpr(i_wt)%lsrain(i)  >   zero) then

        if (l_mr_physics) then

          ! Mixing ratio formulation
          wtrac_mp_cpr(i_wt)%qrain(i) =                                        &
               wtrac_mp_cpr(i_wt)%lsrain(i) * timestep / rhodz_dry(i)

        else ! l_mr_physics

          ! Specific humidity formulation
          wtrac_mp_cpr(i_wt)%qrain(i) =                                        &
               wtrac_mp_cpr(i_wt)%lsrain(i) * timestep / rhodz_moist(i)

        end if  ! l_mr_physics

      else    ! rainrate > 0

        wtrac_mp_cpr(i_wt)%qrain(i) = zero

      end if  ! rainrate > 0
    end do    ! ip
  end do      ! water tracers
end if        ! .not. l_mcr_qrain


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_init_wtrac

end module lsp_init_wtrac_mod
