! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module lsp_settle_wtrac_mod

use um_types, only: real_lsprec

implicit none

! Description:
!   Update water tracers for droplet settling.
!
! Method:
!   Code written for i_fix_mphys_drop_settle == first_fix or second_fix
!   options in lsp_settle.
!   Note, fractionation will occur here for water isotopes, due to
!   evaporation occuring.
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private :: ModuleName = 'LSP_SETTLE_WTRAC_MOD'
contains

! Subroutine Interface:
subroutine lsp_settle_wtrac(points, cfliq, rhor, dhi, q, qcl, t,               &
                            droplet_flux, wtrac_mp_cpr_old, wtrac_mp_cpr)

use science_fixes_mod,       only: i_fix_mphys_drop_settle, first_fix,         &
                                   second_fix
use lsprec_mod,              only: mprog_min, zero, one
use wtrac_calc_ratio_mod,    only: wtrac_calc_ratio_fn
use free_tracers_inputs_mod, only: n_wtrac
use wtrac_mphys_mod,         only: mp_cpr_old_wtrac_type, mp_cpr_wtrac_type

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer, intent(in) :: points      ! No. of points

real (kind=real_lsprec), intent(in) ::                                         &
  cfliq(points),                                                               &
                        ! Liquid cloud fraction at start of timestep
  rhor(points),                                                                &
                        ! 1 / Air density / m3 kg-1
  dhi(points),                                                                 &
                        ! Timestep / layer thickness / s m-1
  q(points),                                                                   &
                        ! Water vapour after settling
  qcl(points),                                                                 &
                        ! Liquid condensate after settling
  t(points),                                                                   &
                        ! Temperature after settling
  droplet_flux(points)
                        ! Droplet flux out of layer (as just calculated in
                        !   lsp_settle)

! Structure containing normal water fields before call to lsp_settle
type(mp_cpr_old_wtrac_type), intent(in out) :: wtrac_mp_cpr_old

! Structure containing water tracer fields
type(mp_cpr_wtrac_type), intent(in out) :: wtrac_mp_cpr(n_wtrac)

! Local variables
integer :: i, m, i_wt         ! Loop counters
integer :: ni                 ! No. of working points
integer :: idx(points)        ! Index for working points

real(kind=real_lsprec) :: flux_into_cloud      ! Flux of tracer in liquid
                                               ! water from layer above
                                               ! falling into
                                               ! cloud (kg m-2 s-1)
real(kind=real_lsprec) :: flux_into_clear_sky  ! Flux of tracer in liquid
                                               ! from layer above falling into
                                               ! clear sky (kg m-2 s-1)
real(kind=real_lsprec) :: droplet_flux_wtrac   ! Water tracer flux out of
                                               ! layer (kg m-2 s-1)
real(kind=real_lsprec) :: dqcl_wtrac           ! Change in tracer qcl this
                                               ! timestep (kg kg-1)
real(kind=real_lsprec) :: dq_wtrac             ! Change in tracer q this
                                               ! timestep (kg kg-1)

real(kind=real_lsprec) :: qcl_ratio            ! Ratio of water tracer
                                               ! to normal water qcl

real(kind=real_lsprec) :: t_mean(points)       ! Mean T during settling

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_SETTLE_WTRAC'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


! Set up working points
ni = 0

select case(i_fix_mphys_drop_settle)
case (first_fix)
  do i = 1, points
    if (wtrac_mp_cpr_old%qcl(i) > mprog_min) then
      ni = ni + 1
      idx(ni) = i
    end if
  end do

case (second_fix)
  do i = 1, points
    if (wtrac_mp_cpr_old%qcl(i) > mprog_min .or.                               &
        wtrac_mp_cpr_old%droplet_flux(i) > zero) then
      ni = ni + 1
      idx(ni) = i
    end if
  end do

  !case DEFAULT
  ! Not needed as check for setting in wtrac_check_setup

end select

! Calculate mean T during settling - for use with isotopes
do i = 1, points
  t_mean(i) = 0.5 * ( t(i) + wtrac_mp_cpr_old%t(i) )
end do

do i_wt = 1, n_wtrac

  do m = 1, ni   ! Only loop over work points
    i = idx(m)

    ! Calculate water tracer flux out of layer
    qcl_ratio = wtrac_calc_ratio_fn(i_wt, wtrac_mp_cpr(i_wt)%qcl(i),           &
                                    wtrac_mp_cpr_old%qcl(i))
    droplet_flux_wtrac = droplet_flux(i) * qcl_ratio

    ! Water tracer flux into layer split between cloud and clear sky
    flux_into_cloud     = cfliq(i) * wtrac_mp_cpr(i_wt)%droplet_flux(i)
    flux_into_clear_sky = (one - cfliq(i)) * wtrac_mp_cpr(i_wt)%droplet_flux(i)

    ! Change in qcl = flux into cloudy region - flux out
    dqcl_wtrac = (flux_into_cloud - droplet_flux_wtrac) * dhi(i)*rhor(i)

    ! Change in q due to qcl falling into clear sky and evaporating
    dq_wtrac   =  flux_into_clear_sky * dhi(i) * rhor(i)

    if (qcl(i) == zero) then
      ! If qcl + dqcl drops below mprog_abs after settling it will have been
      !  set to zero in lsp_settle, so do similarly for water tracers

      if (i_fix_mphys_drop_settle == second_fix) then
        droplet_flux_wtrac = droplet_flux_wtrac                                &
           + ( wtrac_mp_cpr(i_wt)%qcl(i) + dqcl_wtrac ) / ( dhi(i) * rhor(i) )
      end if

      wtrac_mp_cpr(i_wt)%qcl(i) = zero
    else if ( wtrac_mp_cpr(i_wt)%qcl(i) + dqcl_wtrac < 0.0 ) then
      ! Check that there is enough water tracer liquid, if not reduce flux
      ! out of current level consistent with reduced increment

      if (i_fix_mphys_drop_settle == second_fix) then
        droplet_flux_wtrac = droplet_flux_wtrac                                &
           + ( wtrac_mp_cpr(i_wt)%qcl(i) + dqcl_wtrac ) / ( dhi(i) * rhor(i) )
      end if
      wtrac_mp_cpr(i_wt)%qcl(i) = zero

    else
      ! Update water tracer qcl
      wtrac_mp_cpr(i_wt)%qcl(i) = wtrac_mp_cpr(i_wt)%qcl(i) + dqcl_wtrac
    end if

    ! Set droplet_flux to flux out of layer
    wtrac_mp_cpr(i_wt)%droplet_flux(i) = droplet_flux_wtrac

    ! Update water tracer q
    wtrac_mp_cpr(i_wt)%q(i) = wtrac_mp_cpr(i_wt)%q(i) + dq_wtrac

  end do    ! m
end do      ! i_wt

! If water tracers are isotopes, carry out equilibrium fractionation between
! q and qcl here due to evaporation

! Store current q and qcl values
do i = 1, points
  wtrac_mp_cpr_old%q(i)   = q(i)
  wtrac_mp_cpr_old%qcl(i) = qcl(i)
  wtrac_mp_cpr_old%t(i)   = t(i)
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_settle_wtrac

end module lsp_settle_wtrac_mod
