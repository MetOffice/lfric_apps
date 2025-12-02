! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module lsp_sedim_eulexp_wtrac_mod

use um_types, only: real_lsprec

implicit none

! Description:
!
!   First order Eulerian advection scheme for water tracer sedimentation
!   with an exponential-based limiter to ensure stability for CFL>1
!   used for ice/snow and rain (as in lsp_sedim_eulexp)
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private ::                                        &
                              ModuleName = 'LSP_SEDIM_EULEXP_WTRAC_MOD'
contains

! Subroutine Interface:
subroutine lsp_sedim_eulexp_wtrac(points, rho, rhor, dhi, dhir, t,             &
                            fallspeed_thislayer, q, water_type,                &
                            wtrac_mp_cpr_old, wtrac_mp_cpr)

use lsprec_mod,              only: zero, one
use free_tracers_inputs_mod, only: n_wtrac
use wtrac_mphys_mod,         only: mp_cpr_old_wtrac_type, mp_cpr_wtrac_type

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

integer, intent(in) :: points      ! No. of points

real (kind=real_lsprec), intent(in) ::                                         &
  rho(points),                                                                 &
                        ! Air density / kg m-3
  rhor(points),                                                                &
                        ! 1 / Air density / m3 kg-1
  dhi(points),                                                                 &
                        ! Timestep / layer thickness / s m-1
  dhir(points),                                                                &
                        ! 1 / dhi / m s-1
  t(points),                                                                   &
                        ! Temperature after call to lsp_fall (this is only
                        ! passed to this routine to keep it up to date in
                        ! wtrac_mp_cpr_old)
  q(points),                                                                   &
                        ! 'Normal' water field after fall calculation / kg kg-1
                        ! either qcl or qrain
  fallspeed_thislayer(points)
                        ! Ice or rain fall speed leaving the current layer /
                        !   m s-1

character( * ), intent(in) :: water_type   ! 'ice' or 'rain'

! Structure containing normal water fields before call to lsp_settle
type(mp_cpr_old_wtrac_type), intent(in out) :: wtrac_mp_cpr_old

! Structure containing water tracer fields
type(mp_cpr_wtrac_type), intent(in out) :: wtrac_mp_cpr(n_wtrac)

! Local variables
integer :: i, i_wt         ! Loop counters

real(kind=real_lsprec) :: q_old(points)  ! 'Normal' water value before
                                         ! sedimentation

real(kind=real_lsprec) :: mixratio_thislayer  ! Water tracer for the layer

real(kind=real_lsprec) :: flux_fromabove      ! Water tracer flux from above

real(kind=real_lsprec) :: expfactor           ! Fraction of ice that remains
                                              ! in layer

real(kind=real_lsprec) :: flux_out_wtrac      ! Water tracer snow fall rate
                                              ! out of layer


integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_SEDIM_EULEXP_WTRAC'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Water tracer fields
if (water_type == 'ice') then
  do i = 1, points
    q_old(i) = wtrac_mp_cpr_old%qcf(i)
    ! Store current qcf values for use in future routines
    wtrac_mp_cpr_old%qcf(i) = q(i)
    wtrac_mp_cpr_old%t(i)   = t(i)
  end do
else if (water_type == 'rain') then
  do i = 1, points
    q_old(i) = wtrac_mp_cpr_old%qrain(i)
    ! Store current qrain values for use in future routines
    wtrac_mp_cpr_old%qrain(i) = q(i)
    wtrac_mp_cpr_old%t(i)     = t(i)
  end do
end if

do i = 1, points

  do i_wt = 1, n_wtrac

    ! Set up water tracer fields
    if (water_type == 'ice') then
      mixratio_thislayer = wtrac_mp_cpr(i_wt)%qcf(i)
      flux_fromabove     = wtrac_mp_cpr(i_wt)%lssnow(i)
    else if (water_type == 'rain') then
      mixratio_thislayer = wtrac_mp_cpr(i_wt)%qrain(i)
      flux_fromabove     = wtrac_mp_cpr(i_wt)%lsrain(i)
    end if

    if (fallspeed_thislayer(i) > zero) then

      ! Fraction of ice or rain that remains in the same layer
      expfactor = exp(-one*fallspeed_thislayer(i)*dhi(i))

      ! Water tracer flux out of this layer
      flux_out_wtrac =  flux_fromabove + dhir(i) *                             &
                      (rho(i)*mixratio_thislayer                               &
                         - flux_fromabove/fallspeed_thislayer(i))              &
                      * (one-expfactor)

      ! Update water tracer field for the layer
      mixratio_thislayer = flux_fromabove * rhor(i)                            &
                           / fallspeed_thislayer(i)                            &
                           * (one-expfactor)                                   &
                           +  mixratio_thislayer*expfactor

    else

      ! No fall of ice out of layer (so no ice initially in this layer)
      flux_out_wtrac = zero
      mixratio_thislayer = flux_fromabove*rhor(i)*dhi(i)

    end if

    ! Update main water tracer structure
    if (water_type == 'ice') then
      wtrac_mp_cpr(i_wt)%qcf(i)       = mixratio_thislayer
      wtrac_mp_cpr(i_wt)%lssnow(i)    = flux_fromabove
      wtrac_mp_cpr(i_wt)%snowratet(i) = wtrac_mp_cpr(i_wt)%snowratet(i)        &
                                       + flux_out_wtrac
    else if (water_type == 'rain') then
      wtrac_mp_cpr(i_wt)%qrain(i)     = mixratio_thislayer
      wtrac_mp_cpr(i_wt)%lsrain(i)    = flux_fromabove
      wtrac_mp_cpr(i_wt)%rainratet(i) = wtrac_mp_cpr(i_wt)%rainratet(i)        &
                                       + flux_out_wtrac
    end if

  end do   ! i_wt
end do     ! i


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_sedim_eulexp_wtrac

end module lsp_sedim_eulexp_wtrac_mod

