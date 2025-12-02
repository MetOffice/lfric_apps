! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Autoconversion of snow.
module lsp_snow_autoc_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_SNOW_AUTOC_MOD'

contains

subroutine lsp_snow_autoc(                                                     &
  points, timestep,                                                            &
                                          ! Number of points and tstep
  qcf_cry, qcf_agg, t, cttemp,                                                 &
                                          ! Water contents and temp
  m0, t_scaling, cry_nofall, agg_nofall, qcf0,                                 &
                                          ! Parametrization information
  ptransfer,                                                                   &
                                          ! Mass transfer diagnostic
  one_over_tsi                                                                 &
                                          ! 1/(timestep*iterations)
  )

use lsprec_mod, only: zero, one

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,             only: real_lsprec

  ! Dr Hook Modules
use yomhook,          only: lhook, dr_hook
use parkind1,         only: jprb, jpim

implicit none

! Purpose:
!   Update cloud prognostics as a result of the autoconversion of
!   cloud ice to snow aggregates

!  Method:
!   Transfer some mass from ice crystals to snow aggregates depending
!   on the temperature and cloud top temperature.
!   Simple explicit Kessler type param. of autoconversion
!   with the autoconversion limit and rate set to emulate the
!   split-ice scheme.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Description of Code:
!   FORTRAN 77  + common extensions also in Fortran90.
!   This code is written to UMDP3 version 6 programming standards.

!   Documentation: UMDP 26.


! Subroutine Arguments

integer, intent(in) ::                                                         &
  points
                        ! Number of points to process

real (kind=real_lsprec), intent(in) ::                                         &
  timestep,                                                                    &
                        ! Timestep / s
  m0,                                                                          &
                        ! Seed ice mass / kg kg-1
  t_scaling,                                                                   &
                        ! Scaling temperature / K
  qcf0,                                                                        &
                        ! Prescribed ice content / kg kg-1
  t(points),                                                                   &
                        ! Temperature / K
  cttemp(points),                                                              &
                        ! Cloud top temperature / K
  cry_nofall(points),                                                          &
                        ! Fraction of qcf_cry that is not falliing out
  agg_nofall(points),                                                          &
                        ! Fraction of qcf_agg that is not falliing out
  one_over_tsi
                        ! 1/(timestep*iterations)

real (kind=real_lsprec), intent(in out) ::                                     &
  qcf_cry(points),                                                             &
                        ! Ice water content of ice crystals / kg kg-1
  qcf_agg(points),                                                             &
                        ! Ice water content of snow aggs. / kg kg-1
  ptransfer(points)
                        ! Autoconversion rate / kg kg-1 s-1

! Local Variables

integer ::                                                                     &
  i


real (kind=real_lsprec) ::                                                     &
  dpr,                                                                         &
                        ! Transfer amount from ice to snow / kg kg-1
  qcfautolim,                                                                  &
                        ! Autoconversion limit / kg kg-1
  qcfautorate,                                                                 &
                        ! Rate of transfer / s-1
  qc
                        ! Ice remaining after autoconversion
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_SNOW_AUTOC'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
do i = 1, points

  if (qcf_cry(i)*cry_nofall(i) > m0) then

        !-----------------------------------------------
        ! Set autoconversion limit to emulate split-ice scheme
        !-----------------------------------------------
    qcfautolim = (qcf_agg(i)+qcf_cry(i))                                       &
               *max(exp(-t_scaling * max((t(i)-cttemp(i)),zero)                &
               *max(qcf_agg(i)*agg_nofall(i) + qcf_cry(i)*cry_nofall(i),       &
                    zero)*qcf0) , zero)

        !-----------------------------------------------
        ! Set rate to emulate spilt-ice scheme, i.e. infinite
        !-----------------------------------------------
    qcfautorate = one/timestep

    qc  = min(qcfautolim , qcf_cry(i))
    dpr = min(qcfautorate * timestep * (qcf_cry(i)-qc), qcf_cry(i) - qc)

        !-----------------------------------------------
        ! Store process rate (kg kg-1 s-1)
        !-----------------------------------------------
    ptransfer(i) = ptransfer(i) + dpr * one_over_tsi

        !-----------------------------------------------
        ! Update ice/snow variables
        !-----------------------------------------------
    qcf_cry(i) = qcf_cry(i) - dpr
    qcf_agg(i) = qcf_agg(i) + dpr

        ! No cloud fraction updating is needed

  end if  ! qcf_cry > 0

end do  ! Points

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_snow_autoc
end module lsp_snow_autoc_mod
