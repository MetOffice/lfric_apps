! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module environ_wtrac_mod

use um_types, only: real_umphys

implicit none

! Description:
!   Calculate the effect of convection on the water tracers in the large
!   scale atmosphere
!
! Method:
!  The water tracers in the large scale atmosphere are updated for the
!  effect of compensating subsidence and detrainment in the same way as normal
!  water, but with one difference. The difference is due to the implied phase
!  change that occurs in environ_6a, where the condensate increments are
!  partitioned between the liquid and ice phase - allowing the condensate
!  detrained by convection to change phase smoothly from liquid to ice,
!  depending on temperature of environment, rather than being detrained with
!  the same phase as the in the convective plume, which changes abruptly.
!  This is modelled as a normal phase change for water tracers,
!  so it is based on the amount of water changing phase.
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private :: ModuleName = 'ENVIRON_WTRAC_MOD'

contains

! Subroutine Interface:
subroutine environ_wtrac(k, npnts, ni, n_wtrac,                                &
                         idx, timestep, a_smth, deltak,                        &
                         amdetk, flxbydpk, ekp14, delpk, delpkp1,              &
                         qclek_temp, qcfek_temp,                               &
                         dqclek_frz, dqcfek_mlt, dqclekp1_frz, dqcfekp1_mlt,   &
                         qrk_wtrac, wtrac_p, wtrac_e)

use wtrac_conv_mod,         only: conv_e_wtrac_type, conv_p_wtrac_type
use wtrac_calc_ratio_mod,   only: wtrac_calc_ratio_fn
use wtrac_move_phase_mod,   only: wtrac_move_phase

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
integer,intent(in) :: k          ! present model layer
integer,intent(in) :: npnts      ! Number of points
integer,intent(in) :: ni         ! No. of work points
integer,intent(in) :: n_wtrac    ! Number of water tracer variables

integer,intent(in) :: idx(npnts) ! Indices of work points

real(kind=real_umphys),intent(in) :: timestep   ! Timestep
real(kind=real_umphys),intent(in) :: a_smth
                                 ! Parameter determining the weighting
                                 ! between the forced
                                 ! detrainment increments at k and k-1

real(kind=real_umphys),intent(in) :: deltak(npnts)
                                 ! Parcel forced detrainment rate in
                                 ! layer k multiplied by layer thickness
real(kind=real_umphys),intent(in) :: amdetk(npnts)
                                 ! Mixing detrainment coefficient at level k
                                 ! multiplied by appropriate layer thickness
real(kind=real_umphys),intent(in) :: flxbydpk(npnts)
                                 ! mass flux divided by layer thickness (1/s)
real(kind=real_umphys),intent(in) :: ekp14(npnts)
                                 ! Entrainment coefficient at level k+1/4
                                 ! multiplied by appropriate layer thickness
real(kind=real_umphys),intent(in) :: delpk(npnts)
                                 ! pressure difference across layer k (Pa)
real(kind=real_umphys),intent(in) :: delpkp1(npnts)
                                 ! pressure difference across layer k+1 (Pa)

real(kind=real_umphys),intent(in) :: qclek_temp(npnts)
                                 ! qclek updated in environ but before phase
                                 ! change
real(kind=real_umphys),intent(in) :: qcfek_temp(npnts)
                                 ! qcfek updated in environ but before phase
                                 ! change

real(kind=real_umphys),intent(in) :: dqclek_frz(npnts)
                                 ! Freezing rate of liquid condensate at k
real(kind=real_umphys),intent(in) :: dqcfek_mlt(npnts)
                                 ! Melting rate of ice condensate at level k
real(kind=real_umphys),intent(in) :: dqclekp1_frz(npnts)
                                 ! Freezing rate of liquid condensate at k+1
real(kind=real_umphys),intent(in) :: dqcfekp1_mlt(npnts)
                                 ! Melting rate of ice condensate at k+1

real(kind=real_umphys),intent(in) :: qrk_wtrac(npnts,n_wtrac)
                                 ! Water tracer specific humidity of forced
                                 ! detrained parcel in layer k (kg/kg)

type(conv_p_wtrac_type), intent(in) :: wtrac_p(n_wtrac)
                                 ! Structure containing parcel
                                 ! water tracer fields

type(conv_e_wtrac_type), intent(in out) :: wtrac_e(n_wtrac)
                                 ! Structure containing environment
                                 ! water tracer fields

! Local variables
integer :: i, m, i_wt

real(kind=real_umphys) :: tmp_fd_dqek_wtrac
                                 ! Water tracer forced detrainment q inc
                                 ! across levels k and k+1
real(kind=real_umphys) :: tmp_fd_dqclek_wtrac
                                 ! Water tracer forced detrainment qcl inc
                                 ! across levels k and k+1
real(kind=real_umphys) :: tmp_fd_dqcfek_wtrac
                                 ! Water tracer forced detrainment qcf inc
                                 ! across levels k and k+1
real(kind=real_umphys) :: tmp_dqclek_wtrac ! Storage space for qcl rate
real(kind=real_umphys) :: tmp_dqcfek_wtrac ! Storage space for qcf rate
real(kind=real_umphys) :: qclek_temp_wtrac(npnts, n_wtrac)
                                 ! Water tracer qclek updated for change
real(kind=real_umphys) :: qcfek_temp_wtrac(npnts, n_wtrac)
                                 ! Water tracer qclek updated for change
real(kind=real_umphys) :: qchange_wtrac(npnts)
                                 ! Water tracer change

! Ratios of water tracer to normal water for forced detrainment:
real(kind=real_umphys) :: ratio_qclk(npnts,n_wtrac)       ! qcl at k
real(kind=real_umphys) :: ratio_qcfk(npnts,n_wtrac)       ! qcf at k

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='ENVIRON_WTRAC'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Initialise ratio fields
do i_wt = 1, n_wtrac
  do i = 1, npnts
    ratio_qclk(i,i_wt) = 0.0
    ratio_qcfk(i,i_wt) = 0.0
  end do
end do

do i_wt = 1, n_wtrac
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)

    ! Smoothed forced detrainment term (assuming dq_detcond=0)
    tmp_fd_dqek_wtrac   = deltak(i) * (1.0-amdetk(i))                          &
                  * ( qrk_wtrac(i,i_wt) - wtrac_e(i_wt)%q(i,k) )

    ! q increment at level k (assuming dq_detcond = 0)
    wtrac_e(i_wt)%dqbydt(i,k) = wtrac_e(i_wt)%dqbydt(i,k) + flxbydpk(i)        &
    !           Compensating subsidence
                * ( (1.0-amdetk(i))*(1.0-deltak(i))*(1.0+ekp14(i))             &
                * ( wtrac_e(i_wt)%q(i,k+1) - wtrac_e(i_wt)%q(i,k) )            &
    !           Smoothed forced detrainment
                + a_smth*tmp_fd_dqek_wtrac                                     &
    !           Mixing detrainment
                + amdetk(i)                                                    &
                * ( wtrac_p(i_wt)%q(i,k) - wtrac_e(i_wt)%q(i,k) ) )

    ! q increment at level k+1 due to smoothed  forced detrainment
    wtrac_e(i_wt)%dqbydt(i,k+1) = flxbydpk(i) * (1.0-a_smth)                   &
                  * delpk(i)/delpkp1(i)*tmp_fd_dqek_wtrac
  end do ! m
end do  ! i_wt

do i_wt = 1, n_wtrac
!DIR$ IVDEP
  do m=1, ni
    i = idx(m)

    ! -------------------------------------------------------------------
    ! QCL calculations
    ! -------------------------------------------------------------------

    ! Smoothed forced detrainment term for water tracer qcl at level k
    tmp_fd_dqclek_wtrac = deltak(i) * (1.0-amdetk(i))                          &
                * ( wtrac_p(i_wt)%qcl(i,k) - wtrac_e(i_wt)%qcl(i,k) )

    ! Water tracer qcl increment at level k
    tmp_dqclek_wtrac  = flxbydpk(i)                                            &
                ! Compensating subsidence
                * ( (1.0-amdetk(i))*(1.0-deltak(i))*(1.0+ekp14(i))             &
                * ( wtrac_e(i_wt)%qcl(i,k+1) -  wtrac_e(i_wt)%qcl(i,k) )       &
                ! Smoothed forced detrainment
                + a_smth*tmp_fd_dqclek_wtrac                                   &
                ! Mixing detrainment
                + amdetk(i)                                                    &
                * ( wtrac_p(i_wt)%qcl(i,k) - wtrac_e(i_wt)%qcl(i,k) ) )

    ! Created temporary updated field
    qclek_temp_wtrac(i,i_wt) = wtrac_e(i_wt)%qcl(i,k) +                        &
                                  tmp_dqclek_wtrac*timestep

    ! Update dqclbydt field
    wtrac_e(i_wt)%dqclbydt(i,k) = wtrac_e(i_wt)%dqclbydt(i,k) +                &
                                   tmp_dqclek_wtrac

    ! Water tracer qcl increment at level k+1
    wtrac_e(i_wt)%dqclbydt(i,k+1) = wtrac_e(i_wt)%dqclbydt(i,k+1)              &
                + flxbydpk(i) * (1.0-a_smth)                                   &
                * delpk(i)/delpkp1(i)*tmp_fd_dqclek_wtrac

    ! -------------------------------------------------------------------
    ! QCF calculations
    ! -------------------------------------------------------------------

    ! Smoothed forced detrainment term for water tracer qcf
    tmp_fd_dqcfek_wtrac = deltak(i) * (1.0-amdetk(i))                          &
                * ( wtrac_p(i_wt)%qcf(i,k) - wtrac_e(i_wt)%qcf(i,k) )

    ! Water tracer qcf increment at level k
    tmp_dqcfek_wtrac = flxbydpk(i)                                             &
                ! Compensating subsidence
                * ( (1.0-amdetk(i))*(1.0-deltak(i))*(1.0+ekp14(i))             &
                * ( wtrac_e(i_wt)%qcf(i,k+1) - wtrac_e(i_wt)%qcf(i,k) )        &
                ! Smoothed forced detrainment
                + a_smth*tmp_fd_dqcfek_wtrac                                   &
                ! Mixing detrainment
                + amdetk(i)                                                    &
                * ( wtrac_p(i_wt)%qcf(i,k) - wtrac_e(i_wt)%qcf(i,k) ) )

    ! Created temporary updated field
    qcfek_temp_wtrac(i,i_wt) = wtrac_e(i_wt)%qcf(i,k) +                        &
                                   tmp_dqcfek_wtrac*timestep

    ! Update dqcfbydt field
    wtrac_e(i_wt)%dqcfbydt(i,k) = wtrac_e(i_wt)%dqcfbydt(i,k) +                &
                                    tmp_dqcfek_wtrac

    ! Water tracer qcf increment at level k+1
    wtrac_e(i_wt)%dqcfbydt(i,k+1) = wtrac_e(i_wt)%dqcfbydt(i,k+1)              &
                + flxbydpk(i) * (1.0-a_smth)                                   &
                    * delpk(i)/delpkp1(i)*tmp_fd_dqcfek_wtrac

  end do ! m
end do ! i_wt

! Ratio of water tracer to normal water for phase change
! It is not clear what input values to use here.  Currently using environment
! values updated for subsidence and detrainment but before the increment is
! partitioned into liquid and ice phases (i.e. before the phase change).
! (Note that qclek and qclpk can be zero with dqclek_frz > 0, so can't
! simply use parcel or env values.)

do i_wt = 1, n_wtrac
  do m = 1, ni
    i = idx(m)
    ratio_qclk(i,i_wt) = wtrac_calc_ratio_fn(i_wt, qclek_temp_wtrac(i,i_wt),   &
                                             qclek_temp(i))
    ratio_qcfk(i,i_wt) = wtrac_calc_ratio_fn(i_wt, qcfek_temp_wtrac(i,i_wt),   &
                                             qcfek_temp(i))
  end do
end do


! Include changes due to phase changes (caused by partitioning the environment
! increment due to subsidence and detrainment between the liquid and ice
! phases in environ_6a)

do i_wt = 1, n_wtrac
  ! Freezing at level k
  call wtrac_move_phase(npnts, ni, idx, .false., dqclek_frz,                   &
                  ratio_qclk(:,i_wt), ratio_qcfk(:,i_wt),                      &
                  wtrac_e(i_wt)%dqclbydt(:,k), wtrac_e(i_wt)%dqcfbydt(:,k),    &
                  qchange_wtrac)

  ! Melting at level k
  call wtrac_move_phase(npnts, ni, idx, .false., dqcfek_mlt,                   &
                  ratio_qcfk(:,i_wt), ratio_qclk(:,i_wt),                      &
                  wtrac_e(i_wt)%dqcfbydt(:,k), wtrac_e(i_wt)%dqclbydt(:,k),    &
                  qchange_wtrac)

  ! Freezing at level k+1
  call wtrac_move_phase(npnts, ni, idx, .false., dqclekp1_frz,                 &
                  ratio_qclk(:,i_wt), ratio_qcfk(:,i_wt),                      &
                  wtrac_e(i_wt)%dqclbydt(:,k+1),                               &
                  wtrac_e(i_wt)%dqcfbydt(:,k+1), qchange_wtrac)

  ! Melting at level k+1
  call wtrac_move_phase(npnts, ni, idx, .false., dqcfekp1_mlt,                 &
                  ratio_qcfk(:,i_wt), ratio_qclk(:,i_wt),                      &
                  wtrac_e(i_wt)%dqcfbydt(:,k+1),                               &
                  wtrac_e(i_wt)%dqclbydt(:,k+1), qchange_wtrac)
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine environ_wtrac

end module environ_wtrac_mod

