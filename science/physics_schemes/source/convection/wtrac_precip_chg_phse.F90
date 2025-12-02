! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module wtrac_precip_chg_phse_mod

use um_types, only: real_umphys

implicit none

! Description:
!   Update water tracers for precipitation changing phase in convection
!   scheme (as calculated in crs_frzl for downdraughts and chg_phse for
!   environment or below cloud base).
!
! Method:
!   Calculates the amount of water changing phase and then updates the water
!   tracers.
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private ::                                        &
                           ModuleName = 'WTRAC_PRECIP_CHG_PHSE_MOD'
contains

! Subroutine Interface:
subroutine wtrac_precip_chg_phse(npnts, n_wtrac, rain, snow, wtrac_conv_old,   &
                                 wtrac_dd2, wtrac_ev)

use wtrac_conv_mod,          only: conv_dd2_wtrac_type, conv_ev_wtrac_type
use wtrac_conv_store_mod,    only: conv_old_wtrac_type
use wtrac_all_phase_chg_mod, only: wtrac_all_phase_chg

use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguments
integer, intent(in) :: npnts    ! Vector length
integer, intent(in) :: n_wtrac  ! No. of water tracers

real(kind=real_umphys), intent(in) :: rain(npnts)
                                ! Amount of rain descending from k-1
                                ! to k-2 after phase change (kg/m**2/s)
real(kind=real_umphys), intent(in) :: snow(npnts)
                                ! Amount of snow descending from k-1
                                ! to k-2 after phase change (kg/m**2/s)

type(conv_old_wtrac_type), intent(in out) :: wtrac_conv_old
                                ! Water values before phase change

! If called for downdraught calculation (i.e. after crs_frl) use:
type(conv_dd2_wtrac_type), optional, intent(in out) :: wtrac_dd2(n_wtrac)
                                ! Structure containing 2nd compression
                                ! water tracer arrays used in
                                ! downdraught calculations
! Or if called for environment calculation (i.e. after chg_phse) use:
type(conv_ev_wtrac_type), optional, intent(in out) :: wtrac_ev(n_wtrac)
                                ! Structure containing water tracer fields
                                ! used in phase change and evap of precip
                                ! in environment or below cloud base


! Local variables
integer :: i, i_wt       ! Loop counters

! Working arrays
real(kind=real_umphys) :: rain_wtrac(npnts,n_wtrac)  ! Water tracer rainfall
                                                     ! (kg/m**2/s)
real(kind=real_umphys) :: snow_wtrac(npnts,n_wtrac)  ! Water tracer snowfall
                                                     ! (kg/m**2/s)

logical :: l_downd   ! True if downdraught calculations

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='WTRAC_PRECIP_CHG_PHSE'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

if (present(wtrac_dd2)) then
  ! Downdraught
  l_downd = .true.

  ! Set working arrays
  do i_wt = 1, n_wtrac
    do i = 1, npnts
      ! Downdraught
      rain_wtrac(i,i_wt) = wtrac_dd2(i_wt)%rain(i)
      snow_wtrac(i,i_wt) = wtrac_dd2(i_wt)%snow(i)
    end do
  end do

else
  ! Environment
  l_downd = .false.

  ! Set working arrays
  do i_wt = 1, n_wtrac
    do i = 1, npnts
      ! Environment
      rain_wtrac(i,i_wt) = wtrac_ev(i_wt)%rain_env(i)
      snow_wtrac(i,i_wt) = wtrac_ev(i_wt)%snow_env(i)
    end do
  end do

end if   ! l_downd


!-----------------------------------------------------------------------
! Update water tracers for phase change
!-----------------------------------------------------------------------

call wtrac_all_phase_chg(npnts, n_wtrac,                                       &
                        wtrac_conv_old%rain, wtrac_conv_old%snow,              &
                        wtrac_conv_old%precip_frez,                            &
                        rain, snow, 'rai', 'sno', 'two_way',                   &
                        rain_wtrac, snow_wtrac)

! Update fields in water tracer structure
if (l_downd) then
  ! Downdraught
  do i_wt = 1, n_wtrac
    do i = 1, npnts
      wtrac_dd2(i_wt)%rain(i) = rain_wtrac(i,i_wt)
      wtrac_dd2(i_wt)%snow(i) = snow_wtrac(i,i_wt)
    end do
  end do

else
  ! Environment
  do i_wt = 1, n_wtrac
    do i = 1, npnts
      wtrac_ev(i_wt)%rain_env(i) = rain_wtrac(i,i_wt)
      wtrac_ev(i_wt)%snow_env(i) = snow_wtrac(i,i_wt)
    end do
  end do

end if  ! l_downd

! Store current water values
do i = 1, npnts
  wtrac_conv_old%rain(i) = rain(i)
  wtrac_conv_old%snow(i) = snow(i)
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine wtrac_precip_chg_phse

end module wtrac_precip_chg_phse_mod
