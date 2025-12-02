! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module cor_engy_wtrac_mod

use um_types, only: real_umphys

implicit none

! Description:
!   Update water tracers for conservation correction applied in
!   convection scheme (cor_engy_mod-6a)
! Note, water tracers only work with the 6A convection scheme.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
!   This file belongs in section: Water_Tracers
!
! Code description:
!   Language: Fortran 90.
!   This code is written to UM programming standards version 8.3.

character(len=*), parameter, private :: ModuleName='COR_ENGY_WTRAC_MOD'

contains

subroutine cor_engy_wtrac(npnts, nconv, nlev, n_wtrac, index1, timestep, q,    &
                          dqbydt, dmvscale, wtrac_e)

use wtrac_conv_mod,           only: conv_e_wtrac_type
use wtrac_calc_ratio_mod,     only: wtrac_calc_ratio_fn

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none
!
! ------------------------------------------------------------------------------
! Subroutine arguments

integer,intent(in) :: npnts         ! Full vector length
integer,intent(in) :: nconv         ! Number of convecting points
integer,intent(in) :: nlev          ! Number of model levels for calculations
integer,intent(in) :: n_wtrac       ! Number of water tracers

integer,intent(in) :: index1(npnts) ! index of points with convection

real(kind=real_umphys), intent(in) :: timestep
                                            ! Timestep
real(kind=real_umphys),intent(in) :: q(npnts,nlev)
                                            ! Specific humidity on theta levels
                                            ! in kg/kg
real(kind=real_umphys),intent(in) :: dqbydt(npnts,nlev)
                                            ! Increment to specific
                                            ! water vapour (kg/kg/s)
real(kind=real_umphys),intent(in) :: dmvscale(nconv)
                                            ! scaling factor used to correct
                                            !the water increments
type(conv_e_wtrac_type), intent(in out) :: wtrac_e(n_wtrac)
                                            ! Structure containing water
                                            ! tracer fields

!---------------------------------------------------------------------
! Local variables
!---------------------------------------------------------------------

integer :: i,j,k,i_wt                  ! loop counters

real(kind=real_umphys) :: q_new        ! Updated q
real(kind=real_umphys) :: q_wtrac_new  ! Updated q_wtrac
real(kind=real_umphys) :: ratio_wt     ! Ratio of water tracer to water

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='COR_ENGY_WTRAC'

!----------------------------------------------------------------------
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

do i_wt = 1, n_wtrac
  do k = 1, nlev
    do j = 1, nconv
      i = index1(j)

      ! Temporarily update q fields to calculate ratio
      q_new       = q(i,k) + dqbydt(i,k) * timestep
      q_wtrac_new = wtrac_e(i_wt)%q(i,k) + wtrac_e(i_wt)%dqbydt(i,k) * timestep
      ratio_wt    = wtrac_calc_ratio_fn(i_wt, q_wtrac_new, q_new)

      ! Add correction to water tracer rate of change field
      wtrac_e(i_wt)%dqbydt(i,k) = wtrac_e(i_wt)%dqbydt(i,k)  +                 &
                                  ratio_wt * (dmvscale(j)-1.0)*dqbydt(i,k)
    end do  ! nconv
  end do ! nlev
end do ! n_wtrac

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine cor_engy_wtrac

end module cor_engy_wtrac_mod
