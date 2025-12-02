! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module correct_small_q_conv_mod

use um_types, only: real_umphys

implicit none
! Description:
!   Correct negative/very small water vapour amounts at end of convection
!
! Method:
!   Apply an artificial upwards flux from k-1 level to ensure water vapour
!   remains above a minimum value in the column.
!
! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

character(len=*), parameter, private :: ModuleName = 'CORRECT_SMALL_Q_CONV_MOD'

contains

! Subroutine Interface:
subroutine correct_small_q_conv(ni, npnts, nlev, n_wtrac, index1, timestep,    &
                                qmin, r2rho_th, dr_across_th, q,               &
                                field_to_check, conv_type, dqbydt, wtrac_e)

use wtrac_conv_mod,           only: l_wtrac_conv, conv_e_wtrac_type
use wtrac_calc_ratio_mod,     only: wtrac_calc_ratio_fn

use umPrintMgr, only: umPrint, ummessage, printstatus, prstatus_normal

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

! Subroutine arguements
integer, intent(in) :: ni        ! No. of working points
integer, intent(in) :: npnts     ! No. of points
integer, intent(in) :: nlev      ! No. of vertical levels
integer, intent(in) :: n_wtrac   ! No. of water tracers

integer, intent(in) :: index1(ni)  ! Working point indices

real(kind=real_umphys), intent(in) :: timestep   ! Timestep
real(kind=real_umphys), intent(in) :: qmin       ! Minimum value of q
real(kind=real_umphys), intent(in) :: r2rho_th(npnts,nlev)
                                             ! radius**2 density for
                                             ! theta lev (kg/m)
real(kind=real_umphys), intent(in) :: dr_across_th(npnts,nlev)
                                             ! Thickness of theta levels (m)

real(kind=real_umphys), intent(in) :: q(npnts,nlev)
                                             ! Mixing ratio (kg/kg)

character(len=5), intent(in)       :: field_to_check
                                             ! Which field is being tested?
                                             ! ('water' or 'wtrac')

character(len=*), intent(in)       :: conv_type
                                             ! Type of convection
                                             ! ('deep', 'shal' or 'mid')

real(kind=real_umphys), intent(in out) :: dqbydt(npnts,nlev)
                                             ! Increments to q due to
                                             ! convection (kg/kg/s)

type(conv_e_wtrac_type), optional, intent(in out) :: wtrac_e(n_wtrac)
                                               ! Structure containing water
                                               ! tracer fields

! Local variables
integer :: i, j, k, i_wt                        ! Loop counters

real(kind=real_umphys) :: ratio_wt            ! Water tracer ratio
real(kind=real_umphys) :: qminincolumn(ni)    ! Min value for q in column(kg/kg)
real(kind=real_umphys) :: temp1(ni)           ! Temporary array
real(kind=real_umphys) :: delta_dqbydt(ni,nlev)  ! Store change in dqbydt

logical                ::  flag_correct      ! True if correction has been made

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CORRECT_SMALL_Q_CONV'

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Find minimum value in column
do j = 1,ni
  i=index1(j)
  qminincolumn(j) = q(i,nlev)
end do
do k = 1,nlev-1
  do j = 1,ni
    i=index1(j)
    if (q(i,k)  <   qminincolumn(j)) then
      qminincolumn(j) = q(i,k)
    end if
  end do
end do

! Ensure Q does not go below global allowed minimum (QMIN)

do j = 1,ni
  qminincolumn(j)=max(qmin,qminincolumn(j))
end do

! Apply an artificial upwards flux from k-1 level to ensure Q
! remains above minimum value in the column.

flag_correct = .false.

do k = nlev,2,-1                   ! Note, working downwards
  do j = 1,ni
    i=index1(j)
    delta_dqbydt(j,k) = 0.0
    if (dqbydt(i,k) /= 0.0) then
      temp1(j)=q(i,k) + dqbydt(i,k) * timestep
      if (temp1(j)  <   qminincolumn(j)) then

        ! Store change in dqbydt for water tracer use
        delta_dqbydt(j,k) = (qminincolumn(j) - q(i,k)) / timestep              &
                                  - dqbydt(i,k)

        dqbydt(i,k-1) = dqbydt(i,k-1) -                                        &
            ((qminincolumn(j) - q(i,k)) / timestep-dqbydt(i,k))                &
             * (r2rho_th(i,k)*dr_across_th(i,k))                               &
             / (r2rho_th(i,k-1)*dr_across_th(i,k-1))

        dqbydt(i,k) = (qminincolumn(j) - q(i,k)) / timestep

        flag_correct = .true.
      end if
    end if
  end do ! ni loop
end do  ! nlev

if (l_wtrac_conv .and. field_to_check == 'water' .and. flag_correct) then
  ! Update water tracers for changes to water by applying the same
  ! artificial vertical flux
  do k = nlev,2,-1                     ! Note, working downwards
    do j = 1,ni
      i=index1(j)
      if (abs(delta_dqbydt(j,k)) > 0.0 ) then
        do i_wt = 1, n_wtrac

          ratio_wt = wtrac_calc_ratio_fn(i_wt,                                 &
                                         wtrac_e(i_wt)%q(i,k-1), q(i,k-1))

          wtrac_e(i_wt)%dqbydt(i,k-1) =  wtrac_e(i_wt)%dqbydt(i,k-1)           &
                       - (ratio_wt * delta_dqbydt(j,k))                        &
                         * (r2rho_th(i,k)*dr_across_th(i,k))                   &
                         / (r2rho_th(i,k-1)*dr_across_th(i,k-1))

          wtrac_e(i_wt)%dqbydt(i,k) = wtrac_e(i_wt)%dqbydt(i,k)                &
                                     +  (ratio_wt * delta_dqbydt(j,k))

        end do  ! n_wtrac
      end if
    end do ! ni loop
  end do  ! nlev
end if

! Check for negative q at k=1
if (field_to_check == 'water') then

  k=1
  do j = 1,ni
    i=index1(j)
    temp1(j)=q(i,k) + dqbydt(i,k) * timestep
    if (temp1(j)  < qminincolumn(j) .and.                                      &
                printstatus >= prstatus_normal ) then
      write(umMessage,'(a,a27,i0,a10,g26.18,a7,g26.18)')                       &
      conv_type, ' convection, negative q, i:', i,' ,q after ',temp1(j),       &
      ' dq/dt ',dqbydt(i,k)
      call umPrint(umMessage,src='correct_small_q_conv')
    end if
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine correct_small_q_conv

end module correct_small_q_conv_mod
