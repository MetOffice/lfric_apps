! -- begin wtrac_move_phase.h --
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers


integer, intent(in) :: npnts           ! Number of points
integer, intent(in) :: ni              ! Number of work points
integer, intent(in) :: idx(ni)         ! Indices of work points

logical, intent(in) :: l_limit_chg     ! If T, limit the amount of
                                       ! water tracer changing phase so that
                                       ! it does exceed the amount of
                                       ! the source

real(kind=field_kind), intent(in) :: qchange(npnts)
                                       ! Amount of water changing phase (kg/kg)

real(kind=field_kind), intent(in) :: qratio1_wtrac(npnts)
                                       ! Ratio of water tracer to normal
                                       ! water for phase type 1

real(kind=field_kind), intent(in) :: qratio2_wtrac(npnts)
                                       ! Ratio of water tracer to normal
                                       ! water for phase type 2

real(kind=field_kind), intent(in out) :: q1_wtrac(npnts)
                                       ! Water tracer for phase type 1 (kg/kg)
real(kind=field_kind), intent(in out) :: q2_wtrac(npnts)
                                       ! Water tracer for phase type 2 (kg/kg)

real(kind=field_kind),intent(out) :: qchange_wtrac(npnts)
                                       ! Change in water tracer (phase type 2)
                                       ! due to phase change (kg/kg)


! Local variables
integer :: i, m             ! Loop counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

qchange_wtrac(:) = 0.0

! Update water tracer for phase change

if (l_limit_chg .and. l_wtrac_limit_chg) then
  ! Ensure that the amount of water tracer changing phase is not greater than
  ! the source amount.
  do m=1, ni

    i = idx(m)

    if (qchange(i) > 0.0) then   ! q1 is the source
      qchange_wtrac(i) = min(qratio1_wtrac(i) * qchange(i),q1_wtrac(i))
    else                         ! q2 is the source
      qchange_wtrac(i) = max(qratio2_wtrac(i) * qchange(i),-q2_wtrac(i))
    end if

    q1_wtrac(i) = q1_wtrac(i) - qchange_wtrac(i)
    q2_wtrac(i) = q2_wtrac(i) + qchange_wtrac(i)

  end do

else  ! There are some occasions when there should be no limit enforced here:
  ! a) when called from wtrac_pc2_bl which does its own limit enforcement and
  ! b) when called from environ_wtrac as q_wtrac is actually dqbydt_wtrac

  do m=1, ni

    i = idx(m)

    if (qchange(i) > 0.0) then   ! q1 is the source
      qchange_wtrac(i) = qratio1_wtrac(i) * qchange(i)
    else                         ! q2 is the source
      qchange_wtrac(i) = qratio2_wtrac(i) * qchange(i)
    end if

    q1_wtrac(i) = q1_wtrac(i) - qchange_wtrac(i)
    q2_wtrac(i) = q2_wtrac(i) + qchange_wtrac(i)

  end do

end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
