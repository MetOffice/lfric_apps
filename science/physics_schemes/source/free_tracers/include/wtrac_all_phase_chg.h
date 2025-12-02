! -- begin wtrac_all_phase_chg.h --
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers
!

integer, intent(in) :: npnts   ! No. of points
integer, intent(in) :: n_wtrac ! No. of water tracers

real(kind=field_kind), intent(in) :: q1_old(npnts)
                               ! Water field 1 before phase change
real(kind=field_kind), intent(in) :: q2_old(npnts)
                               ! Water field 2 before phase change
real(kind=field_kind), intent(in) :: q_change(npnts)
                               ! Amount of water changing phase
                               ! (+ve for water moving type 1 -> 2)

! The following two fields are currently unused, but will be used in the future
real(kind=field_kind), intent(in) :: q1(npnts)
                               ! Water field 1 after phase change
real(kind=field_kind), intent(in) :: q2(npnts)
                               ! Water field 1 after phase change

character(len=3), intent(in) :: q1_type ! Water field 1 type ('vap','liq' etc.)
character(len=3), intent(in) :: q2_type ! Water field 2 type

character(len=7), intent(in) :: direction ! Direction of phase change
                                          ! ('one_way' or 'two_way')

real(kind=field_kind), intent(in out) :: q1_wtrac(npnts,n_wtrac)
                               ! Water tracer field 1
real(kind=field_kind), intent(in out) :: q2_wtrac(npnts,n_wtrac)
                               ! Water tracer field 2

! The following logical is currently unused but will be used in the future
logical, intent(in),optional :: l_no_fract     ! No fractionation for isotopes

logical, intent(in),optional :: l_limit_chg    ! If T, limit phase change to
                                               ! ensure the source water
                                               ! tracer value stays >=0

real(kind=field_kind), intent(out),optional :: qchange_wtrac(npnts,n_wtrac)
                                               ! Change in water tracer

! Local variables
integer :: i, i_wt, m     ! Loop counters
integer :: ni             ! No. of work points

integer :: idx(npnts)     ! Index of work points

real(kind=field_kind) :: q1_ratio(npnts)
                          ! Ratio of water tracer q1 to 'normal' q1
real(kind=field_kind) :: q2_ratio(npnts)
                          ! Ratio of water tracer q2 to 'normal' q2
real(kind=field_kind) :: qchange_wtrac1(npnts,n_wtrac)
                          ! Change in water tracer

logical :: l_limit_chg1   ! Local version of l_limit_chg

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

! End of header
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Set up l_limit_chg1.  When T this ensures that the amount of water tracer
! changing phase is not greater than the source amount, which is the default
! setting.  However, there are some occasions when there should be no limit
! enforced:
! a) when called from wtrac_pc2_bl which does its own limit enforcement and
! Note, currently this check is not enforced in wtrac_move_phase as more
! work needs to be done to see if this check is necessary or not.

if (present(l_limit_chg)) then
  l_limit_chg1 = l_limit_chg
else
  l_limit_chg1 = .true.
end if

! Set up working points index
ni = 0
do i = 1, npnts

  if (abs(q_change(i)) > 0.0 ) then
    ni = ni + 1
    idx(ni) = i
  end if

  ! Initialise fields
  q1_ratio(i) = 0.0
  q2_ratio(i) = 0.0
  do i_wt = 1, n_wtrac
    qchange_wtrac1(i,i_wt) = 0.0
  end do

end do

if (ni > 0) then
  do i_wt = 1, n_wtrac

    ! Calculate ratios of water tracer to 'normal' water prior to phase change
    do m = 1, ni
      i = idx(m)
      q1_ratio(i) = wtrac_calc_ratio_fn(i_wt, q1_wtrac(i,i_wt), q1_old(i))

      if (direction == 'two_way') then
        ! Only needed if phase change is in both directions
        q2_ratio(i) = wtrac_calc_ratio_fn(i_wt, q2_wtrac(i,i_wt), q2_old(i))
      end if

    end do

    ! Update water tracers for phase change
    call wtrac_move_phase(npnts, ni, idx, l_limit_chg1,                        &
                     q_change, q1_ratio, q2_ratio,                             &
                     q1_wtrac(:,i_wt), q2_wtrac(:,i_wt),                       &
                     qchange_wtrac=qchange_wtrac1(:,i_wt))

  end do   ! i_wt
end if  ! ni > 0

if (present(qchange_wtrac)) then
  qchange_wtrac(:,:) = qchange_wtrac1(:,:)
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
!-- end wtrac_all_phase_chg.h --
