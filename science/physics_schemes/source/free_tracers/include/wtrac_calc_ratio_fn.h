! -- begin wtrac_calc_ratio.h --
! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner:  Please refer to the UM file CodeOwners.txt
! This file belongs in section: Water_Tracers

integer, intent(in) :: i_wt   ! Water tracer number

real(kind=field_kind), intent(in) :: q_wtrac    ! Water tracer (kg/kg)
real(kind=field_kind), intent(in) :: q          ! Normal water (kg/kg)

logical, optional, intent(in)     :: l_rate_in  ! If True, then the q input
                                                ! is a rate of change
real(kind=field_kind) :: q_ratio                ! Output ratio
real(kind=field_kind) :: q_ratio_std            ! Standard ratio
real(kind=field_kind) :: min_q_ratio_use        ! Min q used in ratio calc

logical :: l_check_neg_ratio                    ! If True, do not allow
                                                ! negative ratios if
                                                ! l_wtrac_no_neg_ratio=T


if (present(l_rate_in)) then
  ! If the input is a rate of change, rather than concentration,
  ! then allow negative ratios.
  l_check_neg_ratio = .false.
else
  ! Default is that fields are concentrations and therefore negative
  ! ratios are set to the standard ratio if l_wtrac_no_neg_ratio = T
  l_check_neg_ratio = l_wtrac_no_neg_ratio
end if

q_ratio_std = real(wtrac_info(i_wt)%standard_ratio, kind=field_kind)
min_q_ratio_use = real(min_q_ratio, kind=field_kind)

if (abs(q) > min_q_ratio_use) then
  q_ratio = q_wtrac / q
else
  q_ratio = q_ratio_std
end if

if (l_check_neg_ratio) then
  if (q_ratio < 0.0) q_ratio = q_ratio_std
end if

! -- end wtrac_calc_ratio_fn.h --
