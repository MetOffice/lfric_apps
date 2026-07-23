! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module test_unstable_mod

implicit none

contains

! Subroutine to test whether a model-level might possibly be unstable
! enough to trigger convection
subroutine test_unstable( virt_temp_1, virt_temp_2,                            &
                          lb_z, ub_z, height_1, height_2,                      &
                          lb_p, ub_p, pressure_1, pressure_2,                  &
                          lb_v, ub_v, q_vap,  lb_l, ub_l, q_cl,                &
                          lb_r, ub_r, q_rain, lb_f, ub_f, q_cf,                &
                          lb_s, ub_s, q_snow, lb_g, ub_g,q_graup,              &
                          l_init_poss )

use comorph_constants_mod, only: real_hmprec, nx_full, ny_full,                &
                                 gravity
use dry_adiabat_mod, only: dry_adiabat_2d

implicit none

! Virtual temperature at current and next level
real(kind=real_hmprec), intent(in) :: virt_temp_1 ( nx_full, ny_full )
real(kind=real_hmprec), intent(in) :: virt_temp_2 ( nx_full, ny_full )

! Height of current and next level
integer, intent(in) :: lb_z(2), ub_z(2)
real(kind=real_hmprec), intent(in) :: height_1                                 &
                                      ( lb_z(1):ub_z(1), lb_z(2):ub_z(2) )
real(kind=real_hmprec), intent(in) :: height_2                                 &
                                      ( lb_z(1):ub_z(1), lb_z(2):ub_z(2) )

! Pressure of current and next level
integer, intent(in) :: lb_p(2), ub_p(2)
real(kind=real_hmprec), intent(in) :: pressure_1                               &
                                      ( lb_p(1):ub_p(1), lb_p(2):ub_p(2) )
real(kind=real_hmprec), intent(in) :: pressure_2                               &
                                      ( lb_p(1):ub_p(1), lb_p(2):ub_p(2) )

! Water species mixing ratios at current level
integer, intent(in) :: lb_v(2), ub_v(2)
real(kind=real_hmprec), intent(in) :: q_vap                                    &
                                      ( lb_v(1):ub_v(1), lb_v(2):ub_v(2) )
integer, intent(in) :: lb_l(2), ub_l(2)
real(kind=real_hmprec), intent(in) :: q_cl                                     &
                                      ( lb_l(1):ub_l(1), lb_l(2):ub_l(2) )
integer, intent(in) :: lb_r(2), ub_r(2)
real(kind=real_hmprec), intent(in) :: q_rain                                   &
                                      ( lb_r(1):ub_r(1), lb_r(2):ub_r(2) )
integer, intent(in) :: lb_f(2), ub_f(2)
real(kind=real_hmprec), intent(in) :: q_cf                                     &
                                      ( lb_f(1):ub_f(1), lb_f(2):ub_f(2) )
integer, intent(in) :: lb_s(2), ub_s(2)
real(kind=real_hmprec), intent(in) :: q_snow                                   &
                                      ( lb_s(1):ub_s(1), lb_s(2):ub_s(2) )
integer, intent(in) :: lb_g(2), ub_g(2)
real(kind=real_hmprec), intent(in) :: q_graup                                  &
                                      ( lb_g(1):ub_g(1), lb_g(2):ub_g(2) )

! Flag for whether convective initiation might be possible from current level
logical, intent(in out) :: l_init_poss ( nx_full, ny_full )

! Virtual temperature for ascent from current to next level
real(kind=real_hmprec) :: virt_temp_test ( nx_full, ny_full )

! Dry static stability
real(kind=real_hmprec) :: Nsq_dry ( nx_full, ny_full )

! Parameters converted to host-model precision
real(kind=real_hmprec), parameter :: zero_p = 0.0_real_hmprec
real(kind=real_hmprec) :: gravity_p


! Loop counters
integer :: i, j


! Convert constants to same precision as the full 2-D fields
gravity_p = real( gravity, real_hmprec )

! Calculate virtual temperature of level k ascended from current to next level
do j = 1, ny_full
  do i = 1, nx_full
    virt_temp_test(i,j) = virt_temp_1(i,j)
  end do
end do
call dry_adiabat_2d( lb_p, ub_p, pressure_1, pressure_2,                       &
                     lb_v, ub_v, q_vap,  lb_l, ub_l, q_cl,                     &
                     lb_r, ub_r, q_rain, lb_f, ub_f, q_cf,                     &
                     lb_s, ub_s, q_snow, lb_g, ub_g,q_graup,                   &
                     virt_temp_test )

! Compute dry static stability
do j = 1, ny_full
  do i = 1, nx_full
    Nsq_dry(i,j) = ( gravity_p / virt_temp_2(i,j) )                            &
                 * ( virt_temp_2(i,j) - virt_temp_test(i,j) )                  &
                 / ( height_2(i,j) - height_1(i,j) )
  end do
end do

! No gradient adjustment used in triggering, so the mask only needs
! to include points where the stratification is actually unstable
do j = 1, ny_full
  do i = 1, nx_full
    l_init_poss(i,j) = l_init_poss(i,j) .or. Nsq_dry(i,j) < zero_p
  end do
end do


return
end subroutine test_unstable

end module test_unstable_mod
