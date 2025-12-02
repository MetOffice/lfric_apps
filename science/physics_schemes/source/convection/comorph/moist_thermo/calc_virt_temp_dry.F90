! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_virt_temp_dry_mod

implicit none

contains


! Subroutine to compute dry virtual temperature
! = p / ( R_dry rho_dry )
! = T ( 1 + R_vap/R_dry q_vap )

! Version for compressed arrays
subroutine calc_virt_temp_dry( n_points,                                       &
                               temperature, q_vap,                             &
                               virt_temp_dry )

use comorph_constants_mod, only: real_cvprec, R_dry, R_vap

implicit none

! Number of points
integer, intent(in) :: n_points

! Temperature
real(kind=real_cvprec), intent(in) :: temperature(n_points)

! Water vapour mixing ratio
real(kind=real_cvprec), intent(in) :: q_vap(n_points)

! Output dry virtual temperature
! = p / ( R_dry rho_dry )
real(kind=real_cvprec), intent(out) :: virt_temp_dry(n_points)

! Ratio of gas constants R_vap/R_dry
real(kind=real_cvprec) :: R_vap_on_R_dry

real(kind=real_cvprec), parameter :: one = 1.0_real_cvprec

! Loop counter
integer :: i

R_vap_on_R_dry = R_vap / R_dry

do i = 1, n_points
  virt_temp_dry(i) = temperature(i)                                            &
                   * ( one + R_vap_on_R_dry * q_vap(i) )
end do

return
end subroutine calc_virt_temp_dry



! Version for full 3-D arrays
subroutine calc_virt_temp_dry_3d( lb_T, ub_T, temperature,                     &
                                  lb_q, ub_q, q_vap,                           &
                                  virt_temp_dry )

use comorph_constants_mod, only: real_hmprec, R_dry, R_vap,                    &
                     nx_full, ny_full, k_bot_conv, k_top_conv

implicit none

! Note: the array arguments below may have halos, and in here
! we don't care about the halos.
! We need to pass in the lower bound of each array to use in
! the declarations, to ensure the do loops pick out the same
! indices of the arrays regardless of whether or not they
! have halos.  The lower and upper bounds are passed in through the
! argument list in the lb_* and ub_* integer arrays.  Those storing the
! bounds for 3D arrays must have 3 elements; one for each dimension
! of the array.

! Temperature
integer, intent(in) :: lb_T(3), ub_T(3)
real(kind=real_hmprec), intent(in) :: temperature                              &
                         ( lb_T(1):ub_T(1), lb_T(2):ub_T(2), lb_T(3):ub_T(3) )

! Water vapour mixing ratio
integer, intent(in) :: lb_q(3), ub_q(3)
real(kind=real_hmprec), intent(in) :: q_vap                                    &
                         ( lb_q(1):ub_q(1), lb_q(2):ub_q(2), lb_q(3):ub_q(3) )

! Output dry virtual temperature
! = p / ( R_dry rho_dry )
real(kind=real_hmprec), intent(out) :: virt_temp_dry                           &
                     ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Ratio of gas constants R_vap/R_dry
real(kind=real_hmprec) :: R_vap_on_R_dry

real(kind=real_hmprec), parameter :: one = 1.0_real_hmprec

! Loop counters
integer :: i, j, k

R_vap_on_R_dry = real(R_vap,real_hmprec)                                       &
               / real(R_dry,real_hmprec)

!$OMP PARALLEL do DEFAULT(none) SCHEDULE(STATIC) private( i, j, k )            &
!$OMP SHARED( nx_full, ny_full, k_bot_conv, k_top_conv,                        &
!$OMP         virt_temp_dry, temperature, q_vap, R_vap_on_R_dry )
do k = k_bot_conv, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full
      virt_temp_dry(i,j,k) = temperature(i,j,k)                                &
                     * ( one + R_vap_on_R_dry * q_vap(i,j,k) )
    end do
  end do
end do
!$OMP end PARALLEL do

return
end subroutine calc_virt_temp_dry_3d


end module calc_virt_temp_dry_mod
