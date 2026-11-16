! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module interp_virt_temp_mod

implicit none

contains


! Subroutine to interpolate the virtual temperature profile from
! full-levels to half-levels.  Attempts to account for kinks in the
! profile, e.g. at inversions.  This calculation is done on the full
! 3D fields, before the compression onto convecting points, because
! the kink calculations are vertically non-local.
subroutine interp_virt_temp( lb_hf, ub_hf, height_full,                        &
                             lb_hh, ub_hh, height_half,                        &
                             lb_pf, ub_pf, pressure_full,                      &
                             lb_ph, ub_ph, pressure_half,                      &
                             virt_temp_full, virt_temp_half )

use comorph_constants_mod, only: real_hmprec,                                  &
                                 nx_full, ny_full, k_bot_conv, k_top_conv,     &
                                 R_dry, cp_dry

implicit none

! Heights of full-levels
integer, intent(in) :: lb_hf(3)
integer, intent(in) :: ub_hf(3)
real(kind=real_hmprec), intent(in) :: height_full                              &
             ( lb_hf(1):ub_hf(1), lb_hf(2):ub_hf(2), lb_hf(3):ub_hf(3) )

! Heights of half-levels
integer, intent(in) :: lb_hh(3)
integer, intent(in) :: ub_hh(3)
real(kind=real_hmprec), intent(in) :: height_half                              &
             ( lb_hh(1):ub_hh(1), lb_hh(2):ub_hh(2), lb_hh(3):ub_hh(3) )

! Pressures of full-levels
integer, intent(in) :: lb_pf(3)
integer, intent(in) :: ub_pf(3)
real(kind=real_hmprec), intent(in) :: pressure_full                            &
             ( lb_pf(1):ub_pf(1), lb_pf(2):ub_pf(2), lb_pf(3):ub_pf(3) )

! Pressures of half-levels
integer, intent(in) :: lb_ph(3)
integer, intent(in) :: ub_ph(3)
real(kind=real_hmprec), intent(in) :: pressure_half                            &
             ( lb_ph(1):ub_ph(1), lb_ph(2):ub_ph(2), lb_ph(3):ub_ph(3) )

! Virtual temperature on full-levels
real(kind=real_hmprec), intent(in) :: virt_temp_full                           &
             ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Virtual temperature on half-levels
real(kind=real_hmprec), intent(out) :: virt_temp_half                          &
             ( nx_full, ny_full, k_bot_conv:k_top_conv+1 )

! Virtual potential temperature used in interpolation
real(kind=real_hmprec) :: theta_v_full                                         &
             ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Vertical gradient of virtual temperature (on half-levels)
real(kind=real_hmprec) :: dtvdz

! R_dry / cp_dry
real(kind=real_hmprec) :: Rd_over_cpd

! Loop counters
integer :: i, j, k


! Store physics constants in host-model precision
Rd_over_cpd = real( R_dry / cp_dry, real_hmprec )

!$OMP PARALLEL DEFAULT(NONE)                                                   &
!$OMP PRIVATE( i, j, k, dtvdz )                                                &
!$OMP SHARED( nx_full, ny_full, k_bot_conv, k_top_conv,                        &
!$OMP         virt_temp_full, virt_temp_half, height_full, height_half,        &
!$OMP         pressure_full, pressure_half, theta_v_full, Rd_over_cpd )

  ! Dry adiabat conserves potential temperature
  ! Scale virtual temperature by pressure^(-R/cp)
!$OMP DO SCHEDULE(STATIC)
do k = k_bot_conv, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full
      theta_v_full(i,j,k) = virt_temp_full(i,j,k)                              &
                          * ( pressure_full(i,j,k)**(-Rd_over_cpd) )
    end do
  end do
end do
!$OMP END DO

! Set values at top and bottom by just copying the values from the nearest
! full level, since we can't interpolate (these values will only be used
! if convection tries to go beyond the top / bottom, in which case something
! has already gone wrong!)
!$OMP DO SCHEDULE(STATIC)
do j = 1, ny_full
  do i = 1, nx_full
    virt_temp_half(i,j,k_bot_conv)   = theta_v_full(i,j,k_bot_conv)
    virt_temp_half(i,j,k_top_conv+1) = theta_v_full(i,j,k_top_conv)
  end do
end do
!$OMP END DO NOWAIT

! Set all other half-levels using linear interpolation initially
!$OMP DO SCHEDULE(STATIC)
do k = k_bot_conv+1, k_top_conv
  do j = 1, ny_full
    do i = 1, nx_full
      ! Compute vertical gradients, centred on half-levels
      dtvdz = ( theta_v_full(i,j,k) - theta_v_full(i,j,k-1) )                  &
            / (  height_full(i,j,k) -  height_full(i,j,k-1) )
      ! Linear interpolation onto half-levels
      virt_temp_half(i,j,k) = theta_v_full(i,j,k-1)                            &
                            + ( height_half(i,j,k) - height_full(i,j,k-1) )    &
                              * dtvdz
    end do
  end do
end do
!$OMP END DO

! Dry adiabat conserves potential temperature
! Scale by pressure^(R/cp) to retrieve virtual temperature
!$OMP DO SCHEDULE(STATIC)
do k = k_bot_conv, k_top_conv+1
  do j = 1, ny_full
    do i = 1, nx_full
      virt_temp_half(i,j,k) = virt_temp_half(i,j,k)                            &
                            * ( pressure_half(i,j,k)**(Rd_over_cpd) )
    end do
  end do
end do
!$OMP END DO NOWAIT

!$OMP END PARALLEL


return
end subroutine interp_virt_temp

end module interp_virt_temp_mod
