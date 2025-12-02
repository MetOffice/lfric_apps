! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module calc_wp_below_t_mod

! Purpose: calculates any water path in units of kg m-2

! This is similar in nature to the calculations performed
! for Section 30 (Climate diagnostics), the code of which can
! be found in vert_eng_massq (within energy_correction)
! However, the code in calc_wp is much cheaper than vert_eng_massq
! as it only performs one integral for one liquid species only.
! Therefore calc_wp is recommended for simple calculations involving
! one or two moist variables.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_calc_wp

use atm_fields_bounds_mod, only: tdims

! Dr Hook modules
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='CALC_WP_BELOW_T_MOD'

contains

subroutine calc_wp_below_t( wc, rhodz, t, t_thresh, wp_tot )

implicit none

!-------------------------------
! Subroutine arguments
!-------------------------------

real(kind=real_umphys), intent(in)::                                           &
                      wc(tdims%i_start : tdims%i_end,                          &
                         tdims%j_start : tdims%j_end,                          &
                                     1 : tdims%k_end)
! Any water mixing ratio (q, qcl, qcf, qrain, qgraup)- units kg/kg

real(kind=real_umphys), intent(in)::                                           &
                   rhodz(tdims%i_start : tdims%i_end,                          &
                         tdims%j_start : tdims%j_end,                          &
                                     1 : tdims%k_end)
! Air density * depth of model level (kg m-2)

real(kind=real_umphys), intent(in)::                                           &
                   t(tdims%i_start : tdims%i_end,                              &
                     tdims%j_start : tdims%j_end,                              &
                                 1 : tdims%k_end)
! Temperature (K)

real(kind=real_umphys), intent(in):: t_thresh
! Temperature threshold (K)

real(kind=real_umphys), intent(out) ::                                         &
                     wp_tot(tdims%i_start : tdims%i_end,                       &
                            tdims%j_start : tdims%j_end)

!-------------------
! Local variables
!-------------------

integer :: i, j, k ! Loop indices

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CALC_WP_BELOW_T'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL DEFAULT(none) private(i,j,k)                                    &
!$OMP SHARED(tdims,wp_tot,wc,t,rhodz,t_thresh)

!$OMP do SCHEDULE(STATIC)
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end
    wp_tot(i,j) = 0.0
  end do ! i
end do ! j
!$OMP end do

do k = 1, tdims%k_end
!$OMP do SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end
    do i = tdims%i_start, tdims%i_end

      if (wc(i,j,k) > 0.0 .and. t(i,j,k) <= t_thresh) then
         ! only include if water content is positive and temperature
         ! below threshold
        wp_tot(i,j) = wp_tot(i,j) + (rhodz(i,j,k) * wc(i,j,k) )
      end if

    end do ! i
  end do ! j
!$OMP end do
end do ! k
!$OMP end PARALLEL

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine calc_wp_below_t

end module calc_wp_below_t_mod

