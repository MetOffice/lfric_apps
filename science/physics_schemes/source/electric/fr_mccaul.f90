! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module fr_mccaul_mod

! Purpose: Calculates flash rate based on McCaul et al (2009), Weather and
! forecasting.

! Full paper reference:
! McCaul, E. W., Goodman, S. J., LaCasse, K. M., Cecil, D. J. 2009.
! Forecasting Lightning Threat Using Cloud-Resolving Model Simulations.
! Weather and Forecasting, Volume 24, pp 709-729.
! doi: http://dx.doi.org/10.1175/2008WAF2222152.1

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Electric

use atm_fields_bounds_mod,  only: tdims, wdims_s
use electric_constants_mod, only: minus15, mccaul_r1, mccaul_r2, kg_to_g,      &
                                  fivemin2sec, zero
use electric_inputs_mod,    only: k1, k2

! Dr Hook modules
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='FR_MCCAUL_MOD'

contains

subroutine fr_mccaul(storm_field, t, w, qgraup, gwp, tiwp, flash, fr1_mc,      &
                     fr2_mc )

implicit none

logical, intent(in) :: storm_field( tdims%i_start : tdims%i_end,               &
                                    tdims%j_start : tdims%j_end )
! Defines where storms exist in model

real(kind=real_umphys), intent(in) :: w(                                       &
                                    wdims_s%i_start : wdims_s%i_end,           &
                                    wdims_s%j_start : wdims_s%j_end,           &
                                    wdims_s%k_start : wdims_s%k_end )
! Vertical velocity (m/s)

real(kind=real_umphys), intent(in) :: t(                                       &
                                    tdims%i_start : tdims%i_end,               &
                                    tdims%j_start : tdims%j_end,               &
                                                1 : tdims%k_end )
! Temperature (K)
real(kind=real_umphys), intent(in) :: qgraup(                                  &
                                    tdims%i_start : tdims%i_end,               &
                                    tdims%j_start : tdims%j_end,               &
                                                1 : tdims%k_end )
! Graupel mixing ratio (kg/kg)

real(kind=real_umphys), intent(in) :: gwp(                                     &
                                    tdims%i_start : tdims%i_end,               &
                                    tdims%j_start : tdims%j_end)
real(kind=real_umphys), intent(in) :: tiwp(                                    &
                                    tdims%i_start : tdims%i_end,               &
                                    tdims%j_start : tdims%j_end)
! Graupel and total ice water paths (kg m-2)

real(kind=real_umphys), intent(in out) :: flash(                               &
                              tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end )

real(kind=real_umphys), intent(in out) :: fr1_mc(                              &
                                tdims%i_start : tdims%i_end,                   &
                              tdims%j_start : tdims%j_end )

real(kind=real_umphys), intent(in out) :: fr2_mc(                              &
                                tdims%i_start : tdims%i_end,                   &
                              tdims%j_start : tdims%j_end )

! Flash rates (flash 1 and flash2 for diagnostics). Units: s-1


! local variables

integer :: m15l, upper_limit, lower_limit
integer :: i, j, k
real(kind=real_umphys) :: flash1
real(kind=real_umphys) :: flash2

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='FR_MCCAUL'

!==================================================================
! Start the subroutine
!==================================================================

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Define upper and lower limits for minus 15 Celius calculation
upper_limit = tdims%k_end   - 2
lower_limit = tdims%k_start + 1
! initialise m15l variable
m15l=-999

!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none) private(i,j,flash1,           &
!$OMP flash2,k,m15l)                                                           &
!$OMP SHARED(tdims,storm_field,t,upper_limit,lower_limit,w,qgraup,k1,          &
!$OMP k2,gwp,tiwp,fr1_mc,fr2_mc,flash)
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end

    ! Initialise flash1 and flash2

    flash1 = zero
    flash2 = zero

    ! Only perform calculation if this point is defined as a storm point
    if (storm_field(i,j)) then

      ! First find the minus 15 level
      do k = tdims%k_end, 1, -1

        if (t(i,j,k) > minus15) then

          m15l = k+1 ! define minus 15 level
          exit ! leave the loop early to avoid overwriting

        end if

      end do ! k

      ! Determine flash rates

      if (m15l > upper_limit .or. m15l < lower_limit) then
        ! Set flash1 to be zero if you cannot define the
        ! -15C level in the model.

        flash1 = zero

      else ! m15l

        flash1 = max(zero, k1 * w(i,j, m15l) * ( kg_to_g * qgraup(i,j, m15l) ) )

      end if ! m15l

      flash2 = max(zero, k2 * ( gwp(i,j) + tiwp(i,j) ) )

      flash(i,j) = fivemin2sec * ( (mccaul_r1 * flash1) + (mccaul_r2 * flash2) )

      ! Determine the value of each component for diagnostics
      ! (these will have been initialised to zero in electric_init).

      fr1_mc(i,j) = fivemin2sec * flash1
      fr2_mc(i,j) = fivemin2sec * flash2

    end if ! storm_field

  end do ! i
end do ! j
!$OMP end PARALLEL do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine fr_mccaul

end module fr_mccaul_mod
