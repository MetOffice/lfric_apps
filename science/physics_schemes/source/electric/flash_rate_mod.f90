! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module flash_rate_mod

! Purpose: Calls subroutines to calculate thunderstorm flash rate
!          dependent on the method requested by the user.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Electric

use atm_fields_bounds_mod,  only: wdims_s, tdims
use electric_inputs_mod,    only: electric_method, em_gwp, em_mccaul
use electric_constants_mod, only: minus5

use fr_gwp_mod,             only: fr_gwp
use fr_mccaul_mod,          only: fr_mccaul
use calc_wp_below_t_mod,    only: calc_wp_below_t

! Dr Hook modules
use yomhook,                only: lhook, dr_hook
use parkind1,               only: jprb, jpim

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='FLASH_RATE_MOD'

contains

subroutine flash_rate( storm_field, qgraup, t, w, rhodz, gwp, tiwp, flash,     &
                       fr1_mc, fr2_mc )

implicit none

real(kind=real_umphys), intent(in) :: qgraup(                                  &
                                tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end,                   &
                                            1 : tdims%k_end      )

real(kind=real_umphys), intent(in) :: t(                                       &
                                tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end,                   &
                                            1 : tdims%k_end      )

real(kind=real_umphys), intent(in) :: w(                                       &
                              wdims_s%i_start : wdims_s%i_end,                 &
                              wdims_s%j_start : wdims_s%j_end,                 &
                              wdims_s%k_start : wdims_s%k_end    )

real(kind=real_umphys), intent(in) :: rhodz(                                   &
                                tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end,                   &
                                            1 : tdims%k_end      )

real(kind=real_umphys), intent(in) :: gwp(                                     &
                                tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end      )

real(kind=real_umphys), intent(in) :: tiwp(                                    &
                                tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end      )

real(kind=real_umphys), intent(in out) :: flash(                               &
                                tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end      )

real(kind=real_umphys), intent(in out) :: fr1_mc(                              &
                                tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end      )

real(kind=real_umphys), intent(in out) :: fr2_mc(                              &
                                tdims%i_start : tdims%i_end,                   &
                                tdims%j_start : tdims%j_end      )

logical, intent(in) :: storm_field( tdims%i_start : tdims%i_end,               &
                                    tdims%j_start : tdims%j_end      )

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='FLASH_RATE'

!-----------------------------------------------------------
! Local variables
!-----------------------------------------------------------
real(kind=real_umphys) :: gwp_m5(                                              &
                 tdims%i_start : tdims%i_end,                                  &
                 tdims%j_start : tdims%j_end      )


!==================================================================
! Start the subroutine
!==================================================================

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!====================================================================
! Determine flash rate according to electric scheme
!====================================================================

select case (electric_method)

case (em_gwp) ! GWP relation

  ! First generate gwp below -5 C threshold
  call calc_wp_below_t(qgraup, rhodz, t, minus5, gwp_m5 )

  ! Next call the flash rate
  call fr_gwp(storm_field, gwp_m5, flash)

case (em_mccaul) ! McCaul et al (2009), Weather and Forecasting

  call fr_mccaul( storm_field, t, w, qgraup, gwp, tiwp, flash, fr1_mc, fr2_mc)

case DEFAULT ! GWP relation
             ! (this is your failsafe option in case electric method
             !  is not defined)

  ! First generate gwp below -5 C threshold
  call calc_wp_below_t(qgraup, rhodz, t, minus5, gwp_m5 )

  ! Next call the flash rate
  call fr_gwp(storm_field, gwp_m5, flash)

end select

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine flash_rate

end module flash_rate_mod
