! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module fr_gwp_mod

! Purpose: Calculates flash rate using a simple graupel water path
!          relationship

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Electric

use atm_fields_bounds_mod,  only: tdims
use electric_inputs_mod,    only: g1, g2
use cderived_mod,           only: delta_lambda, delta_phi
use level_heights_mod,      only: r_theta_levels
use trignometric_mod,       only: cos_theta_latitude

! Dr Hook modules
use yomhook,                only: lhook, dr_hook
use parkind1,               only: jprb, jpim

use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='FR_GWP_MOD'

contains

subroutine fr_gwp(storm_field, gwp, flash)

implicit none

real(kind=real_umphys), intent(in)  :: gwp(                                    &
                            tdims%i_start : tdims%i_end,                       &
                            tdims%j_start : tdims%j_end )
! Graupel water path (kg m-2)

real(kind=real_umphys), intent(in out) :: flash(                               &
                              tdims%i_start : tdims%i_end,                     &
                              tdims%j_start : tdims%j_end )
! Flash rate (s-1)

logical, intent(in) :: storm_field(tdims%i_start : tdims%i_end,                &
                                   tdims%j_start : tdims%j_end )

! Local variables
real(kind=real_umphys) :: gbs  ! grid box size (square metres)

real(kind=real_umphys), pointer :: xx_cos_theta_latitude(:,:)

integer :: i, j
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='FR_GWP'

!==================================================================
! Start the subroutine
!==================================================================

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Grid box size variable dependent on latitude

! Calculate the trignometric constant dependent
! on dynamical core

xx_cos_theta_latitude => cos_theta_latitude

! Perform calculate of grid box size

!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none) private(i,j,gbs)              &
!$OMP SHARED(tdims,r_theta_levels,delta_lambda,delta_phi,                      &
!$OMP xx_cos_theta_latitude,storm_field,flash,gwp,g1,g2)
do j = tdims%j_start, tdims%j_end
  do i = tdims%i_start, tdims%i_end

    if ( storm_field(i,j) ) then
      ! Determine grid box size
      gbs = r_theta_levels(i,j,1) * delta_lambda                               &
          * r_theta_levels(i,j,1) * delta_phi                                  &
          * xx_cos_theta_latitude(i,j)
      flash(i,j) = max( 0.0, ( gbs * g1 * gwp(i,j) ) - g2 )
    end if

  end do
end do
!$OMP end PARALLEL do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine fr_gwp

end module fr_gwp_mod
