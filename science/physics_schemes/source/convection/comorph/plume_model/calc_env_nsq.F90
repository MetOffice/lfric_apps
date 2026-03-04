! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_env_nsq_mod

implicit none

contains

! Subroutine to calculate the environment dry static stability N^2
! at the current half-level step within the convective plume integration
subroutine calc_env_nsq( n_points, n_points_env, l_to_full_level,              &
                         env_prev_virt_temp, env_next_virt_temp,               &
                         grid_prev, grid_next, env_k_fields,                   &
                         Nsq_dry )

use comorph_constants_mod, only: real_cvprec, gravity
use grid_type_mod, only: n_grid, i_height, i_pressure
use fields_type_mod, only: n_fields, i_temperature, i_q_vap,                   &
                           i_qc_first, i_qc_last
use dry_adiabat_mod, only: dry_adiabat

implicit none

! Number of points
integer, intent(in) :: n_points

! Length of list of points in re-used super-arrays
integer, intent(in) :: n_points_env

! Flag for first or 2nd half of the level-step
! (whether it ends on a full-level or half-level)
logical, intent(in) :: l_to_full_level

! Environment virtual temperatures at start and end of level-step
real(kind=real_cvprec), intent(in) :: env_prev_virt_temp(n_points)
real(kind=real_cvprec), intent(in) :: env_next_virt_temp(n_points)

! Heights and pressures at start and end of level-step
real(kind=real_cvprec), intent(in) :: grid_prev ( n_points_env, n_grid )
real(kind=real_cvprec), intent(in) :: grid_next ( n_points_env, n_grid )

! Primary fields at the nearest full level k
real(kind=real_cvprec), intent(in) :: env_k_fields ( n_points_env, n_fields )

! Output environment dry static stability
real(kind=real_cvprec), intent(out) :: Nsq_dry(n_points)

! Virtual temperature of test parcel used to compute Nsq
real(kind=real_cvprec) :: virt_temp_test(n_points)

! Loop counter
integer :: ic


! Calculate static stability of the environment...
do ic = 1, n_points
  virt_temp_test(ic) = env_k_fields(ic,i_temperature)
end do
if ( l_to_full_level ) then
  ! First half of level-step; env_k_fields are valid at "next"

  ! Calculate temperature of environment air at k if subsided to prev
  call dry_adiabat( n_points, n_points_env,                                    &
                    grid_next(:,i_pressure), grid_prev(:,i_pressure),          &
                    env_k_fields(:,i_q_vap),                                   &
                    env_k_fields(:,i_qc_first:i_qc_last),                      &
                    virt_temp_test )

  do ic = 1, n_points
    ! Convert to virtual temperature; we already have both T and Tv at the
    ! current level, so just scale by their ratio.
    virt_temp_test(ic) = virt_temp_test(ic)                                    &
      * env_next_virt_temp(ic) / env_k_fields(ic,i_temperature)
    ! Compute N^2
    Nsq_dry(ic) = ( gravity / env_prev_virt_temp(ic) )                         &
                * ( virt_temp_test(ic) - env_prev_virt_temp(ic) )              &
                / ( grid_next(ic,i_height) - grid_prev(ic,i_height) )
  end do

else  ! ( .NOT. l_to_full_level )
  ! Second half of level-step; env_k_fields are valid at "prev"

  ! Calculate temperature of environment air at k if lifted to next
  call dry_adiabat( n_points, n_points_env,                                    &
                    grid_prev(:,i_pressure), grid_next(:,i_pressure),          &
                    env_k_fields(:,i_q_vap),                                   &
                    env_k_fields(:,i_qc_first:i_qc_last),                      &
                    virt_temp_test )

  do ic = 1, n_points
    ! Convert to virtual temperature; we already have both T and Tv at the
    ! current level, so just scale by their ratio.
    virt_temp_test(ic) = virt_temp_test(ic)                                    &
      * env_prev_virt_temp(ic) / env_k_fields(ic,i_temperature)
    ! Compute N^2
    Nsq_dry(ic) = ( gravity / env_next_virt_temp(ic) )                         &
                * ( env_next_virt_temp(ic) - virt_temp_test(ic) )              &
                / ( grid_next(ic,i_height) - grid_prev(ic,i_height) )
  end do

end if  ! ( .NOT. l_to_full_level )


return
end subroutine calc_env_nsq

end module calc_env_nsq_mod
