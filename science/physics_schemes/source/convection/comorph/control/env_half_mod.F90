! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module env_half_mod

implicit none

! Module stores addresses in super-arrays containing
! environment fields which are required at the next half-level
! (layer interface) in the convection parcel model.

! Number of fields needed at half-level
integer, parameter :: n_env_half = 2

! Indices of the fields in the super-array
integer, parameter :: i_wind_w_half = 1     ! Vertical velocity
integer, parameter :: i_virt_temp = 2  ! Virtual temperature


contains


!----------------------------------------------------------------
! Subroutine to interpolate the required fields onto half-level
!----------------------------------------------------------------
subroutine env_half_interp( l_last_level, k_full, n_points_super, cmpr,        &
                            height_full, height_k, height_half,                &
                            wind_w, virt_temp,                                 &
                            wind_w_k, virt_temp_k,                             &
                            env_half )

use cmpr_type_mod, only: cmpr_type
use comorph_constants_mod, only: nx_full, ny_full, k_bot_conv, k_top_conv,     &
                                 real_cvprec, real_hmprec, one
use compress_mod, only: compress

implicit none

! Flag for reached the final model-level
logical, intent(in) :: l_last_level

! k-index of neighbouring full model-level needed for interpolating to k_half
integer, intent(in) :: k_full

! Array dimension of the env_half super-array
integer, intent(in) :: n_points_super

! Structure containing compression list info
type(cmpr_type), intent(in) :: cmpr

! Full 3-D array of model-level heights
! (needed so we can extract values at k_full for interpolation)
real(kind=real_hmprec), pointer, intent(in) :: height_full(:,:,:)

! Already compressed heights of the current level and half-level
real(kind=real_cvprec), intent(in) :: height_k                                 &
                                      ( cmpr%n_points )
real(kind=real_cvprec), intent(in) :: height_half                              &
                                      ( cmpr%n_points )


! Full 3-D array of vertical velocity
real(kind=real_hmprec), pointer, intent(in) :: wind_w(:,:,:)
! Latest virtual temperature 3D array
real(kind=real_hmprec), intent(in) :: virt_temp                                &
       ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Note: the array dimensions of virt_temp are known, because it
! is a local array within comorph.
! But wind_w is passed in from outside and may
! or may not have halos / BC points that we don't want to mess
! with by mistake.  Hence it is passed in as a pointer,
! to ensure the full array bounds specs are passed in with it.

! Vertical velocity already compressed onto convecting points at level k
real(kind=real_cvprec), intent(in) :: wind_w_k                                 &
                                      ( cmpr%n_points )
! Virtual temperature  already compressed onto convecting points at level k
real(kind=real_cvprec), intent(in) :: virt_temp_k                              &
                                      ( cmpr%n_points )

! Super-array to contain fields interpolated to the half-level
real(kind=real_cvprec), intent(out) :: env_half                                &
                                  ( n_points_super, n_env_half )

! Work array stores interpolation weight
real(kind=real_cvprec) :: weight ( cmpr%n_points )

! Work array for compressing stuff
real(kind=real_cvprec) :: work_cmpr ( cmpr%n_points )

! Array lower and upper bounds
integer :: lb(3), ub(3)

! Loop counter
integer :: ic


! TEMPORARY CODE TO PRESERVE KGO:
! If in the last level but interpolating to find Tv at prev, we should
! enter the ELSE branch below and interpolate between k and k-dk;
! this code wrongly sets Tv at prev equal to Tv(k); fix this soon...
!IF ( k_full > k_top_conv .OR. k_full < k_bot_conv ) THEN
if ( l_last_level ) then
  ! If trying to interpolate beyond the first or last full level,
  ! just assume fields are constant
  ! beyond and copy the fields at the current level, as the
  ! required fields may not exist at level k_full, making
  ! interpolation onto k_half impossible.

  do ic = 1, cmpr%n_points
    env_half(ic,i_wind_w_half) = wind_w_k(ic)
  end do
  do ic = 1, cmpr%n_points
    env_half(ic,i_virt_temp) = virt_temp_k(ic)
  end do

else  ! ( k_full <= k_top_conv .AND. k_full >== k_bot_conv )
  ! If not at the top or bottom, interpolate as usual...

  ! Compress height of next model-level into the work array
  lb = lbound(height_full)
  ub = ubound(height_full)
  call compress( cmpr, lb(1:2), ub(1:2), height_full(:,:,k_full), work_cmpr )

  ! Calculate interpolation weight
  do ic = 1, cmpr%n_points
    weight(ic) = ( height_half(ic) - height_k(ic) )                            &
               / ( work_cmpr(ic) - height_k(ic) )
  end do

  ! Compress vertical velocities from the next full model-level
  lb = lbound(wind_w)
  ub = ubound(wind_w)
  call compress( cmpr, lb(1:2), ub(1:2), wind_w(:,:,k_full), work_cmpr )
  ! Interpolate compressed vertical velocities onto half-level
  do ic = 1, cmpr%n_points
    env_half(ic,i_wind_w_half)                                                 &
      = (one-weight(ic)) * wind_w_k(ic)                                        &
      +      weight(ic)  * work_cmpr(ic)
  end do

  ! Compress virtual temperatures from the next full model-level
  lb = [1,1,k_bot_conv]
  ub = [nx_full,ny_full,k_top_conv]
  call compress( cmpr, lb(1:2), ub(1:2), virt_temp(:,:,k_full), work_cmpr )
  ! Interpolate compressed virtual temperatures onto half-level
  do ic = 1, cmpr%n_points
    env_half(ic,i_virt_temp)                                                   &
      = (one-weight(ic)) * virt_temp_k(ic)                                     &
      +      weight(ic)  * work_cmpr(ic)
  end do

end if  ! ( k_full <= k_top_conv .AND. k_full >== k_bot_conv )


return
end subroutine env_half_interp


end module env_half_mod
