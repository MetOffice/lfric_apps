! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module update_edge_virt_temp_mod

implicit none

contains

! Subroutine to update the virtual temperature at the edge of the in-parcel
! assumed PDF.  Relaxes the edge Tv towards the environment Tv.
! Care needs to be taken over exactly how to do this; e.g. the relaxation
! can cause very strange PDF shapes if the PDF is too narrow, or the core
! was closer to the environment value than the edge was.
! The final value of edge_virt_temp is crucial for determining the
! power of the assumed power-law PDF shape at the next level, so this
! calculation has a significant impact on the detrainment rate profile.
subroutine update_edge_virt_temp( n_points, max_points, l_down,                &
                                  core_mean_ratio,                             &
                                  sublevs, i_next, par_conv_super )

use comorph_constants_mod, only: real_cvprec, one
use sublevs_mod, only: max_sublevs, n_sublev_vars,                             &
                       j_delta_tv, j_env_tv, j_mean_buoy, j_core_buoy
use parcel_type_mod, only: n_par, i_edge_virt_temp

implicit none

! Number of points
integer, intent(in) :: n_points
! Super-array dimension
integer, intent(in) :: max_points

! Flag for downdraft vs updraft
logical, intent(in) :: l_down

! Previous ratio of the parcel core buoyancy over the parcel mean buoyancy
! (sets the assumed PDF shape)
real(kind=real_cvprec), intent(in) :: core_mean_ratio(n_points)

! Super-array storing the parcel buoyancies and other properties
! at up to 4 sub-level heights within the current level-step:
! a) Previous model-level interface
! b) Next model-level interface
! c) Height where parcel mean properties first hit saturation
! d) Height where parcel core first hit saturation
real(kind=real_cvprec), intent(in) :: sublevs                                  &
                                     ( n_points, n_sublev_vars, max_sublevs )

! Address of next model-level interface in sublevs
integer, intent(in) :: i_next(n_points)

! Super-array containing miscellaneous parcel properties
real(kind=real_cvprec), intent(in out) :: par_conv_super                       &
                                          ( max_points, n_par )

! Environment virtual temperature at end of half-level-step, after subsidence
real(kind=real_cvprec) :: env_next_virt_temp(n_points)

! Loop counter
integer :: ic


! Calculate final parcel edge virtual temperature,
! used for recalculating the core-mean-ratio at the next level.
!
! cmr = ( core_buoy - edge_buoy ) / ( mean_buoy - edge_buoy )
! where buoy = Tv_par - (Tv_env + dTv_sub)
!
! => edge_buoy ( 1 - cmr ) = core_buoy - mean_buoy cmr
!                          = core_buoy - mean_buoy + mean_buoy ( 1 - cmr )
! => edge_buoy = mean_buoy - ( core_buoy - mean_buoy ) / ( cmr - 1 )
!
do ic = 1, n_points
  ! Set updated Tv env = Tv_env_next + delta_Tv
  env_next_virt_temp(ic) = sublevs(ic,j_env_tv,i_next(ic))                     &
                         + sublevs(ic,j_delta_tv,i_next(ic))
  ! Set updated parel edge Tv = edge_buoy + updated env Tv
  par_conv_super(ic,i_edge_virt_temp)                                          &
    = ( sublevs(ic,j_mean_buoy,i_next(ic)) + env_next_virt_temp(ic) )          &
    - ( sublevs(ic,j_core_buoy,i_next(ic))                                     &
      - sublevs(ic,j_mean_buoy,i_next(ic)) ) / ( core_mean_ratio(ic) - one )
end do

! If the edge Tv has pulled away from the environment Tv (become buoyant
! for updrafts, negatively buoyant for downdrafts), adjust the edge Tv
! back to the environment Tv, changing the shape of the PDF.
if ( l_down ) then
  do ic = 1, n_points
    par_conv_super(ic,i_edge_virt_temp)                                        &
      = max( par_conv_super(ic,i_edge_virt_temp), env_next_virt_temp(ic) )
  end do
else
  do ic = 1, n_points
    par_conv_super(ic,i_edge_virt_temp)                                        &
      = min( par_conv_super(ic,i_edge_virt_temp), env_next_virt_temp(ic) )
  end do
end if


return
end subroutine update_edge_virt_temp

end module update_edge_virt_temp_mod

