! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module set_edge_virt_temp_mod

implicit none

contains

! Subroutine to set the parcel edge virtual temperature, consistent
! with the assumed-PDF of Tv used for detrainment.
subroutine set_edge_virt_temp( n_points, n_buoy_vars, max_buoy_heights,        &
                               i_next, buoyancy_super, l_down,                 &
                               env_next_virt_temp, core_mean_ratio,            &
                               edge_virt_temp )

use comorph_constants_mod, only: real_cvprec, one
use buoyancy_mod, only: i_mean_buoy, i_core_buoy

implicit none

! Number of points
integer, intent(in) :: n_points

! Dimensions of super-array containing buoyancies at different sub-level steps
integer, intent(in) :: n_buoy_vars
integer, intent(in) :: max_buoy_heights

! Index of next model-level interface in the buoyancy super-array
integer, intent(in) :: i_next(n_points)

! Buoyancy super-array
real(kind=real_cvprec), intent(in) :: buoyancy_super                           &
                    ( n_points, n_buoy_vars, max_buoy_heights )

! Flag for downdrafts versus updrafts
logical, intent(in) :: l_down

! Environment virtual temperature at next model-level interface
real(kind=real_cvprec), intent(in) :: env_next_virt_temp(n_points)

! Ratio of core buoyancy parcel-mean buoyancy
real(kind=real_cvprec), intent(in) :: core_mean_ratio(n_points)

! Virtual temperature at the parcel edge
real(kind=real_cvprec), intent(out) :: edge_virt_temp(n_points)

! Loop counter
integer :: ic


! Set new parcel edge virtual temperature at next model-level interface
do ic = 1, n_points
  ! (Tv_core - Tv_edge) = core_rat (Tv_mean - Tv_edge)
  ! => Tv_edge (core_rat - 1) = Tv_mean core_rat - Tv_core
  !                           = Tv_mean (core_rat - 1) - (Tv_core - Tv_mean)
  ! Tv_edge = Tv_mean - (Tv_core - Tv_mean) / (core_rat - 1)
  edge_virt_temp(ic) = ( buoyancy_super(ic,i_mean_buoy,i_next(ic))             &
                       + env_next_virt_temp(ic) )                              &
                     - ( buoyancy_super(ic,i_core_buoy,i_next(ic))             &
                       - buoyancy_super(ic,i_mean_buoy,i_next(ic)) )           &
                     / ( core_mean_ratio(ic) - one )
end do

! If the edge Tv has pulled away from the environment Tv (become buoyant
! for updrafts, negatively buoyant for downdrafts), adjust the edge Tv
! back to the environment Tv, changing the shape of the PDF.
if ( l_down ) then
  do ic = 1, n_points
    edge_virt_temp(ic) = max( edge_virt_temp(ic), env_next_virt_temp(ic) )
  end do
else
  do ic = 1, n_points
    edge_virt_temp(ic) = min( edge_virt_temp(ic), env_next_virt_temp(ic) )
  end do
end if


return
end subroutine set_edge_virt_temp

end module set_edge_virt_temp_mod
