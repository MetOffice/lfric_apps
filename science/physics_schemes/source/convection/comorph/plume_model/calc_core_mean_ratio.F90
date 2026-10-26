! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_core_mean_ratio_mod

implicit none

contains

! Subroutine to calculate the ratio of core buoyancy over mean buoyancy,
! used to contruct the in-parcel assumed-PDF
subroutine calc_core_mean_ratio( n_points, n_points_super, l_down,             &
                                 par_mean_virt_temp, par_core_virt_temp,       &
                                 par_conv_super, core_mean_ratio )

use comorph_constants_mod, only: real_cvprec, min_cmr, max_cmr, zero
use parcel_type_mod, only: n_par, i_edge_virt_temp

implicit none

! Number of points
integer, intent(in) :: n_points
! Length of super-array
integer, intent(in) :: n_points_super

! Flag for downdraft versus updraft
logical, intent(in) :: l_down

! Parcel mean virtual temperature at next model-level interface
real(kind=real_cvprec), intent(in) :: par_mean_virt_temp(n_points)

! Parcel core virtual temperature at next model-level interface
real(kind=real_cvprec), intent(in) :: par_core_virt_temp(n_points)

! Parcel properties super-array
real(kind=real_cvprec), intent(in) :: par_conv_super ( n_points_super, n_par )

! Ratio of core buoyancy over mean buoyancy
real(kind=real_cvprec), intent(out) :: core_mean_ratio(n_points)

! Parcel mean and core buoyancy, relative to parcel edge
real(kind=real_cvprec) :: mean_buoy
real(kind=real_cvprec) :: core_buoy

! Loop counter
integer :: ic


! Let cmr = ( Tv_core - Tv_edge )
!         / ( Tv_mean - Tv_edge )
! where Tv_edge is stored in edge_virt_temp

! Calculate mean and core buoyancies, relative to the previous parcel
! edge virtual temperature
do ic = 1, n_points
  ! Initialise to max value
  core_mean_ratio(ic) = max_cmr
end do

! Safety checks depend on whether updraft or downdraft
if ( l_down ) then
  ! Downdraft: require core and mean to be negatively buoyant
  do ic = 1, n_points
    mean_buoy = par_mean_virt_temp(ic) - par_conv_super(ic,i_edge_virt_temp)
    core_buoy = par_core_virt_temp(ic) - par_conv_super(ic,i_edge_virt_temp)
    if ( mean_buoy < zero .and. core_buoy < zero ) then
      core_mean_ratio(ic) = core_buoy / mean_buoy
    end if
  end do
else
  ! Updraft: require core and mean to be positively buoyant
  do ic = 1, n_points
    mean_buoy = par_mean_virt_temp(ic) - par_conv_super(ic,i_edge_virt_temp)
    core_buoy = par_core_virt_temp(ic) - par_conv_super(ic,i_edge_virt_temp)
    if ( mean_buoy > zero .and. core_buoy > zero ) then
      core_mean_ratio(ic) = core_buoy / mean_buoy
    end if
  end do
end if

! Don't let this go outside specified bounds,
! for safety in the detrainment calculations.
do ic = 1, n_points
  core_mean_ratio(ic) = max( min( core_mean_ratio(ic), max_cmr ), min_cmr )
end do


return
end subroutine calc_core_mean_ratio

end module calc_core_mean_ratio_mod
