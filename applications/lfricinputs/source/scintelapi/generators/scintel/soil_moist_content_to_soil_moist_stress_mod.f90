! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module soil_moist_content_to_soil_moist_stress_mod
!
! This module contains generator that converts soil moisture content
! in kg per m2 to soil moisture stress
!

use dependency_graph_mod, only: dependency_graph

implicit none

private

public :: soil_moist_content_to_soil_moist_stress, soil_moist_stress

contains

subroutine  soil_moist_content_to_soil_moist_stress(dep_graph)
!
! This generator calculates a linear combination of the two input fields in a
! dependency graph and store result in the output field.
!

use gen_io_check_mod,       only: gen_io_check
use field_mod,              only: field_type, field_proxy_type
use log_mod,                only: log_event, log_scratch_space, LOG_LEVEL_ERROR
use fs_continuity_mod,      only: W3
use function_space_mod,     only: function_space_type

use lfricinp_surface_parameters_mod, only: dzsoil, sm_levels
use lfricinp_physics_constants_mod,  only: density_h2o

implicit none

!
! Argument definitions:
!
! Dependency graph to be processed
class(dependency_graph), intent(in out) :: dep_graph

!
! Local variables
!
! Field pointers to use
type(field_type), pointer :: field_soil_moist_content => null(),               &
                             field_soil_moist_content_crit => null(),          &
                             field_soil_moist_content_wilt => null(),          &
                             field_soil_moist_stress => null()

! Field proxies for accessing data
type(field_proxy_type)    :: field_proxy_soil_moist_content,                   &
                             field_proxy_soil_moist_content_crit,              &
                             field_proxy_soil_moist_content_wilt,              &
                             field_proxy_soil_moist_stress

type(function_space_type), pointer :: field_func_space => null()

!
! Local integers
integer :: i, j, l, size_horisontal, ndata

!
! Perform some initial input checks
!
call gen_io_check(                                                             &
                  dep_graph=dep_graph,                                         &
                  input_field_no=3,                                            &
                  input_field_fs=[W3, W3, W3],                                 &
                  output_field_no=1,                                           &
                  output_field_fs=[W3]                                         &
                 )
!
! Done with initial input checks
!

! Set up field pointers
field_soil_moist_content => dep_graph % input_field(1) % field_ptr
field_soil_moist_content_crit => dep_graph % input_field(2) % field_ptr
field_soil_moist_content_wilt => dep_graph % input_field(3) % field_ptr
field_soil_moist_stress => dep_graph % output_field(1) % field_ptr

! Set up field proxies
field_proxy_soil_moist_content = field_soil_moist_content % get_proxy()
field_proxy_soil_moist_content_crit =                                          &
                                    field_soil_moist_content_crit % get_proxy()
field_proxy_soil_moist_content_wilt =                                          &
                                    field_soil_moist_content_wilt % get_proxy()
field_proxy_soil_moist_stress = field_soil_moist_stress % get_proxy()

! Set some size variables
size_horisontal = size(field_proxy_soil_moist_content_crit%data)
field_func_space => field_soil_moist_stress % get_function_space()
ndata = field_func_space % get_ndata()
nullify(field_func_space)

if ( ndata /= sm_levels ) then
  write(log_scratch_space,'(A)')                                               &
       'error: NDATA for soil moisture field is not equal to sm_levels'
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

! Calculate soil moisture stress from soil moisture content
l = 0
do i = 1, size_horisontal
  do j = 1, ndata
    l = l + 1
    field_proxy_soil_moist_stress % data(l) =                                  &
              soil_moist_stress(field_proxy_soil_moist_content % data(l),      &
                                field_proxy_soil_moist_content_crit % data(i), &
                                field_proxy_soil_moist_content_wilt % data(i), &
                                dzsoil(j), density_h2o)
  end do
end do

! Nullify field pointers
nullify(field_soil_moist_content)
nullify(field_soil_moist_content_crit)
nullify(field_soil_moist_content_wilt)
nullify(field_soil_moist_stress)

end subroutine soil_moist_content_to_soil_moist_stress


function soil_moist_stress(soil_moist_content, soil_moist_content_crit,        &
                           soil_moist_content_wilt, dz, rho_water)             &
                          result (sm_stress)

use constants_def_mod, only: r_def, rmdi
use mdi_mod, only: is_rmdi

implicit none

! Arguments
real(kind=r_def), intent(in) :: soil_moist_content, soil_moist_content_crit,   &
                                soil_moist_content_wilt, dz, rho_water

! The result
real(kind=r_def) :: sm_stress, tiny_real

tiny_real = tiny(1.0_r_def)

if ( is_rmdi(soil_moist_content)      .or.                                     &
     is_rmdi(soil_moist_content_crit) .or.                                     &
     is_rmdi(soil_moist_content_wilt) .or.                                     &
     abs(soil_moist_content_crit - soil_moist_content_wilt) < tiny_real ) then

  sm_stress = rmdi

else

  sm_stress = (soil_moist_content / (dz * rho_water) -                         &
               soil_moist_content_wilt) /                                      &
              (soil_moist_content_crit - soil_moist_content_wilt)

end if

end function soil_moist_stress

end module soil_moist_content_to_soil_moist_stress_mod
