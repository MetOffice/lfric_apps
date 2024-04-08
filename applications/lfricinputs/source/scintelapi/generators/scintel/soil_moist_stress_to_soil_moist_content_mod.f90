! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE soil_moist_stress_to_soil_moist_content_mod
!
! This module contains generator that converts soil moisture stress
! to soil moisture content in kg per m2
!

USE dependency_graph_mod, ONLY: dependency_graph

IMPLICIT NONE

PRIVATE

PUBLIC :: soil_moist_stress_to_soil_moist_content, soil_moist_content

CONTAINS

SUBROUTINE  soil_moist_stress_to_soil_moist_content(dep_graph)
!
! This generator calculates a linear combination of the two input fields in a
! dependency graph and store result in the output field.
!

USE gen_io_check_mod,       ONLY: gen_io_check
USE field_mod,              ONLY: field_type, field_proxy_type
USE log_mod,                ONLY: log_event, log_scratch_space, LOG_LEVEL_ERROR
USE fs_continuity_mod,      ONLY: W3
USE function_space_mod,     ONLY: function_space_type

USE lfricinp_surface_parameters_mod, ONLY: dzsoil, sm_levels
USE lfricinp_physics_constants_mod,  ONLY: density_h2o

IMPLICIT NONE

!
! Argument definitions:
!
! Dependency graph to be processed
CLASS(dependency_graph), INTENT(IN OUT) :: dep_graph

!
! Local variables
!
! Field pointers to use
TYPE(field_type), POINTER :: field_soil_moist_stress => NULL(),                &
                             field_soil_moist_content_crit => NULL(),          &
                             field_soil_moist_content_wilt => NULL(),          &
                             field_soil_moist_content => NULL()

! Field proxies for accessing data
TYPE(field_proxy_type)    :: field_proxy_soil_moist_stress,                    &
                             field_proxy_soil_moist_content_crit,              &
                             field_proxy_soil_moist_content_wilt,              &
                             field_proxy_soil_moist_content

TYPE(function_space_type), POINTER :: field_func_space => NULL()

!
! Local integers
INTEGER :: i, j, l, size_horisontal, ndata

!
! Perform some initial input checks
!
CALL gen_io_check(                                                             &
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
field_soil_moist_stress => dep_graph % input_field(1) % field_ptr
field_soil_moist_content_crit => dep_graph % input_field(2) % field_ptr
field_soil_moist_content_wilt => dep_graph % input_field(3) % field_ptr
field_soil_moist_content => dep_graph % output_field(1) % field_ptr

! Set up field proxies
field_proxy_soil_moist_stress = field_soil_moist_stress % get_proxy()
field_proxy_soil_moist_content_crit =                                          &
                                    field_soil_moist_content_crit % get_proxy()
field_proxy_soil_moist_content_wilt =                                          &
                                    field_soil_moist_content_wilt % get_proxy()
field_proxy_soil_moist_content = field_soil_moist_content % get_proxy()

! Set some size variables
size_horisontal = SIZE(field_proxy_soil_moist_content_crit%data)
field_func_space => field_soil_moist_content % get_function_space()
ndata = field_func_space % get_ndata()
NULLIFY(field_func_space)

IF ( ndata /= sm_levels ) THEN
  WRITE(log_scratch_space,'(A)')                                               &
       'ERROR: NDATA for soil moisture field is not equal to sm_levels'
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

! Calculate soil moisture stress from soil moisture content
l = 0
DO i = 1, size_horisontal
  DO j = 1, ndata
    l = l + 1
    field_proxy_soil_moist_content % data(l) =                                 &
              soil_moist_content(field_proxy_soil_moist_stress % data(l),      &
                                field_proxy_soil_moist_content_crit % data(i), &
                                field_proxy_soil_moist_content_wilt % data(i), &
                                dzsoil(j), density_h2o)
  END DO
END DO

! Nullify field pointers
NULLIFY(field_soil_moist_stress)
NULLIFY(field_soil_moist_content_crit)
NULLIFY(field_soil_moist_content_wilt)
NULLIFY(field_soil_moist_content)

END SUBROUTINE soil_moist_stress_to_soil_moist_content


FUNCTION soil_moist_content(soil_moist_stress, soil_moist_content_crit,        &
                            soil_moist_content_wilt, dz, rho_water)            &
                           RESULT (sm_content)

USE constants_def_mod, ONLY: r_def, rmdi
USE mdi_mod, ONLY: is_rmdi

IMPLICIT NONE

! Arguments
REAL(KIND=r_def), INTENT(IN) :: soil_moist_stress, soil_moist_content_crit,    &
                                soil_moist_content_wilt, dz, rho_water
! The result
REAL(KIND=r_def) :: sm_content, tiny_real

tiny_real = TINY(1.0_r_def)

IF ( is_rmdi(soil_moist_stress)       .OR.                                     &
     is_rmdi(soil_moist_content_crit) .OR.                                     &
     is_rmdi(soil_moist_content_wilt) ) THEN

  sm_content = rmdi

ELSE


  sm_content = soil_moist_content_wilt + soil_moist_stress *                   &
               (soil_moist_content_crit - soil_moist_content_wilt)
  sm_content = sm_content * dz * rho_water

END IF

END FUNCTION soil_moist_content

END MODULE soil_moist_stress_to_soil_moist_content_mod
