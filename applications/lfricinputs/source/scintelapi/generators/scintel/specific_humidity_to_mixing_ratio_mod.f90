! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE specific_humidity_to_mixing_ratio_mod
!
! This module contains a generator to convert specific humidities and a total
! density into the equivalent mixing ratios and dry density.
!

USE dependency_graph_mod, ONLY: dependency_graph

IMPLICIT NONE

CONTAINS

SUBROUTINE specific_humidity_to_mixing_ratio(dep_graph)
!
! This generator multiplies the input field in a dependency graph by a constant
! value and save the result in the output field of the dependency graph.
!

USE gen_io_check_mod,       ONLY: gen_io_check
USE field_mod,              ONLY: field_type
USE map_um_lbc_inputs_alg_mod, ONLY: map_um_lbc_inputs

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
TYPE(field_type), POINTER :: field_q => NULL(), &
                             field_qcl => NULL(), &
                             field_qcf => NULL(), &
                             field_qrain => NULL(), &
                             field_rho_r2 => NULL(), &
                             field_m_v => NULL(), &
                             field_m_cl => NULL(), &
                             field_m_ci => NULL(), &
                             field_m_r => NULL(), &
                             field_rho => NULL()

!
! Perform some initial input checks
!
CALL gen_io_check(                                                             &
                  dep_graph=dep_graph,                                         &
                  input_field_no=5,                                            &
                  output_field_no=5,                                           &
                  parameter_no=0                                               &
                 )
!
! Done with initial field checks
!


! Convert specific humidities to mixing ratios through lfric module.
field_m_v => dep_graph % output_field(1) % field_ptr
field_m_cl => dep_graph % output_field(2) % field_ptr
field_m_ci => dep_graph % output_field(3) % field_ptr
field_m_r => dep_graph % output_field(4) % field_ptr
field_rho => dep_graph % output_field(5) % field_ptr

field_q => dep_graph % input_field(1) % field_ptr
field_qcl => dep_graph % input_field(2) % field_ptr
field_qcf => dep_graph % input_field(3) % field_ptr
field_qrain => dep_graph % input_field(4) % field_ptr
field_rho_r2 => dep_graph % input_field(5) % field_ptr

CALL map_um_lbc_inputs (field_q, field_qcl, field_qcf, field_qrain, &
                        field_rho_r2, field_m_v, field_m_cl, &
                        field_m_ci, field_m_r, field_rho)

! Nullify field pointers
NULLIFY(field_q, field_qcl, field_qcf, field_qrain, field_rho_r2, field_m_v, &
        field_m_ci, field_m_cl, field_m_r, field_rho)

END SUBROUTINE specific_humidity_to_mixing_ratio

END MODULE specific_humidity_to_mixing_ratio_mod
