! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module specific_humidity_to_mixing_ratio_mod
!
! This module contains a generator to convert specific humidities and a total
! density into the equivalent mixing ratios and dry density.
!

use dependency_graph_mod, only: dependency_graph

implicit none

contains

subroutine specific_humidity_to_mixing_ratio(dep_graph)
!
! This generator multiplies the input field in a dependency graph by a constant
! value and save the result in the output field of the dependency graph.
!

use gen_io_check_mod,       only: gen_io_check
use field_mod,              only: field_type
use map_um_lbc_inputs_alg_mod, only: map_um_lbc_inputs

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
type(field_type), pointer :: field_q => null(), &
                             field_qcl => null(), &
                             field_qcf => null(), &
                             field_qrain => null(), &
                             field_rho_r2 => null(), &
                             field_m_v => null(), &
                             field_m_cl => null(), &
                             field_m_s => null(), &
                             field_m_r => null(), &
                             field_rho => null()

!
! Perform some initial input checks
!
call gen_io_check(                                                             &
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
field_m_s => dep_graph % output_field(3) % field_ptr
field_m_r => dep_graph % output_field(4) % field_ptr
field_rho => dep_graph % output_field(5) % field_ptr

field_q => dep_graph % input_field(1) % field_ptr
field_qcl => dep_graph % input_field(2) % field_ptr
field_qcf => dep_graph % input_field(3) % field_ptr
field_qrain => dep_graph % input_field(4) % field_ptr
field_rho_r2 => dep_graph % input_field(5) % field_ptr

call map_um_lbc_inputs (field_q, field_qcl, field_qcf, field_qrain, &
                        field_rho_r2, field_m_v, field_m_cl, &
                        field_m_s, field_m_r, field_rho)

! Nullify field pointers
nullify(field_q, field_qcl, field_qcf, field_qrain, field_rho_r2, field_m_v, &
        field_m_s, field_m_cl, field_m_r, field_rho)

end subroutine specific_humidity_to_mixing_ratio

end module specific_humidity_to_mixing_ratio_mod
