! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module constants_def_mod
!
! This is a bespoke module for global parameters used in the API code. The
! intention is to have a separate module from the LFRic constants_mod module,
! which can be accessed for both LFRic and API specific global parameters.
!

! Use some constants from LFRic infrastructure
use constants_mod, only: lfric_r_def => r_def,                                 &
                         lfric_i_def => i_def,                                 &
                         real_missing_data_indicator => rmdi

implicit none

! Input string size definitions
integer, parameter :: field_kind_name_len = 20
integer, parameter :: field_name_len = 64
integer, parameter :: gen_id_len = 40
integer, parameter :: genpar_len = 1000
integer, parameter :: field_dim_len = 2
integer, parameter :: file_name_len = 512

! Array sizes
integer, parameter :: field_id_list_max_size = 5
integer, parameter :: max_no_dependency_graphs = 50
integer, parameter :: max_no_fields = 100

! Global parameters used in the API
character,         parameter :: empty_string = ''
integer,           parameter :: r_def = lfric_r_def
integer,           parameter :: i_def = lfric_i_def
real(kind=r_def),  parameter :: rmdi = real_missing_data_indicator

end module constants_def_mod
