! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE constants_def_mod
!
! This is a bespoke module for global parameters used in the API code. The
! intention is to have a separate module from the LFRic constants_mod module,
! which can be accessed for both LFRic and API specific global parameters.
!

! Use some constants from LFRic infrastructure
USE constants_mod, ONLY: lfric_r_def => r_def,                                 &
                         lfric_i_def => i_def,                                 &
                         real_missing_data_indicator => rmdi

IMPLICIT NONE

! Input string size definitions
INTEGER, PARAMETER :: field_kind_name_len = 20
INTEGER, PARAMETER :: field_name_len = 64
INTEGER, PARAMETER :: gen_id_len = 40
INTEGER, PARAMETER :: genpar_len = 1000
INTEGER, PARAMETER :: field_dim_len = 2
INTEGER, PARAMETER :: file_name_len = 512

! Array sizes
INTEGER, PARAMETER :: field_id_list_max_size = 5
INTEGER, PARAMETER :: max_no_dependency_graphs = 50
INTEGER, PARAMETER :: max_no_fields = 100

! Global parameters used in the API
CHARACTER,         PARAMETER :: empty_string = ''
INTEGER,           PARAMETER :: r_def = lfric_r_def
INTEGER,           PARAMETER :: i_def = lfric_i_def
REAL(KIND=r_def),  PARAMETER :: rmdi = real_missing_data_indicator

END MODULE constants_def_mod
