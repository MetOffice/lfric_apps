! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_um_parameters_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY: int64

IMPLICIT NONE

! Filename length
INTEGER, PUBLIC, PARAMETER :: fnamelen = 512

! Message length
INTEGER, PARAMETER :: msglen = 512

! UM definitions of REAL and INTEGER in fieldsfiles
INTEGER, PUBLIC, PARAMETER :: um_real64  = SELECTED_REAL_KIND(15,307)
INTEGER, PUBLIC, PARAMETER :: um_integer64 = SELECTED_INT_KIND(15)
INTEGER, PUBLIC, PARAMETER :: um_real32  = SELECTED_REAL_KIND(6,37)
INTEGER, PUBLIC, PARAMETER :: um_integer32 = SELECTED_INT_KIND(9)

! Missing data, real and integer
REAL(KIND=um_real64), PARAMETER       :: um_rmdi     = -32768.0*32768.0
INTEGER(KIND=um_integer64), PARAMETER :: um_imdi     = -32768
INTEGER(KIND=um_integer32), PARAMETER :: um_imdi_32  = um_imdi


! Meaningful parameter names for real constants header
! East-West   grid spacing in degrees
INTEGER(KIND=int64), PARAMETER, PUBLIC :: rh_deltaEW         = 1
! North-South grid spacing in degrees
INTEGER(KIND=int64), PARAMETER, PUBLIC :: rh_deltaNS         = 2
! Latitude  of first p point in degrees
INTEGER(KIND=int64), PARAMETER, PUBLIC :: rh_baselat         = 3
! Longitude of first p point in degrees
INTEGER(KIND=int64), PARAMETER, PUBLIC :: rh_baselong        = 4
! Latitude  of rotated N pole in degrees
INTEGER(KIND=int64), PARAMETER, PUBLIC :: rh_polelat         = 5
! Longitude of rotated N pole in degrees
INTEGER(KIND=int64), PARAMETER, PUBLIC :: rh_polelong        = 6
! Height of top theta level (m)
INTEGER(KIND=int64), PARAMETER, PUBLIC :: rh_model_top       =16

! Meaningful parameter names for integer constants header
! No. of points E-W
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ih_row_length      = 6
! No. of points N-S
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ih_rows            = 7
! No. of model levels (0=surface)
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ih_model_levels    = 8
! No. of model levels with moisture
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ih_wet_levels      = 9
! No. of deep soil temperature levels
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ih_soilT_levels    = 10
! No. of cloud levels
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ih_cloud_levels    = 11
! No. of tracer levels
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ih_tracer_levels   = 12
! No. of boundary layer levels
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ih_boundary_levels = 13
! No. of field types
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ih_N_types         = 15
! Height generation method
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ih_height_gen      = 17
! First rho level at which height is constant
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ih_1_c_rho_level   = 24
! No. of land points
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ih_land_points     = 25
! No. of ozone levels
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ih_ozone_levels    = 26
! No. of deep soil moisture levels
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ih_soilQ_levels    = 28
! Number of convective cloud levels
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ih_convect_levels  = 34

! Meaningful parameter names for level dependent constants
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ldc_eta_theta  = 1
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ldc_eta_rho    = 2
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ldc_Zsea_theta = 5
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ldc_C_theta    = 6
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ldc_Zsea_rho   = 7
INTEGER(KIND=int64), PARAMETER, PUBLIC :: ldc_C_rho      = 8

END MODULE lfricinp_um_parameters_mod
