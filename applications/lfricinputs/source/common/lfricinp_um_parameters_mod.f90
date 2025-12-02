! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_um_parameters_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only: int64

implicit none

! Filename length
integer, public, parameter :: fnamelen = 512

! Message length
integer, parameter :: msglen = 512

! UM definitions of real and integer in fieldsfiles
integer, public, parameter :: um_real64  = selected_real_kind(15,307)
integer, public, parameter :: um_integer64 = selected_int_kind(15)
integer, public, parameter :: um_real32  = selected_real_kind(6,37)
integer, public, parameter :: um_integer32 = selected_int_kind(9)

! Missing data, real and integer
real(kind=um_real64), parameter       :: um_rmdi     = -32768.0_um_real64*32768.0_um_real64
integer(kind=um_integer64), parameter :: um_imdi     = -32768
integer(kind=um_integer32), parameter :: um_imdi_32  = um_imdi


! Meaningful parameter names for real constants header
! East-West   grid spacing in degrees
integer(kind=int64), parameter, public :: rh_deltaEW         = 1
! North-South grid spacing in degrees
integer(kind=int64), parameter, public :: rh_deltaNS         = 2
! Latitude  of first p point in degrees
integer(kind=int64), parameter, public :: rh_baselat         = 3
! Longitude of first p point in degrees
integer(kind=int64), parameter, public :: rh_baselong        = 4
! Latitude  of rotated N pole in degrees
integer(kind=int64), parameter, public :: rh_polelat         = 5
! Longitude of rotated N pole in degrees
integer(kind=int64), parameter, public :: rh_polelong        = 6
! Height of top theta level (m)
integer(kind=int64), parameter, public :: rh_model_top       =16

! Meaningful parameter names for integer constants header
! No. of points E-W
integer(kind=int64), parameter, public :: ih_row_length      = 6
! No. of points N-S
integer(kind=int64), parameter, public :: ih_rows            = 7
! No. of model levels (0=surface)
integer(kind=int64), parameter, public :: ih_model_levels    = 8
! No. of model levels with moisture
integer(kind=int64), parameter, public :: ih_wet_levels      = 9
! No. of deep soil temperature levels
integer(kind=int64), parameter, public :: ih_soilT_levels    = 10
! No. of cloud levels
integer(kind=int64), parameter, public :: ih_cloud_levels    = 11
! No. of tracer levels
integer(kind=int64), parameter, public :: ih_tracer_levels   = 12
! No. of boundary layer levels
integer(kind=int64), parameter, public :: ih_boundary_levels = 13
! No. of field types
integer(kind=int64), parameter, public :: ih_N_types         = 15
! Height generation method
integer(kind=int64), parameter, public :: ih_height_gen      = 17
! First rho level at which height is constant
integer(kind=int64), parameter, public :: ih_1_c_rho_level   = 24
! No. of land points
integer(kind=int64), parameter, public :: ih_land_points     = 25
! No. of ozone levels
integer(kind=int64), parameter, public :: ih_ozone_levels    = 26
! No. of deep soil moisture levels
integer(kind=int64), parameter, public :: ih_soilQ_levels    = 28
! Number of convective cloud levels
integer(kind=int64), parameter, public :: ih_convect_levels  = 34

! Meaningful parameter names for level dependent constants
integer(kind=int64), parameter, public :: ldc_eta_theta  = 1
integer(kind=int64), parameter, public :: ldc_eta_rho    = 2
integer(kind=int64), parameter, public :: ldc_Zsea_theta = 5
integer(kind=int64), parameter, public :: ldc_C_theta    = 6
integer(kind=int64), parameter, public :: ldc_Zsea_rho   = 7
integer(kind=int64), parameter, public :: ldc_C_rho      = 8

end module lfricinp_um_parameters_mod
