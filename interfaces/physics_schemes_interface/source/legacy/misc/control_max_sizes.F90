! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc

! CMAXSIZE maximum sizes for dimensioning arrays
! of model constants whose sizes are configuration dependent. This
! allows constants to be read in from a namelist file and maintain
! the flexibility of dynamic allocation for primary variables. The
! maximum sizes should agree with the maximum sizes implicit in the
! front-end User Interface.

module Control_Max_Sizes

use Atmos_Max_Sizes, only: model_levels_max
implicit none

! Max levels in boundary layer
integer, parameter :: max_bl_levels = model_levels_max

! Max no. of levels for pvort output
integer, parameter :: max_req_thpv_levs = model_levels_max

! Max no. 1-2-1 rows in polar filter
integer, parameter ::  max_121_rows =  32
! 0 is used for horizontal diffusion pointer

! Max no. of levels (from top) to apply upper level diffusion
integer, parameter ::  max_updiff_levels = 10

! Max size of any sponge zones
integer, parameter ::  max_sponge_width = 10

! Max size of look-up tables for searches
integer, parameter ::  max_look = 2048

! Max no. of atmos interface areas
integer, parameter :: max_n_intf_a =  8

! Max no. of points in LBC
integer, parameter :: max_intf_lbcrow_length = 5000
integer, parameter :: max_intf_lbcrows = 5000

! Max no. of atmos interface levels
integer, parameter :: max_intf_levels = model_levels_max

! Maximum number of physics segments
integer, parameter :: max_no_of_segs = 200

end module Control_Max_Sizes
