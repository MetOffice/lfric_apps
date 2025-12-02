! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module profiles_mod

! Purpose:
!   Contains profiles and associated variables for idealised configurations.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Idealised
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP3 programming standards.

use missing_data_mod, only: rmdi, imdi

implicit none

save

! Parameters
! ======================
! Types of vertical coordinate
integer, parameter :: eta_coord     = 1
integer, parameter :: z_coord       = 2
integer, parameter :: p_coord       = 3

! Types of fields
integer, parameter :: theta_dry     = 10
integer, parameter :: temperature   = 11
integer, parameter :: dtheta_dz     = 12
integer, parameter :: Brunt_Vaisala = 13
integer, parameter :: mixing_ratio  = 20
integer, parameter :: rel_humidity  = 21

! Maximum size of arrays that can be read through namelist
! currently limited to reduced namelist size
integer, parameter :: num_prof_max  = 100  ! max vertical profile values
integer, parameter :: num_time_max  = 100  ! max times
integer, parameter :: num_data_max  = 100  ! may need to be increased as should
                                           ! really be num_time_max*num_prof_max
                                           ! which is 10000, a very big number

! Profiles
! ========

! Derived type for carrying time-varying profile data
type :: varying_profile
  integer           :: n_times
  integer           :: n_heights
  integer           :: coord_type
  integer           :: field_type
  real              :: timescale
  real, allocatable :: tsec(:)
  real, allocatable :: height(:)
  real, allocatable :: vprof(:,:)
end type varying_profile

! Initial Data Profiles
! =====================
type(varying_profile) :: theta_init ! dry potential temperature
type(varying_profile) :: mv_init    ! vapour mixing ratio
type(varying_profile) :: u_init     ! xi1-component of wind
type(varying_profile) :: v_init     ! xi2-component of wind
type(varying_profile) :: o3_init    ! ozone

type(varying_profile) :: tracer_init ! free tracer initial profile - recon only

! Forcing Profiles
! ================

! Geostrophic forcing profiles
type(varying_profile) :: u_geostrophic
type(varying_profile) :: v_geostrophic

! Newton relaxation profiles
type(varying_profile) :: theta_relax ! dry potential temperature
type(varying_profile) :: mv_relax    ! vapour mixing ratio
type(varying_profile) :: u_relax     ! xi1-component of wind
type(varying_profile) :: v_relax     ! xi2-component of wind


! Increment profiles
type(varying_profile) :: theta_inc
type(varying_profile) :: mv_inc
type(varying_profile) :: u_inc
type(varying_profile) :: v_inc

! Subsidence
type(varying_profile) :: w_subs

! Data Read from Namelist
! =======================

! Relaxation
integer :: num_theta_relax_heights = imdi
integer :: num_theta_relax_times   = imdi
integer :: theta_relax_field_type  = imdi
integer :: num_mv_relax_heights    = imdi
integer :: num_mv_relax_times      = imdi
integer :: num_uv_relax_heights    = imdi
integer :: num_uv_relax_times      = imdi

! Timescale (in seconds) on which fields are relaxed
real    :: theta_relax_timescale   = rmdi
real    :: mv_relax_timescale      = rmdi
real    :: uv_relax_timescale      = rmdi


real    :: theta_relax_height(num_prof_max) = rmdi
real    :: theta_relax_time(num_time_max)   = rmdi
real    :: theta_relax_data(num_data_max)   = rmdi
real    :: mv_relax_height(num_prof_max)    = rmdi
real    :: mv_relax_time(num_time_max)      = rmdi
real    :: mv_relax_data(num_data_max)      = rmdi
real    :: uv_relax_height(num_prof_max)    = rmdi
real    :: uv_relax_time(num_time_max)      = rmdi
real    :: u_relax_data(num_data_max)       = rmdi
real    :: v_relax_data(num_data_max)       = rmdi

! Increment
integer :: num_theta_inc_times   = imdi
integer :: num_theta_inc_heights = imdi
integer :: theta_inc_field_type  = imdi
integer :: num_mv_inc_times      = imdi
integer :: num_mv_inc_heights    = imdi
integer :: num_uv_inc_times      = imdi
integer :: num_uv_inc_heights    = imdi
integer :: num_w_force_times     = imdi
integer :: num_w_force_heights   = imdi

real    :: theta_inc_height(num_prof_max) = rmdi
real    :: theta_inc_time(num_time_max)   = rmdi
real    :: theta_inc_data(num_data_max)   = rmdi
real    :: mv_inc_height(num_prof_max)    = rmdi
real    :: mv_inc_time(num_time_max)      = rmdi
real    :: mv_inc_data(num_data_max)      = rmdi
real    :: uv_inc_height(num_prof_max)    = rmdi
real    :: uv_inc_time(num_time_max)      = rmdi
real    :: u_inc_data(num_data_max)       = rmdi
real    :: v_inc_data(num_data_max)       = rmdi
real    :: w_force_height(num_prof_max)   = rmdi
real    :: w_force_time(num_time_max)     = rmdi
real    :: w_force_data(num_data_max)     = rmdi

! For time and vertical height varying geostrophic forcing
integer :: num_uv_geo_times = imdi
integer :: num_uv_geo_heights = imdi

real    :: uv_geo_height(num_prof_max) = rmdi
real    :: uv_geo_time(num_time_max)   = rmdi
real    :: u_geo_data(num_data_max)    = rmdi
real    :: v_geo_data(num_data_max)    = rmdi


! Hydrostatic Reference Pressure
! ==============================

! Flag to indicate that at least one profile is specified against pressure
logical                   :: l_p_profile = .false.
real, allocatable, target :: p_prof_theta(:)
real, allocatable, target :: p_prof_rho(:)

end module profiles_mod
