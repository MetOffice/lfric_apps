! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

#if !defined(LFRIC)
module standalone_test_mod

implicit none

! Size of diagnostics super-array
! (declared in module so that all contained subroutines can
!  see it)
integer, parameter :: n_diags_max_3d = 120
integer, parameter :: n_diags_max_5d = 10
integer, parameter :: n_diags_max_2d = 8
integer, parameter :: n_diags_max_4d = 2

! Unit number to use for writing the output diagnostics file
integer, parameter :: unit_diag = 10

contains


!----------------------------------------------------------------
! Standalone single-timestep test-run of the comorph convection
! scheme.
!----------------------------------------------------------------
! Sets up idealised initial profiles, sets diagnostic
! requests, calls comorph_ctl, and then writes the output
! diagnostics to a text file.
! This subroutine needs to be called by a "shell" program which
! sets the array sizes in the CoMorph constants module.
subroutine standalone_test()

use set_test_profiles_mod, only: set_test_profiles
use comorph_ctl_mod, only: comorph_ctl
use fields_type_mod, only: fields_type
use grid_type_mod, only: grid_type
use turb_type_mod, only: turb_type, n_turb
use cloudfracs_type_mod, only: cloudfracs_type, n_cloudfracs
use comorph_diags_type_mod, only: comorph_diags_type
use diag_type_mod, only: diag_list_type
use subregion_mod, only: n_regions

use comorph_constants_mod, only: real_hmprec,                                  &
                                 nx_full, ny_full, k_top_conv, k_top_init,     &
                                 n_tracers, n_updraft_types,                   &
                                 n_conv_layers_diag,                           &
                                 l_cv_cloudfrac, l_cv_graup,                   &
                                 l_turb_par_gen, l_par_core,                   &
                                 l_calc_cape, l_calc_mfw_cape, l_calc_ccb_cct

implicit none

! Halo size on pressure
integer, parameter :: hp = 1
! Halo size on q fields
integer, parameter :: hq = 1

! Model grid fields
real(kind=real_hmprec), target :: height_full                                  &
        ( nx_full, ny_full, 0:k_top_conv )
real(kind=real_hmprec), target :: height_half                                  &
        ( nx_full, ny_full, 1:k_top_conv+1 )
real(kind=real_hmprec), target :: pressure_full                                &
        ( 1-hp:nx_full+hp, 1-hp:ny_full+hp, 0:k_top_conv )
real(kind=real_hmprec), target :: pressure_half                                &
        ( 1-hp:nx_full+hp, 1-hp:ny_full+hp, 1:k_top_conv+1 )
real(kind=real_hmprec), target :: rho_dry                                      &
        ( nx_full, ny_full, 1:k_top_conv )
real(kind=real_hmprec), target :: r_surf(nx_full,ny_full)

! Turbulence fields
real(kind=real_hmprec), target :: turb_super                                   &
        ( nx_full, ny_full, 1:k_top_init+1, n_turb )
real(kind=real_hmprec), target :: turb_len                                     &
        ( nx_full, ny_full, 1:k_top_init )
real(kind=real_hmprec), target :: par_radius_amp                               &
        ( nx_full, ny_full )
real(kind=real_hmprec), target :: z_bl_top                                     &
        ( nx_full, ny_full )

! Rain / graupel fraction
real(kind=real_hmprec), target :: precfrac                                     &
        ( nx_full, ny_full, 1:k_top_conv )

! Convective cloud fraction and water content
real(kind=real_hmprec), target :: conv_cloud                                   &
        ( nx_full, ny_full, 1:k_top_conv, 3 )

! Environment profiles
! (declared here using super-arrays, but split into a handful
!  of separate arrays so that different halos can be applied
!  to each, to check that the handling of halos in CoMorph
!  works OK)
real(kind=real_hmprec), target :: winds_super                                  &
        ( nx_full, ny_full, 1:k_top_conv, 6 )
real(kind=real_hmprec), target :: temperature                                  &
        ( nx_full, ny_full, 0:k_top_conv, 2 )
real(kind=real_hmprec), target :: q_super                                      &
        ( 1-hq:nx_full+hq, 1-hq:ny_full+hq, 0:k_top_conv, 12 )
real(kind=real_hmprec), target :: cf_super                                     &
        ( nx_full, ny_full, 1:k_top_conv, 2*n_cloudfracs )
real(kind=real_hmprec), target :: tracer                                       &
        ( nx_full, ny_full, 1:k_top_conv, 2*n_tracers )

! Space for diagnostics
real(kind=real_hmprec) :: diags_super_2d                                       &
                          ( nx_full, ny_full,                                  &
                            n_diags_max_2d )
real(kind=real_hmprec) :: diags_super_3d                                       &
                          ( nx_full, ny_full, 1:k_top_conv,                    &
                            n_diags_max_3d )
real(kind=real_hmprec) :: diags_super_4d                                       &
                          ( nx_full, ny_full,                                  &
                            n_conv_layers_diag, n_updraft_types,               &
                            n_diags_max_4d )
real(kind=real_hmprec) :: diags_super_5d                                       &
                          ( nx_full, ny_full, 1:k_top_conv,                    &
                            n_conv_layers_diag, n_updraft_types,               &
                            n_diags_max_5d )

! Filename for output
character(len=100) :: filename

! Inputs to comorph
type(fields_type) :: fields_n
type(fields_type) :: fields_np1
type(grid_type) :: grid
type(turb_type) :: turb
type(cloudfracs_type) :: cloudfracs
type(comorph_diags_type) :: comorph_diags
logical :: l_tracer
integer :: n_segments

! List of pointers to active diagnostic structures,
! to access meta-data.
type(diag_list_type) :: diags_list_2d(n_diags_max_2d)
type(diag_list_type) :: diags_list_3d(n_diags_max_3d)
type(diag_list_type) :: diags_list_4d(n_diags_max_4d)
type(diag_list_type) :: diags_list_5d(n_diags_max_5d)

character(len=2) :: num

integer :: lb(5), ub(5)

! Loop counters
integer :: i_field, i_region, i_diag,                                          &
           n_diags_2d, n_diags_3d, n_diags_4d, n_diags_5d, n


! Set flag for whether there are any tracers
l_tracer = n_tracers > 0

! Assign pointers for fields input to comorph_ctl

grid % height_full => height_full
grid % height_half => height_half
grid % pressure_full => pressure_full
grid % pressure_half => pressure_half
grid % rho_dry => rho_dry
grid % r_surf => r_surf

if ( l_turb_par_gen ) then
  turb % w_var    => turb_super(:,:,:,1)
  turb % f_templ  => turb_super(:,:,:,2)
  turb % f_q_tot  => turb_super(:,:,:,3)
  turb % f_wind_u => turb_super(:,:,:,4)
  turb % f_wind_v => turb_super(:,:,:,5)
end if
turb % lengthscale => turb_len
turb % par_radius_amp => par_radius_amp
turb % z_bl_top => z_bl_top

if ( .not. l_cv_cloudfrac ) then
  lb(1:4) = lbound(cf_super)
  ub(1:4) = ubound(cf_super)
  cloudfracs % frac_liq( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )               &
                                                         => cf_super(:,:,:,1)
  cloudfracs % frac_ice( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )               &
                                                         => cf_super(:,:,:,2)
  cloudfracs % frac_bulk( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )              &
                                                         => cf_super(:,:,:,3)
end if
cloudfracs % frac_precip => precfrac

cloudfracs % frac_liq_conv  => conv_cloud(:,:,:,1)
cloudfracs % frac_bulk_conv => conv_cloud(:,:,:,2)
cloudfracs % q_cl_conv      => conv_cloud(:,:,:,3)

lb(1:4) = lbound(winds_super)
ub(1:4) = ubound(winds_super)
fields_n % wind_u( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                     &
                                                     => winds_super(:,:,:,1)
fields_n % wind_v( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                     &
                                                     => winds_super(:,:,:,2)
fields_n % wind_w( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                     &
                                                     => winds_super(:,:,:,3)
fields_np1 % wind_u( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                   &
                                                       => winds_super(:,:,:,4)
fields_np1 % wind_v( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                   &
                                                       => winds_super(:,:,:,5)
fields_np1 % wind_w( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                   &
                                                       => winds_super(:,:,:,6)

lb(1:4) = lbound(temperature)
ub(1:4) = ubound(temperature)
fields_n % temperature( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                &
                                                     => temperature(:,:,:,1)
fields_np1 % temperature( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )              &
                                                       => temperature(:,:,:,2)

lb(1:4) = lbound(q_super)
ub(1:4) = ubound(q_super)
fields_n % q_vap   ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                   &
                                                         => q_super(:,:,:,1)
fields_n % q_cl    ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                   &
                                                         => q_super(:,:,:,2)
fields_n % q_rain  ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                   &
                                                         => q_super(:,:,:,3)
fields_n % q_cf    ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                   &
                                                         => q_super(:,:,:,4)
fields_n % q_snow  ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                   &
                                                         => q_super(:,:,:,5)
fields_n % q_graup ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                   &
                                                         => q_super(:,:,:,6)
fields_np1 % q_vap   ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                 &
                                                           => q_super(:,:,:,7)
fields_np1 % q_cl    ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                 &
                                                           => q_super(:,:,:,8)
fields_np1 % q_rain  ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                 &
                                                           => q_super(:,:,:,9)
fields_np1 % q_cf    ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                 &
                                                           => q_super(:,:,:,10)
fields_np1 % q_snow  ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                 &
                                                           => q_super(:,:,:,11)
fields_np1 % q_graup ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                 &
                                                           => q_super(:,:,:,12)

lb(1:4) = lbound(cf_super)
ub(1:4) = ubound(cf_super)
fields_n % cf_liq( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                     &
                                                        => cf_super(:,:,:,1)
fields_n % cf_ice( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                     &
                                                        => cf_super(:,:,:,2)
fields_n % cf_bulk( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                    &
                                                        => cf_super(:,:,:,3)
fields_np1 % cf_liq( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                   &
                                                          => cf_super(:,:,:,4)
fields_np1 % cf_ice( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                   &
                                                          => cf_super(:,:,:,5)
fields_np1 % cf_bulk( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )                  &
                                                          => cf_super(:,:,:,6)

if ( l_tracer ) then
  allocate( fields_n % tracers(n_tracers) )
  allocate( fields_np1 % tracers(n_tracers) )
  lb(1:4) = lbound(tracer)
  ub(1:4) = ubound(tracer)
  do i_field = 1, n_tracers
    fields_n % tracers(i_field)%pt( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )    &
                                          => tracer(:,:,:,i_field)
    fields_np1 % tracers(i_field)%pt( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )  &
                                            => tracer(:,:,:,n_tracers+i_field)
  end do
end if
! (Note: CoMorph does not use tracers in fields_n, but
!  we add them here so that we can reset the tracers when
!  calling CoMorph multiple times to check reproducibility
!  under different decompositions).


! Setup idealised heights, pressures and initial profiles
call set_test_profiles( l_tracer, grid, turb, cloudfracs,                      &
                        fields_n )


! Set CoMorph 2-D diagnostic requests...
i_diag = 0

! CAPE and mass-flux-weighted CAPE
call assign_diag_2d( diags_super_2d, diags_list_2d, i_diag,                    &
                     comorph_diags % updraft_diags_2d % cape )
call assign_diag_2d( diags_super_2d, diags_list_2d, i_diag,                    &
                     comorph_diags % updraft_diags_2d % mfw_cape )
! Updraft cloud-top and cloud-base properties
call assign_diag_2d( diags_super_2d, diags_list_2d, i_diag,                    &
                     comorph_diags % updraft_diags_2d % term_height )
call assign_diag_2d( diags_super_2d, diags_list_2d, i_diag,                    &
                     comorph_diags % updraft_diags_2d % ccb_height )
call assign_diag_2d( diags_super_2d, diags_list_2d, i_diag,                    &
                     comorph_diags % updraft_diags_2d % ccb_area )
call assign_diag_2d( diags_super_2d, diags_list_2d, i_diag,                    &
                     comorph_diags % updraft_diags_2d % ccb_massflux_d )
call assign_diag_2d( diags_super_2d, diags_list_2d, i_diag,                    &
                     comorph_diags % updraft_diags_2d % ccb_mean_tv_ex )
call assign_diag_2d( diags_super_2d, diags_list_2d, i_diag,                    &
                     comorph_diags % updraft_diags_2d % ccb_core_tv_ex )
! Make sure flags are set to calculate these inside comorph
l_calc_cape = .true.
l_calc_mfw_cape = .true.
l_calc_ccb_cct = .true.

! Count number of diagnostics assigned
n_diags_2d = i_diag


! Set CoMorph 3-D diagnostic requests...
i_diag = 0

! Increments to primary fields
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
                     comorph_diags % fields_incr % wind_u )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
                     comorph_diags % fields_incr % wind_v )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
                     comorph_diags % fields_incr % wind_w )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
                     comorph_diags % fields_incr % temperature )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
                     comorph_diags % fields_incr % q_vap )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
                     comorph_diags % fields_incr % q_cl )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
                     comorph_diags % fields_incr % q_rain )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
                     comorph_diags % fields_incr % q_cf )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
                     comorph_diags % fields_incr % q_graup )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
                     comorph_diags % fields_incr % cf_liq )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
                     comorph_diags % fields_incr % cf_ice )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
                     comorph_diags % fields_incr % cf_bulk )

! Miscellaneous
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
                     comorph_diags % pressure_incr_env )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
                     comorph_diags % layer_mass )

! Updraft properties
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % par % massflux_d )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % par % radius )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % par % mean % wind_u )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % par % mean % wind_v )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % par % mean % wind_w )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % par % mean % temperature )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % par % mean % q_vap )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % par % mean % q_cl )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % par % mean % q_rain )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % par % mean % q_cf )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % par % mean % q_graup )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % par % mean % virt_temp_excess )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % par % mean % rel_hum_liq )
if ( l_par_core ) then
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
         comorph_diags % updraft % par % core % temperature )
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
         comorph_diags % updraft % par % core % q_vap )
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
         comorph_diags % updraft % par % core % q_cl )
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
         comorph_diags % updraft % par % core % q_rain )
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
         comorph_diags % updraft % par % core % q_cf )
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
         comorph_diags % updraft % par % core % q_graup )
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
         comorph_diags % updraft % par % core % virt_temp_excess)
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
         comorph_diags % updraft % par % core % rel_hum_liq )
end if

! Downdraft properties
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % par % massflux_d )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % par % radius )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % par % mean % wind_u )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % par % mean % wind_v )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % par % mean % wind_w )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % par % mean % temperature )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % par % mean % q_vap )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % par % mean % virt_temp_excess )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % par % mean % rel_hum_liq )

if ( l_par_core ) then
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
         comorph_diags % dndraft % par % core % temperature )
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
         comorph_diags % dndraft % par % core % q_vap )
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
         comorph_diags % dndraft % par % core % virt_temp_excess)
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
         comorph_diags % dndraft % par % core % rel_hum_liq )
end if

! Initiation mass-sources
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % gen % massflux_d )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % gen % massflux_d )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % gen % mean % virt_temp_excess )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % gen % mean % virt_temp_excess )
if ( l_par_core ) then
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
       comorph_diags % updraft % gen % core % virt_temp_excess )
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
       comorph_diags % dndraft % gen % core % virt_temp_excess )
end if

! Updraft entrained and detrained mass and total-water
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % plume_model % ent_mass_d )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % plume_model % det_mass_d )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % plume_model % ent_fields%q_tot )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % plume_model % det_fields%q_tot )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % plume_model % ent_fields % virt_temp_excess )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % plume_model % det_fields % virt_temp_excess )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % plume_model % ent_fields % rel_hum_liq )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % plume_model % det_fields % rel_hum_liq )


! Same for downdrafts
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % plume_model % ent_mass_d )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % plume_model % det_mass_d )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % plume_model % ent_fields%q_tot )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % plume_model % det_fields%q_tot )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % plume_model % ent_fields % virt_temp_excess )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % plume_model % det_fields % virt_temp_excess )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % plume_model % ent_fields % rel_hum_liq )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % dndraft % plume_model % det_fields % rel_hum_liq )


! Updraft precip fall-speeds
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % plume_model % moist_proc                      &
                     % diags_cl % fall_speed )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % plume_model % moist_proc                      &
                     % diags_rain % fall_speed )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % plume_model % moist_proc                      &
                     % diags_cf % fall_speed )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % updraft % plume_model % moist_proc                      &
                     % diags_graup % fall_speed )

! Environment sub-region properties
do i_region = 1, n_regions
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
       comorph_diags % genesis_diags % subregion_diags(i_region) % frac )
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
       comorph_diags % genesis_diags % subregion_diags(i_region) % fields      &
                     % temperature )
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
       comorph_diags % genesis_diags % subregion_diags(i_region) % fields      &
                     % q_vap )
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
       comorph_diags % genesis_diags % subregion_diags(i_region) % fields      &
                     % q_cl )
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
       comorph_diags % genesis_diags % subregion_diags(i_region) % fields      &
                     % virt_temp_excess )
  call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                  &
       comorph_diags % genesis_diags % subregion_diags(i_region) % fields      &
                     % rel_hum_liq )
end do

! Turbulence-based parcel perturbations and radius
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % turb_fields_pert % wind_u )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % turb_fields_pert % wind_v )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % turb_fields_pert % wind_w )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % turb_fields_pert % temperature )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % turb_fields_pert % q_vap )
call assign_diag_3d( diags_super_3d, diags_list_3d, i_diag,                    &
       comorph_diags % turb_radius )

! Count number of diagnostics assigned
n_diags_3d = i_diag


! Set CoMorph 4-D diagnostic requests...
i_diag = 0

! CAPE and mass-flux-weighted CAPE
call assign_diag_4d( n_updraft_types, diags_super_4d, diags_list_4d, i_diag,   &
                     comorph_diags % updraft_diags_2d % cape )
call assign_diag_4d( n_updraft_types, diags_super_4d, diags_list_4d, i_diag,   &
                     comorph_diags % updraft_diags_2d % mfw_cape )

! Count number of diagnostics assigned
n_diags_4d = i_diag


! Set CoMorph 5-D diagnostic requests...
i_diag = 0

! Updraft properties
call assign_diag_5d( n_updraft_types, diags_super_5d, diags_list_5d, i_diag,   &
       comorph_diags % updraft % par % massflux_d )
call assign_diag_5d( n_updraft_types, diags_super_5d, diags_list_5d, i_diag,   &
       comorph_diags % updraft % par % radius )
call assign_diag_5d( n_updraft_types, diags_super_5d, diags_list_5d, i_diag,   &
       comorph_diags % updraft % par % mean % virt_temp_excess )
call assign_diag_5d( n_updraft_types, diags_super_5d, diags_list_5d, i_diag,   &
       comorph_diags % updraft % par % mean % rel_hum_liq )
if ( l_par_core ) then
  call assign_diag_5d( n_updraft_types, diags_super_5d, diags_list_5d, i_diag, &
         comorph_diags % updraft % par % core % virt_temp_excess)
  call assign_diag_5d( n_updraft_types, diags_super_5d, diags_list_5d, i_diag, &
         comorph_diags % updraft % par % core % rel_hum_liq )
end if
call assign_diag_5d( n_updraft_types, diags_super_5d, diags_list_5d, i_diag,   &
       comorph_diags % updraft % plume_model % core_mean_ratio )

! Count number of diagnostics assigned
n_diags_5d = i_diag


! Loop to call comorph and output diagnostics multiple times
! with different decompositions, to check we get the same answer
! every time...
do n = 1, 4

  ! Set inputs differently on different calls;
  ! simulates 2 timesteps, each with 2 calls to convection,
  ! and different segment size each time.
  if ( n==1 ) then
    n_segments = 1                ! All columns at once
    l_tracer = .false.
  else if ( n==2 ) then
    n_segments = nx_full*ny_full  ! One column at a time
    l_tracer = n_tracers > 0
  else if ( n==3 ) then
    n_segments = 3                ! 3 segments
    l_tracer = .false.
  else if ( n==4 ) then
    n_segments = 1                ! All columns at once
    l_tracer = n_tracers > 0
  end if

  ! Initialise the _np1 fields equal to the _n fields
  ! (use the super-arrays that fields structures point at)
  winds_super(:,:,:,4:6) = winds_super(:,:,:,1:3)
  temperature(:,:,:,2)   = temperature(:,:,:,1)
  q_super(:,:,:,7:12)    = q_super(:,:,:,1:6)
  cf_super(:,:,:,4:6)    = cf_super(:,:,:,1:3)
  tracer(:,:,:,n_tracers+1:2*n_tracers)                                        &
                         = tracer(:,:,:,1:n_tracers)


  ! Call comorph!
  call comorph_ctl( l_tracer, n_segments,                                      &
                    grid, turb, cloudfracs, fields_n, fields_np1,              &
                    comorph_diags )


  ! Reset tracer switch
  l_tracer = n_tracers > 0


  ! Print diagnostic outputs to a python module file...

  write(filename,"(A,I1,A)") "standalone_test_out_", n, ".py"

  open(unit=unit_diag,file=filename)

  write(unit_diag,"(A)") "import numpy as np"
  write(unit_diag,"(A)") ""

  ! Write out CoMorph diagnostics
  do i_diag = 1, n_diags_2d
    lb(1:2) = [1,1]
    ub(1:2) = [nx_full,ny_full]
    call write_diag_2d( diags_list_2d(i_diag)%pt % diag_name, lb(1:2), ub(1:2),&
                        diags_super_2d(:,:,i_diag) )
  end do
  do i_diag = 1, n_diags_3d
    lb(1:3) = [1,1,1]
    ub(1:3) = [nx_full,ny_full,k_top_conv]
    call write_diag_3d( diags_list_3d(i_diag)%pt % diag_name, lb(1:3), ub(1:3),&
                        diags_super_3d(:,:,:,i_diag) )
  end do
  do i_diag = 1, n_diags_4d
    lb(1:4) = [1,1,1,1]
    ub(1:4) = [nx_full,ny_full,n_conv_layers_diag,n_updraft_types]
    call write_diag_4d( diags_list_4d(i_diag)%pt % diag_name, lb(1:4), ub(1:4),&
                        n_updraft_types, diags_super_4d(:,:,:,:,i_diag) )
  end do
  do i_diag = 1, n_diags_5d
    lb(1:5) = [1,1,1,1,1]
    ub(1:5) = [nx_full,ny_full,k_top_conv,n_conv_layers_diag,n_updraft_types]
    call write_diag_5d( diags_list_5d(i_diag)%pt % diag_name, lb(1:5), ub(1:5),&
                        n_updraft_types, diags_super_5d(:,:,:,:,:,i_diag) )
  end do

  ! Write out heights and pressures
  lb(1:3) = lbound(grid % height_full)
  ub(1:3) = ubound(grid % height_full)
  call write_diag_3d( "height_full", lb(1:3), ub(1:3), grid % height_full )
  lb(1:3) = lbound(grid % pressure_full)
  ub(1:3) = ubound(grid % pressure_full)
  call write_diag_3d( "pressure_full", lb(1:3), ub(1:3), grid % pressure_full )
  lb(1:3) = lbound(grid % height_half)
  ub(1:3) = ubound(grid % height_half)
  call write_diag_3d( "height_half", lb(1:3), ub(1:3), grid % height_half )
  lb(1:3) = lbound(grid % pressure_half)
  ub(1:3) = ubound(grid % pressure_half)
  call write_diag_3d( "pressure_half", lb(1:3), ub(1:3), grid % pressure_half )

  ! Write out initial fields
  lb(1:3) = lbound(fields_n % wind_u)
  ub(1:3) = ubound(fields_n % wind_u)
  call write_diag_3d( "wind_u_n", lb(1:3), ub(1:3), fields_n % wind_u )
  call write_diag_3d( "wind_v_n", lb(1:3), ub(1:3), fields_n % wind_v )
  call write_diag_3d( "wind_w_n", lb(1:3), ub(1:3), fields_n % wind_w )
  lb(1:3) = lbound(fields_n % temperature)
  ub(1:3) = ubound(fields_n % temperature)
  call write_diag_3d( "temperature_n", lb(1:3), ub(1:3),                       &
                      fields_n % temperature )
  lb(1:3) = lbound(fields_n % q_vap)
  ub(1:3) = ubound(fields_n % q_vap)
  call write_diag_3d( "q_vap_n", lb(1:3), ub(1:3), fields_n % q_vap )
  call write_diag_3d( "q_cl_n", lb(1:3), ub(1:3), fields_n % q_cl )
  call write_diag_3d( "q_rain_n", lb(1:3), ub(1:3), fields_n % q_rain )
  call write_diag_3d( "q_cf_n", lb(1:3), ub(1:3), fields_n % q_cf )
  if ( l_cv_graup )                                                            &
    call write_diag_3d( "q_graup_n", lb(1:3), ub(1:3), fields_n % q_graup )
  if ( l_cv_cloudfrac ) then
    lb(1:3) = lbound(fields_n % cf_liq)
    ub(1:3) = ubound(fields_n % cf_liq)
    call write_diag_3d( "cf_liq_n", lb(1:3), ub(1:3), fields_n % cf_liq )
    call write_diag_3d( "cf_ice_n", lb(1:3), ub(1:3), fields_n % cf_ice )
    call write_diag_3d( "cf_bulk_n", lb(1:3), ub(1:3), fields_n % cf_bulk )
  end if
  if ( l_tracer ) then
    lb(1:3) = lbound(fields_n % tracers(1)%pt)
    ub(1:3) = ubound(fields_n % tracers(1)%pt)
    do i_field = 1, n_tracers
      write(num,"(I2)") i_field
      call write_diag_3d( "tracer_" // trim(adjustl(num)) // "_n",             &
                          lb(1:3), ub(1:3), fields_n % tracers(i_field)%pt )
    end do
  end if

  ! Write out final fields
  lb(1:3) = lbound(fields_np1 % wind_u)
  ub(1:3) = ubound(fields_np1 % wind_u)
  call write_diag_3d( "wind_u_np1", lb(1:3), ub(1:3), fields_np1 % wind_u )
  call write_diag_3d( "wind_v_np1", lb(1:3), ub(1:3), fields_np1 % wind_v )
  call write_diag_3d( "wind_w_np1", lb(1:3), ub(1:3), fields_np1 % wind_w )
  lb(1:3) = lbound(fields_np1 % temperature)
  ub(1:3) = ubound(fields_np1 % temperature)
  call write_diag_3d( "temperature_np1", lb(1:3), ub(1:3),                     &
                      fields_np1 % temperature)
  lb(1:3) = lbound(fields_np1 % q_vap)
  ub(1:3) = ubound(fields_np1 % q_vap)
  call write_diag_3d( "q_vap_np1", lb(1:3), ub(1:3), fields_np1 % q_vap )
  call write_diag_3d( "q_cl_np1", lb(1:3), ub(1:3), fields_np1 % q_cl )
  call write_diag_3d( "q_rain_np1", lb(1:3), ub(1:3), fields_np1 % q_rain )
  call write_diag_3d( "q_cf_np1", lb(1:3), ub(1:3), fields_np1 % q_cf )
  if ( l_cv_graup )                                                            &
    call write_diag_3d( "q_graup_np1", lb(1:3), ub(1:3), fields_np1 % q_graup )
  if ( l_cv_cloudfrac ) then
    lb(1:3) = lbound(fields_np1 % cf_liq)
    ub(1:3) = ubound(fields_np1 % cf_liq)
    call write_diag_3d( "cf_liq_np1", lb(1:3), ub(1:3), fields_np1 % cf_liq )
    call write_diag_3d( "cf_ice_np1", lb(1:3), ub(1:3), fields_np1 % cf_ice )
    call write_diag_3d( "cf_bulk_np1", lb(1:3), ub(1:3), fields_np1 % cf_bulk )
  end if
  if ( l_tracer ) then
    lb(1:3) = lbound(fields_np1 % tracers(1)%pt)
    ub(1:3) = ubound(fields_np1 % tracers(1)%pt)
    do i_field = 1, n_tracers
      write(num,"(I2)") i_field
      call write_diag_3d( "tracer_" // trim(adjustl(num)) // "_np1",           &
                          lb(1:3), ub(1:3), fields_np1 % tracers(i_field)%pt )
    end do
  end if

  lb(1:3) = lbound( cloudfracs % frac_liq_conv )
  ub(1:3) = ubound( cloudfracs % frac_liq_conv )
  call write_diag_3d( "frac_liq_conv", lb(1:3), ub(1:3),                       &
                      cloudfracs % frac_liq_conv )
  call write_diag_3d( "frac_bulk_conv", lb(1:3), ub(1:3),                      &
                      cloudfracs % frac_bulk_conv )
  call write_diag_3d( "q_cl_conv", lb(1:3), ub(1:3),                           &
                      cloudfracs % q_cl_conv )


  close(unit_diag)

end do  ! n = 1, 3

end subroutine standalone_test


!----------------------------------------------------------------
! Routine to switch on a 2-D diagnostic and assign its pointer to
! the next address in the diags super-array
!----------------------------------------------------------------
subroutine assign_diag_2d( diags_super_2d, diags_list_2d, i_diag, diag )

use diag_type_mod, only: diag_type, diag_list_type
use comorph_constants_mod, only: real_hmprec, nx_full, ny_full

implicit none

! Array to store diagnostics
real(kind=real_hmprec), target, intent(in) :: diags_super_2d                   &
                          ( nx_full, ny_full,                                  &
                            n_diags_max_2d )
! List of diagnostic structures
type(diag_list_type), intent(in out) :: diags_list_2d(n_diags_max_2d)

! Current value of diagnostic counter
integer, intent(in out) :: i_diag

! CoMorph diagnostic structure to be assigned
type(diag_type), target, intent(in out) :: diag

! Set flag to request diag in 2-D in x,y
diag % request % x_y = .true.

! Increment diagnostic counter
i_diag = i_diag + 1

! Point element of diag meta-data pointer list at structure
diags_list_2d(i_diag)%pt => diag

! Assign pointer to the next address in the diags super-array
diag % field_2d => diags_super_2d(:,:,i_diag)

return
end subroutine assign_diag_2d


!----------------------------------------------------------------
! Routine to switch on a 3-D diagnostic and assign its pointer to
! the next address in the diags super-array
!----------------------------------------------------------------
subroutine assign_diag_3d( diags_super_3d, diags_list_3d, i_diag, diag )

use diag_type_mod, only: diag_type, diag_list_type
use comorph_constants_mod, only: real_hmprec, nx_full, ny_full, k_top_conv

implicit none

! Array to store diagnostics
real(kind=real_hmprec), target, intent(in) :: diags_super_3d                   &
                          ( nx_full, ny_full, 1:k_top_conv,                    &
                            n_diags_max_3d )
! List of diagnostic structures
type(diag_list_type), intent(in out) :: diags_list_3d(n_diags_max_3d)

! Current value of diagnostic counter
integer, intent(in out) :: i_diag

! CoMorph diagnostic structure to be assigned
type(diag_type), target, intent(in out) :: diag

! Set flag to request diag in 3-D in x,y,z
diag % request % x_y_z = .true.

! Increment diagnostic counter
i_diag = i_diag + 1

! Point element of diag meta-data pointer list at structure
diags_list_3d(i_diag)%pt => diag

! Assign pointer to the next address in the diags super-array
diag % field_3d => diags_super_3d(:,:,:,i_diag)

return
end subroutine assign_diag_3d


! Same again for a 4-D diagnostic
subroutine assign_diag_4d( n_conv_types, diags_super_4d, diags_list_4d,        &
                           i_diag, diag )

use diag_type_mod, only: diag_type, diag_list_type
use comorph_constants_mod, only: real_hmprec, nx_full, ny_full,                &
                                 n_conv_layers_diag

implicit none

! Number of updraft or downdraft types
integer, intent(in) :: n_conv_types

! Array to store diagnostics
real(kind=real_hmprec), target, intent(in) :: diags_super_4d                   &
                          ( nx_full, ny_full,                                  &
                            n_conv_layers_diag, n_conv_types,                  &
                            n_diags_max_4d )

! List of diagnostic structures
type(diag_list_type), intent(in out) :: diags_list_4d(n_diags_max_4d)

! Current value of diagnostic counter
integer, intent(in out) :: i_diag

! CoMorph diagnostic structure to be assigned
type(diag_type), target, intent(in out) :: diag

! Set flag to request diag in 4-D in x,y,lay,typ
diag % request % x_y_lay_typ = .true.

! Increment diagnostic counter
i_diag = i_diag + 1

! Point element of diag meta-data pointer list at structure
diags_list_4d(i_diag)%pt => diag

! Assign pointer to the next address in the diags super-array
diag % field_4d => diags_super_4d(:,:,:,:,i_diag)

return
end subroutine assign_diag_4d


! Same again for a 5-D diagnostic
subroutine assign_diag_5d( n_conv_types, diags_super_5d, diags_list_5d,        &
                           i_diag, diag )

use diag_type_mod, only: diag_type, diag_list_type
use comorph_constants_mod, only: real_hmprec, nx_full, ny_full, k_top_conv,    &
                                 n_conv_layers_diag

implicit none

! Number of updraft or downdraft types
integer, intent(in) :: n_conv_types

! Array to store diagnostics
real(kind=real_hmprec), target, intent(in) :: diags_super_5d                   &
                          ( nx_full, ny_full, 1:k_top_conv,                    &
                            n_conv_layers_diag, n_conv_types,                  &
                            n_diags_max_5d )

! List of diagnostic structures
type(diag_list_type), intent(in out) :: diags_list_5d(n_diags_max_5d)

! Current value of diagnostic counter
integer, intent(in out) :: i_diag

! CoMorph diagnostic structure to be assigned
type(diag_type), target, intent(in out) :: diag

! Set flag to request diag in 5-D in x,y,z,lay,typ
diag % request % x_y_z_lay_typ = .true.

! Increment diagnostic counter
i_diag = i_diag + 1

! Point element of diag meta-data pointer list at structure
diags_list_5d(i_diag)%pt => diag

! Assign pointer to the next address in the diags super-array
diag % field_5d => diags_super_5d(:,:,:,:,:,i_diag)

return
end subroutine assign_diag_5d


!----------------------------------------------------------------
! Routine to write a single 2-D diagnostic field to an output
! text file, in python numpy array format.
!----------------------------------------------------------------
subroutine write_diag_2d( diag_name, lb, ub, field )

use comorph_constants_mod, only: real_hmprec, nx_full, ny_full

implicit none

character(len=*), intent(in) :: diag_name

integer, intent(in) :: lb(2), ub(2)
real(kind=real_hmprec), intent(in) :: field                                    &
                                      ( lb(1):ub(1), lb(2):ub(2) )

! Write format for diagnostics
character(len=30) :: write_fmt

integer :: i, j

write(write_fmt,"(A4,I3,A16)") "(A1,",nx_full,"(ES18.10,A1),A2)"

write(unit_diag,"(A)") trim(adjustl(diag_name)) // " = np.array(["
do j = 1, ny_full
  write(unit_diag,write_fmt)                                                   &
    "[", ( field(i,j), ",", i=1,nx_full ), "],"
end do
write(unit_diag,"(A)") "])"

return
end subroutine write_diag_2d


!----------------------------------------------------------------
! Routine to write a single 3-D diagnostic field to an output
! text file, in python numpy array format.
!----------------------------------------------------------------
subroutine write_diag_3d( diag_name, lb, ub, field )

use comorph_constants_mod, only: real_hmprec, nx_full, ny_full, k_top_conv

implicit none

character(len=*), intent(in) :: diag_name

integer, intent(in) :: lb(3), ub(3)
real(kind=real_hmprec), intent(in) :: field                                    &
                                      ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )

! Write format for diagnostics
character(len=30) :: write_fmt

integer :: i, j, k

write(write_fmt,"(A4,I3,A16)") "(A1,",nx_full,"(ES18.10,A1),A2)"

write(unit_diag,"(A)") trim(adjustl(diag_name)) // " = np.array(["
do k = 1, k_top_conv
  write(unit_diag,"(A)") "["
  do j = 1, ny_full
    write(unit_diag,write_fmt)                                                 &
      "[", ( field(i,j,k), ",", i=1,nx_full ), "],"
  end do
  write(unit_diag,"(A)") "],"
end do
write(unit_diag,"(A)") "])"

return
end subroutine write_diag_3d


! Same again for a 4-D diagnostic
subroutine write_diag_4d( diag_name, lb, ub, n_conv_types, field )

use comorph_constants_mod, only: real_hmprec, nx_full, ny_full,                &
                                 n_conv_layers_diag

implicit none

character(len=*), intent(in) :: diag_name

integer, intent(in) :: lb(4), ub(4)
integer, intent(in) :: n_conv_types
real(kind=real_hmprec), intent(in) :: field                                    &
           ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3), lb(4):ub(4) )

! Write format for diagnostics
character(len=30) :: write_fmt

integer :: i, j, i_layr, i_type

write(write_fmt,"(A4,I3,A16)") "(A1,",nx_full,"(ES18.10,A1),A2)"

write(unit_diag,"(A)") trim(adjustl(diag_name)) // " = np.array(["
do i_type = 1, n_conv_types
  write(unit_diag,"(A)") "["
  do i_layr = 1, n_conv_layers_diag
    write(unit_diag,"(A)") "["
    do j = 1, ny_full
      write(unit_diag,write_fmt)                                               &
        "[", ( field(i,j,i_layr,i_type), ",", i=1,nx_full ), "],"
    end do
    write(unit_diag,"(A)") "],"
  end do
  write(unit_diag,"(A)") "],"
end do
write(unit_diag,"(A)") "])"

return
end subroutine write_diag_4d


! Same again for a 5-D diagnostic
subroutine write_diag_5d( diag_name, lb, ub, n_conv_types, field )

use comorph_constants_mod, only: real_hmprec, nx_full, ny_full, k_top_conv,    &
                                 n_conv_layers_diag

implicit none

character(len=*), intent(in) :: diag_name

integer, intent(in) :: lb(5), ub(5)
integer, intent(in) :: n_conv_types
real(kind=real_hmprec), intent(in) :: field                                    &
           ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3), lb(4):ub(4), lb(5):ub(5) )

! Write format for diagnostics
character(len=30) :: write_fmt

integer :: i, j, k, i_layr, i_type

write(write_fmt,"(A4,I3,A16)") "(A1,",nx_full,"(ES18.10,A1),A2)"

write(unit_diag,"(A)") trim(adjustl(diag_name)) // " = np.array(["
do i_type = 1, n_conv_types
  write(unit_diag,"(A)") "["
  do i_layr = 1, n_conv_layers_diag
    write(unit_diag,"(A)") "["
    do k = 1, k_top_conv
      write(unit_diag,"(A)") "["
      do j = 1, ny_full
        write(unit_diag,write_fmt)                                             &
          "[", ( field(i,j,k,i_layr,i_type), ",", i=1,nx_full ), "],"
      end do
      write(unit_diag,"(A)") "],"
    end do
    write(unit_diag,"(A)") "],"
  end do
  write(unit_diag,"(A)") "],"
end do
write(unit_diag,"(A)") "])"

return
end subroutine write_diag_5d


end module standalone_test_mod
#endif
