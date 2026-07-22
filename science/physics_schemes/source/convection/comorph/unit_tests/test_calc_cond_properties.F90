! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

#if !defined(LFRIC)
! Unit test for the CoMorph microphysics routine calc_cond_properties,
! which is part of BOTEMS (Back Of The Envelope Microphysics Scheme).
! calc_cond_properties performs an implicit solution for the fall-speed,
! number concentration and mixing-ratio of a given hydrometeor
! species, accounting for the fall-out of the hydrometeor during the
! current timestep.  This unit test checks that the implicit solution
! works, and bit-reproduces when vectorised vs doing one grid-point at
! a time.
program test_calc_cond_properties

use comorph_constants_mod, only: real_cvprec, params_rain, q_activate, zero,   &
                                 nx_full, ny_full, k_top_conv
use set_dependent_constants_mod, only: set_dependent_constants

use calc_cond_properties_mod, only: calc_cond_properties,                      &
                                    calc_cond_properties_cmpr
use fall_speed_mod, only: fall_speed

implicit none

! Number of points
integer, parameter :: n_points = 200

! Inputs to calc_cond_properties:
! Reference temperature used for T-dependent number concentration
real(kind=real_cvprec) :: ref_temp(n_points)
! Initial mixing-ratio of rain before fall-out
real(kind=real_cvprec) :: q_rain(n_points)
! Dry and wet density
real(kind=real_cvprec) :: rho_dry(n_points)
real(kind=real_cvprec) :: rho_wet(n_points)
! Factor dt / lz, where dt = timestep,
!and lz = vertical thickness of parcel / layer.
real(kind=real_cvprec) :: dt_over_lz(n_points)

! Outputs from calc_cond_properties:
! Implicitly-solved rain mixing-ratio after fall-out.
real(kind=real_cvprec) :: q_loc_rain(n_points)
! Number concentration of rain.
real(kind=real_cvprec) :: n_rain(n_points)
! Radius of rain droplets consistent with the solved mixing-ratio and number.
real(kind=real_cvprec) :: r_rain(n_points)
! Implicitly-solved fall-speed of the rain.
real(kind=real_cvprec) :: wf_rain(n_points)
! Coefficients for exchange of vapour and heat between the
! rain-drops and the surrounding air, accounting for ventilation effects.
real(kind=real_cvprec) :: kq_rain(n_points)
real(kind=real_cvprec) :: kt_rain(n_points)

! Fall-speed computed from the output solved rain droplet size,
! to check consistency of the implicit solution.
real(kind=real_cvprec) :: wf_after(n_points)

! Work variable used to set initial q_rain
real(kind=real_cvprec) :: frac

! Compression indices
integer :: nc
integer :: index_ic(n_points)

! Loop counter
integer :: ic


! These need to be set before calling set_dependent_constants,
! even though not used in this test:
nx_full = 1
ny_full = 1
k_top_conv = 1

! Setup constants
call set_dependent_constants()


! Invent some inputs:

! Set reference temperature to 0oC
ref_temp(:) = 273.0

! Set the initial q_rain mixing-ratio...
do ic = 1, n_points
  ! Set fraction of the way through the points in the array
  frac = real(ic,real_cvprec)/real(n_points,real_cvprec)
  if (frac <= 0.6) then
    ! First 3/5ths of the points; exponential function to test
    ! very large range of orders of magnitude for initial q_rain.
    q_rain(ic) = exp( -40.0 * frac )
  else if (frac <= 0.8) then
    ! Next 1/5th of the points; set to fixed extremely small value.
    q_rain(ic) = q_activate
  else
    ! last 1/5th; set to zero to test gather/scatter onto points with
    ! non-zero qrain.
    q_rain(ic) = zero
  end if
end do

! Set densities to 1.0
rho_dry(:) = 1.0_real_cvprec
rho_wet(:) = 1.0_real_cvprec

! Set dt/lz to a very large number; this pushes the implicit solve
! to its limits, since any rain present will want to very-nearly
! entirely fall out during the timestep.
dt_over_lz(:) = 10000.0_real_cvprec

! Set compression list of points where q_rain is nonzero
nc = 0
do ic = 1, n_points
  if ( q_rain(ic) > zero ) then
    nc = nc + 1
    index_ic(nc) = ic
  end if
end do

! Compute rain properties; includes implicit solution for the
! fall-speed
call calc_cond_properties_cmpr( n_points, nc, index_ic,                        &
                                params_rain, ref_temp, q_rain,                 &
                                rho_dry, rho_wet,                              &
                                dt_over_lz,                                    &
                                q_loc_rain, n_rain,                            &
                                r_rain, wf_rain,                               &
                                kq_rain, kt_rain )

! Calculate fall-speed using the solved consistent
! hydrometeor radius
call fall_speed( n_points,                                                     &
                 params_rain % area_coef, params_rain % rho,                   &
                 rho_wet, r_rain, wf_after )


! Write outputs to a file, including fall-speed before and after
! for comparison (the implicit solution is meant to find
! wf_rain such that wf_after = wf_rain).
open( 10, file="test_calc_cond_properties_out1.txt" )
do ic = 1, n_points
  write(10,"(7ES18.10)") q_rain(ic), q_loc_rain(ic), r_rain(ic),               &
                         kq_rain(ic), kt_rain(ic), wf_rain(ic), wf_after(ic)
end do
close(10)


! Repeat the above calculations and output, but calling the
! routines one grid-point at a time.  The output from this should
! be identical to that from the above.  If it differs, there
! must be a bug in the vectorisation / compression.

do ic = 1, n_points
  if ( q_rain(ic) > zero ) then
    call calc_cond_properties( 1,                                              &
                               params_rain, ref_temp, q_rain(ic),              &
                               rho_dry(ic), rho_wet(ic),                       &
                               dt_over_lz(ic),                                 &
                               q_loc_rain(ic), n_rain(ic),                     &
                               r_rain(ic), wf_rain(ic),                        &
                               kq_rain(ic), kt_rain(ic) )
    call fall_speed( 1,                                                        &
                     params_rain % area_coef, params_rain % rho,               &
                     rho_wet(ic), r_rain(ic), wf_after(ic) )
  end if
end do

open( 10, file="test_calc_cond_properties_out2.txt" )
do ic = 1, n_points
  write(10,"(7ES18.10)") q_rain(ic), q_loc_rain(ic), r_rain(ic),               &
                         kq_rain(ic), kt_rain(ic), wf_rain(ic), wf_after(ic)
end do
close(10)


end program test_calc_cond_properties
#endif
