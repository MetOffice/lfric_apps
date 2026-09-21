
#if !defined(LFRIC)
program test_phase_change_solve

use cmpr_type_mod, only: cmpr_type
use comorph_constants_mod, only: real_cvprec, zero, name_length,               &
                                 n_cond_species, n_cond_species_liq,           &
                                 nx_full, ny_full, k_bot_conv, k_top_conv,     &
                                 k_top_init, &
                                 melt_temp, R_dry, R_vap, &
                                 cp_dry, cp_vap, cp_liq, cp_ice, &
                                 L_con_ref, L_fus_ref, rho_liq, rho_ice
use set_dependent_constants_mod, only: set_dependent_constants
use linear_qs_mod, only: n_linear_qs_fields
use moist_proc_diags_type_mod, only: moist_proc_diags_type
use phase_change_solve_mod, only: phase_change_solve

implicit none

integer, parameter :: n_points = 1


! Number of points where each condensed water species is non-zero
integer, allocatable :: nc(:)
! Indices of those points
integer, allocatable :: index_ic(:,:)

type(cmpr_type) :: cmpr
integer :: k
character(len=name_length) :: call_string

! Super-array containing qsat at a reference temperature,
! and dqsat/dT, for linearised qsat calculations
real(kind=real_cvprec) :: linear_qs ( n_points, n_linear_qs_fields )

! Height interval for this step; = the vertical distance
! between the current point and the previous point where
! prev_temp is defined.
! If integrating downwards (as is conventional for
! Eulerian calculations), it should be negative.
real(kind=real_cvprec) :: delta_z(n_points)

! Time interval for converting process rates to increments.
! For Eulerian calculations, this is the model timestep length,
! but for Lagrangian ascents, it is the time taken for the
! parcel to rise over the height interval delta_z,
! so delta_t = delta_z/wind_w
real(kind=real_cvprec) :: delta_t(n_points)

! Vertical wind velocity
! For Euelerian calculations, this is the vertical wind-speed.
! For Lagrangian ascents, it is the vertical velocity of the
! parcel relative to the environment, such that
! wind_w = delta_z/delta_t
real(kind=real_cvprec) :: wind_w(n_points)

! Parcel temperature at the previous model-level
real(kind=real_cvprec) :: prev_temp(n_points)

! Vapour exchange coefficient for each hydrometeor species
real(kind=real_cvprec), allocatable :: kq_cond(:,:)
! Heat exchange coefficient for each hydrometeor species
real(kind=real_cvprec), allocatable :: kt_cond(:,:)
! These have intent inout because they occasionally need to be
! limited for numerical safety reasons.  Also we scale them
! by the timestep delta_t in this routine.

! Fall-speed of each hydrometeor species
real(kind=real_cvprec), allocatable :: wf_cond(:,:)

! Total freezing increment onto each ice hydrometeor species
! (includes homogeneous and heterogeneous freezing and riming)
! Needed for the hydrometeor surface heat budget, important for
! determining the melting rate
real(kind=real_cvprec), allocatable :: dq_frz_cond(:,:)

! Local mixing ratio of each condensed water species,
! implicitly accounting for fall-out from current level / parcel
real(kind=real_cvprec), allocatable :: q_loc_cond(:,:)

! Total available mixing ratio of each condensed water species
! (includes amount that falls through during this step, which
!  maybe considerably larger than the amount actually present
!  at a given instant).
! These are the values updated here; fall-out is calculated
! after this routine.
real(kind=real_cvprec), allocatable :: q_cond(:,:)

! Total heat capacity incremented by phase-changes
real(kind=real_cvprec) :: cp_tot(n_points)

! Parcel air temperature and water vapour mixing ratio
real(kind=real_cvprec) :: temperature(n_points)
real(kind=real_cvprec) :: q_vap(n_points)

! Master switch for diagnostics
logical, parameter :: l_diags = .false.
! Structure storing diagnostics switches and meta-data
type(moist_proc_diags_type) :: moist_proc_diags
integer, parameter :: n_diags = 1
real(kind=real_cvprec) :: diags_super( n_points, n_diags )

integer :: ic, i_cond



! Not used here, but need setting to avoid error trap
nx_full = n_points
ny_full = 1
k_bot_conv = 1
k_top_conv = 10
k_top_init = 9

! Set thermodynamics constants
melt_temp = real( 273.15, real_cvprec )
R_dry = real( 287.05, real_cvprec )
R_vap = real( 287.05/0.62198, real_cvprec )
cp_dry = real( 1005.0, real_cvprec )
cp_vap = zero
cp_liq = zero
cp_ice = zero
L_con_ref = real( 2.501e6, real_cvprec )
L_fus_ref = real( 0.334e6, real_cvprec )
rho_liq = real( 1000.0, real_cvprec )
rho_ice = real( 917.0, real_cvprec )

! Setup constants
call set_dependent_constants()

call_string = "test_phase_change_solve"

allocate( nc ( n_cond_species ) )
allocate( index_ic ( n_points, n_cond_species ) )
allocate( cmpr%index_i(n_points) )
allocate( cmpr%index_j(n_points) )
cmpr%n_points = n_points
allocate( kq_cond ( n_points, n_cond_species ) )
allocate( kt_cond ( n_points, n_cond_species ) )
allocate( wf_cond ( n_points, n_cond_species ) )
allocate( dq_frz_cond ( n_points, n_cond_species_liq+1 : n_cond_species ) )
allocate( q_loc_cond ( n_points, n_cond_species ) )
allocate( q_cond ( n_points, n_cond_species ) )

cmpr%index_i(1) = 1
cmpr%index_j(1) = 1
linear_qs(1,:) = [273.14948,  4.49374365E-3,  4.49329196E-3, &
                              3.28750612E-4,  3.726164E-4]
delta_z(1) = 86.81665
delta_t(1) = 347.2666
wind_w(1) = 0.25
prev_temp(1) = 273.67719
kq_cond(1,:) = [0.17745258,  7.7014578E-5,  1.61309999E-5,  5.90306045E-5]
kt_cond(1,:) = [0.1204145,  5.54143517E-5,  1.09460952E-5,  5.06376382E-5]
wf_cond(1,:) = [2.1684817E-3,  1.8552991,  1.82522414E-3,  4.1219993]
dq_frz_cond(1,:) = [0.,  1.01640635E-5]
q_loc_cond(1,:) = [4.22776575E-5, 2.29610287E-5, 3.23526983E-9, 8.08190598E-5]
q_cond(1,:) = [3.59028309E-5,  2.79465403E-5,  3.23648353E-9,  1.59455813E-4]
cp_tot(1) = 1005.0
temperature(1) = 272.83563
q_vap(1) = 4.62010596E-3

do i_cond = 1, n_cond_species
  nc(i_cond) = 0.0
  do ic = 1, n_points
    if ( kq_cond(ic,i_cond) > zero ) then
      nc(i_cond) = nc(i_cond) + 1
      index_ic(nc(i_cond),i_cond) = ic
    end if
  end do
end do

call phase_change_solve( n_points, n_points,                                   &
             nc, index_ic, cmpr, k, call_string, linear_qs,                    &
             delta_z, delta_t, wind_w, prev_temp,                              &
             kq_cond, kt_cond, wf_cond, dq_frz_cond,                           &
             q_loc_cond, q_cond, cp_tot, temperature, q_vap,                   &
             l_diags, moist_proc_diags, n_points, n_diags, diags_super )

end program test_phase_change_solve
#endif
