! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

#if !defined(LFRIC)
! Unit test for the CoMorph moist processes subroutine.

! This program does an idealised moist parcel ascent and descent,
! calling the subroutine moist_proc to do all the phase-change
! and microphysical calculations.
!
! Bit-reproducibility on different decompositions is tested for
! by doing the calculations twice; once calling the subroutine
! on all columns at once in vectorised form, and once calling
! it at one grid-column at a time.  Any difference between
! the two out files will indicate a bug in the vectorisation
! code within moist_proc.

module test_moist_proc_sizes
implicit none

! Number of points
integer, parameter :: n_points = 5
! Number of levels
integer, parameter :: n_levels = 50

end module test_moist_proc_sizes



program test_moist_proc

use test_moist_proc_sizes, only: n_points, n_levels
use comorph_constants_mod, only: real_cvprec, n_cond_species, cond_params,     &
                                 zero, one, half, gravity, cp_dry, R_dry, pi,  &
                                 name_length, i_cond_cl, i_cond_cf, melt_temp, &
                                 nx_full, ny_full, k_top_conv
use set_dependent_constants_mod, only: set_dependent_constants
use linear_qs_mod, only: n_linear_qs_fields, i_ref_temp,                       &
                         i_qsat_liq_ref, i_qsat_ice_ref,                       &
                         i_dqsatdT_liq, i_dqsatdt_ice,                         &
                         linear_qs_set_ref
use moist_proc_diags_type_mod, only: moist_proc_diags_type,                    &
                                     moist_proc_diags_assign
use diag_type_mod, only: diag_list_type, dom_type
use cmpr_type_mod, only: cmpr_type, cmpr_alloc

use set_qsat_mod, only: set_qsat_liq, set_qsat_ice
use dry_adiabat_mod, only: dry_adiabat
use calc_rho_dry_mod, only: calc_rho_dry
use set_cp_tot_mod, only: set_cp_tot
use sat_adjust_mod, only: sat_adjust
use moist_proc_mod, only: moist_proc

implicit none

! Model grid
real(kind=real_cvprec) :: height_full(n_points,0:n_levels)
real(kind=real_cvprec) :: height_half(n_points,n_levels)
real(kind=real_cvprec) :: pressure_full(n_points,0:n_levels)
real(kind=real_cvprec) :: pressure_half(n_points,n_levels)

! Environment profile
real(kind=real_cvprec) :: env_wind_u(n_points,n_levels)
real(kind=real_cvprec) :: env_wind_v(n_points,n_levels)
real(kind=real_cvprec) :: env_wind_w(n_points,n_levels)
real(kind=real_cvprec) :: env_temperature(n_points,0:n_levels)

! Output parcel profiles
real(kind=real_cvprec) :: out_wind_u(n_points,n_levels)
real(kind=real_cvprec) :: out_wind_v(n_points,n_levels)
real(kind=real_cvprec) :: out_wind_w(n_points,n_levels)
real(kind=real_cvprec) :: out_temperature(n_points,n_levels)
real(kind=real_cvprec) :: out_q_vap(n_points,n_levels)
real(kind=real_cvprec), allocatable :: out_q_cond(:,:,:)
real(kind=real_cvprec), allocatable :: out_flux_cond(:,:,:)

! Super-saturation and relative humidity diagnostics
real(kind=real_cvprec) :: ss_liq(n_points,n_levels)
real(kind=real_cvprec) :: ss_ice(n_points,n_levels)
real(kind=real_cvprec) :: rh_liq(n_points,n_levels)
real(kind=real_cvprec) :: rh_ice(n_points,n_levels)

! Local parcel properties
real(kind=real_cvprec) :: par_wind_u(n_points)
real(kind=real_cvprec) :: par_wind_v(n_points)
real(kind=real_cvprec) :: par_wind_w(n_points)
real(kind=real_cvprec) :: par_temperature(n_points)
real(kind=real_cvprec) :: par_q_vap(n_points)
real(kind=real_cvprec), allocatable :: par_q_cond(:,:)
real(kind=real_cvprec), allocatable :: par_q_cond_1(:)

! Re-occurring factor dt / (rho_dry * lz)
real(kind=real_cvprec) :: dt_over_rhod_lz(n_points)

! Other inputs to moist_proc
real(kind=real_cvprec) :: linear_qs(n_points,n_linear_qs_fields)
real(kind=real_cvprec) :: linear_qs_1(n_linear_qs_fields)
real(kind=real_cvprec) :: delta_z(n_points)
real(kind=real_cvprec) :: delta_t(n_points)
real(kind=real_cvprec) :: vert_len(n_points)
real(kind=real_cvprec) :: wind_w_excess(n_points)
real(kind=real_cvprec), allocatable :: flux_cond(:,:)
real(kind=real_cvprec), allocatable :: flux_cond_1(:)

! Master switch for diagnostics; hardwired on here.
logical, parameter :: l_diags = .true.
! Structure storing diagnostics switches and meta-data
type(moist_proc_diags_type) :: moist_proc_diags

! Super-array for diagnostics passed out of moist_proc
real(kind=real_cvprec), allocatable :: par_diags_super(:,:)
real(kind=real_cvprec), allocatable :: par_diags_super_1(:)
! Super-array for 2-D profile diagnostics
real(kind=real_cvprec), target, allocatable ::                                 &
                                       out_diags_super(:,:,:)
! Diags list object for setting up the diags infrastructure
type(diag_list_type), allocatable :: diag_list(:)
type(dom_type) :: doms
logical :: l_count_diags
integer :: i_diag, i_super, n_diags
character(len=name_length) :: parent_name

! Work variables
real(kind=real_cvprec) :: lapse_rate
real(kind=real_cvprec) :: interp
real(kind=real_cvprec) :: qs(n_points)

! Top of model height
real(kind=real_cvprec), parameter :: z_top = 20000.0_real_cvprec
! Height-scale used to set an idealised wind profile:
real(kind=real_cvprec), parameter :: z_scl_u = 500.0_real_cvprec
! Height-scale used to set an idealised temperature profile:
real(kind=real_cvprec), parameter :: z_scl_t = 5000.0_real_cvprec
! Parcel initial Relative Humidity
real(kind=real_cvprec), parameter :: rh_init = 0.8_real_cvprec

! Parcel ascent / descent rate
real(kind=real_cvprec) :: ascent_rate(n_points)

! Flag for turning on ice seeding
logical, parameter :: l_ice_seed = .true.

! Flag passed into sat_adjust to tell it not to update
! q_vap and q_cl when calculating a saturated reference T
logical, parameter :: l_update_q_false = .false.

! Parcel total heat-capacity for calculation of saturated reference temperature
real(kind=real_cvprec) :: cp_tot(n_points)

! Filename for output
character(len=name_length) :: filename

! Stuff passed into moist_proc for error reporting
type(cmpr_type) :: cmpr
character(len=name_length) :: call_string

integer :: k_first, k_last, dk, k_half_prev, k_half_next

! Loop counters
integer :: ic, k, i_cond, n, m


! These need to be set before calling set_dependent_constants,
! even though not used in this test:
nx_full = 1
ny_full = 1
k_top_conv = 1

! Setup constants
call set_dependent_constants()

! Allocate arrays whose dimensions depend on n_cond_species,
! which is set in the above call
allocate( out_q_cond(n_points,n_levels,n_cond_species) )
allocate( out_flux_cond(n_points,n_levels,n_cond_species) )
allocate( par_q_cond(n_points,n_cond_species) )
allocate( par_q_cond_1(n_cond_species) )
allocate( flux_cond(n_points,n_cond_species) )
allocate( flux_cond_1(n_cond_species) )

! Setup a basic state...

! Set model-level heights:
! heights increase as a quadratic formula of model-level index.
do k = 0, n_levels
  do ic = 1, n_points
    ! Full levels use function of k
    height_full(ic,k) = z_top * ( real(k,real_cvprec)                          &
                                 /real(n_levels,real_cvprec) )**2
  end do
end do
do k = 1, n_levels
  do ic = 1, n_points
    ! Half-levels use same function but with k-1/2
    height_half(ic,k) = z_top * ( (real(k,real_cvprec)-half)                   &
                                 /real(n_levels,real_cvprec) )**2
  end do
end do

! Set made-up wind profiles
do k = 1, n_levels
  do ic = 1, n_points
    ! Log profile of u, zero for v,w
    env_wind_u(ic,k) = 20.0_real_cvprec                                        &
                     * log((height_full(ic,k)+z_scl_u)/z_scl_u)
    env_wind_v(ic,k) = zero
    env_wind_w(ic,k) = zero
  end do
end do

! Set temperature profile...
do ic = 1, n_points
  env_temperature(ic,0) = 300.0_real_cvprec
end do
do k = 1, n_levels
  do ic = 1, n_points
    ! Set temperature lapse-rate; stable in lower troposphere,
    ! gradually tending towards dry-adiabatic in the upper troposphere
    lapse_rate = (gravity/cp_dry)                                              &
       * ( one - half / ( one + height_full(ic,k)/z_scl_t ) )
    ! Calc T at next level based on the lapse-rate
    env_temperature(ic,k) = env_temperature(ic,k-1)                            &
       - lapse_rate * ( height_full(ic,k) - height_full(ic,k-1) )
  end do
end do

! Set pressure approximately in hydrostatic balance:
! dp/dz = -rho g = -p/(R Tv) g
do ic = 1, n_points
  pressure_full(ic,0) = 100000.0_real_cvprec
end do
do k = 1, n_levels
  do ic = 1, n_points
    pressure_full(ic,k) = pressure_full(ic,k-1)                                &
      - ( height_full(ic,k) - height_full(ic,k-1) )                            &
        * gravity * pressure_full(ic,k-1)                                      &
        / ( R_dry * env_temperature(ic,k-1) )
  end do
end do
! Interpolate pressure onto half-levels
do k = 1, n_levels
  do ic = 1, n_points
    interp = ( height_half(ic,k) - height_full(ic,k-1) )                       &
           / ( height_full(ic,k) - height_full(ic,k-1) )
    pressure_half(ic,k) = (one-interp) * pressure_full(ic,k-1)                 &
                        +      interp  * pressure_full(ic,k)
  end do
end do

! Set parcel ascent / descent rate to vary over the points
do ic = 1, n_points
  ! Quadratic variation:
  ascent_rate(ic) = 10.0_real_cvprec                                           &
   * ( real(ic,real_cvprec) / real(n_points,real_cvprec) )**2
end do



! Setup moist_proc diagnostics...

! Turn diagnostics on for liquid cloud
moist_proc_diags % diags_cl % dq_fall % request%x_y_z = .true.
moist_proc_diags % diags_cl % dq_cond % request%x_y_z = .true.
moist_proc_diags % diags_cl % dq_frzmlt % request%x_y_z = .true.
moist_proc_diags % diags_cl % dq_col % request%x_y_z = .true.
moist_proc_diags % diags_cl % fall_speed % request%x_y_z = .true.
moist_proc_diags % diags_cl % cond_temp % request%x_y_z = .true.
! Set other species to have the same diags turned on
moist_proc_diags % diags_rain  = moist_proc_diags % diags_cl
moist_proc_diags % diags_cf    = moist_proc_diags % diags_cl
moist_proc_diags % diags_graup = moist_proc_diags % diags_cl
! For rain, also turn on autoconversion increment
moist_proc_diags % diags_rain % dq_aut % request%x_y_z = .true.

! Count number of active diags
parent_name = "par"
doms%x_y_z = .true.
allocate( diag_list(1) )
i_diag = 0
i_super = 0
l_count_diags = .true.
call moist_proc_diags_assign( parent_name, l_count_diags, doms,                &
                              moist_proc_diags, diag_list,                     &
                              i_diag, i_super )

! Allocate diags list to full size and assign contained pointers
deallocate( diag_list )
allocate( diag_list( moist_proc_diags % n_diags ) )
i_diag = 0
i_super = 0
l_count_diags = .false.
call moist_proc_diags_assign( parent_name, l_count_diags, doms,                &
                              moist_proc_diags, diag_list,                     &
                              i_diag, i_super )


! Allocate diag arrays
n_diags = moist_proc_diags % n_diags
allocate( par_diags_super ( n_points, n_diags ) )
allocate( par_diags_super_1 ( n_diags ) )
allocate( out_diags_super ( n_points, n_levels, n_diags ) )

call cmpr_alloc( cmpr, n_points )
cmpr % index_j(1:n_points) = 1



! Do everything twice to check bit-reproducibility on
! different decompositions
do n = 1, 2

  ! Setup stuff for error reporting inside moist_proc
  write(call_string,"(A,I2)") "test_moist_proc call n = ", n
  if ( n == 1 ) then
    cmpr % n_points = n_points
    do ic = 1, n_points
      cmpr % index_i(ic) = ic
    end do
  else if ( n == 2 ) then
    cmpr % n_points = 1
  end if

  ! Initialise parcel temperature and winds using 1st model-level
  do ic = 1, n_points
    par_wind_u(ic) = env_wind_u(ic,1)
    par_wind_v(ic) = env_wind_v(ic,1)
    par_wind_w(ic) = env_wind_w(ic,1)
    par_temperature(ic) = env_temperature(ic,1)
  end do
  ! Set initial condensate to zero
  do i_cond = 1, n_cond_species
    do ic = 1, n_points
      par_q_cond(ic,i_cond) = zero
    end do
  end do
  ! Set initial vapour to fixed RH
  call set_qsat_liq( n_points,                                                 &
                     par_temperature, pressure_full(:,1), par_q_vap )
  do ic = 1, n_points
    par_q_vap(ic) = rh_init * par_q_vap(ic)
  end do


  ! Adjust temperature down to first half-level
  call dry_adiabat( n_points, n_points,                                        &
                    pressure_full(:,1), pressure_half(:,1),                    &
                    par_q_vap, par_q_cond, par_temperature )

  ! Set parcel vertical length-scale
  do ic = 1, n_points
    vert_len(ic) = 1000.0_real_cvprec
  end do

  ! Set lowest model-level in output parcel property profiles
  do ic = 1, n_points
    out_wind_u(ic,1) = par_wind_u(ic)
    out_wind_v(ic,1) = par_wind_v(ic)
    out_wind_w(ic,1) = par_wind_w(ic)
    out_temperature(ic,1) = par_temperature(ic)
    out_q_vap(ic,1) = par_q_vap(ic)
  end do
  do i_cond = 1, n_cond_species
    do ic = 1, n_points
      out_q_cond(ic,1,i_cond) = par_q_cond(ic,i_cond)
    end do
  end do
  ! Set diagnostics to zero in lowest model-level
  do ic = 1, n_points
    ss_liq(ic,1) = zero
    rh_liq(ic,1) = zero
    ss_ice(ic,1) = zero
    rh_ice(ic,1) = zero
  end do
  do i_cond = 1, n_cond_species
    do ic = 1, n_points
      out_flux_cond(ic,n_levels,i_cond) = zero
    end do
  end do
  do i_diag = 1, moist_proc_diags % n_diags
    do ic = 1, n_points
      out_diags_super(ic,1,i_diag) = zero
    end do
  end do


  ! Do 2 sweeps, one up and one down
  do m = 1, 2


    if ( m == 1 ) then
      k_first = 1
      k_last = n_levels - 1
      dk = 1
      ! Set parcel ascent rate
      do ic = 1, n_points
        wind_w_excess(ic) = ascent_rate(ic)
      end do
    else
      k_first = n_levels - 1
      k_last = 1
      dk = -1
      ! Set parcel descent rate
      do ic = 1, n_points
        wind_w_excess(ic) = -ascent_rate(ic)
      end do
    end if


    ! Loop over levels
    do k = k_first, k_last, dk

      if ( m == 1 ) then
        k_half_prev = k
        k_half_next = k + 1
      else
        k_half_prev = k + 1
        k_half_next = k
      end if


      ! Set height intervals etc
      do ic = 1, n_points
        delta_z(ic) = height_half(ic,k_half_next)                              &
                    - height_half(ic,k_half_prev)
        delta_t(ic) = delta_z(ic) / wind_w_excess(ic)
      end do

      if ( m == 1 ) then
        ! Set fall-in flux of each species to zero on ascent
        do i_cond = 1, n_cond_species
          do ic = 1, n_points
            flux_cond(ic,i_cond) = zero
          end do
        end do
        if ( l_ice_seed ) then
          ! If using ice seeding, set a small non-zero flux of
          ! ice when below freezing.
          do ic = 1, n_points
            if ( par_temperature(ic) < melt_temp ) then
              flux_cond(ic,i_cond_cf) = 1.0e-6_real_cvprec
            end if
          end do
        end if
      else
        ! Set fall-in flux using the ascent fall-out fluxes
        ! in the descent
        do i_cond = 1, n_cond_species
          do ic = 1, n_points
            flux_cond(ic,i_cond) = out_flux_cond(ic,k,i_cond)
          end do
        end do
      end if

      if ( n == 1 ) then
        ! Vectorised form

        ! Lift parcel to next half-level
        call dry_adiabat( n_points, n_points,                                  &
                          pressure_half(:,k_half_prev),                        &
                          pressure_half(:,k_half_next),                        &
                          par_q_vap, par_q_cond, par_temperature )

        ! Calculate factor delta_t / ( rho_dry vert_len )
        call calc_rho_dry( n_points,                                           &
                           par_temperature, par_q_vap,                         &
                           pressure_half(:,k_half_next),                       &
                           dt_over_rhod_lz )
        do ic = 1, n_points
          dt_over_rhod_lz(ic) = delta_t(ic)                                    &
                        / ( dt_over_rhod_lz(ic) * vert_len(ic) )
        end do


        ! Set reference temperature and qsat
        do ic = 1, n_points
          linear_qs(ic,i_ref_temp) = par_temperature(ic)
        end do
        ! Adjust the reference temperature to saturation
        call set_cp_tot( n_points, n_points,                                   &
                         par_q_vap, par_q_cond, cp_tot )
        call sat_adjust( n_points, l_update_q_false,                           &
                         pressure_half(:,k_half_next),                         &
                         linear_qs(:,i_ref_temp),                              &
                         par_q_vap, par_q_cond(:,i_cond_cl), cp_tot )
        ! Call routine to compute qsat and dqsat/dT at reference temperature
        call linear_qs_set_ref( n_points,                                      &
                                pressure_half(:,k_half_next),                  &
                                linear_qs )

        ! Initialise diags to zero
        do i_diag = 1, n_diags
          do ic = 1, n_points
            par_diags_super(ic,i_diag) = zero
          end do
        end do

        ! Call moist_proc
        call moist_proc( n_points, n_points, linear_qs,                        &
               delta_z, delta_t, vert_len, wind_w_excess,                      &
               dt_over_rhod_lz, pressure_half(:,k_half_next),                  &
               out_temperature(:,k_half_prev),                                 &
               env_wind_u(:,k), env_wind_v(:,k),                               &
               env_wind_w(:,k), env_temperature(:,k),                          &
               par_wind_u, par_wind_v, par_wind_w,                             &
               par_temperature, par_q_vap, par_q_cond,                         &
               flux_cond, cmpr, k, call_string, l_diags,                       &
               moist_proc_diags, n_points, n_diags, par_diags_super )

      else  ! ( n == 1 )
        ! Single-column form

        do ic = 1, n_points

          cmpr % index_i(1) = ic

          ! Fields contained in super-arrays; need to make a copy
          ! containing only the fields for the current point,
          ! as subsetting the arrays in the argument list to moist_proc
          ! to only pass in the current point triggers a warning about
          ! creating an array temporary, since the array sections are
          ! not contiguous in memory.
          do i_cond = 1, n_cond_species
            flux_cond_1(i_cond) = flux_cond(ic,i_cond)
            par_q_cond_1(i_cond) = par_q_cond(ic,i_cond)
          end do

          ! Lift parcel to next half-level
          call dry_adiabat( 1, 1,                                              &
                            pressure_half(ic,k_half_prev),                     &
                            pressure_half(ic,k_half_next),                     &
                            par_q_vap(ic), par_q_cond_1,                       &
                            par_temperature(ic) )

          ! Set dry-density
          call calc_rho_dry( 1,                                                &
                             par_temperature(ic), par_q_vap(ic),               &
                             pressure_half(ic,k_half_next),                    &
                             dt_over_rhod_lz(ic) )
          dt_over_rhod_lz(ic) = delta_t(ic)                                    &
                        / ( dt_over_rhod_lz(ic) * vert_len(ic) )

          ! Set reference temperature and qsat
          linear_qs_1(i_ref_temp) = par_temperature(ic)
          call set_cp_tot( 1, 1,                                               &
                           par_q_vap(ic), par_q_cond_1, cp_tot(ic) )
          call sat_adjust( 1, l_update_q_false,                                &
                           pressure_half(ic,k_half_next),                      &
                           linear_qs_1(i_ref_temp),                            &
                           par_q_vap(ic), par_q_cond_1(i_cond_cl), cp_tot(ic) )
          call linear_qs_set_ref( 1,                                           &
                                  pressure_half(ic,k_half_next),               &
                                  linear_qs_1 )

          ! Initialise diags to zero
          do i_diag = 1, n_diags
            par_diags_super_1(i_diag) = zero
          end do

          ! Call moist_proc
          call moist_proc( 1, 1, linear_qs_1,                                  &
                 delta_z(ic), delta_t(ic), vert_len(ic),                       &
                 wind_w_excess(ic),                                            &
                 dt_over_rhod_lz(ic),                                          &
                 pressure_half(ic,k_half_next),                                &
                 out_temperature(ic,k_half_prev),                              &
                 env_wind_u(ic,k), env_wind_v(ic,k),                           &
                 env_wind_w(ic,k), env_temperature(ic,k),                      &
                 par_wind_u(ic), par_wind_v(ic),                               &
                 par_wind_w(ic),                                               &
                 par_temperature(ic), par_q_vap(ic),                           &
                 par_q_cond_1,                                                 &
                 flux_cond_1, cmpr, k, call_string, l_diags,                   &
                 moist_proc_diags, 1, n_diags, par_diags_super_1 )

          ! Copy the updated single-point super-array fields
          ! back into their original super-arrays
          do i_cond = 1, n_cond_species
            flux_cond(ic,i_cond) = flux_cond_1(i_cond)
            par_q_cond(ic,i_cond) = par_q_cond_1(i_cond)
          end do
          do i_cond = 1, n_linear_qs_fields
            linear_qs(ic,i_cond) = linear_qs_1(i_cond)
          end do
          do i_diag = 1, moist_proc_diags % n_diags
            par_diags_super(ic,i_diag) =par_diags_super_1(i_diag)
          end do

        end do  ! ic = 1, n_points

      end if  ! ( n == 1 )

      ! Copy parcel properties into output arrays
      do ic = 1, n_points
        out_wind_u(ic,k_half_next) = par_wind_u(ic)
        out_wind_v(ic,k_half_next) = par_wind_v(ic)
        out_wind_w(ic,k_half_next) = par_wind_w(ic)
        out_temperature(ic,k_half_next) = par_temperature(ic)
        out_q_vap(ic,k_half_next) = par_q_vap(ic)
      end do
      do i_cond = 1, n_cond_species
        do ic = 1, n_points
          out_q_cond(ic,k_half_next,i_cond)                                    &
            = par_q_cond(ic,i_cond)
          out_flux_cond(ic,k,i_cond) = flux_cond(ic,i_cond)
        end do
      end do
      ! Copy diagnostics
      do i_diag = 1, moist_proc_diags % n_diags
        do ic = 1, n_points
          out_diags_super(ic,k_half_next,i_diag)                               &
            = par_diags_super(ic,i_diag)
        end do
      end do

      ! Supersaturation and relative humidity diagnostics
      call set_qsat_liq( n_points,                                             &
                         par_temperature, pressure_half(:,k_half_next), qs )
      do ic = 1, n_points
        ss_liq(ic,k_half_next) = par_q_vap(ic) - qs(ic)
        rh_liq(ic,k_half_next) = 100.0_real_cvprec                             &
                                 * par_q_vap(ic) / qs(ic)
      end do
      call set_qsat_ice( n_points,                                             &
                         par_temperature, pressure_half(:,k_half_next), qs )
      do ic = 1, n_points
        ss_ice(ic,k_half_next) = par_q_vap(ic) - qs(ic)
        rh_ice(ic,k_half_next) = 100.0_real_cvprec                             &
                                 * par_q_vap(ic) / qs(ic)
      end do

    end do


    ! Write out fields into a text file in python numpy module
    ! format, so that the data can be loaded directly into python...

    if ( m == 1 ) then
      write(filename,"(A,I1,A)") "test_moist_proc_up_out", n, ".py"
    else
      write(filename,"(A,I1,A)") "test_moist_proc_dn_out", n, ".py"
    end if

    open(unit=10,file=filename)

    write(10,"(A)") "import numpy as np"
    write(10,"(A)") ""


    call write_diag( "height_full", height_full(:,1:) )
    call write_diag( "pressure_full", pressure_full(:,1:) )
    call write_diag( "height_half", height_half )
    call write_diag( "pressure_half", pressure_half )

    call write_diag( "env_wind_u", env_wind_u )
    call write_diag( "env_wind_v", env_wind_v )
    call write_diag( "env_wind_w", env_wind_w )
    call write_diag( "env_temperature", env_temperature )

    call write_diag( "par_wind_u", out_wind_u )
    call write_diag( "par_wind_v", out_wind_v )
    call write_diag( "par_wind_w", out_wind_w )
    call write_diag( "par_temperature", out_temperature )
    call write_diag( "par_q_vap", out_q_vap )
    do i_cond = 1, n_cond_species
      call write_diag( "par_q_" //                                             &
                        trim(adjustl(                                          &
                        cond_params(i_cond)%pt % cond_name )),                 &
                        out_q_cond(:,:,i_cond) )
    end do

    do i_cond = 1, n_cond_species
      call write_diag( "flux_" //                                              &
                        trim(adjustl(                                          &
                        cond_params(i_cond)%pt % cond_name )),                 &
                        out_flux_cond(:,:,i_cond) )
    end do

    call write_diag( "ss_liq", ss_liq )
    call write_diag( "rh_liq", rh_liq )
    call write_diag( "ss_ice", ss_ice )
    call write_diag( "rh_ice", rh_ice )

    do i_diag = 1, moist_proc_diags % n_diags
      i_super = moist_proc_diags % list(i_diag)%pt % i_super
      call write_diag(                                                         &
             moist_proc_diags % list(i_diag)%pt % diag_name,                   &
             out_diags_super(:,:,i_super) )
    end do

    close(10)


  end do  ! m = 1, 2


end do  ! n = 1, 2


end program test_moist_proc



! Subroutine to write out diagnostics
subroutine write_diag( diag_name, field )

use comorph_constants_mod, only: real_cvprec
use test_moist_proc_sizes, only: n_points, n_levels

implicit none

character(len=*), intent(in) :: diag_name

real(kind=real_cvprec), intent(in) :: field(n_points,n_levels)

! Write format for diagnostics
character(len=30) :: write_fmt

integer :: ic, k

write(write_fmt,"(A4,I3,A16)") "(A1,",n_points,"(ES18.10,A1),A2)"

write(10,"(A)") trim(adjustl(diag_name)) // " = np.array(["
do k = 1, n_levels
  write(10,write_fmt)                                                          &
      "[", ( field(ic,k), ",", ic=1,n_points ), "],"
end do
write(10,"(A)") "])"

return
end subroutine write_diag
#endif
