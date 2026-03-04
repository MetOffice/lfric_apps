! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module set_det_mod

implicit none

contains

! Routine to compute the amount of mass-flux detrained during the current
! level-step, and find the mean properties of the detrained and
! non-detrained parcels.  Detrainement is solved implicitly w.r.t.
! the environment virtual temperature, which is itself modified by
! the compensating subsidence from the non-detrained parcel (and so
! depends on the detrainment rate...)
subroutine set_det( n_points, max_points, n_points_res, n_fields_tot,          &
                    l_down, l_last_level, l_to_full_level,                     &
                    cmpr, k, draft_string,                                     &
                    sublevs, i_next, i_next_max,                               &
                    core_mean_ratio,                                           &
                    par_mean_w_drag, par_core_w_drag, res_source_fields,       &
                    par_next_mean_fields, par_next_core_fields,                &
                    par_prev_super, par_conv_super,                            &
                    det_mass_d, det_fields )

use comorph_constants_mod, only: real_cvprec, zero, one, two, half,            &
                                 i_check_bad_values_cmpr, i_check_bad_none,    &
                                 l_par_core, name_length,                      &
                                 min_float, min_delta, sqrt_min_delta, newline

use cmpr_type_mod, only: cmpr_type
use fields_type_mod, only: field_positive, field_names, i_wind_w,              &
                           i_qc_first, i_qc_last
use sublevs_mod, only: max_sublevs, n_sublev_vars, i_prev,                     &
                       j_height, j_mean_buoy, j_core_buoy, j_delta_tv,         &
                       j_massflux_d, j_env_tv, j_env_w,                        &
                       j_mean_wex, j_core_wex
use parcel_type_mod, only: n_par, i_massflux_d

use solve_detrainment_mod, only: solve_detrainment
use wind_w_eqn_mod, only: wind_w_eqn
use check_bad_values_mod, only: check_bad_values_cmpr
use raise_error_mod, only: raise_fatal

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of points in compressed parcel super-arrays
! (to avoid repeatedly deallocating and reallocating these,
!  they are dimensioned with the biggest size they will need,
!  which will often be bigger than the number of points here)
integer, intent(in) :: max_points

! Size of the resolved-scale source terms super-array
integer, intent(in) :: n_points_res

! Number of fields (including tracers if applicable)
integer, intent(in) :: n_fields_tot

! Flag for downdrafts versus updrafts
logical, intent(in) :: l_down

! Flag for whether the current model-level is the last allowed
! level (top model-level for updrafts, bottom model-level for
! downdrafts); all mass is forced to detrain if this is true
logical, intent(in) :: l_last_level

! Flag for half-level ascent from half-level to full-level
! (as opposed to full-level to half-level)
logical, intent(in) :: l_to_full_level

! Stuff used to make error messages more informative:
! Structure storing compression indices of the current points
type(cmpr_type), intent(in) :: cmpr
! Current model-level index
integer, intent(in) :: k
! String identifying what sort of draft (updraft, downdraft, etc)
character(len=name_length), intent(in) :: draft_string

! Super-array storing the parcel buoyancies and other properties
! at up to 4 sub-level heights within the current level-step:
! a) Previous model-level interface
! b) Next model-level interface
! c) Height where parcel mean properties first hit saturation
! d) Height where parcel core first hit saturation
real(kind=real_cvprec), intent(in out) :: sublevs                              &
                     ( n_points, n_sublev_vars, max_sublevs )

! Address of next model-level interface in sublevs
integer, intent(in out) :: i_next(n_points)

! Max sub-level height address being used (max value in i_next)
integer, intent(in out) :: i_next_max

! Previous ratio of the parcel core buoyancy over the parcel
! mean buoyancy.  When detrainment occurs, it is calculated
! so-as to approximately preserve this ratio.
real(kind=real_cvprec), intent(in) :: core_mean_ratio(n_points)

! Parcel mean and core vertical velocity drag / s-1
real(kind=real_cvprec), intent(in) :: par_mean_w_drag(n_points)
real(kind=real_cvprec), intent(in) :: par_core_w_drag(n_points)

! Resolved-scale source terms
real(kind=real_cvprec), intent(in out) :: res_source_fields                    &
                                          ( n_points_res, n_fields_tot )

! Parcel mean properties; IN  after lifting and entrainment
!                         OUT after detrainment as well
real(kind=real_cvprec), intent(in out) :: par_next_mean_fields                 &
                                      ( max_points, n_fields_tot )
! Parcel core properties after lifting and entrainment
real(kind=real_cvprec), intent(in out) :: par_next_core_fields                 &
                                      ( max_points, n_fields_tot )

! Super-arrays containing...
! Mass-flux at start of level-step (before entrainment)
real(kind=real_cvprec), intent(in) :: par_prev_super                           &
                                      ( n_points, n_par )
! Mass-flux; IN  after entrainment, but before detrainment
!            OUT after detrainment as well
real(kind=real_cvprec), intent(in out) :: par_conv_super                       &
                                      ( max_points, n_par )

! Output detrained mass and detrained parcel properties
real(kind=real_cvprec), intent(out) :: det_mass_d(n_points)
real(kind=real_cvprec), intent(out) :: det_fields                              &
                                      ( n_points, n_fields_tot )

! Mass-flux interpolated to sub-levels
! (only accounting for entrainment)
real(kind=real_cvprec) :: massflux_1
real(kind=real_cvprec) :: massflux_2
! Interpolation weight used to calculate the above
real(kind=real_cvprec) :: interp

! Non-detrained fraction of mass-flux
real(kind=real_cvprec) :: frac(n_points)

! Position of edge of convection (beyond which air is
! detrained) in dimensionless variable space
real(kind=real_cvprec) :: x_edge(n_points)

! Power of the assumed power-law PDF at which detrainment occured
real(kind=real_cvprec) :: power(n_points)

! Stored values of frac and x_edge from each sub-level step
real(kind=real_cvprec) :: frac_sublevs                                         &
                          ( n_points, max_sublevs )
real(kind=real_cvprec) :: x_edge_sublevs                                       &
                          ( n_points, max_sublevs )

! Integral of parcel mean and core buoyancies over sub-levels
real(kind=real_cvprec) :: mean_buoydz(n_points)
real(kind=real_cvprec) :: core_buoydz(n_points)

! Parcel mean total-water mixing-ratio
real(kind=real_cvprec) :: q_tot(n_points)

! Weight for calculating properties of the detrained air
real(kind=real_cvprec) :: chi_m1(n_points)

! phi_core - phi_mean term in calculation of detrained and
! non-detrained air properties
real(kind=real_cvprec) :: core_m_mean(n_points)

! Store for new value of x_edge after detrainment
real(kind=real_cvprec) :: x_edge_new

! Vertical velocity increment due to buoyancy and drag
real(kind=real_cvprec) :: wind_w_inc

! List of points where doing iterative detrainment solve
integer :: nc
integer :: index_ic(n_points)

! Description of where we are in the code, for error messages
character(len=name_length) :: where_string
! Field name for error messages
character(len=name_length) :: field_name
! Flag for checking positivity (input to bad value checker)
logical :: l_positive

! Numerical tolerance used in safety-check to avoid negative values
real(kind=real_cvprec) :: tolerance

! Numerical tolerance for triggering iterative detrainment calculation
real(kind=real_cvprec), parameter :: safety_thresh                             &
                              = 10.0_real_cvprec * min_delta

! Loop counters
integer :: ic, i_field, i_lev

character(len=*), parameter :: routinename = "SET_DET"


! Check inputs for bad values (NaN, Inf, etc)
if ( i_check_bad_values_cmpr > i_check_bad_none ) then
  ! Check parcel mean properties
  where_string = "On input to set_det call for "                //             &
                  trim(adjustl(draft_string))   // "; "         //             &
                 "par_next_mean_fields"
  do i_field = 1, n_fields_tot
    call check_bad_values_cmpr( cmpr, k, par_next_mean_fields(:,i_field),      &
                                where_string, field_names(i_field),            &
                                field_positive(i_field) )
  end do
end if

if ( l_par_core ) then

  !--------------------------------------------------------------
  ! 1) Adaptive detrainment based on separate parcel mean
  !    and core virtual temperature excess
  !--------------------------------------------------------------

  ! Formula for detrainment is
  !
  ! frac_det = 1 - ( (P+1)/(P+2)
  !                  Tv'_core / ( Tv'_core - Tv'_mean ) )^(P+1)
  !
  ! where P is the power of an assumed power-law PDF of Tv'
  ! within the parcel.  This power is related to the ratio
  ! of the core perturbation to the mean perturbation by:
  !
  ! power = Tv'_core / Tv'_mean - 2
  !
  ! The formula gives:
  ! Tv'_core <= 0                  :    frac_det = 1
  ! Tv'_mean => 1/(P+2) Tv'_core   :    frac_det = 0
  !
  ! i.e.
  ! If the core (the extreme of the pdf) is no longer buoyant,
  ! detrain everything.
  ! If the mean of the pdf is sufficiently closer to the core
  ! than to zero, detrain nothing.

  ! This formula is applied at the end of the current
  ! level-step, and at the height where the parcel first hits saturation
  ! (when this falls between the above heights).
  ! The max detrainment out of all of them is selected.


  ! Check for bad values in inputs
  if ( i_check_bad_values_cmpr > i_check_bad_none ) then

    ! Check parcel core properties
    where_string = "On input to set_det call for "                //           &
                    trim(adjustl(draft_string))   // "; "         //           &
                   "par_next_core_fields"
    do i_field = 1, n_fields_tot
      call check_bad_values_cmpr( cmpr, k, par_next_core_fields(:,i_field),    &
                                  where_string, field_names(i_field),          &
                                  field_positive(i_field) )
    end do

    where_string = "On input to set_det call for "            //               &
                   trim(adjustl(draft_string))

    field_name = "core_mean_ratio"
    l_positive = .true.
    call check_bad_values_cmpr( cmpr, k, core_mean_ratio,                      &
                                where_string,                                  &
                                field_name, l_positive )

  end if

  do ic = 1, n_points
    ! Initialise power: core buoyancy / mean buoyancy = P + 2
    power(ic) = core_mean_ratio(ic) - two
    ! For safety, it is not recommended to let the power of the PDF fall
    ! too close to -1 (this is imposed by not allowing core_mean_ratio
    ! to go below 1+eps, via the parameter min_cmr in comorph_constants_mod).
    ! The PDF is still well-defined with a negative power, but as the power
    ! approaches -1 it implies a very large extrapolation of the parcel
    ! edge properties below the mean, which may be dangerous
    ! (can lead to negative q_vap in the detrained air).

    ! Initialise non-detrained fraction to one
    frac(ic) = one
    x_edge(ic) = one
    frac_sublevs(ic,i_prev) = one
    x_edge_sublevs(ic,i_prev) = one
  end do

  ! Now apply the detrainment formula using the mean and core
  ! buoyancies at up to 3 heights contained in the current
  ! model-level interval:
  !
  ! next_height
  ! saturation height for the parcel mean
  ! saturation height for the parcel core
  !
  ! Out of these 3 heights, the maximum detrainment rate
  ! calculated over all of them is used.


  ! Reverse the sign of the buoyancies for downdrafts, so that
  ! we always consider positive values for non-detrained air
  if ( l_down ) then
    do i_lev = 1, i_next_max
      do ic = 1, n_points
        sublevs(ic,j_mean_buoy,i_lev) = -sublevs(ic,j_mean_buoy,i_lev)
        sublevs(ic,j_core_buoy,i_lev) = -sublevs(ic,j_core_buoy,i_lev)
        sublevs(ic,j_delta_tv,i_lev)  = -sublevs(ic,j_delta_tv,i_lev)
      end do
    end do
  end if

  ! Loop through the sub-level heights.
  ! Note: starting from i_lev=2, as no need to calculate
  ! detrainment at "prev" (the start of this level-step).
  ! TEMPORARY CODE (looping backwards for now to force bit-reproducibility
  ! with previous version, but really need to loop forwards to ensure
  ! all sub-level fields are populated correctly).
  do i_lev = i_next_max, i_prev+1, -1

    ! CODE TEMPORARILY COMMENTED-OUT TO PRESERVE KGO
    ! Force full detrainment for downdrafts hitting the model-bottom;
    ! done here to ensure all the sublev diagnostics calculated after
    ! this loop are consistent with this.
    ! Note: if updrafts hit the model-lid, we just raise a fatal error
    ! at the end of this subroutine.
    !IF ( l_last_level .AND. l_down ) THEN
    !  DO ic = 1, n_points
    !    IF ( i_lev == i_next(ic) ) THEN
    !      frac(ic) = zero
    !      x_edge(ic) = zero
    !    END IF
    !  END DO
    !END IF

    ! Find points where this sub-level is in use, and where
    ! some detrainment should occur there.
    nc = 0
    do ic = 1, n_points
      if ( i_lev <= i_next(ic) .and. frac(ic) > zero ) then
        ! If level in use and not already fully detrained...

        if ( sublevs(ic,j_mean_buoy,i_lev)                                     &
          >= sublevs(ic,j_core_buoy,i_lev) * (one-sqrt_min_delta) ) then
          ! If mean buoyancy exceeds the core buoyancy
          ! (or is within a numerical tolerance of it), then the
          ! full partial detrainment calculations aren't safe,
          ! since they assume core_buoy > mean_buoy.
          ! If accounting for delta_tv implies detrainment
          ! should occur nonetheless, use simple solution
          ! for limit as Tv'_core -> Tv'_mean
          ! (i.e. ignore core_buoy and assume the PDF is a
          !  delta-function with buoyancy mean_buoy)
          if ( sublevs(ic,j_mean_buoy,i_lev) <= zero ) then
            ! Full detrainment when no longer buoyant
            x_edge(ic) = zero
            frac(ic) = zero
          else if ( sublevs(ic,j_mean_buoy,i_lev)                              &
             - frac(ic) * sublevs(ic,j_delta_tv,i_lev) <= zero ) then
            ! Partial detrainment when parcel only loses buoyancy
            ! due to the subsidence term
            frac(ic) = sublevs(ic,j_mean_buoy,i_lev)                           &
                     / sublevs(ic,j_delta_tv,i_lev)
            x_edge(ic) = frac(ic)**(one/(power(ic)+one))
          end if

          ! Otherwise, we do have a well-defined PDF width with
          ! core_buoy > mean_buoy...
        else  ! core_buoy > mean_buoy

          if ( sublevs(ic,j_core_buoy,i_lev) <= zero ) then
            ! If the parcel core is no longer buoyant, we must
            ! have full detrainment; just set x_edge and frac to
            ! zero consistent with this:
            x_edge(ic) = zero
            frac(ic) = zero

          else if ( sublevs(ic,j_delta_tv,i_lev)                               &
                    <= safety_thresh * sublevs(ic,j_core_buoy,i_lev) ) then
            ! If the subsidence increment term is not positive
            ! (or is negligible compared to the buoyancy),
            ! just use the explicit solution, ignoring the subsidence term.
            ! Without that term, we have:
            ! Tv'(x) = Tv'_core
            !        - (p+2)/(p+1) (Tv'_core - Tv'_mean) x
            !        = 0
            ! => x = (p+1)/(p+2) Tv'_core / (Tv'_core - Tv'_mean)
            x_edge_new = ((power(ic)+one)/(power(ic)+two))                     &
                       * sublevs(ic,j_core_buoy,i_lev)                         &
                       / ( sublevs(ic,j_core_buoy,i_lev)                       &
                         - sublevs(ic,j_mean_buoy,i_lev) )
            ! Use new edge only if it implies more detrainment than
            ! we've already got (x_edge_new < x_edge)
            if ( x_edge_new < x_edge(ic) ) then
              x_edge(ic) = x_edge_new
              frac(ic)   = x_edge_new**(power(ic)+one)
            end if

            ! Otherwise, partial detrainment will occur if the
            ! buoyancy at the current edge is negative, and needs
            ! to be solved iteratively to account for the subsidence term...
            ! We have:
            ! Tv'(x) = Tv'_core
            !        - (p+2)/(p+1) (Tv'_core - Tv'_mean) x
            !        - frac dTv_k
            ! Detrain when this is < 0:
            ! (actually only detrain if it goes negative by
            !  a small numerical tolerance, for safety;
            !  otherwise, just due to rounding errors, the 2
            !  initial guesses in solve_detrainment can fail to
            !  bracket the root).
          else if ( sublevs(ic,j_core_buoy,i_lev)                              &
                  - ((power(ic)+two)/(power(ic)+one))                          &
                    * ( sublevs(ic,j_core_buoy,i_lev)                          &
                      - sublevs(ic,j_mean_buoy,i_lev) ) * x_edge(ic)           &
                  - frac(ic) * sublevs(ic,j_delta_tv,i_lev)                    &
                < -safety_thresh * sublevs(ic,j_core_buoy,i_lev) ) then
            ! Store indices of points where doing iterative
            ! partial detrainment
            nc = nc + 1
            index_ic(nc) = ic
          end if

        end if  ! core_buoy > mean_buoy

      end if  ! ( i_lev <= i_next(ic) .AND. frac(ic) > zero )
    end do  ! ic = 1, n_points

    if ( nc > 0 ) then
      ! If any points need iterative detrainment calculation
      ! at this sub-level...

      ! Call routine to solve for the non-detrained fraction
      call solve_detrainment( n_points, nc, index_ic,                          &
                              sublevs(:,j_mean_buoy,i_lev),                    &
                              sublevs(:,j_core_buoy,i_lev),                    &
                              power, sublevs(:,j_delta_tv,i_lev),              &
                              x_edge, frac )

    end if

    do ic = 1, n_points

      ! Safety check; force frac to be zero (full detrainment) if it is too
      ! close to zero to numerically represent the non-detrained parcel
      if ( frac(ic) < min_delta .or.                                           &
           frac(ic) * par_conv_super(ic,i_massflux_d) < min_float ) then
        x_edge(ic) = zero
        frac(ic)   = zero
      end if

      ! Store non-detrained fraction and PDF edge at sub-level steps
      frac_sublevs(ic,i_lev)   = frac(ic)
      x_edge_sublevs(ic,i_lev) = x_edge(ic)

    end do  ! ic = 1, n_points

  end do  ! i_lev == i_prev+1, i_next_max


  ! Restore the sign of the buoyancies for downdrafts
  if ( l_down ) then
    do i_lev = 1, i_next_max
      do ic = 1, n_points
        sublevs(ic,j_mean_buoy,i_lev) = -sublevs(ic,j_mean_buoy,i_lev)
        sublevs(ic,j_core_buoy,i_lev) = -sublevs(ic,j_core_buoy,i_lev)
        sublevs(ic,j_delta_tv,i_lev)  = -sublevs(ic,j_delta_tv,i_lev)
      end do
    end do
  end if

  ! Check for bad values in x_edge and frac
  if ( i_check_bad_values_cmpr > i_check_bad_none ) then
    where_string = "set_det call for "                        //               &
                   trim(adjustl(draft_string))

    field_name = "x_edge"
    l_positive = .true.
    call check_bad_values_cmpr( cmpr, k, x_edge, where_string,                 &
                                field_name, l_positive )

    field_name = "frac (non-detrained fraction)"
    l_positive = .true.
    call check_bad_values_cmpr( cmpr, k, frac, where_string,                   &
                                field_name, l_positive )

  end if


  ! Update sub-level fields...

  do i_lev = i_prev+1, i_next_max
    do ic = 1, n_points
      ! Store mass-fluxes after detrainment
      sublevs(ic,j_massflux_d,i_lev) = sublevs(ic,j_massflux_d,i_lev)          &
                                     * frac_sublevs(ic,i_lev)

      ! Calculate the final solved env Tv subsidence increment
      sublevs(ic,j_delta_tv,i_lev) = sublevs(ic,j_delta_tv,i_lev)              &
                                   * frac_sublevs(ic,i_lev)
    end do
  end do

  ! Update parcel vertical velocities on sub-levels, integrating
  ! the buoyancies accounting for the subsidence increment to environment
  ! virtual temperature calculated above
  do ic = 1, n_points
    mean_buoydz(ic) = zero
    core_buoydz(ic) = zero
  end do
  do i_lev = i_prev+1, i_next_max
    ! Mean w-excess
    call wind_w_eqn( n_points,                                                 &
           sublevs(:,j_height,i_prev),                                         &
           sublevs(:,j_height,i_lev-1),    sublevs(:,j_height,i_lev),          &
           sublevs(:,j_mean_buoy,i_lev-1), sublevs(:,j_mean_buoy,i_lev),       &
           sublevs(:,j_delta_tv,i_lev-1),  sublevs(:,j_delta_tv,i_lev),        &
           sublevs(:,j_env_tv,i_lev-1),    sublevs(:,j_env_tv,i_lev),          &
           par_mean_w_drag, mean_buoydz, sublevs(:,j_mean_wex,i_lev) )
    ! Core w-excess
    call wind_w_eqn( n_points,                                                 &
           sublevs(:,j_height,i_prev),                                         &
           sublevs(:,j_height,i_lev-1),    sublevs(:,j_height,i_lev),          &
           sublevs(:,j_core_buoy,i_lev-1), sublevs(:,j_core_buoy,i_lev),       &
           sublevs(:,j_delta_tv,i_lev-1),  sublevs(:,j_delta_tv,i_lev),        &
           sublevs(:,j_env_tv,i_lev-1),    sublevs(:,j_env_tv,i_lev),          &
           par_core_w_drag, core_buoydz, sublevs(:,j_core_wex,i_lev) )
  end do

  ! Calculate parcel total-water mixing ratio
  ! (note q_cl currently stores q_vap+q_cl, so we don't separately add q_vap).
  do ic = 1, n_points
    q_tot(ic) = zero
  end do
  do i_field = i_qc_first, i_qc_last
    do ic = 1, n_points
      q_tot(ic) = q_tot(ic) + par_next_mean_fields(ic,i_field)
    end do
  end do

  ! Update parcel vertical velocities and resolved-scale source terms
  do ic = 1, n_points

    ! Store parcel mean w before buoyancy and drag
    wind_w_inc = par_next_mean_fields(ic,i_wind_w)

    ! Convert latest w-excesses to absolute momentum
    ! (parcel fields are currently in conserved variable form, so
    !  wind_w actually stores momentum per unit dry-mass = wind_w (1 + q_tot) )
    par_next_mean_fields(ic,i_wind_w) = ( sublevs(ic,j_env_w,i_next(ic))       &
                                        + sublevs(ic,j_mean_wex,i_next(ic)) )  &
                                      * ( one + q_tot(ic) )
    par_next_core_fields(ic,i_wind_w) = ( sublevs(ic,j_env_w,i_next(ic))       &
                                        + sublevs(ic,j_core_wex,i_next(ic)) )  &
                                      * ( one + q_tot(ic) )

    ! Complete calculation of parcel momentum increment from buoyancy and drag
    wind_w_inc = par_next_mean_fields(ic,i_wind_w) - wind_w_inc

    ! Add contribution to the resolved momentum source terms due to the
    ! reaction force from the buoyancy and drag on the parcel
    res_source_fields(ic,i_wind_w) = res_source_fields(ic,i_wind_w)            &
                              - wind_w_inc * par_conv_super(ic,i_massflux_d)

  end do

  ! Adjust the mean buoyancy and w-excess to account for removal of the
  ! negatively buoyant detrained air from the parcel
  ! ( core_buoy - mean_buoy_new ) = x_edge ( core_buoy - mean_buoy_old )
  ! => mean_buoy_new = core_buoy - x_edge ( core_buoy - mean_buoy_old )
  !                  = mean_buoy_old + ( 1 - x_edge )
  !                                    ( core_buoy - mean_buoy_old )
  do ic = 1, n_points
    do i_lev = i_prev+1, i_next(ic)
      ! CODE TEMPORARILY COMMENTED OUT TO PRESERVE KGO; TO BE ADDED SOON:
      ! Calculation omitted for buoyancy on sub-levels; this is an oversight
      ! and impacts the convective cloud amounts and CAPE contributions computed
      ! from the buoyancy at saturation height.  To be fixed soon (KGO change).
      !sublevs(ic,j_mean_buoy,i_lev) = sublevs(ic,j_mean_buoy,i_lev)           &
      !                              + ( one - x_edge_sublevs(ic,i_lev) )      &
      !  * ( sublevs(ic,j_core_buoy,i_lev)  - sublevs(ic,j_mean_buoy,i_lev) )
      sublevs(ic,j_mean_wex,i_lev)  = sublevs(ic,j_mean_wex,i_lev)             &
                                    + ( one - x_edge_sublevs(ic,i_lev) )       &
        * ( sublevs(ic,j_core_wex,i_lev)   - sublevs(ic,j_mean_wex,i_lev)  )
    end do  ! i_lev = i_prev+1, i_next(ic)
  end do  ! ic = 1, n_points


  ! x_edge now stores
  ! (P+1)/(P+2) Tv'_core / ( Tv'_core - Tv'_mean )
  !
  ! which is equal to
  ! (P+1)/(P+2) ( phi_core - phi_crit ) / ( phi_core - phi_mean )
  !  = ( phi_core - phi_crit ) / ( phi_core - phi_edge )
  !
  ! where phi_crit is the value of each field phi below-which
  ! air is detrained, and above-which air is retained.
  !
  ! The Cumulative Distribution Function is given by
  ! CDF = 1 - frac = 1 - x_edge^(P+1)
  ! This sets the fraction of mass to detrain...

  ! Set detrained and remaining non-detrained air properties
  ! based on the mean over their respective parts of the
  ! assumed-PDF...

  ! The mean properties of the non-detrained air corresponds
  ! to the integral from x = 0 to x_edge, which yields:
  !
  ! phi_nd = phi_core - x_edge ( phi_core - phi_mean )
  !
  ! To ensure we get exactly the right answer in the limit of
  ! very small detrainment (x_edge -> 1), we rearrange this to:
  !
  ! phi_nd = phi_mean + (1-x_edge) ( phi_core - phi_mean )

  ! The mean properties of the detrained air corresponds
  ! to the integral from x = x_edge to 1, which yields:
  !
  ! phi_det = phi_core - ( 1 - x_edge^(p+2) ) / ( 1 - x_edge^(p+1) )
  !                      ( phi_core - phi_mean )
  !
  ! Labelling the term chi = ( 1 - x_edge^(p+2) ) / ( 1 - x_edge^(p+1) ),
  ! we rearrange this to:
  !
  ! phi_det = phi_mean - (chi-1) ( phi_core - phi_mean )
  !
  ! These formulae ensure that:
  ! a) In the limit of full detrainment (x_edge=0, chi=1):
  !    - the non-detrained air has the parcel core properties.
  !    - the detrained air has the parcel mean properties.
  ! b) In limit of no detrainment (x_edge=1, chi=(P+2)/(P+1)):
  !    - the non-detrained air has the parcel mean properties.
  !    - the detrained air has the edge properties implied by
  !      the assumed power-law PDF:
  !      ( phi_mean - phi_edge ) = ( phi_core - phi_mean ) / (P+1)
  ! c) The mass-weighted mean of the detrained and non-detrained
  !    properties equals the parcel mean before detrainment:
  !
  !      x_edge^(p+1) phi_nd + ( 1 - x_edge^(p+1) ) phi_det = phi_mean

  ! Next, we calculate chi-1...
  !
  ! chi-1 = ( 1 - x_edge^(p+2) ) / ( 1 - x_edge^(p+1) ) - 1
  !       = x_edge^(p+1) ( 1 - x_edge ) / ( 1 - x_edge^(p+1) )
  !       = frac ( 1 - x_edge ) / ( 1 - frac )
  !
  ! Note that, while this formula does have a finite value in the
  ! limit x_edge -> 1, calculating it using this would yield a
  ! div-by-zero in that limit.  Therefore, we use a Taylor expansion
  ! of the above in the event that x_edge is near 1.

  ! Calculate chi-1
  do ic = 1, n_points
    if ( one-x_edge(ic) < sqrt_min_delta ) then
      ! If x_edge is too close to 1.0, use Taylor expansion:
      ! Let y = 1 - x_edge
      ! => chi-1 = (1-y)^(p+1) y / ( 1 - (1-y)^(p+1) )
      !
      ! Applying l'Hopital's rule, the limit as y -> 0
      ! is found to be
      ! chi-1(y=0) = 1 / (p+1)                     (1)
      !
      ! Differentiating chi-1 w.r.t. y and rearranging:
      !
      ! d/dy(chi-1) = (1-y)^p ( -(p+2) y - (1-y)^(p+2) + 1 )
      !             / ( 1 - (1-y)^(p+1) )^2
      !
      ! Applying l'Hopital's rule again (requires 2 differentiations
      ! to get non-zero numerator and denominator), the limit as y -> 0
      ! is found to be:
      ! d/dy(chi-1)(y=0) = -1/2 (p+2)/(p+1)        (2)
      !
      ! Combining (1) and (2) and factorising, a 1st order Taylor
      ! expansion in y for chi-1 is:
      ! chi-1(y) ~ 1/(p+1) ( 1 - 1/2 (p+2) y )

      chi_m1(ic) = (one/(power(ic)+one))                                       &
                 * ( one - half*(power(ic)+two) * (one-x_edge(ic)) )
    else
      ! Otherwise use full formula
      chi_m1(ic) = frac(ic) * ( one - x_edge(ic) ) / ( one - frac(ic) )
    end if
  end do

  ! Loop over fields
  do i_field = 1, n_fields_tot

    ! Precalculate the term phi_core - phi_mean, which appears
    ! in the formula for both the detrained and non-detrained properties
    do ic = 1, n_points
      core_m_mean(ic) = ( par_next_core_fields(ic,i_field)                     &
                        - par_next_mean_fields(ic,i_field) )
    end do

    if ( field_positive(i_field) ) then
      ! For fields that aren't supposed to have negative values...
      do ic = 1, n_points
        ! Impose a limit on the core minus mean term to ensure
        ! the edge of the PDF (and potentially the detrained air)
        ! doesn't go negative.
        ! The core and mean values should already be positive on input
        ! to this routine, but the edge properties are extrapolated
        ! from these and so could go negative without this check.
        ! We have:
        !
        ! phi_edge = phi_mean - ( (p+2)/(p+1) - 1 ) ( phi_core - phi_mean )
        !          = phi_mean - 1/(p+1) ( phi_core - phi_mean )
        !
        ! Therefore, to ensure phi_edge >= 0, we must have:
        !
        ! ( phi_core - phi_mean ) <= (p+1) phi_mean

        ! Set a numerical tolerance to avoid still getting negative edge
        ! values due to rounding errors in the subsequent calculations.
        tolerance = safety_thresh
        ! Note safety_thresh is set assuming the precision of the variables
        ! scales with EPSILON (min_delta) times the value.  However, some
        ! compilers allow values less than TINY (min_float) to be stored
        ! with lower precision than this.  Therefore, if the value of the
        ! field is less than min_float (but not zero), we need to scale the
        ! tolerance up in this check.
        if ( par_next_mean_fields(ic,i_field) < min_float .and.                &
             par_next_mean_fields(ic,i_field) > zero ) then
          tolerance = min( tolerance * ( min_float                             &
                                       / par_next_mean_fields(ic,i_field) ),   &
                           one )
        end if

        ! Apply limit to phi_core - phi_mean, with numerical tolerance
        core_m_mean(ic) = min( core_m_mean(ic),                                &
                               par_next_mean_fields(ic,i_field)                &
                               * (power(ic)+one) * (one-tolerance) )

      end do
    end if  ! ( field_positive(i_field) )

    ! Set detrained air properties:
    do ic = 1, n_points
      ! phi_det = phi_mean - (chi-1) ( phi_core - phi_mean )
      det_fields(ic,i_field) = par_next_mean_fields(ic,i_field)                &
                             - chi_m1(ic) * core_m_mean(ic)
    end do
    ! Set non-detrained air properties:
    do ic = 1, n_points
      ! phi_nd = phi_mean + (1-x_edge) ( phi_core - phi_mean )
      par_next_mean_fields(ic,i_field) = par_next_mean_fields(ic,i_field)      &
                                       + (one - x_edge(ic)) * core_m_mean(ic)
    end do

  end do  ! i_field = 1, n_fields_tot


else  ! ( l_par_core )


  !--------------------------------------------------------------
  ! 2) Detrainment rate just based on fractional rate of change
  !    of parcel mean virtual temperature excess with height
  !--------------------------------------------------------------

  ! Set the properties of detrained air;
  ! currently just copies of the mean parcel properties:
  do i_field = 1, n_fields_tot
    do ic = 1, n_points
      det_fields(ic,i_field) = par_next_mean_fields(ic,i_field)
    end do
  end do
  ! Without the parcel core, the parcel properties are assumed
  ! to be homogeneous, so detrainment doesn't modify the parcel
  ! mean fields

  ! Non-detrained fraction of mass is parameterised as:
  !
  ! frac_nondet = (M_prev/M_next) (Tv'_next/Tv'_prev)
  !
  ! where M_next is the mass-flux at the next level accounting
  ! for entrainment but not detrainment.
  !
  ! This is designed to yield:
  ! - Full detrainment if Tv'_next goes to zero
  ! - Detrainment = Entrainment if Tv'_next = Tv'_prev
  !
  ! The calculation below is complicated by the fact that we
  ! account for kinks in the profile between _prev and _next
  ! due to estimated parcel properties at the saturation level
  ! (if the parcel just started releasing latent heat between
  ! prev and next)...

  ! Initialise non-detrained fraction to unity
  do ic = 1, n_points
    frac(ic) = one
  end do

  ! Loop over the sub-levels stored in sublevs
  do i_lev = 2, i_next_max

    do ic = 1, n_points
      ! If the current point / sub-level is in use
      if ( i_lev <= i_next(ic) ) then

        ! If no longer buoyant, detrain all
        if ( ( sublevs(ic,j_mean_buoy,i_lev-1) <= zero .and. (.not. l_down) )  &
        .or. ( sublevs(ic,j_mean_buoy,i_lev-1) >= zero .and. l_down ) ) then
          frac(ic) = zero
        else
          ! Find mass-flux at base and top of the current
          ! sub-level, accounting for entrainment only
          interp = (sublevs(ic,j_height,i_lev-1) - sublevs(ic,j_height,i_prev))&
              / (sublevs(ic,j_height,i_next(ic)) - sublevs(ic,j_height,i_prev))
          massflux_1 = (one-interp) * par_prev_super(ic,i_massflux_d)          &
                     +      interp  * par_conv_super(ic,i_massflux_d)
          interp = (sublevs(ic,j_height,i_lev) - sublevs(ic,j_height,i_prev))  &
            / (sublevs(ic,j_height,i_next(ic)) - sublevs(ic,j_height,i_prev))
          massflux_2 = (one-interp) * par_prev_super(ic,i_massflux_d)          &
                     +      interp  * par_conv_super(ic,i_massflux_d)

          ! Scale frac by detrainment in current sub-level
          frac(ic) = frac(ic) * min( one, max( zero,                           &
                      (massflux_1/massflux_2)                                  &
                    * ( sublevs(ic,j_mean_buoy,i_lev)                          &
                      / sublevs(ic,j_mean_buoy,i_lev-1) )                      &
                                             ) )
        end if

      end if  ! ( i_lev <= i_next(ic) )
    end do  ! ic = 1, n_points

  end do  ! i_lev = 2, i_next_max

end if  ! ( l_par_core )


! Set detrained mass:
do ic = 1, n_points
  ! Scale by mass-flux to get detrained mass
  det_mass_d(ic) = ( one - frac(ic) ) * par_conv_super(ic,i_massflux_d)
end do

! Subtract detrained mass from the mass-flux!
do ic = 1, n_points
  par_conv_super(ic,i_massflux_d) = par_conv_super(ic,i_massflux_d)            &
                                  - det_mass_d(ic)
end do


if ( l_last_level .and. (.not. l_to_full_level) ) then
  ! If this is the last level (and the 2nd half-level step)...
  if ( l_down ) then

    ! For downdrafts hitting the surface, force them to fully detrain here.
    do ic = 1, n_points
      if ( par_conv_super(ic,i_massflux_d) > zero ) then
        ! Set detrained mass equal to total massflux and reset massflux to zero
        det_mass_d(ic) = det_mass_d(ic) + par_conv_super(ic,i_massflux_d)
        par_conv_super(ic,i_massflux_d) = zero
        ! Set detrained air the same as the original in-parcel mean
        do i_field = 1, n_fields_tot
          det_fields(ic,i_field)                                               &
            = (one-frac(ic)) * det_fields(ic,i_field)                          &
            +      frac(ic)  * par_next_mean_fields(ic,i_field)
        end do
      end if
    end do

  else

    ! If this is an updraft, this implies convection has gone to
    ! the top of the model, so raise a fatal error.
    do ic = 1, n_points
      if ( par_conv_super(ic,i_massflux_d) > zero ) then
        ! Note it is expected that downdrafts will sometimes hit the
        ! surface, but updrafts hitting the lid is probably not OK.
        call raise_fatal( routinename,                                         &
             "Convective updraft trying to go beyond the "                  // &
             "uppermost convection level."                         //newline// &
             "Either something has gone very wrong, or the "                // &
             "highest allowed convection level has just "          //newline// &
             "been set too low (try increasing k_top_conv)" )
      end if
    end do

  end if
end if  ! ( l_last_level .AND. (.NOT. l_to_full_level) )


! Check outputs for bad values (NaN, Inf, etc)
if ( i_check_bad_values_cmpr > i_check_bad_none ) then

  ! Check parcel mean properties
  where_string = "On output from set_det call for "             //             &
                  trim(adjustl(draft_string))   // "; "         //             &
                 "par_next_mean_fields"
  do i_field = 1, n_fields_tot
    call check_bad_values_cmpr( cmpr, k, par_next_mean_fields(:,i_field),      &
                                where_string, field_names(i_field),            &
                                field_positive(i_field) )
  end do

  ! Check detrained properties
  where_string = "On output from set_det call for "             //             &
                  trim(adjustl(draft_string))   // "; "         //             &
                 "det_fields"
  do i_field = 1, n_fields_tot
    call check_bad_values_cmpr( cmpr, k, det_fields(:,i_field),                &
                                where_string, field_names(i_field),            &
                                field_positive(i_field) )
  end do

  ! Check detrained mass
  field_name = "det_mass_d"
  l_positive = .true.
  call check_bad_values_cmpr( cmpr, k, det_mass_d, where_string,               &
                              field_name, l_positive )
  ! Check mass-flux
  field_name = "next_massflux_d"
  l_positive = .true.
  call check_bad_values_cmpr( cmpr, k, par_conv_super(:,i_massflux_d),         &
                              where_string,                                    &
                              field_name, l_positive )

end if


return
end subroutine set_det


end module set_det_mod
