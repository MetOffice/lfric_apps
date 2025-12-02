! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module fall_in_mod

implicit none

contains


! Routine to compute fall-in of hydrometeors from the environment
! into the convective parcel
! (converts the local fall-flux into an increase in parcel
!  mixing ratio).
! Calculation of fall-out is done later by a separate routine.
! Also computes the change in parcel winds and temperature
! due to momentum and heat transported into the parcel by the
! falling-in hydrometeors.
subroutine fall_in( n_points, nc, index_ic,                                    &
                    cp_cond, dt_over_rhod_lz, flux_cond,                       &
                    fallin_wind_u, fallin_wind_v, fallin_wind_w,               &
                    fallin_temp,                                               &
                    cp_tot, q_tot, wind_u, wind_v, wind_w,                     &
                    temperature, q_cond,                                       &
                    l_diags, i_cond, moist_proc_diags,                         &
                    n_points_diag, n_diags, diags_super )

use comorph_constants_mod, only: real_cvprec, indi_thresh, one
use moist_proc_diags_type_mod, only: moist_proc_diags_type

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of points where fall-in flux is nonzero
integer, intent(in) :: nc

! Indices of points where fall-in flux is nonzero
integer, intent(in) :: index_ic(nc)

! Specific heat capacity of the current condensed water species
real(kind=real_cvprec), intent(in) :: cp_cond

! Factor delta_t / ( rho_dry vert_len ) used to convert fall-flux
! to mixing ratio increment
real(kind=real_cvprec), intent(in) :: dt_over_rhod_lz(n_points)

! Inward flux of the current hydrometeor species / kg m-2 s-1
real(kind=real_cvprec), intent(in out) :: flux_cond(n_points)
! (this is an input, but needs intent(inout) as we convert it
!  to an increment)

! Winds and temperature of the air from which the falling-in
! hydrometeros have fallen
real(kind=real_cvprec), intent(in) :: fallin_wind_u(n_points)
real(kind=real_cvprec), intent(in) :: fallin_wind_v(n_points)
real(kind=real_cvprec), intent(in) :: fallin_wind_w(n_points)
real(kind=real_cvprec), intent(in) :: fallin_temp(n_points)

! Total heat capacity per unit dry-mass of the air
! (modified by the addition of falling-in precip)
real(kind=real_cvprec), intent(in out) :: cp_tot(n_points)

! Total-water mixing-ratio = q_vap + sum q_cond
! (used to calc change in winds due to precip momentum transfer)
real(kind=real_cvprec), intent(in out) :: q_tot(n_points)

! Wind velocity / ms-1
real(kind=real_cvprec), intent(in out) :: wind_u(n_points)
real(kind=real_cvprec), intent(in out) :: wind_v(n_points)
real(kind=real_cvprec), intent(in out) :: wind_w(n_points)

! Temperature
real(kind=real_cvprec), intent(in out) :: temperature(n_points)

! Mixing ratio of the falling-in hydrometeor species
real(kind=real_cvprec), intent(in out) :: q_cond(n_points)

! Master switch for whether or not to calculate any diagnostics
logical, intent(in) :: l_diags
! Super-array address of current species,
! needed to address the correct diagnostic array field
integer, intent(in) :: i_cond
! Structure containing diagnostic flags etc.
type(moist_proc_diags_type), intent(in) :: moist_proc_diags
! Number of points in the diagnostics super-array
integer, intent(in) :: n_points_diag
! Total number of diagnostics in the super-array
integer, intent(in) :: n_diags
! Super-array to store all the diagnostics
real(kind=real_cvprec), intent(in out) :: diags_super                          &
                                         ( n_points_diag, n_diags )


! Loop counter
integer :: ic, ic2, i_super


! If non-zero fall-in flux at most points, just do the
! calculations at all points
if ( real(nc,real_cvprec) > indi_thresh * real(n_points,real_cvprec) ) then

  ! Convert fall-in flux to a local mixing-ratio increment
  do ic = 1, n_points
    flux_cond(ic) = flux_cond(ic) * dt_over_rhod_lz(ic)
  end do

  ! Increment hydrometeor mixing ratio
  do ic = 1, n_points
    q_cond(ic) = q_cond(ic) + flux_cond(ic)
  end do

  ! Increment heat capacity and total-water mixing-ratio
  do ic = 1, n_points
    cp_tot(ic) = cp_tot(ic) + cp_cond * flux_cond(ic)
  end do
  do ic = 1, n_points
    q_tot(ic) = q_tot(ic) + flux_cond(ic)
  end do

  ! Increment winds
  !
  ! new_wt_mass = old_wet_mass + fallin_mass
  !
  ! new_wind = ( old_wet_mass * old_wind
  !            + fallin_mass * fallin_wind ) / new_wet_mass
  !
  !          = ( new_wet_mass * old_wind
  !            + fallin_mass * ( fallin_wind
  !                            - old_wind ) ) / new_wet_mass
  !
  !          = old_wind + fallin_mass
  !              * ( fallin_wind - old_wind ) / new_wet_mass
  !
  ! Here we've normalised by parcel dry-mass, so replace
  ! fallin_mass with fall-in mixing ratio increment,
  ! and replace new_wet_mass by wet-mass/dry-mass = 1 + q_tot:
  !
  do ic = 1, n_points
    wind_u(ic) = wind_u(ic) + flux_cond(ic)                                    &
         * ( fallin_wind_u(ic) - wind_u(ic) ) / (one+q_tot(ic))
  end do
  do ic = 1, n_points
    wind_v(ic) = wind_v(ic) + flux_cond(ic)                                    &
         * ( fallin_wind_v(ic) - wind_v(ic) ) / (one+q_tot(ic))
  end do
  do ic = 1, n_points
    wind_w(ic) = wind_w(ic) + flux_cond(ic)                                    &
         * ( fallin_wind_w(ic) - wind_w(ic) ) / (one+q_tot(ic))
  end do

  ! Increment temperature
  do ic = 1, n_points
    temperature(ic) = temperature(ic) + flux_cond(ic) *cp_cond                 &
         * ( fallin_temp(ic) - temperature(ic) ) / cp_tot(ic)
  end do

  ! Save diagnostic of increment due to sedimentation
  if ( l_diags ) then
    if ( moist_proc_diags % diags_cond(i_cond)%pt                              &
         % dq_fall % flag ) then
      ! Extract super-array address
      i_super = moist_proc_diags % diags_cond(i_cond)%pt                       &
                % dq_fall % i_super
      ! Save fall-in increment stored in flux_cond
      do ic = 1, n_points
        diags_super(ic,i_super) = flux_cond(ic)
      end do
    end if
  end if


  ! If non-zero fall-in flux at minority of points, use the
  ! stored indices to do calculations at only those points
else !( real(nc,real_cvprec) > indi_thresh * real(n_points,real_cvprec) )

  ! Loop over only points with nonzero fall-in flux
  do ic2 = 1, nc
    ic = index_ic(ic2)

    ! Do exactly the same calculations as above...

    flux_cond(ic) = flux_cond(ic) * dt_over_rhod_lz(ic)

    q_cond(ic) = q_cond(ic) + flux_cond(ic)

    cp_tot(ic) = cp_tot(ic) + cp_cond * flux_cond(ic)
    q_tot(ic) = q_tot(ic) + flux_cond(ic)

    wind_u(ic) = wind_u(ic) + flux_cond(ic)                                    &
         * ( fallin_wind_u(ic) - wind_u(ic) ) / (one+q_tot(ic))
    wind_v(ic) = wind_v(ic) + flux_cond(ic)                                    &
         * ( fallin_wind_v(ic) - wind_v(ic) ) / (one+q_tot(ic))
    wind_w(ic) = wind_w(ic) + flux_cond(ic)                                    &
         * ( fallin_wind_w(ic) - wind_w(ic) ) / (one+q_tot(ic))

    temperature(ic) = temperature(ic) + flux_cond(ic) *cp_cond                 &
         * ( fallin_temp(ic) - temperature(ic) ) / cp_tot(ic)

  end do

  if ( l_diags ) then
    if ( moist_proc_diags % diags_cond(i_cond)%pt                              &
         % dq_fall % flag ) then
      i_super = moist_proc_diags % diags_cond(i_cond)%pt                       &
                % dq_fall % i_super
      do ic2 = 1, nc
        ic = index_ic(ic2)
        diags_super(ic,i_super) = flux_cond(ic)
      end do
    end if
  end if

end if !( real(nc,real_cvprec) > indi_thresh * real(n_points,real_cvprec) )


return
end subroutine fall_in

end module fall_in_mod
