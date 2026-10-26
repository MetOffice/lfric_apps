! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module ice_nucleation_mod

implicit none

contains

! Subroutine to do ice nucleation (spontaneous freezing of
! liquid water into ice)
! Check for any liquid condensed water occuring
! below the homogeneous freezing threshold, and freezes it
! all out instantly if present.
! Then applies a small freezing rate due to heterogeneous
! nucleation if between the heterogeneous nucleation
! threshold and the homogeneous threshold
subroutine ice_nucleation( n_points,                                           &
                           nc_liq, index_ic_liq,                               &
                           nc_ice, index_ic_ice,                               &
                           ref_temp, q_liq, q_ice,                             &
                           temperature, cp_tot,                                &
                           dq_frz_tot, l_diags,                                &
                           i_liq, i_ice, moist_proc_diags,                     &
                           n_points_diag, n_diags, diags_super )

use comorph_constants_mod, only: real_cvprec, zero,                            &
                                 homnuc_temp, hetnuc_temp,                     &
                                 q_activate
use moist_proc_diags_type_mod, only: moist_proc_diags_type

use lat_heat_mod, only: lat_heat_incr, i_phase_change_frz

implicit none

! Number of points
integer, intent(in) :: n_points

! Points where the liquid species exists
integer, intent(in out) :: nc_liq
integer, intent(in out) :: index_ic_liq(n_points)
! (this list gets altered if liquid is completely removed
!  from any points by freezing)

! Points where the ice species exists
integer, intent(in out) :: nc_ice
integer, intent(in out) :: index_ic_ice(n_points)
! (this list get altered if new ice is formed at a point where
!  there wasn't any before)

! Reference temperature used freezing threshold check
real(kind=real_cvprec), intent(in) :: ref_temp(n_points)

! Mixing ratio of a liquid water species
real(kind=real_cvprec), intent(in out) :: q_liq(n_points)

! Mxing ratio of the ice species which q_liq freezes into
real(kind=real_cvprec), intent(in out) :: q_ice(n_points)

! Air temperature to be updated with latent heat from freezing
real(kind=real_cvprec), intent(in out) :: temperature(n_points)
! Total heat capacity per unit dry-mass
real(kind=real_cvprec), intent(in out) :: cp_tot(n_points)

! Total rate of freezing onto the ice species
! (may already include contributions from freezing of other
!  liquid species)
real(kind=real_cvprec), intent(in out) :: dq_frz_tot(n_points)

! Master switch for whether or not to calculate any diagnostics
logical, intent(in) :: l_diags
! Super-array addresses of the liquid and ice species,
! needed to address the correct diagnostic array fields
integer, intent(in) :: i_liq, i_ice
! Structure containing diagnostic flags etc.
type(moist_proc_diags_type), intent(in) :: moist_proc_diags
! Number of points in the diagnostics super-array
integer, intent(in) :: n_points_diag
! Total number of diagnostics in the super-array
integer, intent(in) :: n_diags
! Super-array to store all the diagnostics
real(kind=real_cvprec), intent(in out) :: diags_super                          &
                                         ( n_points_diag, n_diags )

! Amount of mixing ratio to be frozen by this routine
real(kind=real_cvprec) :: dq_frz(n_points)

! Number of points (and their indices) where freezing occurs
integer :: nc_frz
integer :: index_ic_frz(n_points)

! Flag for whether all liquid has been frozen at any points
logical :: l_full_frz
! Flag for whether new ice added where none yet exists
logical :: l_frz_where_no_ice
! Temporary store for number of points when rejigging list
integer :: nc_tmp

! Loop counters
integer :: ic, ic2, i_super


! Initialise flag for fully removing liquid from any grid-points
l_full_frz = .false.
! Initialise flag for adding ice at new grid-points
l_frz_where_no_ice = .false.

! First, do homogeneous freezing...

! Find points where reference temperature below homogeneous freezing threshold
nc_frz = 0
do ic2 = 1, nc_liq
  ic = index_ic_liq(ic2)
  ! Initialise freezing increment to zero at all liquid points
  dq_frz(ic) = zero
  if ( ref_temp(ic) < homnuc_temp ) then
    nc_frz = nc_frz + 1
    index_ic_frz(nc_frz) = ic
    ! Set flag where about to add ice at a grid-point that doesn't have any
    if ( .not. q_ice(ic) > zero )  l_frz_where_no_ice = .true.
  end if
end do

if ( nc_frz > 0 ) then
  ! If any points below homnuc_temp...

  do ic2 = 1, nc_frz
    ic = index_ic_frz(ic2)
    ! Store increment, and transfer all liquid to ice
    dq_frz(ic) = q_liq(ic)
    q_ice(ic) = q_ice(ic) + dq_frz(ic)
    q_liq(ic) = zero
    ! Set flag to indicate total removal of liquid
    l_full_frz = .true.
    ! Increment total rate if freezing onto the ice
    dq_frz_tot(ic) = dq_frz_tot(ic) + dq_frz(ic)
  end do

  ! Increment air temperature and heat capacity
  ! with latent heat of fusion
  call lat_heat_incr( n_points, nc_frz, i_phase_change_frz,                    &
                      cp_tot, temperature,                                     &
                      index_ic=index_ic_frz, dq=dq_frz )

end if  ! ( nc_frz > 0 )


! Now do heterogeneous nucleation...

! Wherever we have liquid present below the heterogeneous nucleation
! threshold, force q_ice to be nonzero.  Note that this nucleation value
! q_activate is set to be the smallest possible floating point number
! using TINY; there is no point adjusting T, q, qcl as the increment
! will be far smaller than the floating point precision.  It is only
! there to add nucleating points to the ice compression list so that
! the ice subsequently grows via other processes
! (vapour deposition, riming, etc)
do ic2 = 1, nc_liq
  ic = index_ic_liq(ic2)
  if ( ref_temp(ic) <= hetnuc_temp .and.                                       &
       ( .not. q_ice(ic) > zero ) ) then
    q_ice(ic) = q_activate
    dq_frz(ic) = q_activate
    l_frz_where_no_ice = .true.
  end if
end do


! Store diagnostics, if requested...
if ( l_diags ) then
  if ( nc_frz > 0 .or. l_frz_where_no_ice ) then
    ! If any freezing was done...

    ! Diagnostic of freezing/melting increment to q_liq
    if ( moist_proc_diags % diags_cond(i_liq)%pt % dq_frzmlt % flag ) then
      ! Extract super-array address
      i_super = moist_proc_diags % diags_cond(i_liq)%pt % dq_frzmlt % i_super
      ! Increment the diagnostic
      do ic2 = 1, nc_liq
        ic = index_ic_liq(ic2)
        diags_super(ic,i_super) = diags_super(ic,i_super) - dq_frz(ic)
      end do
    end if

    ! Diagnostic of freezing/melting increment to q_ice
    if ( moist_proc_diags % diags_cond(i_ice)%pt % dq_frzmlt % flag ) then
      ! Extract super-array address
      i_super = moist_proc_diags % diags_cond(i_ice)%pt % dq_frzmlt % i_super
      ! Increment the diagnostic
      do ic2 = 1, nc_liq
        ic = index_ic_liq(ic2)
        diags_super(ic,i_super) = diags_super(ic,i_super) + dq_frz(ic)
      end do
    end if

  end if  ! ( nc_frz > 0 )
end if  ! ( l_diags )

! If full freezing was done, need to regenerate the list of
! points containing liquid to exclude fully-frozen points
if ( l_full_frz ) then
  nc_tmp = nc_liq
  nc_liq = 0
  do ic2 = 1, nc_tmp
    ic = index_ic_liq(ic2)
    if ( q_liq(ic) > zero ) then
      nc_liq = nc_liq + 1
      index_ic_liq(nc_liq) = ic
    end if
  end do
end if

! If new ice added where none yet existed, need to regenerate
! the list of points containing ice to include new points
if ( l_frz_where_no_ice ) then
  nc_ice = 0
  do ic = 1, n_points
    if ( q_ice(ic) > zero ) then
      nc_ice = nc_ice + 1
      index_ic_ice(nc_ice) = ic
    end if
  end do
end if


return
end subroutine ice_nucleation

end module ice_nucleation_mod
