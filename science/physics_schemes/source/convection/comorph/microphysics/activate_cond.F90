! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module activate_cond_mod

implicit none

contains

! Subroutine to activate ncondensation of water vapour to form
! new liquid cloud on CCN
subroutine activate_cond( n_points, nc, index_ic,                              &
                          qsat_liq_ref, dqsatdT_liq, ref_temp,                 &
                          temperature, q_vap, q_cond )

use comorph_constants_mod, only: real_cvprec, zero, q_activate

implicit none

! Number of points
integer, intent(in) :: n_points
! Number of points where condensed water field is non-zero
integer, intent(in out) :: nc
! Indices of non-zero points
integer, intent(in out) :: index_ic(n_points)
! (these are inout as activation may add new points containing
!  liquid water where there was none before)

! qsat w.r.t. liquid at reference temperature and dqsat/dT
real(kind=real_cvprec), intent(in) :: qsat_liq_ref(n_points)
real(kind=real_cvprec), intent(in) :: dqsatdT_liq(n_points)
real(kind=real_cvprec), intent(in) :: ref_temp(n_points)

! Air temperature
real(kind=real_cvprec), intent(in) :: temperature(n_points)

! Water vapour mixing ratio
real(kind=real_cvprec), intent(in) :: q_vap(n_points)

! Condensed water mixing ratio to be initialised to tiny
! non-zero value at activated points; the actual condensation
! rate will be calculated later.
real(kind=real_cvprec), intent(in out) :: q_cond(n_points)

! Flag for activation done at any points
logical :: l_any_activate


! Loop counter
integer :: ic

! The rest of the code after this routine will only do anything
! at points where there is already non-zero condensate.
! To allow condensation to initiate, a minute amount of
! liquid cloud is seeded here.
! Note we don't fully condense to adjust the air to saturation here;
! this is because other processes maybe acting (e.g. freezing or melting)
! which could change the RH during the current timestep, so at this
! point we don't know how much condensation we need to do to adjust to
! saturation.  All the microphysical processes that affect T and q_vap
! (including condensation / evaporation of liquid-cloud) are solved
! simultaneously in phase_change_solve, so the adjustment to saturation
! is done in there.
! If the air actually ends up sub-saturated, the minute
! seeded condensate added here will just be re-evaporated without
! consequence.
! If we make this into a 2-moment scheme in the future,
! the activated number concentration should also be set here,
! based on the rate of increase of saturation and the
! modelled aerosols spectrum.

l_any_activate = .false.
do ic = 1, n_points
  ! If zero condensed water mixing ratio present
  if ( q_cond(ic) == zero ) then
    ! If at or above liquid saturation point
    if ( q_vap(ic) >= qsat_liq_ref(ic)                                         &
                    + dqsatdT_liq(ic)                                          &
                      * (temperature(ic) - ref_temp(ic))  ) then
      q_cond(ic) = q_activate
      ! Note: technically we have added water from nowhere
      ! here; however, q_activate is set to the smallest
      ! possible non-zero floating-point number, using tiny.
      ! This is almost certainly smaller than q_vap * epsilon,
      ! so that subtracting it from the water vapour would
      ! make no difference to the floating-point
      ! representation of q_vap.  Therefore, it will make
      ! no difference to conservation in the model.
      l_any_activate = .true.
    end if
  end if
end do

! If any activation was done, recreate the list of points
! containing this condensate species to include the newly
! activated points
if ( l_any_activate ) then
  nc = 0
  do ic = 1, n_points
    if ( q_cond(ic) > zero ) then
      nc = nc + 1
      index_ic(nc) = ic
    end if
  end do
end if


return
end subroutine activate_cond

end module activate_cond_mod
