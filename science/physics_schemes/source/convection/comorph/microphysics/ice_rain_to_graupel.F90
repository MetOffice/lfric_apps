! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module ice_rain_to_graupel_mod

implicit none

contains

! Subroutine to convert collided ice and rain into graupel.
! This routine assumes that colection of rain by ice-cloud has
! already been computed, with the collided rain's mixing ratio
! added onto the ice.  It transfers the collided rain's mixing
! ratio into graupel instead, and additionally calculates the
! amount of ice-cloud involved in the collisions; this is also
! converted to graupel.
subroutine ice_rain_to_graupel( n_points,                                      &
                                nc_col_cf, index_ic_col_cf,                    &
                                nc_graup, index_ic_graup,                      &
                                dq_col_cf,                                     &
                                dq_col_rain_tot, q_rain,                       &
                                kq_rain, kt_rain,                              &
                                n_cf, n_rain,                                  &
                                q_loc_cf, q_loc_rain,                          &
                                q_cf, q_graup,                                 &
                                dq_frz_cf, dq_frz_graup,                       &
                                kq_graup, kt_graup )


use comorph_constants_mod, only: real_cvprec, zero

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of points where ice-cloud and rain have collided
integer, intent(in) :: nc_col_cf
! Indices of those points
integer, intent(in) :: index_ic_col_cf(nc_col_cf)

! Number of points where graupel is present
integer, intent(in out) :: nc_graup
! Indices of those points
integer, intent(in out) :: index_ic_graup(n_points)
! Note: these are updated in the case where rain-ice collisions
! produce new graupel at points where there was none before.

! Amount of rain mixing ratio collected by ice-cloud
real(kind=real_cvprec), intent(in) :: dq_col_cf(n_points)

! Total increment to rain due to collection by other species
real(kind=real_cvprec), intent(in) :: dq_col_rain_tot(n_points)
! Rain mixing ratio
real(kind=real_cvprec), intent(in) :: q_rain(n_points)

! Vapour and heat exchange coefficients for rain
! (needed for applying minimum limit to graupel coefs)
real(kind=real_cvprec), intent(in) :: kq_rain(n_points)
real(kind=real_cvprec), intent(in) :: kt_rain(n_points)

! Number concentration / local mixing ratio of ice-cloud and rain
real(kind=real_cvprec), intent(in) :: n_cf(n_points)
real(kind=real_cvprec), intent(in) :: n_rain(n_points)
real(kind=real_cvprec), intent(in) :: q_loc_cf(n_points)
real(kind=real_cvprec), intent(in) :: q_loc_rain(n_points)

! Ice-cloud and graupel mixing ratios
real(kind=real_cvprec), intent(in out) :: q_cf(n_points)
real(kind=real_cvprec), intent(in out) :: q_graup(n_points)

! Total amount of liquid frozen onto ice-cloud and graupel
real(kind=real_cvprec), intent(in out) :: dq_frz_cf(n_points)
real(kind=real_cvprec), intent(in out) :: dq_frz_graup(n_points)

! Vapour and heat exchange coefficients for graupel
! (need to be updated where graupel created where none before)
real(kind=real_cvprec), intent(in out) :: kq_graup(n_points)
real(kind=real_cvprec), intent(in out) :: kt_graup(n_points)

! Mixing ratio of ice-cloud involved in collisions with rain
real(kind=real_cvprec) :: dq_col_cf_cmpr(nc_col_cf)

! Flag for adding new graupel at points where there wasn't any
logical :: l_added_where_none

! Fraction of rain collected by ice-cloud
real(kind=real_cvprec) :: frac_rain

! Loop counters
integer :: ic, ic2


! Check whether we are about to add graupel at points
! where its mixing ratio is currently zero
l_added_where_none = .false.
over_nc_col_cf: do ic2 = 1, nc_col_cf
  ic = index_ic_col_cf(ic2)
  if ( q_graup(ic) == zero ) then
    l_added_where_none = .true.
    exit over_nc_col_cf
  end if
end do over_nc_col_cf

do ic2 = 1, nc_col_cf
  ic = index_ic_col_cf(ic2)
  ! Add rain collected by ice onto graupel instead.
  q_cf(ic)    = q_cf(ic)    - dq_col_cf(ic)
  q_graup(ic) = q_graup(ic) + dq_col_cf(ic)
  ! Transfer amount freezing on ice onto graupel
  dq_frz_cf(ic)    = dq_frz_cf(ic)    - dq_col_cf(ic)
  dq_frz_graup(ic) = dq_frz_graup(ic) + dq_col_cf(ic)
end do

! If adding graupel to any points for the first time,
! recompute the graupel compression list to include
! the new points
if ( l_added_where_none ) then
  nc_graup = 0
  do ic = 1, n_points
    if ( q_graup(ic) > 0 ) then
      nc_graup = nc_graup + 1
      index_ic_graup(nc_graup) = ic
    end if
  end do
end if

do ic2 = 1, nc_col_cf
  ic = index_ic_col_cf(ic2)
  ! We initially assume that all collisions between liquid and
  ! ice produce ice.  In fact, some of it may remain molten on
  ! sleety particles.  This is accounted for by potentially
  ! remelting some of it later (melting is treated implicitly
  ! within phase_change_solve).
  ! To get the appropriate remelting, a heat exchange
  ! coefficient for the graupel is required.  However, if we have
  ! just created graupel via rain-ice collisions at points where
  ! there was little or no graupel at the start of the step,
  ! the explicit heat exchange coef for graupel will be at or
  ! near zero, leading to remelting of all the graupel created.
  ! => Apply a minimum limit to the heat exchange coef for
  !    graupel, equal to the heat exchange coef for the collected
  !    rain * the proportion of rain that was collected by ice.
  ! Calculate fraction of rain which collided with ice;
  !  = rain collected by ice
  !  / ( rain left now + rain already removed by collision
  frac_rain = dq_col_cf(ic)                                                    &
            / ( q_rain(ic) + dq_col_rain_tot(ic) )
  ! Apply minimum limit to exchange coefs
  kt_graup(ic) = max( kt_graup(ic), kt_rain(ic) * frac_rain )
  kq_graup(ic) = max( kq_graup(ic), kq_rain(ic) * frac_rain )
end do

do ic2 = 1, nc_col_cf
  ic = index_ic_col_cf(ic2)
  ! Calculate the amount of ice involved in the
  ! collisions with rain.
  ! Formula is identical to the amount of rain involved (already
  ! computed in collision_rate and stored in dq_col_cf)
  ! except scaled by n_rain*q_loc_cf instead of n_cf*q_loc_rain:
  dq_col_cf_cmpr(ic2) = dq_col_cf(ic)                                          &
                      * ( n_rain(ic) * q_loc_cf(ic) )                          &
                      / ( n_cf(ic) * q_loc_rain(ic) )

end do

! Avoid removing more ice than is present
do ic2 = 1, nc_col_cf
  ic = index_ic_col_cf(ic2)
  dq_col_cf_cmpr(ic2) = min( dq_col_cf_cmpr(ic2), q_cf(ic) )
end do

do ic2 = 1, nc_col_cf
  ic = index_ic_col_cf(ic2)
  ! Decrement the ice and increment graupel
  q_cf(ic)    = q_cf(ic)    - dq_col_cf_cmpr(ic2)
  q_graup(ic) = q_graup(ic) + dq_col_cf_cmpr(ic2)
end do

return
end subroutine ice_rain_to_graupel

end module ice_rain_to_graupel_mod
