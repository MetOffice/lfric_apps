! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_sat_height_mod

implicit none

contains

! Subroutine to calculate an interpolated height at which
! the parcel first reached saturation, based on the parcel
! properties before and after a level-step
subroutine calc_sat_height(                                                    &
             n_points, n_points_prev, n_points_next,                           &
             n_points_buoy, n_sat, index_ic_sat,                               &
             max_buoy_heights, n_buoy_vars, i_buoy,                            &
             prev_pressure, par_prev_fields, par_next_fields,                  &
             env_prev_virt_temp, env_next_virt_temp,                           &
             linear_qs_super, i_next, i_sat, buoyancy_super,                   &
             i_core_sat )

use comorph_constants_mod, only: real_cvprec, zero, one,                       &
                     n_cond_species, i_cond_cl
use fields_type_mod, only: i_temperature, i_q_vap,                             &
                           i_qc_first, i_qc_last
use set_cp_tot_mod, only: set_cp_tot
use lat_heat_mod, only: lat_heat_incr, i_phase_change_evp
use set_qsat_mod, only: set_qsat_liq
use calc_virt_temp_mod, only: calc_virt_temp
use linear_qs_mod, only: n_linear_qs_fields, i_ref_temp,                       &
                         i_qsat_liq_ref, i_dqsatdt_liq
use grid_type_mod, only: i_height
use buoyancy_mod, only: i_prev

implicit none

! Number of points
integer, intent(in) :: n_points
! Dimensions of par_prev super-array
integer, intent(in) :: n_points_prev
! Dimensions of par_next super-array
integer, intent(in) :: n_points_next
! Dimensions of buoyancy super-array
integer, intent(in) :: n_points_buoy

! Number of points which have just hit saturation
integer, intent(in) :: n_sat
! Indices of those points
integer, intent(in) :: index_ic_sat(n_sat)

! Dimensions of the buoyancy super-array (may differ between
! different calls to this routine)
integer, intent(in) :: max_buoy_heights
integer, intent(in) :: n_buoy_vars
! Address of buoyancy fields in the super-array
integer, intent(in) :: i_buoy

! Pressure at previous level
real(kind=real_cvprec), intent(in) :: prev_pressure(n_points)

! Super-arrays containing parcel properties at the previous and
! next levels
real(kind=real_cvprec), intent(in) :: par_prev_fields                          &
                       ( n_points_prev, i_temperature:i_qc_last )
real(kind=real_cvprec), intent(in) :: par_next_fields                          &
                       ( n_points_next, i_temperature:i_qc_last )

! Environment virtual temperatures at previous and next levels
real(kind=real_cvprec), intent(in) :: env_prev_virt_temp                       &
                                      ( n_points )
real(kind=real_cvprec), intent(in) :: env_next_virt_temp                       &
                                      ( n_points )

! Super-array containing qsat at a reference temperature,
! and dqsat/dT, for linearised qsat calculations.
! Valid at the "next" level
real(kind=real_cvprec), intent(in) :: linear_qs_super                          &
                            ( n_points, n_linear_qs_fields )

! Address of next model-level interface in buoyancy_super
integer, intent(in out) :: i_next(n_points)
! Address of the newly found saturation height within the
! buoyancy super-array
integer, intent(in out) :: i_sat(n_points)

! Super-array storing the buoyancies of the parcel mean and core
! at up to 4 sub-level heights within the current level-step:
! a) Start of the level-step
! b) End of the level-step
! c) Height where parcel mean properties first hit saturation
! d) Height where parcel core first hit saturation
real(kind=real_cvprec), intent(in out) :: buoyancy_super                       &
                ( n_points_buoy, n_buoy_vars, max_buoy_heights )

! Address of the previously calculated parcel core saturation
! height within the buoyancy super-array
! (only passed in if this is a parcel mean ascent whose
!  saturation height can be modified by a separate core ascent).
integer, optional, intent(in out) :: i_core_sat(n_points)

! Parcel heat capacity
real(kind=real_cvprec) :: cp_tot(n_sat)

! Parcel liquid-water-temperature, q_vap+q_cl, and pressure,
! compressed onto saturated points
real(kind=real_cvprec) :: par_templ(n_sat)
real(kind=real_cvprec) :: par_q_vl(n_sat)
real(kind=real_cvprec) :: par_pressure(n_sat)
! Compressed condensed water species for virtual temperature calc
real(kind=real_cvprec) :: par_q_cond(n_sat,n_cond_species)


! Parcel virtual temperature and supersaturation we would've had
! without any condensation
real(kind=real_cvprec) :: prev_virt_temp_nocond(n_sat)
real(kind=real_cvprec) :: next_virt_temp_nocond(n_sat)
real(kind=real_cvprec) :: prev_ss(n_sat)
real(kind=real_cvprec) :: next_ss(n_sat)

! Height and parcel virtual temperature at saturation point
real(kind=real_cvprec) :: sat_height(n_sat)
real(kind=real_cvprec) :: sat_virt_temp(n_sat)

! Interpolation weight
real(kind=real_cvprec) :: interp

! Points where saturation height distinct from existing
! core saturation height
integer :: n_new
integer :: index_ic_new(n_sat)


! Loop counters
integer :: ic, ic2, i_cond, i_field, i_lev


!----------------------------------------------------------------
! 1) Compute previous parcel virtual temperature and
!    supersaturation if condensate were evaporated
!----------------------------------------------------------------

! Compress previous parcel properties onto saturated points
do ic2 = 1, n_sat
  ic = index_ic_sat(ic2)
  par_templ(ic2) = par_prev_fields(ic,i_temperature)
  par_q_vl(ic2)  = par_prev_fields(ic,i_q_vap)
  par_pressure(ic2) = prev_pressure(ic)
end do
do i_cond = 1, n_cond_species
  do ic2 = 1, n_sat
    par_q_cond(ic2,i_cond) = par_prev_fields(index_ic_sat(ic2),                &
                                             i_qc_first-1+i_cond)
  end do
end do

! Calculate total heat capacity of previous parcel
call set_cp_tot( n_sat, n_sat, par_q_vl, par_q_cond, cp_tot )

! Calculate liquid-water temperature of the parcel
call lat_heat_incr( n_sat, n_sat, i_phase_change_evp,                          &
                    cp_tot, par_templ,                                         &
                    dq=par_q_cond(:,i_cond_cl) )

! Add the liquid cloud to the vapour
do ic2 = 1, n_sat
  par_q_vl(ic2) = par_q_vl(ic2) + par_q_cond(ic2,i_cond_cl)
  par_q_cond(ic2,i_cond_cl) = zero
end do

! Calculate virtual temperature
call calc_virt_temp( n_sat, n_sat,                                             &
                     par_templ, par_q_vl, par_q_cond,                          &
                     prev_virt_temp_nocond )

! Calculate supersaturation
call set_qsat_liq( n_sat, par_templ, par_pressure, prev_ss )
do ic2 = 1, n_sat
  prev_ss(ic2) = par_q_vl(ic2) - prev_ss(ic2)
end do


!----------------------------------------------------------------
! 2) Compute next parcel virtual temperature and
!    supersaturation if condensate were evaporated
!----------------------------------------------------------------

! Compress next parcel properties onto saturated points
do ic2 = 1, n_sat
  ic = index_ic_sat(ic2)
  par_templ(ic2) = par_next_fields(ic,i_temperature)
  par_q_vl(ic2)  = par_next_fields(ic,i_q_vap)
end do
do i_cond = 1, n_cond_species
  do ic2 = 1, n_sat
    par_q_cond(ic2,i_cond) = par_next_fields(index_ic_sat(ic2),                &
                                             i_qc_first-1+i_cond)
  end do
end do

! Calculate total heat capacity of next parcel
call set_cp_tot( n_sat, n_sat, par_q_vl, par_q_cond, cp_tot )

! Calculate liquid-water temperature of the parcel
call lat_heat_incr( n_sat, n_sat, i_phase_change_evp,                          &
                    cp_tot, par_templ,                                         &
                    dq=par_q_cond(:,i_cond_cl) )

! Add the liquid cloud to the vapour
do ic2 = 1, n_sat
  par_q_vl(ic2) = par_q_vl(ic2) + par_q_cond(ic2,i_cond_cl)
  par_q_cond(ic2,i_cond_cl) = zero
end do

! Calculate virtual temperature
call calc_virt_temp( n_sat, n_sat,                                             &
                     par_templ, par_q_vl, par_q_cond,                          &
                     next_virt_temp_nocond )

! Calculate supersaturation, using linearised qsat passed in
do ic2 = 1, n_sat
  ic = index_ic_sat(ic2)
  next_ss(ic2) = par_q_vl(ic2)                                                 &
     - ( linear_qs_super(ic,i_qsat_liq_ref)                                    &
       + linear_qs_super(ic,i_dqsatdt_liq)                                     &
         * ( par_templ(ic2) - linear_qs_super(ic,i_ref_temp) )  )
end do


!----------------------------------------------------------------
! 3) Interpolate to find height where parcel RH reached 100%
!----------------------------------------------------------------

do ic2 = 1, n_sat
  ic = index_ic_sat(ic2)

  ! Fraction of the way through the level where SS is zero
  interp = -prev_ss(ic2) / ( next_ss(ic2) - prev_ss(ic2) )

  ! Interpolate height
  sat_height(ic2)                                                              &
         = (one-interp) * buoyancy_super(ic,i_height,i_prev)                   &
         +      interp  * buoyancy_super(ic,i_height,i_next(ic))

  ! Interpolate to estimate parcel virtual temperature at the
  ! saturation height
  sat_virt_temp(ic2)                                                           &
         = (one-interp) * prev_virt_temp_nocond(ic2)                           &
         +      interp  * next_virt_temp_nocond(ic2)

end do

do ic2 = 1, n_sat
  ! Overwrite with just prev or next value in case where
  ! SS doesn't change sign between prev and next
  ! (should only happen occasionally due to rounding errors)
  if ( prev_ss(ic2) * next_ss(ic2) >= zero ) then
    ic = index_ic_sat(ic2)
    if ( abs(next_ss(ic2)) < abs(prev_ss(ic2)) ) then
      sat_height(ic2)    = buoyancy_super(ic,i_height,i_next(ic))
      sat_virt_temp(ic2) = next_virt_temp_nocond(ic2)
    else
      sat_height(ic2)    = buoyancy_super(ic,i_height,i_prev)
      sat_virt_temp(ic2) = prev_virt_temp_nocond(ic2)
    end if
  end if
end do

! Even when SS does change sign, rounding errors can still very
! occasionally cause sat_height to fall a tiny bit outside its
! allowed range; force it to be between prev height and next height!
do ic2 = 1, n_sat
  ic = index_ic_sat(ic2)
  sat_height(ic2) = max( min( sat_height(ic2),                                 &
                              buoyancy_super(ic,i_height,i_next(ic)) ),        &
                         buoyancy_super(ic,i_height,i_prev) )
end do

! Modify calculation of saturation height in case where the
! parcel core hit saturation before the mean, or remained
! saturated for longer than the mean.
! While the core is saturated, the mean must also contain
! saturated air and so is considered to have the same saturation
! height as the core.
if ( present(i_core_sat) ) then

  do ic2 = 1, n_sat
    ic = index_ic_sat(ic2)
    ! If the core hit saturation at this level
    if ( i_core_sat(ic) > 0 ) then

      ! If going from subsaturated to saturated
      if ( next_ss(ic2) > prev_ss(ic2) ) then
        ! If the core hit saturation before the mean
        if ( buoyancy_super(ic,i_height,i_core_sat(ic))                        &
             <= sat_height(ic2) ) then
          ! Reset saturation height for the mean same as core
          i_sat(ic) = i_core_sat(ic)
        end if
        ! If going from saturated to subsaturated
      else
        ! If the core became subsaturated after the mean
        if ( buoyancy_super(ic,i_height,i_core_sat(ic))                        &
             >= sat_height(ic2) ) then
          ! Reset saturation height for the mean same as core
          i_sat(ic) = i_core_sat(ic)
        end if
      end if

      if ( i_sat(ic) == i_core_sat(ic) ) then
        ! If sat height has been reset to that of the core,
        ! re-interpolate virt temp at sat_height
        interp = ( buoyancy_super(ic,i_height,i_core_sat(ic))                  &
                 - buoyancy_super(ic,i_height,i_prev) )                        &
               / ( buoyancy_super(ic,i_height,i_next(ic))                      &
                 - buoyancy_super(ic,i_height,i_prev) )
        sat_virt_temp(ic2)                                                     &
               = (one-interp) * prev_virt_temp_nocond(ic2)                     &
               +      interp  * next_virt_temp_nocond(ic2)
        sat_height(ic2)                                                        &
               = buoyancy_super(ic,i_height,i_core_sat(ic))
      end if

    end if  ! ( i_core_sat(ic) > 0 )
  end do  ! ic2 = 1, n_sat

  ! Find points where new saturation height found and not reset
  ! to the core saturation height
  n_new = 0
  do ic2 = 1, n_sat
    ic = index_ic_sat(ic2)
    if ( i_sat(ic) == 0 ) then
      n_new = n_new + 1
      index_ic_new(n_new) = ic
      ! Compress sat_height onto "new" points
      sat_height(n_new) = sat_height(ic2)
    end if
  end do

else  ! ( present(i_core_sat) )

  ! No separate core ascent; all sat points included in the list
  n_new = n_sat
  do ic2 = 1, n_sat
    index_ic_new(ic2) = index_ic_sat(ic2)
  end do

end if    ! ( present(i_core_sat) )


! Find level address of saturation height in buoyancy_super
do ic2 = 1, n_new
  ic = index_ic_new(ic2)
  ! Find i_sat
  i_sat(ic) = i_prev + 1
  if ( i_next(ic) > i_prev + 1 ) then
    do i_lev = i_prev+1, i_next(ic)-1
      if ( sat_height(ic2) >=                                                  &
           buoyancy_super(ic,i_height,i_lev) )                                 &
         i_sat(ic) = i_lev + 1
    end do
  end if
  ! Shift all the data beyond this forward by 1 sub-level
  do i_lev = i_next(ic), i_sat(ic), -1
    do i_field = 1, n_buoy_vars
      buoyancy_super(ic,i_field,i_lev+1)                                       &
        = buoyancy_super(ic,i_field,i_lev)
    end do
  end do
  ! Store height of the newly inserted saturation level
  buoyancy_super(ic,i_height,i_sat(ic)) = sat_height(ic2)
end do

! Shift address of next model-level interface
do ic2 = 1, n_new
  ic = index_ic_new(ic2)
  i_next(ic) = i_next(ic) + 1
end do
! Also shift address of core saturation height, if used
if ( present(i_core_sat) ) then
  do ic2 = 1, n_new
    ic = index_ic_new(ic2)
    if ( i_core_sat(ic) >= i_sat(ic) )                                         &
                  i_core_sat(ic) = i_core_sat(ic) + 1
  end do
end if

! Interpolate to get environment virt_temp at sat_height,
! to find the virt_temp excess at sat_height.
do ic2 = 1, n_sat
  ic = index_ic_sat(ic2)
  ! Interpolate sat between prev and next
  interp = ( buoyancy_super(ic,i_height,i_sat(ic))                             &
           - buoyancy_super(ic,i_height,i_prev) )                              &
         / ( buoyancy_super(ic,i_height,i_next(ic))                            &
           - buoyancy_super(ic,i_height,i_prev) )
  buoyancy_super(ic,i_buoy,i_sat(ic)) = sat_virt_temp(ic2)                     &
    - ( (one-interp) * env_prev_virt_temp(ic)                                  &
      +      interp  * env_next_virt_temp(ic) )
end do


return
end subroutine calc_sat_height


end module calc_sat_height_mod
