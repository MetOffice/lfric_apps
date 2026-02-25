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
             n_points, n_points_sublevs, l_mean_with_core, l_down, j_buoy,     &
             prev_ss, next_ss, prev_tvl, next_tvl,                             &
             par_prev_q_cl, par_next_q_cl,                                     &
             i_next, i_sat, sublevs,                                           &
             i_core_sat )

use comorph_constants_mod, only: real_cvprec, zero, one, sqrt_min_float
use sublevs_mod, only: max_sublevs, n_sublev_vars, i_prev,                     &
                       j_height, j_env_tv, j_mean_buoy, j_core_buoy, j_delta_tv

implicit none

! Number of points
integer, intent(in) :: n_points
! Dimensions of sub-levels super-array
integer, intent(in) :: n_points_sublevs

! Flag for parcel mean ascent with an accompanying core
logical, intent(in) :: l_mean_with_core

! Flag for going down (downdrafts)
logical, intent(in) :: l_down

! Address of buoyancy fields in the sub-levels super-array
integer, intent(in) :: j_buoy

! Parcel supersaturation and virtual temperature we would've had
! without any condensation, at previous and next level
real(kind=real_cvprec), intent(in) :: prev_ss(n_points)
real(kind=real_cvprec), intent(in) :: next_ss(n_points)
real(kind=real_cvprec), intent(in) :: prev_tvl(n_points)
real(kind=real_cvprec), intent(in) :: next_tvl(n_points)

! Parcel liquid-cloud content at prev and next
real(kind=real_cvprec), intent(in) :: par_prev_q_cl(n_points)
real(kind=real_cvprec), intent(in) :: par_next_q_cl(n_points)

! Address of next model-level interface in sublevs
integer, intent(in out) :: i_next(n_points)
! Address of the newly found saturation height within the
! sub-levels super-array
integer, intent(in out) :: i_sat(n_points)

! Super-array storing parcel buoyancies and other properties
! at up to 4 sub-level heights within the current level-step:
! a) Start of the level-step
! b) End of the level-step
! c) Height where parcel mean properties first hit saturation
! d) Height where parcel core first hit saturation
real(kind=real_cvprec), intent(in out) :: sublevs                              &
                ( n_points_sublevs, n_sublev_vars, max_sublevs )

! Address of the previously calculated parcel core saturation
! height within the sublevs super-array
! (only passed in if this is a parcel mean ascent whose
!  saturation height can be modified by a separate core ascent).
integer, optional, intent(in out) :: i_core_sat(n_points)

! Number of points which have just hit saturation
integer :: n_sat
! Indices of those points
integer :: index_ic_sat(n_points)

! Points where saturation height distinct from existing
! core saturation height
integer :: n_new
integer :: index_ic_new(n_points)

! Height and parcel virtual temperatures at saturation height
real(kind=real_cvprec) :: sat_height(n_points)
real(kind=real_cvprec) :: sat_par_virt_temp(n_points)

! Interpolation weight
real(kind=real_cvprec) :: interp

! Loop counters
integer :: ic, ic2, i_field, i_lev


!------------------------------------------------------------------------------
! 1) Find saturation height properties for the current ascent
!------------------------------------------------------------------------------

! Find list of points which just crossed from subsaturated to
! saturated, or vice-versa
n_sat = 0
do ic = 1, n_points
  if ( par_prev_q_cl(ic) > zero .neqv. par_next_q_cl(ic) > zero ) then
    n_sat = n_sat + 1
    index_ic_sat(n_sat) = ic
  end if
end do

if ( n_sat > 0 ) then
  ! If any points just crossed saturation...

  ! Interpolate to find height where parcel RH reached 100%
  do ic2 = 1, n_sat
    ic = index_ic_sat(ic2)

    ! Fraction of the way through the level where SS is zero
    interp = -prev_ss(ic) / ( next_ss(ic) - prev_ss(ic) )

    ! Interpolate height
    sat_height(ic) = (one-interp) * sublevs(ic,j_height,i_prev)                &
                   +      interp  * sublevs(ic,j_height,i_next(ic))

    ! Interpolate to estimate parcel virtual temperature at the
    ! saturation height
    sat_par_virt_temp(ic) = (one-interp) * prev_tvl(ic)                        &
                          +      interp  * next_tvl(ic)
  end do

  do ic2 = 1, n_sat
    ic = index_ic_sat(ic2)
    ! Overwrite with just prev or next value in case where
    ! SS doesn't change sign between prev and next
    ! (should only happen occasionally due to rounding errors)
    if ( prev_ss(ic) * next_ss(ic) >= zero ) then
      if ( abs(next_ss(ic)) < abs(prev_ss(ic)) ) then
        sat_height(ic)        = sublevs(ic,j_height,i_next(ic))
        sat_par_virt_temp(ic) = next_tvl(ic)
      else
        sat_height(ic)        = sublevs(ic,j_height,i_prev)
        sat_par_virt_temp(ic) = prev_tvl(ic)
      end if
    end if
  end do

  ! Rounding errors can very occasionally cause sat_height to fall a tiny bit
  ! outside its allowed range; force it to be between
  ! prev height and next height!
  if ( l_down ) then
    ! Downdrafts: want    next_height < sat_height < prev_height
    do ic2 = 1, n_sat
      ic = index_ic_sat(ic2)
      sat_height(ic) = min( max( sat_height(ic),                               &
                                 sublevs(ic,j_height,i_next(ic)) ),            &
                            sublevs(ic,j_height,i_prev) )
    end do
  else
    ! Updrafts:   want    prev_height < sat_height < next_height
    do ic2 = 1, n_sat
      ic = index_ic_sat(ic2)
      sat_height(ic) = max( min( sat_height(ic),                               &
                                 sublevs(ic,j_height,i_next(ic)) ),            &
                            sublevs(ic,j_height,i_prev) )
    end do
  end if

  ! Modify calculation of saturation height in case where the
  ! parcel core hit saturation before the mean, or remained
  ! saturated for longer than the mean.
  ! While the core is saturated, the mean must also contain
  ! saturated air and so is considered to have the same saturation
  ! height as the core.
  if ( l_mean_with_core .and. present(i_core_sat) ) then
    do ic2 = 1, n_sat
      ic = index_ic_sat(ic2)
      if ( i_core_sat(ic) > 0 ) then
        ! If the core hit saturation at this level

        if ( l_down ) then
          ! Going down
          if ( next_ss(ic) > prev_ss(ic) ) then
            if ( sublevs(ic,j_height,i_core_sat(ic)) >= sat_height(ic) ) then
              ! Core reached saturation first
              i_sat(ic) = i_core_sat(ic)
            end if
          else
            if ( sublevs(ic,j_height,i_core_sat(ic)) <= sat_height(ic) ) then
              ! Core stayed saturated for longer
              i_sat(ic) = i_core_sat(ic)
            end if
          end if
        else  ! ( .NOT. l_down )
          ! Going up
          if ( next_ss(ic) > prev_ss(ic) ) then
            if ( sublevs(ic,j_height,i_core_sat(ic)) <= sat_height(ic) ) then
              ! Core reached saturation first
              i_sat(ic) = i_core_sat(ic)
            end if
          else
            if ( sublevs(ic,j_height,i_core_sat(ic)) >= sat_height(ic) ) then
              ! Core stayed saturated for longer
              i_sat(ic) = i_core_sat(ic)
            end if
          end if
        end if  ! ( .NOT. l_down )

        if ( i_sat(ic) == i_core_sat(ic) ) then
          ! If sat height has been reset to that of the core,
          ! re-interpolate virt temp at sat_height
          interp = ( sublevs(ic,j_height,i_core_sat(ic))                       &
                   - sublevs(ic,j_height,i_prev) )                             &
                 / ( sublevs(ic,j_height,i_next(ic))                           &
                   - sublevs(ic,j_height,i_prev) )
          sat_par_virt_temp(ic) = (one-interp) * prev_tvl(ic)                  &
                                +      interp  * next_tvl(ic)
          sat_height(ic) = sublevs(ic,j_height,i_core_sat(ic))
        end if

      end if  ! ( i_core_sat(ic) > 0 )
    end do  ! ic2 = 1, n_sat

    ! Find points where new saturation height found and NOT reset
    ! to the core saturation height
    n_new = 0
    do ic2 = 1, n_sat
      ic = index_ic_sat(ic2)
      if ( i_sat(ic) == 0 ) then
        n_new = n_new + 1
        index_ic_new(n_new) = ic
      end if
    end do

  else  ! ( .NOT. ( l_mean_with_core .AND. PRESENT(i_core_sat) ) )

    ! No separate core ascent; all sat points included in the list
    n_new = n_sat
    do ic2 = 1, n_sat
      index_ic_new(ic2) = index_ic_sat(ic2)
    end do

  end if  ! ( .NOT. ( l_mean_with_core .AND. PRESENT(i_core_sat) ) )


  !----------------------------------------------------------------------------
  ! 2) Store the required data in the sublevs super-array
  !----------------------------------------------------------------------------

  ! Find level address of saturation height in the sublevs array
  ! and insert the data at saturation height
  ! Height comparisons have opposite sign for updrafts and downdrafts
  if ( l_down ) then
    do ic2 = 1, n_new
      ic = index_ic_new(ic2)
      i_sat(ic) = i_prev + 1
      if ( i_next(ic) > i_prev + 1 ) then
        do i_lev = i_prev+1, i_next(ic)-1
          if ( sat_height(ic) <= sublevs(ic,j_height,i_lev) ) then
            i_sat(ic) = i_lev + 1
          end if
        end do
      end if
    end do
  else
    do ic2 = 1, n_new
      ic = index_ic_new(ic2)
      i_sat(ic) = i_prev + 1
      if ( i_next(ic) > i_prev + 1 ) then
        do i_lev = i_prev+1, i_next(ic)-1
          if ( sat_height(ic) >= sublevs(ic,j_height,i_lev) ) then
            i_sat(ic) = i_lev + 1
          end if
        end do
      end if
    end do
  end if

  do ic2 = 1, n_new
    ic = index_ic_new(ic2)

    ! Shift all the data beyond i_sat forward by 1 sub-level
    do i_lev = i_next(ic), i_sat(ic), -1
      do i_field = 1, n_sublev_vars
        sublevs(ic,i_field,i_lev+1) = sublevs(ic,i_field,i_lev)
      end do
    end do

    ! Shift address of next model-level interface
    i_next(ic) = i_next(ic) + 1

    ! Linearly interpolate all sublev fields to sat_height
    ! (except for j_height and the buoyancies and vertical velocity excesses
    !  which are set subsequently...)
    interp = ( sat_height(ic)                  - sublevs(ic,j_height,i_prev) ) &
           / ( sublevs(ic,j_height,i_next(ic)) - sublevs(ic,j_height,i_prev) )
    do i_field = j_height+1, j_mean_buoy-1
      sublevs(ic,i_field,i_sat(ic))                                            &
        = (one-interp) * sublevs(ic,i_field,i_prev)                            &
        +      interp  * sublevs(ic,i_field,i_next(ic))
    end do

    ! TEMPORARY CODE TO PRESERVE KGO; TO BE REMOVED SOON:
    ! delta_tv is assumed constant over whole model-level;
    ! remove bit-level changes due to interpolation
    sublevs(ic,j_delta_tv,i_sat(ic)) = sublevs(ic,j_delta_tv,i_prev)

    ! Store the found saturation height
    sublevs(ic,j_height,i_sat(ic)) = sat_height(ic)

  end do  ! ic2 = 1, n_new

  ! Also shift address of core saturation height, if used
  if ( l_mean_with_core .and. present(i_core_sat) ) then
    do ic2 = 1, n_new
      ic = index_ic_new(ic2)
      if ( i_core_sat(ic) >= i_sat(ic) )  i_core_sat(ic) = i_core_sat(ic) + 1
    end do
  end if

  ! Set buoyancy at saturation height using interpolated env Tv
  do ic2 = 1, n_sat
    ic = index_ic_sat(ic2)
    !sublevs(ic,j_buoy,i_sat(ic)) = sat_par_virt_temp(ic)                      &
    !                             - sublevs(ic,j_env_tv,i_sat(ic))
    ! TEMPORARY CODE TO PRESERVE KGO
    ! (re-do interpolation of env Tv to sat height; using the already
    !  calculated value changes answers for CCE high optimisation)
    interp = ( sublevs(ic,j_height,i_sat(ic))                                  &
             - sublevs(ic,j_height,i_prev) )                                   &
           / ( sublevs(ic,j_height,i_next(ic))                                 &
             - sublevs(ic,j_height,i_prev) )
    sublevs(ic,j_buoy,i_sat(ic)) = sat_par_virt_temp(ic)                       &
      - ( (one-interp) * sublevs(ic,j_env_tv,i_prev)                           &
        +      interp  * sublevs(ic,j_env_tv,i_next(ic)) )
  end do

end if  ! ( n_sat > 0 )


! If this is a mean ascent following a core ascent...
if ( l_mean_with_core .and. present(i_core_sat) ) then
  ! We need to calculate the mean's virt_temp excess at the
  ! core's saturation height, and vice-versa...

  ! Find points where the core hit saturation (and the mean's
  ! saturation height is not defined at the same height)
  n_sat = 0
  do ic = 1, n_points
    if ( i_core_sat(ic) > 0 .and. ( .not. i_sat(ic) == i_core_sat(ic) ) ) then
      n_sat = n_sat + 1
      index_ic_sat(n_sat) = ic
    end if
  end do

  ! If any points
  if ( n_sat > 0 ) then
    ! Interpolate mean buoyancy at core sat_height.
    do ic2 = 1, n_sat
      ic = index_ic_sat(ic2)
      interp = ( sublevs(ic,j_height,i_core_sat(ic))                           &
               - sublevs(ic,j_height,i_core_sat(ic)-1) )                       &
             / ( sublevs(ic,j_height,i_core_sat(ic)+1)                         &
               - sublevs(ic,j_height,i_core_sat(ic)-1) )
      sublevs(ic,j_mean_buoy,i_core_sat(ic))                                   &
        = (one-interp) * sublevs(ic,j_mean_buoy,i_core_sat(ic)-1)              &
        +      interp  * sublevs(ic,j_mean_buoy,i_core_sat(ic)+1)
    end do
  end if

  ! Find points where the mean hit saturation (and the core's
  ! saturation height is not defined at the same height)
  n_sat = 0
  do ic = 1, n_points
    if ( i_sat(ic) > 0 .and. ( .not. i_sat(ic) == i_core_sat(ic) ) ) then
      n_sat = n_sat + 1
      index_ic_sat(n_sat) = ic
    end if
  end do

  ! If any points
  if ( n_sat > 0 ) then
    ! Interpolate core buoyancy at mean sat_height.
    do ic2 = 1, n_sat
      ic = index_ic_sat(ic2)
      interp = ( sublevs(ic,j_height,i_sat(ic))                                &
               - sublevs(ic,j_height,i_sat(ic)-1) )                            &
             / ( sublevs(ic,j_height,i_sat(ic)+1)                              &
               - sublevs(ic,j_height,i_sat(ic)-1) )
      sublevs(ic,j_core_buoy,i_sat(ic))                                        &
        = (one-interp) * sublevs(ic,j_core_buoy,i_sat(ic)-1)                   &
        +      interp  * sublevs(ic,j_core_buoy,i_sat(ic)+1)
    end do
  end if

end if  ! ( l_mean_with_core .AND. PRESENT(i_core_sat) )


return
end subroutine calc_sat_height


end module calc_sat_height_mod
