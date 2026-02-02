! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_init_mass_mod

implicit none

contains

! Subroutine to calculate the mass initiating from a locally
! unstable model-level / sub-grid region of the grid-box.
subroutine calc_init_mass( n_points, nc, index_ic, n_conv_types, i_region,     &
                           i_next, i_sat, sublevs,                             &
                           layer_mass_k, next_height, virt_temp_next_cmpr,     &
                           frac_r_t, init_mass_t, Nsq )

use comorph_constants_mod, only: real_cvprec, min_delta, min_float,            &
                                 zero, one, two,                               &
                                 third, sqrt_min_delta, newline,               &
                                 gravity, par_gen_mass_fac
use sublevs_mod, only: max_sublevs, n_sublev_vars, i_prev,                     &
                       j_height, j_mean_buoy
use subregion_mod, only: n_regions
use raise_error_mod, only: raise_fatal

implicit none

! Total number of points being processed in init_mass_moist_frac
integer, intent(in) :: n_points

! Number of points in the current sub-grid region
integer, intent(in) :: nc

! Indices of the above points, for compression / decompression
integer, intent(in) :: index_ic(nc)

! Number of convection types
integer, intent(in) :: n_conv_types

! Index of the current sub-grid region (dry,liq,mph,icr)
integer, intent(in) :: i_region

! Address of next level and saturation height in sublevs
integer, intent(in) :: i_next(nc)
integer, intent(in) :: i_sat(nc)

! Parcel sub-levels super-array output by parcel_dyn
real(kind=real_cvprec), intent(in) :: sublevs                                  &
                                      ( nc, n_sublev_vars, max_sublevs )

! Layer mass on all points
real(kind=real_cvprec), intent(in) :: layer_mass_k(n_points)

! Height of next model-level interface
real(kind=real_cvprec), intent(in) :: next_height(n_points)

! Grid-mean virtual temperature at next, compressed onto region points
real(kind=real_cvprec), intent(in) :: virt_temp_next_cmpr(nc)

! Triggering area fraction of each convection type in each sub-grid region
real(kind=real_cvprec), intent(in) :: frac_r_t                                 &
                                      ( n_points, n_regions, n_conv_types )

! Initiating dry-mass from current level / sub-grid region
! kg m-2 s-1 (for each convection types)
real(kind=real_cvprec), intent(out) :: init_mass_t                             &
                                       ( n_points, n_conv_types )

! Moist static stability from test ascent
real(kind=real_cvprec), intent(out) :: Nsq(nc)

! Moist static stability after crossing saturation threshold
real(kind=real_cvprec) :: Nsq_next

! Linear ramp to apply to mass-sources from unstable layers that
! are only borderline vertically-resolved
real(kind=real_cvprec) :: ramp

! Work sub-level index in the sub-levels super-array
integer :: i_lev(nc)

! Loop counters
integer :: ic, ic2, i_type

! Minimum acceptable precision when computing height interval
real(kind=real_cvprec), parameter :: min_precision                             &
                   = 100.0_real_cvprec * min_delta

character(len=*), parameter :: routinename                                     &
                               = "CALC_INIT_MASS"


! Calculate moist static stability from the ascent, from
! level k to the next sub-level step output in the buoyancy
! super-array:
!            Nsq_moist = -(g/Tv) dTv'/dz
!
! Note that care must be taken regarding loss of precision
! when calculating the vertical gradient of Tv' between
! adjacent sub-level steps.  If the 2nd sub-level height
! is so close to the first that the height difference
! delta_z will have insufficient precision, try the next
! sub-level interval instead.
do ic2 = 1, nc
  i_lev(ic2) = i_prev
  if ( abs( sublevs(ic2,j_height,i_prev+1) - sublevs(ic2,j_height,i_prev) )    &
     < min_precision * sublevs(ic2,j_height,i_prev) ) then
    i_lev(ic2) = i_lev(ic2) + 1

    ! NOTE: The above may yield horific nonsense if the
    ! vertical resolution ever gets so fine that the whole
    ! half-level height interval implies too much loss of
    ! precision in delta_z to safely compute the gradient.
    ! Raise a fatal error if this ever happens:
    if ( i_lev(ic2) >= i_next(ic2) ) then
      call raise_fatal( routinename,                                           &
             "Failed to find a sub-level step with a large "  //               &
             "enough height interval delta_z to "    //newline//               &
             "safely compute Nsq from the vertical gradient " //               &
             "of buoyancy!  You either need to "     //newline//               &
             "increase the precision of the convection "      //               &
             "scheme variables, or use a coarser "   //newline//               &
             "vertical grid." )
    end if
  end if
end do

! Do the actual calculation of Nsq
do ic2 = 1, nc
  Nsq(ic2) = -( gravity / virt_temp_next_cmpr(ic2) )                           &
           * ( sublevs(ic2,j_mean_buoy,i_lev(ic2)+1)                           &
             - sublevs(ic2,j_mean_buoy,i_lev(ic2)) )                           &
           / ( sublevs(ic2,j_height,i_lev(ic2)+1)                              &
             - sublevs(ic2,j_height,i_lev(ic2)) )
end do

do i_type = 1, n_conv_types
  ! For each convection type...

  ! Calculate initiating mass as a function of Nsq
  do ic2 = 1, nc
    ic = index_ic(ic2)
    init_mass_t(ic2,i_type)                                                    &
      = par_gen_mass_fac * layer_mass_k(ic) * frac_r_t(ic,i_region,i_type)     &
        * sqrt( max( -Nsq(ic2), zero ) )

    ! Correct for the case where the parcel crosses the saturation threshold
    ! during the test ascent, and ceases to be unstable after it has crossed...
    if ( i_sat(ic2)==i_next(ic2)-1 .and. init_mass_t(ic2,i_type) > zero ) then
      ! If the test ascent crossed its saturation height, and was
      ! unstable before crossing the threshold
      if ( abs( sublevs(ic2,j_height,i_next(ic2))                              &
              - sublevs(ic2,j_height,i_sat(ic2)) )                             &
         > min_precision * sublevs(ic2,j_height,i_sat(ic2)) ) then
        ! If distance between saturation height and next full-level is
        ! large enough for safe division

        ! Calculate the static stability between the saturation point
        ! and the next model-level
        Nsq_next = -( gravity / virt_temp_next_cmpr(ic2) )                     &
                 * ( sublevs(ic2,j_mean_buoy,i_next(ic2))                      &
                   - sublevs(ic2,j_mean_buoy,i_sat(ic2)) )                     &
                 / ( sublevs(ic2,j_height,i_next(ic2))                         &
                   - sublevs(ic2,j_height,i_sat(ic2)) )

        if ( Nsq_next > zero ) then
          ! If test ascent was stable beyond the saturation threshold...

          ! This is mainly expected to happen for downdrafts triggering
          ! from within liquid cloud, but starting with too small
          ! liquid water content to remain saturated at the next level below.

          ! If the parcel ran out of water to evaporate before even reaching
          ! the next model-level interface, we expect the downdraft to
          ! occur entirely within the current model-level.
          ! Such a feature cannot be vertically resolved, so we suppress
          ! initiation in this case.
          ! To yield smooth behaviour, if the parcel ran out of water between
          ! the next model-level interface and the next full level,
          ! we apply a linear ramp to the initiating mass
          ramp = ( sublevs(ic2,j_height,i_sat(ic2))  - next_height(ic) )       &
               / ( sublevs(ic2,j_height,i_next(ic2)) - next_height(ic) )
          ramp = max( min( ramp, one ), zero )
          ! Ramp goes to zero for sat height < next interface height
          !              one  for sat height > next full-level height
          init_mass_t(ic2,i_type) = init_mass_t(ic2,i_type) * ramp

        end if  ! ( Nsq_next > zero )

      end if  ! Safe to divide
    end if  ! ( i_sat(ic2)==i_next(ic2)-1 .AND. Nsq(ic2) < zero )

  end do  ! ic2 = 1, nc

end do  ! i_type = 1, n_conv_types


return
end subroutine calc_init_mass

end module calc_init_mass_mod
