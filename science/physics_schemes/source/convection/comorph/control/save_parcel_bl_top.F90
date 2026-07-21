! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module save_parcel_bl_top_mod

implicit none

contains


! Subroutine to save the mass-flux and in-parcel mean properties
! of the parcel at the first level above the boundary-layer top;
! for use in homogenizing the resolved-scale source terms
! below the BL-top
subroutine save_parcel_bl_top( n_fields_tot, k, dk, grid, turb,                &
                               par_conv, par_bl_top )

use comorph_constants_mod, only: zero
use cmpr_type_mod, only: cmpr_alloc
use grid_type_mod, only: grid_type
use turb_type_mod, only: turb_type
use parcel_type_mod, only: parcel_type, i_massflux_d
use set_l_within_bl_mod, only: set_l_within_bl

implicit none

! Number of primary fields
integer, intent(in) :: n_fields_tot

! Current full model-level
integer, intent(in) :: k
! Model-level increment (+1 for updrafts, -1 for downdrafts)
integer, intent(in) :: dk

! Structure containing grid fields
type(grid_type), intent(in) :: grid

! Structure containing the boundary-layer height
type(turb_type), intent(in) :: turb

! Current parcel properties
type(parcel_type), intent(in) :: par_conv

! Saved parcel properties at the boundary-layer top
type(parcel_type), intent(in out) :: par_bl_top

! Flags for whether current and next levels are within the boundary-layer
logical :: l_within_bl_k    ( par_conv % cmpr % n_points )
logical :: l_within_bl_kpdk ( par_conv % cmpr % n_points )

! Points where parcel has just crossed the BL-top
integer :: nc
integer :: index_ic ( par_conv % cmpr % n_points )

! Array bounds
integer :: lb_h(3), ub_h(3), lb_t(2), ub_t(2)

! Loop counters
integer :: ic, ic0, i_field


! Check whether the current and next full levels are within the
! boundary-layer
lb_h = lbound( grid % height_half )
ub_h = ubound( grid % height_half )
lb_t = lbound( turb % z_bl_top )
ub_t = ubound( turb % z_bl_top )
call set_l_within_bl( par_conv % cmpr, k,                                      &
                      lb_h, ub_h, grid % height_half,                          &
                      lb_t, ub_t, turb % z_bl_top,                             &
                      l_within_bl_k )
call set_l_within_bl( par_conv % cmpr, k+dk,                                   &
                      lb_h, ub_h, grid % height_half,                          &
                      lb_t, ub_t, turb % z_bl_top,                             &
                      l_within_bl_kpdk )

! Find points where the parcel just crossed the BL-top
nc = 0
do ic = 1, par_conv % cmpr % n_points
  ! If current and next level are either side of the BL-top
  if ( ( l_within_bl_kpdk(ic) .neqv. l_within_bl_k(ic) )                       &
       ! TEMPORARY CODE TO PRESERVE KGO; TO BE REMOVED SOON
       ! (only include points where mass-flux +ive at end of level-step)
       .and. par_conv % par_super(ic,i_massflux_d) > zero  ) then
    ! Increment counter and save address of this point
    nc = nc + 1
    index_ic(nc) = ic
  end if
end do

! If any points just crossed the BL-top
if ( nc > 0 ) then

  ! Save number of points that have reached the boundary-layer
  ! top at this level
  par_bl_top % cmpr % n_points = nc
  ! Allocate arrays
  call cmpr_alloc( par_bl_top % cmpr, nc )
  ! Note: we don't need all the parcel fields here; only
  ! mass-flux and in-parcel means, so doing a
  ! custom allocate instead of calling parcel_alloc:
  allocate( par_bl_top % par_super(nc,i_massflux_d:i_massflux_d))
  allocate( par_bl_top % mean_super(nc,n_fields_tot) )

  ! Copy i,j indices of points
  do ic0 = 1, nc
    ic = index_ic(ic0)
    par_bl_top % cmpr % index_i(ic0)                                           &
      = par_conv % cmpr % index_i(ic)
    par_bl_top % cmpr % index_j(ic0)                                           &
      = par_conv % cmpr % index_j(ic)
  end do

  ! Copy mass-flux
  do ic0 = 1, nc
    par_bl_top % par_super(ic0,i_massflux_d)                                   &
      = par_conv % par_super(index_ic(ic0),i_massflux_d)
  end do

  ! Copy parcel properties
  do i_field = 1, n_fields_tot
    do ic0 = 1, nc
      par_bl_top % mean_super( ic0, i_field )                                  &
        = par_conv % mean_super( index_ic(ic0), i_field )
    end do
  end do

end if  ! ( nc > 0 )


return
end subroutine save_parcel_bl_top

end module save_parcel_bl_top_mod
