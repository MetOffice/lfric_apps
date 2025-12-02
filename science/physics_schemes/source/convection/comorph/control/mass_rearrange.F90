! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module mass_rearrange_mod

implicit none


contains

! Subroutine to perform column vertical mass rearrangement
! (i.e. apply the "compensatig subsidence" term in the
!  resolved-scale increment calculation).
! Scatters the properties of the current model-level k
! onto the appropriate levels using a cell-averaged
! vertical advection algorithm.
subroutine mass_rearrange( n_points_super, n_fields_tot,                       &
                           ij_first, ij_last, index_ic,                        &
                           cmpr_any, lb_p, ub_p, pressure,                     &
                           fields_k_super, layer_mass, k, k2,                  &
                           layer_mass_k, layer_mass_k2_added,                  &
                           fields_cmpr, comorph_diags )

use comorph_constants_mod, only: real_cvprec, real_hmprec, zero,               &
                                 k_bot_conv, k_top_conv, nx_full, ny_full,     &
                                 newline, min_delta, i_check_bad_values_3d,    &
                                 i_check_bad_warn, i_check_bad_fatal
use cmpr_type_mod, only: cmpr_type, cmpr_alloc, cmpr_dealloc,                  &
                         cmpr_copy
use fields_type_mod, only: ragged_super_type,                                  &
                           fields_k_pressure_adjust
use comorph_diags_type_mod, only: comorph_diags_type
use mass_rearrange_calc_mod, only: mass_rearrange_calc
use calc_pressure_incr_diag_mod, only: calc_pressure_incr_diag
use raise_error_mod, only: raise_warning, raise_fatal

implicit none

! Max number of convecting points (dimensions super-array)
integer, intent(in) :: n_points_super
! Number of primary fields to update
integer, intent(in) :: n_fields_tot

! First and last ij-indices on the current segment
integer, intent(in) :: ij_first
integer, intent(in) :: ij_last
! cmpr_any compression list addresses scattered into a
! common array for referencing from res_source compression list
integer, intent(in) :: index_ic( ij_first : ij_last,                           &
                                 k_bot_conv-1 : k_top_conv+1 )

! Compression indices for current level k
type(cmpr_type), intent(in) :: cmpr_any

! 3-D array of pressure on thermodynamic levels
! (and lower and upper bounds of the array, in case it has halos)
integer, intent(in) :: lb_p(3), ub_p(3)
real(kind=real_hmprec), intent(in) :: pressure                                 &
                          ( lb_p(1):ub_p(1), lb_p(2):ub_p(2), lb_p(3):ub_p(3) )

! In: Compressed fields on level k, updated by
! entrainment and detrainment, but not by vertical rearrangement.
! Gets modified as a work array in this routine
real(kind=real_cvprec), intent(in out) :: fields_k_super                       &
                                ( n_points_super, n_fields_tot )

! Full 3-D array of layer-masses
! (valid for both start and end of the convection step)
real(kind=real_hmprec), intent(in) :: layer_mass                               &
                     ( nx_full, ny_full, k_bot_conv:k_top_conv )

! Current model-level
integer, intent(in) :: k
! Model-level that fields from level k are being scattered into
integer, intent(in out) :: k2(ij_first:ij_last)

! In: layer-mass on level k after removing entrainment and
! adding detrainment.
! Out: gets decremented, as is used as a work variable to
! keep track of how much of the mass from level k we have not
! yet assigned to levels k2.
real(kind=real_cvprec), intent(in out) :: layer_mass_k                         &
                                         (n_points_super)

! Amount of dry-mass so-far added to model-level k2
real(kind=real_cvprec), intent(in out) :: layer_mass_k2_added                  &
                                         (ij_first:ij_last)

! Separate compressed fields on each level
type(ragged_super_type), intent(in out) :: fields_cmpr                         &
                                         (k_bot_conv:k_top_conv)

! Structure containing diagnostic flags and meta-data
type(comorph_diags_type), intent(in out) :: comorph_diags

! Structure storing compression list info for points
! where mass needs to be added
type(cmpr_type) :: cmpr_add

! Number of points (and indices of points) where level k2 is now
! fully assigned, and we need to increment k2 to the next level
integer :: n_points_incr_k2
integer :: index_ic_incr_k2( cmpr_any % n_points )

! Fraction of mass on level k2 to come from level k
! (used as the weight for calculating new mean field values
!  on level k2)
real(kind=real_cvprec) :: frac_k2_from_k( cmpr_any % n_points )

! Layer-mass at level k2, gathered onto cmpr_any points
real(kind=real_cvprec) :: layer_mass_k2( cmpr_any % n_points )

! Pressure on layers k and k2, on the cmpr_any compression list
! (needed for dry adiabatic adjustment of temperature)
real(kind=real_cvprec) :: pressure_k( cmpr_any % n_points )
real(kind=real_cvprec) :: pressure_k2( cmpr_any % n_points )

! Index list for referencing compressed fields on common grid
integer :: ij( cmpr_any % n_points )

! Integer index, used for compressing arrays in situ
integer :: ic_first

real(kind=real_cvprec) :: small_number

! Character strings for adding info to errors
character(len=16) :: str1, str2

! Lower and upper bounds of Lagrangian pressure increment diagnostic array
integer :: lb_d(3), ub_d(3)

! Loop counters
integer :: i, j, ic, ic_k2, ic_incr, i_field

character(len=*), parameter :: routinename = "MASS_REARRANGE"


! Set threshold for decting mass errors that are unlikely to be
! explainable just from precision-level rounding errors
small_number = 10.0_real_cvprec                                                &
               * sqrt(real(k-k_bot_conv+1,real_cvprec))                        &
               * min_delta
! Expected error likely scales with sqrt of number of levels
! over-which the error can've accumulated
! (assuming each level-step gives a random rounding error)

! Find ij indices for referencing the k2 and index_ic arrays
do ic = 1, cmpr_any % n_points
  i = cmpr_any % index_i(ic)
  j = cmpr_any % index_j(ic)
  ij(ic) = nx_full*(j-1)+i
end do

! If this is the first model-level in a convecting layer,
! initialise k2 and mass added at k2
do ic = 1, cmpr_any % n_points
  if ( index_ic(ij(ic),k) > 0 .and.                                            &
       index_ic(ij(ic),k-1) == 0 ) then
    k2(ij(ic)) = k
    layer_mass_k2_added(ij(ic)) = zero
  end if
end do

! Initialise compression indices for points where mass needs to
! be added;
! to start with, this is all convecting points on level k
call cmpr_alloc( cmpr_add, cmpr_any % n_points )
call cmpr_copy( cmpr_any, cmpr_add )

do ic = 1, cmpr_add % n_points
  i = cmpr_add % index_i(ic)
  j = cmpr_add % index_j(ic)
  ! Gather layer-masses from model-level k2
  layer_mass_k2(ic) = real(layer_mass(i,j,k2(ij(ic))),real_cvprec)
  ! Gather pressures at k and k2
  pressure_k(ic)  = real(pressure(i,j,k),real_cvprec)
  pressure_k2(ic) = real(pressure(i,j,k2(ij(ic))),real_cvprec)
end do

! Adjust fields due to the pressure change from k to k2
call fields_k_pressure_adjust( cmpr_add % n_points,                            &
                               n_points_super, n_fields_tot,                   &
                               fields_k_super,                                 &
                               pressure_k, pressure_k2 )

! Compute mass from k to add to k2, and find points
! where the k2 level counter needs to be incremented
call mass_rearrange_calc( cmpr_add, n_points_incr_k2,                          &
                          index_ic_incr_k2,                                    &
                          ij_first, ij_last, ij,                               &
                          layer_mass_k2, layer_mass_k2_added,                  &
                          layer_mass_k, frac_k2_from_k )

! Increment properties at level k2 consistent with adding
! the mass fraction frac_k2_from_k from level k onto level k2
do i_field = 1, n_fields_tot
  do ic = 1, cmpr_add % n_points
    ic_k2 = index_ic(ij(ic),k2(ij(ic)))
    fields_cmpr(k2(ij(ic))) % super(ic_k2,i_field)                             &
      = fields_cmpr(k2(ij(ic))) % super(ic_k2,i_field)                         &
      + fields_k_super(ic,i_field) * frac_k2_from_k(ic)
  end do
end do

! If requested, calculate contribution to mean source-layer
! pressure at k2 from level k
if ( comorph_diags % pressure_incr_env % flag ) then
  lb_d = lbound( comorph_diags % pressure_incr_env % field_3d )
  ub_d = ubound( comorph_diags % pressure_incr_env % field_3d )
  call calc_pressure_incr_diag( cmpr_add, lb_p, ub_p, lb_d, ub_d,              &
                                ij_first, ij_last, k, k2, ij,                  &
                                frac_k2_from_k, pressure,                      &
                                comorph_diags % pressure_incr_env % field_3d )
end if


! If any points still need to add mass to the next level up
do while ( n_points_incr_k2 > 0 )


  ! Compress fields onto points where k2 needs to be incremented.
  ! Since the new set of points is guaranteed to be a subset of
  ! the existing points, and we no longer need the uncompressed
  ! data, we can compress in-situ, avoiding the
  ! need to copy into new arrays...

  ! Only any compression to do if the new list is shorter than
  ! the existing one
  if ( n_points_incr_k2 < cmpr_add % n_points ) then

    ! Find the first compression index which differs
    ic_first = 0
    over_n_points_incr_k2: do ic_incr = 1, n_points_incr_k2
      if ( .not. index_ic_incr_k2(ic_incr) == ic_incr ) then
        ic_first = ic_incr
        exit over_n_points_incr_k2
      end if
    end do over_n_points_incr_k2

    ! Only any compression to do if at least one index differs
    if ( .not. ic_first == 0 ) then

      ! Compress unassigned layer-mass at k
      do ic_incr = ic_first, n_points_incr_k2
        ic = index_ic_incr_k2(ic_incr)
        layer_mass_k(ic_incr) = layer_mass_k(ic)
      end do

      ! Compress fields at k
      do i_field = 1, n_fields_tot
        do ic_incr = ic_first, n_points_incr_k2
          ic = index_ic_incr_k2(ic_incr)
          fields_k_super(ic_incr,i_field)                                      &
            = fields_k_super(ic,i_field)
        end do
      end do

      ! Compress i,j indices
      do ic_incr = ic_first, n_points_incr_k2
        ic = index_ic_incr_k2(ic_incr)
        cmpr_add % index_i(ic_incr) = cmpr_add % index_i(ic)
        cmpr_add % index_j(ic_incr) = cmpr_add % index_j(ic)
        ij(ic_incr) = ij(ic)
      end do

    end if

    ! Set number of points to act on to new smaller value
    cmpr_add % n_points = n_points_incr_k2

  end if  ! ( n_points_incr_k2 < cmpr_add % n_points )

  ! Compress pressure at k2 into the pressure_k array;
  ! the fields in fields_k_super are already adjusted consistent
  ! with the pressure at the old k2.  Next we'll adjust them
  ! consistent with the pressure change from old k2 to new k2.
  do ic_incr = 1, n_points_incr_k2
    ic = index_ic_incr_k2(ic_incr)
    pressure_k(ic_incr) = pressure_k2(ic)
  end do

  ! At all points in the updated list...
  do ic_incr = 1, cmpr_add % n_points
    ! Increment k2
    k2(ij(ic_incr)) = k2(ij(ic_incr)) + 1
    ! Reset counter for mass added at k2 to zero
    layer_mass_k2_added(ij(ic_incr)) = zero
  end do

  ! Check that the new (i,j,k2) points are convecting points
  do ic = 1, cmpr_add % n_points
    ! If this index is not labelled as a convecting point...
    if ( index_ic(ij(ic),k2(ij(ic))) == 0 ) then
      i = cmpr_add % index_i(ic)
      j = cmpr_add % index_j(ic)

      ! Rounding error in floating point arithmetic can cause
      ! very tiny amounts of mass to be left over and spuriously
      ! added to the level above the convection top.

      ! When this happens, shift k2 down to the level below just
      ! to avoid out-of-bounds error.
      k2(ij(ic)) = k2(ij(ic)) - 1

      ! Warn or fatal error if we're trying to add a greater mass fraction
      ! than is expected from rounding errors, as this suggests
      ! something has gone very wrong!
      if ( layer_mass_k(ic) > small_number                                     &
            * real(layer_mass(i,j,k2(ij(ic))),real_cvprec) ) then
        write(str1,"(A1,I3,A1,I3,A1,I3,A1)")                                   &
          "(", i, ",", j, ",", k2(ij(ic)), ")"
        write(str2,"(ES14.6)")                                                 &
          layer_mass_k(ic) / real(layer_mass(i,j,k2(ij(ic))),real_cvprec)
        select case(i_check_bad_values_3d)
        case (i_check_bad_fatal)
          call raise_fatal( routinename,                                       &
               "CoMorph's vertical mass rearrangement calculation " //         &
               "is trying to add mass to a grid point which is "    //newline//&
               "not labelled as a convecting point, and the mass "  //         &
               "error exceeds the tolerance accounting for "        //newline//&
               "expected floating point rounding errors. "          //         &
               "Coordinates (i,j,k) = " // str1                     //newline//&
               "Fractional mass error = " // str2 )
        case (i_check_bad_warn)
          call raise_warning( routinename,                                     &
               "CoMorph's vertical mass rearrangement calculation " //         &
               "is trying to add mass to a grid point which is "    //newline//&
               "not labelled as a convecting point, and the mass "  //         &
               "error exceeds the tolerance accounting for "        //newline//&
               "expected floating point rounding errors. "          //         &
               "Coordinates (i,j,k) = " // str1                     //newline//&
               "Fractional mass error = " // str2 )
        end select

      end if

      ! Reset layer_mass_k to zero (i.e. just ditch the surplus mass).
      layer_mass_k(ic) = zero

    end if
  end do

  do ic = 1, cmpr_add % n_points
    i = cmpr_add % index_i(ic)
    j = cmpr_add % index_j(ic)
    ! Gather layer-masses from new model-level k2
    layer_mass_k2(ic) = real(layer_mass(i,j,k2(ij(ic))),real_cvprec)
    ! Gather pressure at the new k2 level
    pressure_k2(ic) = real(pressure(i,j,k2(ij(ic))), real_cvprec)
  end do

  ! Adjust fields due to pressure change from old k2 to new k2
  call fields_k_pressure_adjust( cmpr_add%n_points,                            &
                                 n_points_super, n_fields_tot,                 &
                                 fields_k_super,                               &
                                 pressure_k, pressure_k2 )

  ! Compute mass from k to add to the new k2, and find points
  ! where the k2 level counter needs to be incremented again
  call mass_rearrange_calc( cmpr_add, n_points_incr_k2,                        &
                            index_ic_incr_k2,                                  &
                            ij_first, ij_last, ij,                             &
                            layer_mass_k2, layer_mass_k2_added,                &
                            layer_mass_k, frac_k2_from_k )

  ! Increment properties at layer k2 consistent with adding
  ! the mass fraction frac_k2_from_k from level k onto level k2
  do i_field = 1, n_fields_tot
    do ic = 1, cmpr_add % n_points
      ic_k2 = index_ic(ij(ic),k2(ij(ic)))
      fields_cmpr(k2(ij(ic))) % super(ic_k2,i_field)                           &
        = fields_cmpr(k2(ij(ic))) % super(ic_k2,i_field)                       &
        + fields_k_super(ic,i_field) * frac_k2_from_k(ic)
    end do
  end do

  ! If requested, calculate contribution to mean source-layer
  ! pressure at k2 from level k
  if ( comorph_diags % pressure_incr_env % flag ) then
    call calc_pressure_incr_diag( cmpr_add, lb_p, ub_p, lb_d, ub_d,            &
                                  ij_first, ij_last, k, k2, ij,                &
                                  frac_k2_from_k, pressure,                    &
                                  comorph_diags % pressure_incr_env % field_3d)
  end if


end do  ! while ( n_points_incr_k2 > 0 )


! Deallocate compression list info ready for next level
call cmpr_dealloc( cmpr_add )


return
end subroutine mass_rearrange

end module mass_rearrange_mod
