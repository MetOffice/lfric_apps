! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module res_source_type_mod

use comorph_constants_mod, only: real_cvprec
use cmpr_type_mod, only: cmpr_type

implicit none

! Module contains the type definition for a derived type
! structure containing the convective resolved-scale
! source terms, stored in compression lists of points where
! convection is active.
!
! Also contains subroutines to allocate and deallocate the fields

!----------------------------------------------------------------
! Structure to store resolved-scale sources/sinks of mass and
! primary variables on a level, in compressed arrays
!----------------------------------------------------------------
type :: res_source_type

  ! Compression list properties
  type(cmpr_type) :: cmpr

  ! Super-array containing entrainment and detrainment of dry-mass
  ! on the current level / kg m-2 s-1
  ! (net injection of dry-mass by convection per s
  !  per unit surface area is det minus ent).
  real(kind=real_cvprec), allocatable :: res_super(:,:)

  ! Super-array containing resolved-scale source terms for each
  ! primary field, and tracers if applicable.
  ! These store the sum of contributions from
  ! detrainment and entrainment
  !   det * phi_det  -  ent * phi_ent
  ! plus sources / sinks due to non-mass-exchanging processes
  ! (pressure gradient forces, precipitation fall).
  real(kind=real_cvprec), allocatable :: fields_super(:,:)

  ! Super-array containing contributions to diagnosed convective
  ! cloud fractions and condensate mixing-ratios
  real(kind=real_cvprec), allocatable :: convcloud_super(:,:)

end type res_source_type


! Number of fields in the res_super array
integer, parameter :: n_res = 2

! Addresses of entrainment & detrainment in the res_super array
integer, parameter :: i_ent = 1
integer, parameter :: i_det = 2


contains


!----------------------------------------------------------------
! Subroutine to allocate fields in a res_source_type structure
!----------------------------------------------------------------
subroutine res_source_alloc( l_tracer, n_points, res_source )

use fields_type_mod, only: n_fields
use cloudfracs_type_mod, only: n_convcloud
use comorph_constants_mod, only: n_tracers
use cmpr_type_mod, only: cmpr_alloc

implicit none

! Flag for whether the structure needs to store tracers
logical, intent(in) :: l_tracer

! Size to allocate compression arrays to
integer, intent(in) :: n_points

! res_source structure to be allocated
type(res_source_type), intent(in out) :: res_source

! Total number of fields, including tracers if applicable
integer :: n_fields_tot

! Allocate space for compression list indices
call cmpr_alloc( res_source % cmpr, n_points )

! Allocate array for entrainment and detrainment
allocate( res_source % res_super( n_points, n_res ) )

! Allocate super-array for resolved-scale source term for
! each primary field and tracer (if applicable)
n_fields_tot = n_fields
if ( l_tracer )  n_fields_tot = n_fields_tot + n_tracers
allocate( res_source % fields_super( n_points, n_fields_tot ) )

! Allocate array for convective cloud fields
allocate( res_source % convcloud_super( n_points, n_convcloud ) )

return
end subroutine res_source_alloc


!----------------------------------------------------------------
! Subroutine to deallocate fields in a res_source_type structure
!----------------------------------------------------------------
subroutine res_source_dealloc( res_source )

use cmpr_type_mod, only: cmpr_dealloc

implicit none

! Structure containing fields to be deallocated
type(res_source_type), intent(in out) :: res_source

! Deallocate convective cloud super-array
deallocate( res_source % convcloud_super )

! Deallocate source term super-array
deallocate( res_source % fields_super )

! Deallocate entrainment, detrainment and cloud array
deallocate( res_source % res_super )

! Deallocate compression list indices
call cmpr_dealloc( res_source % cmpr )

return
end subroutine res_source_dealloc


!----------------------------------------------------------------
! Subroutine to initialise the fields to zero
!----------------------------------------------------------------
subroutine res_source_init_zero( l_tracer, res_source )

use comorph_constants_mod, only: n_tracers, zero
use fields_type_mod, only: n_fields
use cloudfracs_type_mod, only: n_convcloud

implicit none

! Flag for whether the fields structure includes tracers
logical, intent(in) :: l_tracer

! Structure containing fields to be initialised to zero
type(res_source_type), intent(in out) :: res_source

! Total number of fields, including tracers
integer :: n_fields_tot

! Loop counters
integer :: ic, i_field

! Set number of fields in res_source % fields_super
n_fields_tot = n_fields
if ( l_tracer )  n_fields_tot = n_fields_tot + n_tracers

! Initialise entrainement, detrainment and convective cloud
do i_field = 1, n_res
  do ic = 1, res_source % cmpr % n_points
    res_source % res_super(ic,i_field) = zero
  end do
end do

! Initialise the fields super-array
do i_field = 1, n_fields_tot
  do ic = 1, res_source % cmpr % n_points
    res_source % fields_super(ic,i_field) = zero
  end do
end do

! Initialise the convective cloud fields to zero
if ( n_convcloud > 0 ) then
  do i_field = 1, n_convcloud
    do ic = 1, res_source % cmpr % n_points
      res_source % convcloud_super(ic,i_field) = zero
    end do
  end do
end if

return
end subroutine res_source_init_zero


!----------------------------------------------------------------
! Subroutine to expand resolved-scale source-terms within the existing
! arrays, to put them on a new, larger compression list
!----------------------------------------------------------------
subroutine res_source_expand( l_tracer, n_points_new, index_ic, res_source )

use comorph_constants_mod, only: n_tracers, zero

use fields_type_mod, only: n_fields
use cloudfracs_type_mod, only: n_convcloud

implicit none

! Flag for whether the structure needs to store tracers
logical, intent(in) :: l_tracer

! Number of points in the new larger compression list
integer, intent(in) :: n_points_new

! Resolved-scale source-term structure whose fields need to be expanded
type(res_source_type), intent(in out) :: res_source

! Indices for transferring data from the old to the new compression list
integer, intent(in) :: index_ic( res_source % cmpr % n_points )

! Total number of fields, including tracers if needed
integer :: n_fields_tot

! Number of points in the new compression list which will be vacant,
! and their indices
integer :: n_vacant
integer :: index_vacant(n_points_new)

! Loop counters
integer :: ic, ic2, ic_first, i_field

! Set number of fields in res_source % fields_super
n_fields_tot = n_fields
if ( l_tracer )  n_fields_tot = n_fields_tot + n_tracers

! See if any points need to be moved
ic_first = 0
over_n_points: do ic = 1, res_source % cmpr % n_points
  ! Search until we find at least one point whose index changes
  ! then exit the loop.
  if ( .not. index_ic(ic) == ic ) then
    ! Save the index of the first point to move
    ic_first = ic
    exit over_n_points
  end if
end do over_n_points

! If any points need to be moved
if ( ic_first > 0 ) then

  ! For each super-array, loop through the points that need
  ! to be moved (only from ic_first onwards) and use index_ic to
  ! transfer the data to its new location.
  ! Note we need to loop backwards through the points to avoid
  ! overwriting data that hasn't been transfered yet.

  ! Fields stored in the res super-array
  do i_field = 1, n_res
    do ic = res_source % cmpr % n_points, ic_first, -1
      res_source % res_super(index_ic(ic),i_field)                             &
        = res_source % res_super(ic,i_field)
    end do
  end do

  ! Fields stored in the fields super-array
  do i_field = 1, n_fields_tot
    do ic = res_source % cmpr % n_points, ic_first, -1
      res_source % fields_super(index_ic(ic),i_field)                          &
        = res_source % fields_super(ic,i_field)
    end do
  end do

  ! Convective cloud fields
  if ( n_convcloud > 0 ) then
    do i_field = 1, n_convcloud
      do ic = res_source % cmpr % n_points, ic_first, -1
        res_source % convcloud_super(index_ic(ic),i_field)                     &
          = res_source % convcloud_super(ic,i_field)
      end do
    end do
  end if

end if  ! ( ic_first > 0 )

! Find the indices of points left vacant by the expansion.
n_vacant = 0
do ic2 = 1, n_points_new
  index_vacant(ic2) = 0
end do
do ic = 1, res_source % cmpr % n_points
  index_vacant(index_ic(ic)) = 1
end do
do ic2 = 1, n_points_new
  if ( index_vacant(ic2) == 0 ) then
    n_vacant = n_vacant + 1
    index_vacant(n_vacant) = ic2
  end if
end do

! If any vacant points
if ( n_vacant > 0 ) then
  ! Zero all data in the vacant points, for safety.  They will presumably
  ! be populated with new data after the call to this routine.

  do i_field = 1, n_res
    do ic = 1, n_vacant
      res_source % res_super(index_vacant(ic),i_field) = zero
    end do
  end do
  do i_field = 1, n_fields_tot
    do ic = 1, n_vacant
      res_source % fields_super(index_vacant(ic),i_field) = zero
    end do
  end do
  if ( n_convcloud > 0 ) then
    do i_field = 1, n_convcloud
      do ic = 1, n_vacant
        res_source % convcloud_super(index_vacant(ic),i_field) = zero
      end do
    end do
  end if

end if  ! ( n_vacant > 0 )

return
end subroutine res_source_expand


!----------------------------------------------------------------
! Subroutine to combine resolved-scale source terms into
! a merged structure
!----------------------------------------------------------------
subroutine res_source_combine( l_tracer, index_ic,                             &
                               res_source_a, res_source_m )

use comorph_constants_mod, only: n_tracers
use fields_type_mod, only: n_fields
use cloudfracs_type_mod, only: n_convcloud

implicit none

! Flag for whether the fields structure includes tracers
logical, intent(in) :: l_tracer

! Input resolved-scale sourc terms to be added on
type(res_source_type), intent(in) :: res_source_a

! Output combined resolved-scale source terms
type(res_source_type), intent(in out) :: res_source_m

! Index list for referencing the _m compression list from _a
integer, intent(in) :: index_ic( res_source_a % cmpr % n_points )

! Total number of fields, including tracers
integer :: n_fields_tot

! Loop counters
integer :: ic, i_field


! Set number of fields in res_source % fields
n_fields_tot = n_fields
if ( l_tracer )  n_fields_tot = n_fields_tot + n_tracers

! Need to use the input index list index_ic to transfer
! data from the _a compression list to the _m compression list.

! Entrainment, detrainment and cloud
do i_field = 1, n_res
  do ic = 1, res_source_a % cmpr % n_points
    res_source_m % res_super( index_ic(ic), i_field )                          &
      = res_source_m % res_super( index_ic(ic), i_field )                      &
      + res_source_a % res_super( ic, i_field )
  end do
end do

! Source terms for primary fields
do i_field = 1, n_fields_tot
  do ic = 1, res_source_a % cmpr % n_points
    res_source_m % fields_super( index_ic(ic), i_field )                       &
      = res_source_m % fields_super( index_ic(ic), i_field )                   &
      + res_source_a % fields_super( ic, i_field )
  end do
end do

! Convective cloud fields
if ( n_convcloud > 0 ) then
  do i_field = 1, n_convcloud
    do ic = 1, res_source_a % cmpr % n_points
      res_source_m % convcloud_super( index_ic(ic), i_field )                  &
        = res_source_m % convcloud_super( index_ic(ic), i_field )              &
        + res_source_a % convcloud_super( ic, i_field )
    end do
  end do
end if

return
end subroutine res_source_combine


!----------------------------------------------------------------
! Subroutine to check for bad values in resolved-scale sources
!----------------------------------------------------------------
subroutine res_source_check_bad_values( res_source, n_fields_tot,              &
                                        k, where_string )

use comorph_constants_mod, only: name_length
use fields_type_mod, only: field_names
use cloudfracs_type_mod, only: n_convcloud, convcloud_names
use check_bad_values_mod, only: check_bad_values_cmpr

implicit none

! res_source structure whose fields are to be checked
type(res_source_type), intent(in) :: res_source

! Total number of primary fields in super-array (incl tracers)
integer, intent(in) :: n_fields_tot

! Current model-level index
integer, intent(in) :: k

! Character string describing where in the convection scheme
! we are, for constructing error message if bad value found.
character(len=name_length), intent(in) :: where_string

! Name of individual field
character(len=name_length) :: field_name
! Flag for whether field is positive-only
logical :: l_positive

! Loop counter
integer :: i_field


! All source terms can be positive or negative
l_positive = .false.

! Check entrained and detrained mass fields
field_name = "ent_mass_d"
call check_bad_values_cmpr( res_source % cmpr, k,                              &
                            res_source % res_super(:,i_ent),                   &
                            where_string, field_name, l_positive)
field_name = "det_mass_d"
call check_bad_values_cmpr( res_source % cmpr, k,                              &
                            res_source % res_super(:,i_det),                   &
                            where_string, field_name, l_positive)

! Check source terms for primary fields
do i_field = 1, n_fields_tot
  call check_bad_values_cmpr( res_source % cmpr, k,                            &
                              res_source % fields_super                        &
                                           (:,i_field),                        &
                              where_string, field_names(i_field),              &
                              l_positive )
end do

if ( n_convcloud > 0 ) then
  ! Check convective cloud fields...

  ! These must all be positive
  l_positive = .true.

  do i_field = 1, n_convcloud
    call check_bad_values_cmpr( res_source % cmpr, k,                          &
                                res_source % convcloud_super                   &
                                             (:,i_field),                      &
                                where_string,                                  &
                                convcloud_names(i_field),                      &
                                l_positive )
  end do

end if


return
end subroutine res_source_check_bad_values



end module res_source_type_mod
