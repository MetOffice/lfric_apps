! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module parcel_type_mod

use cmpr_type_mod, only: cmpr_type
use comorph_constants_mod, only: real_cvprec

implicit none

! Module contains the type definition for a derived type
! structure containing the properties of a convective parcel,
! stored in compression lists of points where convection is
! occuring.
!
! Also contains subroutines that carry out various actions on
! the fields contained in parcel_type structures; eg allocate
! and deallocate the fields, or combine 2 parcels together.


!----------------------------------------------------------------
! Structure for convective parcel properties in compression lists
!----------------------------------------------------------------
! If adding a new field to the parcel, don't forget to also
! add it to parcel_diags_type in parcel_diags_type_mod
! so that it can be output as a diagnostic.
type :: parcel_type

  ! Structure containing compression list info
  type(cmpr_type) :: cmpr

  ! Super-array for main parcel properties
  ! (mass-flux, buoyancy, etc)
  real(kind=real_cvprec), allocatable :: par_super(:,:)

  ! Super-array containing in-plume means of primary fields
  ! (and tracers) carried by the parcel
  real(kind=real_cvprec), allocatable :: mean_super(:,:)
  ! These will be converted back and forth between actual
  ! temperature, winds etc and versions which are
  ! conserved following dry-mass
  ! (cp_tot*temperature, momentum per unit dry-mass).

  ! Super-array containing most intense core properties
  ! of the parcel
  real(kind=real_cvprec), allocatable :: core_super(:,:)

end type parcel_type


!----------------------------------------------------------------
! Addresses of fields in the compressed par_super array
!----------------------------------------------------------------

! Number of fields in the par_super array
integer, parameter :: n_par = 3

! Indices of each parcel field in the super-array
! Note: the code in parcel_combine (in this module)
! assumes that i_massflux_d = 1 and that the addresses
! of all the other fields are positive.

! Flux of dry-mass per unit surface area / kg m-2 s-1
integer, parameter :: i_massflux_d = 1

! Radius of parcel / m
integer, parameter :: i_radius = 2

! Virtual temperature at parcel edge (for constructing buoyancy PDF)
integer, parameter :: i_edge_virt_temp = 3


contains


!----------------------------------------------------------------
! Subroutine to allocate the fields in a parcel_type structure
!----------------------------------------------------------------
subroutine parcel_alloc( l_tracer, n_points, parcel )

use comorph_constants_mod, only: n_tracers, l_par_core
use fields_type_mod, only: n_fields
use cmpr_type_mod, only: cmpr_alloc

implicit none

! Flag for whether the structure needs to store tracers
logical, intent(in) :: l_tracer

! Size to allocate compression arrays to
integer, intent(in) :: n_points

! Parcel structure to be allocated
type(parcel_type), intent(in out) :: parcel

! Total number of mean-fields
! (depends on whether including tracers)
integer :: n_fields_tot

! Allocate space for compression list indices
call cmpr_alloc( parcel % cmpr, n_points )

! Allocate space for main parcel properties
allocate( parcel % par_super( n_points, n_par ) )

! Allocate space for in-parcel mean fields
n_fields_tot = n_fields
if ( l_tracer )  n_fields_tot = n_fields_tot + n_tracers
allocate( parcel % mean_super( n_points, n_fields_tot ) )

! Also allocate array for parcel core fields if needed
if ( l_par_core ) then
  allocate( parcel % core_super( n_points, n_fields_tot ) )
else
  ! Minimal allocation if not needed
  allocate( parcel % core_super(1,1) )
end if

return
end subroutine parcel_alloc


!----------------------------------------------------------------
! Subroutine to deallocate the fields in a parcel_type structure
!----------------------------------------------------------------
subroutine parcel_dealloc( parcel )

use cmpr_type_mod, only: cmpr_dealloc

implicit none

type(parcel_type), intent(in out) :: parcel

! Deallocate parcel core fields
deallocate( parcel % core_super )

! Deallocate in-parcel mean fields
deallocate( parcel % mean_super )

! Deallocate main parcel properties
deallocate( parcel % par_super )

! Deallocate compression list indices
call cmpr_dealloc( parcel % cmpr )

return
end subroutine parcel_dealloc


!----------------------------------------------------------------
! Subroutine to copy data from one parcel structure
! to another in various different ways
!----------------------------------------------------------------
subroutine parcel_copy( l_tracer, parcel, parcel_new, index_ic )

use comorph_constants_mod, only: n_tracers, l_par_core
use fields_type_mod, only: n_fields

implicit none

! Flag for whether to include tracers
logical, intent(in) :: l_tracer

! Parcel structure from which data is to be copied
type(parcel_type), intent(in) :: parcel
! Parcel structure to copy data into
type(parcel_type), intent(in out) :: parcel_new

! List of compression indices to copy from
integer, optional, intent(in) :: index_ic                                      &
                                 ( parcel_new % cmpr % n_points )

! Total number of fields, including tracers
integer :: n_fields_tot

! Loop counters
integer :: ic, i_field


! Set number of fields in parcel % mean
n_fields_tot = n_fields
if ( l_tracer )  n_fields_tot = n_fields_tot + n_tracers

! If sub-compression index list is input, use that
if ( present(index_ic) ) then

  do ic = 1, parcel_new % cmpr % n_points
    parcel_new % cmpr % index_i(ic)                                            &
      = parcel % cmpr % index_i(index_ic(ic))
    parcel_new % cmpr % index_j(ic)                                            &
      = parcel % cmpr % index_j(index_ic(ic))
  end do
  do i_field = 1, n_par
    do ic = 1, parcel_new % cmpr % n_points
      parcel_new % par_super(ic,i_field)                                       &
        = parcel % par_super(index_ic(ic),i_field)
    end do
  end do
  do i_field = 1, n_fields_tot
    do ic = 1, parcel_new % cmpr % n_points
      parcel_new % mean_super(ic,i_field)                                      &
        = parcel % mean_super(index_ic(ic),i_field)
    end do
  end do
  if ( l_par_core ) then
    do i_field = 1, n_fields_tot
      do ic = 1, parcel_new % cmpr % n_points
        parcel_new % core_super(ic,i_field)                                    &
          = parcel % core_super(index_ic(ic),i_field)
      end do
    end do
  end if

  ! No sub-compression indices given...
else

  ! Straight copy of a compressed parcel
  parcel_new % cmpr % n_points = parcel % cmpr % n_points
  do ic = 1, parcel_new % cmpr % n_points
    parcel_new % cmpr % index_i(ic)                                            &
      = parcel % cmpr % index_i(ic)
    parcel_new % cmpr % index_j(ic)                                            &
      = parcel % cmpr % index_j(ic)
  end do
  do i_field = 1, n_par
    do ic = 1, parcel_new % cmpr % n_points
      parcel_new % par_super(ic,i_field)                                       &
        = parcel % par_super(ic,i_field)
    end do
  end do
  do i_field = 1, n_fields_tot
    do ic = 1, parcel_new % cmpr % n_points
      parcel_new % mean_super(ic,i_field)                                      &
        = parcel % mean_super(ic,i_field)
    end do
  end do
  if ( l_par_core ) then
    do i_field = 1, n_fields_tot
      do ic = 1, parcel_new % cmpr % n_points
        parcel_new % core_super(ic,i_field)                                    &
          = parcel % core_super(ic,i_field)
      end do
    end do
  end if

end if


return
end subroutine parcel_copy


!----------------------------------------------------------------
! Subroutine to compress parcel properties within the existing
! arrays, overwriting to make a desired subset of data contiguous
!----------------------------------------------------------------
subroutine parcel_compress( l_tracer, parcel )

use comorph_constants_mod, only: n_tracers, zero, l_par_core
use fields_type_mod, only: n_fields

implicit none

! Flag for whether the structure needs to store tracers
logical, intent(in) :: l_tracer

! Parcel structure whose fields need to be compressed
type(parcel_type), intent(in out) :: parcel

! Count for number of points still containing active parcels
integer :: n_points_new

! List of indices of non-terminated points
integer, allocatable :: index_ic(:)

! Total number of fields, including tracers if needed
integer :: n_fields_tot

! Loop counters
integer :: ic, ic_first, i_field


! Set number of fields in parcel % mean
n_fields_tot = n_fields
if ( l_tracer )  n_fields_tot = n_fields_tot + n_tracers

! See if any points are not to be included
ic_first = 0
over_n_points: do ic = 1, parcel % cmpr % n_points
  ! Search until we find at least one terminated point,
  ! then exit the loop.
  if ( .not. parcel % par_super(ic,i_massflux_d) > zero ) then
    ! Save the index of the first terminated point
    ic_first = ic
    exit over_n_points
  end if
end do over_n_points

! If any points have terminated
if ( ic_first > 0 ) then

  ! Set number of unterminated points already found
  n_points_new = ic_first - 1

  ! If there might be any more unterminated points after ic_first
  if ( ic_first < parcel % cmpr % n_points ) then

    ! Store a list of indices of points which haven't terminated
    ! (only need to find indices of unterminated points that are
    !  after the first terminated point ic_first, as all earlier
    !  unterminated points can remain where they are).
    allocate( index_ic( ic_first : parcel % cmpr % n_points ) )
    do ic = ic_first+1, parcel % cmpr % n_points
      if ( parcel % par_super(ic,i_massflux_d) > zero ) then
        n_points_new = n_points_new + 1
        index_ic(n_points_new) = ic
      end if
    end do

  end if

  ! Set new number of active points
  parcel % cmpr % n_points = n_points_new

  ! If only a fraction of points have terminated, and
  ! not all at the end of the list (and therefore some
  ! data needs to be rearranged)
  if ( n_points_new >= ic_first ) then

    ! Compress the arrays in-situ, leaving junk in the
    ! elements with ic > n_points_new...

    ! Note that ic_first (saved earlier) is the index
    ! of the first terminated point in the existing list,
    ! and therefore corresponds to the index of the first
    ! point that needs to be moved to form the
    ! recompressed list.

    ! Indices of points in the full 2-D grid
    do ic = ic_first, n_points_new
      parcel % cmpr % index_i(ic) = parcel % cmpr % index_i(index_ic(ic))
      parcel % cmpr % index_j(ic) = parcel % cmpr % index_j(index_ic(ic))
    end do

    ! Fields stored in the parcel super-array
    do i_field = 1, n_par
      do ic = ic_first, n_points_new
        parcel % par_super(ic,i_field)                                         &
          = parcel % par_super(index_ic(ic),i_field)
      end do
    end do

    ! Fields stored in the parcel mean fields super-array
    do i_field = 1, n_fields_tot
      do ic = ic_first, n_points_new
        parcel % mean_super(ic,i_field)                                        &
          = parcel % mean_super(index_ic(ic),i_field)
      end do
    end do

    ! Fields stored in the parcel core fields super-array
    if ( l_par_core ) then
      do i_field = 1, n_fields_tot
        do ic = ic_first, n_points_new
          parcel % core_super(ic,i_field)                                      &
            = parcel % core_super(index_ic(ic),i_field)
        end do
      end do
    end if

  end if  ! ( n_points_new >= ic_first )

end if  ! ( ic_first > 0 )


return
end subroutine parcel_compress


!----------------------------------------------------------------
! Subroutine to expand parcel properties within the existing
! arrays, to put them on a new, larger compression list
!----------------------------------------------------------------
subroutine parcel_expand( l_tracer, n_points_new, index_ic, parcel )

use comorph_constants_mod, only: n_tracers, zero, l_par_core
use fields_type_mod, only: n_fields

implicit none

! Flag for whether the structure needs to store tracers
logical, intent(in) :: l_tracer

! Number of points in the new larger compression list
integer, intent(in) :: n_points_new

! Parcel structure whose fields need to be exanded
type(parcel_type), intent(in out) :: parcel

! Indices for transferring data from the old to the new compression list
integer, intent(in) :: index_ic( parcel % cmpr % n_points )

! Total number of fields, including tracers if needed
integer :: n_fields_tot

! Number of points in the new compression list which will be vacant,
! and their indices
integer :: n_vacant
integer :: index_vacant(n_points_new)

! Loop counters
integer :: ic, ic2, ic_first, i_field

! Set number of fields in parcel % mean
n_fields_tot = n_fields
if ( l_tracer )  n_fields_tot = n_fields_tot + n_tracers

! See if any points need to be moved
ic_first = 0
over_n_points: do ic = 1, parcel % cmpr % n_points
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

  ! Fields stored in the parcel super-array
  do i_field = 1, n_par
    do ic = parcel % cmpr % n_points, ic_first, -1
      parcel % par_super(index_ic(ic),i_field)                                 &
        = parcel % par_super(ic,i_field)
    end do
  end do

  ! Fields stored in the parcel mean fields super-array
  do i_field = 1, n_fields_tot
    do ic = parcel % cmpr % n_points, ic_first, -1
      parcel % mean_super(index_ic(ic),i_field)                                &
        = parcel % mean_super(ic,i_field)
    end do
  end do

  ! Fields stored in the parcel core fields super-array
  if ( l_par_core ) then
    do i_field = 1, n_fields_tot
      do ic = parcel % cmpr % n_points, ic_first, -1
        parcel % core_super(index_ic(ic),i_field)                              &
          = parcel % core_super(ic,i_field)
      end do
    end do
  end if

end if  ! ( ic_first > 0 )

! Find the indices of points left vacant by the expansion.
n_vacant = 0
do ic2 = 1, n_points_new
  index_vacant(ic2) = 0
end do
do ic = 1, parcel % cmpr % n_points
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

  do i_field = 1, n_par
    do ic = 1, n_vacant
      parcel % par_super(index_vacant(ic),i_field) = zero
    end do
  end do
  do i_field = 1, n_fields_tot
    do ic = 1, n_vacant
      parcel % mean_super(index_vacant(ic),i_field) = zero
    end do
  end do
  if ( l_par_core ) then
    do i_field = 1, n_fields_tot
      do ic = 1, n_vacant
        parcel % core_super(index_vacant(ic),i_field) = zero
      end do
    end do
  end if

end if  ! ( n_vacant > 0 )

return
end subroutine parcel_expand


!---------------------------------------------------------------
! Subroutine to initialise parcel properties to zero
!---------------------------------------------------------------
subroutine parcel_init_zero( l_tracer, parcel )

use comorph_constants_mod, only: n_tracers, zero, l_par_core
use fields_type_mod, only: n_fields

implicit none

! Flag for whether the structure contains tracers
logical, intent(in) :: l_tracer

! Structure containing fields to be initialised to zero
type(parcel_type), intent(in out) :: parcel

! Total number of fields, including tracers if needed
integer :: n_fields_tot

! Loop counters
integer :: ic, i_field

! Set number of fields in parcel % mean
n_fields_tot = n_fields
if ( l_tracer )  n_fields_tot = n_fields_tot + n_tracers

! Initialise parcel fields to zero
do i_field = 1, n_par
  do ic = 1, parcel % cmpr % n_points
    parcel % par_super(ic,i_field) = zero
  end do
end do

! Initialise parcel mean-fields to zero
do i_field = 1, n_fields_tot
  do ic = 1, parcel % cmpr % n_points
    parcel % mean_super(ic,i_field) = zero
  end do
end do

! Initialise parcel core fields to zero
if ( l_par_core ) then
  do i_field = 1, n_fields_tot
    do ic = 1, parcel % cmpr % n_points
      parcel % core_super(ic,i_field) = zero
    end do
  end do
end if

return
end subroutine parcel_init_zero


!----------------------------------------------------------------
! Subroutine to combine parcel properties into a merged parcel
!----------------------------------------------------------------
subroutine parcel_combine( l_tracer, l_down, index_ic,                         &
                           parcel_a, parcel_m )

use comorph_constants_mod, only: n_tracers, real_cvprec, zero, one,            &
                                 l_par_core
use fields_type_mod, only: n_fields, i_temperature, i_q_vap,                   &
                           i_qc_first, i_qc_last
use calc_virt_temp_mod, only: calc_virt_temp

implicit none

! Flag for whether tracers are included
logical, intent(in) :: l_tracer

! Flag for downdraft versus updraft
logical, intent(in) :: l_down

! Input properties of one of the parcels to combine
type(parcel_type), intent(in) :: parcel_a

! in:  properties of the other parcel to combine
! out: combined merged parcel properties
type(parcel_type), intent(in out) :: parcel_m

! Index list for referencing the parcel_m compression list
! from parcel_a
integer, intent(in) :: index_ic( parcel_a % cmpr % n_points )

! Total number of fields, including tracers if needed
integer :: n_fields_tot

! Weight to apply to parcel_a's properties
real(kind=real_cvprec) :: weight_a( parcel_a % cmpr % n_points )
! Weight to apply to the existing properties of parcel m
real(kind=real_cvprec) :: weight_m( parcel_a % cmpr % n_points )

! Virtual temperature of the cores of parcels a and m
real(kind=real_cvprec) :: core_a_virt_temp                                     &
                          ( parcel_a % cmpr % n_points )
real(kind=real_cvprec) :: core_m_virt_temp                                     &
                          ( parcel_m % cmpr % n_points )

! Store for indices of points where parcel core to be overwritten
integer :: n_core
integer :: index_ic_core( parcel_a % cmpr % n_points )

! Loop counters
integer :: ic, ic2, ic3, i_field


! Set number of fields in parcel % mean_super
n_fields_tot = n_fields
if ( l_tracer )  n_fields_tot = n_fields_tot + n_tracers

! Need to use the input index list index_ic to transfer
! data from the parcel_a compression list to the
! parcel_m compression list...

! Add the massflux of parcel a onto parcel_m
do ic = 1, parcel_a % cmpr % n_points
  ic2 = index_ic(ic)
  parcel_m % par_super(ic2,i_massflux_d)                                       &
    = parcel_m % par_super(ic2,i_massflux_d)                                   &
    + parcel_a % par_super(ic,i_massflux_d)
end do

! Calculate weight to apply to parcel_a
do ic = 1, parcel_a % cmpr % n_points
  weight_a(ic) = parcel_a % par_super(ic,i_massflux_d)                         &
               / parcel_m % par_super(index_ic(ic),i_massflux_d)
end do
! NOTE: to achieve bit-reproducibility on different
! proc decompositions, it is vital that at points where
! the input parcel_m had zero mass-flux, (and hence
! massflux_m = massflux_a, now that we've added a to m),
! the weight stored in weight_a comes out at exactly 1.0.
! However, on some compilers, dividing one number by
! another which is numerically identical is not guaranteed
! to yield exactly 1.0.
! Hence need an explicit check to ensure this happens:
do ic = 1, parcel_a % cmpr % n_points
  if ( parcel_a % par_super(ic,i_massflux_d)                                   &
    == parcel_m % par_super(index_ic(ic),i_massflux_d) )                       &
      weight_a(ic) = one
end do

! Weight to apply to existing parcel_m properties is 1 - weight_a
do ic = 1, parcel_a % cmpr % n_points
  weight_m(ic) = one - weight_a(ic)
end do

! Take mass-flux-weighted-mean of the parcel radii
do ic = 1, parcel_a % cmpr % n_points
  ic2 = index_ic(ic)
  parcel_m % par_super(ic2,i_radius)                                           &
    = weight_m(ic) * parcel_m % par_super(ic2,i_radius)                        &
    + weight_a(ic) * parcel_a % par_super(ic,i_radius)
end do

! Parcel edge virtual temperature;
! for updrafts, use the lower of the 2 values; for downdrafts, use the higher
! (also always use the new value from parcel_a at points where there is
!  nothing in parcel_m yet, as indicated by edge Tv = 0.0).
if ( l_down ) then
  do ic = 1, parcel_a % cmpr % n_points
    ic2 = index_ic(ic)
    if ( parcel_a % par_super(ic,i_edge_virt_temp)                             &
       > parcel_m % par_super(ic2,i_edge_virt_temp) .or.                       &
         parcel_m % par_super(ic2,i_edge_virt_temp) == zero ) then
      parcel_m % par_super(ic2,i_edge_virt_temp)                               &
        = parcel_a % par_super(ic,i_edge_virt_temp)
    end if
  end do
else
  do ic = 1, parcel_a % cmpr % n_points
    ic2 = index_ic(ic)
    if ( parcel_a % par_super(ic,i_edge_virt_temp)                             &
       < parcel_m % par_super(ic2,i_edge_virt_temp) .or.                       &
         parcel_m % par_super(ic2,i_edge_virt_temp) == zero ) then
      parcel_m % par_super(ic2,i_edge_virt_temp)                               &
        = parcel_a % par_super(ic,i_edge_virt_temp)
    end if
  end do
end if

! Parcel mean properties; add mass-flux-weighted contribution
do i_field = 1, n_fields_tot
  do ic = 1, parcel_a % cmpr % n_points
    ic2 = index_ic(ic)
    parcel_m % mean_super(ic2,i_field)                                         &
      = weight_m(ic) * parcel_m % mean_super(ic2,i_field)                      &
      + weight_a(ic) * parcel_a % mean_super(ic,i_field)
  end do
end do

! Parcel core properties...
if ( l_par_core ) then
  ! These are combined by choosing the max
  ! between the 2 parcels, rather than averaging them together.

  ! Choose properties from parcel with the higher buoyancy.
  ! NOTE: this calculation relies on the fact that the parcel
  ! core properties are not in conserved variable form at this
  ! point, whereas the parcel mean properties are.

  ! Calculate core virtual temperature of parcels a and m
  call calc_virt_temp( parcel_a % cmpr % n_points,                             &
                       size(parcel_a % cmpr % index_i),                        &
                       parcel_a % core_super(:,i_temperature),                 &
                       parcel_a % core_super(:,i_q_vap),                       &
                       parcel_a % core_super(:,i_qc_first:i_qc_last),          &
                       core_a_virt_temp )
  call calc_virt_temp( parcel_m % cmpr % n_points,                             &
                       size(parcel_m % cmpr % index_i),                        &
                       parcel_m % core_super(:,i_temperature),                 &
                       parcel_m % core_super(:,i_q_vap),                       &
                       parcel_m % core_super(:,i_qc_first:i_qc_last),          &
                       core_m_virt_temp )

  ! Find points where parcel_a has a higher virtual temperature
  ! than parcel_m (or lower in the case of downdrafts), or
  ! where nothing has been added to parcel_m yet (as indicated
  ! by zero virtual temperature).
  n_core = 0
  if ( l_down ) then
    do ic = 1, parcel_a % cmpr % n_points
      ic2 = index_ic(ic)
      if ( core_a_virt_temp(ic) <= core_m_virt_temp(ic2) .or.                  &
           core_m_virt_temp(ic2) == zero ) then
        n_core = n_core + 1
        index_ic_core(n_core) = ic
      end if
    end do
  else
    do ic = 1, parcel_a % cmpr % n_points
      ic2 = index_ic(ic)
      if ( core_a_virt_temp(ic) >= core_m_virt_temp(ic2) .or.                  &
           core_m_virt_temp(ic2) == zero ) then
        n_core = n_core + 1
        index_ic_core(n_core) = ic
      end if
    end do
  end if

  ! Overwrite parcel core properties in m with a at points where
  ! a is more buoyant than m
  if ( n_core > 0 ) then
    do i_field = 1, n_fields_tot
      do ic3 = 1, n_core
        ic = index_ic_core(ic3)
        ic2 = index_ic(ic)
        parcel_m % core_super(ic2,i_field)                                     &
          =  parcel_a % core_super(ic,i_field)
      end do
    end do
  end if

end if  ! ( l_par_core )


return
end subroutine parcel_combine


!----------------------------------------------------------------
! Subroutine to check for bad values in parcel properties
!----------------------------------------------------------------
subroutine parcel_check_bad_values( parcel, n_fields_tot, k,                   &
                                    where_string )

use comorph_constants_mod, only: name_length, l_par_core
use fields_type_mod, only: field_names, field_positive
use check_bad_values_mod, only: check_bad_values_cmpr

implicit none

! Parcel whose fields are to be checked
type(parcel_type), intent(in) :: parcel

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


! Check mass-flux, parcel radius, and edge virtual temperature
l_positive = .true.
field_name = "massflux_d"
call check_bad_values_cmpr( parcel % cmpr, k,                                  &
                            parcel % par_super(:,i_massflux_d),                &
                            where_string, field_name, l_positive)
field_name = "radius"
call check_bad_values_cmpr( parcel % cmpr, k,                                  &
                            parcel % par_super(:,i_radius),                    &
                            where_string, field_name, l_positive)
field_name = "edge_virt_temp"
call check_bad_values_cmpr( parcel % cmpr, k,                                  &
                            parcel % par_super(:,i_edge_virt_temp),            &
                            where_string, field_name, l_positive)

! Check mean primary fields
do i_field = 1, n_fields_tot
  field_name = "mean_" // trim(adjustl(field_names(i_field)))
  call check_bad_values_cmpr( parcel % cmpr, k,                                &
                              parcel % mean_super(:,i_field),                  &
                              where_string, field_name,                        &
                              field_positive(i_field) )
end do

! Check parcel core fields if used
if ( l_par_core ) then
  do i_field = 1, n_fields_tot
    field_name = "core_" // trim(adjustl(field_names(i_field)))
    call check_bad_values_cmpr( parcel % cmpr, k,                              &
                                parcel % core_super(:,i_field),                &
                                where_string, field_name,                      &
                                field_positive(i_field) )
  end do
end if


return
end subroutine parcel_check_bad_values


end module parcel_type_mod
