! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module cmpr_type_mod

implicit none

!----------------------------------------------------------------
! Structure to store compression list information within various
! other structure types that store compressed data
!----------------------------------------------------------------
type :: cmpr_type
  ! Number of compression list points on current model-level
  integer :: n_points
  ! i,j indices of the points
  integer, allocatable :: index_i(:)
  integer, allocatable :: index_j(:)
end type cmpr_type


contains


!----------------------------------------------------------------
! Subroutine to allocate the compression index lists
!----------------------------------------------------------------
subroutine cmpr_alloc( cmpr, n_points )
implicit none
type(cmpr_type), intent(in out) :: cmpr
integer, intent(in) :: n_points

! Note: the allocated size of the compression index arrays is
! passed in as an argument here and is not always the same as the
! compression list length stored in the structure
! (cmpr % n_points).
! This is because, once allocated, it is often advantageous to
! keep re-using the same structure for different compression
! lists of different lengths, and save having to repeatedly
! deallocate and reallocate the arrays
! (provided the initial allocation uses the max size they will
!  ever need).
allocate( cmpr % index_i( n_points ) )
allocate( cmpr % index_j( n_points ) )

return
end subroutine cmpr_alloc


!----------------------------------------------------------------
! Subroutine to deallocate the compression index lists
!----------------------------------------------------------------
subroutine cmpr_dealloc( cmpr )
implicit none
type(cmpr_type), intent(in out) :: cmpr

deallocate( cmpr % index_j )
deallocate( cmpr % index_i )

return
end subroutine cmpr_dealloc


!----------------------------------------------------------------
! Subroutine to copy one compression list object into another
!----------------------------------------------------------------
subroutine cmpr_copy( cmpr_in, cmpr_out )
implicit none
type(cmpr_type), intent(in) :: cmpr_in
type(cmpr_type), intent(in out) :: cmpr_out
integer :: ic

cmpr_out % n_points = cmpr_in % n_points
if ( cmpr_out % n_points > 0 ) then
  do ic = 1, cmpr_out % n_points
    cmpr_out % index_i(ic) = cmpr_in % index_i(ic)
    cmpr_out % index_j(ic) = cmpr_in % index_j(ic)
  end do
end if

return
end subroutine cmpr_copy


!----------------------------------------------------------------
! Subroutine to merge 2 compression lists
!----------------------------------------------------------------
subroutine cmpr_merge( cmpr_a, cmpr_b, cmpr_m,                                 &
                       index_ic_a, index_ic_b )

use comorph_constants_mod, only: nx_full

implicit none

! Compression list objects to be merged
type(cmpr_type), intent(in) :: cmpr_a
type(cmpr_type), intent(in) :: cmpr_b

! Output compression list object containing all the points that
! are in either a or b, in order
type(cmpr_type), intent(in out) :: cmpr_m
! Note: this routine assumes the index arrays cmpr_m % index_i/j
! are already allocated to the required size.

! Indices for transfering data from the cmpr_a / cmpr_b
! compression lists into the merged compression list
integer, intent(out) :: index_ic_a( cmpr_a % n_points )
integer, intent(out) :: index_ic_b( cmpr_b % n_points )

! Flag to continue while loop
logical :: l_continue

! Loop counters
integer :: i, i_a, i_b, i_m, n_a, n_b


if ( cmpr_a % n_points == 0 .and. cmpr_b % n_points == 0 ) then
  ! Both lists are empty

  cmpr_m % n_points = 0

else

  if ( cmpr_a % n_points == 0 ) then
    ! List a is empty; just copy list b

    call cmpr_copy( cmpr_b, cmpr_m )

  else if ( cmpr_b % n_points == 0 ) then
    ! List b is empty; just copy list a

    call cmpr_copy( cmpr_a, cmpr_m )

  else
    ! Neither list is empty, so we have actual work to do...

    ! Initialise counters through lists a,b and the merged list
    i_a = 1
    i_b = 1
    i_m = 1
    l_continue = .true.

    ! Keep looping until we reach the end of one of the lists
    do while ( l_continue )

      ! Calculate position of points in the array
      n_a = ( cmpr_a % index_j(i_a) - 1 )*nx_full                              &
          + cmpr_a % index_i(i_a)
      n_b = ( cmpr_b % index_j(i_b) - 1 )*nx_full                              &
          + cmpr_b % index_i(i_b)

      ! Add the point with the smaller index to the next
      ! element of the merged list
      if ( n_a < n_b ) then
        cmpr_m % index_i(i_m) = cmpr_a % index_i(i_a)
        cmpr_m % index_j(i_m) = cmpr_a % index_j(i_a)
        index_ic_a(i_a) = i_m
        i_a = i_a + 1
      else if (  n_b < n_a ) then
        cmpr_m % index_i(i_m) = cmpr_b % index_i(i_b)
        cmpr_m % index_j(i_m) = cmpr_b % index_j(i_b)
        index_ic_b(i_b) = i_m
        i_b = i_b + 1
        ! If both lists contain the same point (n_a=n_b),
        ! just add it once and advance the counters for both lists
      else
        cmpr_m % index_i(i_m) = cmpr_a % index_i(i_a)
        cmpr_m % index_j(i_m) = cmpr_a % index_j(i_a)
        index_ic_a(i_a) = i_m
        index_ic_b(i_b) = i_m
        i_a = i_a + 1
        i_b = i_b + 1
      end if
      ! Increment counter for merged list
      i_m = i_m + 1

      ! If reached end of list a
      if ( i_a > cmpr_a % n_points ) then
        ! Set flag to terminate while loop
        l_continue = .false.
        ! Mop up any remaining points in list b
        if ( i_b <= cmpr_b % n_points ) then
          do i = i_b, cmpr_b % n_points
            cmpr_m % index_i(i_m) = cmpr_b % index_i(i)
            cmpr_m % index_j(i_m) = cmpr_b % index_j(i)
            index_ic_b(i) = i_m
            i_m = i_m + 1
          end do
        end if
        ! If reached end of list b
      else if ( i_b > cmpr_b % n_points ) then
        ! Set flag to terminate while loop
        l_continue = .false.
        ! Mop up any remaining points in list a
        if ( i_a <= cmpr_a % n_points ) then
          do i = i_a, cmpr_a % n_points
            cmpr_m % index_i(i_m) = cmpr_a % index_i(i)
            cmpr_m % index_j(i_m) = cmpr_a % index_j(i)
            index_ic_a(i) = i_m
            i_m = i_m + 1
          end do
        end if
      end if

    end do  ! while

    ! Save number of points in the merged list
    cmpr_m % n_points = i_m - 1

  end if  ! Neither list is empty

end if  ! Not case when both lists are empty

return
end subroutine cmpr_merge


end module cmpr_type_mod
