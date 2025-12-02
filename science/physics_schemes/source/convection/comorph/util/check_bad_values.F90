! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module check_bad_values_mod

implicit none

contains

! Subroutines to check fields for bad values
! (e.g. NaN, Inf, or negative values for positive-only fields),
! for debugging purposes.

!----------------------------------------------------------------
! Bad value check for full 3-D fields
!----------------------------------------------------------------
subroutine check_bad_values_3d( lb, ub, field, where_string,                   &
                                field_name, l_positive, l_half, l_init )

use comorph_constants_mod, only: real_hmprec, real_cvprec, name_length,        &
                                 nx_full, ny_full, k_bot_conv, k_top_conv,     &
                                 k_top_init, i_check_bad_values_3d
use cmpr_type_mod, only: cmpr_type, cmpr_alloc, cmpr_dealloc
use compress_mod, only: compress

implicit none

! Array lower and upper bounds (in case it has halos)
integer, intent(in) :: lb(3), ub(3)

! Full 3-D field to be checked
real(kind=real_hmprec), intent(in) :: field                                    &
                                      ( lb(1):ub(1), lb(2):ub(2), lb(3):ub(3) )

! String indicating where in the code this check is being done
character(len=name_length), intent(in) :: where_string

! Name of the field (for error message)
character(len=name_length), intent(in) :: field_name

! Flag for whether the field is positive-only
! (in which case we check for negative values)
logical, intent(in) :: l_positive

! Flag for the field being checked is on half-levels; in this case,
! the uppermost model-level used is k_top_conv+1, since this is the
! uppermost model-level interface.  Optional and assumed false
! if not present
logical, optional, intent(in) :: l_half

! Flag for the field being checked is only needed on convection
! initiation levels, and so only needs to be checked up to
! level k_top_init instead of k_top_conv
logical, optional, intent(in) :: l_init

! Field values converted to comorph native precision
real(kind=real_cvprec) :: field_cmpr(nx_full*ny_full)

! Compression indices of points, to pass into 1D bad value check routine
type(cmpr_type) :: cmpr

! Loop counters
integer :: i, j, k, k_top, ic


! Set uppermost model-level to check...
k_top = k_top_conv
! For fields only used in the convection initiation calculation,
! only need to check up to k_top_init
if ( present(l_init) ) then
  if ( l_init )  k_top = k_top_init
end if
! Uppermost level to check is 1 level higher for fields on
! half-levels, to ensure we include the uppermost model-level interface.
if ( present(l_half) ) then
  if ( l_half )  k_top = k_top + 1
end if

! Set compression indices for all points on the horizontal domain
cmpr%n_points = nx_full * ny_full
call cmpr_alloc( cmpr, cmpr%n_points )
do j = 1, ny_full
  do i = 1, nx_full
    ic = nx_full*(j-1) + i
    cmpr%index_i(ic) = i
    cmpr%index_j(ic) = j
  end do
end do

! Loop over model-levels
do k = k_bot_conv, k_top

  ! Gather field values and convert to comorph native precision
  call compress( cmpr, lb(1:2), ub(1:2), field(:,:,k), field_cmpr )

  ! Call routine to check for bad values at native precision.
  ! Note: crucially this will catch instances where the value in field
  ! was representable in the host-model at 64-bit, but becomes garbage
  ! when converted to comorph native precision at 32-bit.
  call check_bad_values_cmpr( cmpr, k, field_cmpr, where_string,               &
                              field_name, l_positive,                          &
                              i_check_bad=i_check_bad_values_3d )

end do

! Deallocate compression indices
call cmpr_dealloc( cmpr )


return
end subroutine check_bad_values_3d


!----------------------------------------------------------------
! Bad value check for compressed 1-D fields
!----------------------------------------------------------------
subroutine check_bad_values_cmpr( cmpr, k, field, where_string,                &
                                  field_name, l_positive, i_check_bad )

use comorph_constants_mod, only: real_cvprec, zero, name_length, newline,      &
                                 i_check_bad_values_cmpr
use cmpr_type_mod, only: cmpr_type

implicit none

! Structure containing compression list indices
type(cmpr_type), intent(in) :: cmpr

! Current model-level
integer, intent(in) :: k

! Compressed 1-D field to be checked
real(kind=real_cvprec), intent(in) :: field(cmpr%n_points)

! String indicating where in the code this check is being done
character(len=name_length), intent(in) :: where_string

! Name of the field (for error message)
character(len=name_length), intent(in) :: field_name

! Flag for whether the field is positive-only
! (in which case we check for negative values)
logical, intent(in) :: l_positive

! Optionally override default setting for whether to do warning or fatal error
integer, optional, intent(in) :: i_check_bad

! Flags set true where bad value found
logical :: l_bad(cmpr%n_points)

! Min and max allowed values, outside of which bad value is diagnosed
real(kind=real_cvprec) :: min_val
real(kind=real_cvprec) :: max_val

real(kind=real_cvprec), parameter :: large_number = huge(field)

! Warning vs fatal error switch to use
integer :: i_check_bad_use

! Loop counter
integer :: ic


! Set max allowed value to the max possible floating point real
! (so really just checking for Infinities)
max_val = large_number
if ( l_positive ) then
  ! If the input is supposed to be positive-only, set minimum
  ! allowed value to zero
  min_val = zero
else
  ! Otherwise set to largest possible negative to just check for -Inf
  min_val = -large_number
end if

! Loop over points
do ic = 1, cmpr%n_points

  ! Bad value test; should trip if NaN or Inf or if outside range
  l_bad(ic) =  ( .not. ( field(ic) >= min_val .and.                            &
                         field(ic) <= max_val ) )

end do

if ( present(i_check_bad) ) then
  ! If override for warning vs fatal error switch is input, use that
  i_check_bad_use = i_check_bad
else
  ! By default, this routine uses centrally-held setting for compressed fields
  i_check_bad_use = i_check_bad_values_cmpr
end if

! Call routine to print coordinates and values where bad values found
call check_bad_values_print(                                                   &
           cmpr, k, l_bad, i_check_bad_use,                                    &
           "(" // trim(adjustl(where_string)) // ")"               //newline// &
           "Bad value in field: " // trim(adjustl(field_name)),                &
           field1=field )


return
end subroutine check_bad_values_cmpr


!----------------------------------------------------------------
! Routine to print flagged bad values / raise warning or error
!----------------------------------------------------------------
subroutine check_bad_values_print( cmpr, k, l_bad, i_check_bad, info_str,      &
                                   field1, field2, field3, field4 )

use comorph_constants_mod, only: real_cvprec, newline,                         &
                                 i_check_bad_warn, i_check_bad_fatal
use cmpr_type_mod, only: cmpr_type
use raise_error_mod, only: raise_fatal, raise_warning

implicit none

! Structure storing i,j indices of points in the compression list
type(cmpr_type), intent(in) :: cmpr

! Current model-level index
integer, intent(in) :: k

! Logical flags set true where value is bad and needs printing
logical, intent(in) :: l_bad(cmpr%n_points)

! Integer switch indicating whether to print a warning or raise a fatal error
integer, intent(in) :: i_check_bad

! Character string containing any preceding info to be printed
character(len=*), intent(in) :: info_str

! Up to 4 fields whose values may optionally be printed at bad value points
real(kind=real_cvprec), optional, intent(in) :: field1(cmpr%n_points)
real(kind=real_cvprec), optional, intent(in) :: field2(cmpr%n_points)
real(kind=real_cvprec), optional, intent(in) :: field3(cmpr%n_points)
real(kind=real_cvprec), optional, intent(in) :: field4(cmpr%n_points)

! Max length of string containing coordinates and values of bad values found
integer, parameter :: str_len = 256
! Max length of string for a single grid-point
integer, parameter :: pnt_len = 64
! Number of characters to use to print real values
integer, parameter :: val_dec = 3            ! Number of decimal places
integer, parameter :: val_len = val_dec + 8  ! Full length needed for ES format

! String containing coordinates and values of any bad values found
character(len=str_len) :: bad_values
! String containing coordinates and value for current bad value
character(len=pnt_len) :: current_bad
! String for the value
character(len=val_len) :: val_str
! Format for writing real value
character(len=val_len) :: val_fmt

! Counter for current write position in bad_values string
integer :: n
! Number of characters occupied by current bad value string
integer :: dn

! Loop counter
integer :: ic

character(len=*), parameter :: routinename                                     &
                               = "CHECK_BAD_VALUES_PRINT"


! Initialise
n = 0
bad_values(1:str_len) = " "
! Creat write format for printing real values
write(val_fmt,"(A3,I0,A1,I0,A1)") "(ES", val_len, ".", val_dec, ")"

! Loop over points
over_n_points: do ic = 1, cmpr%n_points

  if ( l_bad(ic) ) then
    ! If this point has a bad value...

    current_bad(1:pnt_len) = " "

    ! Write coordinates of current point
    write(current_bad,"(3(A1,I0),A2)") "(", cmpr%index_i(ic), ",",             &
                                            cmpr%index_j(ic), ",",             &
                                            k, "):"
    current_bad = adjustl(current_bad)

    ! Write field values
    if ( present(field1) ) then
      write(val_str,val_fmt) field1(ic)
      current_bad = trim(current_bad) // val_str
    end if
    if ( present(field2) ) then
      write(val_str,val_fmt) field2(ic)
      current_bad = trim(current_bad) // val_str
    end if
    if ( present(field3) ) then
      write(val_str,val_fmt) field3(ic)
      current_bad = trim(current_bad) // val_str
    end if
    if ( present(field4) ) then
      write(val_str,val_fmt) field4(ic)
      current_bad = trim(current_bad) // val_str
    end if

    ! Find number of characters needed to print the above
    dn = len(trim(current_bad)) + 2  ! + 2 blanks to space out

    if ( n + dn <= str_len ) then
      ! If enough space to add this bad value to the error string, copy it in
      bad_values(n+1:n+dn) = current_bad(1:dn)
    else
      ! No more space in the string; exit
      exit over_n_points
    end if

    ! Increment counter for space used up so far
    n = n + dn

  end if  ! Bad value tests

end do over_n_points

! Error if any bad values found:
if ( n > 0 ) then
  select case(i_check_bad)
  case (i_check_bad_fatal)
    ! Fatal error: kill the run.
    call raise_fatal( routinename, trim(adjustl(info_str))         //newline// &
           "In at least the following (i,j,k) coordinates: "  //               &
           trim(adjustl(bad_values)) )
  case (i_check_bad_warn)
    ! Not fatal error: just print a warning.
    call raise_warning( routinename, trim(adjustl(info_str))       //newline// &
           "In at least the following (i,j,k) coordinates: "  //               &
           trim(adjustl(bad_values)) )
  end select
end if

return
end subroutine check_bad_values_print


end module check_bad_values_mod
