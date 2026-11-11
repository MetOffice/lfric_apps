! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

module chk_opts_mod
! Description:
! This module is the interface to a utility subroutine
! that checks namelist input variables for sensible
! values based on developer input strings or an integer array
!
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc

use errormessagelength_mod, only: errormessagelength
use umprintmgr,             only: umprint, ummessage, newline
use ereport_mod,            only: ereport

implicit none

private
public :: chk_var, def_src

integer, parameter :: max_n_conditions = 15
integer, parameter :: condition_len    = 20
integer, parameter :: str_len          = 100

character(len=*), parameter :: modname ='CHK_OPTS_MOD'

character(len=2) :: operators(max_n_conditions)
character(len=2) :: operator_string

character(len=errormessagelength) :: errstr
character(len=condition_len) :: values  (max_n_conditions)
character(len=condition_len) :: lbounds (max_n_conditions)
character(len=condition_len) :: ubounds (max_n_conditions)
character(len=condition_len) :: tmp_str

character(len=str_len) :: srcname
character(len=str_len) :: loc_str
character(len=str_len) :: loc_name
character(len=str_len) :: loc_src

character(len=8)  :: result_chk

logical :: pass_chk
integer :: nconditions

! Setup inteface to allow user to call same routinename in order to check
! supported datatypes. In order to pass, the variable must pass any of
! the code owner specified checks, these checks should be consistent with
! the UM rose-metadata
interface chk_var
  module procedure chk_int_str   ! Integers (values/ranges)
  module procedure chk_real_str  ! Reals    (values/ranges)
  module procedure chk_int_arr   ! Integers (values only)
end interface

character(len=str_len) :: def_src = '' ! Calling routine if specified

contains

subroutine chk_int_str (var, var_name, string, report_pass, cmessage)

! Checks an integer variable against a list of criteria
! passed as a string in the form '[Check_A, Check_B, ... ]'.
!
! Recognised checks
!  * ==Integer (This is assumed if no operator is provided)
!  * <Integer
!  * >Integer
!  * >=Integer
!  * <=Integer
!  * Integer lower bound : Integer upper bound
!
! As this set of checks should be for valid values of the variable
! under test, the defensive check is passed if any of the criteria
! are met. If none of the specified checks are met, ereport will be
! called to abort.
!
! Two optional arguments are provided:
! report_pass - (logical) value indicates if any of the checks passed
!               Ereport is not called with this option.
! cmessage    - If Chk_var calls ereport, the variable name and
!               value will be included, along with some information.
!               if a cmessage argument is provided, the message
!               will be appended to chk_var's ereport message.

implicit none

integer,          intent(in) :: var       ! The integer variable to be checked
character(len=*), intent(in) :: var_name  ! The name of the test variable
character(len=*), intent(in) :: string    ! String containg test criteria

! Return flag for pass/failure of test criteria. If this given in argument
! list then ereport will not be called
logical,          intent(in out), optional :: report_pass

! Additonal comments to af to ereport in event of failure of tests
character(len=*), intent(in), optional :: cmessage


! Local variables
character(len=*), parameter :: routinename ='CHK_INT_STR'

integer :: i, index_val, icode
integer :: condition
integer :: lowerbound
integer :: upperbound

!==========================================================
if (.not. present(cmessage)) then
  errstr = ''
else
  errstr = newline//trim(adjustl(cmessage))
end if

if (def_src == '') then
  loc_src = ''
  srcname = modname//':'//routinename
else
  loc_src = 'Problem performing defensive checks ('//                          &
            trim(adjustl(def_src))//')'//newline
  srcname = modname//':'//routinename//' called by '//trim(adjustl(def_src))
end if

loc_str  = string
pass_chk   = .false.
icode      = 0

result_chk = 'FAIL'

! Strip off brackets
index_val = index(loc_str, '[')
loc_str(index_val:index_val) = ''
index_val = index(loc_str, ']')
loc_str(index_val:index_val) = ''

loc_str  = trim(adjustl(loc_str))
loc_name = trim(adjustl(var_name))

! Check for reals in the condition string
if ( index(loc_str,'.') /= 0) then
  write(errstr,'(A,I0,A)')                                            newline//&
    trim(loc_src)//'Attempting to test an integer variable, '//trim(loc_name)//&
    ', (',var,')'//newline//'against a condition involving a non-integer, (' //&
    trim(loc_str) //                                                           &
    ')'
  icode = 10
  call ereport(trim(srcname), icode, errstr)
end if

! Separate the condition string into separate tests
call split_conditions()

! Check that there are some valid tests
if (nconditions == 0) then
  write(errstr,'(A,I0,A)')                                            newline//&
    trim(loc_src)//'No valid checks, '//string//                               &
    ', provided for integer variable, '//trim(loc_name)//', (',var,')'
  icode = 30
  call ereport(trim(srcname), icode, errstr)
end if


! Loop over the conditions, converting them to integers so they can be tested
do i=1, nconditions

  operator_string  = operators(i)

  if (trim(operator_string) == ':') then

    ! For range tests which have upper and lower bounds
    read(unit=lbounds(i),fmt=*, iostat=icode) lowerbound
    if (icode /= 0) then
      icode = 40
      write(errstr,'(A,I0,A)')                                        newline//&
        trim(loc_src)//'Cannot convert lower bound string "'//                 &
        trim(lbounds(i))//'" into an integer when'//                  newline//&
        'running defensive checks on '//trim(loc_name)//', (',var,').'
      call ereport(trim(srcname), icode, errstr)
    end if

    read(unit=ubounds(i),fmt=*, iostat=icode) upperbound
    if (icode /= 0) then
      icode = 50
      write(errstr,'(A,I0,A)')                                        newline//&
        trim(loc_src)//'Cannot convert upper bound string "'//                 &
        trim(ubounds(i))//'" into an integer when'//                  newline//&
        'running defensive checks on '//trim(loc_name)//', (',var,').'
      call ereport(trim(srcname), icode, errstr)
    end if

  else

    read(unit=values(i),fmt=*, iostat=icode) condition
    if (icode /= 0) then
      icode = 60
      write(errstr,'(A,I0,A)')                                        newline//&
        trim(loc_src)//'Cannot convert string "'//                             &
        trim(values(i))//'" into an integer when running '//          newline//&
        'defensive checks on '//trim(loc_name)//', (',var,').'
      call ereport(trim(srcname), icode, errstr)
    end if

  end if


  ! Now do the test
  select case (trim(operator_string))
  case ('<')
    if ( var < condition )  pass_chk = .true.
  case ('>')
    if ( var > condition )  pass_chk = .true.
  case ('>=')
    if ( var >= condition ) pass_chk = .true.
  case ('<=')
    if ( var <= condition ) pass_chk = .true.
  case (':')
    if ( var >= lowerbound .and. var <= upperbound ) pass_chk = .true.
  case ('==')
    if ( var == condition ) pass_chk = .true.

  case DEFAULT
    write(errstr,'(A)')                                               newline//&
      trim(loc_src)//'Unrecognised operator, ('//trim(operator_string)//       &
      '), when checking'//newline//'integer variable, '//trim(loc_name)//      &
      ', (',var,')'
    icode = 70
    call ereport(trim(srcname), icode, errstr)
    exit
  end select

  if (pass_chk) exit

end do

if (present(report_pass)) then
  report_pass = pass_chk
  if (pass_chk) result_chk = 'pass'
  write(ummessage,'(A,I0,A)')                                                  &
    'Defensive checks on '//trim(loc_name)//' (',var,') against '//            &
    trim(adjustl(string))//' .... '//trim(result_chk)
  call umprint(ummessage, src=trim(srcname))

else if (.not. pass_chk) then

  write(errstr,'(A,I0,A)')                                            newline//&
    trim(loc_name)//', (',var,'), has failed to pass any defensive checks, ' //&
    newline//trim(adjustl(string))//'.'//trim(errstr)
  icode = 80
  call ereport(trim(srcname), icode, errstr)

end if

return
end subroutine chk_int_str



!==============================================================================
!==============================================================================



subroutine chk_real_str (var, var_name, string, report_pass, cmessage)

! Checks a real variable against a list of criteria
! passed as a string in the form '[Check_A, Check_B, ... ]'.
!
! Recognised checks
!  * ==Real (This is assumed if no operator is provided)
!  * <Real
!  * >Real
!  * >=Real
!  * <=Real
!  * Real lower bound : Real upper bound
!
! As this set of checks should be for valid values of the variable
! under test, the defensive check is passed if any of the criteria
! are met. If none of the specified checks are met, ereport will
! called to abort.
!
! Two optional arguments are provided:
! report_pass - (logical) value indicates if any of the checks passed
!               Ereport is not called with this option.
! cmessage    - If Chk_var calls ereport, the variable name and
!               value will be included, along with some information.
!               if a cmessage argument is provided, the message
!               will be appended to chk_var's ereport message.

implicit none

real,             intent(in) :: var        ! The real variable to be checked
character(len=*), intent(in) :: var_name   ! The name of the test variable
character(len=*), intent(in) :: string     ! String containg test criteria

! Return flag for pass/failure of test criteria. If this given in argument
! list then ereport will not be called
logical,          intent(in out), optional :: report_pass

! Additonal comments to af to ereport in event of failure of tests
character(len=*), intent(in), optional :: cmessage

! Local variables

character(len=*), parameter :: routinename ='CHK_REAL_STR'

integer :: i, index_val, icode
real :: condition
real :: lowerbound
real :: upperbound

!==========================================================
if (.not. present(cmessage)) then
  errstr = ''
else
  errstr = newline//trim(adjustl(cmessage))
end if

if (def_src =='') then
  loc_src = ''
  srcname = modname//':'//routinename
else
  loc_src = 'Problem performing defensive checks ('//                          &
            trim(adjustl(def_src))//')'//newline
  srcname = modname//':'//routinename//' called by '//trim(adjustl(def_src))
end if

loc_str    = string
pass_chk   = .false.
icode      = 0
result_chk = 'FAIL'

! Strip off brackets
index_val = index(loc_str, '[')
loc_str(index_val:index_val) = ''
index_val = index(loc_str, ']')
loc_str(index_val:index_val) = ''

write(tmp_str,'(E10.3)') var

loc_str  = adjustl(loc_str)
loc_name = adjustl(var_name)
tmp_str  = adjustl(tmp_str)

call split_conditions ()

if (nconditions == 0) then
  write(errstr,'(A,I0,A)')                                            newline//&
    trim(loc_src)//'No valid checks, '//trim(adjustl(string))//                &
    ', provided for real variable, '//trim(loc_name)//', (',var,')'
  icode = 20
  call ereport(trim(srcname), icode, errstr)
end if


do i=1, nconditions

  operator_string = operators(i)

  if (trim(adjustl(operator_string)) == ':') then
    if ( index(lbounds(i),'.') == 0) then

      write(errstr, '(A)')                                            newline//&
        trim(loc_src)//'You are attempting to test a real variable, '//        &
        trim(loc_name)//'('//trim(tmp_str)//'), '//                   newline//&
        'with a non-real lower bound, ('//trim(lbounds(i))//')'

      icode = 30
      call ereport(trim(srcname), icode, errstr)
    end if

    read(unit=lbounds(i),fmt=*, iostat=icode) lowerbound
    if (icode /= 0) then
      icode = 40
      write(errstr,'(A)')                                             newline//&
        trim(loc_src)//'Cannot convert lower bound string "'//                 &
        trim(lbounds(i))//'" into a real when'//                      newline//&
        'running defensive checks on '//trim(loc_name)//', ('//                &
        trim(tmp_str)//').'
      call ereport(trim(srcname), icode, errstr)
    end if


    if ( index(ubounds(i),'.') == 0) then
      write(errstr,'(A)')                                             newline//&
        trim(loc_src)//'You are attempting to test a real variable, '//        &
        trim(loc_name)//'('//trim(tmp_str)//'), '//                   newline//&
        'with a non-real upper bound, ('//trim(ubounds(i))//')'

      icode = 50
      call ereport(trim(srcname), icode, errstr)
    end if

    read(unit=ubounds(i),fmt=*, iostat=icode) upperbound
    if (icode /= 0) then
      icode = 60
      write(errstr,'(A)')                                             newline//&
        trim(loc_src)//'Cannot convert upper bound string "'//                 &
        trim(ubounds(i))//'" into a real when'//                      newline//&
        'running defensive checks on '//trim(loc_name)//', ('//                &
        trim(tmp_str)//').'
      call ereport(trim(srcname), icode, errstr)
    end if

  else
    if ( index(values(i),'.') == 0) then
      write(errstr,'(A)')                                             newline//&
        trim(loc_src)//'You are attempting to test real variable, '//          &
        trim(loc_name)//', ('//trim(tmp_str)//'), '//                 newline//&
        'with a non-real, ('//trim(values(i))//')'
      icode = 70
      call ereport(trim(srcname), icode, errstr)
    end if

    read(unit=values(i),fmt=*, iostat=icode) condition
    if (icode /= 0) then
      icode = 80
      write(errstr,'(A)')                                             newline//&
        trim(loc_src)//'Cannot convert string "'//                             &
        trim(values(i))//'" into a real when running '//              newline//&
        'defensive checks on '//trim(adjustl(var_name))//', ('//               &
        trim(tmp_str)//').'
      call ereport(trim(srcname), icode, errstr)
    end if

  end if


  select case (trim(operator_string))
  case ('<')
    if ( var <  condition ) pass_chk = .true.
  case ('>')
    if ( var >  condition ) pass_chk = .true.
  case ('>=')
    if ( var >= condition ) pass_chk = .true.
  case ('<=')
    if ( var <= condition ) pass_chk = .true.
  case (':')
    if ( var >= lowerbound .and. var <= upperbound ) pass_chk = .true.
  case ('==')
    if ( var == condition ) pass_chk = .true.

  case DEFAULT
    write(errstr,'(A)')                                                        &
      trim(loc_src)//'Unrecognised operator, ('//trim(operator_string)//')'
    icode = 90
    call ereport(trim(srcname), icode, errstr)
  end select

  if (pass_chk) exit

end do

if (present(report_pass)) then
  report_pass = pass_chk

  if (pass_chk) result_chk = 'pass'
  write(ummessage,'(A)')                                                       &
    'Defensive checks on '//trim(loc_name)//' ('//trim(tmp_str)//') against '//&
    trim(adjustl(string))//' .... '//trim(result_chk)
  call umprint(ummessage,src=trim(srcname))

else if (.not. pass_chk) then
  write(errstr,'(A)')                                                 newline//&
    trim(loc_src)//trim(loc_name)//', ('//trim(tmp_str)//                      &
    '), has failed to pass any '//newline//'defensive checks, '//trim(string)//&
    '.'//trim(errstr)
  icode = 100
  call ereport(trim(srcname), icode, errstr)
end if

return
end subroutine chk_real_str



!==============================================================================
!==============================================================================



subroutine chk_int_arr (var, var_name, intarr, report_pass, cmessage)

! Checks an integer variable against a list of criteria
! passed as an integer array in the form '[Check_A, Check_B, ... ]'.
!
! This routine only checks if the integer variable exists in the
! the array of valid values
!
! As this set of checks should be for valid values of the variable
! under test, the defensive check is passed if any of the criteria
! are met. If none of the specified checks are met, ereport will
! called to abort.
!
! Two optional arguments are provided:
! report_pass - (logical) value indicates if any of the checks passed
!               Ereport is not called with this option.
! cmessage    - If chk_var calls ereport, the variable name and
!               value will be included, along with some information.
!               if a cmessage argument is provided, the message
!               will be appended to chk_var's ereport message.

implicit none

integer,          intent(in) :: var       ! The integer variable to be checked
character(len=*), intent(in) :: var_name  ! The name of the test variable
integer,          intent(in) :: IntArr(:) ! Valid integer values

! Return flag for pass/failure of test criteria. If this given in argument
! list then ereport will not be called
logical,          intent(in out), optional :: report_pass

! Additonal comments to af to ereport in event of failure of tests
character(len=*), intent(in),    optional :: cmessage


! Local variables
character(len=*), parameter :: routinename ='CHK_INT_ARR'
character(len=str_len) :: StrIntArr ! Arbitrarily sized local temporary string.
character(len=str_len) :: StrTest ! Arbitrarily sized local test string.

logical :: pass_chk
integer :: i
integer :: icode
integer :: ioerr

!==========================================================
if (.not. present(cmessage)) then
  errstr = ''
else
  errstr = newline//trim(adjustl(cmessage))
end if

if ( def_src == '' ) then
  loc_src = ''
  srcname = modname//':'//routinename
else
  loc_src = 'Problem performing defensive checks ('//                          &
            trim(adjustl(def_src))//')'//newline
  srcname = modname//':'//routinename//' called by '//trim(adjustl(def_src))
end if

result_chk = 'FAIL'
pass_chk   = .false.
loc_name   = adjustl(var_name)

if (any(var == IntArr(:))) pass_chk = .true.

! Create a list of our array elements in the form of a single string.
! Clearly we must have at least one element.
! Append any subsequent elements, if there are any, separated by commas.
StrIntArr=''
do i = 1, size(IntArr)

  ! Make sure we don't try to append more data than we have space for!
  ! If we don't have enough space, add a message to that effect and exit.
  write(StrTest,'(I0)',iostat=ioerr) IntArr(i)
  if (((len_trim(StrIntArr) + len_trim(StrTest) + 1) > str_len) .or.           &
        ioerr /= 0) then
    StrIntArr(str_len-19:str_len) = " (details truncated)"
    exit
  end if
  if (i > 1) then
    StrIntArr=trim(StrIntArr)//','//trim(StrTest)
  else
    StrIntArr=trim(StrTest)
  end if

end do

if (present(report_pass)) then

  report_pass = pass_chk
  if (pass_chk) result_chk = 'pass'

  write(ummessage,'(A,I0,A)')                                                  &
    'Defensive checks on '//trim(loc_name)//' (',var,                          &
    '), against ['//trim(StrIntArr)//'] .... '//trim(result_chk)
  call umprint(ummessage,src=trim(srcname))

else if (.not. pass_chk) then

  write(errstr,'(A,I0,A)')                                            newline//&
    trim(loc_src)//trim(loc_name)//',(',var,                                   &
    '), has failed defensive checks, ['//trim(StrIntArr)//'].'//trim(errstr)
  icode = 10
  call ereport(trim(srcname), icode, errstr)

end if

return
end subroutine chk_int_arr



!==============================================================================
!==============================================================================



subroutine split_conditions ()

! Support routine for chk_var, the provided string is separated and
! used to populate the module arrays:
!
! * operators
! * values
! * lbounds
! * ubounds
!
! String inputs are separated based on "," with the following recognised
! operators
!  * ==
!  * <
!  * >
!  * >=
!  * <=
!  * :
! If none of these are found, equals (==) is assumed


implicit none

character(len=condition_len) :: condition_str

integer :: i, index_val
logical :: last_condition

! Take a local copy of the input string and
! initialise checking arrays
operators(:) = ''
values(1:max_n_conditions) = '                    '
lbounds(1:max_n_conditions) = '                    '
ubounds(1:max_n_conditions) = '                    '
nconditions  = 0

if (len(trim(loc_str)) == 0) return

! 1) Grab the first character,
!    to see if it is > or <, if so check for an >= or <=.
! 2) Check for range character ":"
! 3) Assume the check is for equals, i.e. "=="
do i=1, max_n_conditions

  index_val = index(loc_str, ',')

  ! Check if this is the last condition, i.e. there
  ! should be no more "," to separate the checks.
  last_condition = (index_val == 0)

  if (last_condition) then
    condition_str  = trim(adjustl(loc_str))
  else
    condition_str  = trim(adjustl(loc_str(1:index_val-1)))

    loc_str(1:index_val) = ''
    loc_str = trim(adjustl(loc_str))
  end if

  if ( condition_str(1:1) == '<' .or. condition_str(1:1) == '>' ) then

    if (condition_str(2:2) == '=') then
      operators(i) = condition_str(1:2)
      values(i)    = trim(adjustl(condition_str(3:)))
    else
      operators(i) = condition_str(1:1)
      values(i)    = trim(adjustl(condition_str(2:)))
    end if

  else if ( index(condition_str,':') /= 0 ) then
    operators(i) = ':'
    index_val  = index(condition_str, ':')
    lbounds(i) = trim(adjustl(condition_str(:index_val-1)))
    ubounds(i) = trim(adjustl(condition_str(index_val+1:)))
  else
    operators(i) = '=='
    if ( index(condition_str,'==') /= 0 ) then
      values(i) = trim(adjustl(condition_str(3:)))
    else
      values(i) = trim(adjustl(condition_str))
    end if
  end if

  nconditions = nconditions + 1
  if (last_condition) exit
end do

return
end subroutine split_conditions


end module chk_opts_mod
