! *********************************COPYRIGHT************************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *********************************COPYRIGHT************************************
! Some of the content of this file has been produced with the assistance of
! Anthropic Claude Opus 5 (Claude Code).
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: atmos_service_vera
!
! Part of the Vera scheme for diagnosing visual range.
!
! This module contains the Vera precision definitions.
!
! There is no executable code here, just some definitions for
! data types.
!

module vera_kind_mod

! grab the real_umphys data type
use um_types, only: real_umphys

implicit none

! Description:
!   This module contains the Vera precision definitions.
!
! Method:
!   There is no executable code here, just some definitions for
!   data types - integer and real number precisions and logicals.
!
!   For more detail, please refer to the Vera user guide.
!
! Code description:
!   Language: Fortran 2003
!   This code is written to UMDP3 standards.

character (len=*), parameter, private :: ModuleName='VERA_KIND_MOD'

private

!=============================================================================
! some real, integer and logical definitions
!=============================================================================

! single precision real numbers
integer, parameter, public :: vera_sp  = selected_real_kind( 6, 37 )

! double precision real numbers
integer, parameter, public :: vera_dp  = selected_real_kind( 15, 307 )

! quadruple precision real numbers
integer, parameter, public :: vera_qp  = selected_real_kind( 30, 291 )

! 8 bit integer numbers have range
! [-128, 127]
integer, parameter, public ::  vera_i8 = selected_int_kind( 2 )

! 16 bit integer numbers have range
! [-32768, 32767]
integer, parameter, public :: vera_i16 = selected_int_kind( 4 )

! 32 bit integer numbers have range
! [-2147483648, 2147483647]
integer, parameter, public :: vera_i32 = selected_int_kind( 8 )

! 64 bit integer numbers have range
! [-9223372036854775808, -9223372036854775807]
integer, parameter, public :: vera_i64 = selected_int_kind( 16 )

! define a logical type
logical (kind=4), public   :: vera_logical_flag4

!=============================================================================
! UM real, integer and logical definitions - use the compiler defaults
!=============================================================================

! define the real precision to use for calling UM subroutines,
! e.g. qsat_wat
integer, parameter, public :: um_r = real_umphys

! define the integer precision to use for calling UM subroutines,
! e.g. qsat_wat
integer            :: um_integer_example
integer, parameter, public :: um_i = kind( um_integer_example )

! define the default logical type used in the UM
logical            :: um_logical_example
integer, parameter, public :: um_l = kind( um_logical_example )

!=============================================================================
! MINPACK default real, integer and logical definitions
!=============================================================================

! define the real precision to use for calling MINPACK subroutines
integer, parameter, public :: minpack_r = um_r

! define the integer precision to use for calling MINPACK subroutines
integer, parameter, public :: minpack_i = um_i

! define the logical type to use in MINPACK subroutines
integer, parameter, public :: minpack_l = um_l

!=============================================================================
! Vera default real, integer and logical definitions
!=============================================================================

! define the default real precision to use
integer, parameter, public :: vera_real    = um_r

! define the default integer range to use
integer, parameter, public :: vera_integer = um_i

! define the logical type to use
integer, parameter, public :: vera_logical = um_l

end module vera_kind_mod
