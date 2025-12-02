! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! This module has an underlying dependency on LFRic's logging module

module um2lfric_read_um_file_mod

! Intrinsic modules
use, intrinsic :: iso_c_binding, only: c_bool

! Shumlib modules
use f_shum_file_mod, only: shum_file_type
use f_shum_field_mod, only: shum_field_type

! UM2LFRic modules
use lfricinp_check_shumlib_status_mod, only: shumlib
use lfricinp_um_parameters_mod, only: fnamelen, um_integer64

! LFRic modules
use log_mod,     only: log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR,             &
                       log_scratch_space

implicit none

private

public :: um2lfric_read_um_file, um2lfric_close_um_file, um_input_file

type(shum_file_type), save :: um_input_file

contains

! DEPENDS ON: c_shum_byteswap.o
! This is required to force fcm-make to compile the C code; whilst the built-in
! dependency analyser successfully works out that it needs to compile the
! Fortran side of the byte-swapping code, it requires an explicit statement
! to force it to compile the C part of the byte-swapping code. This is
! currently the approved way of linking Fortran and C in fcm-make.

!-------------------------------------------------------------------------------

subroutine um2lfric_read_um_file(fname)

  implicit none

  character(len=fnamelen), intent(in) :: fname
  character(len=*), parameter :: routinename='um2lfric_read_um_file'

  ! Load the UM file
  call log_event('Loading file '//trim(fname), LOG_LEVEL_INFO)
  call shumlib(routinename//'::open_file', um_input_file%open_file(fname),&
                print_on_success=.true._C_BOOL)

  call shumlib(routinename//'::open_file', um_input_file%read_header(),  &
                print_on_success=.true._C_BOOL)

end subroutine um2lfric_read_um_file

subroutine um2lfric_close_um_file()
  implicit none
  character(len=*), parameter :: routinename='um2lfric_close_um_file'

  ! Close file
  call shumlib(routinename//'::close_file', um_input_file%close_file() )

end subroutine um2lfric_close_um_file

end module um2lfric_read_um_file_mod
