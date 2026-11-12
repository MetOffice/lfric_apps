! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!
!  Reads in a UKCA_RADAER look-up table (namelist format) and copy
!  its contents in the structure given in argument.
!
!
! Subroutine Interface:
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: UKCA_UM
!
module ukca_radaer_lut_read_in

use umprintmgr,             only: newline
use errormessagelength_mod, only: errormessagelength
use ereport_mod,            only: ereport
use parkind1,               only: jpim, jprb
use yomhook,                only: lhook, dr_hook

implicit none

!
! Local copy of the contents of a UKCA look-up table.
! Here, fixed array sizes are needed.
! Note that array indexing starts at 0 in the namelist.
integer, parameter :: npd_x  = 50
integer, parameter :: npd_nr = 50
integer, parameter :: npd_ni = 50

real :: stdev
integer :: n_x
integer :: n_nr
integer :: n_ni
integer :: funit
real :: x_min
real :: x_max
real :: nr_min
real :: nr_max
real :: ni_min
real :: ni_max
real :: ni_c
real :: ukca_absorption(0:npd_x,0:npd_ni,0:npd_nr)
real :: ukca_scattering(0:npd_x,0:npd_ni,0:npd_nr)
real :: ukca_asymmetry( 0:npd_x,0:npd_ni,0:npd_nr)
real :: volume_fraction(0:npd_x)

namelist /ukcanml/                                                             &
          stdev,                                                               &
          n_x, n_nr, n_ni,                                                     &
          x_min, x_max,                                                        &
          nr_min, nr_max,                                                      &
          ni_min, ni_max, ni_c,                                                &
          ukca_absorption, ukca_scattering, ukca_asymmetry,                    &
          volume_fraction

character(len=*), parameter, private :: ModuleName='UKCA_RADAER_LUT_READ_IN'

! Overload interface
interface ukca_radaer_get_lut_index
  module procedure ukca_radaer_get_lut_index_scalar
  module procedure ukca_radaer_get_lut_index_array
end interface

contains


subroutine ukca_radaer_get_lut_index_scalar(Nni, im_m, ni_min, ni_max, ni_c,   &
                                            ni_ind, ni_c_power)
implicit none

integer, intent(in)  :: Nni       ! number of ni values in LUT
real,    intent(in)  :: im_m      ! calculated value of ni
real,    intent(in)  :: ni_min    ! minimum value of ni in the LUT
real,    intent(in)  :: ni_max    ! maximum value of ni in the LUT
real,    intent(in)  :: ni_c   ! parameter required to calculate ni increments
integer, intent(out) :: ni_ind    ! nearest-neighbour index for ni
real,    intent(in), optional :: ni_c_power

real    :: im_m_array(1)
integer :: ni_ind_array(1)

im_m_array(1) = im_m
call ukca_radaer_get_lut_index_array(Nni, im_m_array, ni_min, ni_max, ni_c,    &
                                     ni_ind_array, 1, ni_c_power=ni_c_power)
ni_ind = ni_ind_array(1)

return
end subroutine ukca_radaer_get_lut_index_scalar

subroutine ukca_radaer_get_lut_index_array(Nni, im_m, ni_min, ni_max, ni_c,    &
                                           ni_ind, length, ni_c_power)
use vectlib_mod, only: log_v
implicit none

integer, intent(in)  :: length           ! Array length
integer, intent(in)  :: Nni              ! number of ni values in LUT
real,    intent(in)  :: im_m(1:length)   ! calculated value of ni
real,    intent(in)  :: ni_min           ! minimum value of ni in the LUT
real,    intent(in)  :: ni_max           ! maximum value of ni in the LUT
real,    intent(in)  :: ni_c    ! parameter required to calculate ni increments
integer, intent(out) :: ni_ind(1:length) ! nearest-neighbour index for ni
real,    intent(in), optional :: ni_c_power

! local variables
real :: a, b, incr_ni
real :: ni_c_power_local
real, parameter :: min_ni_c = 0.001 ! Lowest value of ni_c to accept
real, parameter :: max_ni_c = 5.0   ! Highest value of ni_c to accept

real, parameter :: inv_ln_10 = 1.0 / log(10.0)

real :: logs_array_in(1:length)
real :: logs_array_out(1:length)

integer :: icode
character(len=errormessagelength) :: cmessage

!DrHook-related parameters
! note dr_hook commented out here until overhead can be reduced, see UM:#2787
!integer(kind=jpim), parameter :: zhook_in  = 0
!integer(kind=jpim), parameter :: zhook_out = 1
!real(kind=jprb)               :: zhook_handle
character(len=*), parameter :: RoutineName='UKCA_RADAER_GET_LUT_INDEX_ARRAY'

!if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_in, zhook_handle)

icode = 0

! If ni_c > max_ni_c then produce a fatal error
!
if (ni_c > max_ni_c) then
  icode = 1
  cmessage='UKCA RADAER Look-up table'//                   newline//           &
           'NI_C exceeds upper limit'

  !if (lhook) call dr_hook(ModuleName//':'//RoutineName, zhook_out,zhook_handle)
  return
end if

! As ni_c tends towards zero the function below tends towards a linear scale
! With ni_c = 0, a and b become infinity. So we will protect this
! from happening by using a linear scale if ni_c < min_ni_c
!
if (ni_c > min_ni_c) then

  if (present(ni_c_power)) then
    ni_c_power_local = ni_c_power
  else
    ni_c_power_local = 10.0**ni_c
  end if

  a = ni_max / ((ni_c_power_local) - 1.0)
  b = real(Nni) / ni_c

  logs_array_in(:) = (iM_m(:)/a) + 1.0

  call log_v(length, logs_array_in, logs_array_out)

  ni_ind(:) = nint(b * logs_array_out(:) * inv_ln_10 ) + 1

else

  incr_ni   = (ni_max - ni_min) / real(Nni-1)
  ni_ind(:) = nint( (im_m(:) - ni_min) / incr_ni ) + 1

end if

ni_ind(:) = min(Nni, max(1, ni_ind(:)))

if (icode /= 0) then
  call Ereport(RoutineName,icode,cmessage)
end if

!if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine ukca_radaer_get_lut_index_array

end module ukca_radaer_lut_read_in
