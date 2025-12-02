! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module cloudfracs_type_mod

use comorph_constants_mod, only: real_hmprec, name_length

implicit none

save

! Flag indicating whether the variables contained in this
! module have been set
logical :: l_init_cloudfracs_type_mod = .false.


! Type definition to store sub-grid fractional areas for
! cloud and precipitation.  Note that If the cloud fractions
! are prognostic, they are passed in through the fields_type
! structure instead of here.  But if not using prognostic
! cloud fractions, the diagnostic fractions are passed
! in this structure.  CoMorph doesn't yet include the precip
! fraction (rain / graupel fraction) as a primary prognostic
! field, so frac_precip is only passed in through here.
! This type is also used to pass the diagnosed convective cloud
! fraction and properties in and out of CoMorph.  Values passed
! in are valid at start-of-timestep, values passed out are the
! updated end-of-timestep convective cloud fields.
type :: cloudfracs_type

  ! Liquid cloud fraction
  real(kind=real_hmprec), pointer :: frac_liq(:,:,:) => null()

  ! Ice cloud fraction
  real(kind=real_hmprec), pointer :: frac_ice(:,:,:) => null()

  ! Bulk cloud fraction
  real(kind=real_hmprec), pointer :: frac_bulk(:,:,:) => null()

  ! Rain / graupel fraction
  real(kind=real_hmprec), pointer :: frac_precip(:,:,:) => null()

  ! Convective liquid, ice and bulk cloud fraction
  real(kind=real_hmprec), pointer :: frac_liq_conv(:,:,:) => null()
  real(kind=real_hmprec), pointer :: frac_ice_conv(:,:,:) => null()
  real(kind=real_hmprec), pointer :: frac_bulk_conv(:,:,:) => null()

  ! Convective cloud grid-mean liquid and ice mixing-ratios
  real(kind=real_hmprec), pointer :: q_cl_conv(:,:,:) => null()
  real(kind=real_hmprec), pointer :: q_cf_conv(:,:,:) => null()

end type cloudfracs_type


! Addresses of the cloud-fraction fields in their own
! super-array (used when prognostic cloud is off)
integer, parameter :: n_cloudfracs = 4
integer, parameter :: i_frac_liq = 1
integer, parameter :: i_frac_ice = 2
integer, parameter :: i_frac_bulk = 3
integer, parameter :: i_frac_precip = 4

! Addresses of the convective cloud fields in their own
! super-array (these depend on whether using convective cloud)
integer :: n_convcloud = 0
integer :: i_frac_liq_conv = 0
integer :: i_frac_ice_conv = 0
integer :: i_frac_bulk_conv = 0
integer :: i_q_cl_conv = 0
integer :: i_q_cf_conv = 0

! List of the names of each field, stored as character strings
! (used for labelling diagnostics and error print-outs)
character(len=name_length), allocatable :: convcloud_names(:)


contains


! Subroutine to set the addresses of convective cloud fields,
! depending on switches
subroutine cloudfracs_set_addresses()

use comorph_constants_mod, only: i_convcloud, i_convcloud_liqonly,             &
                                 i_convcloud_mph

implicit none

! Set flag to indicate that we don't need to call this routine
! again!
l_init_cloudfracs_type_mod = .true.

! Set addresses for fields in the convective cloud super-array,
! and set strings containing the name of each field in the list
select case( i_convcloud )
case ( i_convcloud_liqonly )
  ! Only using liquid convective cloud
  ! (but we also still output the bulk convective cloud fraction,
  !  which includes ice cloud, for the purpose of calculating
  !  the convective cloud base and top).
  n_convcloud = 3
  i_frac_liq_conv = 1
  i_frac_bulk_conv = 2
  i_q_cl_conv = 3
  allocate( convcloud_names(n_convcloud) )
  convcloud_names(i_frac_liq_conv)  = "frac_liq_conv"
  convcloud_names(i_frac_bulk_conv) = "frac_bulk_conv"
  convcloud_names(i_q_cl_conv)      = "q_cl_conv"
case ( i_convcloud_mph )
  ! Using separate liquid and ice cloud, with variable overlap
  n_convcloud = 5
  i_frac_liq_conv = 1
  i_frac_ice_conv = 2
  i_frac_bulk_conv = 3
  i_q_cl_conv = 4
  i_q_cf_conv = 5
  allocate( convcloud_names(n_convcloud) )
  convcloud_names(i_frac_liq_conv)  = "frac_liq_conv"
  convcloud_names(i_frac_ice_conv)  = "frac_ice_conv"
  convcloud_names(i_frac_bulk_conv) = "frac_bulk_conv"
  convcloud_names(i_q_cl_conv)      = "q_cl_conv"
  convcloud_names(i_q_cf_conv)      = "q_cf_conv"
end select

return
end subroutine cloudfracs_set_addresses


! Subroutine to nullify the pointers when finished
subroutine cloudfracs_nullify( cloudfracs )
implicit none
type(cloudfracs_type), intent(in out) :: cloudfracs

cloudfracs % frac_liq    => null()
cloudfracs % frac_ice    => null()
cloudfracs % frac_bulk   => null()
cloudfracs % frac_precip => null()
cloudfracs % frac_liq_conv  => null()
cloudfracs % frac_ice_conv  => null()
cloudfracs % frac_bulk_conv => null()
cloudfracs % q_cl_conv => null()
cloudfracs % q_cf_conv => null()

return
end subroutine cloudfracs_nullify


! Subroutine to check the above fields for bad values
subroutine cloudfracs_check_bad_values( cloudfracs,                            &
                                        where_string )

use comorph_constants_mod, only: name_length, l_cv_cloudfrac,                  &
                     i_convcloud, i_convcloud_none,                            &
                     i_convcloud_mph
use check_bad_values_mod, only: check_bad_values_3d

implicit none

! Structure containing pointers to cloudfrac fields to check
type(cloudfracs_type), intent(in) :: cloudfracs

! Character string describing where in the convection scheme
! we are, for constructing error message if bad value found.
character(len=name_length), intent(in) :: where_string

! Character strong to store field names
character(len=name_length) :: field_name

! Lower and upper bounds of array
integer :: lb(3), ub(3)

! Flag passed into check_bad_values;
! all fields checked are positive-only so hardwire to true
logical, parameter :: l_positive = .true.


if ( .not. l_cv_cloudfrac ) then
  ! If CoMorph not treating the cloud fractions as prognostics,
  ! then diagnostic cloud fractions are included in here...

  field_name = "frac_liq"
  lb = lbound( cloudfracs % frac_liq )
  ub = ubound( cloudfracs % frac_liq )
  call check_bad_values_3d( lb, ub, cloudfracs % frac_liq,                     &
                            where_string, field_name,                          &
                            l_positive )

  field_name = "frac_ice"
  lb = lbound( cloudfracs % frac_ice )
  ub = ubound( cloudfracs % frac_ice )
  call check_bad_values_3d( lb, ub, cloudfracs % frac_ice,                     &
                            where_string, field_name,                          &
                            l_positive )

  field_name = "frac_bulk"
  lb = lbound( cloudfracs % frac_bulk )
  ub = ubound( cloudfracs % frac_bulk )
  call check_bad_values_3d( lb, ub, cloudfracs % frac_bulk,                    &
                            where_string, field_name,                          &
                            l_positive )

end if

! Check the precip fraction
field_name = "frac_precip"
lb = lbound( cloudfracs % frac_precip )
ub = ubound( cloudfracs % frac_precip )
call check_bad_values_3d( lb, ub, cloudfracs % frac_precip,                    &
                          where_string, field_name,                            &
                          l_positive )

if ( i_convcloud > i_convcloud_none ) then
  ! Check convective cloud fields if they are used

  ! Liquid convective cloud in use...

  field_name = "frac_liq_conv"
  lb = lbound( cloudfracs % frac_liq_conv )
  ub = ubound( cloudfracs % frac_liq_conv )
  call check_bad_values_3d( lb, ub, cloudfracs % frac_liq_conv,                &
                            where_string, field_name,                          &
                            l_positive )

  field_name = "frac_bulk_conv"
  lb = lbound( cloudfracs % frac_bulk_conv )
  ub = ubound( cloudfracs % frac_bulk_conv )
  call check_bad_values_3d( lb, ub, cloudfracs % frac_bulk_conv,               &
                            where_string, field_name,                          &
                            l_positive )

  field_name = "q_cl_conv"
  lb = lbound( cloudfracs % q_cl_conv )
  ub = ubound( cloudfracs % q_cl_conv )
  call check_bad_values_3d( lb, ub, cloudfracs % q_cl_conv,                    &
                            where_string, field_name,                          &
                            l_positive )

  if ( i_convcloud == i_convcloud_mph ) then
    ! Ice and mixed-phase convective cloud in use as well...

    field_name = "frac_ice_conv"
    lb = lbound( cloudfracs % frac_ice_conv )
    ub = ubound( cloudfracs % frac_ice_conv )
    call check_bad_values_3d( lb, ub, cloudfracs % frac_ice_conv,              &
                              where_string, field_name,                        &
                              l_positive )

    field_name = "q_cf_conv"
    lb = lbound( cloudfracs % q_cf_conv )
    ub = ubound( cloudfracs % q_cf_conv )
    call check_bad_values_3d( lb, ub, cloudfracs % q_cf_conv,                  &
                              where_string, field_name,                        &
                              l_positive )

  end if

end if  ! ( i_convcloud > i_convcloud_none )


return
end subroutine cloudfracs_check_bad_values


end module cloudfracs_type_mod
