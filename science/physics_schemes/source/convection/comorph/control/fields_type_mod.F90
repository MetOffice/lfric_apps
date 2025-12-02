! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module fields_type_mod

use comorph_constants_mod, only: real_cvprec, real_hmprec, name_length

implicit none

save

! Flag indicating whether the variables contained in this
! module have been set
logical :: l_init_fields_type_mod = .false.


!----------------------------------------------------------------
! Type definition to store pointers to the fields in a list
!----------------------------------------------------------------
! This can be looped over, to allow shorter code when doing
! exactly the same thing to all the fields
type :: fields_list_type
  ! Points at full fields which are in/out from CoMorph,
  ! so declared with the host-model's precision.
  real(kind=real_hmprec), pointer :: pt(:,:,:) => null()
end type fields_list_type


!----------------------------------------------------------------
! Type definition which groups together the full 3-D fields
!----------------------------------------------------------------
! This just contains pointers to all the primary model fields
! and tracers updated by convection;
! different instances of this structure can be used for
! start-of-timestep, end-of-timestep, increments etc.
! Since this structure is used entirely for variables
! which are passed in and out of comorph to the host model,
! they are declared with the host-model's precision,
! not the convection scheme's precision

! If adding a new primary field, don't forget to also
! add it to fields_diags_type in fields_diags_type_mod
! so that it can be output as a diagnostic.
! You also need to take care to add it to all the subroutines
! below which act upon the primary fields
! (fields_set_addresses, fields_nullify, etc).
type :: fields_type

  ! Wind components
  real(kind=real_hmprec), pointer :: wind_u(:,:,:) => null()
  real(kind=real_hmprec), pointer :: wind_v(:,:,:) => null()
  real(kind=real_hmprec), pointer :: wind_w(:,:,:) => null()

  ! Temperature
  real(kind=real_hmprec), pointer :: temperature(:,:,:) => null()

  ! Water vapour mixing ratio
  real(kind=real_hmprec), pointer :: q_vap(:,:,:) => null()

  ! Hydrometeor mixing ratios
  real(kind=real_hmprec), pointer :: q_cl(:,:,:) => null()
  real(kind=real_hmprec), pointer :: q_rain(:,:,:) => null()
  real(kind=real_hmprec), pointer :: q_cf(:,:,:) => null()
  real(kind=real_hmprec), pointer :: q_snow(:,:,:) => null()
  real(kind=real_hmprec), pointer :: q_graup(:,:,:) => null()

  ! Liquid and ice cloud fractions, and combined cloud fraction
  real(kind=real_hmprec), pointer :: cf_liq(:,:,:) => null()
  real(kind=real_hmprec), pointer :: cf_ice(:,:,:) => null()
  real(kind=real_hmprec), pointer :: cf_bulk(:,:,:) => null()

  ! List of pointers to passive tracer fields for which
  ! convective transport is required
  type(fields_list_type), allocatable :: tracers(:)

  ! List of pointers to all the above fields,
  ! for the purposes of looping over them to do
  ! the same thing to all the fields
  type(fields_list_type), allocatable :: list(:)

end type fields_type


!----------------------------------------------------------------
! Types for defining ragged arrays of single fields, or generic
! super-arrays
!----------------------------------------------------------------
type :: ragged_array_type
  real(kind=real_cvprec), allocatable :: f(:)
end type ragged_array_type
type :: ragged_super_type
  real(kind=real_cvprec), allocatable :: super(:,:)
end type ragged_super_type


!----------------------------------------------------------------
! Addresses of fields in the lists and in compressed super-arrays
!----------------------------------------------------------------
! These are set in fields_set_addresses below, based on which
! fields are actually in use.

! Total number of primary fields in use, excluding tracers.
integer :: n_fields = 0

! Adresses of fields in the super-array
! (need to be set based on switches)

integer :: i_wind_u = 0
integer :: i_wind_v = 0
integer :: i_wind_w = 0

integer :: i_temperature = 0
integer :: i_q_vap = 0
integer :: i_q_cl = 0
integer :: i_q_rain = 0
integer :: i_q_cf = 0
integer :: i_q_snow = 0
integer :: i_q_graup = 0

integer :: i_cf_liq = 0
integer :: i_cf_ice = 0
integer :: i_cf_bulk = 0

! Addresses of tracers
integer, allocatable :: i_tracers(:)

! Addresses of the first and last condensed water species
! mixing ratio fields in the super-array
integer :: i_qc_first = 0
integer :: i_qc_last = 0

! Dummy array for unused pointers to point at, just to be safe.
real(kind=real_hmprec), allocatable, target :: dummy_full(:,:,:)

! List of the names of each field, stored as character strings
! (used for labelling diagnostics and error print-outs)
character(len=name_length), allocatable :: field_names(:)

! List of flags for whether each field is positive-only
! (e.g. mixing-ratios are not allowed to go negative)
logical, allocatable :: field_positive(:)


contains


!----------------------------------------------------------------
! Subroutine to set the address of each field in the super-arrays
! and pointer lists
!----------------------------------------------------------------
subroutine fields_set_addresses()

use comorph_constants_mod, only: l_cv_rain, l_cv_cf, l_cv_snow, l_cv_graup,    &
                                 l_cv_cloudfrac, n_tracers, n_cond_species,    &
                                 cond_params, k_bot_conv, k_top_conv,          &
                                 i_cond_cl, i_cond_rain,                       &
                                 i_cond_cf, i_cond_snow, i_cond_graup,         &
                                 tracer_positive

implicit none

! Character for storing tracer number
character(len=5) :: tr_num

! Loop counter
integer :: i_field


! Set flag to indicate that we don't need to call this routine again!
l_init_fields_type_mod = .true.

! Wind fields:
! u, v, w
i_wind_u = 1
i_wind_v = 2
i_wind_w = 3

! Thermodynamic fields:
! Fields which are always required:
! T, q
i_temperature = 4
i_q_vap = 5

! Number of fields so far:
n_fields = 5

! Condensed water fields...
! Use condensed water super-array addresses stored in constants
i_q_cl = n_fields + i_cond_cl
if ( l_cv_rain )  i_q_rain  = n_fields + i_cond_rain
if ( l_cv_cf )    i_q_cf    = n_fields + i_cond_cf
if ( l_cv_snow )  i_q_snow  = n_fields + i_cond_snow
if ( l_cv_graup ) i_q_graup = n_fields + i_cond_graup
n_fields = n_fields + n_cond_species

! Store address of first and last condensed water
! mixing ratio field
i_qc_first = i_q_cl
i_qc_last = i_qc_first + n_cond_species - 1

! Cloud-fraction fields only required on switch:
! cf_liq, cf_ice, cf_bulk
if ( l_cv_cloudfrac ) then
  i_cf_liq = n_fields + 1
  i_cf_ice = n_fields + 2
  i_cf_bulk = n_fields + 3
  n_fields = n_fields + 3
end if

! Tracers
if ( n_tracers > 0 ) then
  allocate( i_tracers(n_tracers) )
  do i_field = 1, n_tracers
    i_tracers(i_field) = n_fields + i_field
  end do
end if

! Note: tracers are not included in the count n_fields;
! structures which need to include tracers need to be
! dimensioned n_fields + n_tracers.
! This is because some instances of the fields structures don't
! need to include tracers even if tracers are being used.

! Allocate dummy array for unused fields to point at for safety
allocate( dummy_full(1,1,k_bot_conv:k_top_conv) )

! Allocate list of field names and list of flags for fields
! that aren't allowed to have negative values
allocate( field_names(n_fields+n_tracers) )
allocate( field_positive(n_fields+n_tracers) )

! Set the field names...
field_names(i_wind_u) = "wind_u"
field_names(i_wind_v) = "wind_v"
field_names(i_wind_w) = "wind_w"
field_names(i_temperature) = "temperature"
field_names(i_q_vap) = "q_vap"
! Use condensed water species names already set in cond_params
do i_field = 1, n_cond_species
  field_names(i_qc_first+i_field-1) = "q_" //                                  &
              trim(adjustl( cond_params(i_field)%pt % cond_name ))
end do
if ( l_cv_cloudfrac ) then
  field_names(i_cf_liq) = "cf_liq"
  field_names(i_cf_ice) = "cf_ice"
  field_names(i_cf_bulk) = "cf_bulk"
end if
if ( n_tracers > 0 ) then
  do i_field = 1, n_tracers
    write(tr_num,"(I5)") i_field
    field_names(i_tracers(i_field)) = "tracer_" //                             &
                                      trim(adjustl(tr_num))
  end do
end if

! Initialise positive-only flag to false for all fields
do i_field = 1, n_fields+n_tracers
  field_positive(i_field) = .false.
end do
! Set to true for those fields which are positive-only
field_positive(i_temperature) = .true.
field_positive(i_q_vap) = .true.
field_positive(i_qc_first:i_qc_last) = .true.
if ( l_cv_cloudfrac ) field_positive(i_cf_liq:i_cf_bulk) = .true.

if ( n_tracers > 0 ) then
  ! For tracers, this is user-defined, based on input list of logicals.
  ! If tracer_positive isn't allocated on input to comorph, field_positive
  ! is left as false (i.e. allow negative values for all tracers).
  if ( allocated(tracer_positive) ) then
    do i_field = 1, n_tracers
      field_positive(i_tracers(i_field)) = tracer_positive(i_field)
    end do
  end if
end if


return
end subroutine fields_set_addresses


!----------------------------------------------------------------
! Subroutine to nullify the pointers in a fields_type structure
!----------------------------------------------------------------
! This ought to be called before deallocating any of the arrays
! that the contained pointers might be pointing at
subroutine fields_nullify( fields )
use comorph_constants_mod, only: n_tracers
implicit none
type(fields_type), intent(in out) :: fields
integer :: i_tracer

fields % wind_u      => null()
fields % wind_v      => null()
fields % wind_w      => null()
fields % temperature => null()
fields % q_vap       => null()
fields % q_cl        => null()
fields % q_rain      => null()
fields % q_cf        => null()
fields % q_snow      => null()
fields % q_graup     => null()
fields % cf_liq      => null()
fields % cf_ice      => null()
fields % cf_bulk     => null()
if ( allocated( fields % tracers ) ) then
  do i_tracer = 1, n_tracers
    fields % tracers(i_tracer)%pt => null()
  end do
end if

return
end subroutine fields_nullify


!----------------------------------------------------------------
! Subroutine to assign a list of pointers to the full 3-D fields
!----------------------------------------------------------------
! This is used to allow a super do-loop to be applied when
! we want to do exactly the same thing to all the fields
subroutine fields_list_make( l_tracer, fields )

use comorph_constants_mod, only: l_cv_rain, l_cv_cf, l_cv_snow, l_cv_graup,    &
                     l_cv_cloudfrac, n_tracers,                                &
                     nx_full, ny_full, k_bot_conv, k_top_conv
use raise_error_mod, only: raise_fatal

implicit none

! Flag for whether the structure needs to store tracers
logical, intent(in) :: l_tracer

! Fields structure to be pointed at
type(fields_type), intent(in out) :: fields

! Total number of fields including tracers
integer :: n_fields_tot

! Bounds of the input array pointers
integer :: lb(3), ub(3)

! Loop counter
integer :: i_field

character(len=*), parameter :: routinename = "FIELDS_LIST_MAKE"


! Set number of fields
n_fields_tot = n_fields
if ( l_tracer )  n_fields_tot = n_fields_tot + n_tracers

! Allocate list of pointers
allocate( fields%list(n_fields_tot) )

! Assign pointers to fields

! Wind components
fields%list(i_wind_u)%pt => fields % wind_u
fields%list(i_wind_v)%pt => fields % wind_v
fields%list(i_wind_w)%pt => fields % wind_w

! Thermodynamic fields that are always required
fields%list(i_temperature)%pt => fields % temperature
fields%list(i_q_vap)%pt       => fields % q_vap
fields%list(i_q_cl)%pt        => fields % q_cl

! Optional water species fields if in use
! (if not in use, point the pointers at a dummy array for safety)
if ( l_cv_rain ) then
  fields%list(i_q_rain)%pt  => fields % q_rain
else
  fields % q_rain           => dummy_full
end if
if ( l_cv_cf ) then
  fields%list(i_q_cf)%pt    => fields % q_cf
else
  fields % q_cf             => dummy_full
end if
if ( l_cv_snow ) then
  fields%list(i_q_snow)%pt  => fields % q_snow
else
  fields % q_snow           => dummy_full
end if
if ( l_cv_graup ) then
  fields%list(i_q_graup)%pt => fields % q_graup
else
  fields % q_graup          => dummy_full
end if

if ( l_cv_cloudfrac ) then
  ! Cloud fractions
  fields%list(i_cf_liq)%pt  => fields % cf_liq
  fields%list(i_cf_ice)%pt  => fields % cf_ice
  fields%list(i_cf_bulk)%pt => fields % cf_bulk
else
  fields % cf_liq  => dummy_full
  fields % cf_ice  => dummy_full
  fields % cf_bulk => dummy_full
end if

! Tracers
if ( l_tracer .and. n_tracers > 0 ) then
  do i_field = 1, n_tracers
    fields%list(i_tracers(i_field)) % pt                                       &
      => fields % tracers(i_field) % pt
  end do
end if

! Now loop over all the fields and check that the contained
! pointers all actually point at arrays with the required extent
do i_field = 1, n_fields_tot
  if ( .not. associated( fields%list(i_field)%pt ) ) then
    call raise_fatal( routinename,                                             &
           "Required input primary field " //                                  &
           trim(adjustl(field_names(i_field))) //                              &
           "has not been assigned to any memory." )
  else
    lb = lbound( fields%list(i_field)%pt )
    ub = ubound( fields%list(i_field)%pt )
    if ( .not. ( lb(1) <= 1 .and. ub(1) >= nx_full                             &
           .and. lb(2) <= 1 .and. ub(2) >= ny_full                             &
           .and. lb(3) <= k_bot_conv .and. ub(3) >= k_top_conv ) ) then
      call raise_fatal( routinename,                                           &
             "Required input primary field " //                                &
             trim(adjustl(field_names(i_field))) //                            &
             "has insufficient extent / the wrong shape." )
    end if
  end if
end do


return
end subroutine fields_list_make


!----------------------------------------------------------------
! Subroutine to clear the above list of pointers
!----------------------------------------------------------------
subroutine fields_list_clear( fields )

implicit none

type(fields_type), intent(in out) :: fields

integer :: i_field

! Nullify the pointers
do i_field = 1, size(fields%list)
  fields%list(i_field)%pt => null()
end do

! Deallocate the list
deallocate( fields%list )

return
end subroutine fields_list_clear


!----------------------------------------------------------------
! Subroutine to convert the fields into variables which are
! conserved following dry-mass (at constant pressure)
!----------------------------------------------------------------
subroutine fields_k_conserved_vars( n_points, n_points_super,                  &
                                    n_fields_tot, l_reverse,                   &
                                    fields_k_super )

use comorph_constants_mod, only: real_cvprec, l_cv_cloudfrac, one

use calc_q_tot_mod, only: calc_q_tot
use set_cp_tot_mod, only: set_cp_tot
use calc_virt_temp_dry_mod, only: calc_virt_temp_dry

implicit none

! Array dimensions:

! Number of points where we actually want to do something
integer, intent(in) :: n_points
! Number of points in the compressed super-array
! (to avoid having to repeatedly deallocate and reallocate when
!  re-using the same array, it maybe bigger than needed here).
integer, intent(in) :: n_points_super
! Number of fields in the super-array
integer, intent(in) :: n_fields_tot

! Flag indicating whether to convert from normal form to
! conserved form, or vice versa
logical, intent(in) :: l_reverse

! Structure containing fields to be converted into conserved form
real(kind=real_cvprec), intent(in out) :: fields_k_super                       &
                                ( n_points_super, n_fields_tot )

! Work array storing various conversion factors
real(kind=real_cvprec) :: factor(n_points)

! Loop counter
integer :: ic, i_field


! Calculate dry-mass to wet-mass conversion factor = 1 + q_tot
call calc_q_tot( n_points, n_points_super,                                     &
                 fields_k_super(:,i_q_vap),                                    &
                 fields_k_super(:,i_qc_first:i_qc_last),                       &
                 factor )
do ic = 1, n_points
  factor(ic) = one + factor(ic)
end do

! If converting from momentum per unit dry-mass back to winds,
! take reciprocal of the conversion factor
if ( l_reverse ) then
  do ic = 1, n_points
    factor(ic) = one / factor(ic)
  end do
end if

! Convert between winds and momentum per unit dry-mass
! Loop over the 3 wind components in the super-array:
do i_field = i_wind_u, i_wind_w
  do ic = 1, n_points
    fields_k_super(ic,i_field) = fields_k_super(ic,i_field)                    &
                               * factor(ic)
  end do
end do

! Note: no conversion needed for the water species mixing ratios,
! as these are already conserved following dry-mass.

! Note: temperature is needed in its normal form for the
! conversion of cloud fractions to / from conserved form.
! Therefore when converting to conserved variables,
! need to convert cloud fractions before temperature,
! but when converting FROM conserved variables,
! need to convert temperature first.

if ( l_reverse ) then
  ! Converting from conserved form to normal:

  ! Calculate the total heat capacity of the parcel,
  ! per unit dry-mass
  call set_cp_tot( n_points, n_points_super,                                   &
                   fields_k_super(:,i_q_vap),                                  &
                   fields_k_super(:,i_qc_first:i_qc_last),                     &
                   factor )

  ! Convert from enthalpy per unit dry-mass to temperature
  do ic = 1, n_points
    fields_k_super(ic,i_temperature)                                           &
      = fields_k_super(ic,i_temperature) / factor(ic)
  end do

  ! If using cloud fractions...
  if ( l_cv_cloudfrac ) then

    ! Calculate parcel volume per unit dry-mass * pressure / R_dry
    !  = p / ( R_dry rho_dry )
    !  = T ( 1 + (R_vap/R_dry) q_vap )  "dry virtual temperature"
    call calc_virt_temp_dry( n_points,                                         &
                             fields_k_super(:,i_temperature),                  &
                             fields_k_super(:,i_q_vap),                        &
                             factor )

    ! Converting back to actual cloud fraction so take reciprocal
    do ic = 1, n_points
      factor(ic) = one / factor(ic)
    end do

    ! Convert from CF * Tv_dry to CF
    ! Loop over the 3 cloud-fraction fields in the super-array:
    do i_field = i_cf_liq, i_cf_bulk
      do ic = 1, n_points
        fields_k_super(ic,i_field) = fields_k_super(ic,i_field)                &
                                   * factor(ic)
      end do
    end do

  end if

else  ! ( l_reverse )
  ! Converting from normal to conserved form:

  if ( l_cv_cloudfrac ) then

    ! Calculate parcel volume per unit dry-mass * pressure / R_dry
    !  = p / ( R_dry rho_dry )
    !  = T ( 1 + (R_vap/R_dry) q_vap )  "dry virtual temperature"
    call calc_virt_temp_dry( n_points,                                         &
                             fields_k_super(:,i_temperature),                  &
                             fields_k_super(:,i_q_vap),                        &
                             factor )

    ! Convert from CF to CF * Tv_dry
    ! Loop over the 3 cloud-fraction fields in the super-array:
    do i_field = i_cf_liq, i_cf_bulk
      do ic = 1, n_points
        fields_k_super(ic,i_field) = fields_k_super(ic,i_field)                &
                                   * factor(ic)
      end do
    end do

  end if

  ! Calculate the total heat capacity of the parcel,
  ! per unit dry-mass
  call set_cp_tot( n_points, n_points_super,                                   &
                   fields_k_super(:,i_q_vap),                                  &
                   fields_k_super(:,i_qc_first:i_qc_last),                     &
                   factor )

  ! Convert from temperature to enthalpy per unit dry-mass
  do ic = 1, n_points
    fields_k_super(ic,i_temperature)                                           &
      = fields_k_super(ic,i_temperature) * factor(ic)
  end do

end if  ! ( l_reverse )

return
end subroutine fields_k_conserved_vars


!----------------------------------------------------------------
! Subroutine to adjust the fields under a change in pressure
!----------------------------------------------------------------
subroutine fields_k_pressure_adjust( n_points, n_points_super,                 &
                                     n_fields_tot,                             &
                                     fields_k_super,                           &
                                     pressure_1, pressure_2 )

use comorph_constants_mod, only: l_cv_cloudfrac, real_cvprec
use dry_adiabat_mod, only: dry_adiabat

implicit none

! Number of points
integer, intent(in) :: n_points

! Size of super-array
integer, intent(in) :: n_points_super

! Number of fields in the super-array
integer, intent(in) :: n_fields_tot

! Super-array containing fields to be adjusted
real(kind=real_cvprec), intent(in out) :: fields_k_super                       &
                                ( n_points_super, n_fields_tot )

! Pressure before and after the adjustment
real(kind=real_cvprec), intent(in) :: pressure_1(n_points)
real(kind=real_cvprec), intent(in) :: pressure_2(n_points)

! Work array storing exner ratio for adjusting temperature
real(kind=real_cvprec) :: exner_ratio(n_points)

! Loop counters
integer :: ic, i_field

! Compute exner ratio between the pressures at k and k2
call dry_adiabat( n_points, n_points_super,                                    &
                  pressure_1, pressure_2,                                      &
                  fields_k_super(:,i_q_vap),                                   &
                  fields_k_super(:,i_qc_first:i_qc_last),                      &
                  exner_ratio )

! Calculate temperature the air at k would have if moved to k2
do ic = 1, n_points
  fields_k_super(ic,i_temperature)                                             &
    = fields_k_super(ic,i_temperature) * exner_ratio(ic)
end do

! If using cloud-fractions, these currently store
! cloud volume per unit dry-mass ~ Tv_dry * CF.
! Therefore, these need to be adjusted in the same way
! as temperature
if ( l_cv_cloudfrac ) then
  do i_field = i_cf_liq, i_cf_bulk
    do ic = 1, n_points
      fields_k_super(ic,i_field)                                               &
        = fields_k_super(ic,i_field) * exner_ratio(ic)
    end do
  end do
end if


return
end subroutine fields_k_pressure_adjust


end module fields_type_mod
