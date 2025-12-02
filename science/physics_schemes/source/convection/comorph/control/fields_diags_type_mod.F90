! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module fields_diags_type_mod

use diag_type_mod, only: diag_type, diag_list_type

implicit none

!----------------------------------------------------------------
! Type definition for a structure to store diagnostics for
! each primary field
! (eg resolved-scale increments, or parcel properties)
!----------------------------------------------------------------
! This follows the structure of fields_type in fields_type_mod
type :: fields_diags_type

  ! Number of diagnostics in this structure requested for output
  integer :: n_diags = 0

  ! Number of primary field diagnostics requested (excludes
  ! derived fields cp_tot, q_tot which aren't in fields_type)
  integer :: n_diags_fields = 0

  ! List of pointers to active diagnostic structures
  ! (points to subset of the main list stored in the
  !  comorph_diags_type structure)
  type(diag_list_type), pointer :: list(:) => null()

  ! Wind components
  type(diag_type) :: wind_u
  type(diag_type) :: wind_v
  type(diag_type) :: wind_w

  ! Temperature
  type(diag_type) :: temperature

  ! Water vapour mixing ratio
  type(diag_type) :: q_vap

  ! Hydrometeor mixing ratios
  type(diag_type) :: q_cl
  type(diag_type) :: q_rain
  type(diag_type) :: q_cf
  type(diag_type) :: q_snow
  type(diag_type) :: q_graup

  ! Liquid and ice cloud fractions, and combined cloud fraction
  type(diag_type) :: cf_liq
  type(diag_type) :: cf_ice
  type(diag_type) :: cf_bulk

  ! Tracers
  type(diag_type), allocatable :: tracers(:)

  ! Extra derived fields; some are needed to compute means over
  ! multiple parcels

  ! Total-water mixing-ratio
  type(diag_type) :: q_tot
  ! Total heat capacity per unit dry-mass
  type(diag_type) :: cp_tot
  ! Virtual temperature
  type(diag_type) :: virt_temp
  ! Virtual temperature excess (minus grid-mean)
  type(diag_type) :: virt_temp_excess
  ! Saturation water-vapour mixing ratio (w.r.t. liquid water)
  type(diag_type) :: qsat_liq
  ! Relative humidity (w.r.t. liquid water)
  type(diag_type) :: rel_hum_liq

end type fields_diags_type


contains


!----------------------------------------------------------------
! Subroutine to count the number of active diagnostics in a
! fields_diags structure and store a list of pointers to them
!----------------------------------------------------------------
! Note: this routine gets called twice
!   1st call (l_count_diags = .true.):
!     Check whether each diag is requested and count how many
!   2nd call (l_count_diags = .false.):
!     Set other properties for requested diags, and assign
!     pointers from active diags list.
subroutine fields_diags_assign( parent_name, l_count_diags, doms,              &
                                fields_diags, parent_list,                     &
                                parent_i_diag, i_super, l_mean )

use comorph_constants_mod, only: n_tracers, name_length
use diag_type_mod, only: diag_list_type, diag_assign, dom_type
use fields_type_mod, only: i_wind_u, i_wind_v, i_wind_w,                       &
                           i_temperature, i_q_vap,                             &
                           i_q_cl, i_q_rain,                                   &
                           i_q_cf, i_q_snow, i_q_graup,                        &
                           i_cf_liq, i_cf_ice, i_cf_bulk,                      &
                           i_tracers, n_fields, field_names

implicit none

! Character string to prepend onto all the diagnostic names
! in here; e.g. "updraft" or "downdraft"
character(len=name_length), intent(in) :: parent_name

! Flag for 1st call to this routine, in which we just check
! whether each diag is requested and count up number requested
logical, intent(in) :: l_count_diags

! Structure containing the allowed domain profiles for the diags
type(dom_type), intent(in) :: doms

! The main diagnostics structure, containing diag structures
! for all diagnostics, including those not in use
type(fields_diags_type), target, intent(in out) :: fields_diags

! List of active diagnostics in the parent structure
type(diag_list_type), target, intent(in out) :: parent_list(:)

! Counter for requested diagnostics in the parent structure
integer, intent(in out) :: parent_i_diag

! Counter for addresses of fields in compressed diags super-array
integer, intent(in out) :: i_super

! Flag can be input to indicate that means of any requested
! fields need to be computed, averaging over multiple layers
! or types of convection.  When this is the case, q_tot
! and cp_tot are needed for computing the conserved means, even
! if not requested for output.
logical, optional, intent(in) :: l_mean

! Local list of pointers pointing to all the diag structures,
! to loop over them
type(diag_list_type) :: list_all(n_fields+n_tracers)

! Total number of fields (primary fields + tracers)
integer :: n_fields_tot

! Name for each diagnostic
character(len=name_length) :: diag_name

! Counters
integer :: i_diag, i_field


! Set fields_diags list to point at a section of the parent list
if ( l_count_diags ) then
  ! Not used in the first call, but needs to be assigned
  ! to avoid error when passing it into diag_assign
  fields_diags % list => parent_list
else
  ! Assign to subset address in parent:
  fields_diags % list => parent_list                                           &
    ( parent_i_diag + 1 :                                                      &
      parent_i_diag + fields_diags % n_diags )
end if

! Setup a list of pointers to all the diag structures contained
! in fields_diags, so we can loop over them
list_all(i_wind_u)%pt      => fields_diags%wind_u
list_all(i_wind_v)%pt      => fields_diags%wind_v
list_all(i_wind_w)%pt      => fields_diags%wind_w
list_all(i_temperature)%pt => fields_diags%temperature
list_all(i_q_vap)%pt       => fields_diags%q_vap
list_all(i_q_cl)%pt        => fields_diags%q_cl
if ( i_q_rain>0 )  list_all(i_q_rain)%pt  => fields_diags%q_rain
if ( i_q_cf>0 )    list_all(i_q_cf)%pt    => fields_diags%q_cf
if ( i_q_snow>0 )  list_all(i_q_snow)%pt  => fields_diags%q_snow
if ( i_q_graup>0 ) list_all(i_q_graup)%pt => fields_diags%q_graup
if ( i_cf_liq>0 )  list_all(i_cf_liq)%pt  => fields_diags%cf_liq
if ( i_cf_ice>0 )  list_all(i_cf_ice)%pt  => fields_diags%cf_ice
if ( i_cf_bulk>0 ) list_all(i_cf_bulk)%pt => fields_diags%cf_bulk
if ( n_tracers > 0 .and.                                                       &
     allocated( fields_diags % tracers ) ) then
  n_fields_tot = n_fields + n_tracers
  do i_field = 1, n_tracers
    list_all(i_tracers(i_field))%pt =>                                         &
                                    fields_diags%tracers(i_field)
  end do
else
  n_fields_tot = n_fields
end if


! Initialise local requested diagnostics counter to zero
i_diag = 0

! Process each diagnostic within the fields_diags structure...
do i_field = 1, n_fields_tot

  ! Construct full name of the diagnostic using field_names
  diag_name = trim(adjustl(parent_name)) // "_" //                             &
         trim(adjustl(field_names(i_field)))

  ! Process diagnostic request and assign meta-data
  call diag_assign( diag_name, l_count_diags, doms, list_all(i_field)%pt,      &
         fields_diags%list, i_diag, i_super, i_field=i_field )

end do

! Save number of primary fields and tracers requested
if ( l_count_diags )  fields_diags % n_diags_fields = i_diag

! Extra derived fields; some are needed to compute means over
! multiple parcels
diag_name = trim(adjustl(parent_name)) // "_q_tot"
call diag_assign( diag_name, l_count_diags, doms, fields_diags%q_tot,          &
                  fields_diags%list, i_diag, i_super )
diag_name = trim(adjustl(parent_name)) // "_cp_tot"
call diag_assign( diag_name, l_count_diags, doms, fields_diags%cp_tot,         &
                  fields_diags%list, i_diag, i_super )
diag_name = trim(adjustl(parent_name)) // "_virt_temp"
call diag_assign( diag_name, l_count_diags, doms, fields_diags%virt_temp,      &
                  fields_diags%list, i_diag, i_super )
diag_name = trim(adjustl(parent_name)) // "_virt_temp_excess"
call diag_assign( diag_name, l_count_diags, doms, fields_diags                 &
                                       % virt_temp_excess,                     &
                  fields_diags%list, i_diag, i_super )
diag_name = trim(adjustl(parent_name)) // "_qsat_liq"
call diag_assign( diag_name, l_count_diags, doms, fields_diags%qsat_liq,       &
                  fields_diags%list, i_diag, i_super )
diag_name = trim(adjustl(parent_name)) // "_rel_hum_liq"
call diag_assign( diag_name, l_count_diags, doms, fields_diags%rel_hum_liq,    &
                  fields_diags%list, i_diag, i_super )


! Increment parent active diagnostics counter
parent_i_diag = parent_i_diag + i_diag


! If this is the initial call to this routine to count the
! number of requested diags
if ( l_count_diags ) then

  ! Set number of active diags in this fields_diags structure
  fields_diags % n_diags = i_diag

  ! Nullify the list pointer ready to be reassigned next call
  fields_diags % list => null()

  ! Set flag true for diags needed for calculating other diags...

  if ( present(l_mean) ) then
    if ( l_mean ) then
      ! Fields needed for computing conserved means of other
      ! fields, even if not requested...

      ! If any cloud-fractions requested, then the temperature
      ! and water-vapour are also needed for computing the
      ! conversion factor virt_temp_dry.
      if ( fields_diags % cf_liq % flag .or.                                   &
           fields_diags % cf_ice % flag .or.                                   &
           fields_diags % cf_bulk % flag ) then
        fields_diags % temperature % flag = .true.
        fields_diags % q_vap % flag = .true.
      end if

      ! If any wind fields were requested, then the
      ! total-water is also needed.
      if ( fields_diags % wind_u % flag .or.                                   &
           fields_diags % wind_v % flag .or.                                   &
           fields_diags % wind_w % flag ) then
        fields_diags % q_tot % flag = .true.
      end if

      ! If temperature requested,
      ! then the total heat capacity is also needed.
      if ( fields_diags % temperature % flag ) then
        fields_diags % cp_tot % flag = .true.
      end if

      ! Note: the 2nd sweep with l_count_diags=.false. will now
      ! account for flag when setting super-array addresses
      ! in diag_assign.

    end if  ! ( l_mean )
  end if  ! ( present(l_mean) )

end if  ! ( l_count_diags )


return
end subroutine fields_diags_assign


!----------------------------------------------------------------
! Subroutine to copy primary field properties into a super-array
! for diagnostic output
!----------------------------------------------------------------
subroutine fields_diags_copy( n_points, n_points_fields, n_points_diag,        &
                              n_diags, n_fields_tot,                           &
                              l_conserved_form,                                &
                              fields_diags, fields_super,                      &
                              env_virt_temp, pressure,                         &
                              diags_super )

use comorph_constants_mod, only: real_cvprec
use fields_type_mod, only: i_temperature, i_q_vap,                             &
                           i_qc_first, i_qc_last
use calc_q_tot_mod, only: calc_q_tot
use set_cp_tot_mod, only: set_cp_tot
use calc_virt_temp_mod, only: calc_virt_temp
use set_qsat_mod, only: set_qsat_liq

implicit none

! Number of points
integer, intent(in) :: n_points

! Size of the fields super-array (maybe bigger than needed,
! to save having to deallocate and reallocate)
integer, intent(in) :: n_points_fields

! Size of the diagnostics super-array
integer, intent(in) :: n_points_diag

! Number of diagnostics in diags_super
! (maybe larger than number of fields diags to copy here, as the
!  super-array may also include other diags not processed here)
integer, intent(in) :: n_diags

! Number of fields in fields_super
integer, intent(in) :: n_fields_tot

! Flag for whether or not the input primary fields (fields_super)
! are in conserved variable form
logical, intent(in) :: l_conserved_form

! Structure containing meta-data for fields diags
type(fields_diags_type), intent(in) :: fields_diags

! Input super-array containing primary fields to be copied.
! NOTE: the diagnostics system assumes these are in
! conserved variable form at this point!
real(kind=real_cvprec), intent(in) :: fields_super                             &
                                ( n_points_fields, n_fields_tot )

! Environment virtual temperature; only used if virtual
! temperature excess diagnostic is requested
real(kind=real_cvprec), intent(in) :: env_virt_temp(n_points)

! Ambient air pressure; used for calculating qsat if requested
real(kind=real_cvprec), intent(in) :: pressure(n_points)

! Diagnostics super-array to copy fields into
real(kind=real_cvprec), intent(in out) :: diags_super                          &
                                ( n_points_diag, n_diags )

! Work array for calculations
real(kind=real_cvprec) :: work(n_points)

! Loop counters
integer :: ic, i_req, i_diag, i_field

! Loop over requested primary fields
do i_req = 1, fields_diags % n_diags_fields
  ! Extract diag-array index and fields-array index
  i_diag  = fields_diags % list(i_req)%pt % i_super
  i_field = fields_diags % list(i_req)%pt % i_field
  ! Copy field
  do ic = 1, n_points
    diags_super(ic,i_diag) = fields_super(ic,i_field)
  end do
end do

! Temperature and q_vap might be required even if not requested
if ( fields_diags % temperature % flag ) then
  if ( .not. fields_diags % temperature % l_req ) then
    i_diag = fields_diags % temperature % i_super
    do ic = 1, n_points
      diags_super(ic,i_diag) = fields_super(ic,i_temperature)
    end do
  end if
end if
if ( fields_diags % q_vap % flag ) then
  if ( .not. fields_diags % q_vap % l_req ) then
    i_diag = fields_diags % q_vap % i_super
    do ic = 1, n_points
      diags_super(ic,i_diag) = fields_super(ic,i_q_vap)
    end do
  end if
end if

! Calculate total-water mixing-ratio if required
if ( fields_diags % q_tot % flag ) then
  i_diag = fields_diags % q_tot % i_super
  call calc_q_tot( n_points, n_points_fields,                                  &
                   fields_super(:,i_q_vap),                                    &
                   fields_super(:,i_qc_first:i_qc_last),                       &
                   diags_super(:,i_diag) )
end if

! Calculate total heat capacity if required
if ( fields_diags % cp_tot % flag ) then
  i_diag = fields_diags % cp_tot % i_super
  call set_cp_tot( n_points, n_points_fields,                                  &
                   fields_super(:,i_q_vap),                                    &
                   fields_super(:,i_qc_first:i_qc_last),                       &
                   diags_super(:,i_diag) )
end if


! Derived diagnostics that depend on temperature...
if ( fields_diags % virt_temp % flag .or.                                      &
     fields_diags % virt_temp_excess % flag .or.                               &
     fields_diags % qsat_liq % flag .or.                                       &
     fields_diags % rel_hum_liq % flag ) then

  ! If primary fields are in conserved variable form,
  ! and we've requested a diagnostic that needs T,
  ! convert cp*T to just T
  if ( l_conserved_form ) then
    call set_cp_tot( n_points, n_points_fields,                                &
                     fields_super(:,i_q_vap),                                  &
                     fields_super(:,i_qc_first:i_qc_last),                     &
                     work )
    do ic = 1, n_points
      work(ic) = fields_super(ic,i_temperature) / work(ic)
    end do
  else
    do ic = 1, n_points
      work(ic) = fields_super(ic,i_temperature)
    end do
  end if

  ! Virtual temperature diagnostic:
  if ( fields_diags % virt_temp % flag ) then
    i_diag = fields_diags % virt_temp % i_super
    call calc_virt_temp( n_points, n_points_fields,                            &
                         work, fields_super(:,i_q_vap),                        &
                         fields_super(:,i_qc_first:i_qc_last),                 &
                         diags_super(:,i_diag) )
  end if

  ! Virtual temperature excess diagnostic:
  if ( fields_diags % virt_temp_excess % flag ) then
    i_diag = fields_diags % virt_temp_excess % i_super
    call calc_virt_temp( n_points, n_points_fields,                            &
                         work, fields_super(:,i_q_vap),                        &
                         fields_super(:,i_qc_first:i_qc_last),                 &
                         diags_super(:,i_diag) )
    do ic = 1, n_points
      diags_super(ic,i_diag) = diags_super(ic,i_diag)                          &
                             - env_virt_temp(ic)
    end do
  end if

  ! Saturation mixing ratio w.r.t. liquid water
  if ( fields_diags % qsat_liq % flag ) then
    i_diag = fields_diags % qsat_liq % i_super
    call set_qsat_liq( n_points, work, pressure, diags_super(:,i_diag) )
  end if

  ! Relative humidity w.r.t. liquid water
  if ( fields_diags % rel_hum_liq % flag ) then
    i_diag = fields_diags % rel_hum_liq % i_super
    call set_qsat_liq( n_points, work, pressure, diags_super(:,i_diag) )
    do ic = 1, n_points
      diags_super(ic,i_diag) = fields_super(ic,i_q_vap)                        &
                             / diags_super(ic,i_diag)
    end do
  end if

end if


return
end subroutine fields_diags_copy


!----------------------------------------------------------------
! Subroutine to convert compressed primary field diagnostics
! between normal form and conserved form
!----------------------------------------------------------------
! This routine tediously duplicates the contents of
! fields_k_conserved_vars in fields_type_mod, but with a key
! difference that it needs to check whether each field is
! actually needed, rather than always converting everything.
subroutine fields_diags_conserved_vars(                                        &
             n_points, n_points_diag, n_diags_super, l_reverse,                &
             fields_diags, diags_super )

use comorph_constants_mod, only: real_cvprec, one
use calc_virt_temp_dry_mod, only: calc_virt_temp_dry

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of points in the super-array
! (maybe larger than n_points due to reuse of work arrays)
integer, intent(in) :: n_points_diag

! Number of diagnostics in the super-array
integer, intent(in) :: n_diags_super

! Flag indicating whether to convert from normal form to
! conserved form, or vice versa
logical, intent(in) :: l_reverse

! Structure containing meta-data for diags
type(fields_diags_type), intent(in) :: fields_diags

! Super-array containing fields to be converted
real(kind=real_cvprec), intent(in out) :: diags_super                          &
                               ( n_points_diag, n_diags_super )

! Conversion factor
real(kind=real_cvprec) :: factor(n_points)

! Super-array addresses
integer :: i_q_tot
integer :: i_wind
integer :: i_cp_tot
integer :: i_temp
integer :: i_q_vap
integer :: i_cf

! Loop counter
integer :: ic


! If any wind components requested
if ( fields_diags % wind_u % flag .or.                                         &
     fields_diags % wind_v % flag .or.                                         &
     fields_diags % wind_w % flag ) then
  ! Calculate conversion factor (1 + q_tot)
  i_q_tot = fields_diags % q_tot % i_super
  do ic = 1, n_points
    factor(ic) = one + diags_super(ic,i_q_tot)
  end do
  ! If converting from momentum per unit dry-mass back to winds,
  ! take reciprocal of the conversion factor
  if ( l_reverse ) then
    do ic = 1, n_points
      factor(ic) = one / factor(ic)
    end do
  end if
  ! Convert between momentum per unit dry-mass / winds
  if ( fields_diags % wind_u % flag ) then
    i_wind = fields_diags % wind_u % i_super
    do ic = 1, n_points
      diags_super(ic,i_wind) = diags_super(ic,i_wind)*factor(ic)
    end do
  end if
  if ( fields_diags % wind_v % flag ) then
    i_wind = fields_diags % wind_v  % i_super
    do ic = 1, n_points
      diags_super(ic,i_wind) = diags_super(ic,i_wind)*factor(ic)
    end do
  end if
  if ( fields_diags % wind_w % flag ) then
    i_wind = fields_diags % wind_w  % i_super
    do ic = 1, n_points
      diags_super(ic,i_wind) = diags_super(ic,i_wind)*factor(ic)
    end do
  end if
end if

! Note: temperature is needed in its normal form for the
! conversion of cloud fractions to / from conserved form.
! Therefore when converting to conserved variables,
! need to convert cloud fractions before temperature,
! but when converting FROM conserved variables,
! need to convert temperature first.

if ( l_reverse ) then
  ! Converting from conserved form to normal:

  ! If temperature requested
  if ( fields_diags % temperature % flag ) then
    ! Conversion factor is 1 / cp_tot
    i_cp_tot = fields_diags % cp_tot % i_super
    ! Convert enthalpy per unit dry-mas back to temperature
    i_temp = fields_diags % temperature % i_super
    do ic = 1, n_points
      diags_super(ic,i_temp) = diags_super(ic,i_temp)                          &
                             / diags_super(ic,i_cp_tot)
    end do
  end if

  ! If any cloud-fractions requested
  if ( fields_diags % cf_liq % flag .or.                                       &
       fields_diags % cf_ice % flag .or.                                       &
       fields_diags % cf_bulk % flag ) then
    ! Calculate 1/Tv_dry
    i_temp = fields_diags % temperature % i_super
    i_q_vap = fields_diags % q_vap % i_super
    call calc_virt_temp_dry( n_points,                                         &
                             diags_super(:,i_temp),                            &
                             diags_super(:,i_q_vap),                           &
                             factor )
    do ic = 1, n_points
      factor(ic) = one / factor(ic)
    end do
    ! Convert Tv_dry*CF back to CF
    if ( fields_diags % cf_liq % flag ) then
      i_cf = fields_diags % cf_liq % i_super
      do ic = 1, n_points
        diags_super(ic,i_cf) = diags_super(ic,i_cf) * factor(ic)
      end do
    end if
    if ( fields_diags % cf_ice % flag ) then
      i_cf = fields_diags % cf_ice % i_super
      do ic = 1, n_points
        diags_super(ic,i_cf) = diags_super(ic,i_cf) * factor(ic)
      end do
    end if
    if ( fields_diags % cf_bulk % flag ) then
      i_cf = fields_diags % cf_bulk % i_super
      do ic = 1, n_points
        diags_super(ic,i_cf) = diags_super(ic,i_cf) * factor(ic)
      end do
    end if
  end if

else  ! ( l_reverse )
  ! Converting from normal to conserved form:

  ! If any cloud-fractions requested
  if ( fields_diags % cf_liq % flag .or.                                       &
       fields_diags % cf_ice % flag .or.                                       &
       fields_diags % cf_bulk % flag ) then
    ! Calculate Tv_dry
    i_temp = fields_diags % temperature % i_super
    i_q_vap = fields_diags % q_vap % i_super
    call calc_virt_temp_dry( n_points,                                         &
                             diags_super(:,i_temp),                            &
                             diags_super(:,i_q_vap),                           &
                             factor )
    ! Convert CF to Tv_dry*CF
    if ( fields_diags % cf_liq % flag ) then
      i_cf = fields_diags % cf_liq % i_super
      do ic = 1, n_points
        diags_super(ic,i_cf) = diags_super(ic,i_cf) * factor(ic)
      end do
    end if
    if ( fields_diags % cf_ice % flag ) then
      i_cf = fields_diags % cf_ice % i_super
      do ic = 1, n_points
        diags_super(ic,i_cf) = diags_super(ic,i_cf) * factor(ic)
      end do
    end if
    if ( fields_diags % cf_bulk % flag ) then
      i_cf = fields_diags % cf_bulk % i_super
      do ic = 1, n_points
        diags_super(ic,i_cf) = diags_super(ic,i_cf) * factor(ic)
      end do
    end if
  end if

  ! If temperature requested
  if ( fields_diags % temperature % flag ) then
    ! Conversion factor is cp_tot
    i_cp_tot = fields_diags % cp_tot % i_super
    ! Convert enthalpy per unit dry-mas back to temperature
    i_temp = fields_diags % temperature % i_super
    do ic = 1, n_points
      diags_super(ic,i_temp) = diags_super(ic,i_temp)                          &
                             * diags_super(ic,i_cp_tot)
    end do
  end if

end if  ! ( l_reverse )


return
end subroutine fields_diags_conserved_vars


end module fields_diags_type_mod
