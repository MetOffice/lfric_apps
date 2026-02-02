! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

#if !defined(LFRIC)
module assign_tracers_mod

use um_types, only: real_umphys

implicit none

contains

! Subroutine to assign a list of pointers to each tracer field
! for which convective transport is required in CoMorph.
!
! This mirrors the copying of the tracer fields into a
! a super-array in SL-advection and in the 6A convection scheme.
! But here, no data is actually copied; the tracer fields can be
! looped over generically by looping over the pointers in the
! list.  This avoids wasting large memory and CPU resources
! by pointlessly copying all the tracer fields, as is done
! elsewhere in the UM.
!
! Note: the number of tracer fields assigned here MUST be kept
! consistent with the total number of tracers counted in the
! routine tracer_total_var, as that sets the length of the
! list of pointers used.
subroutine assign_tracers(                                                     &
             aerosol, dust_div1, dust_div2,                                    &
             dust_div3, dust_div4, dust_div5, dust_div6,                       &
             so2, so4_aitken, so4_accu, so4_diss,                              &
             dms, nh3, soot_new, soot_aged, soot_cld,                          &
             bmass_new, bmass_aged, bmass_cld,                                 &
             ocff_new, ocff_aged, ocff_cld, nitr_acc, nitr_diss,               &
             co2, free_tracers, tracer_ukca, ozone_tracer,                     &
             fields_np1 )

use nlsizes_namelist_mod, only: tr_vars, tr_ukca, tr_ukca_actv
use ukca_scavenging_mod, only:  tracer_info
use atm_fields_bounds_mod, only: tdims_s
use run_aerosol_mod, only: l_sulpc_so2, l_sulpc_dms, l_sulpc_nh3,              &
                           l_soot, l_ocff, l_biomass, l_nitrate
use cv_run_mod, only: l_murk_conv
use rad_input_mod, only: l_use_cariolle
use carbon_options_mod, only: l_co2_interactive
use dust_parameters_mod, only: l_dust, l_twobin_dust

use comorph_constants_mod, only: n_tracers, newline
use fields_type_mod, only: fields_type
use raise_error_mod, only: raise_fatal

implicit none

! Aerosol tracer (called "murk" in atm_step).
real(kind=real_umphys), target, intent(in) :: aerosol                          &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)

! Dust fields
real(kind=real_umphys), target, intent(in) :: dust_div1                        &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: dust_div2                        &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: dust_div3                        &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: dust_div4                        &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: dust_div5                        &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: dust_div6                        &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)

! Chemistry species
real(kind=real_umphys), target, intent(in) :: so2                              &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: so4_aitken                       &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: so4_accu                         &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: so4_diss                         &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: dms                              &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: nh3                              &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: soot_new                         &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: soot_aged                        &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: soot_cld                         &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: bmass_new                        &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: bmass_aged                       &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: bmass_cld                        &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: ocff_new                         &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: ocff_aged                        &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: ocff_cld                         &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: nitr_acc                         &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: nitr_diss                        &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)
real(kind=real_umphys), target, intent(in) :: co2                              &
                             (tdims_s%i_start:tdims_s%i_end,                   &
                              tdims_s%j_start:tdims_s%j_end,                   &
                              tdims_s%k_start:tdims_s%k_end)

! Free tracers
real(kind=real_umphys), target, intent(in) :: free_tracers                     &
                       (tdims_s%i_start:tdims_s%i_end,                         &
                        tdims_s%j_start:tdims_s%j_end,                         &
                        tdims_s%k_start:tdims_s%k_end, tr_vars)

! UKCA tracers
real(kind=real_umphys), target, intent(in) :: tracer_ukca                      &
                       (tdims_s%i_start:tdims_s%i_end,                         &
                        tdims_s%j_start:tdims_s%j_end,                         &
                        tdims_s%k_start:tdims_s%k_end, tr_ukca)

! Ozone tracer
real(kind=real_umphys), target, intent(in) :: ozone_tracer                     &
                       (tdims_s%i_start:tdims_s%i_end,                         &
                        tdims_s%j_start:tdims_s%j_end,                         &
                        tdims_s%k_start:tdims_s%k_end)

! Derived-type structure containing pointers
type(fields_type), intent(in out) :: fields_np1

! Array bounds, for bounds-specification when assigning pointer
integer :: lb(4)

! Counters for tracer fields
integer :: i_tracer, i_field

character(len=*), parameter :: routinename = "ASSIGN_TRACERS"


! Initialise counter for tracer fields
i_tracer = 0

! Dust scheme
if (l_dust) then
  fields_np1 % tracers(i_tracer+1)%pt => dust_div1
  fields_np1 % tracers(i_tracer+2)%pt => dust_div2
  if (l_twobin_dust) then
    i_tracer = i_tracer + 2
  else
    fields_np1 % tracers(i_tracer+3)%pt => dust_div3
    fields_np1 % tracers(i_tracer+4)%pt => dust_div4
    fields_np1 % tracers(i_tracer+5)%pt => dust_div5
    fields_np1 % tracers(i_tracer+6)%pt => dust_div6
    i_tracer = i_tracer + 6
  end if
end if

! Sulphur cycle
if (l_sulpc_so2) then
  fields_np1 % tracers(i_tracer+1)%pt => so2
  fields_np1 % tracers(i_tracer+2)%pt => so4_aitken
  fields_np1 % tracers(i_tracer+3)%pt => so4_accu
  fields_np1 % tracers(i_tracer+4)%pt => so4_diss
  i_tracer = i_tracer + 4
  if (l_sulpc_dms) then
    i_tracer = i_tracer + 1
    fields_np1 % tracers(i_tracer)%pt => dms
  end if
  if (l_sulpc_nh3) then
    i_tracer = i_tracer + 1
    fields_np1 % tracers(i_tracer)%pt => nh3
  end if
end if

! Soot
if (l_soot) then
  fields_np1 % tracers(i_tracer+1)%pt => soot_new
  fields_np1 % tracers(i_tracer+2)%pt => soot_aged
  fields_np1 % tracers(i_tracer+3)%pt => soot_cld
  i_tracer = i_tracer + 3
end if

! Biomass burning aerosol
if (l_biomass) then
  fields_np1 % tracers(i_tracer+1)%pt => bmass_new
  fields_np1 % tracers(i_tracer+2)%pt => bmass_aged
  fields_np1 % tracers(i_tracer+3)%pt => bmass_cld
  i_tracer = i_tracer + 3
end if

! Interactive CO2
if (l_co2_interactive) then
  i_tracer = i_tracer + 1
  fields_np1 % tracers(i_tracer)%pt => co2
end if

! Organic Carbon from Fossil Fuels
if (l_ocff) then
  fields_np1 % tracers(i_tracer+1)%pt => ocff_new
  fields_np1 % tracers(i_tracer+2)%pt => ocff_aged
  fields_np1 % tracers(i_tracer+3)%pt => ocff_cld
  i_tracer = i_tracer + 3
end if

! Nitrate
if (l_nitrate) then
  fields_np1 % tracers(i_tracer+1)%pt => nitr_acc
  fields_np1 % tracers(i_tracer+2)%pt => nitr_diss
  i_tracer = i_tracer + 2
end if

! Interactive ozone
if (l_use_cariolle) then
  i_tracer = i_tracer + 1
  fields_np1 % tracers(i_tracer)%pt => ozone_tracer
end if

! "Murk" (simplified aerosol scheme with just one prognostic)
if (l_murk_conv) then
  i_tracer = i_tracer + 1
  fields_np1 % tracers(i_tracer)%pt => aerosol
end if

! Free tracers
if ( tr_vars > 0 ) then
  lb = lbound(free_tracers)
  do i_field = 1, tr_vars
    i_tracer = i_tracer + 1
    fields_np1 % tracers(i_tracer)%pt( lb(1):, lb(2):, lb(3): )                &
      => free_tracers(:,:,:,i_field)
    ! Note: this pointer assignment is technically referencing
    ! an array section, rather than a full array.  In fortran,
    ! this results in the bounds specs being lost and
    ! resetting the array lower-bounds to 1 in the pointer.
    ! Therefore, the above use of bounds specification is needed
    ! in order to retain the knowledge of the haloes in the
    ! array pointer.
  end do
end if

! UKCA fields
if ( tr_ukca_actv > 0 ) then

  tracer_info%i_ukca_first = i_tracer + 1

  lb = lbound(tracer_ukca)
  do i_field = 1, tr_ukca
    i_tracer = i_tracer + 1
    fields_np1 % tracers(i_tracer)%pt( lb(1):, lb(2):, lb(3): )                &
      => tracer_ukca(:,:,:,i_field)
    ! Use of bounds specification needed to retain the array
    ! lower bounds in case of haloes; see the note above for
    ! the free tracers super-array.
  end do

  tracer_info%i_ukca_last = i_tracer

end if

! Check the final tally is correct!
! Note: if i_tracer exceeded n_tracers above, the model will have
! already referenced fields_np1 % tracers out of bounds, causing
! either segmentation fault or "unexpected behaviour".
if ( .not. i_tracer == n_tracers ) then
  call raise_fatal( routinename,                                               &
         "Number of tracer fields assigned does not equal "   //               &
         "the number of tracer fields counted in "   //newline//               &
         "tracer_total_var." )
end if


return
end subroutine assign_tracers


end module assign_tracers_mod
#endif
