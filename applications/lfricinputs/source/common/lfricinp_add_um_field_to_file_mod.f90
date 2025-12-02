! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_add_um_field_to_file_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only: real64, int64

! lfricinp modules
use lfricinp_stashmaster_mod,          only: get_stashmaster_item, grid,   &
                                             p_points, u_points, v_points, &
                                             ozone_points,                 &
                                             ppfc,                         &
                                             sm_lbvc => lbvc,              &
                                             cfff,                         &
                                             levelt,                       &
                                             rho_levels, theta_levels,     &
                                             single_level,                 &
                                             cfll,                         &
                                             datat
use lfricinp_um_parameters_mod,        only: um_imdi, um_rmdi,             &
                                             rh_polelat, rh_polelong,      &
                                             ldc_zsea_theta, ldc_zsea_rho, &
                                             ldc_c_theta, ldc_c_rho,       &
                                             rh_deltaEW, rh_deltaNS,       &
                                             ih_model_levels
use lfricinp_grid_type_mod,            only: lfricinp_grid_type
use lfricinp_um_level_codes_mod,       only: lfricinp_get_first_level_num
use lfricinp_check_shumlib_status_mod, only: shumlib


! lfric modules
use log_mod, only : log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR, &
                    log_scratch_space

! shumlib modules
use f_shum_fieldsfile_mod, only: f_shum_fixed_length_header_len
use f_shum_file_mod, only: shum_file_type
use f_shum_field_mod, only: shum_field_type
use f_shum_lookup_indices_mod, only: &
                                    lbyr, lbmon, lbdat, lbhr, lbmin, lbday, lbsec, lbyrd, lbmond, lbdatd, &
    lbhrd, lbmind, lbdayd, lbsecd, lbtim, lbft, lbcode, lbhem,    &
    lbrow, lbnpt, lbpack, lbrel, lbfc, lbcfc, lbproc, lbvc, lbrvc, &
    lbtyp, lblev, lbrsvd1, lbrsvd2,        &
    lbrsvd3, lbrsvd4, lbsrce, lbuser1, lbuser4, &
    lbuser7, bulev, bhulev, brsvd3, brsvd4, bdatum, bacc, blev,  &
    brlev, bhlev, bhrlev, bplat, bplon, bgor, bzy, bdy, bzx, bdx, bmdi,   &
    bmks

use f_shum_fixed_length_header_indices_mod, only:                        &
                                                                        vert_coord_type, horiz_grid_type, &
 dataset_type, run_identifier, calendar, projection_number, model_version, &
 grid_staggering, sub_model, t1_year, t1_month, t1_day, t1_hour, t1_minute, &
 t1_second, t2_year, t2_month, t2_day, t2_hour, t2_minute, t2_second, &
 t3_year, t3_month, t3_day, t3_hour, t3_minute, t3_second


implicit none

! Lookup lengths
integer(kind=int64), parameter :: len_int_lookup = 45
integer(kind=int64), parameter :: len_real_lookup = 19

private

public :: lfricinp_add_um_field_to_file

contains

subroutine lfricinp_add_um_field_to_file(um_file, stashcode, level_index, &
                                         um_grid, lbtim_from_conf, lbproc_from_conf)
! Description:
!  Adds a field object of given stashcode and level index to the um file
!  object. Uses metadata from um_file headers as well as the um_grid
!  object in order to populate the field metadata. Allocates space for
!  field data.
implicit none

type(shum_file_type), intent(in out) :: um_file
type(lfricinp_grid_type), intent(in) :: um_grid

integer(kind=int64), intent(in) :: stashcode
integer(kind=int64), intent(in) :: level_index ! Level index in field array
integer(kind=int64), intent(in) :: lbtim_from_conf  ! 'lbtim' and 'lbproc' metadata
integer(kind=int64), intent(in) :: lbproc_from_conf ! values passed in from the conf file

character(len=*), parameter :: routinename = 'lfricinp_add_um_field_to_file'

integer(kind=int64) :: grid_type_code ! Stashmaster grid type code
integer(kind=int64) :: level_code     ! Stashmaster level code

integer(kind=int64) :: level_number, i_field

integer(kind=int64) :: lookup_int(len_int_lookup) = um_imdi
real(kind=real64)   :: lookup_real_tmp(len_real_lookup+len_int_lookup) = um_rmdi
type(shum_field_type) :: temp_field

! Access file headers
integer(kind=int64) :: fixed_length_header(f_shum_fixed_length_header_len)
integer(kind=int64), allocatable :: integer_constants(:)
real(kind=real64), allocatable :: real_constants(:)
real(kind=real64), allocatable :: level_dep_c(:,:)

! Get a copy of file headers
call shumlib(routinename//'::get_fixed_length_header', &
     um_file % get_fixed_length_header(fixed_length_header))
call shumlib(routinename//'::get_integer_constants', &
     um_file % get_integer_constants(integer_constants))
call shumlib(routinename//'::get_real_constants', &
     um_file % get_real_constants(real_constants))
call shumlib(routinename// &
     '::get_level_dependent_constants', &
     um_file % get_level_dependent_constants(level_dep_c))

! Populate lookup
! Validity Time
lookup_int(lbyr) = fixed_length_header(t2_year)
lookup_int(lbmon) = fixed_length_header(t2_month)
lookup_int(lbdat) = fixed_length_header(t2_day)
lookup_int(lbhr) = fixed_length_header(t2_hour)
lookup_int(lbmin) = fixed_length_header(t2_minute)
lookup_int(lbsec) = fixed_length_header(t2_second)
! Data time
lookup_int(lbyrd) = fixed_length_header(t1_year)
lookup_int(lbmond) = fixed_length_header(t1_month)
lookup_int(lbdatd) = fixed_length_header(t1_day)
lookup_int(lbhrd) = fixed_length_header(t1_hour)
lookup_int(lbmind) = fixed_length_header(t1_minute)
lookup_int(lbsecd) = fixed_length_header(t1_second)
lookup_int(lbtim) = lbtim_from_conf
lookup_int(lbft) = 0 ! Forecast time, difference between datatime and
                     ! validity time

lookup_int(lbcode) = 1 ! harcode to lat/long non-rotated
lookup_int(lbhem) = 0 ! Hardcode to global
! Determine grid type
grid_type_code = get_stashmaster_item(stashcode, grid)

! Set horiz grid
select case(grid_type_code)
case(p_points, ozone_points)
  lookup_int(lbrow) = um_grid%num_p_points_y
  lookup_int(lbnpt) = um_grid%num_p_points_x
  ! "Zeroth" start lat/lon so subtract one grid spacing
  lookup_real_tmp(bzy) = um_grid%p_origin_y - um_grid%spacing_y
  lookup_real_tmp(bzx) = um_grid%p_origin_x - um_grid%spacing_x
case(u_points)
  lookup_int(lbrow) = um_grid%num_u_points_y
  lookup_int(lbnpt) = um_grid%num_u_points_x
  ! "Zeroth" start lat/lon so subtract one grid spacing
  lookup_real_tmp(bzy) = um_grid%u_origin_y - um_grid%spacing_y
  lookup_real_tmp(bzx) = um_grid%u_origin_x - um_grid%spacing_x
case(v_points)
  lookup_int(lbrow) = um_grid%num_v_points_y
  lookup_int(lbnpt) = um_grid%num_v_points_x
  ! "Zeroth" start lat/lon so subtract one grid spacing
  lookup_real_tmp(bzy) = um_grid%v_origin_y - um_grid%spacing_y
  lookup_real_tmp(bzx) = um_grid%v_origin_x - um_grid%spacing_x
case DEFAULT
  write(log_scratch_space, '(A,I0,A)') &
       "Grid type code ", grid_type_code, " not supported"
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end select
lookup_real_tmp(bplat) = real_constants(rh_polelat)
lookup_real_tmp(bplon) = real_constants(rh_polelong)
lookup_real_tmp(bdx) = real_constants(rh_deltaEW)
lookup_real_tmp(bdy) = real_constants(rh_deltaNS)


lookup_int(lbpack) = 2 ! Currently only support 32 bit packing
lookup_int(lbrel) = 3 ! UM version 8.1 onwards
lookup_int(lbfc) = get_stashmaster_item(stashcode, ppfc) ! field code
lookup_int(lbcfc) = 0  ! Always set to 0 for UM files
lookup_int(lbproc) = lbproc_from_conf
lookup_int(lbvc) = get_stashmaster_item(stashcode, sm_lbvc) !vertical coord type
lookup_int(lbrvc) = 0 ! Always zero in UM
lookup_int(lbtyp) = get_stashmaster_item(stashcode, cfff) !Fieldsfile field type

! The level index of the field array is indexed from 1. for fields
! with a "zeroth" level on the surface then level_number = level_index - 1
! but for fields without a "zeroth" level then level_number = level_index
if ( lfricinp_get_first_level_num(stashcode) == 0 )  then
  level_number = level_index - 1
else
  level_number = level_index
end if

level_code = get_stashmaster_item(stashcode, levelt)

if (level_code == single_level) then
  ! Single levels get special code from stashmaster
  lookup_int(lblev) = get_stashmaster_item(stashcode, cfll)
else
  if (level_number == 0) then
    ! Levels on surface set to 9999
    lookup_int(lblev) = 9999
  else
    lookup_int(lblev) = level_number
  end if
end if

! Reserved slots (unused)
lookup_int(lbrsvd1) = 0
lookup_int(lbrsvd2) = 0
lookup_int(lbrsvd3) = 0
! Ensemble number - set to 0 for deterministic
lookup_int(lbrsvd4) = 0

! Model version id - UM identifier is 1111
lookup_int(lbsrce) = fixed_length_header(model_version) * 10000 &
                     + 1111
! Data type
lookup_int(lbuser1) = get_stashmaster_item(stashcode, datat)
! Stashcode
lookup_int(lbuser4) = stashcode
! Internal model number: atmosphere
lookup_int(lbuser7) = 1

! Datum value, always 1
lookup_real_tmp(bdatum) = 0.0_real64

! Reserved for future use
lookup_real_tmp(brsvd3) = 0.0_real64
lookup_real_tmp(brsvd4) = 0.0_real64

! Check vertical coordinate type
if (lookup_int(lbvc) >= 126 .and. lookup_int(lbvc) <= 139 &
     .or. lookup_int(lbvc) == 5) then
  ! Special codes inc single level, set to 0.0
  lookup_real_tmp(blev)=0.0_real64
  lookup_real_tmp(bhlev)=0.0_real64
  lookup_real_tmp(brlev)=0.0_real64
  lookup_real_tmp(bhrlev)=0.0_real64
  lookup_real_tmp(bulev)=0.0_real64
  lookup_real_tmp(bhulev)=0.0_real64
else if (lookup_int(lbvc) == 65) then ! Standard hybrid height levels
  ! height of model level k above mean sea level is
  !       z(i,j,k) = Zsea(k) + C(k)*Zorog(i,j)
  ! bulev,bhulev      zsea,C of upper layer boundary
  ! blev ,bhlev       zsea,C of level
  ! brlev,bhrlev      zsea,C of lower level boundary
  ! The level here can refer to either a theta or rho level, with
  ! layer boundaries defined by surrounding rho or theta levels.
  if (level_code == theta_levels) then ! theta level (& w)

    ! When referencing theta arrays need to add 1 to the level number as level
    ! numbers start at 0 but array start at 1
    if (level_number == integer_constants(ih_model_levels)) then ! top level
      lookup_real_tmp(bulev) =  level_dep_c(level_number+1, ldc_zsea_theta)    &
                                * 2.0_real64 - level_dep_c(level_number, ldc_zsea_rho)
      lookup_real_tmp(bhulev)=   level_dep_c(level_number+1, ldc_c_theta)      &
                                * 2.0_real64 - level_dep_c(level_number,ldc_c_rho)
    else
      lookup_real_tmp(bulev) = level_dep_c(level_number+1, ldc_zsea_rho)
      lookup_real_tmp(bhulev)= level_dep_c(level_number+1, ldc_c_rho)
    end if                             ! top level

    lookup_real_tmp(blev) = level_dep_c(level_number+1, ldc_zsea_theta)
    lookup_real_tmp(bhlev)= level_dep_c(level_number+1, ldc_c_theta)

    if (level_number == 0) then        ! Zeroth level
      lookup_real_tmp(brlev) = 0.0_real64     ! zsea at/below surface
      lookup_real_tmp(bhrlev)= 1.0_real64   ! C    at/below surface
    else
      lookup_real_tmp(brlev) = level_dep_c(level_number, ldc_zsea_rho)
      lookup_real_tmp(bhrlev)=    level_dep_c(level_number, ldc_c_rho)
    end if                              ! bottom level

  else if (level_code == rho_levels) then ! rho level (u,v)

    ! if exner above top level
    if (level_number >  integer_constants(ih_model_levels)) then
      lookup_real_tmp(bulev) = level_dep_c(level_number+1, ldc_zsea_theta)     &
                               * 2.0_real64 - level_dep_c(level_number-1, ldc_zsea_rho)
      lookup_real_tmp(bhulev)= level_dep_c(level_number+1, ldc_c_theta)        &
                               * 2.0_real64 - level_dep_c(level_number-1, ldc_c_rho)
      lookup_real_tmp(blev) = lookup_real_tmp(bulev)
      lookup_real_tmp(bhlev)= lookup_real_tmp(bhulev)
    else
      lookup_real_tmp(bulev) = level_dep_c(level_number+1, ldc_zsea_theta)
      lookup_real_tmp(bhulev)=    level_dep_c(level_number+1, ldc_c_theta)
      lookup_real_tmp(blev)  = level_dep_c(level_number, ldc_zsea_rho)
      lookup_real_tmp(bhlev) =    level_dep_c(level_number, ldc_c_rho)
    end if

    lookup_real_tmp(brlev) = level_dep_c(level_number, ldc_zsea_theta)
    lookup_real_tmp(bhrlev)=    level_dep_c(level_number, ldc_c_theta)
  end if
else
   write(log_scratch_space, '(A,I0,A)') &
       "Vertical coord type ", lookup_int(lbvc), " not supported"
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

! Missing data indicator
lookup_real_tmp(bmdi) = um_rmdi
! MKS scaling factor (unity as model uses SI units throughout)
lookup_real_tmp(bmks) = 1.0_real64


! Set lookup
call shumlib(routinename//'::set_lookup',       &
             temp_field%set_lookup(lookup_int, &
             ! Lookup real is dimensioned using full lookup size so that
             ! index parameters match those use in shumlib
             lookup_real_tmp(len_int_lookup+1 :)) )

! Add to file
call shumlib(routinename//'::add_field', um_file%add_field(temp_field))
! Current field index will last one added to file
i_field = um_file%num_fields
! Allocate data array using sizes specified in lookup
select case (get_stashmaster_item(stashcode, datat)) ! Which datatype?
case(1) ! Real
  allocate(um_file%fields(i_field)%rdata(lookup_int(lbnpt), &
                                         lookup_int(lbrow)))
case(2) ! Integer
  allocate(um_file%fields(i_field)%idata(lookup_int(lbnpt), &
                                         lookup_int(lbrow)))
case DEFAULT  ! Logical and others
  write(log_scratch_space, '(A,I0,A)') "Datatype: ", &
       get_stashmaster_item(stashcode, datat), " not supported"
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end select


end subroutine lfricinp_add_um_field_to_file

end module lfricinp_add_um_field_to_file_mod
