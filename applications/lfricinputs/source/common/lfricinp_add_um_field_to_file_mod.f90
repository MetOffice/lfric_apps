! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_add_um_field_to_file_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY: real64, int64

! lfricinp modules
USE lfricinp_stashmaster_mod, ONLY: get_stashmaster_item, grid, &
    p_points, u_points, v_points, ozone_points, ppfc, sm_lbvc => lbvc, cfff,&
    levelt, rho_levels, theta_levels, single_level, cfll, datat
USE lfricinp_um_parameters_mod, ONLY: um_imdi, um_rmdi, rh_polelat, &
    rh_polelong, ldc_zsea_theta, ldc_zsea_rho, ldc_c_theta, ldc_c_rho, &
    rh_deltaEW, rh_deltaNS, ih_model_levels
USE lfricinp_grid_type_mod, ONLY: lfricinp_grid_type
USE lfricinp_um_level_codes_mod, ONLY: lfricinp_get_first_level_num
USE lfricinp_check_shumlib_status_mod, ONLY: shumlib


! lfric modules
USE log_mod, ONLY : log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR, &
     log_scratch_space

! shumlib modules
USE f_shum_fieldsfile_mod, ONLY: f_shum_fixed_length_header_len
USE f_shum_file_mod, ONLY: shum_file_type
USE f_shum_field_mod, ONLY: shum_field_type
USE f_shum_lookup_indices_mod, ONLY: &
    lbyr, lbmon, lbdat, lbhr, lbmin, lbday, lbsec, lbyrd, lbmond, lbdatd, &
    lbhrd, lbmind, lbdayd, lbsecd, lbtim, lbft, lbcode, lbhem,    &
    lbrow, lbnpt, lbpack, lbrel, lbfc, lbcfc, lbproc, lbvc, lbrvc, &
    lbtyp, lblev, lbrsvd1, lbrsvd2,        &
    lbrsvd3, lbrsvd4, lbsrce, lbuser1, lbuser4, &
    lbuser7, bulev, bhulev, brsvd3, brsvd4, bdatum, bacc, blev,  &
    brlev, bhlev, bhrlev, bplat, bplon, bgor, bzy, bdy, bzx, bdx, bmdi,   &
    bmks

USE f_shum_fixed_length_header_indices_mod, ONLY:                        &
 vert_coord_type, horiz_grid_type, &
 dataset_type, run_identifier, calendar, projection_number, model_version, &
 grid_staggering, sub_model, t1_year, t1_month, t1_day, t1_hour, t1_minute, &
 t1_second, t2_year, t2_month, t2_day, t2_hour, t2_minute, t2_second, &
 t3_year, t3_month, t3_day, t3_hour, t3_minute, t3_second


IMPLICIT NONE

! Lookup lengths
INTEGER(KIND=int64), PARAMETER :: len_int_lookup = 45
INTEGER(KIND=int64), PARAMETER :: len_real_lookup = 19

PRIVATE

PUBLIC :: lfricinp_add_um_field_to_file

CONTAINS

SUBROUTINE lfricinp_add_um_field_to_file(um_file, stashcode, level_index, &
                                         um_grid)
! Description:
!  Adds a field object of given stashcode and level index to the um file
!  object. Uses metadata from um_file headers as well as the um_grid
!  object in order to populate the field metadata. Allocates space for
!  field data.
IMPLICIT NONE

TYPE(shum_file_type), INTENT(IN OUT) :: um_file
TYPE(lfricinp_grid_type), INTENT(IN) :: um_grid

INTEGER(KIND=int64), INTENT(IN) :: stashcode
INTEGER(KIND=int64), INTENT(IN) :: level_index ! Level index in field array

CHARACTER(LEN=*), PARAMETER :: routinename = 'lfricinp_add_um_field_to_file'

INTEGER(KIND=int64) :: grid_type_code ! Stashmaster grid type code
INTEGER(KIND=int64) :: level_code     ! Stashmaster level code

INTEGER(KIND=int64) :: level_number, i_field

INTEGER(KIND=int64) :: lookup_int(len_int_lookup) = um_imdi
REAL(KIND=real64)   :: lookup_real_tmp(len_real_lookup+len_int_lookup) = um_rmdi
TYPE(shum_field_type) :: temp_field

! Access file headers
INTEGER(KIND=int64) :: fixed_length_header(f_shum_fixed_length_header_len)
INTEGER(KIND=int64), ALLOCATABLE :: integer_constants(:)
REAL(KIND=real64), ALLOCATABLE :: real_constants(:)
REAL(KIND=real64), ALLOCATABLE :: level_dep_c(:,:)

! Get a copy of file headers
CALL shumlib(routinename//'::get_fixed_length_header', &
     um_file % get_fixed_length_header(fixed_length_header))
CALL shumlib(routinename//'::get_integer_constants', &
     um_file % get_integer_constants(integer_constants))
CALL shumlib(routinename//'::get_real_constants', &
     um_file % get_real_constants(real_constants))
CALL shumlib(routinename// &
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
lookup_int(lbtim) = 1  ! hardcode to proleptic gregorian
lookup_int(lbft) = 0 ! Forecast time, difference between datatime and
                     ! validity time

lookup_int(lbcode) = 1 ! harcode to lat/long non-rotated
lookup_int(lbhem) = 0 ! Hardcode to global
! Determine grid type
grid_type_code = get_stashmaster_item(stashcode, grid)

! Set horiz grid
SELECT CASE(grid_type_code)
CASE(p_points, ozone_points)
  lookup_int(lbrow) = um_grid%num_p_points_y
  lookup_int(lbnpt) = um_grid%num_p_points_x
  ! "Zeroth" start lat/lon so subtract one grid spacing
  lookup_real_tmp(bzy) = um_grid%p_origin_y - um_grid%spacing_y
  lookup_real_tmp(bzx) = um_grid%p_origin_x - um_grid%spacing_x
CASE(u_points)
  lookup_int(lbrow) = um_grid%num_u_points_y
  lookup_int(lbnpt) = um_grid%num_u_points_x
  ! "Zeroth" start lat/lon so subtract one grid spacing
  lookup_real_tmp(bzy) = um_grid%u_origin_y - um_grid%spacing_y
  lookup_real_tmp(bzx) = um_grid%u_origin_x - um_grid%spacing_x
CASE(v_points)
  lookup_int(lbrow) = um_grid%num_v_points_y
  lookup_int(lbnpt) = um_grid%num_v_points_x
  ! "Zeroth" start lat/lon so subtract one grid spacing
  lookup_real_tmp(bzy) = um_grid%v_origin_y - um_grid%spacing_y
  lookup_real_tmp(bzx) = um_grid%v_origin_x - um_grid%spacing_x
CASE DEFAULT
  WRITE(log_scratch_space, '(A,I0,A)') &
       "Grid type code ", grid_type_code, " not supported"
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END SELECT
lookup_real_tmp(bplat) = real_constants(rh_polelat)
lookup_real_tmp(bplon) = real_constants(rh_polelong)
lookup_real_tmp(bdx) = real_constants(rh_deltaEW)
lookup_real_tmp(bdy) = real_constants(rh_deltaNS)


lookup_int(lbpack) = 2 ! Currently only support 32 bit packing
lookup_int(lbrel) = 3 ! UM version 8.1 onwards
lookup_int(lbfc) = get_stashmaster_item(stashcode, ppfc) ! field code
lookup_int(lbcfc) = 0  ! Always set to 0 for UM files
lookup_int(lbproc) = 0 ! No field processing done
lookup_int(lbvc) = get_stashmaster_item(stashcode, sm_lbvc) !vertical coord type
lookup_int(lbrvc) = 0 ! Always zero in UM
lookup_int(lbtyp) = get_stashmaster_item(stashcode, cfff) !Fieldsfile field type

! The level index of the field array is indexed from 1. for fields
! with a "zeroth" level on the surface then level_number = level_index - 1
! but for fields without a "zeroth" level then level_number = level_index
IF ( lfricinp_get_first_level_num(stashcode) == 0 )  THEN
  level_number = level_index - 1
ELSE
  level_number = level_index
END IF

level_code = get_stashmaster_item(stashcode, levelt)

IF (level_code == single_level) THEN
  ! Single levels get special code from stashmaster
  lookup_int(lblev) = get_stashmaster_item(stashcode, cfll)
ELSE
  IF (level_number == 0) THEN
    ! Levels on surface set to 9999
    lookup_int(lblev) = 9999
  ELSE
    lookup_int(lblev) = level_number
  END IF
END IF

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
lookup_real_tmp(bdatum) = 0.0

! Reserved for future use
lookup_real_tmp(brsvd3) = 0.0
lookup_real_tmp(brsvd4) = 0.0

! Check vertical coordinate type
IF (lookup_int(lbvc) >= 126 .AND. lookup_int(lbvc) <= 139 &
     .OR. lookup_int(lbvc) == 5) THEN
  ! Special codes inc single level, set to 0.0
  lookup_real_tmp(blev)=0.0
  lookup_real_tmp(bhlev)=0.0
  lookup_real_tmp(brlev)=0.0
  lookup_real_tmp(bhrlev)=0.0
  lookup_real_tmp(bulev)=0.0
  lookup_real_tmp(bhulev)=0.0
ELSE IF (lookup_int(lbvc) == 65) THEN ! Standard hybrid height levels
  ! height of model level k above mean sea level is
  !       z(i,j,k) = Zsea(k) + C(k)*Zorog(i,j)
  ! bulev,bhulev      zsea,C of upper layer boundary
  ! blev ,bhlev       zsea,C of level
  ! brlev,bhrlev      zsea,C of lower level boundary
  ! The level here can refer to either a theta or rho level, with
  ! layer boundaries defined by surrounding rho or theta levels.
  IF (level_code == theta_levels) THEN ! theta level (& w)

    ! When referencing theta arrays need to add 1 to the level number as level
    ! numbers start at 0 but array start at 1
    IF (level_number == integer_constants(ih_model_levels)) THEN ! top level
      lookup_real_tmp(bulev) =  level_dep_c(level_number+1, ldc_zsea_theta)    &
                                * 2.0 - level_dep_c(level_number, ldc_zsea_rho)
      lookup_real_tmp(bhulev)=   level_dep_c(level_number+1, ldc_c_theta)      &
                                * 2.0 - level_dep_c(level_number,ldc_c_rho)
    ELSE
      lookup_real_tmp(bulev) = level_dep_c(level_number+1, ldc_zsea_rho)
      lookup_real_tmp(bhulev)= level_dep_c(level_number+1, ldc_c_rho)
    END IF                             ! top level

    lookup_real_tmp(blev) = level_dep_c(level_number+1, ldc_zsea_theta)
    lookup_real_tmp(bhlev)= level_dep_c(level_number+1, ldc_c_theta)

    IF (level_number == 0) THEN        ! Zeroth level
      lookup_real_tmp(brlev) = 0.0     ! zsea at/below surface
      lookup_real_tmp(bhrlev)= 1.0     ! C    at/below surface
    ELSE
      lookup_real_tmp(brlev) = level_dep_c(level_number, ldc_zsea_rho)
      lookup_real_tmp(bhrlev)=    level_dep_c(level_number, ldc_c_rho)
    END IF                              ! bottom level

  ELSE IF (level_code == rho_levels) THEN ! rho level (u,v)

    ! if exner above top level
    IF (level_number >  integer_constants(ih_model_levels)) THEN
      lookup_real_tmp(bulev) = level_dep_c(level_number+1, ldc_zsea_theta)     &
                               * 2.0 - level_dep_c(level_number-1, ldc_zsea_rho)
      lookup_real_tmp(bhulev)= level_dep_c(level_number+1, ldc_c_theta)        &
                               * 2.0 - level_dep_c(level_number-1, ldc_c_rho)
      lookup_real_tmp(blev) = lookup_real_tmp(bulev)
      lookup_real_tmp(bhlev)= lookup_real_tmp(bhulev)
    ELSE
      lookup_real_tmp(bulev) = level_dep_c(level_number+1, ldc_zsea_theta)
      lookup_real_tmp(bhulev)=    level_dep_c(level_number+1, ldc_c_theta)
      lookup_real_tmp(blev)  = level_dep_c(level_number, ldc_zsea_rho)
      lookup_real_tmp(bhlev) =    level_dep_c(level_number, ldc_c_rho)
    END IF

    lookup_real_tmp(brlev) = level_dep_c(level_number, ldc_zsea_theta)
    lookup_real_tmp(bhrlev)=    level_dep_c(level_number, ldc_c_theta)
  END IF
ELSE
   WRITE(log_scratch_space, '(A,I0,A)') &
       "Vertical coord type ", lookup_int(lbvc), " not supported"
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

! Missing data indicator
lookup_real_tmp(bmdi) = um_rmdi
! MKS scaling factor (unity as model uses SI units throughout)
lookup_real_tmp(bmks) = 1.0_real64


! Set lookup
CALL shumlib(routinename//'::set_lookup',       &
             temp_field%set_lookup(lookup_int, &
             ! Lookup real is dimensioned using full lookup size so that
             ! index parameters match those use in shumlib
             lookup_real_tmp(len_int_lookup+1 :)) )

! Add to file
CALL shumlib(routinename//'::add_field', um_file%add_field(temp_field))
! Current field index will last one added to file
i_field = um_file%num_fields
! Allocate data array using sizes specified in lookup
SELECT CASE (get_stashmaster_item(stashcode, datat)) ! Which datatype?
CASE(1) ! Real
  ALLOCATE(um_file%fields(i_field)%rdata(lookup_int(lbnpt), &
                                         lookup_int(lbrow)))
CASE(2) ! Integer
  ALLOCATE(um_file%fields(i_field)%idata(lookup_int(lbnpt), &
                                         lookup_int(lbrow)))
CASE DEFAULT  ! Logical and others
  WRITE(log_scratch_space, '(A,I0,A)') "Datatype: ", &
       get_stashmaster_item(stashcode, datat), " not supported"
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END SELECT


END SUBROUTINE lfricinp_add_um_field_to_file

END MODULE lfricinp_add_um_field_to_file_mod
