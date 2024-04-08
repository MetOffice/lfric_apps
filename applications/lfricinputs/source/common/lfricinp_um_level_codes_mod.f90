! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_um_level_codes_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY: int64

! Shumlib modules
USE f_shum_file_mod, ONLY: shum_file_type

! lfricinputs modules
USE lfricinp_check_shumlib_status_mod, ONLY: shumlib
USE lfricinp_grid_type_mod, ONLY: lfricinp_grid_type

! lfric modules
USE log_mod,  ONLY : log_event, log_scratch_space,         &
                     LOG_LEVEL_ERROR

IMPLICIT NONE

PRIVATE

PUBLIC :: lfricinp_get_num_levels, lfricinp_get_first_level_num, &
     lfricinp_get_last_level_num, lfricinp_get_num_pseudo_levels


CONTAINS

!------------------------------------------------------------------

FUNCTION lfricinp_get_num_levels(um_file, stashcode) RESULT(num_levels)
! Description:
!  Returns the number of levels expected in a field as defined by
!  first and last level codes in the stashmaster

IMPLICIT NONE

TYPE(shum_file_type), INTENT(IN) :: um_file
INTEGER(KIND=int64), INTENT(IN) :: stashcode
! Result
INTEGER(KIND=int64) :: num_levels

num_levels = lfricinp_get_last_level_num(um_file, stashcode) - &
             lfricinp_get_first_level_num(stashcode) + 1

END FUNCTION lfricinp_get_num_levels

!------------------------------------------------------------------

FUNCTION lfricinp_get_first_level_num(stashcode) RESULT(first_level_num)
! Description:
!  Takes stashcode as input and interogates the stashmaster first
!  level code to determine the first/bottom level number for the field

! lfricinp modules
USE lfricinp_stashmaster_mod, ONLY: get_stashmaster_item, levelf

IMPLICIT NONE
INTEGER(KIND=int64), INTENT(IN) :: stashcode
! Result
INTEGER(KIND=int64) :: first_level_num
! Local variables
INTEGER(KIND=int64) :: first_level_code

first_level_code =  get_stashmaster_item(stashcode, levelf)

SELECT CASE(first_level_code)
  CASE (-1) ! Unset /single level
    first_level_num = 1
  CASE (1) ! First atmos level
    first_level_num = 1
  CASE (8) ! First soil level
    first_level_num = 1
  CASE (38,40) ! Zeroth atmos level
    first_level_num = 0
CASE DEFAULT
  WRITE(log_scratch_space, '(A,I0,A)') &
     "First level code ", first_level_code, " not supported"
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END SELECT

END FUNCTION lfricinp_get_first_level_num

!------------------------------------------------------------------

FUNCTION lfricinp_get_last_level_num(um_file, stashcode) RESULT(last_level_num)
! Description:
!  Takes um_file and stashcode as input and interogates the stashmaster
!  last level code to determine the last/top level number for the field

! lfricinp modules
USE lfricinp_stashmaster_mod, ONLY: get_stashmaster_item, levell
USE lfricinp_um_parameters_mod, ONLY: ih_model_levels, ih_soilT_levels
IMPLICIT NONE
TYPE(shum_file_type), INTENT(IN) :: um_file
INTEGER(KIND=int64), INTENT(IN) :: stashcode
! Result
INTEGER(KIND=int64) :: last_level_num
! Local variables
INTEGER(KIND=int64) :: last_level_code
INTEGER(KIND=int64) :: model_levels
CHARACTER(LEN=*), PARAMETER :: routinename='lfricinp_get_last_level_num'


last_level_code =  get_stashmaster_item(stashcode, levell)
! Get model_levels
CALL shumlib(routinename //'::get_integer_constants_by_index', &
     um_file%get_integer_constants_by_index(ih_model_levels,   &
     model_levels))

SELECT CASE(last_level_code)
  CASE (-1) ! Unset /single level
    last_level_num = 1
  CASE (2,3) ! Top atmos level
    last_level_num = model_levels
  CASE (9) ! Last soil level
    CALL shumlib(routinename //'::get_integer_constants_by_index', &
     um_file%get_integer_constants_by_index(ih_soilT_levels,   &
     last_level_num))
  CASE(23) ! Top ozone level
    last_level_num = model_levels
CASE DEFAULT
  WRITE(log_scratch_space, '(A,I0,A)') &
     "Last level code ", last_level_code, " not supported"
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END SELECT

END FUNCTION lfricinp_get_last_level_num

!------------------------------------------------------------------

FUNCTION lfricinp_get_num_pseudo_levels(um_grid, stashcode) &
     RESULT(num_pseudo_levels)
! Description:
!  Returns the number of pseudo levels expected in a field as defined by
!  first and last pseudo level codes in the stashmaster

IMPLICIT NONE

TYPE(lfricinp_grid_type), INTENT(IN):: um_grid
INTEGER(KIND=int64), INTENT(IN) :: stashcode
! Result
INTEGER(KIND=int64) :: num_pseudo_levels

num_pseudo_levels = lfricinp_get_last_pseudo_level_num(um_grid, stashcode) - &
             lfricinp_get_first_pseudo_level_num(stashcode) + 1

END FUNCTION lfricinp_get_num_pseudo_levels

!------------------------------------------------------------------

FUNCTION lfricinp_get_first_pseudo_level_num(stashcode) &
     RESULT(first_pseudo_level_num)
! Description:
!  Takes stashcode as input and interogates the stashmaster first
!  pseudo_level code to determine the first/bottom pseudo_level
!  number for the field

! lfricinp modules
USE lfricinp_stashmaster_mod, ONLY: get_stashmaster_item, pseudf

IMPLICIT NONE
INTEGER(KIND=int64), INTENT(IN) :: stashcode
! Result
INTEGER(KIND=int64) :: first_pseudo_level_num
! Local variables
INTEGER(KIND=int64) :: first_pseudo_level_code

first_pseudo_level_code =  get_stashmaster_item(stashcode, pseudf)

SELECT CASE(first_pseudo_level_code)
  CASE (1) ! Dimension starts at 1
    first_pseudo_level_num = 1
CASE DEFAULT
  WRITE(log_scratch_space, '(A,I0,A)') &
     "First pseudo_level code ", first_pseudo_level_code, " not supported"
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END SELECT

END FUNCTION lfricinp_get_first_pseudo_level_num

!------------------------------------------------------------------

FUNCTION lfricinp_get_last_pseudo_level_num(um_grid, stashcode) &
     RESULT(last_pseudo_level_num)
! Description:
!  Takes um_grid and stashcode as input and interogates the stashmaster
!  last pseudo_level code to determine the last/top pseudo_level number
!  for the field

! lfricinp modules
USE lfricinp_stashmaster_mod, ONLY: get_stashmaster_item, pseudl

IMPLICIT NONE
TYPE(lfricinp_grid_type), INTENT(IN) :: um_grid
INTEGER(KIND=int64), INTENT(IN) :: stashcode
! Result
INTEGER(KIND=int64) :: last_pseudo_level_num
! Local variables
INTEGER(KIND=int64) :: last_pseudo_level_code
CHARACTER(LEN=*), PARAMETER :: routinename = &
     'lfricinp_get_last_pseudo_level_num'

last_pseudo_level_code =  get_stashmaster_item(stashcode, pseudl)

SELECT CASE(last_pseudo_level_code)
CASE (7,9) ! ntypes == ntiles (lfricinputs doesn't support aggregate tile)
  last_pseudo_level_num = um_grid % num_surface_types
CASE (11)
  last_pseudo_level_num = um_grid%num_snow_layers * um_grid%num_surface_types
CASE DEFAULT
  WRITE(log_scratch_space, '(A,I0,A)') &
       "Last pseudo_level code ", last_pseudo_level_code, " not supported"
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END SELECT

END FUNCTION lfricinp_get_last_pseudo_level_num


END MODULE lfricinp_um_level_codes_mod
