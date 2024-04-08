! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_masks_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY: int32, int64, real64
USE, INTRINSIC :: ISO_C_BINDING, ONLY: C_BOOL

IMPLICIT NONE

LOGICAL, ALLOCATABLE :: um_land_mask(:,:), um_maritime_mask(:,:)
LOGICAL, ALLOCATABLE :: lfric_land_mask(:), lfric_maritime_mask(:)

CONTAINS

SUBROUTINE lfricinp_init_masks(stashcode_land_mask)

USE log_mod,        ONLY: log_event, log_scratch_space, LOG_LEVEL_INFO
USE field_mod,      ONLY: lfric_field_type => field_type,                      &
                          lfric_proxy_type => field_proxy_type
USE local_mesh_mod, ONLY: local_mesh_type
USE mesh_mod,       ONLY: mesh_type

USE lfricinp_lfric_driver_mod,         ONLY: local_rank
USE lfricinp_check_shumlib_status_mod, ONLY: shumlib
USE um2lfric_read_um_file_mod,        ONLY: um_input_file
USE lfricinp_ancils_mod,               ONLY: ancil_fields, l_land_area_fraction
!
! shumlib modules
USE f_shum_field_mod, ONLY: shum_field_type

IMPLICIT NONE

! LFRic mesh and local mesh
TYPE(mesh_type), POINTER :: mesh
TYPE(local_mesh_type), POINTER :: local_mesh

! Ancil fields
TYPE(lfric_field_type), POINTER :: ancil_field
TYPE(lfric_proxy_type) :: ancil_field_proxy

! Array of shumlib field objects that will be returned from UM file
TYPE(shum_field_type), ALLOCATABLE  :: um_input_fields(:)

! Local variables
INTEGER(KIND=int64)  :: stashcode_land_mask
! Results of SIZE() are the default integer kind
INTEGER :: dim_1d, dim_2dx, dim_2dy, i
INTEGER(KIND=int32)  :: err
LOGICAL(KIND=C_BOOL) :: true_cbool

true_cbool = LOGICAL(.TRUE., KIND=C_BOOL)

! Get UM land mask
CALL shumlib("um2lfric::find_fields_in_file",                                  &
             um_input_file%find_fields_in_file(um_input_fields,                &
                                               stashcode = stashcode_land_mask,&
                                               lbproc = 0_int64),              &
             ignore_warning=true_cbool, errorstatus=err)

! Check that both the source and destination land area fraction fields are
! available. Source from start dump via shumlib and destination from lfric ancil
IF (err /= 0 .OR. .NOT. l_land_area_fraction) THEN

  IF (err /=0 ) THEN
    log_scratch_space = 'Error encountered when trying to read UM land mask. '
  ELSE IF (.NOT. l_land_area_fraction) THEN
    log_scratch_space = 'l_land_area_fraction set to FALSE. Land fraction '//  &
         'ancil not available. '
  END IF
  log_scratch_space = TRIM(log_scratch_space) //                               &
                      'Post regridding adjustments '//                         &
                      'on land and maritime compressed fields will be ignored.'
  CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

ELSE

  ! Set up LFRic mask dimension size
  call ancil_fields%get_field("land_area_fraction", ancil_field)
  CALL ancil_field%read_field(ancil_field%get_name())
  ancil_field_proxy = ancil_field%get_proxy()
  dim_1d = SIZE(ancil_field_proxy%data)
  !
  ! Get local_mesh LFRic mask is defined on
  mesh => ancil_field%get_mesh()
  local_mesh => mesh%get_local_mesh()
  !
  ! Set up LFRic logical land mask
  ALLOCATE(lfric_land_mask(dim_1d))
  DO i = 1, dim_1d
    lfric_land_mask(i) = .FALSE.
    IF (local_mesh%get_cell_owner(i) == local_rank ) THEN
      IF (ancil_field_proxy%data(i) > 0.0) THEN
        lfric_land_mask(i) = .TRUE.
      END IF
    END IF
  END DO
  !
  ! Set up LFRic logical maritime mask
  ALLOCATE(lfric_maritime_mask(dim_1d))
  DO i = 1, dim_1d
    lfric_maritime_mask(i) = .FALSE.
    IF (local_mesh%get_cell_owner(i) == local_rank ) THEN
      IF (ancil_field_proxy%data(i) < 1.0) THEN
        lfric_maritime_mask(i) = .TRUE.
      END IF
    END IF
  END DO
  !
  ! Nullify LFRic mesh and local_mesh pointers
  NULLIFY(mesh, local_mesh)

  ! Set up UM mask dimensions
  dim_2dx = SIZE(um_input_fields(1)%rdata,DIM=1)
  dim_2dy = SIZE(um_input_fields(1)%rdata,DIM=2)
  !
  ! Set up UM logical land mask
  ALLOCATE(um_land_mask(dim_2dx,dim_2dy))
  um_land_mask = (um_input_fields(1)%rdata > 0.0)
  !
  ! Set up UM logical maritime mask
  ALLOCATE(um_maritime_mask(dim_2dx,dim_2dy))
  um_maritime_mask = (um_input_fields(1)%rdata < 1.0)

END IF

END SUBROUTINE lfricinp_init_masks


SUBROUTINE lfricinp_finalise_masks()

IMPLICIT NONE

IF (ALLOCATED(um_land_mask)) DEALLOCATE(um_land_mask)
IF (ALLOCATED(lfric_land_mask)) DEALLOCATE(lfric_land_mask)
IF (ALLOCATED(um_maritime_mask)) DEALLOCATE(um_maritime_mask)
IF (ALLOCATED(lfric_maritime_mask)) DEALLOCATE(lfric_maritime_mask)

END SUBROUTINE lfricinp_finalise_masks

END MODULE lfricinp_masks_mod
