! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_create_lfric_fields_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int64

! lfric modules
USE constants_mod,                 ONLY : i_def, str_def
USE field_mod,                     ONLY : field_type
USE field_parent_mod,              ONLY :  read_interface, write_interface,    &
                                          field_parent_proxy_type
USE lfric_xios_time_axis_mod,      ONLY : time_axis_type
USE lfric_xios_read_mod,           ONLY : read_field_generic
USE lfric_xios_write_mod,          ONLY : write_field_generic
USE finite_element_config_mod,     ONLY : element_order
USE function_space_mod,            ONLY : function_space_type
USE function_space_collection_mod, ONLY : function_space_collection
USE field_collection_mod,          ONLY : field_collection_type
USE fs_continuity_mod,             ONLY : W3, Wtheta, W2H
USE log_mod,                       ONLY : LOG_LEVEL_INFO, LOG_LEVEL_ERROR,     &
                                          log_event, log_scratch_space
USE mesh_mod,                      ONLY : mesh_type

! lfricinp modules
USE lfricinp_stashmaster_mod,        ONLY: get_stashmaster_item, pseudt
USE lfricinp_stash_to_lfric_map_mod, ONLY: w2h_field, w3_field, w3_field_2d,   &
                                           w3_soil_field, wtheta_field,        &
                                           get_field_name, get_lfric_field_kind
USE lfricinp_grid_type_mod,          ONLY: lfricinp_grid_type
USE lfricinp_um_level_codes_mod,     ONLY: lfricinp_get_num_levels,            &
                                           lfricinp_get_num_pseudo_levels
USE lfricinp_regrid_options_mod,     ONLY: winds_on_w3

! shumlib modules
USE f_shum_file_mod, ONLY: shum_file_type

IMPLICIT NONE
PRIVATE
PUBLIC :: lfricinp_create_lfric_fields

CONTAINS

SUBROUTINE lfricinp_create_lfric_fields( mesh, twod_mesh,              &
                                         field_collection, stash_list, &
                                         um_grid, um_file)
! Description:
!  When given list of stashcodes create lfric field. Use stashcode to determine
!  which vertical level set the field is on and create appropiate function
!  space. Set read and write procedure for each field.


IMPLICIT NONE

TYPE(mesh_type), INTENT(IN), POINTER :: mesh
TYPE(mesh_type), INTENT(IN), POINTER :: twod_mesh

TYPE(field_collection_type), INTENT(IN OUT) :: field_collection

INTEGER(KIND=int64), INTENT(IN)  :: stash_list(:)

TYPE(lfricinp_grid_type), INTENT(IN):: um_grid

TYPE(shum_file_type), INTENT(IN) :: um_file

PROCEDURE(read_interface),  POINTER  :: tmp_read_ptr => NULL()
PROCEDURE(write_interface), POINTER  :: tmp_write_ptr => NULL()
TYPE(function_space_type),  POINTER :: vector_space => NULL()

TYPE(mesh_type), POINTER :: type_mesh => null()

TYPE( field_type )            :: field
INTEGER(KIND=int64)           :: stashcode, lfric_field_kind, i
INTEGER(KIND=i_def)           :: fs_id
! A temporary variable to hold the 64bit output
INTEGER(KIND=int64)           :: ndata_64
CHARACTER(LEN=:), ALLOCATABLE :: field_name
LOGICAL                       :: ndata_first

CALL log_event('Creating lfric finite difference fields ...', LOG_LEVEL_INFO)

IF (element_order > 0) CALL log_event('Finite difference fields requires ' //  &
                                      'lowest order elements', LOG_LEVEL_ERROR)

! Create the field collection
CALL field_collection%initialise(name="lfric_fields", table_len=100)

DO i=1, SIZE(stash_list)

  stashcode = stash_list(i)
  lfric_field_kind = get_lfric_field_kind(stashcode)
  field_name = TRIM(get_field_name(stashcode))

  ! TODO - It would be good to move this case statement into its own
  ! routine it could be used elsewhere to add fields to the collection
  ! without the do loop over stash list
  IF ( .NOT. field_collection%field_exists(field_name) ) THEN

    tmp_read_ptr => read_field_generic
    tmp_write_ptr => write_field_generic

    SELECT CASE (lfric_field_kind)

      CASE(w2h_field) ! Stashcodes that would map to W2h, i.e. winds
        type_mesh => mesh
        fs_id = W2H
        ndata_64 = 1_int64
        ndata_first = .FALSE.

      CASE(w3_field) ! Stashcodes that map to W3/rho
        type_mesh => mesh
        fs_id = W3
        ndata_64 = 1_int64
        ndata_first = .FALSE.

      CASE(wtheta_field) ! Stashcodes that maps to Wtheta
        type_mesh => mesh
        fs_id = Wtheta
        ndata_64 = 1_int64
        ndata_first = .FALSE.

      CASE(w3_field_2d) ! Stash that needs 2D mesh
        type_mesh => twod_mesh
        fs_id = W3
        IF ( get_stashmaster_item(stashcode, pseudt) == 0 ) THEN
          ! Field has no pseudo levels
          ndata_64 = 1_int64
        ELSE
          ! Get number of pseudo levels/ndata
          ndata_64 = lfricinp_get_num_pseudo_levels(um_grid, stashcode)
        END IF
        ndata_first = .FALSE.

      CASE(w3_soil_field) ! Soil fields
        type_mesh => twod_mesh
        fs_id = W3
        ndata_64 = lfricinp_get_num_levels(um_file, stashcode)
        ndata_first = .TRUE.

      CASE DEFAULT
        WRITE(log_scratch_space, '(A,I0,A)')                                   &
           "LFRic field kind code ", lfric_field_kind, " not recognised"
        CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)

    END SELECT

    vector_space => function_space_collection%get_fs(                          &
                                              type_mesh,                       &
                                              element_order,                   &
                                              fs_id,                           &
                                              ndata=INT(ndata_64, KIND=i_def), &
                                              ndata_first=ndata_first          &
                                                    )
    CALL field % initialise(vector_space=vector_space, name=field_name, &
                            halo_depth = type_mesh%get_halo_depth() )
    CALL field % set_read_behaviour(tmp_read_ptr)
    CALL field % set_write_behaviour(tmp_write_ptr)
    CALL log_event("Add "//field_name//" to field collection", LOG_LEVEL_INFO)
    CALL field_collection % add_field(field)

    NULLIFY(tmp_read_ptr, tmp_write_ptr)
    NULLIFY(vector_space)

  END IF

END DO

CALL log_event('... Done creating finite difference fields', LOG_LEVEL_INFO)

END SUBROUTINE lfricinp_create_lfric_fields

END MODULE lfricinp_create_lfric_fields_mod
