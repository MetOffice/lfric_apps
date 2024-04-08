! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_ancils_mod

USE constants_mod,                  ONLY : i_def, l_def
USE log_mod,                        ONLY : log_event,         &
                                           log_scratch_space, &
                                           LOG_LEVEL_INFO
USE field_collection_mod,           ONLY: field_collection_type
USE field_mod,                      ONLY : field_type
USE field_parent_mod,               ONLY : read_interface, &
                                           write_interface
use lfric_xios_read_mod,            ONLY : read_field_generic
use lfric_xios_write_mod,           ONLY : write_field_generic
USE field_collection_mod,           ONLY : field_collection_type
USE mesh_mod,                       ONLY : mesh_type
USE function_space_mod,             ONLY : function_space_type
USE function_space_collection_mod,  ONLY : function_space_collection
USE fs_continuity_mod,              ONLY : W3

IMPLICIT NONE

LOGICAL :: l_land_area_fraction = .FALSE.
! Container for ancil fields
TYPE(field_collection_type) :: ancil_fields

PUBLIC   :: lfricinp_create_ancil_fields,         &
            lfricinp_setup_ancil_field,           &
            l_land_area_fraction

CONTAINS

! Organises fields to be read from ancils into ancil_fields collection

SUBROUTINE lfricinp_create_ancil_fields( ancil_fields, mesh, twod_mesh )

IMPLICIT NONE

TYPE( field_collection_type ), INTENT( OUT )         :: ancil_fields
TYPE( mesh_type ),             INTENT( IN ), POINTER :: mesh
TYPE( mesh_type ),             INTENT( IN ), POINTER :: twod_mesh

! Set up ancil_fields collection
WRITE(log_scratch_space,'(A,A)') "Create ancil fields: "// &
     "Setting up ancil field collection"
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
CALL ancil_fields%initialise(name='ancil_fields', table_len=100)

IF (l_land_area_fraction) THEN
  ! Surface ancils
  CALL log_event("Create land area fraction ancil", LOG_LEVEL_INFO)
  CALL lfricinp_setup_ancil_field("land_area_fraction", ancil_fields, mesh, &
       twod_mesh, twod=.TRUE.)
END IF

END SUBROUTINE lfricinp_create_ancil_fields

!------------------------------------------------------------

! Creates fields to be read into from ancillary files

SUBROUTINE lfricinp_setup_ancil_field( name, ancil_fields, mesh, &
                                       twod_mesh, twod, ndata )

IMPLICIT NONE

CHARACTER(*), INTENT(IN)                    :: name
TYPE(field_collection_type), INTENT(IN OUT) :: ancil_fields
TYPE( mesh_type ), INTENT(IN), POINTER      :: mesh
TYPE( mesh_type ), INTENT(IN), POINTER      :: twod_mesh
LOGICAL(l_def), OPTIONAL, INTENT(IN)        :: twod
INTEGER(i_def), OPTIONAL, INTENT(IN)        :: ndata

! Local variables
TYPE(field_type)           :: new_field
INTEGER(i_def)             :: ndat
INTEGER(i_def), PARAMETER  :: fs_order = 0

! Pointers
TYPE(function_space_type),       POINTER  :: w3_space => NULL()
TYPE(function_space_type),       POINTER  :: twod_space => NULL()
PROCEDURE(read_interface),       POINTER  :: tmp_read_ptr => NULL()
PROCEDURE(write_interface),      POINTER  :: tmp_write_ptr => NULL()
TYPE(field_type),                POINTER  :: tgt_ptr => NULL()

! Set field ndata if argument is present, else leave as default value
IF (PRESENT(ndata)) THEN
  ndat = ndata
ELSE
  ndat = 1
END IF

! Set up function spaces for field initialisation
w3_space   => function_space_collection%get_fs( mesh, fs_order, &
              W3, ndat )
twod_space => function_space_collection%get_fs( twod_mesh, fs_order, &
               W3, ndat )

! Create field
WRITE(log_scratch_space,'(3A,I6)') &
     "Creating new field for ", TRIM(name)
CALL log_event(log_scratch_space,LOG_LEVEL_INFO)
IF (PRESENT(twod)) THEN
  CALL new_field%initialise( twod_space, name=TRIM(name), &
                             halo_depth = twod_mesh%get_halo_depth() )
ELSE
  CALL new_field%initialise( w3_space, name=TRIM(name), &
                             halo_depth = mesh%get_halo_depth() )
END IF

! Add the new field to the field collection
CALL ancil_fields%add_field(new_field)

! Get a field pointer from the collection
call ancil_fields%get_field(name, tgt_ptr)

! Set up field read behaviour for 2D and 3D fields
tmp_read_ptr => read_field_generic
tmp_write_ptr => write_field_generic

! Set field read behaviour for target field
CALL tgt_ptr%set_read_behaviour(tmp_read_ptr)
CALL tgt_ptr%set_write_behaviour(tmp_write_ptr)


! Nullify pointers
NULLIFY(w3_space)
NULLIFY(twod_space)
NULLIFY(tmp_read_ptr)
NULLIFY(tmp_write_ptr)
NULLIFY(tgt_ptr)

END SUBROUTINE lfricinp_setup_ancil_field

END MODULE lfricinp_ancils_mod
