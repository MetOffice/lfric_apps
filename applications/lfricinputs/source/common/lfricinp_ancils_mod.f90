! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_ancils_mod

use constants_mod,                  only : i_def, l_def
use log_mod,                        only : log_event,         &
                                           log_scratch_space, &
                                           LOG_LEVEL_INFO
use field_collection_mod,           only: field_collection_type
use field_mod,                      only : field_type
use field_parent_mod,               only : read_interface, &
                                           write_interface
use lfric_xios_read_mod,            only : read_field_generic
use lfric_xios_write_mod,           only : write_field_generic
use field_collection_mod,           only : field_collection_type
use mesh_mod,                       only : mesh_type
use function_space_mod,             only : function_space_type
use function_space_collection_mod,  only : function_space_collection
use fs_continuity_mod,              only : W3

implicit none

logical :: l_land_area_fraction = .false.
! Container for ancil fields
type(field_collection_type) :: ancil_fields

public   :: lfricinp_create_ancil_fields,         &
            lfricinp_setup_ancil_field,           &
            l_land_area_fraction

contains

! Organises fields to be read from ancils into ancil_fields collection

subroutine lfricinp_create_ancil_fields( ancil_fields, mesh, twod_mesh )

implicit none

type( field_collection_type ), intent( out )         :: ancil_fields
type( mesh_type ),             intent( in ), pointer :: mesh
type( mesh_type ),             intent( in ), pointer :: twod_mesh

! Set up ancil_fields collection
write(log_scratch_space,'(A,A)') "Create ancil fields: "// &
     "Setting up ancil field collection"
call log_event(log_scratch_space, LOG_LEVEL_INFO)
call ancil_fields%initialise(name='ancil_fields', table_len=100)

if (l_land_area_fraction) then
  ! Surface ancils
  call log_event("Create land area fraction ancil", LOG_LEVEL_INFO)
  call lfricinp_setup_ancil_field("land_area_fraction", ancil_fields, mesh, &
       twod_mesh, twod=.true.)
end if

end subroutine lfricinp_create_ancil_fields

!------------------------------------------------------------

! Creates fields to be read into from ancillary files

subroutine lfricinp_setup_ancil_field( name, ancil_fields, mesh, &
                                       twod_mesh, twod, ndata )

implicit none

character(*), intent(in)                    :: name
type(field_collection_type), intent(in out) :: ancil_fields
type( mesh_type ), intent(in), pointer      :: mesh
type( mesh_type ), intent(in), pointer      :: twod_mesh
logical(l_def), optional, intent(in)        :: twod
integer(i_def), optional, intent(in)        :: ndata

! Local variables
type(field_type)           :: new_field
integer(i_def)             :: ndat
integer(i_def), parameter  :: fs_order_h = 0
integer(i_def), parameter  :: fs_order_v = 0

! Pointers
type(function_space_type),       pointer  :: w3_space => null()
type(function_space_type),       pointer  :: twod_space => null()
procedure(read_interface),       pointer  :: tmp_read_ptr => null()
procedure(write_interface),      pointer  :: tmp_write_ptr => null()
type(field_type),                pointer  :: tgt_ptr => null()

! Set field ndata if argument is present, else leave as default value
if (present(ndata)) then
  ndat = ndata
else
  ndat = 1
end if

! Set up function spaces for field initialisation
w3_space   => function_space_collection%get_fs( mesh, fs_order_h, &
              fs_order_v, W3, ndat )
twod_space => function_space_collection%get_fs( twod_mesh, fs_order_h, &
               fs_order_v, W3, ndat )

! Create field
write(log_scratch_space,'(3A,I6)') &
     "Creating new field for ", trim(name)
call log_event(log_scratch_space,LOG_LEVEL_INFO)
if (present(twod)) then
  call new_field%initialise( twod_space, name=trim(name), &
                             halo_depth = twod_mesh%get_halo_depth() )
else
  call new_field%initialise( w3_space, name=trim(name), &
                             halo_depth = mesh%get_halo_depth() )
end if

! Add the new field to the field collection
call ancil_fields%add_field(new_field)

! Get a field pointer from the collection
call ancil_fields%get_field(name, tgt_ptr)

! Set up field read behaviour for 2D and 3D fields
tmp_read_ptr => read_field_generic
tmp_write_ptr => write_field_generic

! Set field read behaviour for target field
call tgt_ptr%set_read_behaviour(tmp_read_ptr)
call tgt_ptr%set_write_behaviour(tmp_write_ptr)


! Nullify pointers
nullify(w3_space)
nullify(twod_space)
nullify(tmp_read_ptr)
nullify(tmp_write_ptr)
nullify(tgt_ptr)

end subroutine lfricinp_setup_ancil_field

end module lfricinp_ancils_mod
