! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_create_lfric_fields_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int64

! lfric modules
use constants_mod,                 only : i_def, str_def
use field_mod,                     only : field_type
use field_parent_mod,              only :  read_interface, write_interface,    &
                                           field_parent_proxy_type
use lfric_xios_time_axis_mod,      only : time_axis_type
use lfric_xios_read_mod,           only : read_field_generic
use lfric_xios_write_mod,          only : write_field_generic
use finite_element_config_mod,     only : element_order_h, element_order_v
use function_space_mod,            only : function_space_type
use function_space_collection_mod, only : function_space_collection
use field_collection_mod,          only : field_collection_type
use fs_continuity_mod,             only : W3, Wtheta, W2H
use log_mod,                       only : LOG_LEVEL_INFO, LOG_LEVEL_ERROR,     &
                                          log_event, log_scratch_space
use mesh_mod,                      only : mesh_type

! lfricinp modules
use lfricinp_stashmaster_mod,        only: get_stashmaster_item, pseudt
use lfricinp_stash_to_lfric_map_mod, only: w2h_field, w3_field, w3_field_2d,   &
                                           w3_soil_field, wtheta_field,        &
                                           get_field_name, get_lfric_field_kind
use lfricinp_grid_type_mod,          only: lfricinp_grid_type
use lfricinp_um_level_codes_mod,     only: lfricinp_get_num_levels,            &
                                           lfricinp_get_num_pseudo_levels
use lfricinp_regrid_options_mod,     only: winds_on_w3

! shumlib modules
use f_shum_file_mod, only: shum_file_type

implicit none
private
public :: lfricinp_create_lfric_fields

contains

subroutine lfricinp_create_lfric_fields( mesh, twod_mesh,              &
                                         field_collection, stash_list, &
                                         um_grid, um_file)
! Description:
!  When given list of stashcodes create lfric field. Use stashcode to determine
!  which vertical level set the field is on and create appropiate function
!  space. Set read and write procedure for each field.


implicit none

type(mesh_type), intent(in), pointer :: mesh
type(mesh_type), intent(in), pointer :: twod_mesh

type(field_collection_type), intent(in out) :: field_collection

integer(kind=int64), intent(in)  :: stash_list(:)

type(lfricinp_grid_type), intent(in):: um_grid

type(shum_file_type), intent(in) :: um_file

procedure(read_interface),  pointer  :: tmp_read_ptr => null()
procedure(write_interface), pointer  :: tmp_write_ptr => null()
type(function_space_type),  pointer :: vector_space => null()

type(mesh_type), pointer :: type_mesh => null()

type( field_type )            :: field
integer(kind=int64)           :: stashcode, lfric_field_kind, i
integer(kind=i_def)           :: fs_id
! A temporary variable to hold the 64bit output
integer(kind=int64)           :: ndata_64
character(len=:), allocatable :: field_name
logical                       :: ndata_first

call log_event('Creating lfric finite difference fields ...', LOG_LEVEL_INFO)

if (element_order_h > 0 .or. element_order_v  > 0)        &
  call log_event('Finite difference fields requires ' //  &
                                      'lowest order elements', LOG_LEVEL_ERROR)

! Create the field collection
call field_collection%initialise(name="lfric_fields", table_len=100)

do i=1, size(stash_list)

  stashcode = stash_list(i)
  lfric_field_kind = get_lfric_field_kind(stashcode)
  field_name = trim(get_field_name(stashcode))

  ! TODO - It would be good to move this case statement into its own
  ! routine it could be used elsewhere to add fields to the collection
  ! without the do loop over stash list
  if ( .not. field_collection%field_exists(field_name) ) then

    tmp_read_ptr => read_field_generic
    tmp_write_ptr => write_field_generic

    select case (lfric_field_kind)

      case(w2h_field) ! Stashcodes that would map to W2h, i.e. winds
        type_mesh => mesh
        fs_id = W2H
        ndata_64 = 1_int64
        ndata_first = .false.

      case(w3_field) ! Stashcodes that map to W3/rho
        type_mesh => mesh
        fs_id = W3
        ndata_64 = 1_int64
        ndata_first = .false.

      case(wtheta_field) ! Stashcodes that maps to Wtheta
        type_mesh => mesh
        fs_id = Wtheta
        ndata_64 = 1_int64
        ndata_first = .false.

      case(w3_field_2d) ! Stash that needs 2D mesh
        type_mesh => twod_mesh
        fs_id = W3
        if ( get_stashmaster_item(stashcode, pseudt) == 0 ) then
          ! Field has no pseudo levels
          ndata_64 = 1_int64
        else
          ! Get number of pseudo levels/ndata
          ndata_64 = lfricinp_get_num_pseudo_levels(um_grid, stashcode)
        end if
        ndata_first = .false.

      case(w3_soil_field) ! Soil fields
        type_mesh => twod_mesh
        fs_id = W3
        ndata_64 = lfricinp_get_num_levels(um_file, stashcode)
        ndata_first = .true.

      case DEFAULT
        write(log_scratch_space, '(A,I0,A)')                                   &
           "LFRic field kind code ", lfric_field_kind, " not recognised"
        call log_event(log_scratch_space, LOG_LEVEL_ERROR)

    end select

    vector_space => function_space_collection%get_fs(                          &
                                              type_mesh,                       &
                                              element_order_h,                 &
                                              element_order_v,                 &
                                              fs_id,                           &
                                              ndata=int(ndata_64, kind=i_def), &
                                              ndata_first=ndata_first          &
                                                    )
    call field % initialise(vector_space=vector_space, name=field_name, &
                            halo_depth = type_mesh%get_halo_depth() )
    call field % set_read_behaviour(tmp_read_ptr)
    call field % set_write_behaviour(tmp_write_ptr)
    call log_event("Add "//field_name//" to field collection", LOG_LEVEL_INFO)
    call field_collection % add_field(field)

    nullify(tmp_read_ptr, tmp_write_ptr)
    nullify(vector_space)

  end if

end do

call log_event('... Done creating finite difference fields', LOG_LEVEL_INFO)

end subroutine lfricinp_create_lfric_fields

end module lfricinp_create_lfric_fields_mod
