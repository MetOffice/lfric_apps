! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module um2lfric_regrid_weights_mod

use log_mod, only: log_event, LOG_LEVEL_ERROR, log_scratch_space, &
                   LOG_LEVEL_INFO
use lfricinp_um_parameters_mod, only: fnamelen
use lfricinp_regrid_weights_type_mod, only: lfricinp_regrid_weights_type
use lfricinp_um_grid_mod, only: um_grid
use um2lfric_namelist_mod, only: um2lfric_config

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int32, int64, real64

implicit none

private

! Regridding Weights
type(lfricinp_regrid_weights_type), public, target ::                          &
                                    grid_p_to_mesh_face_centre_bilinear,       &
                                    grid_p_to_mesh_face_centre_neareststod,    &
                                    grid_u_to_mesh_face_centre_bilinear,       &
                                    grid_v_to_mesh_face_centre_bilinear,       &
                                    grid_p_to_w3_wtheta_map, grid_u_to_w2h_map,&
                                    grid_v_to_w2h_map

public ::  um2lfric_regrid_weightsfile_ctl, get_weights

contains

!---------------------------------------------------------
!> @brief   Fetches the appropriate weights object for a given stashcode
!> @param[in]   stashcode   Integer value of stashcode to return weights for
!> @return      lfricinp_regrid_weights_type object
function get_weights(stashcode) result (weights)

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int64
! UM2LFRic modules
use lfricinp_stashmaster_mod,      only: stashmaster,                  &
                                         p_points, u_points, v_points, &
                                         ozone_points,                 &
                                         land_compressed,              &
                                         p_points_values_over_sea
use lfricinp_regrid_options_mod,   only: interp_method, &
                                         winds_on_w3, &
                                         specify_nearest_neighbour, &
                                         nn_fields

implicit none

! Arguments
integer(kind=int64), intent(in) :: stashcode
! Result
type(lfricinp_regrid_weights_type), pointer :: weights

! Local variables
integer(kind=int64) :: horiz_grid_code = 0
integer(kind=int64) :: i_stash
logical :: unspecified

! Check that the STASHmaster record exists.
if (associated(stashmaster(stashcode) % record)) then
  ! Get grid type code from STASHmaster entry
  horiz_grid_code = stashmaster(stashcode) % record % grid
else
  write(log_scratch_space, '(A,I0)')                             &
       "Unassociated STASHmaster record for stashcode ", stashcode
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)
end if

if (horiz_grid_code == u_points) then

  if (trim(interp_method) == 'copy' .and. .not. winds_on_w3) then
    weights => grid_u_to_w2h_map
  else if (trim(interp_method) == 'bilinear' .and. winds_on_w3) then
    weights => grid_u_to_mesh_face_centre_bilinear
  else
    write(log_scratch_space, '(A)')                                            &
                                'Unsupported interpolation method for U points'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

else if (horiz_grid_code == v_points) then

  if (trim(interp_method) == 'copy' .and. .not. winds_on_w3) then
    weights => grid_v_to_w2h_map
  else if (trim(interp_method) == 'bilinear' .and. winds_on_w3) then
    weights => grid_v_to_mesh_face_centre_bilinear
  else
    write(log_scratch_space, '(A)')                                            &
                                'Unsupported interpolation method for V points'
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if

else if (horiz_grid_code == p_points .or.   &
         horiz_grid_code == ozone_points .or. &
         horiz_grid_code == land_compressed .or. &
         horiz_grid_code == p_points_values_over_sea) then

  ! Check for specified interpolation method for this stashcode
  ! and the use appropriate weights
  unspecified = .true.
  if ( nn_fields > 0 ) then
    do i_stash = 1,nn_fields
      if (stashcode == specify_nearest_neighbour(i_stash))then
        weights => grid_p_to_mesh_face_centre_neareststod
        unspecified = .false.

        write(log_scratch_space, '((A,I4))')                          &
           "Will use nearest neigbour interpolation for stashcode: ", &
           stashcode
        call log_event(log_scratch_space, LOG_LEVEL_INFO)
        exit
      end if
    end do
  end if

  ! If no specific interpolation specified for this stashcode, then
  ! use the default
  if (unspecified)then

   if (trim(interp_method) == 'copy') then
     weights => grid_p_to_w3_wtheta_map
   else if (trim(interp_method) == 'bilinear') then
     weights => grid_p_to_mesh_face_centre_bilinear
   else
     write(log_scratch_space, '(A)')                                           &
                                 'Unsupported interpolation method for P points'
     call log_event(log_scratch_space, LOG_LEVEL_ERROR)
   ENDIF

 end if

else

  write(log_scratch_space, '(2(A,I0))')                                        &
       "Unsupported horizontal grid type code: ",                              &
       horiz_grid_code, " encountered during regrid of stashcode", stashcode
  call log_event(log_scratch_space, LOG_LEVEL_ERROR)

end if

if (.not. allocated(weights%remap_matrix)) then
  call log_event("Attempted to select unallocated weights matrix",             &
                 LOG_LEVEL_ERROR)
end if

end function get_weights

!> @brief Initialise all regridding weights objects
!> @details Takes all SCRIP weights files generated for um2lfric and initialises
!!          objects for each, validating and partitioning based on local mesh.
subroutine um2lfric_regrid_weightsfile_ctl()

implicit none

! P points to face centre bilinear interpolation
call grid_p_to_mesh_face_centre_bilinear % load(                               &
     um2lfric_config%weights_file_p_to_face_centre_bilinear)
call grid_p_to_mesh_face_centre_bilinear % validate_src(                       &
     int(um_grid % num_p_points_x, kind=int32) *                               &
     int(um_grid % num_p_points_y, kind=int32))
call grid_p_to_mesh_face_centre_bilinear % populate_src_address_2D(            &
     int(um_grid % num_p_points_x, kind=int32))
call partition_weights(grid_p_to_mesh_face_centre_bilinear)

! P points to face centre nearest neighbour interpolation
call grid_p_to_mesh_face_centre_neareststod % load(                            &
     um2lfric_config%weights_file_p_to_face_centre_neareststod)
call grid_p_to_mesh_face_centre_neareststod % validate_src(                    &
     int(um_grid % num_p_points_x, kind=int32) *                               &
     int(um_grid % num_p_points_y, kind=int32))
call grid_p_to_mesh_face_centre_neareststod % populate_src_address_2D(         &
     int(um_grid % num_p_points_x, kind=int32))
call partition_weights(grid_p_to_mesh_face_centre_neareststod)

! U points to face centre bilinear interpolation
call grid_u_to_mesh_face_centre_bilinear % load(                               &
     um2lfric_config%weights_file_u_to_face_centre_bilinear)
call grid_u_to_mesh_face_centre_bilinear % validate_src(                       &
     int(um_grid % num_u_points_x, kind=int32) *                               &
     int(um_grid % num_u_points_y, kind=int32))
call grid_u_to_mesh_face_centre_bilinear % populate_src_address_2D(            &
     int(um_grid % num_u_points_x, kind=int32))
call partition_weights(grid_u_to_mesh_face_centre_bilinear)

! V points to face center bilinear interpolation
call grid_v_to_mesh_face_centre_bilinear % load(                               &
     um2lfric_config%weights_file_v_to_face_centre_bilinear)
call grid_v_to_mesh_face_centre_bilinear % validate_src(                       &
     int(um_grid % num_v_points_x, kind=int32) *                               &
     int(um_grid % num_v_points_y, kind=int32))
call grid_v_to_mesh_face_centre_bilinear % populate_src_address_2D(            &
     int(um_grid % num_v_points_x, kind=int32))
call partition_weights(grid_v_to_mesh_face_centre_bilinear)

! Set up P points to W3/Wtheta copy map
call grid_p_to_w3_wtheta_map % load(                            &
     um2lfric_config%weights_file_p_to_face_centre_neareststod)
call grid_p_to_w3_wtheta_map % validate_src(                    &
     int(um_grid % num_p_points_x, kind=int32) *                &
     int(um_grid % num_p_points_y, kind=int32))
call grid_p_to_w3_wtheta_map % populate_src_address_2D(         &
     int(um_grid % num_p_points_x, kind=int32))
call partition_weights(grid_p_to_w3_wtheta_map)

! Create W2H copy maps
call create_w2h_copy_maps()

end subroutine um2lfric_regrid_weightsfile_ctl


subroutine create_w2h_copy_maps()
!
! This routine determines the mappings between LFRic mesh edges and their
! corresponding UM grid U and V points. The basic working assumption is that
! the LFRic mesh (2D foot-print of the 3D mesh) exactly corresponds to the
! UM Arakawa C-grid (ENDGAME). The mappings are of the lfricinputs standard
! regrid weights type. The mappings should purely be used for transfer
! of field data from the UM fields to LFRic field objects across the two
! identical grids. The mappings have no meaning outside of this context.
!
! The mappings are determined as follows:
!
! 1) Given the nearest neighbour UM grid P to LFRic face centre map (i.e.
! regrid weight type) for each UM grid P point index, find the corresponding
! LFRic cell index at the base of the LFRic mesh.
!
! 2) For each reference element W-edge index of the LFRic cell, connect it to
! the corresponding UM U point index which is just west of the UM P point. The
! edge index of the LFRic cell is found using the W2H dofmap. Special care
! needs to be taken to ensure only dofs belonging to the current partition are
! considered. The UM grid U point index is found using the fact that the UM grid
! is a rectangular structured grid in lat-lon coordinates and the basic layout
! of the UM ENDGAME grid.
!
! 3) If the UM P-point is on the extreme eastern boundary, also set up the
! mapping between the UM extreme eastern U point index, and the LFRic
! reference element E-edge index.
!
! 4) Apply a similar procedure to points 2 and 3 above to connect the UM V-point
! indices to the LFRic reference element N and S edge indices.


use function_space_collection_mod, only: function_space_collection
use function_space_mod,            only: function_space_type
use fs_continuity_mod,             only: W2H
use finite_element_config_mod,     only: element_order_h, element_order_v
use reference_element_mod,         only: W, S, E, N
use lfricinp_lfric_driver_mod,     only: mesh
use lfricinp_um_grid_mod,          only: um_grid

implicit none

! Local variables
type(function_space_type), pointer  :: fs_w2h => null()
integer(kind=int32),       pointer  :: map_w2h(:,:) => null()
integer(kind=int32)              :: num_links_u, num_links_v, num_links_p,     &
                                    num_wgts, no_layers, cell_index_um_2d(2),  &
                                    cell_index_lfric, wgt, lnk, base_dof_id_u, &
                                    base_dof_id_v, base_last_dof_owned
integer(kind=int32), allocatable :: src_address_2d(:,:), dst_address(:)
logical                          :: l_new_link

! Deallocate all arrays first
if (allocated(grid_u_to_w2h_map % remap_matrix)) then
  deallocate(grid_u_to_w2h_map % remap_matrix)
end if
if (allocated(grid_u_to_w2h_map % src_address)) then
  deallocate(grid_u_to_w2h_map % src_address)
end if
if (allocated(grid_u_to_w2h_map % dst_address)) then
  deallocate(grid_u_to_w2h_map % dst_address)
end if
if (allocated(grid_u_to_w2h_map % src_address_2d)) then
  deallocate(grid_u_to_w2h_map % src_address_2d)
end if
if (allocated(grid_u_to_w2h_map % dst_address_2d)) then
  deallocate(grid_u_to_w2h_map % dst_address_2d)
ENDIF
!
if (allocated(grid_v_to_w2h_map % remap_matrix)) then
  deallocate(grid_v_to_w2h_map % remap_matrix)
end if
if (allocated(grid_v_to_w2h_map % src_address)) then
  deallocate(grid_v_to_w2h_map % src_address)
end if
if (allocated(grid_v_to_w2h_map % dst_address)) then
  deallocate(grid_v_to_w2h_map % dst_address)
end if
if (allocated(grid_v_to_w2h_map % src_address_2d)) then
  deallocate(grid_v_to_w2h_map % src_address_2d)
end if
if (allocated(grid_v_to_w2h_map % dst_address_2d)) then
  deallocate(grid_v_to_w2h_map % dst_address_2d)
end if

! Set global number of source and destination points
grid_u_to_w2h_map % num_points_src = int(um_grid % num_u_points_x *            &
                                         um_grid % num_u_points_y, kind=int32)
grid_u_to_w2h_map % num_points_dst = grid_u_to_w2h_map % num_points_src +      &
                                     int(um_grid % num_u_points_y, kind=int32)
!
grid_v_to_w2h_map % num_points_src = int(um_grid % num_v_points_x *            &
                                         um_grid % num_v_points_y, kind=int32)
grid_v_to_w2h_map % num_points_dst = grid_v_to_w2h_map % num_points_src

! Get whole W2H dofmap
fs_w2h => function_space_collection % get_fs(mesh, element_order_h, &
                                             element_order_v, W2H)
map_w2h => fs_w2h % get_whole_dofmap()
no_layers = fs_w2h % get_nlayers()
base_last_dof_owned =  (fs_w2h % get_last_dof_owned() / no_layers ) + 1
write(log_scratch_space, '(A,I0)') 'Maximum base dof index = ',                &
                                   base_last_dof_owned
call log_event(log_scratch_space, LOG_LEVEL_INFO)

! Allocate tempory source and destination arrays to rough maximum possible size
num_links_p = grid_p_to_w3_wtheta_map % num_links
!
allocate(src_address_2d(2*num_links_p,2))
allocate(dst_address(2*num_links_p))

! Fill U to W2H remap weights src and dst arrays
src_address_2d(:,:) = 0
dst_address(:)      = 0
num_links_u = 0
do lnk = 1, num_links_p

  cell_index_um_2d(:) = grid_p_to_w3_wtheta_map % src_address_2d(lnk,:)
  cell_index_lfric = grid_p_to_w3_wtheta_map % dst_address(lnk)

  base_dof_id_u = (map_w2h(W,cell_index_lfric) / no_layers) + 1
  l_new_link = (.not. any(dst_address == base_dof_id_u))
  l_new_link = l_new_link .and. (base_dof_id_u <= base_last_dof_owned)
  if (l_new_link) then
    num_links_u = num_links_u + 1
    src_address_2d(num_links_u,1) = cell_index_um_2d(1)
    src_address_2d(num_links_u,2) = cell_index_um_2d(2)
    dst_address(num_links_u)      = base_dof_id_u
  end if

  if (cell_index_um_2d(1) == um_grid % num_p_points_x) then
    base_dof_id_u = (map_w2h(E,cell_index_lfric) / no_layers) + 1
    l_new_link = (.not. any(dst_address == base_dof_id_u))
    l_new_link = l_new_link .and. (base_dof_id_u <= base_last_dof_owned)
    if (l_new_link) then
      num_links_u = num_links_u + 1
      src_address_2d(num_links_u,1) = cell_index_um_2d(1)
      src_address_2d(num_links_u,2) = cell_index_um_2d(2)
      dst_address(num_links_u)      = base_dof_id_u
    end if
  end if

end do
allocate(grid_u_to_w2h_map % src_address_2d(num_links_u,2))
allocate(grid_u_to_w2h_map % dst_address(num_links_u))
grid_u_to_w2h_map % src_address_2d(1:num_links_u,:) =                          &
                                                src_address_2d(1:num_links_u,:)
grid_u_to_w2h_map % dst_address(1:num_links_u) = dst_address(1:num_links_u)

! Fill V to W2H remap weights src and dst arrays
src_address_2d(:,:) = 0
dst_address(:)      = 0
num_links_v = 0
do lnk = 1, num_links_p

  cell_index_um_2d(:) = grid_p_to_w3_wtheta_map % src_address_2d(lnk,:)
  cell_index_lfric = grid_p_to_w3_wtheta_map % dst_address(lnk)

  base_dof_id_v = (map_w2h(S,cell_index_lfric) / no_layers) + 1
  l_new_link = (.not. any(dst_address == base_dof_id_v))
  l_new_link = l_new_link .and. (base_dof_id_v <= base_last_dof_owned)
  if (l_new_link) then
    num_links_v = num_links_v + 1
    src_address_2d(num_links_v,1) = cell_index_um_2d(1)
    src_address_2d(num_links_v,2) = cell_index_um_2d(2)
    dst_address(num_links_v)      = base_dof_id_v
  end if

  if (cell_index_um_2d(2) == um_grid % num_p_points_y) then
    base_dof_id_v = (map_w2h(N,cell_index_lfric) / no_layers) + 1
    l_new_link = (.not. any(dst_address == base_dof_id_v))
    l_new_link = l_new_link .and. (base_dof_id_v <= base_last_dof_owned)
    if (l_new_link) then
      num_links_v = num_links_v + 1
      src_address_2d(num_links_v,1) = cell_index_um_2d(1)
      src_address_2d(num_links_v,2) = cell_index_um_2d(2) + 1
      dst_address(num_links_v)      = base_dof_id_v
    end if
  end if

end do
allocate(grid_v_to_w2h_map % src_address_2d(num_links_v,2))
allocate(grid_v_to_w2h_map % dst_address(num_links_v))
grid_v_to_w2h_map % src_address_2d(1:num_links_v,:) =                          &
                                                src_address_2d(1:num_links_v,:)
grid_v_to_w2h_map % dst_address(1:num_links_v) = dst_address(1:num_links_v)

! Deallocate temporary arrays and nullify pointers
deallocate(src_address_2d, dst_address)
nullify(fs_w2h, map_w2h)

! Set number of links and set num_wgts
grid_u_to_w2h_map % num_links = num_links_u
grid_v_to_w2h_map % num_links = num_links_v
!
num_wgts = 1
grid_u_to_w2h_map % num_wgts = num_wgts
grid_v_to_w2h_map % num_wgts = num_wgts

! Fill U and V regrid weights remap matrices
allocate(grid_u_to_w2h_map % remap_matrix(num_wgts, num_links_u))
allocate(grid_v_to_w2h_map % remap_matrix(num_wgts, num_links_v))
do wgt = 1, num_wgts
  do lnk = 1, num_links_u
    grid_u_to_w2h_map % remap_matrix(wgt, lnk) = 1.0_real64
  end do
  do lnk = 1, num_links_v
    grid_v_to_w2h_map % remap_matrix(wgt, lnk) = 1.0_real64
  end do
end do

end subroutine create_w2h_copy_maps


subroutine partition_weights(weights)

  ! Description: Set local weights indices so that the weights type on each
  ! partition only contains indices relevant to that local_mesh.

  ! Intrinsic modules
  use, intrinsic :: iso_fortran_env, only : int32, real64

  ! LFRic modules
  use local_mesh_mod, only: local_mesh_type

  ! LFRic Inputs modules
  use lfricinp_lfric_driver_mod, only: mesh

  implicit none

  ! Arguments

  type(lfricinp_regrid_weights_type), intent(INOUT) :: weights

  ! Local variables
  type(local_mesh_type), pointer :: local_mesh => null()

  ! Local local_mesh versions of weights arrays
  integer(kind=int32), allocatable :: dst_address_local_tmp(:)
  integer(kind=int32), allocatable :: src_address_local_tmp(:)
  integer(kind=int32), allocatable :: src_address_local_2D_tmp(:,:)
  real(kind=real64),   allocatable :: remap_matrix_local_tmp(:,:)

  integer(kind=int32) :: local_links, i_link, w

  ! Set up pointer to LFRic mesh
  local_mesh => mesh % get_local_mesh()

  ! Initially allocate to be the same size as global arrays
  allocate( dst_address_local_tmp( weights%num_links ) )
  allocate( src_address_local_tmp( weights%num_links ) )
  allocate( src_address_local_2D_tmp( weights%num_links, 2) )
  allocate( remap_matrix_local_tmp( weights%num_wgts, weights%num_links) )
  ! Convert weights from global indices to the local local_meshs id/indices
  ! Use local_links to count how many links on local partition
  local_links = 0
  do i_link = 1, weights%num_links
    ! The call to get_lid_from_gid will return -1 if the global id is not
    ! known to this partition. Also the partition will also contain halo cells.
    ! However, we only want to select cells wholy owned by this partition. So
    ! reject any addresses that are -1 or greater than ncells_2D
    if ( local_mesh % get_lid_from_gid(weights%dst_address(i_link)) /= -1 .and.&
         local_mesh % get_lid_from_gid(weights%dst_address(i_link)) <=         &
           mesh % get_ncells_2d() ) then
      local_links = local_links + 1
      dst_address_local_tmp(local_links) =                                     &
                      local_mesh % get_lid_from_gid(weights%dst_address(i_link))
      ! Weights and source address also need filtering
      src_address_local_tmp(local_links) = weights%src_address(i_link)
      src_address_local_2D_tmp(local_links, 1) =                               &
                                              weights%src_address_2D(i_link, 1)
      src_address_local_2D_tmp(local_links, 2) =                               &
                                              weights%src_address_2D(i_link, 2)
      do w = 1, weights%num_wgts
        remap_matrix_local_tmp(w, local_links) = weights%remap_matrix(w, i_link)
      end do
    end if
  end do

  ! Deallocate and reallocate the arrays in the main weights type using
  ! the local sizes as calculated by the local_links counter in the above loop
  if (allocated(weights%src_address_2D)) deallocate(weights%src_address_2D)
  allocate(weights%src_address_2D(local_links, 2))
  if (allocated(weights%src_address)) deallocate(weights%src_address)
  allocate( weights%src_address( local_links ))
  if (allocated(weights%dst_address)) deallocate(weights%dst_address)
  allocate( weights%dst_address( local_links ))
  if (allocated(weights%remap_matrix)) deallocate(weights%remap_matrix)
  allocate(weights%remap_matrix( weights%num_wgts, local_links))

  weights % num_links = local_links
  ! Copy values from the tmp versions of arrays to the main weights type
  do i_link = 1, weights%num_links
    weights%dst_address(i_link) = dst_address_local_tmp(i_link)
    weights%src_address(i_link) = src_address_local_tmp(i_link)
    weights%src_address_2D(i_link,1) = src_address_local_2D_tmp(i_link,1)
    weights%src_address_2D(i_link,2) = src_address_local_2D_tmp(i_link,2)
    weights%remap_matrix(:,i_link) = remap_matrix_local_tmp(:,i_link)
  end do

  deallocate(dst_address_local_tmp)
  deallocate(src_address_local_tmp)
  deallocate(src_address_local_2D_tmp)
  deallocate(remap_matrix_local_tmp)
  nullify(local_mesh)

end subroutine partition_weights

end module um2lfric_regrid_weights_mod
