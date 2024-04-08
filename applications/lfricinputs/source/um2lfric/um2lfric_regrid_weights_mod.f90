! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE um2lfric_regrid_weights_mod

USE log_mod, ONLY: log_event, LOG_LEVEL_ERROR, log_scratch_space, &
                   LOG_LEVEL_INFO
USE lfricinp_um_parameters_mod, ONLY: fnamelen
USE lfricinp_regrid_weights_type_mod, ONLY: lfricinp_regrid_weights_type
USE lfricinp_um_grid_mod, ONLY: um_grid
USE um2lfric_namelist_mod, ONLY: um2lfric_config

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int32, int64, real64

IMPLICIT NONE

PRIVATE

! Regridding Weights
TYPE(lfricinp_regrid_weights_type), PUBLIC, TARGET ::                          &
                                    grid_p_to_mesh_face_centre_bilinear,       &
                                    grid_p_to_mesh_face_centre_neareststod,    &
                                    grid_u_to_mesh_face_centre_bilinear,       &
                                    grid_v_to_mesh_face_centre_bilinear,       &
                                    grid_p_to_w3_wtheta_map, grid_u_to_w2h_map,&
                                    grid_v_to_w2h_map

PUBLIC ::  um2lfric_regrid_weightsfile_ctl, get_weights

CONTAINS

!---------------------------------------------------------

FUNCTION get_weights(stashcode) RESULT (weights)

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY : int64
! UM2LFRic modules
USE lfricinp_stashmaster_mod, ONLY: stashmaster, p_points, &
   u_points, v_points, ozone_points, land_compressed,           &
   p_points_values_over_sea
USE lfricinp_regrid_options_mod, ONLY: interp_method, winds_on_w3, &
     specify_nearest_neighbour, nn_fields

IMPLICIT NONE

! Arguments
INTEGER(KIND=int64), INTENT(IN) :: stashcode
! Result
TYPE(lfricinp_regrid_weights_type), POINTER :: weights

! Local variables
INTEGER(KIND=int64) :: horiz_grid_code = 0
INTEGER(KIND=int64) :: i_stash
LOGICAL :: unspecified

! Check that the STASHmaster record exists.
IF (ASSOCIATED(stashmaster(stashcode) % record)) THEN
  ! Get grid type code from STASHmaster entry
  horiz_grid_code = stashmaster(stashcode) % record % grid
ELSE
  WRITE(log_scratch_space, '(A,I0)')                             &
       "Unassociated STASHmaster record for stashcode ", stashcode
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
END IF

IF (horiz_grid_code == u_points) THEN

  IF (TRIM(interp_method) == 'copy' .AND. .NOT. winds_on_w3) THEN
    weights => grid_u_to_w2h_map
  ELSE IF (TRIM(interp_method) == 'bilinear' .AND. winds_on_w3) THEN
    weights => grid_u_to_mesh_face_centre_bilinear
  ELSE
    WRITE(log_scratch_space, '(A)')                                            &
                                'Unsupported interpolation method for U points'
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END IF

ELSE IF (horiz_grid_code == v_points) THEN

  IF (TRIM(interp_method) == 'copy' .AND. .NOT. winds_on_w3) THEN
    weights => grid_v_to_w2h_map
  ELSE IF (TRIM(interp_method) == 'bilinear' .AND. winds_on_w3) THEN
    weights => grid_v_to_mesh_face_centre_bilinear
  ELSE
    WRITE(log_scratch_space, '(A)')                                            &
                                'Unsupported interpolation method for V points'
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END IF

ELSE IF (horiz_grid_code == p_points .OR.   &
         horiz_grid_code == ozone_points .OR. &
         horiz_grid_code == land_compressed .OR. &
         horiz_grid_code == p_points_values_over_sea) THEN

  ! Check for specified interpolation method for this stashcode
  ! and the use appropriate weights
  unspecified = .TRUE.
  IF ( nn_fields > 0 ) THEN
    DO i_stash = 1,nn_fields
      IF (stashcode == specify_nearest_neighbour(i_stash))THEN
        weights => grid_p_to_mesh_face_centre_neareststod
        unspecified = .FALSE.

        WRITE(log_scratch_space, '((A,I4))')                          &
           "Will use nearest neigbour interpolation for stashcode: ", &
           stashcode
        CALL log_event(log_scratch_space, LOG_LEVEL_INFO)
        EXIT
      END IF
    END DO
  END IF

  ! If no specific interpolation specified for this stashcode, then
  ! use the default
  IF (unspecified)THEN

   IF (TRIM(interp_method) == 'copy') THEN
     weights => grid_p_to_w3_wtheta_map
   ELSE IF (TRIM(interp_method) == 'bilinear') THEN
     weights => grid_p_to_mesh_face_centre_bilinear
   ELSE
     WRITE(log_scratch_space, '(A)')                                           &
                                 'Unsupported interpolation method for P points'
     CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
   ENDIF

 END IF

ELSE

  WRITE(log_scratch_space, '(2(A,I0))')                                        &
       "Unsupported horizontal grid type code: ",                              &
       horiz_grid_code, " encountered during regrid of stashcode", stashcode
  CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)

END IF

IF (.NOT. ALLOCATED(weights%remap_matrix)) THEN
  CALL log_event("Attempted to select unallocated weights matrix",             &
                 LOG_LEVEL_ERROR)
END IF

END FUNCTION get_weights


SUBROUTINE um2lfric_regrid_weightsfile_ctl()

IMPLICIT NONE

! P points to face centre bilinear interpolation
CALL grid_p_to_mesh_face_centre_bilinear % load(                               &
     um2lfric_config%weights_file_p_to_face_centre_bilinear)
CALL grid_p_to_mesh_face_centre_bilinear % validate_src(                       &
     INT(um_grid % num_p_points_x, KIND=int32) *                               &
     INT(um_grid % num_p_points_y, KIND=int32))
CALL grid_p_to_mesh_face_centre_bilinear % populate_src_address_2D(            &
     INT(um_grid % num_p_points_x, KIND=int32))
CALL partition_weights(grid_p_to_mesh_face_centre_bilinear)

! P points to face centre nearest neighbour interpolation
CALL grid_p_to_mesh_face_centre_neareststod % load(                            &
     um2lfric_config%weights_file_p_to_face_centre_neareststod)
CALL grid_p_to_mesh_face_centre_neareststod % validate_src(                    &
     INT(um_grid % num_p_points_x, KIND=int32) *                               &
     INT(um_grid % num_p_points_y, KIND=int32))
CALL grid_p_to_mesh_face_centre_neareststod % populate_src_address_2D(         &
     INT(um_grid % num_p_points_x, KIND=int32))
CALL partition_weights(grid_p_to_mesh_face_centre_neareststod)

! U points to face centre bilinear interpolation
CALL grid_u_to_mesh_face_centre_bilinear % load(                               &
     um2lfric_config%weights_file_u_to_face_centre_bilinear)
CALL grid_u_to_mesh_face_centre_bilinear % validate_src(                       &
     INT(um_grid % num_u_points_x, KIND=int32) *                               &
     INT(um_grid % num_u_points_y, KIND=int32))
CALL grid_u_to_mesh_face_centre_bilinear % populate_src_address_2D(            &
     INT(um_grid % num_u_points_x, KIND=int32))
CALL partition_weights(grid_u_to_mesh_face_centre_bilinear)

! V points to face center bilinear interpolation
CALL grid_v_to_mesh_face_centre_bilinear % load(                               &
     um2lfric_config%weights_file_v_to_face_centre_bilinear)
CALL grid_v_to_mesh_face_centre_bilinear % validate_src(                       &
     INT(um_grid % num_v_points_x, KIND=int32) *                               &
     INT(um_grid % num_v_points_y, KIND=int32))
CALL grid_v_to_mesh_face_centre_bilinear % populate_src_address_2D(            &
     INT(um_grid % num_v_points_x, KIND=int32))
CALL partition_weights(grid_v_to_mesh_face_centre_bilinear)

! Set up P points to W3/Wtheta copy map
CALL grid_p_to_w3_wtheta_map % load(                            &
     um2lfric_config%weights_file_p_to_face_centre_neareststod)
CALL grid_p_to_w3_wtheta_map % validate_src(                    &
     INT(um_grid % num_p_points_x, KIND=int32) *                &
     INT(um_grid % num_p_points_y, KIND=int32))
CALL grid_p_to_w3_wtheta_map % populate_src_address_2D(         &
     INT(um_grid % num_p_points_x, KIND=int32))
CALL partition_weights(grid_p_to_w3_wtheta_map)

! Create W2H copy maps
CALL create_w2h_copy_maps()

END SUBROUTINE um2lfric_regrid_weightsfile_ctl


SUBROUTINE create_w2h_copy_maps()
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


USE function_space_collection_mod, ONLY: function_space_collection
USE function_space_mod,            ONLY: function_space_type
USE fs_continuity_mod,             ONLY: W2H
USE finite_element_config_mod,     ONLY: element_order
USE reference_element_mod,         ONLY: W, S, E, N
USE lfricinp_lfric_driver_mod,     ONLY: mesh
USE lfricinp_um_grid_mod,          ONLY: um_grid

IMPLICIT NONE

! Local variables
TYPE(function_space_type), POINTER  :: fs_w2h => NULL()
INTEGER(KIND=int32),       POINTER  :: map_w2h(:,:) => NULL()
INTEGER(KIND=int32)              :: num_links_u, num_links_v, num_links_p,     &
                                    num_wgts, no_layers, cell_index_um_2d(2),  &
                                    cell_index_lfric, wgt, lnk, base_dof_id_u, &
                                    base_dof_id_v, base_last_dof_owned
INTEGER(KIND=int32), ALLOCATABLE :: src_address_2d(:,:), dst_address(:)
LOGICAL                          :: l_new_link

! Deallocate all arrays first
IF (ALLOCATED(grid_u_to_w2h_map % remap_matrix)) THEN
  DEALLOCATE(grid_u_to_w2h_map % remap_matrix)
END IF
IF (ALLOCATED(grid_u_to_w2h_map % src_address)) THEN
  DEALLOCATE(grid_u_to_w2h_map % src_address)
END IF
IF (ALLOCATED(grid_u_to_w2h_map % dst_address)) THEN
  DEALLOCATE(grid_u_to_w2h_map % dst_address)
END IF
IF (ALLOCATED(grid_u_to_w2h_map % src_address_2d)) THEN
  DEALLOCATE(grid_u_to_w2h_map % src_address_2d)
END IF
IF (ALLOCATED(grid_u_to_w2h_map % dst_address_2d)) THEN
  DEALLOCATE(grid_u_to_w2h_map % dst_address_2d)
ENDIF
!
IF (ALLOCATED(grid_v_to_w2h_map % remap_matrix)) THEN
  DEALLOCATE(grid_v_to_w2h_map % remap_matrix)
END IF
IF (ALLOCATED(grid_v_to_w2h_map % src_address)) THEN
  DEALLOCATE(grid_v_to_w2h_map % src_address)
END IF
IF (ALLOCATED(grid_v_to_w2h_map % dst_address)) THEN
  DEALLOCATE(grid_v_to_w2h_map % dst_address)
END IF
IF (ALLOCATED(grid_v_to_w2h_map % src_address_2d)) THEN
  DEALLOCATE(grid_v_to_w2h_map % src_address_2d)
END IF
IF (ALLOCATED(grid_v_to_w2h_map % dst_address_2d)) THEN
  DEALLOCATE(grid_v_to_w2h_map % dst_address_2d)
END IF

! Set global number of source and destination points
grid_u_to_w2h_map % num_points_src = INT(um_grid % num_u_points_x *            &
                                         um_grid % num_u_points_y, KIND=int32)
grid_u_to_w2h_map % num_points_dst = grid_u_to_w2h_map % num_points_src +      &
                                     INT(um_grid % num_u_points_y, KIND=int32)
!
grid_v_to_w2h_map % num_points_src = INT(um_grid % num_v_points_x *            &
                                         um_grid % num_v_points_y, KIND=int32)
grid_v_to_w2h_map % num_points_dst = grid_v_to_w2h_map % num_points_src

! Get whole W2H dofmap
fs_w2h => function_space_collection % get_fs(mesh, element_order, W2H)
map_w2h => fs_w2h % get_whole_dofmap()
no_layers = fs_w2h % get_nlayers()
base_last_dof_owned =  (fs_w2h % get_last_dof_owned() / no_layers ) + 1
WRITE(log_scratch_space, '(A,I0)') 'Maximum base dof index = ',                &
                                   base_last_dof_owned
CALL log_event(log_scratch_space, LOG_LEVEL_INFO)

! Allocate tempory source and destination arrays to rough maximum possible size
num_links_p = grid_p_to_w3_wtheta_map % num_links
!
ALLOCATE(src_address_2d(2*num_links_p,2))
ALLOCATE(dst_address(2*num_links_p))

! Fill U to W2H remap weights src and dst arrays
src_address_2d(:,:) = 0
dst_address(:)      = 0
num_links_u = 0
DO lnk = 1, num_links_p

  cell_index_um_2d(:) = grid_p_to_w3_wtheta_map % src_address_2d(lnk,:)
  cell_index_lfric = grid_p_to_w3_wtheta_map % dst_address(lnk)

  base_dof_id_u = (map_w2h(W,cell_index_lfric) / no_layers) + 1
  l_new_link = (.NOT. ANY(dst_address == base_dof_id_u))
  l_new_link = l_new_link .AND. (base_dof_id_u <= base_last_dof_owned)
  IF (l_new_link) THEN
    num_links_u = num_links_u + 1
    src_address_2d(num_links_u,1) = cell_index_um_2d(1)
    src_address_2d(num_links_u,2) = cell_index_um_2d(2)
    dst_address(num_links_u)      = base_dof_id_u
  END IF

  IF (cell_index_um_2d(1) == um_grid % num_p_points_x) THEN
    base_dof_id_u = (map_w2h(E,cell_index_lfric) / no_layers) + 1
    l_new_link = (.NOT. ANY(dst_address == base_dof_id_u))
    l_new_link = l_new_link .AND. (base_dof_id_u <= base_last_dof_owned)
    IF (l_new_link) THEN
      num_links_u = num_links_u + 1
      src_address_2d(num_links_u,1) = cell_index_um_2d(1)
      src_address_2d(num_links_u,2) = cell_index_um_2d(2)
      dst_address(num_links_u)      = base_dof_id_u
    END IF
  END IF

END DO
ALLOCATE(grid_u_to_w2h_map % src_address_2d(num_links_u,2))
ALLOCATE(grid_u_to_w2h_map % dst_address(num_links_u))
grid_u_to_w2h_map % src_address_2d(1:num_links_u,:) =                          &
                                                src_address_2d(1:num_links_u,:)
grid_u_to_w2h_map % dst_address(1:num_links_u) = dst_address(1:num_links_u)

! Fill V to W2H remap weights src and dst arrays
src_address_2d(:,:) = 0
dst_address(:)      = 0
num_links_v = 0
DO lnk = 1, num_links_p

  cell_index_um_2d(:) = grid_p_to_w3_wtheta_map % src_address_2d(lnk,:)
  cell_index_lfric = grid_p_to_w3_wtheta_map % dst_address(lnk)

  base_dof_id_v = (map_w2h(S,cell_index_lfric) / no_layers) + 1
  l_new_link = (.NOT. ANY(dst_address == base_dof_id_v))
  l_new_link = l_new_link .AND. (base_dof_id_v <= base_last_dof_owned)
  IF (l_new_link) THEN
    num_links_v = num_links_v + 1
    src_address_2d(num_links_v,1) = cell_index_um_2d(1)
    src_address_2d(num_links_v,2) = cell_index_um_2d(2)
    dst_address(num_links_v)      = base_dof_id_v
  END IF

  IF (cell_index_um_2d(2) == um_grid % num_p_points_y) THEN
    base_dof_id_v = (map_w2h(N,cell_index_lfric) / no_layers) + 1
    l_new_link = (.NOT. ANY(dst_address == base_dof_id_v))
    l_new_link = l_new_link .AND. (base_dof_id_v <= base_last_dof_owned)
    IF (l_new_link) THEN
      num_links_v = num_links_v + 1
      src_address_2d(num_links_v,1) = cell_index_um_2d(1)
      src_address_2d(num_links_v,2) = cell_index_um_2d(2) + 1
      dst_address(num_links_v)      = base_dof_id_v
    END IF
  END IF

END DO
ALLOCATE(grid_v_to_w2h_map % src_address_2d(num_links_v,2))
ALLOCATE(grid_v_to_w2h_map % dst_address(num_links_v))
grid_v_to_w2h_map % src_address_2d(1:num_links_v,:) =                          &
                                                src_address_2d(1:num_links_v,:)
grid_v_to_w2h_map % dst_address(1:num_links_v) = dst_address(1:num_links_v)

! Deallocate temporary arrays and nullify pointers
DEALLOCATE(src_address_2d, dst_address)
NULLIFY(fs_w2h, map_w2h)

! Set number of links and set num_wgts
grid_u_to_w2h_map % num_links = num_links_u
grid_v_to_w2h_map % num_links = num_links_v
!
num_wgts = 1
grid_u_to_w2h_map % num_wgts = num_wgts
grid_v_to_w2h_map % num_wgts = num_wgts

! Fill U and V regrid weights remap matrices
ALLOCATE(grid_u_to_w2h_map % remap_matrix(num_wgts, num_links_u))
ALLOCATE(grid_v_to_w2h_map % remap_matrix(num_wgts, num_links_v))
DO wgt = 1, num_wgts
  DO lnk = 1, num_links_u
    grid_u_to_w2h_map % remap_matrix(wgt, lnk) = 1.0
  END DO
  DO lnk = 1, num_links_v
    grid_v_to_w2h_map % remap_matrix(wgt, lnk) = 1.0
  END DO
END DO

END SUBROUTINE create_w2h_copy_maps


SUBROUTINE partition_weights(weights)

  ! Description: Set local weights indices so that the weights type on each
  ! partition only contains indices relevant to that local_mesh.

  ! Intrinsic modules
  USE, INTRINSIC :: iso_fortran_env, ONLY : int32, real64

  ! LFRic modules
  USE local_mesh_mod, ONLY: local_mesh_type

  ! LFRic Inputs modules
  USE lfricinp_lfric_driver_mod, ONLY: mesh

  IMPLICIT NONE

  ! Arguments

  TYPE(lfricinp_regrid_weights_type), INTENT(INOUT) :: weights

  ! Local variables
  TYPE(local_mesh_type), POINTER :: local_mesh => NULL()

  ! Local local_mesh versions of weights arrays
  INTEGER(KIND=int32), ALLOCATABLE :: dst_address_local_tmp(:)
  INTEGER(KIND=int32), ALLOCATABLE :: src_address_local_tmp(:)
  INTEGER(KIND=int32), ALLOCATABLE :: src_address_local_2D_tmp(:,:)
  REAL(KIND=real64),   ALLOCATABLE :: remap_matrix_local_tmp(:,:)

  INTEGER(KIND=int32) :: local_links, i_link, w

  ! Set up pointer to LFRic mesh
  local_mesh => mesh % get_local_mesh()

  ! Initially allocate to be the same size as global arrays
  ALLOCATE( dst_address_local_tmp( weights%num_links ) )
  ALLOCATE( src_address_local_tmp( weights%num_links ) )
  ALLOCATE( src_address_local_2D_tmp( weights%num_links, 2) )
  ALLOCATE( remap_matrix_local_tmp( weights%num_wgts, weights%num_links) )
  ! Convert weights from global indices to the local local_meshs id/indices
  ! Use local_links to count how many links on local partition
  local_links = 0
  DO i_link = 1, weights%num_links
    ! The call to get_lid_from_gid will return -1 if the global id is not
    ! known to this partition. Also the partition will also contain halo cells.
    ! However, we only want to select cells wholy owned by this partition. So
    ! reject any addresses that are -1 or greater than ncells_2D
    IF ( local_mesh % get_lid_from_gid(weights%dst_address(i_link)) /= -1 .AND.&
         local_mesh % get_lid_from_gid(weights%dst_address(i_link)) <=         &
           mesh % get_ncells_2d() ) THEN
      local_links = local_links + 1
      dst_address_local_tmp(local_links) =                                     &
                      local_mesh % get_lid_from_gid(weights%dst_address(i_link))
      ! Weights and source address also need filtering
      src_address_local_tmp(local_links) = weights%src_address(i_link)
      src_address_local_2D_tmp(local_links, 1) =                               &
                                              weights%src_address_2D(i_link, 1)
      src_address_local_2D_tmp(local_links, 2) =                               &
                                              weights%src_address_2D(i_link, 2)
      DO w = 1, weights%num_wgts
        remap_matrix_local_tmp(w, local_links) = weights%remap_matrix(w, i_link)
      END DO
    END IF
  END DO

  ! Deallocate and reallocate the arrays in the main weights type using
  ! the local sizes as calculated by the local_links counter in the above loop
  IF (ALLOCATED(weights%src_address_2D)) DEALLOCATE(weights%src_address_2D)
  ALLOCATE(weights%src_address_2D(local_links, 2))
  IF (ALLOCATED(weights%src_address)) DEALLOCATE(weights%src_address)
  ALLOCATE( weights%src_address( local_links ))
  IF (ALLOCATED(weights%dst_address)) DEALLOCATE(weights%dst_address)
  ALLOCATE( weights%dst_address( local_links ))
  IF (ALLOCATED(weights%remap_matrix)) DEALLOCATE(weights%remap_matrix)
  ALLOCATE(weights%remap_matrix( weights%num_wgts, local_links))

  weights % num_links = local_links
  ! Copy values from the tmp versions of arrays to the main weights type
  DO i_link = 1, weights%num_links
    weights%dst_address(i_link) = dst_address_local_tmp(i_link)
    weights%src_address(i_link) = src_address_local_tmp(i_link)
    weights%src_address_2D(i_link,1) = src_address_local_2D_tmp(i_link,1)
    weights%src_address_2D(i_link,2) = src_address_local_2D_tmp(i_link,2)
    weights%remap_matrix(:,i_link) = remap_matrix_local_tmp(:,i_link)
  END DO

  DEALLOCATE(dst_address_local_tmp)
  DEALLOCATE(src_address_local_tmp)
  DEALLOCATE(src_address_local_2D_tmp)
  DEALLOCATE(remap_matrix_local_tmp)
  NULLIFY(local_mesh)

END SUBROUTINE partition_weights

END MODULE um2lfric_regrid_weights_mod
