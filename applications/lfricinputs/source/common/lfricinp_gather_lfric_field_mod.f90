! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE lfricinp_gather_lfric_field_mod

! Intrinsic modules
USE, INTRINSIC :: iso_fortran_env, ONLY: real64, int32, int64

! external libraries
USE mpi

! lfric modules
USE field_mod, ONLY: field_type, field_proxy_type
USE constants_mod, ONLY: r_def, i_def
USE mesh_mod, ONLY: mesh_type
USE function_space_mod, ONLY: function_space_type
USE mpi_mod, ONLY: global_mpi
USE log_mod, ONLY: log_scratch_space, log_event, LOG_LEVEL_INFO, LOG_LEVEL_ERROR
IMPLICIT NONE

PRIVATE
PUBLIC :: lfricinp_gather_lfric_field

CONTAINS

SUBROUTINE lfricinp_gather_lfric_field( lfric_field, global_field_array, comm, &
                                        num_levels, level, twod_mesh )

IMPLICIT NONE
!
! Description:
!  Takes an partitioned lfric field and extracts data from a single level,
!  as specified in argument list, then gathers that data onto rank 0
!  and puts into correct location using the global id (gid) map
!
! Arguments
TYPE(field_type),    INTENT(INOUT) :: lfric_field
REAL(KIND=real64),   INTENT(OUT)   :: global_field_array(:)
INTEGER(KIND=i_def), INTENT(IN)    :: comm
INTEGER(KIND=int64), INTENT(IN)    :: num_levels
INTEGER(KIND=int64), INTENT(IN)    :: level
TYPE(mesh_type),     INTENT(IN), POINTER :: twod_mesh

! Local variables
TYPE(mesh_type), POINTER :: mesh => null()
TYPE(field_proxy_type) :: field_proxy
INTEGER(KIND=int32), ALLOCATABLE :: rank_sizes(:)
INTEGER(KIND=int32), ALLOCATABLE :: displacements(:)
INTEGER(KIND=int32) :: local_rank, total_ranks
INTEGER(KIND=int32) :: local_size_2d, global_size_2d
INTEGER(KIND=int32), PARAMETER   :: rank_0 = 0
INTEGER(KIND=int32) :: err, i
INTEGER(KIND=int64) :: index_3d
!, unit_num

REAL(KIND=real64), ALLOCATABLE :: local_data(:)
REAL(KIND=real64), ALLOCATABLE :: temp_global_data(:)

INTEGER(KIND=int32), ALLOCATABLE :: local_gid_lid_map(:)
INTEGER(KIND=int32), ALLOCATABLE :: global_gid_map(:)

! Get objects
mesh => lfric_field%get_mesh()

field_proxy = lfric_field%get_proxy()

! Get number of layers from function space, as we have W3, Wtheta and W3 2d
! fields which have different number of layers
!nlayers = fs%get_nlayers()
local_rank = global_mpi%get_comm_rank()
total_ranks = global_mpi%get_comm_size()

! The local size of single 2D level, just local domain, no haloes etc
local_size_2d = twod_mesh%get_last_edge_cell()
ALLOCATE(local_data(local_size_2d))

! Copy from 1D array that contains full local 3D field in column order into
! 1D array that contains only a 2D slice of the field
index_3d = level
DO i = 1, local_size_2d
  local_data(i) = field_proxy%data(index_3d)
  index_3d = index_3d + num_levels
END DO

! Gather size of each rank onto rank 0
ALLOCATE(rank_sizes(total_ranks))
CALL mpi_gather(local_size_2d, 1, mpi_integer, rank_sizes, 1, mpi_integer, &
               rank_0, comm, err)
IF (err /= mpi_success) THEN
  CALL log_event('Call to mpi_gather failed in MPI error.', &
       LOG_LEVEL_ERROR )
END IF

! Construct displacements array. This tells receiving rank the start position
! of where it should put the data from each rank.
ALLOCATE(displacements(total_ranks))
! Displacements value begin at zero, this correlates to first element
displacements(1) = 0
global_size_2d = rank_sizes(1)
DO i = 2, total_ranks
  ! Next position is previous position + size of previous buffer
  displacements(i) = displacements(i-1) + rank_sizes(i-1)
  global_size_2d = global_size_2d + rank_sizes(i)
END DO

IF (local_rank == rank_0) THEN
  ! Allocate full global array to receive data
  ALLOCATE(temp_global_data(global_size_2d))
ELSE
  ALLOCATE(temp_global_data(1))
END IF

! Gather data from all ranks onto rank 0
CALL mpi_gatherv(local_data, local_size_2d, mpi_double_precision,              &
     temp_global_data, rank_sizes, displacements, mpi_double_precision,        &
     rank_0, comm, err)
IF (err /= mpi_success) THEN
  CALL log_event('Call to mpi_gatherv failed in MPI error.', &
       LOG_LEVEL_ERROR )
END IF

ALLOCATE(local_gid_lid_map(local_size_2d))
! Get global indices for each local point on 2D level
DO i = 1, local_size_2d
  local_gid_lid_map(i) = mesh%get_gid_from_lid(i)
END DO

IF (local_rank == rank_0) THEN
  ALLOCATE(global_gid_map(global_size_2d))
ELSE
  ALLOCATE(global_gid_map(1))
END IF

! Gather gid map from all ranks onto rank 0
CALL mpi_gatherv(local_gid_lid_map, local_size_2d, mpi_integer, &
       global_gid_map, rank_sizes, displacements, mpi_integer, rank_0, &
       comm, err)
IF (err /= mpi_success) THEN
  CALL log_event('Call to mpi_gatherv failed in MPI error.', &
       LOG_LEVEL_ERROR )
END IF

IF (local_rank == rank_0) THEN
  IF (SIZE(global_field_array, 1) /= SIZE(temp_global_data, 1) ) THEN
    WRITE(log_scratch_space, '(2(A,I0))')                         &
         "Mismatch between array sizes global_field_array ",      &
         SIZE(global_field_array, 1), " and temp_global_data", &
         SIZE(temp_global_data, 1)
    CALL log_event(log_scratch_space, LOG_LEVEL_ERROR)
  END IF
  ! Use the gid map to copy data into correct location in main array
  DO i = 1, global_size_2d
    global_field_array(global_gid_map(i)) = temp_global_data(i)
  END DO
END IF

NULLIFY(mesh)

END SUBROUTINE lfricinp_gather_lfric_field

END MODULE lfricinp_gather_lfric_field_mod
