! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_gather_lfric_field_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only: real64, int32, int64

! external libraries
use mpi

! lfric modules
use field_mod, only: field_type, field_proxy_type
use constants_mod, only: r_def, i_def
use mesh_mod, only: mesh_type
use function_space_mod, only: function_space_type
use lfric_mpi_mod, only: global_mpi, lfric_comm_type
use log_mod, only: log_scratch_space, log_event, log_level_info, log_level_error
implicit none

private
public :: lfricinp_gather_lfric_field

contains

subroutine lfricinp_gather_lfric_field( lfric_field, global_field_array, comm, &
                                        num_levels, level, twod_mesh )

implicit none
!
! Description:
!  Takes an partitioned lfric field and extracts data from a single level,
!  as specified in argument list, then gathers that data onto rank 0
!  and puts into correct location using the global id (gid) map
!
! Arguments
type(field_type),    intent(INOUT) :: lfric_field
real(kind=real64),   intent(out)   :: global_field_array(:)
type(lfric_comm_type), intent(in)  :: comm
integer(kind=int64), intent(in)    :: num_levels
integer(kind=int64), intent(in)    :: level
type(mesh_type),     intent(in), pointer :: twod_mesh

! Local variables
type(mesh_type), pointer :: mesh => null()
type(field_proxy_type) :: field_proxy
integer(kind=int32), allocatable :: rank_sizes(:)
integer(kind=int32), allocatable :: displacements(:)
integer(kind=int32) :: local_rank, total_ranks
integer(kind=int32) :: local_size_2d, global_size_2d
integer(kind=int32), parameter   :: rank_0 = 0
integer(kind=int32) :: err, i
integer(kind=int64) :: index_3d
integer(kind=i_def) :: mpi_comm
!, unit_num

real(kind=real64), allocatable :: local_data(:)
real(kind=real64), allocatable :: temp_global_data(:)

integer(kind=int32), allocatable :: local_gid_lid_map(:)
integer(kind=int32), allocatable :: global_gid_map(:)

! Get objects
mesh => lfric_field%get_mesh()

field_proxy = lfric_field%get_proxy()

! Get number of layers from function space, as we have W3, Wtheta and W3 2d
! fields which have different number of layers
!nlayers = fs%get_nlayers()
local_rank = global_mpi%get_comm_rank()
total_ranks = global_mpi%get_comm_size()

! The local size of single 2D level, just local domain, no halos etc
local_size_2d = twod_mesh%get_last_edge_cell()
allocate(local_data(local_size_2d))

! Copy from 1D array that contains full local 3D field in column order into
! 1D array that contains only a 2D slice of the field
index_3d = level
do i = 1, local_size_2d
  local_data(i) = field_proxy%data(index_3d)
  index_3d = index_3d + num_levels
end do

mpi_comm = comm%get_comm_mpi_val()
! Gather size of each rank onto rank 0
allocate(rank_sizes(total_ranks))
call mpi_gather(local_size_2d, 1, mpi_integer, rank_sizes, 1, mpi_integer, &
               rank_0, mpi_comm, err)
if (err /= mpi_success) then
  call log_event('Call to mpi_gather failed in MPI error.', &
       LOG_LEVEL_ERROR )
end if

! Construct displacements array. This tells receiving rank the start position
! of where it should put the data from each rank.
allocate(displacements(total_ranks))
! Displacements value begin at zero, this correlates to first element
displacements(1) = 0
global_size_2d = rank_sizes(1)
do i = 2, total_ranks
  ! Next position is previous position + size of previous buffer
  displacements(i) = displacements(i-1) + rank_sizes(i-1)
  global_size_2d = global_size_2d + rank_sizes(i)
end do

if (local_rank == rank_0) then
  ! Allocate full global array to receive data
  allocate(temp_global_data(global_size_2d))
else
  allocate(temp_global_data(1))
end if

! Gather data from all ranks onto rank 0
call mpi_gatherv(local_data, local_size_2d, mpi_double_precision,              &
     temp_global_data, rank_sizes, displacements, mpi_double_precision,        &
     rank_0, mpi_comm, err)
if (err /= mpi_success) then
  call log_event('Call to mpi_gatherv failed in MPI error.', &
       LOG_LEVEL_ERROR )
end if

allocate(local_gid_lid_map(local_size_2d))
! Get global indices for each local point on 2D level
do i = 1, local_size_2d
  local_gid_lid_map(i) = mesh%get_gid_from_lid(i)
end do

if (local_rank == rank_0) then
  allocate(global_gid_map(global_size_2d))
else
  allocate(global_gid_map(1))
end if

! Gather gid map from all ranks onto rank 0
call mpi_gatherv(local_gid_lid_map, local_size_2d, mpi_integer, &
       global_gid_map, rank_sizes, displacements, mpi_integer, rank_0, &
       mpi_comm, err)
if (err /= mpi_success) then
  call log_event('Call to mpi_gatherv failed in MPI error.', &
       LOG_LEVEL_ERROR )
end if

if (local_rank == rank_0) then
  if (size(global_field_array, 1) /= size(temp_global_data, 1) ) then
    write(log_scratch_space, '(2(A,I0))')                         &
         "Mismatch between array sizes global_field_array ",      &
         size(global_field_array, 1), " and temp_global_data", &
         size(temp_global_data, 1)
    call log_event(log_scratch_space, LOG_LEVEL_ERROR)
  end if
  ! Use the gid map to copy data into correct location in main array
  do i = 1, global_size_2d
    global_field_array(global_gid_map(i)) = temp_global_data(i)
  end do
end if

nullify(mesh)

end subroutine lfricinp_gather_lfric_field

end module lfricinp_gather_lfric_field_mod
