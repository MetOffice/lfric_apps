! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Decription:
! given an integer's data_start,data_end, splits this range into equal
! chunk_sizes for each thread. If the range is not divisible by the number
! of threads, chunk_size is increased by one for all the threads, starting
! with the last and counting backwards. Also returns start and end as indices
! into a shared array. Procedures can be made thread safe by calling this
! routine to determine the worksharing pattern. They can then be called from
! a parallel region with private variables.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Misc
!
module compute_chunk_size_mod

implicit none
contains

subroutine compute_chunk_size(data_start,data_end,iter_start,iter_end)
!$    use OMP_lib, only: omp_get_num_threads, omp_get_thread_num
! only use omp_lib if compiling with openMP
implicit none
integer, intent(in)  :: data_start!start of range to be split across threads
integer, intent(in)  :: data_end  !end of range to be split across threads
integer, intent(out) :: iter_start ! start for this thread
integer, intent(out) :: iter_end ! end for this thread
! local variables
integer :: thread_id, num_threads, extra
integer :: chunk_size
integer :: data_range  ! range of values from start to end

! outside of parallel region num_threads will be 1 if compiled without omp
! num_threads will be 1
num_threads = 1
!$     num_threads=omp_get_num_threads()

thread_id = 0
!$    thread_id=omp_get_thread_num()

data_range = (data_end - data_start) + 1
chunk_size = data_range/num_threads
extra = data_range - chunk_size*num_threads

iter_start = chunk_size*thread_id + data_start
iter_end = iter_start + chunk_size - 1

! account for data_range not divisible by num_threads
if ( (thread_id >= (num_threads-extra)) .and. (extra /= 0) ) then
  chunk_size = chunk_size + 1
  iter_start = iter_start + thread_id - (num_threads-extra)
  iter_end = iter_end + thread_id - (num_threads-extra) + 1
end if

return
end subroutine compute_chunk_size
end module compute_chunk_size_mod
