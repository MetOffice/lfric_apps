!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Adjoint of the sci_psykal_builtin_light_mod routines
module adj_sci_psykal_builtin_light_mod

  use, intrinsic :: iso_fortran_env, only : real32, real64
  use constants_mod,                 only: r_def, i_def
  use field_mod,                     only: field_type, field_proxy_type

  implicit none

  private
  public :: invoke_adj_copy_field_32_64
  public :: invoke_adj_copy_field_32_32
  public :: invoke_adj_copy_field_64_64

  contains

  !=============================================================================
  ! COPY ROUTINES
  !=============================================================================
  ! Adjoint copy routines necessary due to needing adjoint of pre-existing PSyKAl-lite code,
  ! which will be fixed in PSyclone issue #2674. We may potentially need our own PSyclone
  ! issue which addresses the adjoint of said builtins.
  !=============================================================================
  !> @brief Adjoint copy routine for real32 to real64 field
  !> @param[in,out] fsrce_32  Source real32 field
  !> @param[in,out] fdest_64  Destination real64 field
  subroutine invoke_adj_copy_field_32_64(fsrce_32, fdest_64)

     use omp_lib,            only: omp_get_thread_num
     use omp_lib,            only: omp_get_max_threads
     use mesh_mod,           only: mesh_type
     use field_real32_mod,   only: field_real32_type, &
                                   field_real32_proxy_type
     use field_real64_mod,   only: field_real64_type, &
                                   field_real64_proxy_type

     implicit none

     type(field_real32_type), intent(inout)  :: fsrce_32
     type(field_real64_type), intent(inout)  :: fdest_64

     integer(kind=i_def)             :: df
     integer(kind=i_def)             :: loop0_start, loop0_stop
     type(field_real32_proxy_type)   :: fsrce_32_proxy
     type(field_real64_proxy_type)   :: fdest_64_proxy
     integer(kind=i_def)             :: max_halo_depth_mesh
     type(mesh_type), pointer        :: mesh => null()
     !
     ! Initialise field and/or operator proxies
     !
     fsrce_32_proxy = fsrce_32%get_proxy()
     fdest_64_proxy = fdest_64%get_proxy()
     !
     ! Create a mesh object
     !
     mesh => fdest_64_proxy%vspace%get_mesh()
     max_halo_depth_mesh = mesh%get_halo_depth()
     !
     ! Set-up all of the loop bounds
     !
     loop0_start = 1
     if (fsrce_32_proxy%is_dirty(depth=1)) then
       ! only copy the owned dofs
       loop0_stop = fdest_64_proxy%vspace%get_last_dof_annexed()
     else
       ! copy the 1st halo row as well
       loop0_stop = fdest_64_proxy%vspace%get_last_dof_halo(1)
     end if
     !
     ! Call kernels and communication routines
     !
     !$omp parallel default(shared), private(df)
     !$omp do schedule(static)
     do df=loop0_start,loop0_stop
       fsrce_32_proxy%data(df) = fsrce_32_proxy%data(df) + real(fdest_64_proxy%data(df), real32)
       fdest_64_proxy%data(df) = 0.0_real64
     end do
     !$omp end do
     !$omp end parallel
     !
     ! Set halos dirty/clean for fields modified in the above loop
     !
     call fdest_64_proxy%set_dirty()
     if (.not. fsrce_32_proxy%is_dirty(depth=1)) then
       call fdest_64_proxy%set_clean(1)
     end if
     !

  end subroutine invoke_adj_copy_field_32_64

  !=============================================================================
  !> @brief Adjoint copy routine for real32 to real32 field
  !> @param[in,out] fsrce_32  Source real32 field
  !> @param[in,out] fdest_32  Destination real32 field
  subroutine invoke_adj_copy_field_32_32(fsrce_32, fdest_32)

     use omp_lib,            only: omp_get_thread_num
     use omp_lib,            only: omp_get_max_threads
     use mesh_mod,           only: mesh_type
     use field_real32_mod,   only: field_real32_type, &
                                   field_real32_proxy_type

     implicit none

     type(field_real32_type), intent(inout)  :: fsrce_32
     type(field_real32_type), intent(inout)  :: fdest_32

     integer(kind=i_def)             :: df
     integer(kind=i_def)             :: loop0_start, loop0_stop
     type(field_real32_proxy_type)   :: fsrce_32_proxy
     type(field_real32_proxy_type)   :: fdest_32_proxy
     integer(kind=i_def)             :: max_halo_depth_mesh
     type(mesh_type), pointer        :: mesh => null()
     !
     ! Initialise field and/or operator proxies
     !
     fsrce_32_proxy = fsrce_32%get_proxy()
     fdest_32_proxy = fdest_32%get_proxy()
     !
     ! Create a mesh object
     !
     mesh => fdest_32_proxy%vspace%get_mesh()
     max_halo_depth_mesh = mesh%get_halo_depth()
     !
     ! Set-up all of the loop bounds
     !
     loop0_start = 1
     if (fsrce_32_proxy%is_dirty(depth=1)) then
       ! only copy the owned dofs
       loop0_stop = fdest_32_proxy%vspace%get_last_dof_annexed()
     else
       ! copy the 1st halo row as well
       loop0_stop = fdest_32_proxy%vspace%get_last_dof_halo(1)
     end if
     !
     ! Call kernels and communication routines
     !
     !$omp parallel default(shared), private(df)
     !$omp do schedule(static)
     do df=loop0_start,loop0_stop
       fsrce_32_proxy%data(df) = fsrce_32_proxy%data(df) + real(fdest_32_proxy%data(df), real32)
       fdest_32_proxy%data(df) = 0.0_real32
     end do
     !$omp end do
     !$omp end parallel
     !
     ! Set halos dirty/clean for fields modified in the above loop
     !
     call fdest_32_proxy%set_dirty()
     if (.not. fsrce_32_proxy%is_dirty(depth=1)) then
       call fdest_32_proxy%set_clean(1)
     end if
     !

  end subroutine invoke_adj_copy_field_32_32

  !=============================================================================
  !> @brief Adjoint copy routine for real64 to real64 field
  !> @param[in,out] fsrce_64  Source real64 field
  !> @param[in,out] fdest_64  Destination real64 field
  subroutine invoke_adj_copy_field_64_64(fsrce_64, fdest_64)

     use omp_lib,            only: omp_get_thread_num
     use omp_lib,            only: omp_get_max_threads
     use mesh_mod,           only: mesh_type
     use field_real64_mod,   only: field_real64_type, &
                                   field_real64_proxy_type

     implicit none

     type(field_real64_type), intent(in)     :: fsrce_64
     type(field_real64_type), intent(inout)  :: fdest_64

     integer(kind=i_def)             :: df
     integer(kind=i_def)             :: loop0_start, loop0_stop
     type(field_real64_proxy_type)   :: fsrce_64_proxy
     type(field_real64_proxy_type)   :: fdest_64_proxy
     integer(kind=i_def)             :: max_halo_depth_mesh
     type(mesh_type), pointer        :: mesh => null()
     !
     ! Initialise field and/or operator proxies
     !
     fsrce_64_proxy = fsrce_64%get_proxy()
     fdest_64_proxy = fdest_64%get_proxy()
     !
     ! Create a mesh object
     !
     mesh => fdest_64_proxy%vspace%get_mesh()
     max_halo_depth_mesh = mesh%get_halo_depth()
     !
     ! Set-up all of the loop bounds
     !
     loop0_start = 1
     if (fsrce_64_proxy%is_dirty(depth=1)) then
       ! only copy the owned dofs
       loop0_stop = fdest_64_proxy%vspace%get_last_dof_annexed()
     else
       ! copy the 1st halo row as well
       loop0_stop = fdest_64_proxy%vspace%get_last_dof_halo(1)
     end if
     !
     ! Call kernels and communication routines
     !
     !$omp parallel default(shared), private(df)
     !$omp do schedule(static)
     do df=loop0_start,loop0_stop
       fsrce_64_proxy%data(df) = fsrce_64_proxy%data(df) + real(fdest_64_proxy%data(df), real64)
       fdest_64_proxy%data(df) = 0.0_real64
     end do
     !$omp end do
     !$omp end parallel
     !
     ! Set halos dirty/clean for fields modified in the above loop
     !
     call fdest_64_proxy%set_dirty()
     if (.not. fsrce_64_proxy%is_dirty(depth=1)) then
       call fdest_64_proxy%set_clean(1)
     end if
     !
  end subroutine invoke_adj_copy_field_64_64

end module adj_sci_psykal_builtin_light_mod
