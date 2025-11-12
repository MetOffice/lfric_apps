! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfric2um_main_loop_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only: int64, real64

! lfric2um modules
use lfric2um_namelists_mod,       only: lfric2um_config
use lfricinp_lfric_driver_mod,    only: lfric_fields, local_rank, comm,   &
                                        twod_mesh
use lfric2um_initialise_um_mod,   only: um_output_file
use lfricinp_um_grid_mod,         only: um_grid
use lfric2um_regrid_weights_mod,  only: get_weights
use lfric2um_exner_above_top_mod, only: lfric2um_exner_above_top
use lfric2um_conv_exner_mod,      only: lfric2um_conv_p_exner,            &
                                        lfric2um_conv_exner_p

! lfricinp modules
use lfricinp_stashmaster_mod,          only: get_stashmaster_item, levelt,     &
                                             rho_levels, single_level,         &
                                             stashcode_exner, stashcode_p,     &
                                             stashcode_q, stashcode_theta,     &
                                             theta_levels
use lfricinp_add_um_field_to_file_mod, only: lfricinp_add_um_field_to_file
use lfricinp_um_level_codes_mod,       only: lfricinp_get_num_levels
use lfricinp_check_shumlib_status_mod, only: shumlib
use lfricinp_regrid_weights_type_mod,  only: lfricinp_regrid_weights_type
use lfricinp_stash_to_lfric_map_mod,   only: get_field_name
use lfricinp_gather_lfric_field_mod,   only: lfricinp_gather_lfric_field

! lfric modules
use field_mod,  only: field_type
use log_mod,    only: log_event, log_scratch_space, LOG_LEVEL_INFO

implicit none

private
public :: lfric2um_main_loop

contains

subroutine lfric2um_main_loop()
! Description:
!  Main processing field loop in lfric2um. Loops over requested fields,
!  reads field on all ranks, loops over levels, gather full level data
!  onto rank 0, perform regridding and write to disk.
implicit none

integer(kind=int64) :: i_stash, level, i_field
integer(kind=int64) :: stashcode, num_levels
character(len=*), parameter :: routinename='lfric2um_main_loop'
type(field_type), pointer :: lfric_field
type(lfricinp_regrid_weights_type), pointer :: weights
real(kind=real64), allocatable :: global_field_array(:)
real(kind=real64), allocatable :: q_top_buffer(:,:)      ! Buffer arrays for
real(kind=real64), allocatable :: theta_top_buffer(:,:)  ! computing Exner pressure
real(kind=real64), allocatable :: exner_top_buffer(:,:)  ! above the top level
logical :: l_conv_p_exner

if (local_rank == 0) then
  allocate( q_top_buffer(um_grid%num_p_points_x,um_grid%num_p_points_y) )
  allocate( theta_top_buffer(um_grid%num_p_points_x,um_grid%num_p_points_y) )
  allocate( exner_top_buffer(um_grid%num_p_points_x,um_grid%num_p_points_y) )
end if

! Main loop over requested stashcodes
do i_stash = 1, lfric2um_config%num_fields
  stashcode = lfric2um_config%stash_list(i_stash)
  write(log_scratch_space, '(A,I0)') "Processing stashcode ", stashcode
  call log_event(log_scratch_space, LOG_LEVEL_INFO)
  num_levels = lfricinp_get_num_levels(um_output_file, stashcode)

  ! Select appropriate weights
  weights => get_weights(stashcode)

  ! Get pointer to lfric field + read
  call lfric_fields%get_field(get_field_name(stashcode), lfric_field)
  call lfric_field%read_field("read_"//lfric_field%get_name())

  ! Allocate space for global data, only need full field on rank 0
  if (allocated(global_field_array)) deallocate(global_field_array)
  if (local_rank == 0) then
    allocate(global_field_array(weights%num_points_src))
  else
    allocate(global_field_array(1))
  end if

  l_conv_p_exner = .false.

  ! Loop over number of levels in field
  do level = 1, num_levels
    global_field_array(:) = 0.0_real64
    ! Gather local lfric fields into global array on base rank
    call lfricinp_gather_lfric_field( lfric_field, global_field_array, comm, &
         num_levels, level, twod_mesh )
    if (local_rank == 0 ) then
      ! Create UM field in output file
      call lfricinp_add_um_field_to_file(um_output_file, stashcode, &
           level, um_grid, lfric2um_config%lbtim_list(i_stash),     &
           lfric2um_config%lbproc_list(i_stash))
      ! Adding a 2D UM field to the file increments num_fields by 1 each time
      ! Get the index of the last field added, which is just num_fields
      i_field = um_output_file%num_fields
      call weights%validate_dst(size(um_output_file%fields(i_field)%rdata))
      ! Perform the regridding from lfric to um
      call weights%regrid_src_1d_dst_2d(global_field_array(:), &
           um_output_file%fields(i_field)%rdata(:,:))

      ! Copy the pointers for the top-level fields needed to compute the
      ! Exner pressure at the half level immediately above the model top
      if (stashcode == stashcode_q .and. level == num_levels) then
        q_top_buffer(:,:) = um_output_file%fields(i_field)%rdata(:,:)
      else if (stashcode == stashcode_theta .and. level == num_levels) then
        theta_top_buffer(:,:) = um_output_file%fields(i_field)%rdata(:,:)
      else if (stashcode == stashcode_exner .and. level == num_levels) then
        exner_top_buffer(:,:) = um_output_file%fields(i_field)%rdata(:,:)
        l_conv_p_exner = .false.
      else if (stashcode == stashcode_p .and. level == num_levels) then
        exner_top_buffer(:,:) = um_output_file%fields(i_field)%rdata(:,:)
        l_conv_p_exner = .true.
      end if

      ! Write to file
      call shumlib(routinename//'::write_field',  &
           um_output_file%write_field(i_field))
      ! Unload data from field
      call shumlib(routinename//'::unload_field', &
           um_output_file%unload_field(i_field))
    end if
  end do ! end loop over levels

  ! Add the additional upper level for the Exner pressure
  if (stashcode == stashcode_exner .or. stashcode == stashcode_p) then
    level = num_levels + 1
    if (local_rank == 0) then
      call lfricinp_add_um_field_to_file(um_output_file, stashcode, &
           level, um_grid, lfric2um_config%lbtim_list(i_stash),     &
           lfric2um_config%lbproc_list(i_stash))
      i_field = um_output_file%num_fields
      if (l_conv_p_exner) then
        call lfric2um_conv_p_exner(exner_top_buffer)
      end if
      call lfric2um_exner_above_top(theta_top_buffer, exner_top_buffer, q_top_buffer, &
           um_output_file%fields(i_field)%rdata(:,:))
      if (l_conv_p_exner) then
        call lfric2um_conv_exner_p(um_output_file%fields(i_field)%rdata(:,:))
      end if
      call shumlib(routinename//'::write_field',  &
           um_output_file%write_field(i_field))
      call shumlib(routinename//'::unload_field', &
           um_output_file%unload_field(i_field))
    end if
  end if

  ! Unload field data from memory
  call lfric_field%field_final()
end do  ! end loop over field stashcodes

if (local_rank == 0) then ! Only write UM file with rank 0
  ! Write output header
  call shumlib(routinename//'::write_header', &
           um_output_file%write_header())
  ! Close file
  call shumlib(routinename//'::close_file', &
           um_output_file%close_file())

  deallocate( q_top_buffer )
  deallocate( theta_top_buffer )
  deallocate( exner_top_buffer )

end if

end subroutine lfric2um_main_loop

end module lfric2um_main_loop_mod
