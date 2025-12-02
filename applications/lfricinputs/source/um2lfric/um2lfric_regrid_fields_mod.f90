! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module um2lfric_regrid_fields_mod


implicit none

private

public :: um2lfric_regrid_fields

contains

subroutine um2lfric_regrid_fields(fctime)

! Description:
!
!  This routine performs a level by level regridding of the UM fields to
!  an intermediate regridded field array for each field, which is then
!  unpacked/copied into the corresponding lfric field data array as follows:
!
!  The LFRic field is a single 1D array containing all levels where the field
!  dof ordering in the LFRic array loops over each column before moving onto
!  the next horizontal point.
!
!  In this current UM2LFRic code we use a serial implementation of the LFRic
!  infrastructure. This means that there are no complications from the
!  MPI partitioning or from field stencils.
!
!  Simple example of a 2 layer mesh, with 4 points on each layer and both
!  the target and LFRic field indicies labelled for each point
!
!  Note that the target array indices repeat for each layer but the LFRic
!  field index refers to the full 3D field

!             Layer 1                                  Layer 2
!
!
! x  regridded_1      x  regridded_3       x  regridded_1       x  regridded_3
!    LFRIC_1             LFRIC_5              LFRIC_2              LFRIC_6
!
! x  regridded_2      x  regridded_4       x  regridded_2       x  regridded_4
!    LFRIC_3             LFRIC_7              LFRIC_4              LFRIC_8

! 2D array field containing the regridded data. First dimension corresponds
! to number of points per level and the second dimension to the number levels

! Intrinsic modules
use, intrinsic :: iso_fortran_env, &
                                  only: real64, int32, int64
use, intrinsic :: iso_c_binding,  only: c_bool

! lfricinputs modules
use lfricinp_check_shumlib_status_mod, &
                                  only: shumlib
use lfricinp_regrid_weights_type_mod, &
                                  only: lfricinp_regrid_weights_type
use lfricinp_stash_to_lfric_map_mod, &
                                  only: get_field_name
use lfricinp_lfric_driver_mod,    only: lfric_fields
use lfricinp_regrid_options_mod,  only: regrid_type

! um2lfric modules
use um2lfric_namelist_mod,        only: um2lfric_config
use um2lfric_read_um_file_mod,    only: um_input_file
use um2lfric_regrid_weights_mod,  only: get_weights
use um2lfric_post_process_fields_mod, &
                                  only: um2lfric_post_process_fields
use um2lfric_apply_masked_field_adjustments_mod, &
                                  only: um2lfric_apply_masked_field_adjustments

! shumlib modules
use f_shum_field_mod,             only: shum_field_type

! lfric modules
use field_mod,                    only: lfric_field_type => field_type, &
                                        lfric_proxy_type => field_proxy_type
use function_space_mod,           only: function_space_type
use fs_continuity_mod,            only: W3, Wtheta, W2H
use log_mod,                      only: log_event,       &
                                        LOG_LEVEL_INFO,  &
                                        LOG_LEVEL_ERROR, &
                                        log_scratch_space

implicit none

! Forecast time of a field
real(kind=real64), intent(in) :: fctime

! Array of shumlib field objects that will be returned from UM file
type(shum_field_type), allocatable  :: um_input_fields(:)

! Intermediate target field for regridding
real(kind=real64), allocatable :: regridded_field(:,:)

! Pointers to lfric objects
type(lfric_field_type), pointer :: lfric_field => null()
type(function_space_type), pointer :: lfric_field_fs => null()
type(lfricinp_regrid_weights_type), pointer :: weights => null()

! LFRic field proxy
type(lfric_proxy_type) :: lfric_field_proxy

! Other variables
integer(kind=int64) :: stashcode
integer(kind=int32) :: fs_type
! Iterators
integer :: i_field, level, regridded_index, len_regridded_field, lfric_index
integer :: errorstatus
! Number levels of UM field
integer(kind=int32) :: num_levels, num_dofs_per_level, num_lfric_levels
logical(kind=c_bool) :: true_cbool

true_cbool = logical(.true., kind=c_bool)

!-------------------------------------------------------------------------------
! Initialise lfric fields to zero
!-------------------------------------------------------------------------------
do i_field = 1, um2lfric_config%num_fields
  stashcode = um2lfric_config%stash_list(i_field)
  call lfric_fields % get_field(get_field_name(stashcode), lfric_field)
  lfric_field_proxy = lfric_field % get_proxy()
  lfric_field_proxy % data(:) = 0.0_real64
  lfric_field => null()
end do

!-------------------------------------------------------------------------------
! Main loop over fields for regridding
!-------------------------------------------------------------------------------
write(log_scratch_space, '(A,I0,A)') 'Will process ', &
     um2lfric_config%num_fields, ' fields'
call log_event(log_scratch_space, LOG_LEVEL_INFO)

do i_field = 1, um2lfric_config%num_fields

  stashcode = um2lfric_config%stash_list(i_field)
  write(log_scratch_space, '(A,I0)') 'Processing/regrid STASH code: ',         &
                                     stashcode
  call log_event(log_scratch_space, LOG_LEVEL_INFO)
  write(log_scratch_space,'(A,A)') 'LFRic field name: ',                       &
                                   trim(get_field_name(stashcode))
  call log_event(log_scratch_space, LOG_LEVEL_INFO)

  !-----------------------------------------------------------------------------
  ! Get UM field array and number of UM field levels
  !-----------------------------------------------------------------------------
  call shumlib("um2lfric::find_fields_in_file",                                &
                um_input_file%find_fields_in_file(um_input_fields,             &
                stashcode = stashcode, lbproc = 0_int64,                       &
                fctime = fctime), ignore_warning = true_cbool,                 &
                errorstatus = errorstatus)

  if (errorstatus /= 0) then ! Field has not been found in dump

    write(log_scratch_space,'(A,I0,A)') 'WARNING: stashcode ', stashcode,      &
                                        'not found in input dump. Data set '// &
                                        'to zero in LFRic output field.'
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
    if (allocated(um_input_fields)) deallocate(um_input_fields)

  else ! Field has sucessfully been found

    num_levels = size(um_input_fields)

    !---------------------------------------------------------------------------
    ! Create pointers to lfric field and function space field lives on
    !---------------------------------------------------------------------------
    call lfric_fields % get_field(get_field_name(stashcode), lfric_field)
    lfric_field_fs => lfric_field % get_function_space()

    !---------------------------------------------------------------------------
    ! Allocate and set dimensions of intermediate regridded field data array
    !---------------------------------------------------------------------------
    fs_type = lfric_field % which_function_space()
    if ((fs_type == W3) .or.  (fs_type == W2H)) then
      num_lfric_levels = lfric_field_fs % get_nlayers()
    else if (fs_type == Wtheta) then
      num_lfric_levels = lfric_field_fs % get_nlayers() + 1
    else
      write(log_scratch_space,'(A)') 'Function space for field ' //            &
                                     get_field_name(stashcode) //              &
                                     ' not currently supported in regridding'
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if
    num_dofs_per_level = lfric_field_fs % get_undf() /                         &
                         (num_lfric_levels * lfric_field_fs % get_ndata())
    allocate(regridded_field(num_dofs_per_level, num_levels))

    !---------------------------------------------------------------------------
    ! Get required regridding weights
    !---------------------------------------------------------------------------
    weights => get_weights(stashcode)

    !---------------------------------------------------------------------------
    ! Loop for level by level regridding
    !---------------------------------------------------------------------------
    do level = 1, num_levels

      !-------------------------------------------------------------------------
      ! Perform regridding/copying data from one grid/mesh to another
      !-------------------------------------------------------------------------
      call weights % regrid_src_2d_dst_1d (src=um_input_fields(level)%rdata,   &
                                           dst=regridded_field(:, level))
      !-------------------------------------------------------------------------
      ! Perform post regridding masked field adjustments, if regridding was not
      ! a simple copy of data (as is case for LAMs)
      !-------------------------------------------------------------------------
      if (trim(regrid_type) /= 'lam_to_lam') then
        call um2lfric_apply_masked_field_adjustments(                          &
                                            stashcode,                         &
                                            src=um_input_fields(level)%rdata,  &
                                            dst=regridded_field(:, level))
      end if

    end do ! loop over levels

    ! Tidy up input field memory
    if (allocated(um_input_fields)) deallocate(um_input_fields)

    !---------------------------------------------------------------------------
    ! Do any final post-processing to field, if required
    !---------------------------------------------------------------------------
    call um2lfric_post_process_fields(regridded_field, stashcode)
    ! Update number of levels for field if it changed during post-processing
    num_levels = size(regridded_field, 2)

    !---------------------------------------------------------------------------
    ! Copy the regridded data to the lfric field
    !---------------------------------------------------------------------------
    lfric_field_proxy = lfric_field % get_proxy()
    if (size(regridded_field) /= size(lfric_field_proxy % data)) then
      write(log_scratch_space,'(A,I0,A,I0)')                                   &
                              'Mismatch between lfric field data array size ', &
                              size(lfric_field_proxy % data), ' and ' //       &
                              'regridded array size ',                         &
                              size(regridded_field)
      call log_event(log_scratch_space, LOG_LEVEL_ERROR)
    end if
    write(log_scratch_space,'(A)') 'Fill LFRic field object with regridded data'
    call log_event(log_scratch_space, LOG_LEVEL_INFO)
    do level = 1, num_levels
       len_regridded_field = size(regridded_field, 1)
       ! Split up the data here and insert it into the LFRic field in the
       ! correct place. Loop over all points in the 2D regridded field that
       ! represents a single level. The first lfric array index will match the
       ! current level number
       lfric_index = level
       do regridded_index = 1, len_regridded_field
         lfric_field_proxy % data(lfric_index) =                               &
                                      lfric_field_proxy % data(lfric_index) +  &
                                      regridded_field(regridded_index, level)
         ! Need to step by the total number of levels to get to the lfric
         ! array index that corresponds to the next regridded point
         lfric_index = lfric_index + num_levels
       end do
     end do

    deallocate(regridded_field)

    !---------------------------------------------------------------------------
    ! Nullify pointers before redefining
    !---------------------------------------------------------------------------
    weights => null()
    lfric_field => null()
    lfric_field_fs => null()

  end if

end do ! loop over stashcodes

end subroutine um2lfric_regrid_fields

end module um2lfric_regrid_fields_mod
