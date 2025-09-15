! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
module lfricinp_masks_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only: int32, int64, real64
use, intrinsic :: iso_c_binding, only: c_bool

implicit none

logical, allocatable :: um_land_mask(:,:), um_maritime_mask(:,:)
logical, allocatable :: lfric_land_mask(:), lfric_maritime_mask(:)

contains

!> @brief Initialises UM and LFRic land and maritime masks
!> @param[in] stashcode_land_mask Stashcode for the UM land fraction
subroutine lfricinp_init_masks(stashcode_land_mask)

use log_mod,        only: log_event, log_scratch_space, LOG_LEVEL_INFO
use field_mod,      only: lfric_field_type => field_type,                      &
                          lfric_proxy_type => field_proxy_type
use local_mesh_mod, only: local_mesh_type
use mesh_mod,       only: mesh_type

use lfricinp_lfric_driver_mod,         only: local_rank
use lfricinp_check_shumlib_status_mod, only: shumlib
use um2lfric_read_um_file_mod,        only: um_input_file
use lfricinp_ancils_mod,               only: ancil_fields, l_land_area_fraction
!
! shumlib modules
use f_shum_field_mod, only: shum_field_type

implicit none

! LFRic mesh and local mesh
type(mesh_type), pointer :: mesh
type(local_mesh_type), pointer :: local_mesh

! Ancil fields
type(lfric_field_type), pointer :: ancil_field
type(lfric_proxy_type) :: ancil_field_proxy

! Array of shumlib field objects that will be returned from UM file
type(shum_field_type), allocatable  :: um_input_fields(:)

! Local variables
integer(kind=int64)  :: stashcode_land_mask
! Results of size() are the default integer kind
integer :: dim_1d, dim_2dx, dim_2dy, i
integer(kind=int32)  :: err
logical(kind=c_bool) :: true_cbool

true_cbool = logical(.true., kind=c_bool)

! Get UM land mask
call shumlib("um2lfric::find_fields_in_file",                                  &
             um_input_file%find_fields_in_file(um_input_fields,                &
                                               stashcode = stashcode_land_mask,&
                                               lbproc = 0_int64),              &
             ignore_warning=true_cbool, errorstatus=err)

! Check that both the source and destination land area fraction fields are
! available. Source from start dump via shumlib and destination from lfric ancil
if (err /= 0 .or. .not. l_land_area_fraction) then

  if (err /=0 ) then
    log_scratch_space = 'Error encountered when trying to read UM land mask. '
  else if (.not. l_land_area_fraction) then
    log_scratch_space = 'l_land_area_fraction set to false. Land fraction '//  &
         'ancil not available. '
  end if
  log_scratch_space = trim(log_scratch_space) //                               &
                      'Post regridding adjustments '//                         &
                      'on land and maritime compressed fields will be ignored.'
  call log_event(log_scratch_space, LOG_LEVEL_INFO)

else

  ! Set up LFRic mask dimension size
  call ancil_fields%get_field("land_area_fraction", ancil_field)
  call ancil_field%read_field(ancil_field%get_name())
  ancil_field_proxy = ancil_field%get_proxy()
  dim_1d = size(ancil_field_proxy%data)
  !
  ! Get local_mesh LFRic mask is defined on
  mesh => ancil_field%get_mesh()
  local_mesh => mesh%get_local_mesh()
  !
  ! Set up LFRic logical land mask
  allocate(lfric_land_mask(dim_1d))
  do i = 1, dim_1d
    lfric_land_mask(i) = .false.
    if (local_mesh%get_cell_owner(i) == local_rank ) then
      if (ancil_field_proxy%data(i) > 0.0_real64) then
        lfric_land_mask(i) = .true.
      end if
    end if
  end do
  !
  ! Set up LFRic logical maritime mask
  allocate(lfric_maritime_mask(dim_1d))
  do i = 1, dim_1d
    lfric_maritime_mask(i) = .false.
    if (local_mesh%get_cell_owner(i) == local_rank ) then
      if (ancil_field_proxy%data(i) < 1.0_real64) then
        lfric_maritime_mask(i) = .true.
      end if
    end if
  end do
  !
  ! Nullify LFRic mesh and local_mesh pointers
  nullify(mesh, local_mesh)

  ! Set up UM mask dimensions
  dim_2dx = size(um_input_fields(1)%rdata,dim=1)
  dim_2dy = size(um_input_fields(1)%rdata,dim=2)
  !
  ! Set up UM logical land mask
  allocate(um_land_mask(dim_2dx,dim_2dy))
  um_land_mask = (um_input_fields(1)%rdata > 0.0_real64)
  !
  ! Set up UM logical maritime mask
  allocate(um_maritime_mask(dim_2dx,dim_2dy))
  um_maritime_mask = (um_input_fields(1)%rdata < 1.0_real64)

end if

end subroutine lfricinp_init_masks


subroutine lfricinp_finalise_masks()

implicit none

if (allocated(um_land_mask)) deallocate(um_land_mask)
if (allocated(lfric_land_mask)) deallocate(lfric_land_mask)
if (allocated(um_maritime_mask)) deallocate(um_maritime_mask)
if (allocated(lfric_maritime_mask)) deallocate(lfric_maritime_mask)

end subroutine lfricinp_finalise_masks

end module lfricinp_masks_mod
