! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!> @brief   Module containing um2lfric configuration type
!> @details Hold information on input and output files and fields to be
!!          regridded, including proceedures to read this information from
!!          the relevant namelist.
module um2lfric_namelist_mod

! Intrinsic modules
use, intrinsic :: iso_fortran_env, only : int64, real64

! UM2LFRic modules
use lfricinp_um_parameters_mod,    only: fnamelen, um_imdi

implicit none
private
public :: um2lfric_config, required_lfric_namelists

type :: config
  character(len=fnamelen) :: um_file = 'unset'
  character(len=fnamelen) :: stashmaster_file = 'unset'
  character(len=fnamelen) :: weights_file_p_to_face_centre_bilinear = 'unset'
  character(len=fnamelen) :: weights_file_p_to_face_centre_neareststod = 'unset'
  character(len=fnamelen) :: weights_file_u_to_face_centre_bilinear = 'unset'
  character(len=fnamelen) :: weights_file_v_to_face_centre_bilinear = 'unset'
  integer(kind=int64), allocatable ::  stash_list(:)
  integer(kind=int64) :: num_snow_layers = um_imdi
  integer(kind=int64) :: num_surface_types = um_imdi
  integer(kind=int64) :: num_ice_cats = um_imdi

  integer :: num_fields

  integer :: status = -1
  character(len=512) :: message = 'No namelist read'
  integer :: unit_number

contains

  procedure :: load_namelist

end type config

integer(kind=int64), parameter :: max_stash_list = 999

! Input namelist configuration
type(config) :: um2lfric_config

! Namelist filenames read from command line
character(len=fnamelen), public :: um2lfric_nl_fname

character(*), parameter :: required_lfric_namelists(6) =  &
    ['logging             ', &
     'finite_element      ', &
     'base_mesh           ', &
     'planet              ', &
     'extrusion           ', &
     'io                  ']

contains

subroutine load_namelist(self)

  ! lfricinp modules
  use lfricinp_unit_handler_mod, only: get_free_unit
  use lfricinp_ancils_mod, only: lfricinp_l_land_area_fraction => &
                                 l_land_area_fraction

  ! LFRic modules
  use log_mod,          only: log_event, LOG_LEVEL_ERROR

  implicit none
  class(config) :: self

  ! Local variables
  integer :: i_stash

  ! Namelist variables
  character(len=fnamelen) :: um_file = 'unset'
  character(len=fnamelen) :: stashmaster_file = 'unset'
  character(len=fnamelen) :: weights_file_p_to_face_centre_bilinear = 'unset'
  character(len=fnamelen) :: weights_file_p_to_face_centre_neareststod = 'unset'
  character(len=fnamelen) :: weights_file_u_to_face_centre_bilinear = 'unset'
  character(len=fnamelen) :: weights_file_v_to_face_centre_bilinear = 'unset'
  integer(kind=int64)     :: stash_list(max_stash_list)
  integer(kind=int64)     :: num_snow_layers = um_imdi
  integer(kind=int64)     :: num_surface_types = um_imdi
  integer(kind=int64)     :: num_ice_cats = um_imdi
  logical :: l_land_area_fraction = .false.

  namelist /configure_um2lfric/ um_file,                                       &
                                stashmaster_file,                              &
                                weights_file_p_to_face_centre_bilinear,        &
                                weights_file_p_to_face_centre_neareststod,     &
                                weights_file_u_to_face_centre_bilinear,        &
                                weights_file_v_to_face_centre_bilinear,        &
                                stash_list, num_snow_layers, num_surface_types,&
                                num_ice_cats,                                  &
                                l_land_area_fraction

  stash_list(:) = um_imdi

  self%status = 0
  self%message = 'Reading namelist from ' // trim(um2lfric_nl_fname)

  call get_free_unit(self%unit_number)

  open(unit=self%unit_number, file=um2lfric_nl_fname, iostat=self%status,      &
                              iomsg=self%message)
  if (self%status /= 0) call log_event(self%message, LOG_LEVEL_ERROR)

  read(self%unit_number, nml=configure_um2lfric, iostat=self%status,           &
                         iomsg=self%message)
  if (self%status /= 0) call log_event(self%message, LOG_LEVEL_ERROR)

  if (trim(um_file) == 'unset') then
    self%status = 1
    self%message='UM filename is unset'
    call log_event(self%message, LOG_LEVEL_ERROR)
  end if

 ! Further error checking goes here

  ! Load namelist variables into object
  self%um_file = um_file
  self%stashmaster_file = stashmaster_file
  self%weights_file_p_to_face_centre_bilinear =                                &
                                    weights_file_p_to_face_centre_bilinear
  self%weights_file_p_to_face_centre_neareststod =                             &
                                    weights_file_p_to_face_centre_neareststod
  self%weights_file_u_to_face_centre_bilinear =                                &
                                    weights_file_u_to_face_centre_bilinear
  self%weights_file_v_to_face_centre_bilinear =                                &
                                    weights_file_v_to_face_centre_bilinear
  self%num_snow_layers = num_snow_layers
  self%num_surface_types = num_surface_types
  self%num_ice_cats = num_ice_cats
  ! Pass the um2lfric namelist variable to the lfricinputs module variable
  lfricinp_l_land_area_fraction = l_land_area_fraction

  self%num_fields=0
  ! Count how many fields have been requested
  do i_stash = 1, max_stash_list
    if (stash_list(i_stash) == um_imdi) then
      exit
    else
      self%num_fields = self%num_fields + 1
    end if
  end do

  if (self%num_fields <= 0) then
    call log_event('No fields selected in stash_list namelist variable', &
         LOG_LEVEL_ERROR)
  end if

  ! Can now allocate type variable
  allocate(self%stash_list(self%num_fields))
  self%stash_list(:) = stash_list(1:self%num_fields)

  self%status = 0
  self%message = 'Successfully read namelist from ' // trim(um2lfric_nl_fname)

  close(self%unit_number)

! Broadcasting?

end subroutine load_namelist

end module um2lfric_namelist_mod
