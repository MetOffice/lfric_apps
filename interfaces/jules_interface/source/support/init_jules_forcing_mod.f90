!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Create fields holding JULES site forcing data read from file
module init_jules_forcing_mod

  use constants_mod,            only : i_def, str_def
  use field_collection_mod,     only : field_collection_type
  use gungho_time_axes_mod,     only : gungho_time_axes_type
  use init_ancils_mod,          only : setup_ancil_field
  use log_mod,                  only : log_event, LOG_LEVEL_INFO, &
                                       LOG_LEVEL_ERROR
  use mesh_mod,                 only : mesh_type

  implicit none

  private
  public :: create_jules_forcing_fields, jules_forcing_field_names

  !> Names of the model fields populated from the forcing file. The XIOS
  !> field ids in the forcing file definition are these with "_data" appended.
  character(len=str_def), parameter :: jules_forcing_field_names(7) = &
    [ character(len=str_def) :: 'forcing_wind',   & ! 10 m wind speed (m s-1)
                                'forcing_tair',   & ! 10 m air temperature (K)
                                'forcing_qair',   & ! 10 m specific humidity (kg kg-1)
                                'forcing_psurf',  & ! Surface pressure (Pa)
                                'forcing_lwdown', & ! Downward LW flux (W m-2)
                                'forcing_swdown', & ! Downward SW flux (W m-2)
                                'forcing_rainf' ]   ! Rainfall rate (kg m-2 s-1)

contains

  !> @brief  Create the JULES forcing fields and their time axis
  !> @details The fields are 2D W3 fields in jules_forcing_fields, set from
  !>          the data in the jules_forcing file without time interpolation.
  !>          As for the 'nf' interpolation option of standalone JULES, each
  !>          data time is the start of the period the data applies to, and
  !>          the data is used unchanged until the next data time.
  !> @param[in,out] depository            The depository field collection
  !> @param[in,out] jules_forcing_fields  Collection to hold the forcing fields
  !> @param[in]     mesh                  The current 3d mesh
  !> @param[in]     twod_mesh             The current 2d mesh
  !> @param[in,out] model_axes            Time axes, to hold the forcing axis
  subroutine create_jules_forcing_fields( depository, jules_forcing_fields, &
                                          mesh, twod_mesh, model_axes )

    implicit none

    type(field_collection_type), intent(inout)       :: depository
    type(field_collection_type), intent(inout)       :: jules_forcing_fields
    type(mesh_type),             intent(in), pointer :: mesh
    type(mesh_type),             intent(in), pointer :: twod_mesh
    type(gungho_time_axes_type), intent(inout)       :: model_axes

    integer(i_def) :: i

    call log_event( 'Creating JULES site forcing fields', LOG_LEVEL_INFO )

    if (associated(model_axes%jules_forcing_time_axis)) then
      call log_event( 'Attempt to recreate JULES forcing time axis', &
                      LOG_LEVEL_ERROR )
    end if
    allocate(model_axes%jules_forcing_time_axis)
    call model_axes%jules_forcing_time_axis%initialise(        &
            "jules_forcing_time", file_id="jules_forcing",     &
            interp_flag=.false., window_size=1 )

    do i = 1, size(jules_forcing_field_names)
      call setup_ancil_field( trim(jules_forcing_field_names(i)), depository, &
                              jules_forcing_fields, mesh, twod_mesh,          &
                              twod=.true.,                                    &
                              time_axis=model_axes%jules_forcing_time_axis )
    end do

    call model_axes%jules_forcing_times_list%insert_item( &
                                   model_axes%jules_forcing_time_axis )

  end subroutine create_jules_forcing_fields

end module init_jules_forcing_mod
