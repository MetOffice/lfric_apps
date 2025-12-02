!-----------------------------------------------------------------------------
! (C) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> Module controlling the initialisation and finalisation of the:
!>  1. external communication group (XIOS)
!>  2. internal communication group used by the model instance.
!>
!
! Note:
!  This module is based on the components/driver/source/driver_comm_mod.F90
!  where the init_comm and final_comm routines have been split into external
!  setup (for just XIOS here) and internal setup for the model application.
!  This module is based on r42829 and further details on this update can be
!  found in #3803.
module jedi_lfric_comm_mod

  use constants_mod,         only: i_def
  use halo_comms_mod,        only: initialise_halo_comms, &
                                   finalise_halo_comms
  use lfric_mpi_mod,         only: global_mpi, &
                                   lfric_comm_type

! USE_XIOS flag used for models using the XIOS I/O server
#ifdef USE_XIOS
  use lfric_xios_driver_mod, only: lfric_xios_initialise, lfric_xios_finalise
#endif

  implicit none

  public :: init_external_comm, init_internal_comm, final_external_comm, final_internal_comm
  private

contains

  !> @brief  Initialises the external communicator groups
  !>
  !> @param[in]    program_name     The model name
  !> @param[in]    world_comm       The world communicator  which is used
  !>                                if XIOS is not included.
  !> @param[out]   output_comm      The communicator group available to the
  !>                                application after splitting by external
  !>                                libs (only XIOS at present).
  subroutine init_external_comm( program_name, world_comm, output_comm )

    implicit none

    character(len=*),  intent(in)    :: program_name
    integer(i_def),    intent(in)    :: world_comm
    integer(i_def),    intent(out)   :: output_comm

    ! Local
    logical :: comm_is_split
#ifdef USE_XIOS
    type(lfric_comm_type) :: lfric_comm
#endif

    ! Comm has not been split yet
    comm_is_split = .false.

#ifdef USE_XIOS
    ! Initialise XIOS and get back the split communicator
    call lfric_xios_initialise( program_name, lfric_comm, comm_is_split )
    ! Convert the LFRic communicator back to an mpi communicator
    output_comm = lfric_comm%get_comm_mpi_val()
    comm_is_split = .true.
#endif

    ! If communicator has not been split, set model comm as world comm
    if (.not. comm_is_split) output_comm = world_comm

  end subroutine init_external_comm

  !> @brief  Initialises the communicator used by the model
  !>
  !> @param[in]    model_communicator  The communicator to be used by the model
  subroutine init_internal_comm( model_communicator )

    implicit none

    integer(i_def), intent(in) :: model_communicator
    type(lfric_comm_type) :: lfric_comm

    ! Convert the mpi communicator to an LFRic communicator
    call lfric_comm%set_comm_mpi_val(model_communicator)

    ! Store the MPI communicator for later use
    call global_mpi%initialise( lfric_comm )

    ! Initialise halo functionality
    call initialise_halo_comms( lfric_comm )

  end subroutine init_internal_comm


  !> @brief  Finalises the external communicator groups
  subroutine final_external_comm()

    implicit none

#ifdef USE_XIOS
    ! Finalise XIOS
    call lfric_xios_finalise()
#endif

  end subroutine final_external_comm


  !> @brief  Finalise the lfric model communication group
  subroutine final_internal_comm()

    implicit none

    ! Finalise halo exchange functionality
    call finalise_halo_comms()

    ! Finalise the mpi object
    call global_mpi%finalise()

  end subroutine final_internal_comm

end module jedi_lfric_comm_mod
