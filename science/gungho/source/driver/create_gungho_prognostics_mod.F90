!-----------------------------------------------------------------------------
! (C) Crown copyright 2019 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Create empty fields for use by the gungho model
!> @details Creates the empty prognostic fields that will be later
!>          initialise and used by the gungho model. Field creation
!>          consists of three stages: constructing the field object,
!>          placing it in the depository (so it doesn't go out of scope) and
!>          putting a pointer to the depository version of the field into
!>          a 'prognostic_fields' field collection
module create_gungho_prognostics_mod

  use constants_mod,                  only : i_def, l_def, r_def
  use field_mod,                      only : field_type
  use field_parent_mod,               only : write_interface, read_interface, &
                                             checkpoint_write_interface, &
                                             checkpoint_read_interface
  use field_collection_mod,           only : field_collection_type
  use field_spec_mod,                 only : main_coll_dict, &
                                             adv_coll_dict, &
                                             moist_arr_dict, &
                                             processor_type, &
                                             make_spec
  use field_mapper_mod,               only : field_mapper_type
  use field_maker_mod,                only : field_maker_type
  use finite_element_config_mod,      only : element_order_h, &
                                             element_order_v
  use fs_continuity_mod,              only : W0, W2, W3, Wtheta, W2H, W2V
  use function_space_collection_mod , only : function_space_collection
  use log_mod,                        only : log_event,         &
                                             LOG_LEVEL_INFO, LOG_LEVEL_WARNING
  use mesh_mod,                       only : mesh_type
  use mixed_solver_config_mod,        only : reference_reset_time
  use mr_indices_mod,                 only : nummr, &
                                             mr_names
  use moist_dyn_mod,                  only : num_moist_factors, &
                                             moist_dyn_names
  use initialization_config_mod,      only: init_option,               &
                                            init_option_checkpoint_dump
  use io_config_mod,                  only: checkpoint_read, checkpoint_write
  use transport_config_mod,           only : transport_ageofair
  use clock_mod,                      only : clock_type
  implicit none

  private
  public :: create_gungho_prognostics, process_gungho_prognostics

contains

  !> @brief Iterate over active gungho fields and apply an arbitrary
  !> processor to the field specifiers.
  !> @details To be used by create_gungho_prognostics to create fields and by
  !> gungho_model_mod / initialise_infrastrucure to enable checkpoint fields.
  !> These two operations have to be separate because they have to happen
  !> at different times. The enabling of checkpoint fields has to happen
  !> before the io context closes, and this is too early for field creation.
  !> @param  proc Processor to be applied to selected field specifiers
  subroutine process_gungho_prognostics(proc)
    use field_spec_mod,            only : main => main_coll_dict, &
                                          adv => adv_coll_dict
    implicit none

    class(processor_type) :: proc
    class(clock_type), pointer :: clock
    integer(i_def) :: imr, reference_reset_freq, ord_h, ord_v
    logical(l_def) :: legacy
    logical(l_def) :: checkpoint_flag
    logical(l_def) :: is_empty
    real(r_def)    :: dt

    clock => proc%get_clock()

    ! Get the horizontal and vertical order of the function spaces
    ord_h = element_order_h
    ord_v = element_order_v

    ! enable/disable legacy checkpointing
    legacy = .true.

    call proc%apply(make_spec('theta', main%none, Wtheta, order_h=ord_h, &
                              order_v=ord_v, ckp=.true., legacy=legacy))
    call proc%apply(make_spec('u', main%none, W2, order_h=ord_h, order_v=ord_v,&
                              ckp=.true., legacy=legacy))
    if (.not. legacy) then
      call proc%apply(make_spec('h_u', main%none, W2H, order_h=ord_h, &
                                order_v=ord_v, ckp=.true.))
      call proc%apply(make_spec('v_u', main%none, W2V, order_h=ord_h, &
                                order_v=ord_v, ckp=.true.))
    end if
    call proc%apply(make_spec('rho', main%none, W3, order_h=ord_h, &
                              order_v=ord_v, ckp=.true., legacy=legacy))
    call proc%apply(make_spec('exner', main%none, W3, order_h=ord_h, &
                              order_v=ord_v, ckp=.true., legacy=legacy))

    ! Create reference fields for solver. They are only created and checkpointed if either the first or
    ! final steps are not semi-implicit operator recalculation timesteps
    dt = real(clock%get_seconds_per_step(), r_def)
    reference_reset_freq = nint(reference_reset_time / dt, i_def)
    checkpoint_flag =                                                &
      mod(clock%get_first_step()-1, reference_reset_freq) /= 0 .or.  &
      mod(clock%get_last_step(),    reference_reset_freq) /= 0
    if (checkpoint_read .or. &
         init_option == init_option_checkpoint_dump) then
      ! If the first timestep of this run IS an operator calc timestep, but the
      ! first timestep of the next run IS NOT, then checkpoint_flag
      ! must be false to allow model to start running, as the operator
      ! prognostics will not be in the initial dump
      if (mod(clock%get_first_step()-1, reference_reset_freq) == 0 .and.       &
          mod(clock%get_last_step(),    reference_reset_freq) /= 0) then
        checkpoint_flag = .false.
        if (checkpoint_write) then
          ! Any dump written will be incomplete, and the following run
          ! will need to start on an operator calc timestep - print user a
          ! warning
          call log_event('Danger: start of this run is an operator calc ' //   &
                         'timestep, but start of next run is not. ' //         &
                         'Written dump will be incomplete. Next run ' //       &
                         'must start with an operator calc timestep',          &
                         LOG_LEVEL_WARNING)
        end if
      end if
    end if
    is_empty = .not. checkpoint_flag
    call proc%apply(make_spec('theta_ref', main%none, Wtheta, empty=is_empty, &
                    order_h=ord_h, order_v=ord_v, ckp=checkpoint_flag))
    call proc%apply(make_spec('rho_ref', main%none, W3, empty=is_empty, &
                    order_h=ord_h, order_v=ord_v, ckp=checkpoint_flag))
    call proc%apply(make_spec('exner_ref', main%none, W3,  empty=is_empty, &
                    order_h=ord_h, order_v=ord_v, ckp=checkpoint_flag))

    ! The moisture mixing ratio fields (mr) and moist dynamics fields
    ! (moist_dyn) are always passed into the timestep algorithm, so are
    ! always created here, even when moisture_formulation = 'dry'
    do imr = 1,nummr
      call proc%apply(make_spec(trim(mr_names(imr)), main%none, &
        Wtheta, moist_arr=moist_arr_dict%mr, moist_idx=imr, order_h=ord_h,     &
        order_v=ord_v, ckp=.true., legacy=legacy))
    end do

    ! Auxiliary fields holding moisture-dependent factors for dynamics, including checkpointed versions for
    ! semi-implicit operator recalculations
    do imr = 1, num_moist_factors
      call proc%apply(make_spec(trim(moist_dyn_names(imr)), main%none, &
        Wtheta, moist_arr=moist_arr_dict%moist_dyn, moist_idx=imr, &
        order_h=ord_h, order_v=ord_v))
      call proc%apply(make_spec(trim('moist_dyn_'//trim(moist_dyn_names(imr))//'_ref'), &
        main%none, Wtheta,  empty=is_empty, moist_arr=moist_arr_dict%moist_dyn_ref, &
        moist_idx=imr, order_h=ord_h, order_v=ord_v, ckp=checkpoint_flag))
    end do

    if (transport_ageofair) then
      call proc%apply(make_spec('ageofair', main%none, &
        W3, adv_coll=adv%last_con, order_h=ord_h, order_v=ord_v, ckp=.true., &
        legacy=legacy))
    end if
  end subroutine process_gungho_prognostics

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Create empty fields to be used as prognostics by the gungho model
  !> @param[in]    mesh       The current 3d mesh
  !> @param[in]    twod_mesh  The current 2d mesh
  !> @param[in]    coarse_mesh          The coarse 3d mesh
  !> @param[in]    coarse_twod_mesh     The coarse 2d mesh
  !> @param[in]    mapper     Provides access to the field collections
  !> @param[in]    clock      The model clock
  subroutine create_gungho_prognostics(mesh,                                  &
                                       twod_mesh,                             &
                                       coarse_mesh,                           &
                                       coarse_twod_mesh,                      &
                                       mapper,                                &
                                       clock)

    implicit none

    type(mesh_type), intent(in), pointer      :: mesh
    type(mesh_type), intent(in), pointer      :: twod_mesh
    type(mesh_type), intent(in), pointer      :: coarse_mesh
    type(mesh_type), intent(in), pointer      :: coarse_twod_mesh
    type( field_mapper_type ), intent(in)     :: mapper
    class( clock_type ), intent(in)           :: clock

    type( field_maker_type ) :: creator

    call log_event( 'GungHo: Creating prognostics...', LOG_LEVEL_INFO )

    call creator%init(mesh, twod_mesh, coarse_mesh, coarse_twod_mesh, mapper, clock)

    call process_gungho_prognostics(creator)

  end subroutine create_gungho_prognostics

end module create_gungho_prognostics_mod
