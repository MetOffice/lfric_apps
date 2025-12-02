!-----------------------------------------------------------------------------
! (c) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Create LBC fields.
!> @details Create LBC field collection and add fields.
module create_lbcs_mod

  use constants_mod,              only : i_def, l_def, str_def
  use log_mod,                    only : log_event,             &
                                         log_scratch_space,     &
                                         LOG_LEVEL_INFO,        &
                                         LOG_LEVEL_ERROR
  use field_mod,                  only : field_type
  use field_collection_mod,       only : field_collection_type
  use clock_mod,                  only : clock_type
  use field_mapper_mod,           only : field_mapper_type
  use field_maker_mod,            only : field_maker_type
  use fs_continuity_mod,          only : W0, W2, W3, Wtheta, W2H, W2V
  use function_space_mod,         only : function_space_type
  use mesh_mod,                   only : mesh_type
  use mr_indices_mod,             only : nummr,                      &
                                         mr_names
  use gungho_time_axes_mod,       only : gungho_time_axes_type
  use initialization_config_mod,  only : lbc_option,             &
                                         lbc_option_analytic,    &
                                         lbc_option_gungho_file, &
                                         lbc_option_um2lfric_file
  use aerosol_config_mod,         only : murk_lbc

  implicit none

  public  :: process_lbc_fields, create_lbc_fields

  contains

   !> @brief Iterate over active lbc fields and apply an arbitrary
   !> processor to the field specifiers.
   !> @details To be used by create_lbc_fields to create fields and by
   !> gungho_model_mod / initialise_infrastrucure to enable checkpoint fields.
   !> These two operations have to be separate because they have to happen
   !> at different times. The enabling of checkpoint fields has to happen
   !> before the io context closes, and this is too early for field creation.
   !> @param  proc Processor to be applied to selected field specifiers
   subroutine process_lbc_fields(proc)
      use field_spec_mod,             only : main => main_coll_dict,            &
                                             axis => time_axis_dict,            &
                                             processor_type,                    &
                                             make_spec
      implicit none

      class(processor_type) :: proc

      integer(i_def) :: imr
      character(str_def) :: name
      logical(l_def) :: legacy

      legacy = .true.

      select case( lbc_option )

        case ( lbc_option_analytic )

          call proc%apply(make_spec('lbc_theta', main%lbc, Wtheta, ckp=.true., legacy=legacy))
          call proc%apply(make_spec('lbc_u', main%lbc, W2, ckp=.true., legacy=legacy))
          call proc%apply(make_spec('lbc_h_u', main%lbc, W2H, ckp=.true., legacy=legacy))
          call proc%apply(make_spec('lbc_v_u', main%lbc, W2V, ckp=.true., legacy=legacy))
          call proc%apply(make_spec('lbc_rho', main%lbc, W3, ckp=.true., legacy=legacy))
          call proc%apply(make_spec('lbc_exner', main%lbc, W3, ckp=.true., legacy=legacy))
          call proc%apply(make_spec('boundary_u_diff', main%lbc, W2, ckp=.true., legacy=legacy))
          if (.not. legacy) then
            call proc%apply(make_spec('boundary_h_u_diff', main%lbc, W2H, ckp=.true., legacy=.false.))
            call proc%apply(make_spec('boundary_v_u_diff', main%lbc, W2V, ckp=.true., legacy=.false.))
          end if
          call proc%apply(make_spec('boundary_u_driving', main%lbc, W2, ckp=.true., legacy=legacy))
          if (.not. legacy) then
            call proc%apply(make_spec('boundary_h_u_driving', main%lbc, W2H, ckp=.true., legacy=.false.))
            call proc%apply(make_spec('boundary_v_u_driving', main%lbc, W2V, ckp=.true., legacy=.false.))
          end if
          do imr = 1, nummr
            name = trim('lbc_') // adjustl(mr_names(imr))
            call proc%apply(make_spec(name, main%lbc, Wtheta, ckp=.true., legacy=legacy))
          enddo

        case ( lbc_option_gungho_file )

          !------ Fields updated directly from LBC file-----------------
          call proc%apply(make_spec('lbc_theta', main%lbc, Wtheta, time_axis=axis%lbc))
          call proc%apply(make_spec('lbc_rho', main%lbc, W3, time_axis=axis%lbc))
          call proc%apply(make_spec('lbc_exner', main%lbc, W3, time_axis=axis%lbc))
          call proc%apply(make_spec('lbc_h_u', main%lbc, W2H, time_axis=axis%lbc))
          call proc%apply(make_spec('lbc_v_u', main%lbc, Wtheta, time_axis=axis%lbc))

          !----- Fields derived from the fields in the LBC file---------
          call proc%apply(make_spec('lbc_u', main%lbc, W2))
          call proc%apply(make_spec('boundary_u_diff', main%lbc, W2))
          call proc%apply(make_spec('boundary_u_driving', main%lbc, W2))

        case ( lbc_option_um2lfric_file )
          !------ Fields updated directly from LBC file-----------------
          call proc%apply(make_spec('lbc_theta', main%lbc, Wtheta, time_axis=axis%lbc))
          call proc%apply(make_spec('lbc_rho_r2', main%lbc, W3, time_axis=axis%lbc))
          call proc%apply(make_spec('lbc_h_u', main%lbc, W2H, time_axis=axis%lbc))
          call proc%apply(make_spec('lbc_v_u', main%lbc, Wtheta, time_axis=axis%lbc))

          ! Specific humidities (Initially these are read in, but
          ! in the long run we should move to just use the mixing
          ! ratios)
          call proc%apply(make_spec('lbc_q', main%lbc, Wtheta, time_axis=axis%lbc))
          call proc%apply(make_spec('lbc_qcl', main%lbc, Wtheta, time_axis=axis%lbc))
          call proc%apply(make_spec('lbc_qcf', main%lbc, Wtheta, time_axis=axis%lbc))
          call proc%apply(make_spec('lbc_qrain', main%lbc, Wtheta, time_axis=axis%lbc))
          if (murk_lbc) then
            call proc%apply(make_spec('lbc_murk', main%lbc, Wtheta, time_axis=axis%lbc))
          end if

          !----- Fields derived from the fields in the LBC file---------
          call proc%apply(make_spec('lbc_rho', main%lbc, W3))
          call proc%apply(make_spec('lbc_exner', main%lbc, W3))
          call proc%apply(make_spec('lbc_u', main%lbc, W2))
          call proc%apply(make_spec('boundary_u_diff', main%lbc, W2))
          call proc%apply(make_spec('boundary_u_driving', main%lbc, W2))

          ! Mixing ratios
          call proc%apply(make_spec('lbc_m_v', main%lbc, Wtheta))
          call proc%apply(make_spec('lbc_m_cl', main%lbc, Wtheta))
          call proc%apply(make_spec('lbc_m_s', main%lbc, Wtheta))
          call proc%apply(make_spec('lbc_m_r', main%lbc, Wtheta))

        case default
          call log_event( 'This lbc_option not available', LOG_LEVEL_ERROR )
      end select

    end subroutine process_lbc_fields

    !> @brief   Create and add LBC fields.
    !> @details Create the lateral boundary condition field collection.
    !>          On every timestep these fields will be updated and used by the
    !>          limited area model.
    !> @param[in]    mesh       The current 3d mesh
    !> @param[in]    twod_mesh  The current 2d mesh (not used here)
    !> @param[in]    coarse_mesh      The coarse 3d mesh
    !> @param[in]    coarse_twod_mesh The coarse 2d mesh
    !> @param[in]    mapper     Provides access to the field collections
    !> @param[in]    clock      The model clock
    subroutine create_lbc_fields(mesh, twod_mesh, coarse_mesh, coarse_twod_mesh, mapper, clock)
     implicit none
     type(mesh_type), intent(in), pointer    :: mesh
     type(mesh_type), intent(in), pointer    :: twod_mesh
     type(mesh_type), intent(in), pointer    :: coarse_mesh
     type(mesh_type), intent(in), pointer    :: coarse_twod_mesh
     type(field_mapper_type), intent(in)     :: mapper
     class(clock_type), intent(in)           :: clock

     type(gungho_time_axes_type), pointer    :: gungho_axes

     type(field_maker_type) :: creator

     call log_event('GungHo: Creating lbc fields...', LOG_LEVEL_INFO)

     gungho_axes => mapper%get_gungho_axes()
     if (lbc_option /= lbc_option_analytic) call gungho_axes%make_lbc_time_axis()

     call creator%init(mesh, twod_mesh, coarse_mesh, coarse_twod_mesh, mapper, clock)

     call process_lbc_fields(creator)

     call gungho_axes%save_lbc_time_axis()

    end subroutine create_lbc_fields

end module create_lbcs_mod
