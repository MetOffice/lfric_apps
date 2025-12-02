!-----------------------------------------------------------------------------
! (C) Crown copyright 2025 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
program test_swift_outer_update_tracer

  use constants_mod,                  only : r_tran, r_second
  use ffsl_advective_updates_alg_mod, only : swift_outer_update_tracer
  use fs_continuity_mod,              only : W3
  use function_space_mod,             only : function_space_type
  use function_space_collection_mod,  only : function_space_collection
  use test_helpers_mod,               only : create_test_field_flood_constant, &
                                             print_test_field
  use test_tiny_world_mod,            only : initialise_tiny_world, &
                                             finalise_tiny_world,   &
                                             mesh

  implicit none

  character(:), allocatable :: option

  call initialise_tiny_world()
  call run_test_64()
  call finalise_tiny_world()

contains

  subroutine run_test_64()

    use, intrinsic :: iso_fortran_env, only : real64

    use field_real64_mod, only : field_real64_type
    use test_field_mod,   only : test_field_64_type

    implicit none

    type(function_space_type), pointer :: w3_fs

    type(field_real64_type),  allocatable, target :: post_x, post_y
    type(test_field_64_type), allocatable         :: post_x_test, post_y_test
    type(field_real64_type),  allocatable, target :: tracer_field
    type(test_field_64_type), allocatable         :: tracer_test_field
    type(field_real64_type),  allocatable, target :: dry_mass_next
    type(test_field_64_type), allocatable         :: dry_mass_next_test
    type(field_real64_type),  allocatable, target :: dry_mass_x, dry_mass_y
    type(test_field_64_type), allocatable         :: dry_mass_x_test, &
                                                     dry_mass_y_test
    type(field_real64_type),  allocatable, target :: inc_x, inc_y
    type(test_field_64_type), allocatable         :: inc_x_test, inc_y_test

    w3_fs => function_space_collection%get_fs(mesh, 0, 0, W3)

    call create_test_field_flood_constant( w3_fs, "tracer",                 &
                                           tracer_field, tracer_test_field, &
                                           1.0_real64 )

    call create_test_field_flood_constant( w3_fs, "post x", &
                                           post_x, post_x_test, 0.2_real64 )
    call create_test_field_flood_constant( w3_fs, "post y", &
                                           post_y, post_y_test, 0.4_real64 )

    call create_test_field_flood_constant( w3_fs, "dry mass next",            &
                                           dry_mass_next, dry_mass_next_test, &
                                           7.5_real64 )

    call create_test_field_flood_constant( w3_fs, "dry mass x",         &
                                           dry_mass_x, dry_mass_x_test, &
                                           0.3_real64 )
    call create_test_field_flood_constant( w3_fs, "dry mass y",         &
                                           dry_mass_y, dry_mass_y_test, &
                                           0.5_real64 )

    call create_test_field_flood_constant( w3_fs, "increment x", &
                                           inc_x, inc_x_test, 0.7_real64 )
    call create_test_field_flood_constant( w3_fs, "increment y", &
                                           inc_y, inc_y_test, 0.9_real64 )

    call swift_outer_update_tracer( tracer_field,           &
                                    post_x, post_y,         &
                                    dry_mass_next,          &
                                    dry_mass_x, dry_mass_y, &
                                    inc_x, inc_y,           &
                                    dt=0.5_r_second )

    call print_test_field( tracer_test_field )

  end subroutine run_test_64

end program test_swift_outer_update_tracer
