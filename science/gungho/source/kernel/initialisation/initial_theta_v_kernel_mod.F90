!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Computes the initial theta perturbation field
module initial_theta_v_kernel_mod

  use argument_mod,                  only: arg_type, func_type,        &
                                            GH_FIELD, GH_REAL,         &
                                            GH_WRITE, GH_READ,         &
                                            ANY_SPACE_9, GH_BASIS,     &
                                            ANY_DISCONTINUOUS_SPACE_3, &
                                            CELL_COLUMN, GH_EVALUATOR, &
                                            GH_INTEGER, GH_SCALAR
  use constants_mod,                 only: r_def, i_def, l_def, PI
  use fs_continuity_mod,             only: Wtheta
  use kernel_mod,                    only: kernel_type
  use idealised_config_mod,          only: test_squall_line, test_supercell

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: initial_theta_v_kernel_type
      private
      type(arg_type) :: meta_args(8) = (/                                      &
          arg_type(GH_FIELD,   GH_REAL, GH_WRITE, Wtheta),                     &
          arg_type(GH_FIELD,   GH_REAL, GH_WRITE, Wtheta),                     &
          arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta),                     &
          arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta),                     &
          arg_type(GH_FIELD,   GH_REAL, GH_READ,  Wtheta),                     &
          arg_type(GH_FIELD*3, GH_REAL, GH_READ,  ANY_SPACE_9),                &
          arg_type(GH_FIELD,   GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3),  &
          arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)                            &
      /)
      type(func_type) :: meta_funcs(1) = (/                                    &
            func_type(ANY_SPACE_9, GH_BASIS)                                   &
            /)
      integer :: operates_on = CELL_COLUMN
      integer :: gh_shape = GH_EVALUATOR
  contains
      procedure, nopass :: initial_theta_v_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: initial_theta_v_code

contains

  subroutine initial_theta_v_code(nlayers,                               &
                                  theta_v, exner_wt,                     &
                                  theta_eq, exner_eq,                    &
                                  height_wth,                            &
                                  chi_1, chi_2, chi_3,                   &
                                  panel_id, test,                        &
                                  ndf_wtheta, undf_wtheta, map_wtheta,   &
                                  ndf_chi, undf_chi, map_chi, chi_basis, &
                                  ndf_pid, undf_pid, map_pid             )

    use analytic_temperature_profiles_mod, only : integrate_theta_v, integrate_exner_surf
    use sci_chi_transform_mod,             only : chi2llr

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers, ndf_wtheta, ndf_chi, ndf_pid
    integer(kind=i_def), intent(in) :: undf_wtheta, undf_chi, undf_pid, test

    integer(kind=i_def), dimension(ndf_wtheta), intent(in) :: map_wtheta
    integer(kind=i_def), dimension(ndf_chi),    intent(in) :: map_chi
    integer(kind=i_def), dimension(ndf_pid),    intent(in) :: map_pid

    real(kind=r_def),    dimension(undf_wtheta),           intent(inout) :: theta_v
    real(kind=r_def),    dimension(undf_wtheta),           intent(inout) :: exner_wt
    real(kind=r_def),    dimension(undf_wtheta),           intent(in)    :: theta_eq
    real(kind=r_def),    dimension(undf_wtheta),           intent(in)    :: exner_eq
    real(kind=r_def),    dimension(undf_wtheta),           intent(in)    :: height_wth
    real(kind=r_def),    dimension(undf_chi),              intent(in)    :: chi_1, chi_2, chi_3
    real(kind=r_def),    dimension(undf_pid),              intent(in)    :: panel_id
    real(kind=r_def),    dimension(1,ndf_chi,ndf_wtheta),  intent(in)    :: chi_basis

    ! Internal variables
    real(kind=r_def),   dimension(ndf_chi) :: chi_1_e, chi_2_e, chi_3_e
    real(kind=r_def)                       :: coords(3)
    integer(kind=i_def)                    :: dfc, k, ipanel
    real(kind=r_def)                       :: lon, lat, radius, dlat
    real(kind=r_def)                       :: z_km1, z_k, z_kp1
    integer(kind=i_def)                    :: i, lat_index
    logical(kind=l_def)                    :: lat_index_found

    integer(kind=i_def), parameter :: num_iterations = 10
    integer(kind=i_def), parameter :: num_quad_points = 1200

    real(kind=r_def), dimension(nlayers+1, num_quad_points) :: theta_v_prev
    real(kind=r_def), dimension(nlayers+1, num_quad_points) :: theta_v_tab
    real(kind=r_def), dimension(nlayers+1, num_quad_points) :: exner_tab
    real(kind=r_def), dimension(num_quad_points)            :: lat_points
    real(kind=r_def), dimension(nlayers+1, num_quad_points) :: dtheta_v_dz
    real(kind=r_def), dimension(nlayers+1)                  :: ueq2
    real(kind=r_def), dimension(nlayers+1)                  :: dueq2_dz

    real(kind=r_def) :: z, zS, dzU, Us, Uc, A, B, C

    ipanel = int(panel_id(map_pid(1)), i_def)

    ! ======================================================================== !
    ! Set up coordinates
    ! ======================================================================== !

    ! Get latitude ---------------------------------------------------------
    do dfc = 1, ndf_chi
      chi_1_e(dfc) = chi_1( map_chi(dfc) )
      chi_2_e(dfc) = chi_2( map_chi(dfc) )
      chi_3_e(dfc) = chi_3( map_chi(dfc) )
    end do

    coords(:) = 0.0_r_def
    do dfc = 1, ndf_chi
      coords(1) = coords(1) + chi_1_e(dfc)*chi_basis(1,dfc,1)
      coords(2) = coords(2) + chi_2_e(dfc)*chi_basis(1,dfc,1)
      coords(3) = coords(3) + chi_3_e(dfc)*chi_basis(1,dfc,1)
    end do

    call chi2llr(coords(1), coords(2), coords(3), ipanel, lon, lat, radius)

    ! Make an array of latitude points for numerical integral --------------
    if (lat >= 0.0_r_def) then
      dlat = 0.5_r_def * PI / REAL(num_quad_points-1, r_def)
    else
      dlat = -0.5_r_def * PI / REAL(num_quad_points-1, r_def)
    end if
    lat_points(1) = 0.0_r_def
    lat_index_found = .false.
    do i = 2, num_quad_points
      lat_points(i) = lat_points(i-1) + dlat
      if (lat >= 0.0_r_def .and. lat_points(i) > lat .and. .not. lat_index_found) then
        lat_index = i
        lat_index_found = .true.
      else if (lat < 0.0_r_def .and. lat_points(i) < lat .and. .not. lat_index_found) then
        lat_index = i
        lat_index_found = .true.
      end if
    end do
    if (lat >= 0.0_r_def) then
      lat_points(num_quad_points) = PI/2.0_r_def
    else
      lat_points(num_quad_points) = -PI/2.0_r_def
    end if

    ! ======================================================================== !
    ! Calculate du^2/dz at the equator, for each level
    ! ======================================================================== !
    if ( test == test_squall_line ) then
      zS = 2500.0_r_def
      dzU = 1000.0_r_def
      Us = 12.0_r_def
      Uc = 5.0_r_def
    else if ( test == test_supercell ) then
      zS = 5000.0_r_def
      dzU = 1000.0_r_def
      Us = 30.0_r_def
      Uc = 15.0_r_def
    end if

    A = 0.25_r_def * (-zs**2 + 2.0_r_def*zS*dzU - dzU**2) / (zS * dzU)
    B = 0.5_r_def * (zS + dzU) / dzU
    C = - 0.25_r_def * (zS / dzU)

    do k = 1, nlayers+1
      z = height_wth(map_wtheta(1)+k-1)
      if (z < zS - dzU) then
        ueq2(k) = (Us*z/zS - Uc)**2
      else if (ABS(z - zS) < dzU) then
        ueq2(k) = ((A + B*z/zS + C*(z/zS)**2)*Us - Uc)**2
      else
        ueq2(k) = (Us - Uc)**2
      end if
    end do

    ! Fit quadratic polynomial through theta_v points, and take vertical
    ! derivative at each point
    ! Bottom point, assume linear gradient
    dueq2_dz(1) =                                                              &
      (ueq2(2) - ueq2(1)) / (height_wth(map_wtheta(1)+1) - height_wth(map_wtheta(1)))

    do k = 2, nlayers
      ! Calculate derivative at height_wth(map_wtheta(1)+k)
      z_km1 = height_wth(map_wtheta(1)+k-1)
      z_k = height_wth(map_wtheta(1)+k)
      z_kp1 = height_wth(map_wtheta(1)+k+1)
      dueq2_dz(k) = (                                                          &
        ueq2(k-1) * (z_k - z_kp1) / ((z_km1 - z_k) * (z_km1 - z_kp1))          &
        + ueq2(k) * (2.0_r_def*z_k - z_km1 - z_kp1)                            &
                    / ((z_k - z_km1) * (z_k - z_kp1))                          &
        + ueq2(k+1) * (z_k - z_km1) / ((z_kp1 - z_k) * (z_kp1 - z_km1))        &
      )
    end do
    ! Top point, assume linear gradient
    dueq2_dz(nlayers+1) =                                                      &
      (ueq2(nlayers+1) - ueq2(nlayers))                                        &
      / (height_wth(map_wtheta(1)+nlayers) - height_wth(map_wtheta(1)+nlayers-1))


    ! ======================================================================== !
    ! Calculate theta_v
    ! ======================================================================== !

    ! Initialise theta_v_prev
    do k = 1, nlayers+1
      theta_v_prev(k,:) = theta_eq(map_wtheta(1)+k-1)
    end do

    do i = 1, num_iterations
      ! Fit quadratic polynomial through theta_v points, and take vertical
      ! derivative at each point
      ! Bottom point, assume linear gradient
      dtheta_v_dz(1,:) =                                                       &
        (theta_v_prev(2,:) - theta_v_prev(1,:))                                &
        / (height_wth(map_wtheta(1)+1) - height_wth(map_wtheta(1)))

      do k = 2, nlayers
        ! Calculate derivative at height_wth(map_wtheta(1)+k)
        z_km1 = height_wth(map_wtheta(1)+k-1)
        z_k = height_wth(map_wtheta(1)+k)
        z_kp1 = height_wth(map_wtheta(1)+k+1)
        dtheta_v_dz(k,:) = (                                                        &
          theta_v_prev(k-1,:) * (z_k - z_kp1) / ((z_km1 - z_k) * (z_km1 - z_kp1))   &
          + theta_v_prev(k,:) * (2.0_r_def*z_k - z_km1 - z_kp1)                     &
                      / ((z_k - z_km1) * (z_k - z_kp1))                             &
          + theta_v_prev(k+1,:) * (z_k - z_km1) / ((z_kp1 - z_k) * (z_kp1 - z_km1)) &
        )
      end do
      ! Top point, assume linear gradient
      dtheta_v_dz(nlayers+1,:) =                                               &
        (theta_v_prev(nlayers+1,:) - theta_v_prev(nlayers,:))                  &
        / (height_wth(map_wtheta(1)+nlayers) - height_wth(map_wtheta(1)+nlayers-1))

      call integrate_theta_v(                                                  &
              theta_v_tab, theta_v_prev, dtheta_v_dz,                          &
              theta_eq(map_wtheta(1) : map_wtheta(1)+nlayers), dueq2_dz,       &
              height_wth(map_wtheta(1) : map_wtheta(1)+nlayers),               &
              lat_points, nlayers+1, num_quad_points, test                     &
      )

      theta_v_prev(:,:) = theta_v_tab(:,:)
    end do

    call integrate_exner_surf(                                                 &
            exner_tab, theta_v_tab,                                            &
            exner_eq(map_wtheta(1) : map_wtheta(1)+nlayers),                   &
            height_wth(map_wtheta(1) : map_wtheta(1)+nlayers),                 &
            lat_points, nlayers+1, num_quad_points, test                       &
    )

    ! Determine final value
    do k = 0, nlayers
      theta_v(map_wtheta(1)+k) = theta_v_tab(k+1,lat_index)                    &
          + (theta_v_tab(k+1,lat_index+1) - theta_v_tab(k+1,lat_index))        &
          * (lat - lat_points(lat_index)) / (lat_points(lat_index+1) - lat_points(lat_index))
      exner_wt(map_wtheta(1)+k) = exner_tab(k+1,lat_index)                     &
          + (exner_tab(k+1,lat_index+1) - exner_tab(k+1,lat_index))            &
          * (lat - lat_points(lat_index)) / (lat_points(lat_index+1) - lat_points(lat_index))
    end do


  end subroutine initial_theta_v_code

end module initial_theta_v_kernel_mod
