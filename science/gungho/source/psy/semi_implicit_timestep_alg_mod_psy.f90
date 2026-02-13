module semi_implicit_timestep_alg_mod_psy
  use constants_mod
  use field_mod, only : field_proxy_type, field_proxy_type_10=>field_proxy_type, field_proxy_type_11=>field_proxy_type, &
&field_proxy_type_12=>field_proxy_type, field_proxy_type_1=>field_proxy_type, field_proxy_type_2=>field_proxy_type, &
&field_proxy_type_3=>field_proxy_type, field_proxy_type_4=>field_proxy_type, field_proxy_type_5=>field_proxy_type, &
&field_proxy_type_6=>field_proxy_type, field_proxy_type_7=>field_proxy_type, field_proxy_type_8=>field_proxy_type, &
&field_proxy_type_9=>field_proxy_type, field_type
  use r_tran_field_mod, only : r_tran_field_proxy_type, r_tran_field_type
  use operator_mod, only : operator_proxy_type, operator_type
  implicit none
  public

  contains
  subroutine invoke_0(self_dtheta_cld)
    use mesh_mod, only : mesh_type
    type(field_type), intent(in) :: self_dtheta_cld
    integer(kind=i_def) :: df
    type(mesh_type), pointer :: mesh => null()
    integer(kind=i_def) :: max_halo_depth_mesh
    real(kind=r_def), pointer, dimension(:) :: self_dtheta_cld_data => null()
    type(field_proxy_type) :: self_dtheta_cld_proxy
    integer(kind=i_def) :: loop0_start
    integer(kind=i_def) :: loop0_stop

    ! Initialise field and/or operator proxies
    self_dtheta_cld_proxy = self_dtheta_cld%get_proxy()
    self_dtheta_cld_data => self_dtheta_cld_proxy%data

    ! Create a mesh object
    mesh => self_dtheta_cld_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()

    ! Set-up all of the loop bounds
    loop0_start = 1
    loop0_stop = self_dtheta_cld_proxy%vspace%get_last_dof_halo(1)

    ! Call kernels and communication routines
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop0_start, loop0_stop, 1
      ! Built-in: setval_c (set a real-valued field to a real scalar value)
      self_dtheta_cld_data(df) = 0.0_r_def
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_dtheta_cld_proxy%set_dirty()
    call self_dtheta_cld_proxy%set_clean(1)

  end subroutine invoke_0
  subroutine invoke_copy_init_fields_to_state(self_state, u, self_state_1, theta, self_state_2, rho, self_state_3, exner, &
&self_total_dry_flux)
    use mesh_mod, only : mesh_type
    type(field_type), intent(in) :: self_state
    type(field_type), intent(in) :: u
    type(field_type), intent(in) :: self_state_1
    type(field_type), intent(in) :: theta
    type(field_type), intent(in) :: self_state_2
    type(field_type), intent(in) :: rho
    type(field_type), intent(in) :: self_state_3
    type(field_type), intent(in) :: exner
    type(r_tran_field_type), intent(in) :: self_total_dry_flux
    integer(kind=i_def) :: df
    type(mesh_type), pointer :: mesh => null()
    integer(kind=i_def) :: max_halo_depth_mesh
    real(kind=r_def), pointer, dimension(:) :: self_state_data => null()
    real(kind=r_def), pointer, dimension(:) :: u_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_state_1_data => null()
    real(kind=r_def), pointer, dimension(:) :: theta_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_state_2_data => null()
    real(kind=r_def), pointer, dimension(:) :: rho_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_state_3_data => null()
    real(kind=r_def), pointer, dimension(:) :: exner_data => null()
    real(kind=r_tran), pointer, dimension(:) :: self_total_dry_flux_data => null()
    type(field_proxy_type_1) :: self_state_proxy
    type(field_proxy_type_1) :: u_proxy
    type(field_proxy_type_1) :: self_state_1_proxy
    type(field_proxy_type_1) :: theta_proxy
    type(field_proxy_type_1) :: self_state_2_proxy
    type(field_proxy_type_1) :: rho_proxy
    type(field_proxy_type_1) :: self_state_3_proxy
    type(field_proxy_type_1) :: exner_proxy
    type(r_tran_field_proxy_type) :: self_total_dry_flux_proxy
    integer(kind=i_def) :: loop0_start
    integer(kind=i_def) :: loop0_stop
    integer(kind=i_def) :: loop1_start
    integer(kind=i_def) :: loop1_stop
    integer(kind=i_def) :: loop2_start
    integer(kind=i_def) :: loop2_stop
    integer(kind=i_def) :: loop3_start
    integer(kind=i_def) :: loop3_stop
    integer(kind=i_def) :: loop4_start
    integer(kind=i_def) :: loop4_stop

    ! Initialise field and/or operator proxies
    self_state_proxy = self_state%get_proxy()
    self_state_data => self_state_proxy%data
    u_proxy = u%get_proxy()
    u_data => u_proxy%data
    self_state_1_proxy = self_state_1%get_proxy()
    self_state_1_data => self_state_1_proxy%data
    theta_proxy = theta%get_proxy()
    theta_data => theta_proxy%data
    self_state_2_proxy = self_state_2%get_proxy()
    self_state_2_data => self_state_2_proxy%data
    rho_proxy = rho%get_proxy()
    rho_data => rho_proxy%data
    self_state_3_proxy = self_state_3%get_proxy()
    self_state_3_data => self_state_3_proxy%data
    exner_proxy = exner%get_proxy()
    exner_data => exner_proxy%data
    self_total_dry_flux_proxy = self_total_dry_flux%get_proxy()
    self_total_dry_flux_data => self_total_dry_flux_proxy%data

    ! Create a mesh object
    mesh => self_state_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()

    ! Set-up all of the loop bounds
    loop0_start = 1
    loop0_stop = self_state_proxy%vspace%get_last_dof_annexed()
    loop1_start = 1
    loop1_stop = self_state_1_proxy%vspace%get_last_dof_annexed()
    loop2_start = 1
    loop2_stop = self_state_2_proxy%vspace%get_last_dof_annexed()
    loop3_start = 1
    loop3_stop = self_state_3_proxy%vspace%get_last_dof_annexed()
    loop4_start = 1
    loop4_stop = self_total_dry_flux_proxy%vspace%get_last_dof_halo(1)

    ! Call kernels and communication routines
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop0_start, loop0_stop, 1
      ! Built-in: setval_X (set a real-valued field equal to another such field)
      self_state_data(df) = u_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_state_proxy%set_dirty()
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop1_start, loop1_stop, 1
      ! Built-in: setval_X (set a real-valued field equal to another such field)
      self_state_1_data(df) = theta_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_state_1_proxy%set_dirty()
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop2_start, loop2_stop, 1
      ! Built-in: setval_X (set a real-valued field equal to another such field)
      self_state_2_data(df) = rho_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_state_2_proxy%set_dirty()
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop3_start, loop3_stop, 1
      ! Built-in: setval_X (set a real-valued field equal to another such field)
      self_state_3_data(df) = exner_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_state_3_proxy%set_dirty()
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop4_start, loop4_stop, 1
      ! Built-in: setval_c (set a real-valued field to a real scalar value)
      self_total_dry_flux_data(df) = 0.0_r_tran
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_total_dry_flux_proxy%set_dirty()
    call self_total_dry_flux_proxy%set_clean(1)

  end subroutine invoke_copy_init_fields_to_state
  subroutine invoke_2(theta_ref, self_state, rho_ref, self_state_1, exner_ref, self_state_2)
    use mesh_mod, only : mesh_type
    type(field_type), intent(in) :: theta_ref
    type(field_type), intent(in) :: self_state
    type(field_type), intent(in) :: rho_ref
    type(field_type), intent(in) :: self_state_1
    type(field_type), intent(in) :: exner_ref
    type(field_type), intent(in) :: self_state_2
    integer(kind=i_def) :: df
    type(mesh_type), pointer :: mesh => null()
    integer(kind=i_def) :: max_halo_depth_mesh
    real(kind=r_def), pointer, dimension(:) :: theta_ref_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_state_data => null()
    real(kind=r_def), pointer, dimension(:) :: rho_ref_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_state_1_data => null()
    real(kind=r_def), pointer, dimension(:) :: exner_ref_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_state_2_data => null()
    type(field_proxy_type_2) :: theta_ref_proxy
    type(field_proxy_type_2) :: self_state_proxy
    type(field_proxy_type_2) :: rho_ref_proxy
    type(field_proxy_type_2) :: self_state_1_proxy
    type(field_proxy_type_2) :: exner_ref_proxy
    type(field_proxy_type_2) :: self_state_2_proxy
    integer(kind=i_def) :: loop0_start
    integer(kind=i_def) :: loop0_stop
    integer(kind=i_def) :: loop1_start
    integer(kind=i_def) :: loop1_stop
    integer(kind=i_def) :: loop2_start
    integer(kind=i_def) :: loop2_stop

    ! Initialise field and/or operator proxies
    theta_ref_proxy = theta_ref%get_proxy()
    theta_ref_data => theta_ref_proxy%data
    self_state_proxy = self_state%get_proxy()
    self_state_data => self_state_proxy%data
    rho_ref_proxy = rho_ref%get_proxy()
    rho_ref_data => rho_ref_proxy%data
    self_state_1_proxy = self_state_1%get_proxy()
    self_state_1_data => self_state_1_proxy%data
    exner_ref_proxy = exner_ref%get_proxy()
    exner_ref_data => exner_ref_proxy%data
    self_state_2_proxy = self_state_2%get_proxy()
    self_state_2_data => self_state_2_proxy%data

    ! Create a mesh object
    mesh => theta_ref_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()

    ! Set-up all of the loop bounds
    loop0_start = 1
    loop0_stop = theta_ref_proxy%vspace%get_last_dof_annexed()
    loop1_start = 1
    loop1_stop = rho_ref_proxy%vspace%get_last_dof_annexed()
    loop2_start = 1
    loop2_stop = exner_ref_proxy%vspace%get_last_dof_annexed()

    ! Call kernels and communication routines
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop0_start, loop0_stop, 1
      ! Built-in: setval_X (set a real-valued field equal to another such field)
      theta_ref_data(df) = self_state_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call theta_ref_proxy%set_dirty()
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop1_start, loop1_stop, 1
      ! Built-in: setval_X (set a real-valued field equal to another such field)
      rho_ref_data(df) = self_state_1_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call rho_ref_proxy%set_dirty()
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop2_start, loop2_stop, 1
      ! Built-in: setval_X (set a real-valued field equal to another such field)
      exner_ref_data(df) = self_state_2_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call exner_ref_proxy%set_dirty()

  end subroutine invoke_2
  subroutine invoke_update_from_slow_physics(self_state_after_slow, self_dtheta, self_state_after_slow_1, self_du)
    use mesh_mod, only : mesh_type
    type(field_type), intent(in) :: self_state_after_slow
    type(field_type), intent(in) :: self_dtheta
    type(field_type), intent(in) :: self_state_after_slow_1
    type(field_type), intent(in) :: self_du
    integer(kind=i_def) :: df
    type(mesh_type), pointer :: mesh => null()
    integer(kind=i_def) :: max_halo_depth_mesh
    real(kind=r_def), pointer, dimension(:) :: self_state_after_slow_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_dtheta_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_state_after_slow_1_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_du_data => null()
    type(field_proxy_type_3) :: self_state_after_slow_proxy
    type(field_proxy_type_3) :: self_dtheta_proxy
    type(field_proxy_type_3) :: self_state_after_slow_1_proxy
    type(field_proxy_type_3) :: self_du_proxy
    integer(kind=i_def) :: loop0_start
    integer(kind=i_def) :: loop0_stop
    integer(kind=i_def) :: loop1_start
    integer(kind=i_def) :: loop1_stop

    ! Initialise field and/or operator proxies
    self_state_after_slow_proxy = self_state_after_slow%get_proxy()
    self_state_after_slow_data => self_state_after_slow_proxy%data
    self_dtheta_proxy = self_dtheta%get_proxy()
    self_dtheta_data => self_dtheta_proxy%data
    self_state_after_slow_1_proxy = self_state_after_slow_1%get_proxy()
    self_state_after_slow_1_data => self_state_after_slow_1_proxy%data
    self_du_proxy = self_du%get_proxy()
    self_du_data => self_du_proxy%data

    ! Create a mesh object
    mesh => self_state_after_slow_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()

    ! Set-up all of the loop bounds
    loop0_start = 1
    loop0_stop = self_state_after_slow_proxy%vspace%get_last_dof_annexed()
    loop1_start = 1
    loop1_stop = self_state_after_slow_1_proxy%vspace%get_last_dof_annexed()

    ! Call kernels and communication routines
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop0_start, loop0_stop, 1
      ! Built-in: inc_X_plus_Y (increment a real-valued field)
      self_state_after_slow_data(df) = self_state_after_slow_data(df) + self_dtheta_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_state_after_slow_proxy%set_dirty()
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop1_start, loop1_stop, 1
      ! Built-in: inc_X_plus_Y (increment a real-valued field)
      self_state_after_slow_1_data(df) = self_state_after_slow_1_data(df) + self_du_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_state_after_slow_1_proxy%set_dirty()

  end subroutine invoke_update_from_slow_physics
  subroutine invoke_4(self_wind_prev, self_state)
    use mesh_mod, only : mesh_type
    type(field_type), intent(in) :: self_wind_prev
    type(field_type), intent(in) :: self_state
    integer(kind=i_def) :: df
    type(mesh_type), pointer :: mesh => null()
    integer(kind=i_def) :: max_halo_depth_mesh
    real(kind=r_def), pointer, dimension(:) :: self_wind_prev_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_state_data => null()
    type(field_proxy_type_4) :: self_wind_prev_proxy
    type(field_proxy_type_4) :: self_state_proxy
    integer(kind=i_def) :: loop0_start
    integer(kind=i_def) :: loop0_stop

    ! Initialise field and/or operator proxies
    self_wind_prev_proxy = self_wind_prev%get_proxy()
    self_wind_prev_data => self_wind_prev_proxy%data
    self_state_proxy = self_state%get_proxy()
    self_state_data => self_state_proxy%data

    ! Create a mesh object
    mesh => self_wind_prev_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()

    ! Set-up all of the loop bounds
    loop0_start = 1
    loop0_stop = self_wind_prev_proxy%vspace%get_last_dof_annexed()

    ! Call kernels and communication routines
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop0_start, loop0_stop, 1
      ! Built-in: setval_X (set a real-valued field equal to another such field)
      self_wind_prev_data(df) = self_state_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_wind_prev_proxy%set_dirty()

  end subroutine invoke_4
  subroutine invoke_5(self_advected_state, self_rhs_adv, self_advected_state_1, self_rhs_adv_1, self_advected_state_2, self_du)
    use mesh_mod, only : mesh_type
    type(field_type), intent(in) :: self_advected_state
    type(field_type), intent(in) :: self_rhs_adv
    type(field_type), intent(in) :: self_advected_state_1
    type(field_type), intent(in) :: self_rhs_adv_1
    type(field_type), intent(in) :: self_advected_state_2
    type(field_type), intent(in) :: self_du
    integer(kind=i_def) :: df
    type(mesh_type), pointer :: mesh => null()
    integer(kind=i_def) :: max_halo_depth_mesh
    real(kind=r_def), pointer, dimension(:) :: self_advected_state_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_rhs_adv_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_advected_state_1_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_rhs_adv_1_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_advected_state_2_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_du_data => null()
    type(field_proxy_type_5) :: self_advected_state_proxy
    type(field_proxy_type_5) :: self_rhs_adv_proxy
    type(field_proxy_type_5) :: self_advected_state_1_proxy
    type(field_proxy_type_5) :: self_rhs_adv_1_proxy
    type(field_proxy_type_5) :: self_advected_state_2_proxy
    type(field_proxy_type_5) :: self_du_proxy
    integer(kind=i_def) :: loop0_start
    integer(kind=i_def) :: loop0_stop
    integer(kind=i_def) :: loop1_start
    integer(kind=i_def) :: loop1_stop
    integer(kind=i_def) :: loop2_start
    integer(kind=i_def) :: loop2_stop

    ! Initialise field and/or operator proxies
    self_advected_state_proxy = self_advected_state%get_proxy()
    self_advected_state_data => self_advected_state_proxy%data
    self_rhs_adv_proxy = self_rhs_adv%get_proxy()
    self_rhs_adv_data => self_rhs_adv_proxy%data
    self_advected_state_1_proxy = self_advected_state_1%get_proxy()
    self_advected_state_1_data => self_advected_state_1_proxy%data
    self_rhs_adv_1_proxy = self_rhs_adv_1%get_proxy()
    self_rhs_adv_1_data => self_rhs_adv_1_proxy%data
    self_advected_state_2_proxy = self_advected_state_2%get_proxy()
    self_advected_state_2_data => self_advected_state_2_proxy%data
    self_du_proxy = self_du%get_proxy()
    self_du_data => self_du_proxy%data

    ! Create a mesh object
    mesh => self_advected_state_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()

    ! Set-up all of the loop bounds
    loop0_start = 1
    loop0_stop = self_advected_state_proxy%vspace%get_last_dof_annexed()
    loop1_start = 1
    loop1_stop = self_advected_state_1_proxy%vspace%get_last_dof_annexed()
    loop2_start = 1
    loop2_stop = self_advected_state_2_proxy%vspace%get_last_dof_annexed()

    ! Call kernels and communication routines
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop0_start, loop0_stop, 1
      ! Built-in: inc_X_plus_Y (increment a real-valued field)
      self_advected_state_data(df) = self_advected_state_data(df) + self_rhs_adv_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_advected_state_proxy%set_dirty()
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop1_start, loop1_stop, 1
      ! Built-in: inc_X_plus_Y (increment a real-valued field)
      self_advected_state_1_data(df) = self_advected_state_1_data(df) + self_rhs_adv_1_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_advected_state_1_proxy%set_dirty()
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop2_start, loop2_stop, 1
      ! Built-in: inc_X_plus_Y (increment a real-valued field)
      self_advected_state_2_data(df) = self_advected_state_2_data(df) + self_du_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_advected_state_2_proxy%set_dirty()

  end subroutine invoke_5
  subroutine invoke_6(self_wind_prev, self_state)
    use mesh_mod, only : mesh_type
    type(field_type), intent(in) :: self_wind_prev
    type(field_type), intent(in) :: self_state
    integer(kind=i_def) :: df
    type(mesh_type), pointer :: mesh => null()
    integer(kind=i_def) :: max_halo_depth_mesh
    real(kind=r_def), pointer, dimension(:) :: self_wind_prev_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_state_data => null()
    type(field_proxy_type_6) :: self_wind_prev_proxy
    type(field_proxy_type_6) :: self_state_proxy
    integer(kind=i_def) :: loop0_start
    integer(kind=i_def) :: loop0_stop

    ! Initialise field and/or operator proxies
    self_wind_prev_proxy = self_wind_prev%get_proxy()
    self_wind_prev_data => self_wind_prev_proxy%data
    self_state_proxy = self_state%get_proxy()
    self_state_data => self_state_proxy%data

    ! Create a mesh object
    mesh => self_wind_prev_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()

    ! Set-up all of the loop bounds
    loop0_start = 1
    loop0_stop = self_wind_prev_proxy%vspace%get_last_dof_annexed()

    ! Call kernels and communication routines
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop0_start, loop0_stop, 1
      ! Built-in: setval_X (set a real-valued field equal to another such field)
      self_wind_prev_data(df) = self_state_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_wind_prev_proxy%set_dirty()

  end subroutine invoke_6
  subroutine invoke_7(self_dtheta, self_rhs_adv, mm_wt)
    use mesh_mod, only : mesh_type
    use dg_inc_matrix_vector_kernel_mod, only : dg_inc_matrix_vector_code
    type(field_type), intent(in) :: self_dtheta
    type(field_type), intent(in) :: self_rhs_adv
    type(operator_type), intent(in) :: mm_wt
    integer(kind=i_def) :: df
    integer(kind=i_def) :: cell
    type(mesh_type), pointer :: mesh => null()
    integer(kind=i_def) :: max_halo_depth_mesh
    real(kind=r_def), pointer, dimension(:) :: self_dtheta_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_rhs_adv_data => null()
    real(kind=r_def), pointer, dimension(:,:,:) :: mm_wt_local_stencil => null()
    integer(kind=i_def) :: nlayers_self_rhs_adv
    integer(kind=i_def) :: ndf_aspc1_self_dtheta
    integer(kind=i_def) :: undf_aspc1_self_dtheta
    integer(kind=i_def) :: ndf_aspc1_self_rhs_adv
    integer(kind=i_def) :: undf_aspc1_self_rhs_adv
    integer(kind=i_def) :: ndf_adspc1_self_rhs_adv
    integer(kind=i_def) :: undf_adspc1_self_rhs_adv
    integer(kind=i_def), pointer :: map_adspc1_self_rhs_adv(:,:) => null()
    integer(kind=i_def), pointer :: map_aspc1_self_dtheta(:,:) => null()
    type(field_proxy_type_7) :: self_dtheta_proxy
    type(field_proxy_type_7) :: self_rhs_adv_proxy
    type(operator_proxy_type) :: mm_wt_proxy
    integer(kind=i_def) :: loop0_start
    integer(kind=i_def) :: loop0_stop
    integer(kind=i_def) :: loop1_start
    integer(kind=i_def) :: loop1_stop
    integer(kind=i_def) :: loop2_start
    integer(kind=i_def) :: loop2_stop

    ! Initialise field and/or operator proxies
    self_dtheta_proxy = self_dtheta%get_proxy()
    self_dtheta_data => self_dtheta_proxy%data
    self_rhs_adv_proxy = self_rhs_adv%get_proxy()
    self_rhs_adv_data => self_rhs_adv_proxy%data
    mm_wt_proxy = mm_wt%get_proxy()
    mm_wt_local_stencil => mm_wt_proxy%local_stencil

    ! Initialise number of layers
    nlayers_self_rhs_adv = self_rhs_adv_proxy%vspace%get_nlayers()

    ! Create a mesh object
    mesh => self_dtheta_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()

    ! Look-up dofmaps for each function space
    map_adspc1_self_rhs_adv => self_rhs_adv_proxy%vspace%get_whole_dofmap()
    map_aspc1_self_dtheta => self_dtheta_proxy%vspace%get_whole_dofmap()

    ! Initialise number of DoFs for aspc1_self_dtheta
    ndf_aspc1_self_dtheta = self_dtheta_proxy%vspace%get_ndf()
    undf_aspc1_self_dtheta = self_dtheta_proxy%vspace%get_undf()

    ! Initialise number of DoFs for aspc1_self_rhs_adv
    ndf_aspc1_self_rhs_adv = self_rhs_adv_proxy%vspace%get_ndf()
    undf_aspc1_self_rhs_adv = self_rhs_adv_proxy%vspace%get_undf()

    ! Initialise number of DoFs for adspc1_self_rhs_adv
    ndf_adspc1_self_rhs_adv = self_rhs_adv_proxy%vspace%get_ndf()
    undf_adspc1_self_rhs_adv = self_rhs_adv_proxy%vspace%get_undf()

    ! Set-up all of the loop bounds
    loop0_start = 1
    loop0_stop = self_dtheta_proxy%vspace%get_last_dof_annexed()
    loop1_start = 1
    loop1_stop = self_rhs_adv_proxy%vspace%get_last_dof_halo(1)
    loop2_start = 1
    loop2_stop = mesh%get_last_edge_cell()

    ! Call kernels and communication routines
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop0_start, loop0_stop, 1
      ! Built-in: setval_X (set a real-valued field equal to another such field)
      self_dtheta_data(df) = self_rhs_adv_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_dtheta_proxy%set_dirty()
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop1_start, loop1_stop, 1
      ! Built-in: setval_c (set a real-valued field to a real scalar value)
      self_rhs_adv_data(df) = 0.0_r_def
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_rhs_adv_proxy%set_dirty()
    call self_rhs_adv_proxy%set_clean(1)
    !$omp parallel default(shared) private(cell)
    !$omp do schedule(static)
    do cell = loop2_start, loop2_stop, 1
      call dg_inc_matrix_vector_code(cell, nlayers_self_rhs_adv, self_rhs_adv_data, self_dtheta_data, mm_wt_proxy%ncell_3d, &
&mm_wt_local_stencil, ndf_adspc1_self_rhs_adv, undf_adspc1_self_rhs_adv, map_adspc1_self_rhs_adv(:,cell), ndf_aspc1_self_dtheta, &
&undf_aspc1_self_dtheta, map_aspc1_self_dtheta(:,cell))
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_rhs_adv_proxy%set_dirty()

  end subroutine invoke_7
  subroutine invoke_update_rhs_phys_from_fast_physics(self_rhs_phys, self_dtheta, mm_wt, self_rhs_phys_1, self_du, mm_vel)
    use mesh_mod, only : mesh_type
    use dg_inc_matrix_vector_kernel_mod, only : dg_inc_matrix_vector_code
    use matrix_vector_kernel_mod, only : matrix_vector_code
    use sci_enforce_bc_kernel_mod, only : enforce_bc_code
    type(field_type), intent(in) :: self_rhs_phys
    type(field_type), intent(in) :: self_dtheta
    type(operator_type), intent(in) :: mm_wt
    type(field_type), intent(in) :: self_rhs_phys_1
    type(field_type), intent(in) :: self_du
    type(operator_type), intent(in) :: mm_vel
    integer(kind=i_def) :: cell
    type(mesh_type), pointer :: mesh => null()
    integer(kind=i_def) :: max_halo_depth_mesh
    real(kind=r_def), pointer, dimension(:) :: self_rhs_phys_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_dtheta_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_rhs_phys_1_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_du_data => null()
    real(kind=r_def), pointer, dimension(:,:,:) :: mm_wt_local_stencil => null()
    real(kind=r_def), pointer, dimension(:,:,:) :: mm_vel_local_stencil => null()
    integer(kind=i_def) :: nlayers_self_rhs_phys
    integer(kind=i_def) :: nlayers_self_rhs_phys_1
    integer(kind=i_def) :: colour
    integer(kind=i_def), pointer :: cmap(:,:)
    integer(kind=i_def) :: ncolour
    integer(kind=i_def), allocatable, dimension(:,:) :: last_halo_cell_all_colours
    integer(kind=i_def) :: ndf_adspc1_self_rhs_phys
    integer(kind=i_def) :: undf_adspc1_self_rhs_phys
    integer(kind=i_def) :: ndf_aspc1_self_dtheta
    integer(kind=i_def) :: undf_aspc1_self_dtheta
    integer(kind=i_def) :: ndf_aspc1_self_rhs_phys_1
    integer(kind=i_def) :: undf_aspc1_self_rhs_phys_1
    integer(kind=i_def) :: ndf_aspc2_self_du
    integer(kind=i_def) :: undf_aspc2_self_du
    integer(kind=i_def), pointer :: map_adspc1_self_rhs_phys(:,:) => null()
    integer(kind=i_def), pointer :: map_aspc1_self_dtheta(:,:) => null()
    integer(kind=i_def), pointer :: map_aspc1_self_rhs_phys_1(:,:) => null()
    integer(kind=i_def), pointer :: map_aspc2_self_du(:,:) => null()
    integer(kind=i_def), pointer :: boundary_dofs_self_rhs_phys_1(:,:) => null()
    type(field_proxy_type_8) :: self_rhs_phys_proxy
    type(field_proxy_type_8) :: self_dtheta_proxy
    type(field_proxy_type_8) :: self_rhs_phys_1_proxy
    type(field_proxy_type_8) :: self_du_proxy
    type(operator_proxy_type) :: mm_wt_proxy
    type(operator_proxy_type) :: mm_vel_proxy
    integer(kind=i_def) :: loop0_start
    integer(kind=i_def) :: loop0_stop
    integer(kind=i_def) :: loop1_start
    integer(kind=i_def) :: loop1_stop
    integer(kind=i_def) :: loop2_start
    integer(kind=i_def) :: loop3_start
    integer(kind=i_def) :: loop3_stop
    integer(kind=i_def) :: loop4_start

    ! Initialise field and/or operator proxies
    self_rhs_phys_proxy = self_rhs_phys%get_proxy()
    self_rhs_phys_data => self_rhs_phys_proxy%data
    self_dtheta_proxy = self_dtheta%get_proxy()
    self_dtheta_data => self_dtheta_proxy%data
    mm_wt_proxy = mm_wt%get_proxy()
    mm_wt_local_stencil => mm_wt_proxy%local_stencil
    self_rhs_phys_1_proxy = self_rhs_phys_1%get_proxy()
    self_rhs_phys_1_data => self_rhs_phys_1_proxy%data
    self_du_proxy = self_du%get_proxy()
    self_du_data => self_du_proxy%data
    mm_vel_proxy = mm_vel%get_proxy()
    mm_vel_local_stencil => mm_vel_proxy%local_stencil

    ! Initialise number of layers
    nlayers_self_rhs_phys = self_rhs_phys_proxy%vspace%get_nlayers()
    nlayers_self_rhs_phys_1 = self_rhs_phys_1_proxy%vspace%get_nlayers()

    ! Create a mesh object
    mesh => self_rhs_phys_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()

    ! Get the colourmap
    ncolour = mesh%get_ncolours()
    cmap => mesh%get_colour_map()

    ! Look-up dofmaps for each function space
    map_adspc1_self_rhs_phys => self_rhs_phys_proxy%vspace%get_whole_dofmap()
    map_aspc1_self_dtheta => self_dtheta_proxy%vspace%get_whole_dofmap()
    map_aspc1_self_rhs_phys_1 => self_rhs_phys_1_proxy%vspace%get_whole_dofmap()
    map_aspc2_self_du => self_du_proxy%vspace%get_whole_dofmap()
    boundary_dofs_self_rhs_phys_1 => self_rhs_phys_1_proxy%vspace%get_boundary_dofs()

    ! Initialise number of DoFs for adspc1_self_rhs_phys
    ndf_adspc1_self_rhs_phys = self_rhs_phys_proxy%vspace%get_ndf()
    undf_adspc1_self_rhs_phys = self_rhs_phys_proxy%vspace%get_undf()

    ! Initialise number of DoFs for aspc1_self_dtheta
    ndf_aspc1_self_dtheta = self_dtheta_proxy%vspace%get_ndf()
    undf_aspc1_self_dtheta = self_dtheta_proxy%vspace%get_undf()

    ! Initialise number of DoFs for aspc1_self_rhs_phys_1
    ndf_aspc1_self_rhs_phys_1 = self_rhs_phys_1_proxy%vspace%get_ndf()
    undf_aspc1_self_rhs_phys_1 = self_rhs_phys_1_proxy%vspace%get_undf()

    ! Initialise number of DoFs for aspc2_self_du
    ndf_aspc2_self_du = self_du_proxy%vspace%get_ndf()
    undf_aspc2_self_du = self_du_proxy%vspace%get_undf()
    last_halo_cell_all_colours = mesh%get_last_halo_cell_all_colours()

    ! Set-up all of the loop bounds
    loop0_start = 1
    loop0_stop = mesh%get_last_edge_cell()
    loop1_start = 1
    loop1_stop = ncolour
    loop2_start = 1
    loop3_start = 1
    loop3_stop = ncolour
    loop4_start = 1

    ! Call kernels and communication routines
    !$omp parallel default(shared) private(cell)
    !$omp do schedule(static)
    do cell = loop0_start, loop0_stop, 1
      call dg_inc_matrix_vector_code(cell, nlayers_self_rhs_phys, self_rhs_phys_data, self_dtheta_data, mm_wt_proxy%ncell_3d, &
&mm_wt_local_stencil, ndf_adspc1_self_rhs_phys, undf_adspc1_self_rhs_phys, map_adspc1_self_rhs_phys(:,cell), &
&ndf_aspc1_self_dtheta, undf_aspc1_self_dtheta, map_aspc1_self_dtheta(:,cell))
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_rhs_phys_proxy%set_dirty()
    if (self_du_proxy%is_dirty(depth=1)) then
      call self_du_proxy%halo_exchange(depth=1)
    end if
    do colour = loop1_start, loop1_stop, 1
      !$omp parallel default(shared) private(cell)
      !$omp do schedule(static)
      do cell = loop2_start, last_halo_cell_all_colours(colour,1), 1
        call matrix_vector_code(cmap(colour,cell), nlayers_self_rhs_phys_1, self_rhs_phys_1_data, self_du_data, &
&mm_vel_proxy%ncell_3d, mm_vel_local_stencil, ndf_aspc1_self_rhs_phys_1, undf_aspc1_self_rhs_phys_1, &
&map_aspc1_self_rhs_phys_1(:,cmap(colour,cell)), ndf_aspc2_self_du, undf_aspc2_self_du, map_aspc2_self_du(:,cmap(colour,cell)))
      enddo
      !$omp end do
      !$omp end parallel
    enddo

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_rhs_phys_1_proxy%set_dirty()
    do colour = loop3_start, loop3_stop, 1
      !$omp parallel default(shared) private(cell)
      !$omp do schedule(static)
      do cell = loop4_start, last_halo_cell_all_colours(colour,1), 1
        call enforce_bc_code(nlayers_self_rhs_phys_1, self_rhs_phys_1_data, ndf_aspc1_self_rhs_phys_1, undf_aspc1_self_rhs_phys_1, &
&map_aspc1_self_rhs_phys_1(:,cmap(colour,cell)), boundary_dofs_self_rhs_phys_1)
      enddo
      !$omp end do
      !$omp end parallel
    enddo

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_rhs_phys_1_proxy%set_dirty()

  end subroutine invoke_update_rhs_phys_from_fast_physics
  subroutine invoke_9(self_rhs_np1, self_rhs_np1_1)
    use mesh_mod, only : mesh_type
    type(field_type), intent(in) :: self_rhs_np1
    type(field_type), intent(in) :: self_rhs_np1_1
    integer(kind=i_def) :: df
    type(mesh_type), pointer :: mesh => null()
    integer(kind=i_def) :: max_halo_depth_mesh
    real(kind=r_def), pointer, dimension(:) :: self_rhs_np1_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_rhs_np1_1_data => null()
    type(field_proxy_type_9) :: self_rhs_np1_proxy
    type(field_proxy_type_9) :: self_rhs_np1_1_proxy
    integer(kind=i_def) :: loop0_start
    integer(kind=i_def) :: loop0_stop
    integer(kind=i_def) :: loop1_start
    integer(kind=i_def) :: loop1_stop

    ! Initialise field and/or operator proxies
    self_rhs_np1_proxy = self_rhs_np1%get_proxy()
    self_rhs_np1_data => self_rhs_np1_proxy%data
    self_rhs_np1_1_proxy = self_rhs_np1_1%get_proxy()
    self_rhs_np1_1_data => self_rhs_np1_1_proxy%data

    ! Create a mesh object
    mesh => self_rhs_np1_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()

    ! Set-up all of the loop bounds
    loop0_start = 1
    loop0_stop = self_rhs_np1_proxy%vspace%get_last_dof_halo(1)
    loop1_start = 1
    loop1_stop = self_rhs_np1_1_proxy%vspace%get_last_dof_halo(1)

    ! Call kernels and communication routines
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop0_start, loop0_stop, 1
      ! Built-in: setval_c (set a real-valued field to a real scalar value)
      self_rhs_np1_data(df) = 0.0_r_def
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_rhs_np1_proxy%set_dirty()
    call self_rhs_np1_proxy%set_clean(1)
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop1_start, loop1_stop, 1
      ! Built-in: setval_c (set a real-valued field to a real scalar value)
      self_rhs_np1_1_data(df) = 0.0_r_def
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_rhs_np1_1_proxy%set_dirty()
    call self_rhs_np1_1_proxy%set_clean(1)

  end subroutine invoke_9
  subroutine invoke_10(self_state, self_state_after_slow, self_dtheta, self_state_1, self_state_after_slow_1, self_du)
    use mesh_mod, only : mesh_type
    type(field_type), intent(in) :: self_state
    type(field_type), intent(in) :: self_state_after_slow
    type(field_type), intent(in) :: self_dtheta
    type(field_type), intent(in) :: self_state_1
    type(field_type), intent(in) :: self_state_after_slow_1
    type(field_type), intent(in) :: self_du
    integer(kind=i_def) :: df
    type(mesh_type), pointer :: mesh => null()
    integer(kind=i_def) :: max_halo_depth_mesh
    real(kind=r_def), pointer, dimension(:) :: self_state_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_state_after_slow_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_dtheta_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_state_1_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_state_after_slow_1_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_du_data => null()
    type(field_proxy_type_10) :: self_state_proxy
    type(field_proxy_type_10) :: self_state_after_slow_proxy
    type(field_proxy_type_10) :: self_dtheta_proxy
    type(field_proxy_type_10) :: self_state_1_proxy
    type(field_proxy_type_10) :: self_state_after_slow_1_proxy
    type(field_proxy_type_10) :: self_du_proxy
    integer(kind=i_def) :: loop0_start
    integer(kind=i_def) :: loop0_stop
    integer(kind=i_def) :: loop1_start
    integer(kind=i_def) :: loop1_stop

    ! Initialise field and/or operator proxies
    self_state_proxy = self_state%get_proxy()
    self_state_data => self_state_proxy%data
    self_state_after_slow_proxy = self_state_after_slow%get_proxy()
    self_state_after_slow_data => self_state_after_slow_proxy%data
    self_dtheta_proxy = self_dtheta%get_proxy()
    self_dtheta_data => self_dtheta_proxy%data
    self_state_1_proxy = self_state_1%get_proxy()
    self_state_1_data => self_state_1_proxy%data
    self_state_after_slow_1_proxy = self_state_after_slow_1%get_proxy()
    self_state_after_slow_1_data => self_state_after_slow_1_proxy%data
    self_du_proxy = self_du%get_proxy()
    self_du_data => self_du_proxy%data

    ! Create a mesh object
    mesh => self_state_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()

    ! Set-up all of the loop bounds
    loop0_start = 1
    loop0_stop = self_state_proxy%vspace%get_last_dof_annexed()
    loop1_start = 1
    loop1_stop = self_state_1_proxy%vspace%get_last_dof_annexed()

    ! Call kernels and communication routines
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop0_start, loop0_stop, 1
      ! Built-in: X_plus_Y (add real-valued fields)
      self_state_data(df) = self_state_after_slow_data(df) + self_dtheta_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_state_proxy%set_dirty()
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop1_start, loop1_stop, 1
      ! Built-in: X_plus_Y (add real-valued fields)
      self_state_1_data(df) = self_state_after_slow_1_data(df) + self_du_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call self_state_1_proxy%set_dirty()

  end subroutine invoke_10
  subroutine invoke_11(u, self_state, theta, self_state_1, rho, self_state_2, exner, self_state_3)
    use mesh_mod, only : mesh_type
    type(field_type), intent(in) :: u
    type(field_type), intent(in) :: self_state
    type(field_type), intent(in) :: theta
    type(field_type), intent(in) :: self_state_1
    type(field_type), intent(in) :: rho
    type(field_type), intent(in) :: self_state_2
    type(field_type), intent(in) :: exner
    type(field_type), intent(in) :: self_state_3
    integer(kind=i_def) :: df
    type(mesh_type), pointer :: mesh => null()
    integer(kind=i_def) :: max_halo_depth_mesh
    real(kind=r_def), pointer, dimension(:) :: u_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_state_data => null()
    real(kind=r_def), pointer, dimension(:) :: theta_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_state_1_data => null()
    real(kind=r_def), pointer, dimension(:) :: rho_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_state_2_data => null()
    real(kind=r_def), pointer, dimension(:) :: exner_data => null()
    real(kind=r_def), pointer, dimension(:) :: self_state_3_data => null()
    type(field_proxy_type_11) :: u_proxy
    type(field_proxy_type_11) :: self_state_proxy
    type(field_proxy_type_11) :: theta_proxy
    type(field_proxy_type_11) :: self_state_1_proxy
    type(field_proxy_type_11) :: rho_proxy
    type(field_proxy_type_11) :: self_state_2_proxy
    type(field_proxy_type_11) :: exner_proxy
    type(field_proxy_type_11) :: self_state_3_proxy
    integer(kind=i_def) :: loop0_start
    integer(kind=i_def) :: loop0_stop
    integer(kind=i_def) :: loop1_start
    integer(kind=i_def) :: loop1_stop
    integer(kind=i_def) :: loop2_start
    integer(kind=i_def) :: loop2_stop
    integer(kind=i_def) :: loop3_start
    integer(kind=i_def) :: loop3_stop

    ! Initialise field and/or operator proxies
    u_proxy = u%get_proxy()
    u_data => u_proxy%data
    self_state_proxy = self_state%get_proxy()
    self_state_data => self_state_proxy%data
    theta_proxy = theta%get_proxy()
    theta_data => theta_proxy%data
    self_state_1_proxy = self_state_1%get_proxy()
    self_state_1_data => self_state_1_proxy%data
    rho_proxy = rho%get_proxy()
    rho_data => rho_proxy%data
    self_state_2_proxy = self_state_2%get_proxy()
    self_state_2_data => self_state_2_proxy%data
    exner_proxy = exner%get_proxy()
    exner_data => exner_proxy%data
    self_state_3_proxy = self_state_3%get_proxy()
    self_state_3_data => self_state_3_proxy%data

    ! Create a mesh object
    mesh => u_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()

    ! Set-up all of the loop bounds
    loop0_start = 1
    loop0_stop = u_proxy%vspace%get_last_dof_annexed()
    loop1_start = 1
    loop1_stop = theta_proxy%vspace%get_last_dof_annexed()
    loop2_start = 1
    loop2_stop = rho_proxy%vspace%get_last_dof_annexed()
    loop3_start = 1
    loop3_stop = exner_proxy%vspace%get_last_dof_annexed()

    ! Call kernels and communication routines
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop0_start, loop0_stop, 1
      ! Built-in: setval_X (set a real-valued field equal to another such field)
      u_data(df) = self_state_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call u_proxy%set_dirty()
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop1_start, loop1_stop, 1
      ! Built-in: setval_X (set a real-valued field equal to another such field)
      theta_data(df) = self_state_1_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call theta_proxy%set_dirty()
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop2_start, loop2_stop, 1
      ! Built-in: setval_X (set a real-valued field equal to another such field)
      rho_data(df) = self_state_2_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call rho_proxy%set_dirty()
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop3_start, loop3_stop, 1
      ! Built-in: setval_X (set a real-valued field equal to another such field)
      exner_data(df) = self_state_3_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call exner_proxy%set_dirty()

  end subroutine invoke_11
  subroutine invoke_12(copied_generic_field, single_generic_field)
    use mesh_mod, only : mesh_type
    type(field_type), intent(in) :: copied_generic_field
    type(field_type), intent(in) :: single_generic_field
    integer(kind=i_def) :: df
    type(mesh_type), pointer :: mesh => null()
    integer(kind=i_def) :: max_halo_depth_mesh
    real(kind=r_def), pointer, dimension(:) :: copied_generic_field_data => null()
    real(kind=r_def), pointer, dimension(:) :: single_generic_field_data => null()
    type(field_proxy_type_12) :: copied_generic_field_proxy
    type(field_proxy_type_12) :: single_generic_field_proxy
    integer(kind=i_def) :: loop0_start
    integer(kind=i_def) :: loop0_stop

    ! Initialise field and/or operator proxies
    copied_generic_field_proxy = copied_generic_field%get_proxy()
    copied_generic_field_data => copied_generic_field_proxy%data
    single_generic_field_proxy = single_generic_field%get_proxy()
    single_generic_field_data => single_generic_field_proxy%data

    ! Create a mesh object
    mesh => copied_generic_field_proxy%vspace%get_mesh()
    max_halo_depth_mesh = mesh%get_halo_depth()

    ! Set-up all of the loop bounds
    loop0_start = 1
    loop0_stop = copied_generic_field_proxy%vspace%get_last_dof_annexed()

    ! Call kernels and communication routines
    !$omp parallel default(shared) private(df)
    !$omp do schedule(static)
    do df = loop0_start, loop0_stop, 1
      ! Built-in: setval_X (set a real-valued field equal to another such field)
      copied_generic_field_data(df) = single_generic_field_data(df)
    enddo
    !$omp end do
    !$omp end parallel

    ! Set halos dirty/clean for fields modified in the above loop(s)
    call copied_generic_field_proxy%set_dirty()

  end subroutine invoke_12

end module semi_implicit_timestep_alg_mod_psy
