!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!------------------------------------------------------------------------------
!> @brief Computes the mask to define the normal winds on the boundary.
!> @details The mask is defined using the W2 dofs. The coordinates of the dofs
!!          are compared with the coordinates of the boundary, to within a
!!          prescribed tolerance. To define the mask, the dofs on the boundary
!!          take the value 1, whilst those in the interior/exterior are zero.
!!          The normal winds are used to apply Dirichlet boundary conditions in
!!          the linear solver.
!!          If rim_width_ns < 0 then use periodic boundary conditions in y.
!!          If rim_width_ns >=0 then use W2 Dirichlet boundary conditions in y.
!!          If rim_width_ew < 0 then use periodic boundary conditions in x.
!!          If rim_width_ew >=0 then use W2 Dirichlet boundary conditions in x.
module create_w2bmask_kernel_mod

  use argument_mod,              only : arg_type, func_type, &
                                        GH_SCALAR, GH_FIELD, &
                                        GH_REAL, GH_INTEGER, &
                                        GH_READ, GH_WRITE,   &
                                        GH_BASIS,            &
                                        CELL_COLUMN, GH_EVALUATOR
  use constants_mod,             only : r_def, i_def
  use fs_continuity_mod,         only : W2, Wchi
  use kernel_mod,                only : kernel_type
  use base_mesh_config_mod,      only : geometry,              &
                                        geometry_spherical
  use finite_element_config_mod, only : coord_system, &
                                        coord_system_xyz

  implicit none

  private

  !-------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------

  type, public, extends(kernel_type) :: create_w2bmask_kernel_type
    private
    type(arg_type) :: meta_args(8) = (/                  &
         arg_type(GH_SCALAR,  GH_REAL,    GH_READ),      &
         arg_type(GH_SCALAR,  GH_REAL,    GH_READ),      &
         arg_type(GH_SCALAR,  GH_REAL,    GH_READ),      &
         arg_type(GH_SCALAR,  GH_REAL,    GH_READ),      &
         arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),      &
         arg_type(GH_SCALAR,  GH_INTEGER, GH_READ),      &
         arg_type(GH_FIELD,   GH_REAL,    GH_WRITE, W2), &
         arg_type(GH_FIELD*3, GH_REAL,    GH_READ, Wchi) &
         /)
    type(func_type) :: meta_funcs(1) = (/                &
         func_type(Wchi, GH_BASIS)                       &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: create_w2bmask_code
  end type

  !-------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------
  public :: create_w2bmask_code

contains

!> @param[in] nlayers      Number of layers
!> @param[in] boundary_s   Southern boundary coordinates
!> @param[in] boundary_n   Northern boundary coordinates
!> @param[in] boundary_e   Eastern  boundary coordinates
!> @param[in] boundary_w   Western  boundary coordinates
!> @param[in] rim_width_ns If >= 0 apply LBC in north south
!> @param[in] rim_width_ew If >= 0 apply LBC in east west
!> @param[in,out] u_mask   The mask to define the boundary for W2 (wind, u)
!> @param[in] chi_1        X component of the chi coordinate field
!> @param[in] chi_2        Y component of the chi coordinate field
!> @param[in] chi_3        Z component of the chi coordinate field
!> @param[in] ndf_w2       Number of degrees of freedom per cell for W2
!> @param[in] undf_w2      Total number of degrees of freedom for W2
!> @param[in] map_w2       Dofmap for the cell at the base of the column for W2
!> @param[in] ndf_chi      Number of degrees of freedom per cell for chi
!> @param[in] undf_chi     Number of degrees of freedom for chi
!> @param[in] map_chi      Dofmap for the cell at the base of the column for chi
!> @param[in] chi_basis    Basis functions evaluated at gaussian quadrature points
subroutine create_w2bmask_code( nlayers,      &
                                boundary_s,   &
                                boundary_n,   &
                                boundary_e,   &
                                boundary_w,   &
                                rim_width_ns, &
                                rim_width_ew, &
                                u_mask,       &
                                chi_1,        &
                                chi_2,        &
                                chi_3,        &
                                ndf_w2,       &
                                undf_w2,      &
                                map_w2,       &
                                ndf_chi,      &
                                undf_chi,     &
                                map_chi,      &
                                chi_basis )

  implicit none

  ! Arguments
  integer(kind=i_def),                           intent(in) :: nlayers,       &
                                                               ndf_w2,        &
                                                               ndf_chi,       &
                                                               undf_w2,       &
                                                               undf_chi,      &
                                                               rim_width_ns,  &
                                                               rim_width_ew
  integer(kind=i_def), dimension(ndf_w2),        intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_chi),       intent(in) :: map_chi
  real(kind=r_def),                              intent(in) :: boundary_s,    &
                                                               boundary_n,    &
                                                               boundary_e,    &
                                                               boundary_w
  real(kind=r_def), dimension(undf_w2),          intent(inout) :: u_mask
  real(kind=r_def), dimension(1,ndf_chi,ndf_w2), intent(in)    :: chi_basis
  real(kind=r_def), dimension(undf_chi),         intent(in)    :: chi_1,      &
                                                                  chi_2,      &
                                                                  chi_3

  ! Internal variables
  integer(kind=i_def)                    :: k, df1, df2
  real(kind=r_def), dimension(ndf_chi)   :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def)                       :: x(3),           &
                                            query_value_ns, &
                                            query_value_ew
  real(kind=r_def), parameter            :: tol = 1.0e-15_r_def

  ! Calculate the coordinates for the bottom vertical layer and
  ! copy the resulting mask to other layers
  k = 0

  ! Get the chi values in the kth layer
  ! chi_e (u) = chi ( f(u) + k )
  do df2 = 1, ndf_chi
    chi_1_e(df2) = chi_1( map_chi(df2) + k)
    chi_2_e(df2) = chi_2( map_chi(df2) + k)
    chi_3_e(df2) = chi_3( map_chi(df2) + k)
  end do

  ! Get the x values (x,y,z) for each cell w2 dof, by summing over
  ! W0 to W2 basis functions
  ! x(v) = sum_u chi_e(u) * g(u,v)
  do df1 = 1, ndf_w2
    x(:) = 0.0_r_def
    do df2 = 1,ndf_chi
      x(1) = x(1) + chi_1_e(df2) * chi_basis(1,df2,df1)
      x(2) = x(2) + chi_2_e(df2) * chi_basis(1,df2,df1)
      x(3) = x(3) + chi_3_e(df2) * chi_basis(1,df2,df1)
    end do

    if ( geometry == geometry_spherical .and. &
         coord_system == coord_system_xyz ) then
      ! In alpha beta space - and only pick the face where x(2)>0
      if ( x(2) > 0.0_r_def ) then
        query_value_ns = atan2(x(3),x(2))
        query_value_ew = atan2(x(1),x(2))
      else
      ! Set to an arbitrary value greater than pi/4
        query_value_ns = 10.0_r_def
        query_value_ew = 10.0_r_def
      endif
    else
      ! In cartesian space
      query_value_ns = x(2)
      query_value_ew = x(1)
    end if

    ! Define the limited area domain: 1 on the boundary, 0 elsewhere
    ! if rim_width_ns < 0 then use periodic boundary conditions in y
    ! if rim_width_ns >=0 then use W2 Dirichlet boundary conditions in y
    ! if rim_width_ew < 0 then use periodic boundary conditions in x
    ! if rim_width_ew >=0 then use W2 Dirichlet boundary conditions in x

    ! 111111
    ! 100001
    ! 100001
    ! 111111

    ! Loop over vertical layers
    do k = 0, nlayers-1
      u_mask(map_w2(df1) + k) = 0.0_r_def
    end do

    ! Check if the dof is on the North or South boundaries of the subdomain
    if ( rim_width_ns >= 0_i_def ) then
      ! Want a box, not a noughts and crosses pattern
      if ( query_value_ew > (boundary_w - tol) .and. &
           query_value_ew < (boundary_e + tol) ) then
        if ( (abs(query_value_ns - boundary_s) < tol) .or. &
             (abs(query_value_ns - boundary_n) < tol) ) then
          ! Loop over vertical layers
          do k = 0, nlayers-1
            u_mask(map_w2(df1) + k) = 1.0_r_def
          enddo

        endif
      endif
    endif

    ! Check if the dof is on the East or West boundaries of the subdomain
    if ( rim_width_ew >= 0_i_def ) then
      ! Want a box, not a noughts and crosses pattern
      if ( query_value_ns > (boundary_s - tol) .and. &
           query_value_ns < (boundary_n + tol) ) then
        if ( (abs(query_value_ew - boundary_e) < tol) .or. &
             (abs(query_value_ew - boundary_w) < tol) ) then
          ! Loop over vertical layers
          do k = 0, nlayers-1
            u_mask(map_w2(df1) + k) = 1.0_r_def
          enddo
        endif
      endif
    endif

  enddo

end subroutine create_w2bmask_code

end module create_w2bmask_kernel_mod
