!-----------------------------------------------------------------------------
! (C) Crown copyright 2020 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the mask to define the interior of the limited area region
!!        for W3 space.
!> @details The coordinates of the W3 dofs are compared with the coordinates
!!          of the boundary. The mask takes the value 1 in the interior and
!!          0 on the exterior.
module create_w3mask_kernel_mod

  use argument_mod,              only : arg_type, func_type, &
                                        GH_SCALAR, GH_FIELD, &
                                        GH_READ, GH_WRITE,   &
                                        GH_REAL, GH_BASIS,   &
                                        CELL_COLUMN, GH_EVALUATOR
  use constants_mod,             only : r_def, i_def, l_def
  use fs_continuity_mod,         only : W3, Wchi
  use kernel_mod,                only : kernel_type
  use base_mesh_config_mod,      only : geometry,            &
                                        geometry_spherical
  use finite_element_config_mod, only : coord_system, &
                                        coord_system_xyz

  implicit none

  private

  !-------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------

  type, public, extends(kernel_type) :: create_w3mask_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),       &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),       &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),       &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),       &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, W3),  &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ,  Wchi) &
         /)
    type(func_type) :: meta_funcs(1) = (/              &
         func_type(Wchi, GH_BASIS)                     &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: create_w3mask_code
  end type

  !-------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------
  public :: create_w3mask_code

contains

!> @param[in] nlayers      Number of layers
!> @param[in] boundary_s   Southern boundary coordinates
!> @param[in] boundary_n   Northern boundary coordinates
!> @param[in] boundary_e   Eastern  boundary coordinates
!> @param[in] boundary_w   Western  boundary coordinates
!> @param[in,out] rho_mask The mask to define the interior for W3 (density, rho)
!> @param[in] chi_1        X component of the chi coordinate field
!> @param[in] chi_2        Y component of the chi coordinate field
!> @param[in] chi_3        Z component of the chi coordinate field
!> @param[in] ndf_w3       Number of degrees of freedom per cell for W3
!> @param[in] undf_w3      Total number of degrees of freedom for W3
!> @param[in] map_w3       Dofmap for the cell at the base of the column for W3
!> @param[in] ndf_chi      Number of degrees of freedom per cell for chi
!> @param[in] undf_chi     Number of degrees of freedom for chi
!> @param[in] map_chi      Dofmap for the cell at the base of the column for chi
!> @param[in] chi_basis    Basis functions evaluated at gaussian quadrature point
subroutine create_w3mask_code( nlayers,    &
                               boundary_s, &
                               boundary_n, &
                               boundary_e, &
                               boundary_w, &
                               rho_mask,   &
                               chi_1,      &
                               chi_2,      &
                               chi_3,      &
                               ndf_w3,     &
                               undf_w3,    &
                               map_w3,     &
                               ndf_chi,    &
                               undf_chi,   &
                               map_chi,    &
                               chi_basis )

  implicit none

  ! Arguments
  integer(kind=i_def),                           intent(in) :: nlayers,        &
                                                               ndf_w3,         &
                                                               ndf_chi,        &
                                                               undf_w3,        &
                                                               undf_chi
  integer(kind=i_def), dimension(ndf_w3),        intent(in) :: map_w3
  integer(kind=i_def), dimension(ndf_chi),       intent(in) :: map_chi
  real(kind=r_def),                              intent(in) :: boundary_s,     &
                                                               boundary_n,     &
                                                               boundary_e,     &
                                                               boundary_w
  real(kind=r_def), dimension(undf_w3),          intent(inout)  :: rho_mask
  real(kind=r_def), dimension(1,ndf_chi,ndf_w3), intent(in)     :: chi_basis
  real(kind=r_def), dimension(undf_chi),         intent(in)     :: chi_1,      &
                                                                   chi_2,      &
                                                                   chi_3

  ! Internal variables
  integer(kind=i_def)                    :: k, df1, df2
  real(kind=r_def), dimension(ndf_chi)   :: chi_1_e, chi_2_e, chi_3_e
  real(kind=r_def)                       :: x(3),           &
                                            query_value_ns, &
                                            query_value_ew, &
                                            mask_value
  logical(kind=l_def)                    :: innerbox

  ! Define the mask for the bottom vertical layer and copy to others
  k = 0

  ! Get the chi values in the kth layer
  ! chi_e (u) = chi ( f(u) + k )
  do df2 = 1, ndf_chi
    chi_1_e(df2) = chi_1( map_chi(df2) + k)
    chi_2_e(df2) = chi_2( map_chi(df2) + k)
    chi_3_e(df2) = chi_3( map_chi(df2) + k)
  end do

  ! Get the x values (x,y,z) for each cell centre, by summing over
  ! W0 to W3 basis functions
  ! x(v) = sum_u chi_e(u) * g(u,v)
  do df1 = 1, ndf_w3
    x(:) = 0.0_r_def
    do df2 = 1,ndf_chi
      x(1) = x(1) + chi_1_e(df2) * chi_basis(1,df2,df1)
      x(2) = x(2) + chi_2_e(df2) * chi_basis(1,df2,df1)
      x(3) = x(3) + chi_3_e(df2) * chi_basis(1,df2,df1)
    end do

    ! Change of coordinates for spherical geometry
    if ( geometry == geometry_spherical .and. &
         coord_system == coord_system_xyz ) then
      ! in alpha beta space - and only pick the face where x(2)>0
      if ( x(2) > 0.0_r_def ) then
        query_value_ns = atan2(x(3),x(2))
        query_value_ew = atan2(x(1),x(2))
      else
       ! set to an arbitrary value bigger than pi/4
        query_value_ns = 10.0_r_def
        query_value_ew = 10.0_r_def
      endif
    else
      ! in cartesian space
      query_value_ns = x(2)
      query_value_ew = x(1)
    end if

    ! Define the limited area domain: 1 in interior, 0 in exterior

    ! 00000000000000000
    ! 00---innerbox--00
    ! 00|11111111111|00
    ! 00|11111111111|00
    ! 00|11111111111|00
    ! 00-------------00
    ! 00000000000000000

    ! Is the dof inside the innerbox?
    if ( query_value_ns > boundary_s .and. &
         query_value_ns < boundary_n .and. &
         query_value_ew > boundary_w .and. &
         query_value_ew < boundary_e ) then
      innerbox = .true.
    else
      innerbox = .false.
    endif

    ! Define the mask_value
    mask_value = 1.0_r_def
    if (.not. innerbox) then
      mask_value = 0.0_r_def
    endif

    ! Loop over vertical layers
    do k = 0, nlayers-1
      rho_mask(map_w3(df1) + k) = mask_value
    enddo

  end do

end subroutine create_w3mask_code

end module create_w3mask_kernel_mod
