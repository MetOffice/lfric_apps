!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes the rhs for the initialisation of the wind field.
!>
!> The kernel computes the rhs of the equation u = u0 where u0 is the
!> analytically defined wind field. The analytic wind field is projected onto
!> using Galerkin projection.
!>
module initial_u_kernel_mod

  use argument_mod,            only : arg_type, func_type,       &
                                      GH_FIELD, GH_SCALAR,       &
                                      GH_INC, GH_READ,           &
                                      GH_REAL, ANY_SPACE_9,      &
                                      ANY_DISCONTINUOUS_SPACE_3, &
                                      GH_BASIS, GH_DIFF_BASIS,   &
                                      CELL_COLUMN, GH_QUADRATURE_XYoZ
  use constants_mod,           only : r_def, i_def, PI
  use fs_continuity_mod,       only : W2
  use initial_wind_config_mod, only : profile_sin_uv,                        &
                                      profile, sbr_angle_lat, sbr_angle_lon, &
                                      u0, v0, shear, wavelength
  use kernel_mod,              only : kernel_type
  use log_mod,                 only : log_event, LOG_LEVEL_ERROR

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: initial_u_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                     &
         arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),                        &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ, ANY_SPACE_9),               &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),                            &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ),                            &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)                             &
         /)
    type(func_type) :: meta_funcs(2) = (/                                   &
         func_type(W2,          GH_BASIS),                                  &
         func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS)                    &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: initial_u_code
  end type

!-----------------------------------------------------------------------------
! Contained functions/subroutines
!-----------------------------------------------------------------------------
public :: initial_u_code

contains

  !> @brief Compute the right hand side to initialise the wind field.
  !! @param[in] nlayers Number of layers
  !! @param[in,out] rhs Right hand side field to compute
  !! @param[in] chi_1 1st coordinate field
  !! @param[in] chi_2 2nd coordinate field
  !! @param[in] chi_3 3rd coordinate field
  !! @param[in] panel_id A field giving the ID for mesh panels.
  !! @param[in] time Time (timestep multiplied by dt)
  !! @param[in] domain_max_x Domain maximum x-coordinate.
  !! @param[in] domain_max_y Domain maximum y-coordinate.
  !! @param[in] ndf Number of degrees of freedom per cell for W2
  !! @param[in] undf Total number of degrees of freedom for W2
  !! @param[in] map Dofmap for the cell at the base of the column for W2
  !! @param[in] basis Basis functions for W2 evaluated at gaussian quadrature points
  !! @param[in] ndf_chi Number of degrees of freedom per cell for chi
  !! @param[in] undf_chi Number of unique degrees of freedom for chi
  !! @param[in] map_chi Dofmap for the cell at the base of the column for chi
  !! @param[in] chi_basis Basis functions for Wchi evaluated at
  !!                      gaussian quadrature points
  !! @param[in] chi_diff_basis Differential of the Wchi basis functions
  !!                           evaluated at gaussian quadrature point
  !! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
  !! @param[in] undf_pid Number of unique degrees of freedom for panel_id
  !! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
  !! @param[in] nqp_h Number of quadrature points in the horizontal
  !! @param[in] nqp_v Number of quadrature points in the vertical
  !! @param[in] wqp_h Horizontal quadrature weights
  !! @param[in] wqp_v Vertical quadrature weights
  subroutine initial_u_code(nlayers,                         &
                            rhs,                             &
                            chi_1, chi_2, chi_3,             &
                            panel_id,                        &
                            time,                            &
                            domain_max_x,                    &
                            domain_max_y,                    &
                            ndf, undf, map, basis,           &
                            ndf_chi, undf_chi,               &
                            map_chi, chi_basis,              &
                            chi_diff_basis,                  &
                            ndf_pid, undf_pid, map_pid,      &
                            nqp_h, nqp_v, wqp_h, wqp_v       &
                            )

  use analytic_wind_profiles_mod, only : analytic_wind
  use base_mesh_config_mod,       only : geometry,           &
                                         geometry_spherical
  use chi_transform_mod,          only : chi2llr
  use coordinate_jacobian_mod,    only : coordinate_jacobian
  use coord_transform_mod,        only : sphere2cart_vector

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf, ndf_pid, ndf_chi
  integer(kind=i_def), intent(in) :: undf, undf_pid, undf_chi
  integer(kind=i_def), intent(in) :: nqp_h, nqp_v

  integer(kind=i_def), dimension(ndf),         intent(in) :: map
  integer(kind=i_def), dimension(ndf_chi), intent(in) :: map_chi
  integer(kind=i_def), dimension(ndf_pid),     intent(in) :: map_pid

  real(kind=r_def), intent(in), dimension(3,ndf,    nqp_h,nqp_v) :: basis
  real(kind=r_def), intent(in), dimension(3,ndf_chi,nqp_h,nqp_v) :: chi_diff_basis
  real(kind=r_def), intent(in), dimension(1,ndf_chi,nqp_h,nqp_v) :: chi_basis

  real(kind=r_def), dimension(undf),      intent(inout) :: rhs
  real(kind=r_def), dimension(undf_chi),  intent(in)    :: chi_1, &
                                                           chi_2, &
                                                           chi_3
  real(kind=r_def), dimension(undf_pid),  intent(in)    :: panel_id
  real(kind=r_def), intent(in)                          :: time
  real(kind=r_def), intent(in)                          :: domain_max_x
  real(kind=r_def), intent(in)                          :: domain_max_y
  real(kind=r_def), dimension(nqp_h), intent(in)        :: wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)        :: wqp_v


  integer(kind=i_def), parameter :: n_options = 3
  real(kind=r_def)               :: opt_args(n_options)

  ! Internal variables
  integer(kind=i_def)                          :: df, k, qp1, qp2
  real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
  real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jacobian
  real(kind=r_def), dimension(ndf_chi)     :: chi_1_cell, &
                                                  chi_2_cell, &
                                                  chi_3_cell
  real(kind=r_def), dimension(3)               :: u_physical, u_spherical, coords, llr
  real(kind=r_def)                             :: integrand

  integer(kind=i_def) :: ipanel

  ipanel = int(panel_id(map_pid(1)), i_def)

  if ( geometry == geometry_spherical ) then
    ! Options for spherical domains
    opt_args = [ u0, sbr_angle_lat, sbr_angle_lon ]
  else
    ! Options for Cartesian domains
    if ( profile == profile_sin_uv ) then
      opt_args = [ u0, v0, wavelength ]
    else
      opt_args = [ u0, v0, shear ]
    end if
  end if

  do k = 0, nlayers-1
    do df = 1, ndf_chi
      chi_1_cell(df) = chi_1( map_chi(df) + k)
      chi_2_cell(df) = chi_2( map_chi(df) + k)
      chi_3_cell(df) = chi_3( map_chi(df) + k)
    end do

    call coordinate_jacobian(ndf_chi,        &
                             nqp_h,          &
                             nqp_v,          &
                             chi_1_cell,     &
                             chi_2_cell,     &
                             chi_3_cell,     &
                             ipanel,         &
                             chi_basis,      &
                             chi_diff_basis, &
                             jacobian,       &
                             dj)

    do qp2 = 1, nqp_v
      do qp1 = 1, nqp_h
        ! Compute analytical vector wind in physical space
        coords(:) = 0.0_r_def
        do df = 1, ndf_chi
          coords(1) = coords(1) + chi_1_cell(df)*chi_basis(1,df,qp1,qp2)
          coords(2) = coords(2) + chi_2_cell(df)*chi_basis(1,df,qp1,qp2)
          coords(3) = coords(3) + chi_3_cell(df)*chi_basis(1,df,qp1,qp2)
        end do

        if ( geometry == geometry_spherical ) then

          ! Need to obtain longitude, latitude and radius from position vector
          call chi2llr(coords(1), coords(2), coords(3), &
                       ipanel, llr(1), llr(2), llr(3))

          ! Obtain (lon,lat,r) components of u and then transform to (X,Y,Z) components
          u_spherical = analytic_wind( llr, time, profile, n_options, &
                                       domain_max_x, domain_max_y,    &
                                       opt_args )
          u_physical  = sphere2cart_vector(u_spherical,llr)
        else
          ! For planar geometries, no transformation is needed
          u_physical = analytic_wind( coords, time, profile, n_options, &
                                      domain_max_x, domain_max_y,       &
                                      opt_args )
        end if

        do df = 1, ndf
          integrand = dot_product(matmul(jacobian(:,:,qp1,qp2),&
                                         basis(:,df,qp1,qp2)),u_physical)
          rhs(map(df) + k) = rhs(map(df) + k) &
                           + wqp_h(qp1)*wqp_v(qp2)*integrand
        end do
      end do
    end do
  end do

end subroutine initial_u_code

end module initial_u_kernel_mod
