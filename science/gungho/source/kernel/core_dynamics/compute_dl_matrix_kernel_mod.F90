!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief A 'vertical W2' mass matrix used in the damping layer term.
!>
!> @details The kernel modifies the two rows of the velocity mass matrix
!> corresponding to the degrees of freedom of the vertical component of velocity
!> to account for Rayleigh damping in the absorbing layer region for use in runs
!> with non-flat bottom boundary.
!>
module compute_dl_matrix_kernel_mod

  use argument_mod,              only: arg_type, func_type,       &
                                       GH_OPERATOR, GH_FIELD,     &
                                       GH_READ, GH_WRITE,         &
                                       GH_REAL, ANY_SPACE_9,      &
                                       ANY_DISCONTINUOUS_SPACE_3, &
                                       GH_BASIS, GH_DIFF_BASIS,   &
                                       GH_SCALAR, GH_INTEGER,     &
                                       CELL_COLUMN, GH_QUADRATURE_XYoZ
  use base_mesh_config_mod,      only: geometry, geometry_spherical
  use constants_mod,             only: i_def, r_def, r_second, &
                                       PI, degrees_to_radians
  use chi_transform_mod,         only: chi2llr
  use damping_layer_config_mod,  only: dl_type, dl_type_latitude
  use finite_element_config_mod, only: coord_system,    &
                                       coord_system_xyz
  use fs_continuity_mod,         only: W2
  use kernel_mod,                only: kernel_type
  use coordinate_jacobian_mod,   only: coordinate_jacobian

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------

  type, public, extends(kernel_type) :: compute_dl_matrix_kernel_type
    private
    type(arg_type) :: meta_args(9) = (/                                       &
         arg_type(GH_OPERATOR, GH_REAL, GH_WRITE, W2, W2),                    &
         arg_type(GH_FIELD*3,  GH_REAL, GH_READ,  ANY_SPACE_9),               &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,  ANY_DISCONTINUOUS_SPACE_3), &
         arg_type(GH_SCALAR,   GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR,   GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR,   GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR,   GH_REAL, GH_READ),                             &
         arg_type(GH_SCALAR,   GH_INTEGER, GH_READ),                          &
         arg_type(GH_SCALAR,   GH_REAL, GH_READ)                              &
         /)
    type(func_type) :: meta_funcs(2) = (/                                    &
         func_type(ANY_SPACE_9, GH_BASIS, GH_DIFF_BASIS),                    &
         func_type(W2,          GH_BASIS)                                    &
         /)
        integer :: operates_on = CELL_COLUMN
        integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass :: compute_dl_matrix_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: compute_dl_matrix_code
  public damping_layer_func

contains

  !> @brief Computes the reduced mass matrix for the damping layer term in the momentum equation.
  !!
  !! @param[in] cell     Identifying number of cell
  !! @param[in] nlayers  Number of layers
  !! @param[in] ncell_3d ncell*ndf
  !! @param[in,out] mm   Local stencil or mass matrix
  !! @param[in] chi1     1st coordinate field in Wchi
  !! @param[in] chi2     2nd coordinate field in Wchi
  !! @param[in] chi3     3rd coordinate field in Wchi
  !! @param[in] panel_id Field giving the ID for mesh panels
  !! @param[in] dl_base_height
  !!                     Base height of damping layer
  !! @param[in] dl_strength
  !!                     Strength of damping layer
  !! @param[in] domain_top
  !!                     The model domain height
  !! @param[in] radius   The planet radius
  !! @param[in] element_order The model finite element order
  !! @param[in] dt       The model timestep length
  !! @param[in] ndf_w2   Degrees of freedom per cell
  !! @param[in] basis_w2 Vector basis functions evaluated at quadrature points.
  !! @param[in] ndf_chi  Degrees of freedom per cell for chi field
  !! @param[in] undf_chi Unique degrees of freedom for chi field
  !! @param[in] map_chi  Dofmap for the cell at the base of the column, for the
  !!                     space on which the chi field lives
  !! @param[in] basis_chi Vector basis functions evaluated at quadrature points
  !! @param[in] diff_basis_chi Vector differential basis functions evaluated at
  !!                           quadrature points
  !! @param[in] ndf_pid  Number of degrees of freedom per cell for panel_id
  !! @param[in] undf_pid Number of unique degrees of freedom for panel_id
  !! @param[in] map_pid  Dofmap for the cell at the base of the column for panel_id
  !! @param[in] nqp_h    Number of horizontal quadrature points
  !! @param[in] nqp_v    Number of vertical quadrature points
  !! @param[in] wqp_h    Horizontal quadrature weights
  !! @param[in] wqp_v    Vertical quadrature weights
  subroutine compute_dl_matrix_code(cell, nlayers, ncell_3d,     &
                                    mm, chi1, chi2, chi3,        &
                                    panel_id, dl_base_height,    &
                                    dl_strength, domain_top,     &
                                    radius, element_order, dt,   &
                                    ndf_w2, basis_w2,            &
                                    ndf_chi, undf_chi, map_chi,  &
                                    basis_chi, diff_basis_chi,   &
                                    ndf_pid, undf_pid,  map_pid, &
                                    nqp_h, nqp_v, wqp_h, wqp_v)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: cell, nqp_h, nqp_v
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ncell_3d
    integer(kind=i_def), intent(in)    :: ndf_w2
    integer(kind=i_def), intent(in)    :: ndf_chi
    integer(kind=i_def), intent(in)    :: undf_chi
    integer(kind=i_def), intent(in)    :: ndf_pid, undf_pid
    integer(kind=i_def), intent(in)    :: map_chi(ndf_chi)
    integer(kind=i_def), intent(in)    :: map_pid(ndf_pid)

    real(kind=r_def),    intent(inout) :: mm(ndf_w2,ndf_w2,ncell_3d)
    real(kind=r_def),    intent(in)    :: basis_chi(1,ndf_chi,nqp_h,nqp_v)
    real(kind=r_def),    intent(in)    :: diff_basis_chi(3,ndf_chi,nqp_h,nqp_v)
    real(kind=r_def),    intent(in)    :: chi1(undf_chi)
    real(kind=r_def),    intent(in)    :: chi2(undf_chi)
    real(kind=r_def),    intent(in)    :: chi3(undf_chi)
    real(kind=r_def),    intent(in)    :: panel_id(undf_pid)
    real(kind=r_def),    intent(in)    :: dl_base_height
    real(kind=r_def),    intent(in)    :: dl_strength
    real(kind=r_def),    intent(in)    :: domain_top
    real(kind=r_def),    intent(in)    :: radius
    real(kind=r_second), intent(in)    :: dt
    real(kind=r_def),    intent(in)    :: wqp_h(nqp_h)
    real(kind=r_def),    intent(in)    :: wqp_v(nqp_v)
    real(kind=r_def),    intent(in)    :: basis_w2(3,ndf_w2,nqp_h,nqp_v)
    integer(kind=i_def), intent(in)    :: element_order

    ! Internal variables
    integer(kind=i_def) :: df, df2, dfc, k, ik
    integer(kind=i_def) :: qp1, qp2

    real(kind=r_def), dimension(ndf_chi)         :: chi1_e, chi2_e, chi3_e
    real(kind=r_def)                             :: integrand
    real(kind=r_def)                             :: mu_at_quad
    real(kind=r_def)                             :: chi1_at_quad
    real(kind=r_def)                             :: chi2_at_quad
    real(kind=r_def)                             :: chi3_at_quad
    real(kind=r_def)                             :: lat_at_quad
    real(kind=r_def)                             :: long_at_quad
    real(kind=r_def)                             :: r_at_quad
    real(kind=r_def)                             :: z
    real(kind=r_def), dimension(nqp_h,nqp_v)     :: dj
    real(kind=r_def), dimension(3,3,nqp_h,nqp_v) :: jac

    integer(kind=i_def) :: ipanel

    ipanel = int(panel_id(map_pid(1)), i_def)

    ! Loop over layers: Start from 1 as in this loop k is not an offset
    do k = 1, nlayers
      ik = k + (cell-1)*nlayers
      ! Indirect the chi coord field here
      do df = 1, ndf_chi
        chi1_e(df) = chi1(map_chi(df) + k - 1)
        chi2_e(df) = chi2(map_chi(df) + k - 1)
        chi3_e(df) = chi3(map_chi(df) + k - 1)
      end do
      call coordinate_jacobian(ndf_chi, nqp_h, nqp_v, chi1_e, chi2_e, chi3_e,  &
                               ipanel, basis_chi, diff_basis_chi, jac, dj)

      ! Only use dofs corresponding to vertical part of basis function
      do df2 = 1, ndf_w2
        do df = 1, ndf_w2 ! Mass matrix is not symmetric for damping layer
          mm(df,df2,ik) = 0.0_r_def
          do qp2 = 1, nqp_v
            do qp1 = 1, nqp_h
              chi1_at_quad = 0.0_r_def
              chi2_at_quad = 0.0_r_def
              chi3_at_quad = 0.0_r_def
              do dfc = 1,ndf_chi
                chi1_at_quad = chi1_at_quad + chi1_e(dfc)*basis_chi(1,dfc,qp1,qp2)
                chi2_at_quad = chi2_at_quad + chi2_e(dfc)*basis_chi(1,dfc,qp1,qp2)
                chi3_at_quad = chi3_at_quad + chi3_e(dfc)*basis_chi(1,dfc,qp1,qp2)
              end do

              if (geometry == geometry_spherical) then

                call chi2llr(chi1_at_quad,chi2_at_quad,chi3_at_quad, &
                             ipanel, long_at_quad,lat_at_quad,r_at_quad)
                z=r_at_quad - radius

                if (dl_type == dl_type_latitude) then
                  mu_at_quad = damping_layer_func(z, dl_strength, &
                                                  dl_base_height, domain_top, lat_at_quad)
                else
                  mu_at_quad = damping_layer_func(z, dl_strength, &
                                                  dl_base_height, domain_top, 0.0_r_def)
                end if
              else
                mu_at_quad = damping_layer_func(chi3_at_quad, dl_strength, &
                                                dl_base_height, domain_top, 0.0_r_def)
              end if

              integrand = wqp_h(qp1) * wqp_v(qp2) *                          &
                          dot_product(                                       &
                          matmul(jac(:,:,qp1,qp2),basis_w2(:,df,qp1,qp2)),   &
                          matmul(jac(:,:,qp1,qp2),basis_w2(:,df2,qp1,qp2)) ) &
                          /dj(qp1,qp2)
              ! Only modify dofs corresponding to vertical part of w-basis function (for lowest order: bottom two rows)
              if (df > ndf_w2 - (element_order+2)*(element_order+1)**2) then
                mm(df,df2,ik) = mm(df,df2,ik) +                          &
                                (1.0_r_def + real(dt, r_def)*mu_at_quad) &
                                * integrand
              else
                mm(df,df2,ik) = mm(df,df2,ik) + integrand
              end if
            end do
          end do
        end do
      end do

    end do ! End of k loop

  end subroutine compute_dl_matrix_code

  !> @brief Computes the value of damping layer function at height z.
  !!
  !! @param[in]  height         Physical height on which to compute value of damping layer function.
  !! @param[in]  dl_strength    Damping layer function coefficient and maximum value of damping layer function.
  !! @param[in]  dl_base_height Height above which damping layer is active.
  !! @param[in]  domain_top     Top of the computational domain.
  !! @param[in]  latitude     Latitude at which to computer the damping function
  !! @return     dl_val         Value of damping layer function.
  function damping_layer_func(height, dl_strength, dl_base_height, domain_top, latitude) result (dl_val)

    implicit none

    ! Arguments
    real(kind=r_def),   intent(in)  :: height
    real(kind=r_def),   intent(in)  :: dl_strength
    real(kind=r_def),   intent(in)  :: dl_base_height
    real(kind=r_def),   intent(in)  :: domain_top
    real(kind=r_def),   intent(in)  :: latitude
    real(kind=r_def)                :: dl_val
    real(kind=r_def)                :: height_star
    ! Latitude over which to reduce damping height from dl_base_height
    ! chosen to reduce it to approximately dl_base_height/2
    real(kind=r_def), parameter     :: taper_lat = 50.0_r_def

    if (abs(latitude) < taper_lat*degrees_to_radians) then
      height_star = domain_top + (height-domain_top)*cos(latitude)
    else
      height_star = domain_top + (height-domain_top)*cos(taper_lat*degrees_to_radians)
    end if

    if (height_star >= dl_base_height) then
      dl_val = dl_strength*sin(0.5_r_def*PI*(height_star-dl_base_height) / &
                                            (domain_top-dl_base_height))**2
    else
      dl_val = 0.0_r_def
    end if

  end function damping_layer_func

end module compute_dl_matrix_kernel_mod
