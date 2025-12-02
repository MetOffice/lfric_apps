!-----------------------------------------------------------------------------
! Copyright (c) 2017,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-----------------------------------------------------------------------------
!> @brief Computes the boundary integral part of the projection operator from
!>        the velocity space to the potential temperature space weighted by
!>        the potential temperature gradient.
!>
!> @details The kernel computes the boundary integral part of the projection
!>          operator from the velocity space to the potential temperature space
!>          weighted by the potential temperature gradient.
!>          Computes the projection operator
!>          \f<[<\gamma, flux(\theta*v)\cdot n>\f],
!>          where v is in W2 and gamma is in the potential temperature space.
!>
!> @todo Create unit test for this kernel, see #2935
module weighted_proj_theta2_bd_kernel_mod

  use argument_mod,      only : arg_type, func_type,         &
                                mesh_data_type,              &
                                reference_element_data_type, &
                                GH_OPERATOR, GH_FIELD,       &
                                GH_SCALAR, GH_REAL,          &
                                GH_READ, GH_READWRITE,       &
                                STENCIL, CROSS, GH_BASIS,    &
                                GH_QUADRATURE_face,          &
                                CELL_COLUMN, adjacent_face,  &
                                outward_normals_to_horizontal_faces
  use constants_mod,     only : r_def, i_def, l_def, r_solver
  use cross_product_mod, only : cross_product
  use fs_continuity_mod, only : W2, Wtheta
  use kernel_mod,        only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> PSy layer.
  type, public, extends(kernel_type) :: weighted_proj_theta2_bd_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                        &
         arg_type(GH_OPERATOR, GH_REAL, GH_READWRITE, Wtheta, W2),             &
         arg_type(GH_FIELD,    GH_REAL, GH_READ,      Wtheta, STENCIL(CROSS)), &
         arg_type(GH_SCALAR,   GH_REAL, GH_READ)                               &
         /)
    type(func_type) :: meta_funcs(2) = (/                                      &
         func_type(Wtheta, GH_BASIS),                                          &
         func_type(W2,     GH_BASIS)                                           &
         /)
    type(mesh_data_type) :: meta_mesh(1) = (/                                  &
         mesh_data_type( adjacent_face )                                       &
         /)
    type(reference_element_data_type) :: meta_reference_element(1) = (/        &
         reference_element_data_type( outward_normals_to_horizontal_faces )    &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_face
  contains
    procedure, nopass :: weighted_proj_theta2_bd_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: weighted_proj_theta2_bd_code

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Computes the weighted projection from W2 to Wtheta.
  !>
  !> @param[in] cell Cell number
  !> @param[in] nlayers Number of layers
  !> @param[in] ncell_3d ncell*ndf
  !> @param[in,out] projection Projection operator to compute
  !> @param[in] theta Potential temperature
  !> @param[in] stencil_wtheta_size Number of cells in the Wtheta stencil map
  !> @param[in] stencil_wtheta_map Stencil dofmap for the Wtheta space
  !> @param[in] scalar Real to scale matrix by
  !> @param[in] ndf_wtheta Number of degrees of freedom per cell for Wtheta
  !> @param[in] undf_wtheta Number of unique degrees of freedom for Wtheta
  !> @param[in] map_wtheta Dofmap for the cell at the base of the column for
  !!                       Wtheta
  !> @param[in] wtheta_basis_face Wtheta basis functions evaluated at Gaussian
  !!                              quadrature points on horizontal faces
  !> @param[in] ndf_w2 Number of degrees of freedom per cell for W2
  !> @param[in] w2_basis_face W2 basis functions evaluated at Gaussian
  !!                          quadrature points on horizontal faces
  !> @param[in] nfaces_re_h Number of reference element faces bisected by a
  !!                        horizontal plane
  !> @param[in] outward_normals Vector of normals to the
  !!                            reference element horizontal
  !!                            "outward faces"
  !> @param[in] adjacent_face Vector containing information on neighbouring face
  !!                          index for the current cell
  !> @param[in] nfaces_qr Number of faces in the quadrature rule
  !> @param[in] nqp_f Number of quadrature points on horizontal faces
  !> @param[in] wqp_f Quadrature weights on horizontal faces
  !>
  subroutine weighted_proj_theta2_bd_code( cell, nlayers, ncell_3d,             &
                                           projection, theta,                   &
                                           stencil_wtheta_size,                 &
                                           stencil_wtheta_map,                  &
                                           scalar,                              &
                                           ndf_wtheta, undf_wtheta, map_wtheta, &
                                           wtheta_basis_face,                   &
                                           ndf_w2, w2_basis_face,               &
                                           nfaces_re_h,                         &
                                           outward_normals,                     &
                                           adjacent_face,                       &
                                           nfaces_qr, nqp_f, wqp_f )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: cell, ncell_3d
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: nfaces_qr, nqp_f
    integer(kind=i_def), intent(in) :: ndf_wtheta, undf_wtheta, ndf_w2
    integer(kind=i_def), intent(in) :: nfaces_re_h

    integer(kind=i_def), intent(in) :: stencil_wtheta_size
    integer(kind=i_def), dimension(ndf_wtheta, stencil_wtheta_size), intent(in) :: stencil_wtheta_map

    integer(kind=i_def), dimension(ndf_wtheta), intent(in) :: map_wtheta

    real(kind=r_solver), dimension(ncell_3d,ndf_wtheta,ndf_w2), intent(inout) :: projection
    real(kind=r_solver), dimension(undf_wtheta),                intent(in)    :: theta
    real(kind=r_solver),                                        intent(in)    :: scalar

    real(kind=r_def), dimension(3,ndf_w2,nqp_f,nfaces_qr),     intent(in) :: w2_basis_face
    real(kind=r_def), dimension(1,ndf_wtheta,nqp_f,nfaces_qr), intent(in) :: wtheta_basis_face
    real(kind=r_def), dimension(nqp_f,nfaces_qr),              intent(in) :: wqp_f

    integer(kind=i_def), intent(in) :: adjacent_face(:)

    real(kind=r_def),    intent(in) :: outward_normals(:,:)

    ! Internal variables
    integer(kind=i_def) :: df, k, ik, face, face_next, dft, df2
    integer(kind=i_def) :: qp

    real(kind=r_solver), dimension(ndf_wtheta) :: theta_e, theta_next_e
    real(kind=r_solver) :: theta_at_fquad, theta_next_at_fquad, v_dot_n
    real(kind=r_solver) :: flux_term, theta_av

    logical(kind=l_def), parameter :: upwind = .false.

    real(kind=r_solver), dimension(3,ndf_w2,nqp_f,nfaces_qr)     :: rsol_w2_basis_face
    real(kind=r_solver), dimension(1,ndf_wtheta,nqp_f,nfaces_qr) :: rsol_wtheta_basis_face
    real(kind=r_solver), dimension(nqp_f,nfaces_qr)              :: rsol_wqp_f

    ! If we're near the edge of the regional domain then the
    ! stencil size will be less that 5 so don't do anything here
    ! This should be removed with lfric ticket #2958
    if (stencil_wtheta_size < 5)return

    rsol_w2_basis_face     = real(w2_basis_face, r_solver)
    rsol_wtheta_basis_face = real(wtheta_basis_face, r_solver)
    rsol_wqp_f             = real(wqp_f, r_solver)

    ! Assumes same number of horizontal qp in x and y
    do k = 0, nlayers-1
      ik = k + 1 + (cell-1)*nlayers

      do face = 1, nfaces_re_h
        ! Storing opposite face number on neighbouring cell
        face_next = adjacent_face(face)

        ! Computing theta in adjacent cells
        do df = 1, ndf_wtheta
          theta_e(df)      = theta(stencil_wtheta_map(df, 1) + k )
          theta_next_e(df) = theta(stencil_wtheta_map(df, face+1) + k )
        end do

        do qp = 1, nqp_f

          theta_at_fquad = 0.0_r_solver
          theta_next_at_fquad = 0.0_r_solver
          do df = 1, ndf_wtheta
            theta_at_fquad       = theta_at_fquad      + theta_e(df)     * &
                                   rsol_wtheta_basis_face(1,df,qp,face)
            theta_next_at_fquad  = theta_next_at_fquad + theta_next_e(df)* &
                                   rsol_wtheta_basis_face(1,df,qp,face_next)
          end do
          theta_av = 0.5_r_solver * (theta_at_fquad + theta_next_at_fquad)

          do df2 = 1,ndf_w2
            v_dot_n = dot_product(rsol_w2_basis_face(:,df2,qp,face), &
                                  real(outward_normals(:, face), r_solver))
            flux_term = rsol_wqp_f(qp,face) * theta_av * v_dot_n
            if (upwind) then
              flux_term = flux_term + 0.5_r_solver * abs(v_dot_n) * &
                          (theta_at_fquad - theta_next_at_fquad)
            end if
            do dft = 1,ndf_wtheta
              projection(ik,dft,df2) = projection(ik,dft,df2) &
                                     + rsol_wtheta_basis_face(1,dft,qp,face) &
                                     * flux_term * scalar

            end do ! dft
          end do ! df2

        end do ! qp

      end do ! faces

    end do ! layers

  end subroutine weighted_proj_theta2_bd_code

end module weighted_proj_theta2_bd_kernel_mod
