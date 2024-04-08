!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the facet integral for the rhs of the thermal shallow water
!!        vector invariant momentum equation, for buoyancy in V2 (W3)
!!        and zero surface geopotential.
!!
!> @details The rhs of the thermal shallow water vector invariant momentum equation
!!          consists of three terms:
!!          buoyancy gradient: \f[ \frac{1}{2} \Phi \nabla(s')\f]
!!          geopotential gradient: \f[ \nabla(s' \Phi)\f]
!!          gradient of kinetic energy: \f[ \nabla(\frac{1}{2} u \cdot u) \f]
!!
!!          The gradient of kinetic energy is computed in swe_rhs_alg_mod.
!!          Multiplying the buoyancy and geopotential gradients by the
!!          test function (v) and using the reverse product rule gives
!!
!!          \f[ s v \cdot \nabla \Phi + \frac{1}{2}\Phi v \cdot \nabla s =
!!              - s \Phi \nabla \cdot v + \frac{1}{2} s \nabla \cdot (v \Phi)
!!              + \frac{1}{2} \nabla \cdot (v s \Phi)  \f]
!!
!!          The first two terms are evaluated as a volume integrals and are handled
!!          in swe_buoyancy_gradient_volume_kernel_mod.
!!          The third term can be integrated over each element, then the sum taken.
!!          Applying the divergence theorem, and using that geopotential and
!!          buoyancy are in V2 (W3) and thus discontinuous, gives
!!
!!          \f[ \frac{1}{2} \int [\Phi v]\{s\} \cdot \hat{n} dS\f]
!!
!!          This is evaluated here as 0.5*geopotential*v*normal_vector*average(buoyancy).
!!          Note that this kernel is only suitable for lowest order function
!!          space, zero surface geopotential, and planar geometry.
!!
!> @todo Modify kernel to work for non-zero surface geopotential and
!!       and any geometry, see ticket #3002.

module swe_buoyancy_gradient_facet_kernel_mod

  use argument_mod,             only : arg_type, func_type,         &
                                       mesh_data_type,              &
                                       reference_element_data_type, &
                                       GH_FIELD, GH_READ, GH_INC,   &
                                       GH_BASIS,                    &
                                       GH_DIFF_BASIS, CELL_COLUMN,  &
                                       GH_QUADRATURE_face,          &
                                       adjacent_face,               &
                                       outward_normals_to_faces,    &
                                       normals_to_horizontal_faces, &
                                       GH_REAL, STENCIL, CROSS
  use constants_mod,            only : r_def, i_def
  use cross_product_mod,        only : cross_product
  use fs_continuity_mod,        only : W2, W3
  use kernel_mod,               only : kernel_type

  implicit none

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the Psy layer
  type, public, extends(kernel_type) :: swe_buoyancy_gradient_facet_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                 &
      arg_type(GH_FIELD, GH_REAL, GH_INC,  W2),                         & ! r_u_bd
      arg_type(GH_FIELD, GH_REAL, GH_READ, W3),                         & ! geopot
      arg_type(GH_FIELD, GH_REAL, GH_READ, W3, STENCIL(CROSS))          & ! buoyancy
      /)
    type(func_type) :: meta_funcs(2) = (/                               &
      func_type(W2, GH_BASIS),                                          &
      func_type(W3, GH_BASIS)                                           &
      /)
    type(mesh_data_type) :: meta_mesh(1) = (/                           &
        mesh_data_type( adjacent_face )                                 &
      /)
    type(reference_element_data_type) :: meta_reference_element(2) = (/ &
         reference_element_data_type( normals_to_horizontal_faces ),    &
         reference_element_data_type( outward_normals_to_faces )        &
      /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_face
  contains
    procedure, nopass ::swe_buoyancy_gradient_facet_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public swe_buoyancy_gradient_facet_code
contains

  !> @brief Compute the facet integral for the for the rhs of the thermal
  !!        shallow water momentum equation.
  !> @param[in]     nlayers         Number of layers
  !> @param[in,out] r_u_bd          Right hand side of the momentum equation
  !> @param[in]     geopot          The geopotential
  !> @param[in]     buoyancy        The buoyancy
  !> @param[in]     stencil_size_w3 Size of the W3 stencil (number of cells)
  !> @param[in]     stencil_map_w3  W3 dofmaps for the stencil
  !> @param[in]     ndf_w2          Number of degrees of freedom per cell for w2
  !> @param[in]     undf_w2         Number unique of degrees of freedom  for w2
  !> @param[in]     map_w2          Integer array holding the dofmap for the cell for w2
  !> @param[in]     w2_basis_face   Basis functions evaluated at gaussian quadrature
  !!                                points on horizontal faces
  !> @param[in]     ndf_w3          Number of degrees of freedom per cell for w3
  !> @param[in]     undf_w3         Number unique of degrees of freedom  for w3
  !> @param[in]     map_w3          W3 dofmaps for the cell for w3
  !> @param[in]     w3_basis_face   Basis functions evaluated at gaussian quadrature
  !!                                points on horizontal faces
  !> @param[in]     nfaces_re_h     Number of horizontal faces
  !> @param[in]     nfaces_re       Number of faces
  !> @param[in]     normals_to_horizontal_faces
  !!                                Vector normal to the horizontal faces of the reference element.
  !> @param[in]     outward_normals_to_faces
  !!                                Vector normal to the out faces of the reference element.
  !> @param[in]     adjacent_face   Vector containing information on neighbouring
  !!                                face index for the current cell
  !> @param[in]     nfaces_qrf      Number of faces in the quadrature rule
  !> @param[in]     nqp             Number of quadrature points on each face
  !> @param[in]     wqp             Quadrature weights
  !!
  subroutine swe_buoyancy_gradient_facet_code( nlayers,                      &
                                               r_u_bd, geopot, buoyancy,     &
                                               stencil_size_w3,              &
                                               stencil_map_w3,               &
                                               ndf_w2, undf_w2,              &
                                               map_w2, w2_basis_face,        &
                                               ndf_w3, undf_w3,              &
                                               map_w3, w3_basis_face,        &
                                               nfaces_re_h,                  &
                                               nfaces_re,                    &
                                               normals_to_horizontal_faces,  &
                                               outward_normals_to_faces,     &
                                               adjacent_face,                &
                                               nfaces_qrf, nqp, wqp         )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: nqp
    integer(kind=i_def), intent(in) :: ndf_w2, ndf_w3
    integer(kind=i_def), intent(in) :: undf_w2, undf_w3
    integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
    integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3

    integer(kind=i_def), intent(in) :: stencil_size_w3
    integer(kind=i_def), dimension(ndf_w3, stencil_size_w3), intent(in)  :: stencil_map_w3

    real(kind=r_def), dimension(3,ndf_w2,nqp,4), intent(in) :: w2_basis_face
    real(kind=r_def), dimension(1,ndf_w3,nqp,4), intent(in) :: w3_basis_face

    integer(kind=i_def), intent(in) :: adjacent_face(:)
    real(kind=r_def),    intent(in) :: outward_normals_to_faces(:,:)
    real(kind=r_def),    intent(in) :: normals_to_horizontal_faces(:,:)
    integer(kind=i_def), intent(in) :: nfaces_re_h
    integer(kind=i_def), intent(in) :: nfaces_re
    integer(kind=i_def), intent(in) :: nfaces_qrf

    real(kind=r_def), dimension(undf_w2), intent(inout) :: r_u_bd
    real(kind=r_def), dimension(undf_w3), intent(in)    :: geopot
    real(kind=r_def), dimension(undf_w3), intent(in)    :: buoyancy

    real(kind=r_def), dimension(nqp,4), intent(in)      ::  wqp

    ! Internal variables
    integer(kind=i_def) :: df, face, face_next
    integer(kind=i_def) :: qp

    real(kind=r_def), dimension(ndf_w2) :: ru_bd_e
    real(kind=r_def), dimension(ndf_w3) :: geopot_e
    real(kind=r_def), dimension(ndf_w3) :: buoyancy_e, buoyancy_next_e

    real(kind=r_def) :: v(3)
    real(kind=r_def) :: buoyancy_av
    real(kind=r_def) :: geopot_at_fquad, bdary_term

    ru_bd_e(1:ndf_w2) = 0.0_r_def

    do face = 1, size( adjacent_face, 1 )
      ! Storing opposite face number on neighbouring cell
      face_next = adjacent_face(face)

      ! Computing geopotential in local cell
      do df = 1, ndf_w3
        geopot_e(df) = geopot( map_w3(df) )
      end do

      ! Computing buoyancy in local and adjacent cell
      do df = 1, ndf_w3
        buoyancy_e(df)      = buoyancy( stencil_map_w3(df, 1) )
        buoyancy_next_e(df) = buoyancy( stencil_map_w3(df, face+1) )
      end do

      ! Compute the boundary RHS integrated over one horizontal face
      do qp = 1, nqp
        buoyancy_av = 0.0_r_def
        geopot_at_fquad = 0.0_r_def
        do df = 1, ndf_w3
          buoyancy_av     = buoyancy_av &
                            + 0.5_r_def*(  buoyancy_e(df)*w3_basis_face(1,df,qp,face) &
                                         + buoyancy_next_e(df)*w3_basis_face(1,df,qp,face_next) )
          geopot_at_fquad = geopot_at_fquad + geopot_e(df)*w3_basis_face(1,df,qp,face)
        end do

        do df = 1, ndf_w2
          v  = w2_basis_face(:,df,qp,face)

          bdary_term = 0.5_r_def * dot_product(v, outward_normals_to_faces(:, face)) &
                                 * geopot_at_fquad * buoyancy_av
          ru_bd_e(df) = ru_bd_e(df) + wqp(qp,face) * bdary_term
        end do

      end do ! qp
    end do ! faces

    do df = 1, ndf_w2
      r_u_bd( map_w2(df) ) =  r_u_bd( map_w2(df) ) + ru_bd_e(df)
    end do

  end subroutine swe_buoyancy_gradient_facet_code

end module swe_buoyancy_gradient_facet_kernel_mod
