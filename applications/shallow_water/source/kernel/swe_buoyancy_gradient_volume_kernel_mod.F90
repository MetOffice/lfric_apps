!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Computes the volume integral for the rhs of the thermal shallow water
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
!!          The third term is evaluated as a facet integral and is handled in
!!          swe_buoyancy_gradient_facet_kernel_mod.
!!          The first two terms can be combined to give the volume integral
!!
!!          \f[ -\frac{1}{2} \int s \Phi \nabla \cdot v dV\f]
!!
!!          which is evaluated here.
!!          Note that this kernel is only suitable for lowest order function
!!          space, zero surface geopotential, and planar geometry.
!!
!> @todo Modify kernel to work for non-zero surface geopotential and
!!       and any geometry, see ticket #3002.

module swe_buoyancy_gradient_volume_kernel_mod

  use argument_mod,            only: arg_type, func_type,       &
                                     GH_FIELD, GH_READ, GH_INC, &
                                     GH_BASIS, GH_DIFF_BASIS,   &
                                     GH_QUADRATURE_XYoZ,        &
                                     GH_REAL, CELL_COLUMN
  use constants_mod,           only: r_def, i_def
  use fs_continuity_mod,       only: W2, W3
  use kernel_mod,              only: kernel_type

  implicit none

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: swe_buoyancy_gradient_volume_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/              &
        arg_type(GH_FIELD, GH_REAL, GH_INC,  W2),    &
        arg_type(GH_FIELD, GH_REAL, GH_READ, W3),    &
        arg_type(GH_FIELD, GH_REAL, GH_READ, W3)     &
        /)
    type(func_type) :: meta_funcs(2) = (/            &
        func_type(W2,     GH_DIFF_BASIS),            &
        func_type(W3,     GH_BASIS)                  &
        /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_QUADRATURE_XYoZ
  contains
    procedure, nopass ::swe_buoyancy_gradient_volume_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public swe_buoyancy_gradient_volume_code

contains

!> @brief Compute the volume integral of the buoyancy and geopotential gradients
!!        for the rhs of the thermal shallow water momentum equation.
!! @param[in]     nlayers       Number of layers
!! @param[in,out] r_u           Momentum equation right hand side
!! @param[in]     geopot        Geopotential field
!! @param[in]     buoyancy      Buoyancy field
!! @param[in]     ndf_w2        Number of degrees of freedom per cell for w2
!! @param[in]     undf_w2       Number unique of degrees of freedom  for w2
!! @param[in]     map_w2        Dofmap for the cell at the base of the column for w2
!! @param[in]     w2_diff_basis Differential of the basis functions evaluated at quadrature points
!! @param[in]     ndf_w3        Number of degrees of freedom per cell for w3
!! @param[in]     undf_w3       Number unique of degrees of freedom  for w3
!! @param[in]     map_w3        Dofmap for the cell at the base of the column for w3
!! @param[in]     w3_basis      Basis functions evaluated at gaussian quadrature points
!! @param[in]     nqp_h         Number of quadrature points in the horizontal
!! @param[in]     nqp_v         Number of quadrature points in the vertical
!! @param[in]     wqp_h         Horizontal quadrature weights
!! @param[in]     wqp_v         Vertical quadrature weights
subroutine swe_buoyancy_gradient_volume_code(nlayers, r_u, geopot, buoyancy,     &
                                         ndf_w2, undf_w2, map_w2, w2_diff_basis, &
                                         ndf_w3, undf_w3, map_w3, w3_basis,      &
                                         nqp_h, nqp_v, wqp_h, wqp_v )

  implicit none

  ! Arguments
  integer, intent(in) :: nlayers, nqp_h, nqp_v
  integer, intent(in) :: ndf_w2, ndf_w3
  integer, intent(in) :: undf_w2, undf_w3
  integer, dimension(ndf_w2),  intent(in)  :: map_w2, map_w3

  real(kind=r_def), dimension(1,ndf_w2,nqp_h,nqp_v),   intent(in) :: w2_diff_basis
  real(kind=r_def), dimension(1,ndf_w3,nqp_h,nqp_v),   intent(in) :: w3_basis

  real(kind=r_def), dimension(undf_w2),  intent(inout) :: r_u
  real(kind=r_def), dimension(undf_w3),  intent(in)    :: geopot
  real(kind=r_def), dimension(undf_w3),  intent(in)    :: buoyancy

  real(kind=r_def), dimension(nqp_h), intent(in)      :: wqp_h
  real(kind=r_def), dimension(nqp_v), intent(in)      :: wqp_v

  ! Internal variables
  integer(kind=i_def) :: df, qp1, qp2

  real(kind=r_def), dimension(ndf_w2) :: ru_e
  real(kind=r_def), dimension(ndf_w3) :: geopot_e
  real(kind=r_def), dimension(ndf_w3) :: b_e

  real(kind=r_def) :: geopot_at_quad, b_at_quad, geo_grad_term, dv

  ! Initialise rhs, geopotential and buoyancy at dofs
  ru_e(1:ndf_w2) = 0.0_r_def
  do df = 1, ndf_w3
    geopot_e(df) = geopot( map_w3(df) )
  end do
  do df = 1, ndf_w3
    b_e(df) = buoyancy( map_w3(df) )
  end do
  ! Compute the RHS integrated over one cell
  do qp2 = 1, nqp_v
    do qp1 = 1, nqp_h

      ! Geopotential and Buoyancy
      geopot_at_quad = 0.0_r_def
      b_at_quad = 0.0_r_def
      do df = 1, ndf_w3
        geopot_at_quad = geopot_at_quad + geopot_e(df)*w3_basis(1,df,qp1,qp2)
        b_at_quad      = b_at_quad + b_e(df)*w3_basis(1,df,qp1,qp2)
      end do

      ! Buoyancy
      b_at_quad = 0.0_r_def
      do df = 1, ndf_w3
        b_at_quad = b_at_quad + b_e(df)*w3_basis(1,df,qp1,qp2)
      end do

      ! Combine buoyancy, geopotential, and differential of W2 basis function
      ! Works for lowest order DG space only
      do df = 1, ndf_w2
        dv  = w2_diff_basis(1,df,qp1,qp2)

        geo_grad_term = 0.5_r_def*dv*geopot_at_quad*b_at_quad

        ru_e(df) = ru_e(df) +  wqp_h(qp1)*wqp_v(qp2)*geo_grad_term
      end do
    end do
  end do
  do df = 1, ndf_w2
    r_u( map_w2(df) ) =  r_u( map_w2(df) ) + ru_e(df)
  end do

end subroutine swe_buoyancy_gradient_volume_code

end module swe_buoyancy_gradient_volume_kernel_mod
