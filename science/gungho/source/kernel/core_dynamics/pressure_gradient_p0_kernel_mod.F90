!-----------------------------------------------------------------------------
! (c) Crown copyright 2024 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Calculates the pressure and geopotential gradients for oreder p=0
!!          elements
!> @details Computation of the dynamics forcing terms for the momentum equation
!!          cp*theta_v*grad(exner) + grad(geopotential) where lowest order
!!          elements are being used, i.e. exner, geopotential are a DG0 space
!!          and rhs is a RT0 space

module pressure_gradient_p0_kernel_mod

  use argument_mod,             only : arg_type, func_type,         &
                                       mesh_data_type,              &
                                       GH_FIELD, GH_READ, GH_INC,   &
                                       GH_SCALAR,                   &
                                       GH_REAL, STENCIL, CROSS2D,   &
                                       CELL_COLUMN
  use constants_mod,            only : r_def, i_def
  use fs_continuity_mod,        only : W2, W3, Wtheta
  use kernel_mod,               only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: pressure_gradient_p0_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                     &
         arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),                        &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W3),                        &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W3),                        &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, Wtheta, STENCIL(CROSS2D)),  &
         arg_type(GH_SCALAR,  GH_REAL, GH_READ)                             &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: pressure_gradient_p0_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: pressure_gradient_p0_code

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !> @brief Calculates the pressure and geopotential gradients for oreder p=0
  !!        elements.
  !> @param[in]     nlayers      Number of layers
  !> @param[in,out] rhs          Right-hand side of the momentum equation
  !> @param[in]     exner        Exner pressure
  !> @param[in]     geopotential Geopotential field
  !> @param[in]     theta_v      Virtual potential temperature
  !> @param[in]     smap_wt_size Size of the Wtheta stencil map
  !> @param[in]     smap_wt_max  Max size of the Wtheta stencil map
  !> @param[in]     smap_wt      Wtheta stencil map
  !> @param[in]     cp           Specific heat of dry air at constant pressure
  !> @param[in]     ndf_w2       Number of degrees of freedom per cell for W2
  !> @param[in]     undf_w2      Number of unique degrees of freedom for W2
  !> @param[in]     map_w2       Dofmap for the cell at the base of the column for W2
  !> @param[in]     ndf_w3       Number of degrees of freedom per cell for W3
  !> @param[in]     undf_w3      Number of unique degrees freedom for W3
  !> @param[in]     map_w3       Dofmap for the cell at the base of the column for W3
  !> @param[in]     ndf_wt       Number of degrees of freedom per cell for Wtheta
  !> @param[in]     undf_wt      Number of unique degrees of freedom for Wtheta
  !> @param[in]     map_wt       Dofmap for the cell at the base of the column for
  !!                             Wtheta
  subroutine pressure_gradient_p0_code( nlayers,                  &
                                        rhs, exner, geopotential, &
                                        theta_v,                  &
                                        smap_wt_size,             &
                                        smap_wt_max,              &
                                        smap_wt,                  &
                                        cp,                       &
                                        ndf_w2, undf_w2, map_w2,  &
                                        ndf_w3, undf_w3, map_w3,  &
                                        ndf_wt, undf_wt, map_wt )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w2, ndf_w3, ndf_wt
    integer(kind=i_def), intent(in) :: undf_w2, undf_w3, undf_wt

    integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
    integer(kind=i_def), dimension(ndf_w3), intent(in) :: map_w3
    integer(kind=i_def), dimension(ndf_wt), intent(in) :: map_wt

    integer(kind=i_def),                                    intent(in) :: smap_wt_max
    integer(kind=i_def), dimension(4),                      intent(in) :: smap_wt_size
    integer(kind=i_def), dimension(ndf_wt, smap_wt_max, 4), intent(in) :: smap_wt

    real(kind=r_def), dimension(undf_w2), intent(inout) :: rhs
    real(kind=r_def), dimension(undf_w3), intent(in)    :: exner
    real(kind=r_def), dimension(undf_w3), intent(in)    :: geopotential
    real(kind=r_def), dimension(undf_wt), intent(in)    :: theta_v
    real(kind=r_def),                     intent(in)    :: cp

    ! Internal variables
    integer(kind=i_def) :: df, dfv, k

    real(kind=r_def)               :: theta_av, theta_v_here, theta_v_next
    real(kind=r_def), dimension(6) :: div_v

    div_v = (/ -1.0_r_def, 1.0_r_def, 1.0_r_def, -1.0_r_def, -1.0_r_def, 1.0_r_def /)

    ! Horizontal gradients
    do df = 1, 4
      ! Check if there is a neighbouring cell, if not do not compute any
      ! horizontal terms
      if ( smap_wt_size(df) > 1 ) then
        do k = 0, nlayers-1

          theta_v_here = 0.5_r_def*(theta_v(map_wt(1)+k)       + theta_v(map_wt(2)+k))
          theta_v_next = 0.5_r_def*(theta_v(smap_wt(1,2,df)+k) + theta_v(smap_wt(2,2,df)+k))

          theta_av = 0.5_r_def * (theta_v_here + theta_v_next)

          rhs(map_w2(df)+k) = rhs(map_w2(df)+k) + cp*theta_av*exner(map_w3(1)+k)*div_v(df) &
                                                + geopotential(map_w3(1)+k)*div_v(df)
        end do
      end if
    end do

    ! Vertical gradients
    do df = 1,2
      dfv = df + 4
      do k = 0, nlayers-1
        rhs(map_w2(dfv)+k) = rhs(map_w2(dfv)+k) + cp*theta_v(map_wt(df)+k)*exner(map_w3(1)+k)*div_v(dfv) &
                                                + geopotential(map_w3(1)+k)*div_v(dfv)
      end do
    end do


  end subroutine pressure_gradient_p0_code

end module pressure_gradient_p0_kernel_mod
