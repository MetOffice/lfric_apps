!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Modify the vertical wind to include consistent computation of grid metric.
!> @details Modify the vertical wind profile so that:
!!          w => w + u*(dz/dx - dz/dx_A) + v*(dz/dy - dz/dy_A)
!!          where subscript A indicates computation of terms by the advection
!!          operator.

module consistent_wind_kernel_mod

use argument_mod,      only : arg_type, func_type,     &
                              GH_FIELD, GH_REAL,       &
                              GH_READWRITE, GH_READ,   &
                              GH_BASIS, GH_DIFF_BASIS, &
                              CELL_COLUMN, GH_EVALUATOR
use constants_mod,     only : r_def, i_def, r_tran
use fs_continuity_mod, only : Wtheta, W2, W2v, Wchi
use kernel_mod,        only : kernel_type

implicit none

private
!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: consistent_wind_kernel_type
  private
  type(arg_type) :: meta_args(4) = (/                     &
       arg_type(GH_FIELD, GH_REAL, GH_READWRITE, W2v),    &
       arg_type(GH_FIELD, GH_REAL, GH_READ,      W2),     &
       arg_type(GH_FIELD, GH_REAL, GH_READ,      Wtheta), &
       arg_type(GH_FIELD, GH_REAL, GH_READ,      Wchi)    &
       /)
  type(func_type) :: meta_funcs(2) = (/                   &
       func_type(W2,   GH_BASIS),                         &
       func_type(Wchi, GH_DIFF_BASIS)                     &
       /)
  integer :: operates_on = CELL_COLUMN
  integer :: gh_shape = GH_EVALUATOR
contains
  procedure, nopass :: consistent_wind_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: consistent_wind_code

contains

!> @brief Modify the vertical wind to include consistent computation of grid metric.
!> @param[in]     nlayers         Number of vertical layers
!> @param[in,out] consistent_wind Wind with the vertical component modified to
!!                                include consistent metrics
!> @param[in]     wind            Unmodified wind field
!> @param[in]     theta_metrics   Horizontal component of u.grad(z) computed by
!!                                advection scheme
!> @param[in]     height          Height above the surface
!> @param[in]     ndf_w2v         Number of degrees of freedom per cell for W2v
!> @param[in]     undf_w2v        Total number of degrees of freedom for W2v
!> @param[in]     map_w2v         Dofmap of the wind field
!> @param[in]     ndf_w2          Number of degrees of freedom per cell for W2
!> @param[in]     undf_w2         Total number of degrees of freedom for W2
!> @param[in]     map_w2          Dofmap of the wind field
!> @param[in]     basis_w2        Basis function of the wind space evaluated on
!!                                W2 nodal points
!> @param[in]     ndf_wt          Number of degrees of freedom per cell for Wtheta
!> @param[in]     undf_wt         Total number of degrees of freedom for Wtheta
!> @param[in]     map_wt          Dofmap of the metric field
!> @param[in]     ndf_wx          Number of degrees of freedom per cell for the coordinate space
!> @param[in]     undf_wx         Total number of degrees of freedom for the coordinate space
!> @param[in]     map_wx          Dofmap of the coordinate space
!> @param[in]     diff_basis_wx   Differential of basis function of the coordinate space
!!                                evaluated on W2v nodal points
subroutine consistent_wind_code(nlayers,                   &
                                consistent_wind,           &
                                wind,                      &
                                theta_metrics,             &
                                height,                    &
                                ndf_w2v,                   &
                                undf_w2v,                  &
                                map_w2v,                   &
                                ndf_w2,                    &
                                undf_w2,                   &
                                map_w2,                    &
                                basis_w2,                  &
                                ndf_wt,                    &
                                undf_wt,                   &
                                map_wt,                    &
                                ndf_wx,                    &
                                undf_wx,                   &
                                map_wx,                    &
                                diff_basis_wx )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_wt, undf_wt,   &
                                     ndf_w2, undf_w2,   &
                                     ndf_w2v, undf_w2v, &
                                     ndf_wx, undf_wx

  integer(kind=i_def), dimension(ndf_wt),  intent(in) :: map_wt
  integer(kind=i_def), dimension(ndf_w2),  intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_w2v), intent(in) :: map_w2v
  integer(kind=i_def), dimension(ndf_wx),  intent(in) :: map_wx

  real(kind=r_tran), dimension(undf_w2v), intent(inout) :: consistent_wind
  real(kind=r_tran), dimension(undf_wt),  intent(in)    :: theta_metrics
  real(kind=r_tran), dimension(undf_w2),  intent(in)    :: wind
  real(kind=r_def), dimension(undf_wx),   intent(in)    :: height

  real(kind=r_def), dimension(3,ndf_wx,ndf_w2v), intent(in) :: diff_basis_wx
  real(kind=r_def), dimension(3,ndf_w2,ndf_w2v), intent(in) :: basis_w2

  ! Local variables
  integer(kind=i_def) :: k, df, w2v_df

  real(kind=r_tran) :: dz, dzdx, dzdy
  real(kind=r_tran), dimension(3,0:nlayers) :: u_av

  ! Compute u and v averaged to w points
  u_av = 0.0_r_tran
  do k = 0, nlayers-1
    do df = 1,4
      ! Basis W2 has shape (dim, ndf_w2, ndf_w2v) = (3,6,2)
      u_av(:,k)   = u_av(:,k)   + wind(map_w2(df)+k)*real(basis_w2(:,df,1), r_tran)
      u_av(:,k+1) = u_av(:,k+1) + wind(map_w2(df)+k)*real(basis_w2(:,df,2), r_tran)
    end do
  end do
  u_av(:,1:nlayers-1) = 0.5_r_tran*u_av(:,1:nlayers-1)

  layer_loop: do k = 1, nlayers-1
    ! Compute dz/dx & dz/dy on w/theta points
    w2v_df = 1
    dz = 0.0_r_tran
    dzdx = 0.0_r_tran
    dzdy = 0.0_r_tran
    do df = 1,ndf_wx
      dzdx = dzdx + real(height(map_wx(df)+k)*diff_basis_wx(1,df,w2v_df), r_tran)
      dzdy = dzdy + real(height(map_wx(df)+k)*diff_basis_wx(2,df,w2v_df), r_tran)
      dz   = dz   + real(height(map_wx(df)+k)*diff_basis_wx(3,df,w2v_df), r_tran)
    end do

    consistent_wind(map_w2v(w2v_df)+k) = consistent_wind(map_w2v(w2v_df)+k) &
      + (u_av(1,k)*dzdx + u_av(2,k)*dzdy - theta_metrics(map_wt(1)+k))/dz
  end do layer_loop
  ! Can leave first & last w points unchanged since we want to enforce zero flux
  ! boundary conditions

end subroutine consistent_wind_code

end module consistent_wind_kernel_mod
