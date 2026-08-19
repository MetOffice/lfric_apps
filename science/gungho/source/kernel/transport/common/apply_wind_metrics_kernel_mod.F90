!-------------------------------------------------------------------------------
! (C) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------

!-----------------------------------------------------------------------------!
!> @brief Kernel for applying metric terms to transported wind components.
!> @details Applies geometric and slope-related corrections to winds on W2,
!!          including spherical rotation terms where required.
!-----------------------------------------------------------------------------!
module apply_wind_metrics_kernel_mod

  use argument_mod,           only: arg_type,            &
                                    func_type, GH_BASIS, &
                                    GH_DIFF_BASIS,       &
                                    GH_FIELD, GH_REAL,   &
                                    GH_READ, GH_INC,     &
                                    CELL_COLUMN,         &
                                    STENCIL, CROSS2D,    &
                                    GH_EVALUATOR
  use constants_mod,          only: r_def, i_def, l_def
  use fs_continuity_mod,      only: W2, Wchi
  use kernel_mod,             only: kernel_type
  use base_mesh_config_mod,   only: geometry,           &
                                    geometry_spherical, &
                                    topology,           &
                                    topology_fully_periodic
  !use formulation_config_mod, only: shallow

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  type, public, extends(kernel_type) :: apply_wind_metrics_kernel_type
    private
    type(arg_type) :: meta_args(6) = (/                                     &
         arg_type(GH_FIELD,   GH_REAL, GH_INC,  W2),                        &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W2, STENCIL(CROSS2D)),      &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W2),                        &
         arg_type(GH_FIELD*3, GH_REAL, GH_READ, Wchi),                      &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W2),                        &
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W2)                         &
         /)
    type(func_type) :: meta_funcs(1) = (/         &
         func_type(Wchi, GH_BASIS, GH_DIFF_BASIS) &
         /)
    integer :: operates_on = CELL_COLUMN
    integer :: gh_shape = GH_EVALUATOR
  contains
    procedure, nopass :: apply_wind_metrics_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: apply_wind_metrics_code

contains

!> @brief Applies wind metrics transformation to convert wind components
!! @param[in]     nlayers       Number of layers
!! @param[in,out] metrics_uvw   Output wind with metrics applied in W2
!! @param[in]     uvw           Input vector wind component in W2
!! @param[in]     stencil_sizes Sizes of the cross2d stencil in each direction
!! @param[in]     stencil_max   Maximum stencil size
!! @param[in]     stencil_map   Stencil dofmap for W2
!! @param[in]     da            Area element in W2
!! @param[in]     chi1          First coordinate field in Wchi
!! @param[in]     chi2          Second coordinate field in Wchi
!! @param[in]     chi3          Third coordinate field in Wchi
!! @param[in]     chi1_d        First departure coordinate in W2
!! @param[in]     chi2_d        Second departure coordinate in W2
!! @param[in]     ndf_w2        Number of degrees of freedom per cell for W2
!! @param[in]     undf_w2       Number of unique degrees of freedom for W2
!! @param[in]     map_w2        Dofmap for the cell at the base of the column for W2
!! @param[in]     ndf_wx        Number of degrees of freedom per cell for Wchi
!! @param[in]     undf_wx       Number of unique degrees of freedom for Wchi
!! @param[in]     map_wx        Dofmap for the cell at the base of the column for Wchi
!! @param[in]     basis_wx      Basis functions for Wchi evaluated at W2 nodes
!! @param[in]     diff_basis_wx Differential basis functions for Wchi evaluated
!!                              at W2 nodes
subroutine apply_wind_metrics_code(nlayers,                 &
                                   metrics_uvw, uvw,        &
                                   stencil_sizes,           &
                                   stencil_max,             &
                                   stencil_map,             &
                                   da,                      &
                                   chi1, chi2, chi3,        &
                                   chi1_d, chi2_d,          &
                                   ndf_w2, undf_w2, map_w2, &
                                   ndf_wx, undf_wx, map_wx, &
                                   basis_wx, diff_basis_wx  &
                                  )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2
  integer(kind=i_def), intent(in) :: ndf_wx, undf_wx
  integer(kind=i_def), intent(in) :: stencil_max
  integer(kind=i_def), intent(in) :: stencil_sizes(4)

  integer(kind=i_def), dimension(ndf_w2), intent(in) :: map_w2
  integer(kind=i_def), dimension(ndf_wx), intent(in) :: map_wx

  integer(kind=i_def), intent(in) :: stencil_map(ndf_w2, stencil_max, 4)

  real(kind=r_def), dimension(1,ndf_wx,ndf_w2), intent(in) :: basis_wx
  real(kind=r_def), dimension(3,ndf_wx,ndf_w2), intent(in) :: diff_basis_wx

  real(kind=r_def), dimension(undf_w2), intent(in)    :: uvw
  real(kind=r_def), dimension(undf_w2), intent(inout) :: metrics_uvw
  real(kind=r_def), dimension(undf_w2), intent(in)    :: da
  real(kind=r_def), dimension(undf_wx), intent(in)    :: chi1, chi2, chi3
  real(kind=r_def), dimension(undf_w2), intent(in)    :: chi1_d
  real(kind=r_def), dimension(undf_w2), intent(in)    :: chi2_d

  ! Internal variables
  integer(kind=i_def) :: k, df, nm1, dfx, id, kmin
  integer(kind=i_def) :: d(4), dir(6), av_ids(3,4,6)

  real(kind=r_def) :: u_av, v_av
  real(kind=r_def) :: rot_mat(3,3), chi_uvw(2), chi_d_uvw(2), deltaAD, denom
  real(kind=r_def) :: uvw_lhs(3), uvw_rhs(3)
  real(kind=r_def) :: dx_z, dy_z, dz_z

  ! The rotation metric contains code for shallow atmosphere geometry,
  ! however the existing shallow switch only controls shallow
  ! (constant with height) gravity and not shallow geometry and so this
  ! option is disabled here until it is correctly implemented
  logical(kind=l_def), parameter :: shallow = .false.

  nm1 = nlayers-1

  ! Compute metric term based upon geometry (effectively to turn the delta chi1 &
  ! delta chi2 terms into m)
  if ( geometry == geometry_spherical ) then

    ! Compute averaging ids of u,v,w at u,v,w points
    do df = 1, 4
      d(df) = min(2, stencil_sizes(df) )
    end do

    ! Average u, v, w to West (u) point
    av_ids(1,:,1) = stencil_map(1,1,1)
    av_ids(2,:,1) = (/ stencil_map(2,1,1), stencil_map(4,1,1), stencil_map(2,d(1),1), stencil_map(4,d(1),1) /)
    av_ids(3,:,1) = (/ stencil_map(5,1,1), stencil_map(5,d(1),1), stencil_map(5,1,1)+1, stencil_map(5,d(1),1)+1 /)

    ! Average u, v, w to South (v) point
    av_ids(1,:,2) = (/ stencil_map(1,1,1), stencil_map(3,1,1), stencil_map(1,d(2),2), stencil_map(3,d(2),2) /)
    av_ids(2,:,2) = stencil_map(2,1,1)
    av_ids(3,:,2) = (/ stencil_map(5,1,1), stencil_map(5,d(2),2), stencil_map(5,1,1)+1, stencil_map(5,d(2),2)+1 /)

    ! Average u, v, w to East (u) point
    av_ids(1,:,3) = stencil_map(3,1,1)
    av_ids(2,:,3) = (/ stencil_map(2,1,1), stencil_map(4,1,1), stencil_map(2,d(3),3), stencil_map(4,d(3),3) /)
    av_ids(3,:,3) = (/ stencil_map(5,1,1), stencil_map(5,d(3),3), stencil_map(5,1,1)+1, stencil_map(5,d(3),3)+1 /)

    ! Average u, v, w to North (v) point
    av_ids(1,:,4) = (/ stencil_map(1,1,1), stencil_map(3,1,1), stencil_map(1,d(4),4), stencil_map(3,d(4),4) /)
    av_ids(2,:,4) = stencil_map(4,1,1)
    av_ids(3,:,4) = (/ stencil_map(5,1,1), stencil_map(5,d(4),4), stencil_map(5,1,1)+1, stencil_map(5,d(4),4)+1 /)

    ! Average u, v, w to Bottom (w) point
    av_ids(1,:,5) = (/ stencil_map(1,1,1), stencil_map(3,1,1), stencil_map(1,1,1)-1, stencil_map(3,1,1)-1 /)
    av_ids(2,:,5) = (/ stencil_map(2,1,1), stencil_map(4,1,1), stencil_map(2,1,1)-1, stencil_map(4,1,1)-1 /)
    av_ids(3,:,5) = stencil_map(5,1,1)

    ! Average u, v, w to Top (w) point
    av_ids(1,:,6) = (/ stencil_map(1,1,1), stencil_map(3,1,1), stencil_map(1,1,1)+1, stencil_map(3,1,1)+1 /)
    av_ids(2,:,6) = (/ stencil_map(2,1,1), stencil_map(4,1,1), stencil_map(2,1,1)+1, stencil_map(4,1,1)+1 /)
    av_ids(3,:,6) = stencil_map(6,1,1)

    ! directions of w2 dofs
    dir(:) = (/ 1, 2, 1, 2, 3, 3 /)

    if ( topology == topology_fully_periodic ) then
      ! Global model -> alpha-beta-h coordinates
    else
      ! Limited area model -> lon-lat-h coordinates
      ! Compute deep atmosphere rotation matrix

      ! Loop over all dofs ignoring the top w dof as this will be done by the layer above and the
      ! top of the model dof will be fixed by the boundary conditions
      do df = 1, ndf_w2 - 1
        ! Need rotation matrix on u, v & w points and hence need chi1, chi2, chi1_d, chi2_d w2 points
        chi_uvw = 0.0_r_def
        do dfx = 1, ndf_wx
          chi_uvw(1) = chi_uvw(1) + chi1(map_wx(dfx))*basis_wx(1, dfx, df)
          chi_uvw(2) = chi_uvw(2) + chi2(map_wx(dfx))*basis_wx(1, dfx, df)
        end do

        if ( df == 5 ) then
          kmin = 1
        else
          kmin = 0
        end if

        do k = kmin, nlayers-1
          chi_d_uvw(1) = chi1_d(map_w2(df)+k)
          chi_d_uvw(2) = chi2_d(map_w2(df)+k)

          deltaAD = chi_uvw(1) - chi_d_uvw(1)

          if ( shallow ) then
            ! Shallow atmosphere
            denom = 1.0_r_def/(1.0_r_def + sin(chi_uvw(2))*sin(chi_d_uvw(2)) + cos(chi_uvw(2))*cos(chi_d_uvw(2))*cos(deltaAD))
            rot_mat(1,1) = (cos(chi_uvw(2))*cos(chi_d_uvw(2)) + (1.0_r_def + sin(chi_uvw(2))*sin(chi_d_uvw(2)))*cos(deltaAD))*denom
            rot_mat(1,2) = (sin(chi_uvw(2))+sin(chi_d_uvw(2)))*sin(deltaAD)*denom
            rot_mat(1,3) = 0.0_r_def
            rot_mat(2,1) = -(sin(chi_uvw(2))+sin(chi_d_uvw(2)))*sin(deltaAD)*denom
            rot_mat(2,2) = rot_mat(1,1)
            rot_mat(2,3) = 0.0_r_def
            rot_mat(3,1) = 0.0_r_def
            rot_mat(3,2) = 0.0_r_def
            rot_mat(3,3) = 1.0_r_def
          else
            ! Deep atmosphere
            rot_mat(1,1) = cos(deltaAD)
            rot_mat(1,2) = sin(chi_d_uvw(2))*sin(deltaAD)
            rot_mat(1,3) = -cos(chi_d_uvw(2))*sin(deltaAD)
            rot_mat(2,1) = -sin(chi_uvw(2))*sin(deltaAD)
            rot_mat(2,2) = cos(chi_uvw(2))*cos(chi_d_uvw(2))+sin(chi_uvw(2))*sin(chi_d_uvw(2))*cos(deltaAD)
            rot_mat(2,3) = cos(chi_uvw(2))*sin(chi_d_uvw(2))-sin(chi_uvw(2))*cos(chi_d_uvw(2))*cos(deltaAD)
            rot_mat(3,1) = cos(chi_uvw(2))*sin(deltaAD)
            rot_mat(3,2) = sin(chi_uvw(2))*cos(chi_d_uvw(2))-cos(chi_uvw(2))*sin(chi_d_uvw(2))*cos(deltaAD)
            rot_mat(3,3) = sin(chi_uvw(2))*sin(chi_d_uvw(2))+cos(chi_uvw(2))*cos(chi_d_uvw(2))*cos(deltaAD)
          end if

          ! Compute uvw at this uvw point
          uvw_rhs = 0.0_r_def
          do id = 1,4
            uvw_rhs(1) = uvw_rhs(1) + 0.25_r_def*uvw(av_ids(1,id,df)+k)
            uvw_rhs(2) = uvw_rhs(2) - 0.25_r_def*uvw(av_ids(2,id,df)+k)
            uvw_rhs(3) = uvw_rhs(3) + 0.25_r_def*uvw(av_ids(3,id,df)+k)
          end do
          uvw_lhs = matmul(rot_mat,uvw_rhs)
          uvw_lhs(2) = -uvw_lhs(2)
          metrics_uvw(map_w2(df)+k) = uvw_lhs(dir(df))*da(map_w2(df)+k)
        end do ! nlayers
      end do ! ndf_w2
      metrics_uvw(map_w2(5)) = 0.0_r_def
      metrics_uvw(map_w2(6)+nlayers-1) = 0.0_r_def
    end if ! topology
  else
    do df = 1, 5
      metrics_uvw(map_w2(df):map_w2(df)+nm1) = uvw(map_w2(df):map_w2(df)+nm1) &
                                               *da(map_w2(df):map_w2(df)+nm1)
    end do
  end if

  ! Correction to vertical wind based upon model slope
  metrics_uvw(map_w2(5)) = 0.0_r_def
  do k = 1, nlayers - 1
    dx_z = 0.0_r_def ! delta_x(z) on W point
    dy_z = 0.0_r_def ! delta_y(z) on W point
    dz_z = 0.0_r_def ! delta_z(z) on W point
    do df = 1, ndf_wx
      dx_z = dx_z + chi3(map_wx(df)+k)*diff_basis_wx(1, df, 1)
      dy_z = dy_z + chi3(map_wx(df)+k)*diff_basis_wx(2, df, 1)
      dz_z = dz_z + chi3(map_wx(df)+k)*diff_basis_wx(3, df, 1)
    end do

    u_av = 0.5_r_def * ( metrics_uvw(map_w2(1)+k) + metrics_uvw(map_w2(3)+k))
    v_av = 0.5_r_def * ( metrics_uvw(map_w2(2)+k) + metrics_uvw(map_w2(4)+k) )
    metrics_uvw(map_w2(5)+k) = metrics_uvw(map_w2(5)+k) - u_av*dx_z/dz_z + v_av*dy_z/dz_z
  end do
  metrics_uvw(map_w2(6)+nlayers-1) = 0.0_r_def

end subroutine apply_wind_metrics_code

end module apply_wind_metrics_kernel_mod
