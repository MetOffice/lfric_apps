!-------------------------------------------------------------------------------
! (c) Crown copyright 2026 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Compute advective increment u.grad(z) using the Nivana scheme used in
!!        the ffsl and sl transport schemes.
!> @details Compute the advective increment u.grad(z) = u*dz/dx + v*dz/dy where
!! the horizontal gradients are computed using the same reconstruction method as
!! the transport scheme (in this case Nivana).
!!
!> @note This kernel only works when field is a W3/Wtheta field at lowest order.

module horizontal_cubic_sl_metric_kernel_mod

  use argument_mod,          only: arg_type,                     &
                                   GH_FIELD, GH_REAL,            &
                                   CELL_COLUMN, GH_WRITE,        &
                                   GH_READ, GH_SCALAR,           &
                                   STENCIL, CROSS2D, GH_INTEGER, &
                                   ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,         only: r_def, i_def, l_def
  use fs_continuity_mod,     only: W2H
  use kernel_mod,            only: kernel_type

  implicit none

  private

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the PSy layer
  type, public, extends(kernel_type) :: horizontal_cubic_sl_metric_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                        &
        arg_type(GH_FIELD,  GH_REAL,    GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),  & ! increment
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1,   &
                                                            STENCIL(CROSS2D)), & ! z
        arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2H)                         & ! wind
    /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: horizontal_cubic_sl_metric_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public :: horizontal_cubic_sl_metric_code

contains

  !> @brief Compute advective increment u.grad(z) using the Nivana scheme used in
  !!        the ffsl and sl transport schemes.
  !> @param[in]     nlayers           Number of layers
  !> @param[in,out] increment         Horizontal metric increment
  !> @param[in]     z                 Height field values used in reconstruction
  !> @param[in]     stencil_sizes     Sizes of the branches of the cross stencil
  !> @param[in]     stencil_max       Maximum size of a cross stencil branch
  !> @param[in]     stencil_map       Dofmap for the z stencil
  !> @param[in]     wind              Horizontal wind/departure field on W2H
  !> @param[in]     ndf_wf            Num of DoFs for field per cell
  !> @param[in]     undf_wf           Num of DoFs for this partition for field
  !> @param[in]     map_wf            Map for Wf
  !> @param[in]     ndf_w2h           Num of DoFs for W2H per cell
  !> @param[in]     undf_w2h          Num of DoFs for this partition for W2H
  !> @param[in]     map_w2h           Map for W2H
  subroutine horizontal_cubic_sl_metric_code( nlayers,         &
                                              increment,       &
                                              z,               &
                                              stencil_sizes,   &
                                              stencil_max,     &
                                              stencil_map,     &
                                              wind,            &
                                              ndf_wf,          &
                                              undf_wf,         &
                                              map_wf,          &
                                              ndf_w2h,         &
                                              undf_w2h,        &
                                              map_w2h )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: undf_wf
    integer(kind=i_def), intent(in) :: ndf_wf
    integer(kind=i_def), intent(in) :: undf_w2h
    integer(kind=i_def), intent(in) :: ndf_w2h
    integer(kind=i_def), intent(in) :: stencil_max
    integer(kind=i_def), intent(in) :: stencil_sizes(4)

    ! Arguments: Maps
    integer(kind=i_def), intent(in) :: map_wf(ndf_wf)
    integer(kind=i_def), intent(in) :: map_w2h(ndf_w2h)
    integer(kind=i_def), intent(in) :: stencil_map(ndf_wf,stencil_max,4)

    ! Arguments: Fields
    real(kind=r_def),   intent(inout) :: increment(undf_wf)
    real(kind=r_def),   intent(in)    :: z(undf_wf)
    real(kind=r_def),   intent(in)    :: wind(undf_w2h)

    ! Local scalars
    integer(kind=i_def) :: k, kp, km, d, d2(4), d3(4)

    real(kind=r_def)   :: dzdx, dzdy
    real(kind=r_def)   :: z_l, z_r, up, um, vp, vm

    ! Interpolation coefficients
    real(kind=r_def), parameter :: b0 = -1.0_r_def/6.0_r_def
    real(kind=r_def), parameter :: b1 = 5.0_r_def/6.0_r_def
    real(kind=r_def), parameter :: b2 = 2.0_r_def/6.0_r_def
    real(kind=r_def), parameter :: c1 = 2.0_r_def/6.0_r_def
    real(kind=r_def), parameter :: c2 = 5.0_r_def/6.0_r_def
    real(kind=r_def), parameter :: c3 = -1.0_r_def/6.0_r_def

    ! Ensure that we don't do out of the domain and if there are not enough points
    ! then revert to constant reconstruction
    d2(:) = 1
    d3(:) = 1
    do d = 1, 4
      if ( stencil_sizes(d) == stencil_max ) then
        d2(d) = 2
        d3(d) = 3
      end if
    end do

    do k = 0, nlayers
      km = max(0, k-1)
      kp = min(nlayers-1, k)

      um =  0.5_r_def*(wind(map_w2h(1)+km) + wind(map_w2h(1)+kp))
      up =  0.5_r_def*(wind(map_w2h(3)+km) + wind(map_w2h(3)+kp))
      vm = -0.5_r_def*(wind(map_w2h(2)+km) + wind(map_w2h(2)+kp))
      vp = -0.5_r_def*(wind(map_w2h(4)+km) + wind(map_w2h(4)+kp))

    ! dzdx
    ! Compute upwind Z on the left and right sides of the cell
    if ( um > 0.0_r_def ) then
      z_l = b0*z(stencil_map(1,d3(1),1)+k) + b1*z(stencil_map(1,d2(1),1)+k) + b2*z(stencil_map(1,1,1)+k)
    else
      z_l = c1*z(stencil_map(1,d2(1),1)+k) + c2*z(stencil_map(1,1,1)+k) + c3*z(stencil_map(1,d2(3),3)+k)
    end if
    if ( up > 0.0_r_def ) then
      z_r = b0*z(stencil_map(1,d2(1),1)+k) + b1*z(stencil_map(1,1,1)+k) + b2*z(stencil_map(1,d2(3),3)+k)
    else
      z_r = c1*z(stencil_map(1,1,3)+k) + c2*z(stencil_map(1,d2(3),3)+k) + c3*z(stencil_map(1,d3(3),3)+k)
    end if

    dzdx = (z_r- z_l)

    ! dzdy
    ! Compute upwind Z on the left and right sides of the cell
    if ( vm > 0.0_r_def ) then
      z_l = b0*z(stencil_map(1,d3(2),2)+k) + b1*z(stencil_map(1,d2(2),2)+k) + b2*z(stencil_map(1,1,2)+k)
    else
      z_l = c1*z(stencil_map(1,d2(2),2)+k) + c2*z(stencil_map(1,1,2)+k) + c3*z(stencil_map(1,d2(4),4)+k)
    end if
    if ( vp > 0.0_r_def ) then
      z_r = b0*z(stencil_map(1,d2(2),2)+k) + b1*z(stencil_map(1,1,2)+k) + b2*z(stencil_map(1,d2(4),4)+k)
    else
      z_r = c1*z(stencil_map(1,1,4)+k) + c2*z(stencil_map(1,d2(4),4)+k) + c3*z(stencil_map(1,d3(4),4)+k)
    end if

    dzdy = (z_r- z_l)

    increment(map_wf(1)+k) = 0.5_r_def*(um+up)*dzdx + 0.5_r_def*(vm+vp)*dzdy
  end do

  end subroutine horizontal_cubic_sl_metric_code

end module horizontal_cubic_sl_metric_kernel_mod
