!-----------------------------------------------------------------------------
! (c) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Set the number of cloud layers to the highest cloud layer
module set_n_cloud_layer_kernel_mod

  use argument_mod,      only : arg_type,                  &
                                GH_FIELD, GH_SCALAR,       &
                                GH_REAL, GH_INTEGER,       &
                                GH_READ, GH_WRITE,         &
                                GH_READWRITE, CELL_COLUMN, &
                                ANY_DISCONTINUOUS_SPACE_1
  use fs_continuity_mod, only : Wtheta
  use constants_mod,     only : r_def, i_def
  use kernel_mod,        only : kernel_type

  implicit none

  private

  type, public, extends(kernel_type) :: set_n_cloud_layer_kernel_type
    private
    type(arg_type) :: meta_args(4) = (/ &
      arg_type(GH_FIELD, GH_REAL,     GH_READ, Wtheta ),                    & ! radiative_conv_fraction
      arg_type(GH_FIELD, GH_REAL,     GH_READ, Wtheta ),                    & ! radiative_cloud_fraction
      arg_type(GH_FIELD, GH_INTEGER,  GH_WRITE, ANY_DISCONTINUOUS_SPACE_1), & ! n_cloud_layers
      arg_type(GH_SCALAR, GH_REAL,    GH_READ )                             & ! cloud_fraction_min
      /)

    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: set_n_cloud_layer_code

  end type set_n_cloud_layer_kernel_type

  public :: set_n_cloud_layer_code

contains
  !> @param[in]    nlayers                  Number of layers in W3 2D field
  !> @param[in]    radiative_conv_fraction  Convective cloud fraction
  !> @param[in]    radiative_cloud_fraction Large scale cloud fraction
  !> @param[inout] n_cloud_layers           Number of cloud layers
  !> @param[in]    cloud_fraction_min       Minimum cloud fraction
  !> @param[in]    ndf_wtheta               No. DOFs per cell for Wtheta space
  !> @param[in]    undf_wtheta              No. unique DOFs for Wtheta space
  !> @param[in]    map_wtheta               Dofmap for Wtheta space column base cell
  !> @param[in]    ndf_n_cloud_layers       No. DOFs per cell for W3 2D space
  !> @param[in]    undf_n_cloud_layers      No. unique DOFs for W3 2D space
  !> @param[in]    map_n_cloud_layers       Dofmap for W3 space column base cell
  subroutine set_n_cloud_layer_code( nlayers, &
                                          radiative_conv_fraction, &
                                          radiative_cloud_fraction, &
                                          n_cloud_layers, &
                                          cloud_fraction_min, &
                                          ndf_wtheta, &
                                          undf_wtheta, &
                                          map_wtheta, &
                                          ndf_n_cloud_layers, &
                                          undf_n_cloud_layers, &
                                          map_n_cloud_layers )

    use constants_mod, only: r_def, i_def

    implicit none

    ! Dummy arguments
    integer(i_def), intent(in) :: nlayers
    integer(i_def), intent(in) :: ndf_n_cloud_layers
    integer(i_def), intent(in), dimension(ndf_n_cloud_layers) :: map_n_cloud_layers
    integer(i_def), intent(in) :: ndf_wtheta
    integer(i_def), intent(in), dimension(ndf_wtheta) :: map_wtheta
    integer(i_def), intent(in) :: undf_n_cloud_layers, undf_wtheta
    real(r_def), intent(in) :: cloud_fraction_min
    integer(i_def), intent(inout), dimension(undf_n_cloud_layers) :: n_cloud_layers
    real(r_def), intent(in), dimension(undf_wtheta) :: radiative_conv_fraction
    real(r_def), intent(in), dimension(undf_wtheta) :: radiative_cloud_fraction

    ! Internal constants and variables
    integer(i_def) :: k, j

    ! Set the number of cloud layers to the highest cloud layer
    n_cloud_layers(map_n_cloud_layers(1)) = 0_i_def
    do j=nlayers, 1, -1
      k = map_wtheta(1) + j
      if ( (radiative_cloud_fraction(k) + radiative_conv_fraction(k)) &
        > cloud_fraction_min ) then
        n_cloud_layers(map_n_cloud_layers(1)) = j
        exit
      end if
    end do


  end subroutine set_n_cloud_layer_code

end module set_n_cloud_layer_kernel_mod