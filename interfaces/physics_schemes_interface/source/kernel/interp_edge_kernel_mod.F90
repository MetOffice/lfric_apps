!-----------------------------------------------------------------------------
! (C) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------

!> @brief Interpolates scalar variables from W3 (Wtheta) to W2 (shifted W2) dofs
!> @details Takes all the variables required for BL momentum mixing and
!>          interpolates them from their lowest order W3 dof to W2 dofs
!>          so that wind increments can be calculated in their native space

module interp_edge_kernel_mod

  use kernel_mod,               only: kernel_type
  use argument_mod,             only: arg_type, GH_SCALAR, GH_FIELD, &
                                      GH_REAL, GH_INTEGER, GH_LOGICAL, &
                                      GH_INC, GH_READ, &
                                      ANY_SPACE_1, CELL_COLUMN, &
                                      ANY_DISCONTINUOUS_SPACE_1

  use constants_mod,            only: r_def, i_def, l_def
  use fs_continuity_mod,        only: W2
  use kernel_mod,               only: kernel_type

  implicit none

  private

  !----------------------------------------------------------------------------
  ! Public types
  !----------------------------------------------------------------------------
  !> Kernel metadata type.
  type, public, extends(kernel_type) :: interp_edge_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                         &
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                              &! flag_surface
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                              &! n_interp
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  ANY_DISCONTINUOUS_SPACE_1),  &! data_center
         arg_type(GH_FIELD,  GH_REAL,    GH_INC,   ANY_SPACE_1),                &! data_edge
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,  W2)                          &! w2_rmultiplicity
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: interp_edge_code
  end type interp_edge_kernel_type

  !----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !----------------------------------------------------------------------------
  public :: interp_edge_code

contains

  !> @brief Subroutine to do the re-mapping
  !> @param[in]     nlayers          Number of layers
  !> @param[in]     flag_surface     Flag of interpolating surface variables
  !> @param[in]     n_interp         Number of indices to interpolate
  !> @param[in]     data_center      Input data defined in W3(Wtheta) space
  !> @param[in,out] data_edge        Output data defined in W2(shifted W2) space
  !> @param[in]     w2_rmultiplicity Reciprocal of multiplicity for w2
  !> @param[in]     ndf_in           Number of DOFs for W3(Wtheta) space
  !> @param[in]     undf_in          Number of unique DOFs for W3(Wtheta) space
  !> @param[in]     map_in           Dofmap for W3(Wtheta) space
  !> @param[in]     ndf_out          Number of DOFs for W2(shifted W2) space
  !> @param[in]     undf_out         Number of unique DOFs for W2(shifted W2) space
  !> @param[in]     map_out          Dofmap for W2(shifted W2) space
  !> @param[in]     ndf_w2           Number of DOFs for W2 surface space
  !> @param[in]     undf_w2          Number of unique DOFs for W2 surface space
  !> @param[in]     map_w2           Dofmap for W2 surface space
  subroutine interp_edge_code(nlayers,          &
                              flag_surface,     &
                              n_interp,         &
                              data_center,      &
                              data_edge,        &
                              w2_rmultiplicity, &
                              ndf_in,           &
                              undf_in,          &
                              map_in,           &
                              ndf_out,          &
                              undf_out,         &
                              map_out,          &
                              ndf_w2,           &
                              undf_w2,          &
                              map_w2)

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in) :: nlayers
    logical(kind=l_def), intent(in) :: flag_surface
    integer(kind=i_def), intent(in) :: n_interp

    integer(kind=i_def), intent(in) :: ndf_in, ndf_out, ndf_w2
    integer(kind=i_def), intent(in) :: undf_in, undf_out, undf_w2
    integer(kind=i_def), intent(in) :: map_in(ndf_in)
    integer(kind=i_def), intent(in) :: map_out(ndf_out)
    integer(kind=i_def), intent(in) :: map_w2(ndf_w2)

    real(kind=r_def), dimension(undf_in), intent(in) :: data_center
    real(kind=r_def), dimension(undf_out), intent(inout) :: data_edge
    real(kind=r_def), dimension(undf_w2), intent(in) :: w2_rmultiplicity

    ! Internal variables
    integer :: k, df

    ! Map w3 (wtheta) variables into w2 (shifted w2) space
    if (flag_surface) then
      do df = 1,4
        do k = 0, n_interp-1
          data_edge(map_out(df) + k) = data_edge(map_out(df) + k) +            &
                                      w2_rmultiplicity(map_w2(df)) *           &
                                      data_center(map_in(1) + k)
        end do
      end do
    else
      do df = 1,4
        do k = 0, n_interp-1
          data_edge(map_out(df) + k) = data_edge(map_out(df) + k) +            &
                                      w2_rmultiplicity(map_w2(df) + k) *       &
                                      data_center(map_in(1) + k)
        end do
      end do
    end if

  end subroutine interp_edge_code

end module interp_edge_kernel_mod
