!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Kernel to compute the vertical cubic semi-Lagragian advection of a field
!!        in the vertical direction for the linear model.
!> @Details The 1D vertical advective transport equation for a W3/Wtheta variable
!!          is solved using a cubic semi-Lagragian advection scheme. There are two
!!          parts to the TL advection equation, and the update has the form
!!          f^{n+1} = f^{n} - ls_u dt grad f - u_pert dt grad ls_f
!!          The first two terms on the RHS are solved using the SL scheme, and the
!!          third term is computed using the gradient of the ls field in the
!!          departure cell.

module tl_vertical_cubic_sl_kernel_mod

  use argument_mod,          only : arg_type,                  &
                                    GH_FIELD, GH_SCALAR,       &
                                    GH_REAL, GH_INTEGER,       &
                                    GH_READWRITE, GH_READ,     &
                                    CELL_COLUMN, GH_LOGICAL,   &
                                    ANY_DISCONTINUOUS_SPACE_1, &
                                    ANY_DISCONTINUOUS_SPACE_2
  use fs_continuity_mod,     only : W2v
  use constants_mod,         only : r_tran, i_def, l_def, EPS_R_TRAN
  use kernel_mod,            only : kernel_type
  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed
  !>                                      by the PSy layer.
  type, public, extends(kernel_type) :: tl_vertical_cubic_sl_kernel_type
    private
    type(arg_type) :: meta_args(11) = (/                                           &
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! field
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1), & ! ls_field
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      W2v),                       & ! dep_dist_pert
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! cubic_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! cubic_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! cubic_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! cubic_coef
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! cubic_indices
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! cubic_indices
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! cubic_indices
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2)  & ! cubic_indices
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: tl_vertical_cubic_sl_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: tl_vertical_cubic_sl_code

  contains

  !-------------------------------------------------------------------------------
  !> @details This kernel interpolates the field to the departure point
  !!          using 1d-Cubic-Lagrange interpolation. It then adds the
  !!          contribution from the gradient of the ls field in the
  !!          departure cell.
  !> @param[in]     nlayers         The number of layers
  !> @param[in,out] field           The perturbation field to be advected
  !> @param[in]     ls_field        The ls field
  !> @param[in]     dep_dist_pert   The perturbation wind departure point
  !> @param[in]     cubic_coef      The cubic interpolation coefficients (1-4)
  !> @param[in]     cubic_indices   The cubic interpolation indices (1-4)
  !> @param[in]     ndf_wf          Num Dofs per cell for the field
  !> @param[in]     undf_wf         Num Dofs in this partition for the field
  !> @param[in]     map_wf          Dofmap for the field
  !> @param[in]     ndf_w2          Num Dofs per cell for dep_dist_pert
  !> @param[in]     undf_w2         Num Dofs in this partition for dep_dist_pert
  !> @param[in]     map_w2          Dofmap for dep_dist_pert
  !> @param[in]     ndf_wc          Num Dofs per cell for the coefficients
  !> @param[in]     undf_wc         Num Dofs per cell in this partition
  !!                                for the coefficients
  !> @param[in]     map_wc          Dofmap for the coefficients
  !-------------------------------------------------------------------------------
  subroutine tl_vertical_cubic_sl_code( nlayers,                 &
                                        field,                   &
                                        ls_field,                &
                                        dep_dist_pert,           &
                                        cubic_coef_1,            &
                                        cubic_coef_2,            &
                                        cubic_coef_3,            &
                                        cubic_coef_4,            &
                                        cubic_indices_1,         &
                                        cubic_indices_2,         &
                                        cubic_indices_3,         &
                                        cubic_indices_4,         &
                                        ndf_wf, undf_wf, map_wf, &
                                        ndf_w2, undf_w2, map_w2, &
                                        ndf_wc, undf_wc, map_wc )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ndf_wf
    integer(kind=i_def), intent(in)    :: undf_wf
    integer(kind=i_def), intent(in)    :: ndf_w2
    integer(kind=i_def), intent(in)    :: undf_w2
    integer(kind=i_def), intent(in)    :: ndf_wc
    integer(kind=i_def), intent(in)    :: undf_wc
    integer(kind=i_def), intent(in)    :: map_wf(ndf_wf)
    integer(kind=i_def), intent(in)    :: map_w2(ndf_w2)
    integer(kind=i_def), intent(in)    :: map_wc(ndf_wc)
    real(kind=r_tran),   intent(inout) :: field(undf_wf)
    real(kind=r_tran),   intent(in)    :: ls_field(undf_wf)
    real(kind=r_tran),   intent(in)    :: dep_dist_pert(undf_w2)
    real(kind=r_tran),   intent(in)    :: cubic_coef_1(undf_wc)
    real(kind=r_tran),   intent(in)    :: cubic_coef_2(undf_wc)
    real(kind=r_tran),   intent(in)    :: cubic_coef_3(undf_wc)
    real(kind=r_tran),   intent(in)    :: cubic_coef_4(undf_wc)
    integer(kind=i_def), intent(in)    :: cubic_indices_1(undf_wc)
    integer(kind=i_def), intent(in)    :: cubic_indices_2(undf_wc)
    integer(kind=i_def), intent(in)    :: cubic_indices_3(undf_wc)
    integer(kind=i_def), intent(in)    :: cubic_indices_4(undf_wc)

    ! Local arrays
    real(kind=r_tran) :: field_local(nlayers+ndf_wf-1,4)
    real(kind=r_tran) :: ls_field_local(nlayers+ndf_wf-1,2)
    real(kind=r_tran) :: field_dep(nlayers+ndf_wf-1)
    real(kind=r_tran) :: grad_ls_field(nlayers+ndf_wf-1)
    real(kind=r_tran) :: pert_dist(nlayers+ndf_wf-1)

    ! Indices
    integer(kind=i_def) :: k, nl, wf_idx, wc_idx, w2_idx

    ! nl = nlayers    for w3
    !    = nlayers+1  for wtheta
    nl = nlayers + ndf_wf - 1
    wf_idx = map_wf(1)
    w2_idx = map_w2(1)
    wc_idx = map_wc(1)

    ! Create local arrays
    do k = 1, nl
      field_local(k,1) = field(wf_idx + cubic_indices_1(wc_idx+k-1) - 1)
      field_local(k,2) = field(wf_idx + cubic_indices_2(wc_idx+k-1) - 1)
      field_local(k,3) = field(wf_idx + cubic_indices_3(wc_idx+k-1) - 1)
      field_local(k,4) = field(wf_idx + cubic_indices_4(wc_idx+k-1) - 1)
      ! Require indices 2 and 3 for ls_field as these lie around
      ! the departure point
      ls_field_local(k,1) = ls_field(wf_idx + cubic_indices_2(wc_idx+k-1) - 1)
      ls_field_local(k,2) = ls_field(wf_idx + cubic_indices_3(wc_idx+k-1) - 1)
    end do

    ! Interpolate field
    field_dep(:) = (                                                           &
        cubic_coef_1(wc_idx : wc_idx+nl-1)*field_local(:,1)                    &
        + cubic_coef_2(wc_idx : wc_idx+nl-1)*field_local(:,2)                  &
        + cubic_coef_3(wc_idx : wc_idx+nl-1)*field_local(:,3)                  &
        + cubic_coef_4(wc_idx : wc_idx+nl-1)*field_local(:,4)                  &
    )

    ! Compute gradient of ls_field in departure cell
    grad_ls_field(:) = ls_field_local(:,2)-ls_field_local(:,1)

    ! Compute pert dist based on whether this is W3 or W3theta field
    if (ndf_wf == 1) then
      ! W3 field so require pert dist averaged to W3 point
      pert_dist(1:nlayers) = ( dep_dist_pert(w2_idx : w2_idx + nlayers - 1)    &
                             + dep_dist_pert(w2_idx + 1 : w2_idx + nlayers) ) / 2.0_r_tran
    else
      ! Wtheta field so can use pert dist at W2v point
      pert_dist(1:nl) = dep_dist_pert(w2_idx : w2_idx + nl - 1)
    end if

    ! Put answer back from local array into global field including
    ! the gradient of the ls_field
    field(wf_idx : wf_idx+nl-1) = field_dep(:) - grad_ls_field(:)*pert_dist(:)

  end subroutine tl_vertical_cubic_sl_code

end module tl_vertical_cubic_sl_kernel_mod