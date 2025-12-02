!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
!> @brief Kernel to do a vertical cubic semi-Lagragian advection of a field
!!        in the vertical direction.
!> @Details The 1D vertical advective transport equation for a W3 variable
!!          is solved using a cubic semi-Lagragian advection scheme.

module vertical_cubic_sl_kernel_mod

  use argument_mod,          only : arg_type,                  &
                                    GH_FIELD, GH_SCALAR,       &
                                    GH_REAL, GH_INTEGER,       &
                                    GH_READWRITE, GH_READ,     &
                                    CELL_COLUMN, GH_LOGICAL,   &
                                    ANY_DISCONTINUOUS_SPACE_1, &
                                    ANY_DISCONTINUOUS_SPACE_2
  use constants_mod,         only : r_tran, i_def, l_def, EPS_R_TRAN
  use kernel_mod,            only : kernel_type
  use transport_enumerated_types_mod, only : monotone_none,                    &
                                             vertical_monotone_order_constant, &
                                             vertical_monotone_order_linear,   &
                                             vertical_monotone_order_high
  use sl_support_mod,   only : monotone_cubic_sl
  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed
  !>                                      by the PSy layer.
  type, public, extends(kernel_type) :: vertical_cubic_sl_kernel_type
    private
    type(arg_type) :: meta_args(14) = (/                                            &
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! field
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! cubic_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! cubic_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! cubic_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! cubic_coef
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! cubic_indices
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! cubic_indices
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! cubic_indices
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! cubic_indices
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! linear_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! linear_coef
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                 & ! monotone scheme
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                 & ! monotone order
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                                  & ! log_space
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: vertical_cubic_sl_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: vertical_cubic_sl_code

  contains

  !-------------------------------------------------------------------------------
  !> @details This kernel interpolates the field to the
  !!          departure point using 1d-Cubic-Lagrange interpolation.
  !> @param[in]     nlayers         The number of layers
  !> @param[in,out] field           The field to be advected
  !> @param[in]     cubic_coef      The cubic interpolation coefficients (1-4)
  !> @param[in]     cubic_indices   The cubic interpolation indices (1-4)
  !> @param[in]     linear_coef     The linear interpolation coefficients (1-2,
  !!                                used for monotonicity)
  !> @param[in]     vertical_monotone
  !!                                Enumerator indicating the monotone scheme
  !> @param[in]     vertical_monotone_order
  !!                                Order of the monotone scheme
  !> @param[in]     log_space       Switch to use natural logarithmic space
  !!                                for the SL interpolation
  !> @param[in]     ndf_wf          Num Dofs per cell for the field
  !> @param[in]     undf_wf         Num Dofs in this partition for the field
  !> @param[in]     map_wf          Dofmap for the field
  !> @param[in]     ndf_wc          Num Dofs per cell for the coefficients
  !> @param[in]     undf_wc         Num Dofs per cell in this partition
  !!                                for the coefficients
  !> @param[in]     map_wc          Dofmap for the coefficients
  !-------------------------------------------------------------------------------
  subroutine vertical_cubic_sl_code( nlayers,                 &
                                     field,                   &
                                     cubic_coef_1,            &
                                     cubic_coef_2,            &
                                     cubic_coef_3,            &
                                     cubic_coef_4,            &
                                     cubic_indices_1,         &
                                     cubic_indices_2,         &
                                     cubic_indices_3,         &
                                     cubic_indices_4,         &
                                     linear_coef_1,           &
                                     linear_coef_2,           &
                                     vertical_monotone,       &
                                     vertical_monotone_order, &
                                     log_space,               &
                                     ndf_wf, undf_wf, map_wf, &
                                     ndf_wc, undf_wc, map_wc )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ndf_wf
    integer(kind=i_def), intent(in)    :: undf_wf
    integer(kind=i_def), intent(in)    :: ndf_wc
    integer(kind=i_def), intent(in)    :: undf_wc
    integer(kind=i_def), intent(in)    :: map_wf(ndf_wf)
    integer(kind=i_def), intent(in)    :: map_wc(ndf_wc)
    real(kind=r_tran),   intent(inout) :: field(undf_wf)
    real(kind=r_tran),   intent(in)    :: cubic_coef_1(undf_wc)
    real(kind=r_tran),   intent(in)    :: cubic_coef_2(undf_wc)
    real(kind=r_tran),   intent(in)    :: cubic_coef_3(undf_wc)
    real(kind=r_tran),   intent(in)    :: cubic_coef_4(undf_wc)
    integer(kind=i_def), intent(in)    :: cubic_indices_1(undf_wc)
    integer(kind=i_def), intent(in)    :: cubic_indices_2(undf_wc)
    integer(kind=i_def), intent(in)    :: cubic_indices_3(undf_wc)
    integer(kind=i_def), intent(in)    :: cubic_indices_4(undf_wc)
    real(kind=r_tran),   intent(in)    :: linear_coef_1(undf_wc)
    real(kind=r_tran),   intent(in)    :: linear_coef_2(undf_wc)
    integer(kind=i_def), intent(in)    :: vertical_monotone
    integer(kind=i_def), intent(in)    :: vertical_monotone_order
    logical(kind=l_def), intent(in)    :: log_space

    ! Local arrays
    real(kind=r_tran) :: field_local(nlayers+ndf_wf-1,4)
    real(kind=r_tran) :: field_dep(nlayers+ndf_wf-1)

    real(kind=r_tran), allocatable :: log_field_local(:,:)

    ! Indices
    integer(kind=i_def) :: k, nl, wf_idx, wc_idx

    ! nl = nlayers    for w3
    !    = nlayers+1  for wtheta
    nl = nlayers + ndf_wf - 1
    wf_idx = map_wf(1)
    wc_idx = map_wc(1)

    do k = 1, nl
      field_local(k,1) = field(wf_idx + cubic_indices_1(wc_idx+k-1) - 1)
      field_local(k,2) = field(wf_idx + cubic_indices_2(wc_idx+k-1) - 1)
      field_local(k,3) = field(wf_idx + cubic_indices_3(wc_idx+k-1) - 1)
      field_local(k,4) = field(wf_idx + cubic_indices_4(wc_idx+k-1) - 1)
    end do

    if (log_space) then
      allocate(log_field_local(nl,4))
      log_field_local(:,:) = LOG(MAX(ABS(field_local(:,:)), EPS_R_TRAN))

      ! Do interpolation on log(field) using cubic coefficients and indices
      field_dep(:) = EXP(                                                      &
          cubic_coef_1(wc_idx : wc_idx+nl-1)*log_field_local(:,1)              &
          + cubic_coef_2(wc_idx : wc_idx+nl-1)*log_field_local(:,2)            &
          + cubic_coef_3(wc_idx : wc_idx+nl-1)*log_field_local(:,3)            &
          + cubic_coef_4(wc_idx : wc_idx+nl-1)*log_field_local(:,4)            &
      )

      deallocate(log_field_local)

    else
      ! Interpolate field as is
      field_dep(:) = (                                                         &
          cubic_coef_1(wc_idx : wc_idx+nl-1)*field_local(:,1)                  &
          + cubic_coef_2(wc_idx : wc_idx+nl-1)*field_local(:,2)                &
          + cubic_coef_3(wc_idx : wc_idx+nl-1)*field_local(:,3)                &
          + cubic_coef_4(wc_idx : wc_idx+nl-1)*field_local(:,4)                &
      )
    end if

    ! Enforce monotonicity if required
    if ( vertical_monotone /= monotone_none ) then
      ! Apply monotonicity
      call monotone_cubic_sl(                                                  &
              field_dep, field_local,                                          &
              linear_coef_1(wc_idx : wc_idx+nl-1),                             &
              linear_coef_2(wc_idx : wc_idx+nl-1),                             &
              vertical_monotone, vertical_monotone_order, nl                   &
      )
    end if

    ! Put answer back from local array into global field
    field(wf_idx : wf_idx+nl-1) = field_dep(:)

  end subroutine vertical_cubic_sl_code

end module vertical_cubic_sl_kernel_mod