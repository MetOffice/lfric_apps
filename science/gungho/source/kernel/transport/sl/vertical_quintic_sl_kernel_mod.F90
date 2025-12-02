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

module vertical_quintic_sl_kernel_mod

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
                                             monotone_strict,                  &
                                             monotone_relaxed,                 &
                                             vertical_monotone_order_constant, &
                                             vertical_monotone_order_linear,   &
                                             vertical_monotone_order_high
  use sl_support_mod,   only : monotone_quintic_sl
  implicit none

  private

  !-------------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed
  !>                                      by the PSy layer.
  type, public, extends(kernel_type) :: vertical_quintic_sl_kernel_type
    private
    type(arg_type) :: meta_args(18) = (/                                           &
         arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), & ! field
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! quintic_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! quintic_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! quintic_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! quintic_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! quintic_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! quintic_coef
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! quintic_indices
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! quintic_indices
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! quintic_indices
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! quintic_indices
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! quintic_indices
         arg_type(GH_FIELD,  GH_INTEGER, GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! quintic_indices
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! linear_coef
         arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_2), & ! linear_coef
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                 & ! monotone scheme
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ),                                 & ! monotone order
         arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                                  & ! log_space
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: vertical_quintic_sl_code
  end type

  !-------------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------------
  public :: vertical_quintic_sl_code

  contains

  !-------------------------------------------------------------------------------
  !> @details This kernel interpolates the field to the
  !!          departure point using 1d-Cubic-Lagrange interpolation.
  !> @param[in]     nlayers         The number of layers
  !> @param[in,out] field           The field to be advected
  !> @param[in]     quintic_coef    The quintic interpolation coefficients (1-6)
  !> @param[in]     quintic_indices The quintic interpolation indices (1-6)
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
  !> @param[in]     ndf_wq          Num Dofs per cell for the coefficients
  !> @param[in]     undf_wq         Num Dofs per cell in this partition
  !!                                for the coefficients
  !> @param[in]     map_wq          Dofmap for the coefficients
  !-------------------------------------------------------------------------------
  subroutine vertical_quintic_sl_code( nlayers,                 &
                                       field,                   &
                                       quintic_coef_1,          &
                                       quintic_coef_2,          &
                                       quintic_coef_3,          &
                                       quintic_coef_4,          &
                                       quintic_coef_5,          &
                                       quintic_coef_6,          &
                                       quintic_indices_1,       &
                                       quintic_indices_2,       &
                                       quintic_indices_3,       &
                                       quintic_indices_4,       &
                                       quintic_indices_5,       &
                                       quintic_indices_6,       &
                                       linear_coef_1,           &
                                       linear_coef_2,           &
                                       vertical_monotone,       &
                                       vertical_monotone_order, &
                                       log_space,               &
                                       ndf_wf, undf_wf, map_wf, &
                                       ndf_wq, undf_wq, map_wq )

    implicit none

    ! Arguments
    integer(kind=i_def), intent(in)    :: nlayers
    integer(kind=i_def), intent(in)    :: ndf_wf
    integer(kind=i_def), intent(in)    :: undf_wf
    integer(kind=i_def), intent(in)    :: ndf_wq
    integer(kind=i_def), intent(in)    :: undf_wq
    integer(kind=i_def), intent(in)    :: map_wf(ndf_wf)
    integer(kind=i_def), intent(in)    :: map_wq(ndf_wq)
    real(kind=r_tran),   intent(inout) :: field(undf_wf)
    real(kind=r_tran),   intent(in)    :: quintic_coef_1(undf_wq)
    real(kind=r_tran),   intent(in)    :: quintic_coef_2(undf_wq)
    real(kind=r_tran),   intent(in)    :: quintic_coef_3(undf_wq)
    real(kind=r_tran),   intent(in)    :: quintic_coef_4(undf_wq)
    real(kind=r_tran),   intent(in)    :: quintic_coef_5(undf_wq)
    real(kind=r_tran),   intent(in)    :: quintic_coef_6(undf_wq)
    integer(kind=i_def), intent(in)    :: quintic_indices_1(undf_wq)
    integer(kind=i_def), intent(in)    :: quintic_indices_2(undf_wq)
    integer(kind=i_def), intent(in)    :: quintic_indices_3(undf_wq)
    integer(kind=i_def), intent(in)    :: quintic_indices_4(undf_wq)
    integer(kind=i_def), intent(in)    :: quintic_indices_5(undf_wq)
    integer(kind=i_def), intent(in)    :: quintic_indices_6(undf_wq)
    real(kind=r_tran),   intent(in)    :: linear_coef_1(undf_wq)
    real(kind=r_tran),   intent(in)    :: linear_coef_2(undf_wq)
    integer(kind=i_def), intent(in)    :: vertical_monotone
    integer(kind=i_def), intent(in)    :: vertical_monotone_order
    logical(kind=l_def), intent(in)    :: log_space

    ! Local arrays
    real(kind=r_tran) :: field_local(nlayers+ndf_wf-1,6)
    real(kind=r_tran) :: field_dep(nlayers+ndf_wf-1)

    real(kind=r_tran), allocatable :: log_field_local(:,:)

    ! Indices
    integer(kind=i_def) :: k, nl, wf_idx, wq_idx

    ! nl = nlayers    for w3
    !    = nlayers+1  for wtheta
    nl = nlayers + ndf_wf - 1
    wf_idx = map_wf(1)
    wq_idx = map_wq(1)

    do k = 1, nl
      field_local(k,1) = field(wf_idx + quintic_indices_1(wq_idx+k-1) - 1)
      field_local(k,2) = field(wf_idx + quintic_indices_2(wq_idx+k-1) - 1)
      field_local(k,3) = field(wf_idx + quintic_indices_3(wq_idx+k-1) - 1)
      field_local(k,4) = field(wf_idx + quintic_indices_4(wq_idx+k-1) - 1)
      field_local(k,5) = field(wf_idx + quintic_indices_5(wq_idx+k-1) - 1)
      field_local(k,6) = field(wf_idx + quintic_indices_6(wq_idx+k-1) - 1)
    end do

    if (log_space) then
      allocate(log_field_local(nl,6))
      log_field_local(:,:) = LOG(MAX(ABS(field_local(:,:)), EPS_R_TRAN))

      ! Do interpolation on log(field) using cubic coefficients and indices
      field_dep(:) = EXP(                                                      &
          quintic_coef_1(wq_idx : wq_idx+nl-1)*log_field_local(:,1)            &
          + quintic_coef_2(wq_idx : wq_idx+nl-1)*log_field_local(:,2)          &
          + quintic_coef_3(wq_idx : wq_idx+nl-1)*log_field_local(:,3)          &
          + quintic_coef_4(wq_idx : wq_idx+nl-1)*log_field_local(:,4)          &
          + quintic_coef_5(wq_idx : wq_idx+nl-1)*log_field_local(:,5)          &
          + quintic_coef_6(wq_idx : wq_idx+nl-1)*log_field_local(:,6)          &
      )

      deallocate(log_field_local)

    else
      ! Interpolate field as is
      field_dep(:) = (                                                         &
          quintic_coef_1(wq_idx : wq_idx+nl-1)*field_local(:,1)                &
          + quintic_coef_2(wq_idx : wq_idx+nl-1)*field_local(:,2)              &
          + quintic_coef_3(wq_idx : wq_idx+nl-1)*field_local(:,3)              &
          + quintic_coef_4(wq_idx : wq_idx+nl-1)*field_local(:,4)              &
          + quintic_coef_5(wq_idx : wq_idx+nl-1)*field_local(:,5)              &
          + quintic_coef_6(wq_idx : wq_idx+nl-1)*field_local(:,6)              &
      )
    end if

    ! Enforce monotonicity if required
    if ( vertical_monotone /= monotone_none ) then
      ! Apply monotonicity
      call monotone_quintic_sl(                                                &
              field_dep, field_local,                                          &
              linear_coef_1(wq_idx : wq_idx+nl-1),                             &
              linear_coef_2(wq_idx : wq_idx+nl-1),                             &
              vertical_monotone, vertical_monotone_order, nl                   &
      )
    end if

    ! Put answer back from local array into global field
    field(wf_idx : wf_idx+nl-1) = field_dep(:)

  end subroutine vertical_quintic_sl_code

end module vertical_quintic_sl_kernel_mod