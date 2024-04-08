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
  use transport_enumerated_types_mod, only : vertical_monotone_none,           &
                                             vertical_monotone_strict,         &
                                             vertical_monotone_relaxed,        &
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
    type(arg_type) :: meta_args(18) = (/                                            &
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
  !> @param[in,out] field           The field at time level n interpolated to the departure point
  !> @param[in]     quintic_coef    The quintic interpolation coefficients (1-6)
  !> @param[in]     quintic_indices The quintic interpolation indices (1-6)
  !> @param[in]     linear_coef     The linear interpolation coefficients (1-2, used for monotonicity)
  !> @param[in]     vertical_monotone
  !!                                The monotone scheme
  !> @param[in]     vertical_monotone_order
  !!                                Order of the monotone scheme
  !> @param[in]     log_space       Switch to use natural logarithmic space
  !!                                for the SL interpolation
  !> @param[in]     ndf_wf          The number of degrees of freedom per cell
  !!                                for the field (i.e. w3/wtheta) space
  !> @param[in]     undf_wf         The number of unique degrees of freedom
  !!                                for the field (i.e. w3/wtheta) space
  !> @param[in]     map_wf          The dofmap for the cell at the base of the column
  !!                                for the field (i.e. w3/wtheta) space
  !> @param[in]     ndf_wq          The number of degrees of freedom per cell
  !!                                for the coefficients space
  !> @param[in]     undf_wq         The number of unique degrees of freedom
  !!                                for the coefficients space
  !> @param[in]     map_wq          The dofmap for the cell at the base of the column
  !!                                for the coefficients space
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
    integer(kind=i_def),                     intent(in)    :: nlayers
    integer(kind=i_def),                     intent(in)    :: ndf_wf
    integer(kind=i_def),                     intent(in)    :: undf_wf
    integer(kind=i_def),                     intent(in)    :: ndf_wq
    integer(kind=i_def),                     intent(in)    :: undf_wq
    integer(kind=i_def), dimension(ndf_wf),  intent(in)    :: map_wf
    integer(kind=i_def), dimension(ndf_wq),  intent(in)    :: map_wq
    real(kind=r_tran),   dimension(undf_wf), intent(inout) :: field
    real(kind=r_tran),   dimension(undf_wq), intent(in)    :: quintic_coef_1
    real(kind=r_tran),   dimension(undf_wq), intent(in)    :: quintic_coef_2
    real(kind=r_tran),   dimension(undf_wq), intent(in)    :: quintic_coef_3
    real(kind=r_tran),   dimension(undf_wq), intent(in)    :: quintic_coef_4
    real(kind=r_tran),   dimension(undf_wq), intent(in)    :: quintic_coef_5
    real(kind=r_tran),   dimension(undf_wq), intent(in)    :: quintic_coef_6
    integer(kind=i_def), dimension(undf_wq), intent(in)    :: quintic_indices_1
    integer(kind=i_def), dimension(undf_wq), intent(in)    :: quintic_indices_2
    integer(kind=i_def), dimension(undf_wq), intent(in)    :: quintic_indices_3
    integer(kind=i_def), dimension(undf_wq), intent(in)    :: quintic_indices_4
    integer(kind=i_def), dimension(undf_wq), intent(in)    :: quintic_indices_5
    integer(kind=i_def), dimension(undf_wq), intent(in)    :: quintic_indices_6
    real(kind=r_tran),   dimension(undf_wq), intent(in)    :: linear_coef_1
    real(kind=r_tran),   dimension(undf_wq), intent(in)    :: linear_coef_2
    integer(kind=i_def), intent(in)  :: vertical_monotone,  &
                                        vertical_monotone_order
    logical(kind=l_def), intent(in)  :: log_space

    ! Local arrays
    real(kind=r_tran),   dimension(nlayers+1)   :: f0, fd
    real(kind=r_tran),   dimension(6,nlayers+1) :: qc
    integer(kind=i_def), dimension(6,nlayers+1) :: sq
    real(kind=r_tran),   dimension(2,nlayers+1) :: cl

    ! Indices
    integer(kind=i_def) :: k, nz, nl

    ! nl = nlayers-1  for w3
    !    = nlayers    for wtheta
    nl = nlayers - 1 + (ndf_wf - 1)
    nz = nl + 1

    ! Map global field into 1d-array f0
    ! Coeffs are on a layer-first multidata field
    ! so the indices are: map_md(1) + index*nz + k
    do k = 0, nl
      f0(k+1) = field(map_wf(1)+k)
      qc(1,k+1) = quintic_coef_1(map_wq(1)+k)
      qc(2,k+1) = quintic_coef_2(map_wq(1)+k)
      qc(3,k+1) = quintic_coef_3(map_wq(1)+k)
      qc(4,k+1) = quintic_coef_4(map_wq(1)+k)
      qc(5,k+1) = quintic_coef_5(map_wq(1)+k)
      qc(6,k+1) = quintic_coef_6(map_wq(1)+k)
      sq(1,k+1) = quintic_indices_1(map_wq(1)+k)
      sq(2,k+1) = quintic_indices_2(map_wq(1)+k)
      sq(3,k+1) = quintic_indices_3(map_wq(1)+k)
      sq(4,k+1) = quintic_indices_4(map_wq(1)+k)
      sq(5,k+1) = quintic_indices_5(map_wq(1)+k)
      sq(6,k+1) = quintic_indices_6(map_wq(1)+k)
      cl(1,k+1) = linear_coef_1(map_wq(1)+k)
      cl(2,k+1) = linear_coef_2(map_wq(1)+k)
    end do

    ! Apply log to f0 if required
    ! If using the log_space option, f0 is forced to be positive
    if (log_space) then
      do k = 1, nz
        f0(k) = log(max(EPS_R_TRAN,abs(f0(k))))
      end do
    end if

    ! Do field interpolation using cubic coefficients and indices

    do k = 1, nz
      fd(k) = qc(1,k)*f0(sq(1,k)) + qc(2,k)*f0(sq(2,k)) + &
              qc(3,k)*f0(sq(3,k)) + qc(4,k)*f0(sq(4,k)) + &
              qc(5,k)*f0(sq(5,k)) + qc(6,k)*f0(sq(6,k))
    end do

    ! If using log_space then convert back
    if (log_space) then
      do k = 1, nz
        fd(k) = exp(fd(k))
        f0(k) = field(map_wf(1)+k-1)
      end do
    end if

    ! Enforce monotonicity if required
    if ( vertical_monotone /= vertical_monotone_none ) then
      ! Apply monotonicity
      call monotone_quintic_sl( fd,f0,sq,cl,vertical_monotone, &
                                vertical_monotone_order,1,nz )
    end if

    ! Remap the column answer back to the global data

    do k=0,nl
      field(map_wf(1)+k) = fd(k+1)
    end do

  end subroutine vertical_quintic_sl_code

end module vertical_quintic_sl_kernel_mod