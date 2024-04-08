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
  use transport_enumerated_types_mod, only : vertical_monotone_none,           &
                                             vertical_monotone_strict,         &
                                             vertical_monotone_relaxed,        &
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
  !> @param[in]     nlayers       The number of layers
  !> @param[in,out] field         The field at time level n interpolated to the departure point
  !> @param[in]     cubic_coef    The cubic interpolation coefficients (1-4)
  !> @param[in]     cubic_indices The cubic interpolation indices (1-4)
  !> @param[in]     linear_coef   The linear interpolation coefficients (1-2, used for monotonicity)
  !> @param[in]     vertical_monotone
  !!                              The monotone scheme
  !> @param[in]     vertical_monotone_order
  !!                              Order of the monotone scheme
  !> @param[in]     log_space     Switch to use natural logarithmic space
  !!                              for the SL interpolation
  !> @param[in]     ndf_wf        The number of degrees of freedom per cell
  !!                              for the field (i.e. w3/wtheta) space
  !> @param[in]     undf_wf       The number of unique degrees of freedom
  !!                              for the field (i.e. w3/wtheta) space
  !> @param[in]     map_wf        The dofmap for the cell at the base of the column
  !!                              for the field (i.e. w3/wtheta) space
  !> @param[in]     ndf_wc        The number of degrees of freedom per cell
  !!                              for the coefficients space
  !> @param[in]     undf_wc       The number of unique degrees of freedom
  !!                              for the coefficients space
  !> @param[in]     map_wc        The dofmap for the cell at the base of the column
  !!                              for the coefficients space
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
    integer(kind=i_def),                     intent(in)    :: nlayers
    integer(kind=i_def),                     intent(in)    :: ndf_wf
    integer(kind=i_def),                     intent(in)    :: undf_wf
    integer(kind=i_def),                     intent(in)    :: ndf_wc
    integer(kind=i_def),                     intent(in)    :: undf_wc
    integer(kind=i_def), dimension(ndf_wf),  intent(in)    :: map_wf
    integer(kind=i_def), dimension(ndf_wc),  intent(in)    :: map_wc
    real(kind=r_tran),   dimension(undf_wf), intent(inout) :: field
    real(kind=r_tran),   dimension(undf_wc), intent(in)    :: cubic_coef_1
    real(kind=r_tran),   dimension(undf_wc), intent(in)    :: cubic_coef_2
    real(kind=r_tran),   dimension(undf_wc), intent(in)    :: cubic_coef_3
    real(kind=r_tran),   dimension(undf_wc), intent(in)    :: cubic_coef_4
    integer(kind=i_def), dimension(undf_wc), intent(in)    :: cubic_indices_1
    integer(kind=i_def), dimension(undf_wc), intent(in)    :: cubic_indices_2
    integer(kind=i_def), dimension(undf_wc), intent(in)    :: cubic_indices_3
    integer(kind=i_def), dimension(undf_wc), intent(in)    :: cubic_indices_4
    real(kind=r_tran),   dimension(undf_wc), intent(in)    :: linear_coef_1
    real(kind=r_tran),   dimension(undf_wc), intent(in)    :: linear_coef_2
    integer(kind=i_def), intent(in)  :: vertical_monotone,  &
                                        vertical_monotone_order
    logical(kind=l_def), intent(in)  :: log_space

    ! Local arrays
    real(kind=r_tran),   dimension(nlayers+1)   :: f0, fd
    real(kind=r_tran),   dimension(4,nlayers+1) :: cc
    integer(kind=i_def), dimension(4,nlayers+1) :: sc
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
      cc(1,k+1) = cubic_coef_1(map_wc(1)+k)
      cc(2,k+1) = cubic_coef_2(map_wc(1)+k)
      cc(3,k+1) = cubic_coef_3(map_wc(1)+k)
      cc(4,k+1) = cubic_coef_4(map_wc(1)+k)
      sc(1,k+1) = cubic_indices_1(map_wc(1)+k)
      sc(2,k+1) = cubic_indices_2(map_wc(1)+k)
      sc(3,k+1) = cubic_indices_3(map_wc(1)+k)
      sc(4,k+1) = cubic_indices_4(map_wc(1)+k)
      cl(1,k+1) = linear_coef_1(map_wc(1)+k)
      cl(2,k+1) = linear_coef_2(map_wc(1)+k)
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
      fd(k) = cc(1,k)*f0(sc(1,k)) + cc(2,k)*f0(sc(2,k)) + &
              cc(3,k)*f0(sc(3,k)) + cc(4,k)*f0(sc(4,k))
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
      call monotone_cubic_sl( fd,f0,sc,cl,vertical_monotone, &
                              vertical_monotone_order,1,nz )
    end if

    ! Remap the column answer back to the global data

    do k=0,nl
      field(map_wf(1)+k) = fd(k+1)
    end do

  end subroutine vertical_cubic_sl_code

end module vertical_cubic_sl_kernel_mod