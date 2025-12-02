!-------------------------------------------------------------------------------
!(c) Crown copyright 2023 Met Office. All rights reserved.
!The file LICENCE, distributed with this code, contains details of the terms
!under which the code may be used.
!-------------------------------------------------------------------------------
!> @brief Compute the SKEB2 rotational increments
module skeb_durot_kernel_mod

  use argument_mod,      only: arg_type, GH_FIELD,          &
                               GH_REAL, GH_WRITE, GH_READ,  &
                               CELL_COLUMN, GH_INTEGER,     &
                               GH_SCALAR, STENCIL, CROSS
  use fs_continuity_mod, only: W1, W2
  use constants_mod,     only: r_def, i_def
  use kernel_mod,        only: kernel_type

  implicit none

  !> Kernel metadata for Psyclone
  type, public, extends(kernel_type) :: skeb_durot_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/        &
    arg_type(GH_FIELD, GH_REAL, GH_WRITE, W2), & ! du_rot_w2
    arg_type(GH_FIELD, GH_REAL, GH_READ, W1),  & ! psif_hat
    arg_type(GH_FIELD, GH_REAL, GH_READ,  W2), & ! dx_in_w2
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ ), & ! skeb_level_bottom
    arg_type(GH_SCALAR, GH_INTEGER, GH_READ )  & ! skeb_level_top
    /)
    integer :: operates_on = CELL_COLUMN

  contains
    procedure, nopass :: skeb_durot_code
  end type skeb_durot_kernel_type

  public skeb_durot_code
contains

  !> @brief Calculate SKEB2 rotational increments
  !> @param[in]    nlayers     The number of layers
  !> @param[out]   du_rot_w2   Rotational wind increments
  !> @param[in]    psif_hat    Stremfunction
  !> @param[in]    dx_in_w2    delta-x and w2 dofs
  !> @param[in]    skeb_level_bottom Bottom SKEB level
  !> @param[in]    skeb_level_top    Top SKEB level
  !> @param[in]    ndf_w2      Number of DOFs per cell for w2 space
  !> @param[in]    undf_w2     Number of unique DOFs  for w2 space
  !> @param[in]    map_w2      dofmap for the cell at the base of the column for w2 space
  !> @param[in]    ndf_w1      Number of DOFs per cell for w1 space
  !> @param[in]    undf_w1     Number of unique DOFs  for w1 space
  !> @param[in]    map_w1      dofmap for the cell at the base of the column for w1 space

  subroutine skeb_durot_code(nlayers,           &
                             du_rot_w2,         &
                             psif_hat,          &
                             dx_at_w2,          &
                             skeb_level_bottom, &
                             skeb_level_top,    &
                             ndf_w2,            &
                             undf_w2,           &
                             map_w2,            &
                             ndf_w1,            &
                             undf_w1,           &
                             map_w1)

    implicit none

    !Arguments
    integer(kind=i_def), intent(in) :: nlayers
    integer(kind=i_def), intent(in) :: ndf_w1, ndf_w2
    integer(kind=i_def), intent(in) :: undf_w1, undf_w2
    integer(kind=i_def), intent(in),  dimension(ndf_w1)   :: map_w1
    integer(kind=i_def), intent(in),  dimension(ndf_w2)   :: map_w2

    ! Fields
    real(kind=r_def),    intent(in),    dimension(undf_w1)  :: psif_hat
    real(kind=r_def),    intent(in),    dimension(undf_w2)  :: dx_at_w2
    real(kind=r_def),    intent(inout), dimension(undf_w2)  :: du_rot_w2

    ! Scalars
    integer(kind=i_def), intent(in) :: skeb_level_bottom, skeb_level_top

    integer(kind=i_def) :: k, df, fac, dfp3, dfm1

!
! Here we're calculating:
!
! (du, dv) = (-dpsi/dy, dpsi/dx)
!
! Calculated as a finite difference using the following:
! x-comp = j-1 - j
! dof 3 this is 6 - 7 = df+3 - df+4
! dof 1 this is 5 - 8 = df+4 - df+3(+4)
! y-comp = i+1 - i but because beta = -y we do i - i+1
! dof 4 this is 8 - 7 = df+4 - df+3
! dof 2 this is 5 - 6 = df+3 - df+4

    do df = 1,4

      if (df == 1 .or. df == 4) then
        fac = -1.0_r_def
      else
        fac = 1.0_r_def
      end if

      if (df == 1) then
        dfp3 = df+7
        dfm1 = df+3
      else
        dfp3 = df+3
        dfm1 = df-1
      end if

      do k = skeb_level_bottom, skeb_level_top


        du_rot_w2(map_w2(df)+k-1) = fac*(psif_hat(map_w1(dfp3)+k-1) - &
                                         psif_hat(map_w1(df+4)+k-1)) / &
                                         dx_at_w2(map_w2(dfm1)+k-1)

      end do

    end do

  end subroutine skeb_durot_code

end module skeb_durot_kernel_mod
