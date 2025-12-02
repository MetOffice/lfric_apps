!-------------------------------------------------------------------------------
! Copyright (c) 2024,  Met Office, on behalf of HMSO and Queen's Printer
! For further details please refer to the file LICENCE.original which you
! should have received as part of this distribution.
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Add a random perturbation to a supplied scalar field.
!> @details Samples two arrays of uniformly distributed random numbers, performs
!! a Box-Muller transformation to obtain an array distributed by the Gaussian
!! distribution, and adds it to the column.

module random_perturb_field_kernel_mod

  use argument_mod,               only : arg_type,                              &
                                         GH_FIELD, GH_SCALAR,                   &
                                         GH_REAL, GH_INTEGER,                   &
                                         GH_READ, GH_READWRITE,                 &
                                         CELL_COLUMN,                           &
                                         ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,              only : r_def, i_def, pi
  use kernel_mod,                 only : kernel_type

  use log_mod,                   only : log_event, log_scratch_space,     &
                                        LOG_LEVEL_INFO

  implicit none

  !-----------------------------------------------------------------------------
  ! Public types
  !-----------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !! Psy layer
  type, public, extends(kernel_type) :: random_perturb_field_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                                         &
          arg_type(GH_FIELD, GH_REAL, GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1), &
          arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                              &
          /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass ::random_perturb_field_code
  end type

  !-----------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-----------------------------------------------------------------------------
  public random_perturb_field_code
  contains

  subroutine random_perturb_field_code(nlayers, field, perturb_magnitude,       &
                                       ndf_field, undf_field, map_field)

    implicit none

    integer(kind=i_def), intent(in)     :: nlayers
    integer(kind=i_def), intent(in)     :: ndf_field
    integer(kind=i_def), intent(in)     :: undf_field

    integer(kind=i_def), intent(in),    dimension(ndf_field)  :: map_field
    real(kind=r_def),    intent(inout), dimension(undf_field) :: field

    integer(kind=i_def), intent(in)     :: perturb_magnitude

    integer(kind=i_def) :: k, df, top_df
    real(kind=r_def)    :: pert, u

    ! Assume lowest order scalar field
    df = 1

    ! If the field is Wtheta loop to the very top
    top_df = nlayers + ndf_field - 2

    do k = 0, top_df
      ! Get the uniformly distributed numbers
      call random_number(pert)
      call random_number(u)

      ! Apply the Box-Muller transformation (0 <= u < 1 so this LOG is
      ! safe). Reuse pert for this
      pert = SQRT(-2.0_r_def*LOG(1.0_r_def - u))*COS(2.0_r_def*pi*pert)

      ! Apply the magnitude scaling defined in the namelist options
      pert = pert*10.0_r_def**perturb_magnitude

      ! Add the perturbation
      field(map_field(df) + k) = field(map_field(df) + k) + pert
    end do

  end subroutine random_perturb_field_code

end module random_perturb_field_kernel_mod