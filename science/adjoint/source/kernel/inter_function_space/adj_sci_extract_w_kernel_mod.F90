!-----------------------------------------------------------------------------
! (c) Crown copyright Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Adjoint of the extraction of the vertical component of the wind in w2
!>        onto wtheta.
!>
module adj_sci_extract_w_kernel_mod

  use, intrinsic :: iso_fortran_env, only: real32, real64

  use argument_mod,      only: arg_type,             &
                               GH_FIELD, GH_REAL,    &
                               GH_READWRITE, GH_INC, &
                               CELL_COLUMN
  use constants_mod,     only: i_def
  use fs_continuity_mod, only: W2, Wtheta
  use kernel_mod,        only: kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  !> The type declaration for the kernel. Contains the metadata needed by the
  !> Psy layer.
  !>
  type, public, extends(kernel_type) :: adj_extract_w_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                     &
         arg_type(GH_FIELD, GH_REAL, GH_READWRITE, Wtheta), &
         arg_type(GH_FIELD, GH_REAL, GH_INC,       W2)      &
         /)
    integer :: operates_on = CELL_COLUMN
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: adj_extract_w_code

  ! Generic interface for real32 and real64 types
  interface adj_extract_w_code
    module procedure  &
      adj_extract_w_code_real32, &
      adj_extract_w_code_real64
  end interface


contains

!> @brief The adjoint extraction subroutine which is called directly by the psy layer
!! @param[in] nlayers Integer the number of layers
!! @param[in,out] velocity_w2v  component of u_physics normal to cell top
!! @param[in,out] u_physics Real array winds for physics
!! @param[in] ndf_wth The number of degrees of freedom per cell for wth
!! @param[in] undf_wth The number of unique degrees of freedom for wth
!! @param[in] map_wth Integer array holding the dofmap for the cell at the
!>            base of the column for wth
!! @param[in] ndf_w2 The number of degrees of freedom per cell for w2
!! @param[in] undf_w2 The number of unique degrees of freedom for w2
!! @param[in] map_w2 Integer array holding the dofmap for the cell at the
!>            base of the column for w2
subroutine adj_extract_w_code_real32( nlayers,                    &
                                      velocity_w2v,               &
                                      u_physics,                  &
                                      ndf_wth, undf_wth, map_wth, &
                                      ndf_w2, undf_w2, map_w2     &
                                    )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2

  real(kind=real32),   dimension(undf_wth), intent(inout) :: velocity_w2v
  real(kind=real32),   dimension(undf_w2),  intent(inout) :: u_physics
  integer(kind=i_def), dimension(ndf_wth),  intent(in)    :: map_wth
  integer(kind=i_def), dimension(ndf_w2),   intent(in)    :: map_w2

  ! Internal variables
  integer(kind=i_def) :: k

  do k = nlayers, 0, -1

    u_physics(map_w2(5) + k ) = u_physics(map_w2(5) + k ) &
                              + velocity_w2v(map_wth(1) + k)
    velocity_w2v(map_wth(1) + k) = 0.0_real32

  end do

end subroutine adj_extract_w_code_real32

subroutine adj_extract_w_code_real64( nlayers,                    &
                                      velocity_w2v,               &
                                      u_physics,                  &
                                      ndf_wth, undf_wth, map_wth, &
                                      ndf_w2, undf_w2, map_w2     &
                                    )

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers

  integer(kind=i_def), intent(in) :: ndf_wth, undf_wth
  integer(kind=i_def), intent(in) :: ndf_w2, undf_w2

  real(kind=real64),   dimension(undf_wth), intent(inout) :: velocity_w2v
  real(kind=real64),   dimension(undf_w2),  intent(inout) :: u_physics
  integer(kind=i_def), dimension(ndf_wth),  intent(in)    :: map_wth
  integer(kind=i_def), dimension(ndf_w2),   intent(in)    :: map_w2

  ! Internal variables
  integer(kind=i_def) :: k

  do k = nlayers, 0, -1

    u_physics(map_w2(5) + k ) = u_physics(map_w2(5) + k ) &
                              + velocity_w2v(map_wth(1) + k)
    velocity_w2v(map_wth(1) + k) = 0.0_real64

  end do

end subroutine adj_extract_w_code_real64

end module adj_sci_extract_w_kernel_mod
