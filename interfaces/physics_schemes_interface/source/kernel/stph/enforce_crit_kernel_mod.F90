!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Resets the value of a field if it exceeds a certain critical
!>        value, and also halves the value of any surrounding points
module enforce_crit_kernel_mod

  use argument_mod,  only : arg_type,                  &
                            GH_FIELD, GH_SCALAR,       &
                            GH_REAL, GH_READ,          &
                            GH_WRITE,                  &
                            STENCIL, REGION,           &
                            ANY_DISCONTINUOUS_SPACE_1, &
                            CELL_COLUMN, GH_INTEGER
  use constants_mod, only : i_def, r_def, r_double
  use kernel_mod,    only : kernel_type

  implicit none

  private

  !---------------------------------------------------------------------------
  ! Public types
  !---------------------------------------------------------------------------
  ! The type declaration for the kernel. Contains the metadata needed by the
  ! Psy layer.
  !
  type, public, extends(kernel_type) :: enforce_crit_kernel_type
    private
    type(arg_type) :: meta_args(5) = (/                                                    &
         arg_type(GH_FIELD, GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),                 &
         arg_type(GH_FIELD, GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1, STENCIL(REGION)), &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                                            &
         arg_type(GH_SCALAR, GH_REAL, GH_READ),                                            &
         arg_type(GH_SCALAR, GH_INTEGER, GH_READ)                                          &
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: enforce_crit_code
  end type

  !---------------------------------------------------------------------------
  ! Contained functions/subroutines
  !---------------------------------------------------------------------------
  public :: enforce_crit_code


contains

!> @brief Resets the value of a field to a given level if it exceeds a
!!        certain critical value. Also halves the value of any neighbouring
!!        points
!! @param[in] nlayers Number of layers
!! @param[in,out] field_out Field to be modified
!! @param[in] field_in The unmodified input field
!! @param[in] stencil_size The size of the field_in stencil
!! @param[in] stencil_map  The map for the field_in stencil points
!! @param[in] crit The upper bound
!! @param[in] enforce_val The value to reset the field to
!! @param[in] level2km  The k value of the level 2km above the surface
!! @param[in] ndf Number of degrees of freedom per cell
!! @param[in] undf Total number of degrees of freedom
!! @param[in] map Dofmap for the cell at the base of the column

subroutine enforce_crit_code(nlayers, field_out, field_in,    &
                             stencil_size, stencil_map, crit, &
                             enforce_val, level2km,           &
                             ndf, undf, map)

  implicit none

  ! Arguments
  integer(kind=i_def), intent(in) :: nlayers, ndf, undf, stencil_size, level2km
  integer(kind=i_def), dimension(ndf), intent(in) :: map
  real(kind=r_def), dimension(undf), intent(inout) :: field_out
  real(kind=r_def), dimension(undf), intent(in) :: field_in
  real(kind=r_def), intent(in) :: crit
  real(kind=r_def), intent(in) :: enforce_val
  integer(kind=i_def), dimension(ndf,stencil_size), intent(in) :: stencil_map

  ! Internal variables
  integer(kind=i_def) :: df, k, dfs
  real(kind=r_def) :: logscale,levfac, local_crit

  logscale = 10.0_r_def/level2km

  do df = 1, ndf

    do k = 0, nlayers-1
      levfac = MAX(0.1_r_def, LOG10(logscale*real(k+1,r_def)))
      local_crit = crit/levfac

      ! Clip field
      if (field_in(map(df)+k) > local_crit) then
        field_out(map(df)+k) = enforce_val
      else
        do dfs = 2, stencil_size
          if (field_in(stencil_map(df,dfs)+k) > local_crit) then
            field_out(map(df)+k) = 0.5_r_def*field_in(map(df)+k)
          end if
        end do
      end if
    end do
  end do

end subroutine enforce_crit_code

end module enforce_crit_kernel_mod
