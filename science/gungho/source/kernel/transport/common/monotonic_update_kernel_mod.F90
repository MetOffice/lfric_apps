!-----------------------------------------------------------------------------
! (C) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------

!> @brief Kernel which enforces monotonicity on an advective update.
!> @details Modifies the advective update A(field) such that
!!          \f[ M_\field (field^{n+1} - field^{n}) + \Delta t A(\field) \f]
!!          Returns a monotonic field^{n+1}. The monotonticity is enforced by
!!          ensuring that the the implied field^{n+1} lies within the range
!!          [min(field_s),max(field_s)] for all s in the stencil of values used
!!          to compute field.
module monotonic_update_kernel_mod

use argument_mod,      only : arg_type,                  &
                              GH_FIELD, GH_SCALAR,       &
                              GH_REAL, GH_INTEGER,       &
                              GH_LOGICAL,                &
                              GH_READWRITE, GH_READ,     &
                              ANY_DISCONTINUOUS_SPACE_1, &
                              STENCIL, CROSS, CELL_COLUMN
use constants_mod,     only : r_def, i_def, l_def, r_tran
use kernel_mod,        only : kernel_type

implicit none

private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the PSy layer
type, public, extends(kernel_type) :: monotonic_update_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                                                            &
       arg_type(GH_FIELD,  GH_REAL,    GH_READWRITE, ANY_DISCONTINUOUS_SPACE_1),                 &
       arg_type(GH_FIELD,  GH_REAL,    GH_READ,      ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS)), &
       arg_type(GH_SCALAR, GH_REAL,    GH_READ),                                                 &
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ),                                                 &
       arg_type(GH_SCALAR, GH_LOGICAL, GH_READ)                                                  &
       /)
  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: monotonic_update_code
end type

!-------------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------------
public :: monotonic_update_code

contains

!> @brief Apply a montonic update to a field by clipping it
!> @param[in]     nlayers       Number of layers
!> @param[in,out] adv           Advective update field to apply monotonicity to
!> @param[in]     field         Field to be advected, used to compute the advective output
!> @param[in]     stencil_size  Size of the stencil (number of cells)
!> @param[in]     stencil_map   Dofmaps for the stencil
!> @param[in]     dt            Timestep
!> @param[in]     do_horizontal Do horizontal monotone update
!> @param[in]     do_vertical   Do vertical monotone update
!> @param[in]     ndf           Number of degrees of freedom per cell
!> @param[in]     undf          Number of unique degrees of freedom for the tracer field
!> @param[in]     map           Dofmap for the tracer field
subroutine monotonic_update_code( nlayers,              &
                                  adv,                  &
                                  field,                &
                                  stencil_size,         &
                                  stencil_map,          &
                                  dt,                   &
                                  do_horizontal,        &
                                  do_vertical,          &
                                  ndf,                  &
                                  undf,                 &
                                  map )

  implicit none

  ! Arguments
  integer(kind=i_def),                              intent(in) :: nlayers
  integer(kind=i_def),                              intent(in) :: ndf
  integer(kind=i_def),                              intent(in) :: undf
  integer(kind=i_def),                              intent(in) :: stencil_size
  integer(kind=i_def), dimension(ndf,stencil_size), intent(in) :: stencil_map
  integer(kind=i_def), dimension(ndf),              intent(in) :: map

  logical(kind=l_def), intent(in) :: do_horizontal, do_vertical

  real(kind=r_tran), dimension(undf), intent(inout) :: adv
  real(kind=r_tran), dimension(undf), intent(in)    :: field
  real(kind=r_tran),                  intent(in)    :: dt

  ! Internal variables
  integer(kind=i_def) :: k, cell, p, ijkp, pmin, pmax, nl

  real(kind=r_tran) :: field0, field_max, field_min, a_min, a_max

  ! Upper limit
  ! nl  == nlayers - 1 for W3 spaces
  ! nl  == nlayers for Wtheta spaces
  nl = (nlayers - 1) + (ndf - 1)

  ! Enforce monotonicity of the advective update across the stencil
  do k = 0, nl
    field0    = field(stencil_map(1,1)+k)
    field_max = field(stencil_map(1,1)+k)
    field_min = field(stencil_map(1,1)+k)
    ! Max min across the horzontal stencil (this assumes a cross stencil which
    ! is incorrect for 2D interpolation but should still be a good guess)
    if ( do_horizontal ) then
      do cell = 2,stencil_size
        field_max = max(field_max,field(stencil_map(1,cell)+k))
        field_min = min(field_min,field(stencil_map(1,cell)+k))
      end do
    end if
    if ( do_vertical ) then
      ! Max min across the vertical stencil
      ! Use simple k-1,k+1 limits
      pmin = max(0,  k-1)
      pmax = min(nl, k+1)
      do p = pmin,pmax
        ijkp = stencil_map(1,1) + p
        field_max = max(field_max,field(ijkp))
        field_min = min(field_min,field(ijkp))
      end do
    end if
    a_min = (field0 - field_max)/dt
    a_max = (field0 - field_min)/dt
    adv(map(1)+k) = min(a_max,max(adv(stencil_map(1,1)+k),a_min))
  end do

end subroutine monotonic_update_code

end module monotonic_update_kernel_mod
