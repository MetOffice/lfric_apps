!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief  Contains support functions/routines relating to the Koren scheme
!------------------------------------------------------------------------------
module koren_support_mod

use constants_mod, only : r_tran

implicit none

private

! Public subroutines
public :: interpolate_to_regular_grid

!------------------------------------------------------------------------------
! Contained functions / subroutines
!------------------------------------------------------------------------------
contains

!---------------------------------------------------------------------------------
!> @brief Interpolate 3 regular values from 3 irregular ones for the Koren scheme
!> @details For a cell edge with 2 upwind and 1 downwind values (tracer1,tracer2,tracer3)
!>          with variable cell sizes (dz1,dz2,dz3), this subroutine interpolate three
!>          values (tracer1i,tracer2i,tracer3i) for an equivalent constant
!>          grid-mesh of size h=min(dz1,dz2,dz3), surrounding the same edge, .
!>          using quadratic Lagrange interpolation
!> @param[in]  tracer1 Tracer value for cell 1 (cell 1 is defined as the furthest cell
!>                     from the edge of interest in the upwind direction)
!> @param[in]  tracer2 Tracer value for cell 2 (cell 2 is the cell adjacent to the edge
!>                      in question in the upwind direction)
!> @param[in]  tracer3 Tracer value for cell 3 (cell 3 is the cell adjacent to the edge
!>                     in question in the downwind direction)
!> @param[in]  dz1  Mesh size of cell 1 [tracer1 at centre of dz1]
!> @param[in]  dz2  Mesh size of cell 2 [tracer2 at centre of dz2]
!> @param[in]  dz3  Mesh size of cell 3 [tracer3 at centre of dz3]
!> @param[out] tracer1i Interpolated value for cell 1 of the regular grid
!> @param[out] tracer2i Interpolated value for cell 2 of the regular grid
!> @param[out] tracer3i Interpolated value for cell 3 of the regular grid
!>
!> Explanatory diagram for wind going from cell 1 to cell 3
!>
!>             Irregular grid       ====>   Regular/constant grid
!>            dz1 /= dz2 /= dz3                 grid size h
!>                                            h=min(dz1,dz2,dz3)
!>              +---------+  +
!>              |         |  |
!>              |         |  |               +----------+ +
!>       Cell 1 | tracer1 !  | dz1           |          | |
!>              |         |  |        Cell 1 | tracer1i | | h
!>              |         |  |               |          | |
!>              +---------+  +               +----------+ +
!>              |         |  |               |          | |
!>       Cell 2 | tracer2 |  | dz2    Cell 2 | tracer2i | | h
!>              |         |  |               |          | |
!>              +---------+  + <-----------  +----------+ + <--- Edge of interest
!>              |         |  |               |          | |
!>              |         |  |        Cell 3 | tracer3i | | h
!>              |         |  |               |          | |
!>              |         |  |               +----------+ +
!>       Cell 3 | tracer3 |  | dz3
!>              |         |  |
!>              |         |  |
!>              |         |  |
!>              |         |  |
!>              +---------+  +
!-------------------------------------------------------------------------------
subroutine interpolate_to_regular_grid( tracer1,tracer2,tracer3,dz1,dz2,dz3, &
                                        tracer1i,tracer2i,tracer3i )

  implicit none

  real(kind=r_tran), intent(in)     :: tracer1,tracer2,tracer3,dz1,dz2,dz3
  real(kind=r_tran), intent(out)    :: tracer1i,tracer2i,tracer3i
  real(kind=r_tran)                 :: h,b1,b2,b3,d1,d2,d3,fmin,fmax
  real(kind=r_tran), dimension(3,3) :: c

  h = min(dz1,min(dz2,dz3))
  fmin = min(tracer1,min(tracer2,tracer3))
  fmax = max(tracer1,max(tracer2,tracer3))

  b1 = dz1 + dz2
  b2 = dz2 + dz3
  b3 = 2.0_r_tran*dz2 + dz1

  d1 = b1*(b1+b2)
  d2 = b1*b2
  d3 = (b1+b2)*b2

  c(1,1) = (-3.0_r_tran*h + dz2)*(-3.0_r_tran*h - dz3)/d1
  c(1,2) = (-3.0_r_tran*h + b3 )*( 3.0_r_tran*h + dz3)/d2
  c(1,3) = (-3.0_r_tran*h + b3 )*(-3.0_r_tran*h + dz2)/d3

  c(2,1) = (-h + dz2)*(-h - dz3)/d1
  c(2,2) = (-h + b3 )*( h + dz3)/d2
  c(2,3) = (-h + b3 )*(-h + dz2)/d3

  c(3,1) = (h + dz2)*( h - dz3)/d1
  c(3,2) = (h + b3 )*(-h + dz3)/d2
  c(3,3) = (h + b3 )*( h + dz2)/d3

  tracer1i = c(1,1)*tracer1 + c(1,2)*tracer2 + c(1,3)*tracer3
  tracer2i = c(2,1)*tracer1 + c(2,2)*tracer2 + c(2,3)*tracer3
  tracer3i = c(3,1)*tracer1 + c(3,2)*tracer2 + c(3,3)*tracer3

  tracer1i = max(fmin,min(tracer1i,fmax))
  tracer2i = max(fmin,min(tracer2i,fmax))
  tracer3i = max(fmin,min(tracer3i,fmax))

end subroutine interpolate_to_regular_grid

end module koren_support_mod
