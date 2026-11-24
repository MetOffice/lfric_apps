!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Computes a mask with 1s on the cells on the domain boundaries,
!>          all other cells marked as inactive with 0s
!> @details Makes use of the stencil extents of the W3 field to determine
!>          if the column is on the edge of the domain. This
!>          creates a mask that has 0s everywhere, except for the cells that are
!>          directly next to the domain boundary. i.e. The cell depth of the rim is 1
!>          Requires stencil with depth 1.
module create_edge_mask_kernel_mod

  use argument_mod,              only : arg_type, func_type, &
                                        GH_SCALAR, GH_FIELD, &
                                        GH_READ, GH_WRITE,   &
                                        GH_REAL,             &
                                        CELL_COLUMN,         &
                                        STENCIL, CROSS,  GH_INTEGER, &
                                        ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,             only : r_def, i_def
  use kernel_mod,                only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------

  type, public, extends(kernel_type) :: create_edge_mask_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                       &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE,              & ! edge_mask
                  ANY_DISCONTINUOUS_SPACE_1),                 &
         arg_type(GH_FIELD,   GH_REAL, GH_READ,               & ! mask
                  ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS)), &
         arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)            & ! stencil_depth
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: create_edge_mask_code
  end type

  !-------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------
  public :: create_edge_mask_code

contains

!> @param[in] nlayers      Number of layers
!> @param[in,out] edge_mask The mask to define the domain boundary cells
!> @param[in] mask         Input mask field - this is essentially a dummy
!>                         field to read from and need only have the same
!>                         stencil as the edge_mask
!> @param[in] stencil_size Local stencil size (will be reduced near a boundary)
!> @param[in] stencil_map Stencil map
!> @param[in] stencil_depth Stencil depth
!> @param[in] ndf       Number of degrees of freedom per cell for ANY_DISCONTINUOUS_SPACE_1
!> @param[in] undf      Total number of degrees of freedom for ANY_DISCONTINUOUS_SPACE_1
!> @param[in] map       Dofmap for the cell at the base of the column for ANY_DISCONTINUOUS_SPACE_1
subroutine create_edge_mask_code( nlayers,       &
                                  edge_mask,     &
                                  mask,          &
                                  stencil_size,  &
                                  stencil_map,   &
                                  stencil_depth, &
                                  ndf,           &
                                  undf,          &
                                  map)

  implicit none

  ! Arguments
  integer(kind=i_def),                  intent(in) :: nlayers,     &
                                                      ndf,         &
                                                      undf
  real(kind=r_def), dimension(undf), intent(inout) :: edge_mask
  real(kind=r_def), dimension(undf), intent(in)    :: mask
  integer(kind=i_def), intent(in)                  :: stencil_size
  integer(kind=i_def), intent(in)                  :: stencil_depth
  integer(kind=i_def), dimension(ndf,stencil_size), intent(in) :: stencil_map
  integer(kind=i_def), dimension(ndf),  intent(in) :: map

  ! Internal variables
  integer(kind=i_def)                    :: k, df

  ! The usual size of a cross stencil is 4*stencil_depth+1, but is less than
  ! this near the edges/corners of the domain.
  if (stencil_size < 4*stencil_depth+1)then
    do k=0,nlayers-1
      do df=1,ndf
        edge_mask(stencil_map(df,1)+k) = 1.0_r_def
      end do
    end do
  end if

end subroutine create_edge_mask_code

end module create_edge_mask_kernel_mod
