!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief   Increments an onion_layer mask resulting in a deeper mask.
!> @details Increments by 1 the values in a mask where there are non-zero values
!>          surrounding a cell.  This has the effect of increasing the depth
!>          of the mask and generating an additional 'onion-layer'. The layer
!>          numbering then increases from inner to outer layers.
!>          Requires stencil with depth 1.
module propagate_onion_layers_kernel_mod

  use argument_mod,              only : arg_type, GH_FIELD,   &
                                        GH_READ, GH_WRITE,    &
                                        GH_REAL, CELL_COLUMN, &
                                        STENCIL, CROSS,       &
                                        ANY_DISCONTINUOUS_SPACE_1
  use constants_mod,             only : r_def, i_def
  use kernel_mod,                only : kernel_type

  implicit none

  private

  !-------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------

  type, public, extends(kernel_type) :: propagate_onion_layers_kernel_type
    private
    type(arg_type) :: meta_args(2) = (/                &
         arg_type(GH_FIELD,   GH_REAL, GH_WRITE, ANY_DISCONTINUOUS_SPACE_1),               & ! propagated_onion_layers
         arg_type(GH_FIELD,   GH_REAL, GH_READ, ANY_DISCONTINUOUS_SPACE_1, STENCIL(CROSS)) & ! onion_layers
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: propagate_onion_layers_code
  end type

  !-------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------
  public :: propagate_onion_layers_code

contains

!> @param[in] nlayers           Number of layers
!> @param[in,out] propagated_onion_layers  Updated onion_layers
!> @param[in] onion_layers      Input Onion_layers
!> @param[in] stencil_size      Local stencil size (will be reduced near a boundary)
!> @param[in] stencil_map       Stencil map
!> @param[in] ndf               Number of degrees of freedom per cell
!>                              for ANY_DISCONTINUOUS_SPACE_1
!> @param[in] undf              Total number of degrees of freedom for
!>                              ANY_DISCONTINUOUS_SPACE_1
!> @param[in] map               Dofmap for the cell at the base of the
!>                              column for ANY_DISCONTINUOUS_SPACE_1
subroutine propagate_onion_layers_code( nlayers,                 &
                                        propagated_onion_layers, &
                                        onion_layers,            &
                                        stencil_size,            &
                                        stencil_map,             &
                                        ndf,                     &
                                        undf,                    &
                                        map                      &
                                        )

  implicit none

  ! Arguments
  integer(kind=i_def),                  intent(in) :: nlayers, &
                                                      ndf,     &
                                                      undf
  real(kind=r_def), dimension(undf), intent(inout) :: propagated_onion_layers
  real(kind=r_def), dimension(undf), intent(in)    :: onion_layers
  integer(kind=i_def), intent(in)                  :: stencil_size
  integer(kind=i_def), dimension(ndf,stencil_size), intent(in) :: stencil_map
  integer(kind=i_def), dimension(ndf), intent(in)  :: map
  ! Internal variables
  integer(kind=i_def) :: k, df, i
  real(kind=r_def) :: max_of_surrounding_cells

  max_of_surrounding_cells=0.0_r_def
  do i=2,stencil_size
    max_of_surrounding_cells =  &
       max(max_of_surrounding_cells, onion_layers(stencil_map(1,i)))
  end do
  if (max_of_surrounding_cells > 0)then
    do k=0,nlayers-1
      do df=1,ndf
        propagated_onion_layers(stencil_map(df,1)+k) = &
           onion_layers(stencil_map(df,1)+k) + 1.0_r_def
      end do
    end do
  end if


end subroutine propagate_onion_layers_code

end module propagate_onion_layers_kernel_mod
