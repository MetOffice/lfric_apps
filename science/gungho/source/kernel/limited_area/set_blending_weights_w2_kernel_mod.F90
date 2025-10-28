!-----------------------------------------------------------------------------
! (C) Crown copyright 2022 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Populates the w2 weights field with the blending weights obtained from the
!>        configuration.
!> @details Uses the onion layer values as array index for blending_weights array
!>          to populate a field of blending weights.
!>          The onion layers are on the discontinuous W3 space.  For the W2 dofs,
!>          average values are taken using the weights from either side of the face.
module set_blending_weights_w2_kernel_mod

  use argument_mod,              only : arg_type,             &
                                        GH_SCALAR, GH_FIELD,  &
                                        GH_READ, GH_REAL,     &
                                        GH_INC,               &
                                        GH_INTEGER, GH_BASIS, &
                                        CELL_COLUMN,          &
                                        STENCIL, CROSS2D
  use fs_continuity_mod,         only : W3, W2
  use constants_mod,             only : r_def, i_def
  use kernel_mod,                only : kernel_type
  use reference_element_mod,     only : T, B

  implicit none

  private

  !-------------------------------------------------------------------------
  ! Public types
  !-------------------------------------------------------------------------

  type, public, extends(kernel_type) :: set_blending_weights_w2_kernel_type
    private
    type(arg_type) :: meta_args(3) = (/                                &
         arg_type(GH_FIELD,   GH_REAL, GH_INC, W2),                    & ! weights_field
         arg_type(GH_FIELD,   GH_REAL, GH_READ, W3, STENCIL(CROSS2D)), & ! onion_layers
         arg_type(GH_SCALAR,  GH_INTEGER, GH_READ)                     & ! depth
         /)
    integer :: operates_on = CELL_COLUMN
  contains
    procedure, nopass :: set_blending_weights_w2_code
  end type

  !-------------------------------------------------------------------------
  ! Contained functions/subroutines
  !-------------------------------------------------------------------------
  public :: set_blending_weights_w2_code

contains

!> @param[in] nlayers      Number of layers
!> @param[in,out] weights_field The output weights field
!> @param[in] onion_layers   The onion layer field
!> @param[in] stencil_sizes  Sizes of the branches of the cross stencil
!> @param[in] stencil_max    Maximum branch size of the cross stencil
!> @param[in] stencil_map    The stencil map
!> @param[in] depth          Depth of the weight array
!> @param[in] ndf_out        Number of degrees of freedom for weights_field
!> @param[in] undf_out       Total number of degrees of freedom for weights_field
!> @param[in] map_out        Dofmap for the cell at the base of the column for weights_field
!> @param[in] ndf_in         Number of degrees of freedom for onion_layers
!> @param[in] undf_in        Total number of degrees of freedom for onion_layers
!> @param[in] map_in         Dofmap for the cell at the base of the column for onion_layers
subroutine set_blending_weights_w2_code( nlayers,        &
                                         weights_field,  &
                                         onion_layers,   &
                                         stencil_sizes,  &
                                         stencil_max,    &
                                         stencil_map,    &
                                         depth,          &
                                         ndf_out,        &
                                         undf_out,       &
                                         map_out,        &
                                         ndf_in,         &
                                         undf_in,        &
                                         map_in)

  use boundaries_config_mod,        only : blending_weights, &
                                           blending_weights_w2v

  implicit none

  ! Arguments
  integer(kind=i_def),  intent(in)    :: nlayers
  integer(kind=i_def),  intent(in)    :: ndf_out, undf_out
  integer(kind=i_def),  intent(in)    :: ndf_in, undf_in
  integer(kind=i_def),  intent(in)    :: stencil_max
  real(kind=r_def),     intent(inout) :: weights_field(undf_out)
  real(kind=r_def),     intent(in)    :: onion_layers(undf_in)
  integer(kind=i_def),  intent(in)    :: stencil_sizes(4)
  integer(kind=i_def),  intent(in)    :: stencil_map(ndf_in,stencil_max,4)
  integer(kind=i_def),  intent(in)    :: depth
  integer(kind=i_def),  intent(in)    :: map_out(ndf_out)
  integer(kind=i_def),  intent(in)    :: map_in(ndf_in)

  ! Internal variables
  integer(kind=i_def) :: k, df
  integer(kind=i_def) :: index
  integer(kind=i_def) :: onion_layer

  onion_layer = int(onion_layers(map_in(1)), i_def)

  ! W2 weights should not go right up to the inner region
  if (onion_layer > 0_i_def)then
    index = depth - onion_layer + 1

    ! Vertical dofs first - these are all set to the W2V blending weights
    do k = 0, nlayers-1
      do df = B, T
        weights_field(map_out(df)+k) = blending_weights_w2v(index)
      end do
    end do

    ! Next the horizontal dofs - using W2H blending weights
    do k = 0, nlayers-1
      ! Loop over four sides of cross stencil (for quadratileral cell)
      do df = 1, 4
        ! At very edge of domain, set weight to be 1
        if (stencil_sizes(df) == 1_i_def) then
          weights_field(map_out(df)+k) = 1.0_r_def
        else
          ! Otherwise, take the average weight from neighbouring columns.
          ! NB: in the domain interior, the blending weight is 0.
          ! Therefore we can compute the W2 blending weights by adding
          ! half the weight for each column
          weights_field(map_out(df)+k) = min( 1.0_r_def,                       &
              weights_field(map_out(df)+k)                                     &
              + 0.5_r_def*blending_weights(index)                              &
          )
        end if
      end do
    end do
  end if

end subroutine set_blending_weights_w2_code

end module set_blending_weights_w2_kernel_mod
