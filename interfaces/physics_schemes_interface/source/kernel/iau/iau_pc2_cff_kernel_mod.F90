!-----------------------------------------------------------------------------
! (c) Crown copyright 2023 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
!> @brief Calculate changes to ice cloud fraction

module iau_pc2_cff_kernel_mod
  use argument_mod,      only : arg_type,     &
                                GH_FIELD,     &
                                GH_REAL,      &
                                GH_READWRITE, &
                                CELL_COLUMN
  use constants_mod,     only : r_def, i_def
  use fs_continuity_mod, only : WTHETA
  use kernel_mod,        only : kernel_type

  implicit none

  private

!-------------------------------------------------------------------------------
! Public types
!-------------------------------------------------------------------------------
!> The type declaration for the kernel. Contains the metadata needed by the Psy layer
type, public, extends(kernel_type) :: iau_pc2_cff_kernel_type
  private
  type(arg_type) :: meta_args(5) = (/                  &
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! qcf
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! del_qcf
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! bcf
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA), & ! cfl
        arg_type(GH_FIELD,  GH_REAL, GH_READWRITE,  WTHETA)  & ! cff
       /)

  integer :: operates_on = CELL_COLUMN
contains
  procedure, nopass :: iau_pc2_cff_code
end type
!-------------------------------------------------------------------------
! Contained functions/subroutines
!-------------------------------------------------------------------------
    public :: iau_pc2_cff_code

contains

!> @brief Calculate change in total cloud fraction due to a change in liquid
!> cloud fraction. This depends upon the sign of the change of liquid
!> cloud fraction.
!> @param[in]        nlayers      number of layers
!> @param[in]        qcf          ice content
!> @param[in]        del_qcf      forcing of ice
!> @param[in,out]    bcf          bulk cloud fraction
!> @param[in,out]    cff          frozen cloud fraction
!> @param[in,out]    cfl          liquid cloud fraction
!> @param[in]        ndf_wth      the number of degrees of freedom for a field in wtheta space.
!> @param[in]        undf_wth     the number of unique degrees of freedom for a field in wtheta
!> @param[in]        map_wth      dofmap for the cell at the base of the column

subroutine iau_pc2_cff_code( nlayers,                           &
                             qcf,del_qcf,                       &
                             bcf,cfl,cff,                       &
                             ndf_wth, undf_wth, map_wth         &
                           )

  implicit none

      ! Arguments
      integer(kind=i_def), intent(in)                         :: nlayers
      integer(kind=i_def), intent(in)                         :: ndf_wth
      integer(kind=i_def), intent(in)                         :: undf_wth
      real(kind=r_def), intent(in),       dimension(undf_wth) :: qcf
      real(kind=r_def), intent(in),       dimension(undf_wth) :: del_qcf
      real(kind=r_def), intent(inout),    dimension(undf_wth) :: bcf
      real(kind=r_def), intent(inout),    dimension(undf_wth) :: cff
      real(kind=r_def), intent(inout),    dimension(undf_wth) :: cfl
      integer(kind=i_def), intent(in),    dimension(ndf_wth)  :: map_wth

      ! Internal variables
      integer(kind=i_def)                                 :: k
      real(kind=r_def)                                    :: del_cff
      real(kind=r_def)                                    :: bcf_correction
      real(kind=r_def), parameter                         :: qcf0 = 1.0e-04_r_def

      do k=0,nlayers

        if (del_qcf( map_wth(1) + k ) /= 0.0_r_def) then

          ! Adjust del_cff
          del_cff = del_qcf( map_wth(1) + k ) / &
                  ( qcf( map_wth(1) + k ) + (1.0_r_def - cff( map_wth(1)+k ) ) * qcf0 )
          ! Check that del_cff is sensible
          del_cff = max( min( del_cff, ( 1.0_r_def-cff( map_wth(1)+k ) ) ),   &
                        -cff( map_wth(1)+k )  )
          if ( del_cff > 0.0_r_def .and. cff(map_wth(1) + k) < 1.0_r_def ) then
            bcf_correction = min( del_cff , ( 1.0_r_def - bcf(map_wth(1) + k ) ) )

          else if ( del_cff < 0.0_r_def .and. cff(map_wth(1) + k) > 0.0_r_def ) then
            bcf_correction = max( del_cff, ( cfl( map_wth(1)+ k ) - &
                             bcf( map_wth(1) + k ) ) )
          else
            bcf_correction = 0.0_r_def
          end if

          ! Update bcf
          bcf( map_wth(1) + k )   =   bcf( map_wth(1) + k ) + bcf_correction
          ! Check that the bulk cloud fraction is constrained within 0 and 1
          bcf( map_wth(1) + k ) = max( min( bcf( map_wth(1) + k ), 1.0_r_def), &
          0.0_r_def)

          ! Update cff
          cff( map_wth(1) + k ) = cff ( map_wth(1) + k ) + del_cff
          cff( map_wth(1) + k ) = max(min(cff(map_wth(1) + k), 1.0_r_def), &
          0.0_r_def)
          ! Check that the liquid cloud fraction is constrained within 0 and 1
          cfl( map_wth(1) + k ) = max(min(cfl(map_wth(1) + k), 1.0_r_def), &
          0.0_r_def)

        end if

      end do

end subroutine iau_pc2_cff_code

end module iau_pc2_cff_kernel_mod
