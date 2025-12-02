! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module entdet_res_source_mod

implicit none

contains


! Subroutine to calculate the resolved-scale source or sink of
! each primary variable due to entrainment or detrainment
subroutine entdet_res_source( n_points, n_points_entdet_super,                 &
                              n_points_res_super, n_fields_tot, l_ent,         &
                              entdet_mass_d, entdet_fields,                    &
                              res_source_super, res_source_fields )

use comorph_constants_mod, only: real_cvprec
use res_source_type_mod, only: n_res, i_ent, i_det

implicit none

! Array dimensions:

! Number of points
integer, intent(in) :: n_points
! Number of points in the fields and res_source super-arrays
! (to avoid having to repeatedly deallocate and reallocate when
!  re-using the same array, it maybe bigger than needed here).
integer, intent(in) :: n_points_entdet_super
integer, intent(in) :: n_points_res_super
! Number of fields in the super-array
integer, intent(in) :: n_fields_tot

! Flag for doing calculations for entrainment versus detrainment
logical, intent(in) :: l_ent

! Rate of entrainment or detrainment of dry-mass per s per
! unit surface area.
real(kind=real_cvprec), intent(in) :: entdet_mass_d(n_points)

! Super-array storing primary field values for the
! entrained or detrained air, expressed as variables
! which are conserved following dry-mass at constant pressure
real(kind=real_cvprec), intent(in) :: entdet_fields                            &
                                      ( n_points_entdet_super, n_fields_tot )

! Super-array storing sum of entrainment and detrainment mass sources/sinks
real(kind=real_cvprec), intent(in out) :: res_source_super                     &
                                          ( n_points_res_super, n_res )

! Super-array storing resolved-scale source terms for primary fields
real(kind=real_cvprec), intent(in out) :: res_source_fields                    &
                                          ( n_points_res_super, n_fields_tot )

! Loop counters
integer :: ic, i_field


! Subtract sink terms from the environment due to entrainment,
! or add source terms to the environment due to detrainment.

if ( l_ent ) then
  ! If doing entrainment...

  ! Add to total entrainment
  do ic = 1, n_points
    res_source_super(ic,i_ent) = res_source_super(ic,i_ent) + entdet_mass_d(ic)
  end do

  ! Decrement all the fields in the super-array
  do i_field = 1, n_fields_tot
    do ic = 1, n_points
      res_source_fields(ic,i_field) = res_source_fields(ic,i_field)            &
        - entdet_mass_d(ic) * entdet_fields(ic,i_field)
    end do
  end do

else
  ! If doing detrainment...

  ! Add to total detrainment
  do ic = 1, n_points
    res_source_super(ic,i_det) = res_source_super(ic,i_det) + entdet_mass_d(ic)
  end do

  ! Increment all the fields in the super-array
  do i_field = 1, n_fields_tot
    do ic = 1, n_points
      res_source_fields(ic,i_field) = res_source_fields(ic,i_field)            &
        + entdet_mass_d(ic) * entdet_fields(ic,i_field)
    end do
  end do

end if  ! ( l_ent )


return
end subroutine entdet_res_source


end module entdet_res_source_mod
