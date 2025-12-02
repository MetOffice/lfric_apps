! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module entrain_fields_mod

implicit none

contains

! Subroutine to perform entrainment of the environment
! primary fields into the parcel
subroutine entrain_fields( n_points, n_points_par_super,                       &
                           n_points_ent_super, n_fields_tot,                   &
                           ent_mass_d, massflux_d,                             &
                           ent_fields, par_fields )

use comorph_constants_mod, only: real_cvprec, one, min_float

implicit none

! Number of points
integer, intent(in) :: n_points

! Number of points in the compressed super-arrays for the
! parcel and the entrained air
! (to avoid repeatedly deallocating and reallocating these,
!  they are dimensioned with the biggest size they will need,
!  which will often be bigger than the number of points here)
integer, intent(in) :: n_points_par_super
integer, intent(in) :: n_points_ent_super

! Number of fields in the super-arrays
integer, intent(in) :: n_fields_tot

! Dry-mass entrained from level k
real(kind=real_cvprec), intent(in) :: ent_mass_d(n_points)

! Dry-mass flux
! (assumed to already have the entrained mass added on)
real(kind=real_cvprec), intent(in) :: massflux_d(n_points)

! Super-arrays containing the entrained air properties
! and the parcel properties to be updated
real(kind=real_cvprec), intent(in) :: ent_fields                               &
                         ( n_points_ent_super, n_fields_tot )
real(kind=real_cvprec), intent(in out) :: par_fields                           &
                         ( n_points_par_super, n_fields_tot )


! Weight aplied to entrained air in combined parcel calculation
real(kind=real_cvprec) :: weight(n_points)

! Loop counters
integer :: ic, i_field

! Convert the entrained mass to a weight for the contribution
! of the entrained properties to the combined parcel properties
do ic = 1, n_points
  ! massflux_d already includes the entrained mass, so the fraction
  ! of the mass present that came from the entrained air is:
  weight(ic) = ent_mass_d(ic) / max( massflux_d(ic), min_float )
end do

! Loop over all the fields in the super-array and compute
! new parcel properties after entrainment; this is a weighted mean
! of the existing parcel and the entrained air.
do i_field = 1, n_fields_tot
  do ic = 1, n_points
    par_fields(ic,i_field) = (one-weight(ic)) * par_fields(ic,i_field)         &
                           +      weight(ic)  * ent_fields(ic,i_field)
  end do
end do

return
end subroutine entrain_fields


end module entrain_fields_mod

