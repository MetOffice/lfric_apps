! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: convection_comorph

module calc_kqkt_mod

implicit none

contains


! Subroutine to calculate the coefficients for exchange of water
! vapour and heat between hydrometeors of a given species and the
! surrounding parcel air
subroutine calc_kqkt( n_points, area_coef,                                     &
                      n_cond, r_cond, wf_cond, rho_dry,                        &
                      kq_cond, kt_cond )

use comorph_constants_mod, only: pi, mol_diff_q, mol_diff_t, vent_factor,      &
                     real_cvprec, four

implicit none

! Number of points
integer, intent(in) :: n_points

! Coefficient scaling the surface area of hydrometeors,
! to account for non-spheres.
real(kind=real_cvprec), intent(in) :: area_coef

! Number concentration per unit dry-mass
real(kind=real_cvprec), intent(in) :: n_cond(n_points)
! Radius of the hydrometeor particles
real(kind=real_cvprec), intent(in) :: r_cond(n_points)
! Fall-speed of the hydrometeor particles
real(kind=real_cvprec), intent(in) :: wf_cond(n_points)
! Dry density of parcel air
real(kind=real_cvprec), intent(in) :: rho_dry(n_points)

! Coefficient for exchange of water vapour between parcel vapour
! and hydrometeor condensate / s-1
real(kind=real_cvprec), intent(out) :: kq_cond(n_points)
! Coefficient for exchange of heat between parcel air and the
! hydrometeors / s-1
real(kind=real_cvprec), intent(out) :: kt_cond(n_points)

real(kind=real_cvprec) :: work(n_points)

! Loop counter
integer :: ic

! Exchange coefficients for moisture and heat have almost
! identical formulae...

! First calculate ventilation term, which is the same for both:
do ic = 1, n_points
  work(ic) = vent_factor * r_cond(ic) * wf_cond(ic)
end do

! Add on the molecular diffusivity term,
! which uses different constants for moisture and heat
do ic = 1, n_points
  kq_cond(ic) = four*mol_diff_q + work(ic)
  kt_cond(ic) = four*mol_diff_t + work(ic)
end do

! Calculate remaining multiplying term which is the same for both
do ic = 1, n_points
  work(ic) = area_coef * n_cond(ic) * rho_dry(ic)                              &
             * pi * r_cond(ic)
end do

! Scale the k's to get final answers
do ic = 1, n_points
  kq_cond(ic) = kq_cond(ic) * work(ic)
  kt_cond(ic) = kt_cond(ic) * work(ic)
end do

return
end subroutine calc_kqkt


end module calc_kqkt_mod
