! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Subroutine to do safety / consistency checks on the prognostic
! precipitation fraction.
! This needs to be called at the end of the model timetep, to sort out
! instances where the dynamics has caused precip mass to go non-zero
! at grid-points where precip fraction is currently zero.

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Large_Scale_Precipitation

module lsp_precfrac_checks_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_PRECFRAC_CHECKS_MOD'

contains

subroutine lsp_precfrac_checks( dims, p_rho_levels,                            &
                                q, qcl, qcf, qcf2, qrain, qgraup,              &
                                precfrac )

use atm_fields_bounds_mod, only: array_dims, tdims, pdims_s
use um_types,              only: real_umphys
use mphys_inputs_mod,      only: l_mcr_qcf2, l_mcr_qrain, l_mcr_qgraup,        &
                                 l_subgrid_graupel_frac,                       &
                                 l_improve_precfrac_checks
use pc2_constants_mod,     only: max_in_cloud_qcf
use yomhook,               only: lhook, dr_hook
use parkind1,              only: jprb, jpim

implicit none

! Dimensions of the input arrays (this is passed through the argument
! list to allow this routine to be called on versions of the fields
! with or without processor halos; pass in tdims_l for dims if they have
! extended halos, or just tdims if they have no halos, etc).
! Note that we don't update any halos present, so we only loop over
! tdims in the calculations, regardless of dims.
type(array_dims), intent(in) :: dims

! Pressure on rho-levels, used to compute vertical mean q_total
real(kind=real_umphys), intent(in) :: p_rho_levels                             &
                                      ( pdims_s%i_start:pdims_s%i_end,         &
                                        pdims_s%j_start:pdims_s%j_end,         &
                                        pdims_s%k_start:pdims_s%k_end )

! Water-vapour mass, used to calculate limit on in-fraction precip-mass
real(kind=real_umphys), intent(in) :: q      ( dims%i_start:dims%i_end,        &
                                               dims%j_start:dims%j_end,        &
                                               dims%k_start:dims%k_end )

! Condensed water species masses
real(kind=real_umphys), intent(in) :: qcl    ( dims%i_start:dims%i_end,        &
                                               dims%j_start:dims%j_end,        &
                                               dims%k_start:dims%k_end )
real(kind=real_umphys), intent(in) :: qcf    ( dims%i_start:dims%i_end,        &
                                               dims%j_start:dims%j_end,        &
                                               dims%k_start:dims%k_end )
real(kind=real_umphys), intent(in) :: qcf2   ( dims%i_start:dims%i_end,        &
                                               dims%j_start:dims%j_end,        &
                                               dims%k_start:dims%k_end )
real(kind=real_umphys), intent(in) :: qrain  ( dims%i_start:dims%i_end,        &
                                               dims%j_start:dims%j_end,        &
                                               dims%k_start:dims%k_end )
real(kind=real_umphys), intent(in) :: qgraup ( dims%i_start:dims%i_end,        &
                                               dims%j_start:dims%j_end,        &
                                               dims%k_start:dims%k_end )

! Prognostic precipitation fraction to be checked / updated
real(kind=real_umphys), intent(in out) :: precfrac ( dims%i_start:dims%i_end,  &
                                                     dims%j_start:dims%j_end,  &
                                                     dims%k_start:dims%k_end )

! Max allowed in-precip-fraction condensate mixing-ratio divided by
! precip-fraction
! (max in-fraction condensate = qp_max * PF)
! Used to compute a corresponding minimum allowed precip-fraction
! imposed as a safety check to avoid advection creating instances of
! near-zero fraction but non-zero mixing-ratio
real(kind=real_umphys) :: qp_max ( tdims%i_start:tdims%i_end,                  &
                                   tdims%j_start:tdims%j_end,                  &
                                               1:tdims%k_end )

! Dimensionless factor scaling the max allowed in-precip-fraction condensate
! We impose  qp/pf < qp_max_fac * q_tot * pf
real(kind=real_umphys), parameter :: qp_max_fac = 1000.0_real_umphys
              ! gives 5 g kg-1 when CF = 0.0005 and q_tot = 10 g kg-1

! Smallest non-zero number, used to avoid div-by-zero
real(kind=real_umphys), parameter :: min_float = tiny(qp_max)

! Zero and one stored in native precision
real(kind=real_umphys), parameter :: zero = 0.0_real_umphys
real(kind=real_umphys), parameter :: one  = 1.0_real_umphys

! Loop counters
integer :: i, j, k

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_PRECFRAC_CHECKS'


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)


!$OMP PARALLEL DEFAULT(NONE) PRIVATE( i, j, k )                                &
!$OMP SHARED( l_mcr_qrain, l_mcr_qcf2, l_mcr_qgraup, tdims, qp_max,            &
!$OMP         q, qcl, qcf, qcf2, qrain, qgraup, p_rho_levels,                  &
!$OMP         l_subgrid_graupel_frac, l_improve_precfrac_checks, precfrac )

if ( l_improve_precfrac_checks ) then
  ! Improved version of precip fraction checking; uses variable max
  ! in-rainshaft precip-mass (scales with mean total-water in the column
  ! and with the precip-fraction) instead of a fixed constant.
  ! Also adds checks to ignore negative values on precip-mass, and removes
  ! instances of non-zero precip-fraction where precip-mass has gone to zero.

  ! Compute max in-precip-fraction water content
  ! (pending scaling by precip-fraction).
  ! Make this scale with the grid-mean total-water content
!$OMP DO SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    ! Sum condensate species that are always used (ignoring negative values)
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        qp_max(i,j,k) = max( q(i,j,k),   zero )                                &
                      + max( qcl(i,j,k), zero )                                &
                      + max( qcf(i,j,k), zero )
      end do
    end do
    ! Add on the 3 optional species if used
    if ( l_mcr_qcf2 ) then
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          qp_max(i,j,k) = qp_max(i,j,k) + max( qcf2(i,j,k),   zero )
        end do
      end do
    end if
    if ( l_mcr_qrain ) then
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          qp_max(i,j,k) = qp_max(i,j,k) + max( qrain(i,j,k),  zero )
        end do
      end do
    end if
    if ( l_mcr_qgraup ) then
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          qp_max(i,j,k) = qp_max(i,j,k) + max( qgraup(i,j,k), zero )
        end do
      end do
    end if
    ! Apply dimensionless scaling factor and tiny limit to avoid div-by-zero
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        qp_max(i,j,k) = max( qp_max(i,j,k) * qp_max_fac, min_float )
      end do
    end do
  end do  ! k = 1, tdims%k_end
!$OMP END DO

  ! Setting qp_max to scale with q_total seems to unduly over-restrict
  ! in-precip-fraction graupel at upper-levels, where typical ratios of
  ! qp / q_vap are bigger (especially where deep convection is detraining high
  ! concentrations of graupel into otherwise very dry air).
  ! This check is spuriously inflating precfrac where there
  ! is convectively-generated graupel.
  ! To avoid this, replace qp_max with its vertical mean over the layers
  ! below (just a way of not making it so much smaller higher-up,
  ! without introducing additional ad-hoc thresholds).
!$OMP DO SCHEDULE(STATIC)
  do j = tdims%j_start, tdims%j_end
    do k = 2, tdims%k_end-1
      do i = tdims%i_start, tdims%i_end
        qp_max(i,j,k) = ( qp_max(i,j,k-1)                                      &
                          * ( p_rho_levels(i,j,1) - p_rho_levels(i,j,k) )      &
                        + qp_max(i,j,k)                                        &
                          * ( p_rho_levels(i,j,k) - p_rho_levels(i,j,k+1) )    &
                        ) / ( p_rho_levels(i,j,1) - p_rho_levels(i,j,k+1) )
      end do
    end do
    do i = tdims%i_start, tdims%i_end
      qp_max(i,j,tdims%k_end) = qp_max(i,j,tdims%k_end-1)
    end do
  end do
!$OMP END DO

  ! Impose a minimum limit on the precip-fraction so-as to enforce
  ! a consistent maximum limit on the in-fraction precip water content.
  ! We make the max in-fraction water content scale with precip-fraction,
  ! so-as to mainly affect "noise" present in very small precip-fractions,
  ! while leaving genuine occurences of high water contents alone.
  ! We impose:
  !  qp/pf < qp_max * pf
  ! => pf^2 > qp / qp_max
  ! The checks don't make sense if the condensed water masses
  ! have gone negative, so we treat negative values as zeros here.
  ! Note: this check also removes instances of
  ! frac < 0 or frac > 1 which SL advection can create.

  if ( l_mcr_qgraup .and. l_subgrid_graupel_frac ) then
    ! Graupel included in the prognostic precip fraction as well as rain
!$OMP DO SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          precfrac(i,j,k) = min( max( precfrac(i,j,k),                         &
                                      sqrt( ( max( qrain(i,j,k), zero )        &
                                            + max( qgraup(i,j,k), zero ) )     &
                                          / qp_max(i,j,k) )                    &
                                    ), one )
          if ( .not. ( qrain(i,j,k) > zero .or. qgraup(i,j,k) > zero ) ) then
            ! Reset to zero if no precip mass
            precfrac(i,j,k) = zero
          end if
        end do
      end do
    end do
!$OMP END DO NOWAIT
  else
    ! Only rain included in the prognostic precip fraction
!$OMP DO SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          precfrac(i,j,k) = min( max( precfrac(i,j,k),                         &
                                      sqrt( max( qrain(i,j,k), zero )          &
                                          / qp_max(i,j,k) )                    &
                                    ), one )
          if ( .not. ( qrain(i,j,k) > zero ) ) then
            ! Reset to zero if no precip mass
            precfrac(i,j,k) = zero
          end if
        end do
      end do
    end do
!$OMP END DO NOWAIT
  end if

else  ! ( .NOT. l_improve_precfrac_checks )
  ! Original version of the checks; just imposes a fixed max limit on
  ! in-rainshaft precip-mass (taken from pc2_constants qcf limit),
  ! and keeps precfrac between 0 and 1

  if ( l_subgrid_graupel_frac ) then
    ! Graupel included in the "precip" fraction
!$OMP DO SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          precfrac(i,j,k) = max( precfrac(i,j,k),                              &
                                 ( qrain(i,j,k)                                &
                                 + qgraup(i,j,k) ) / max_in_cloud_qcf )
        end do
      end do
    end do
!$OMP END DO
  else
    ! "Precip" fraction only includes rain
!$OMP DO SCHEDULE(STATIC)
    do k = 1, tdims%k_end
      do j = tdims%j_start, tdims%j_end
        do i = tdims%i_start, tdims%i_end
          precfrac(i,j,k) = max( precfrac(i,j,k),                              &
                                 qrain(i,j,k) / max_in_cloud_qcf )
        end do
      end do
    end do
!$OMP END DO
  end if

  ! Ensure fraction is between 0 and 1
!$OMP DO SCHEDULE(STATIC)
  do k = 1, tdims%k_end
    do j = tdims%j_start, tdims%j_end
      do i = tdims%i_start, tdims%i_end
        precfrac(i,j,k) = max( min( precfrac(i,j,k), one ), zero )
      end do
    end do
  end do
!$OMP END DO

end if  ! ( .NOT. l_improve_precfrac_checks )

!$OMP END PARALLEL


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_precfrac_checks

end module lsp_precfrac_checks_mod
