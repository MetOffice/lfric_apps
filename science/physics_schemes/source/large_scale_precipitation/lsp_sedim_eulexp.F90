! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Microphysics hydrometeor Eulerian sedimentation scheme
module lsp_sedim_eulexp_mod

implicit none

character(len=*), parameter, private :: ModuleName='LSP_SEDIM_EULEXP_MOD'

contains

subroutine lsp_sedim_eulexp(                                                   &
  points,m0,dhi,dhir,rho,rhor,                                                 &
  flux_fromabove, fallspeed_thislayer,                                         &
  mixratio_thislayer, fallspeed_fromabove,                                     &
  total_flux_out)

use lsprec_mod, only: zero, one

! Use in kind for large scale precip, used for compressed variables passed down
! from here
use um_types,             only: real_lsprec

use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim
implicit none

! Description:
!   First order Eulerian advection scheme for hydrometeor sedimentation
!   with an exponential-based limiter to ensure stability for CFL>1
!   used by routine LSP_ICE for ice crystals, snow, rain and graupel.

! Method:
!   Based on method described in Rotstayn (1997)(QJRMS, 123, 1227-1282)
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

! Code Description:
!   Language: FORTRAN 77 + common extensions.
!   This code is written to UMDP3 v6 programming standards.


! Subroutine arguments

    ! Intent (In)
integer, intent(in) :: points       ! number of points to process

    ! Intent (In)
real (kind=real_lsprec), intent(in) ::                                         &
  m0,                                                                          &
                         ! Small mass (kg/kg) defined in c_lspmic
  dhi(points),                                                                 &
                         ! CFL limit (s m-1)
  dhir(points),                                                                &
                         ! 1.0/DHI (m s-1)
  rho(points),                                                                 &
                         ! Air density (kg m-3)
  rhor(points),                                                                &
                         ! 1.0/Rho
  flux_fromabove(points)

    ! Intent (InOut)
real (kind=real_lsprec), intent(in out) ::                                     &
  fallspeed_thislayer(points)

    ! Intent (InOut)
real (kind=real_lsprec), intent(in out) ::                                     &
  mixratio_thislayer(points),                                                  &
  fallspeed_fromabove(points)

    ! Intent (InOut)
real (kind=real_lsprec), intent(in out) ::                                     &
  total_flux_out(points)

! Local variables

real (kind=real_lsprec) ::                                                     &
  mixratio_fromabove,                                                          &
                                 ! Mixing Ratio from above
  flux_out,                                                                    &
                                 ! Temporary flux out of layer
  expfactor                  ! Exponential Factor

integer :: i                    ! Loop counter

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='LSP_SEDIM_EULEXP'


!-----------------------------------------------------------------------

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
do i = 1, points

      ! -------------------------------------------------------------
      ! Adjust fall speeds to make a linear combination of the fall
      ! speed from the layer above. This ensures that the fall speed
      ! in this layer will not be calculated as zero even though
      ! there is mass falling into it from above.
      ! -------------------------------------------------------------

  mixratio_fromabove = flux_fromabove(i)*dhi(i)*rhor(i)

  if (mixratio_thislayer(i) + mixratio_fromabove  >   m0) then

    fallspeed_thislayer(i) =                                                   &
          (fallspeed_thislayer(i)*mixratio_thislayer(i)                        &
          +fallspeed_fromabove(i)*mixratio_fromabove)                          &
          /(mixratio_thislayer(i)+mixratio_fromabove)

  else

    fallspeed_thislayer(i) = zero

  end if

      ! -------------------------------------------------------------
      ! Eulerian solution with exponential limiter
      ! -------------------------------------------------------------

  if (fallspeed_thislayer(i)  >   zero) then

    expfactor = exp(-one*fallspeed_thislayer(i)*dhi(i))

        ! Calculate flux out of this layer

    flux_out = flux_fromabove(i)+dhir(i)                                       &
                 *(rho(i)*mixratio_thislayer(i)-flux_fromabove(i)              &
                 /fallspeed_thislayer(i))*(one-expfactor)

        ! Calculate mass (kg/kg) that remains in this layer

    mixratio_thislayer(i) = flux_fromabove(i)*rhor(i)                          &
                     / fallspeed_thislayer(i)                                  &
                     * (one-expfactor) + mixratio_thislayer(i) * expfactor

  else

        ! No fall out of the layer.
        ! Set MixingRatio to be the amount of mass falling in.
        ! FallSpeed can only be zero if MixingRatio_ThisLayer le M0
        ! and MixingRatio_FromAbove le M0
        ! so this is slightly inconsistent.
    flux_out       = zero
    mixratio_thislayer(i) = flux_fromabove(i)*rhor(i)*dhi(i)

  end if

      ! No need to compute fall speed out of the layer in this method
      ! -------------------------------------------------------------
      !  Total_Flux_Out is the flux out of this layer that falls
      !  into the next layer down
      ! -------------------------------------------------------------

  total_flux_out(i) = total_flux_out(i) + flux_out

      ! -------------------------------------------------------------
      ! Store fall speed in this layer to be the fallspeed from above
      ! for the next layer down
      ! -------------------------------------------------------------

  fallspeed_fromabove(i) = fallspeed_thislayer(i)

end do ! on loop over points

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine lsp_sedim_eulexp
end module lsp_sedim_eulexp_mod
