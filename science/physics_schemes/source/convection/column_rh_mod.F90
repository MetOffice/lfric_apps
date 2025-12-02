! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! To calculate the column integrated relative humidity
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection
module column_rh_mod


use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

! Model level heights from centre of Earth
use level_heights_mod, only:                                                   &
  r_theta_levels             &  ! Radii on theta levels (m)
 ,r_rho_levels                  ! Radii on rho levels (m)

! Number of model levels
use nlsizes_namelist_mod, only: model_levels

use qsat_mod, only: qsat, qsat_mix

use gen_phys_inputs_mod, only: l_mr_physics
use bl_option_mod, only: zero
use um_types, only: r_bl

implicit none

character(len=*), parameter, private :: ModuleName='COLUMN_RH_MOD'

contains

subroutine calc_column_rh(                                                     &
   npnts, nlcl                                                                 &
   , row_length,rows                                                           &
   , index_i, index_j                                                          &
   , p, t, q, rho_only                                                         &
   , z_half, zmax                                                              &
   , qsat_calc, column_rh, column_rh_bl                                        &
   , column_q )

implicit none

integer, intent(in) ::                                                         &
  npnts                    & ! Number of points
, nlcl(npnts)                ! Level of LCL

integer, intent(in) ::                                                         &
  row_length               & ! Local number of points on a row
, rows                       ! Local number of rows in a theta field

integer, intent(in) ::                                                         &
  index_i(npnts)           & ! Column number of unstable points
, index_j(npnts)             ! Row number of unstable points

real(kind=r_bl), intent(in) ::                                                 &
   p(npnts, model_levels)    & ! pressure (Pa)
   , t(npnts, model_levels)  & ! temperature (K)
   , q(npnts, model_levels)    ! water vapour (kg/kg)

real(kind=r_bl), intent(in) ::                                                 &
   rho_only(row_length,rows,model_levels)
                             ! Density (kg/m3)

real(kind=r_bl), intent(in) ::                                                 &
   z_half(npnts, model_levels) ! Height on half levels

real(kind=r_bl), intent(in) ::                                                 &
   zmax                      ! Top of column

real(kind=r_bl), intent(out) ::                                                &
   qsat_calc(npnts, model_levels) ! q saturation value

real(kind=r_bl), intent(out) ::                                                &
   column_rh(npnts)        & ! Column integrated RH (fraction)
   , column_rh_bl(npnts)     ! Column integrated RH up to LCL

real(kind=r_bl), intent(out) ::                                                &
   column_q(npnts)           ! Column integrated q

!------------------------------
! Local variables
!------------------------------

real(kind=r_bl) :: column_qsat          ! column_qsat
real(kind=r_bl) :: temp_mass            ! column_qsat

integer :: k,ii,i,j          ! loop counters

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='CALC_COLUMN_RH'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!$OMP PARALLEL do SCHEDULE(STATIC) DEFAULT(none) private(k)                    &
!$OMP SHARED(model_levels,l_mr_physics,qsat_calc,t,p,npnts)
do  k = 1,model_levels
  if ( l_mr_physics ) then
    call qsat_mix(qsat_calc(:,k),t(:,k),p(:,k),npnts)
  else
    call qsat(qsat_calc(:,k),t(:,k),p(:,k),npnts)
  end if
end do
!$OMP end PARALLEL do

do ii=1,npnts
  i = index_i(ii)
  j = index_j(ii)

  column_q(ii) = zero
  column_qsat  = zero

  do k = 2,model_levels-1

    if (z_half(ii,k) <= zmax) then

      temp_mass = rho_only(i,j,k)                                              &
                *  r_theta_levels(i,j,k) * r_theta_levels(i,j,k)               &
                * (r_rho_levels(i,j,k+1) - r_theta_levels(i,j,k))

      column_q(ii) = column_q(ii) + temp_mass * q(ii,k)
      column_qsat  = column_qsat  + temp_mass * qsat_calc(ii,k)

      if ( k == nlcl(ii) ) then
        column_rh_bl(ii) = column_q(ii)/column_qsat
      end if
    end if

  end do

  column_rh(ii)=column_q(ii)/column_qsat

end do
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return


end subroutine calc_column_rh

end module column_rh_mod
