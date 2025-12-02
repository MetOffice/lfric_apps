! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!  subroutine TRSRCE ----------------------------------------------
!
!  Purpose: Adds source increments to a single level of the aerosol
!  field.
!
!  Suitable for single-column use.
!
!  Programming standard: Unified Model Documentation Paper No 3,
!
!
!  Arguments:---------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: tracer_advection

module trsrce_mod
use um_types, only: real_umphys

implicit none

character(len=*), parameter, private :: ModuleName='TRSRCE_MOD'

contains
subroutine trsrce(                                                             &
 rows, row_length, offx, offy, halo_i, halo_j, r_theta_levels, r_rho_levels,   &
 theta, q, qcl, qcf, exner, rho, tracer, srce,                                 &
 level, timestep, i_hour, i_minute, amp                                        &
)

use planet_constants_mod, only: kappa, c_virtual, pref, cp

use nlsizes_namelist_mod, only: model_levels

use conversions_mod, only: pi
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use errormessagelength_mod, only: errormessagelength

implicit none

integer, intent(in)    :: rows         ! number of P/U rows
integer, intent(in)    :: row_length   !
integer, intent(in)    :: offx         ! EW size of std. halo
integer, intent(in)    :: offy         ! NS size of std. halo
integer, intent(in)    :: halo_i       ! EW extended halo
integer, intent(in)    :: halo_j       ! NS extended halo
integer, intent(in)    :: i_hour       ! Local time hour
integer, intent(in)    :: i_minute     ! Local time minute
integer, intent(in)    :: level        ! level of the tracer

real(kind=real_umphys), intent(in)       :: timestep     ! Timestep in seconds
real(kind=real_umphys), intent(in)       :: amp          ! Amplitude of diurnal
                                       ! variation of emission

real, intent(in) :: r_theta_levels (1-halo_i:row_length+halo_i,                &
                                    1-halo_j:rows+halo_j,0:model_levels)
! height of theta levels from Earth centre
real, intent(in) :: r_rho_levels   (1-halo_i:row_length+halo_i,                &
                                    1-halo_j:rows+halo_j,model_levels)
! height of rho levels from Earth centre

real(kind=real_umphys), intent(in)       ::                                    &
      theta( 1 - offx : row_length + offx,     & ! pot. temperature
             1 - offy : rows + offy,           & !
             model_levels )                                                    &
,     q    ( 1 - halo_i :row_length + halo_i,  & ! Q on theta
             1 - halo_j :rows + halo_j,        & ! levels
             model_levels )                                                    &
,     qcl  ( 1 - halo_i :row_length + halo_i,  & ! Qcl on theta
             1 - halo_j :rows + halo_j,        & ! levels
             model_levels )                                                    &
,     qcf  ( 1 - halo_i :row_length + halo_i,  & ! Qcf on theta
             1 - halo_j :rows + halo_j,        & ! levels
             model_levels )                                                    &
,     exner( 1 - offx : row_length + offx,     & ! exner on rho
             1 - offy : rows + offy,           & ! levels
             model_levels + 1)                                                 &
,     rho  ( 1 - offx : row_length + offx,     & ! density * r * r
             1 - offy : rows + offy,           & ! on rho levels
             model_levels )

real(kind=real_umphys), intent(in)       ::                                    &
      srce( : , : )                              ! tracer source

real(kind=real_umphys), intent(in out)    ::                                   &
      tracer( 1 - offx : ,                     & ! level of tracer
              1 - offy : )                       ! to be updated


! Local, including save'd, storage------------------------------------

!  (a) Scalars effectively expanded to workspace by the Cray (using
!      vector registers).
real(kind=real_umphys)             :: dm
real(kind=real_umphys)             :: ts
real(kind=real_umphys)             :: thetav
                                  ! virtual potential temperature
real(kind=real_umphys)             :: exner_ave     ! an averaged exner term
real(kind=real_umphys)             :: rho_theta     ! rho on theta level
real(kind=real_umphys)             :: rho1          ! } values of rho after the
real(kind=real_umphys)             :: rho2          ! } r-squared factor removed

real(kind=real_umphys), parameter  :: Factor = 1.0
                                  ! Factor to multiply source term
real(kind=real_umphys), parameter  :: TZero  = 12.0 ! Time of maximum emissions

!  (b) Others.
integer          :: i   ! Loop counter
integer          :: j   ! Loop counter

character (len=*), parameter :: RoutineName='TRSRCE'

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

!-----------------------------------------------------------------------
!  Subroutine structure :
!  Loop over field adding source term*timestep/level mass/unit area.
!-----------------------------------------------------------------------

!     Allow for source varying with time of day
ts = 1.0
if (amp > 0.0) then
  ts = 1.0 +                                                                   &
  amp * cos( (real(i_hour) + real(i_minute)/60.0 - TZero)                      &
      * pi/12.0)
end if

ts = ts * timestep * factor

if (level < model_levels) then
!$OMP  PARALLEL do DEFAULT(SHARED) private(i,j, rho1, rho2, DM)                &
!$OMP  SCHEDULE(STATIC)
  do j=1, rows
    do i = 1, row_length
      ! Remove the r squared factor from rho before interpolation
      rho1 =  rho(i,j,level)/(r_rho_levels(i,j,level) *                        &
                              r_rho_levels(i,j,level) )
      rho2 =  rho(i,j,level+1)/(r_rho_levels(i,j,level+1) *                    &
                                r_rho_levels(i,j,level+1) )


      ! DM = density (interpolated on to theta levels) * delta r

      dm = rho2 * (r_theta_levels(i,j,level) -                                 &
                   r_rho_levels(i,j,level) ) +                                 &
           rho1 * (r_rho_levels(i,j,level+1) -                                 &
                   r_theta_levels(i,j,level) )
      !
      ! Special case for lowest layer to get correct mass
      if (level == 1) then
        dm = dm * (r_rho_levels(i,j,2) - r_theta_levels(i,j,0))                &
                / (r_rho_levels(i,j,2) - r_rho_levels(i,j,1))
      end if
      !
      ! Convert DM to DRY density
      dm = dm * (1.0 -q(i,j,level)-qcl(i,j,level)-qcf(i,j,level))
      !
      tracer(i, j) = tracer(i, j) + srce(i, j) * ts/dm
    end do
  end do
!$OMP end PARALLEL do

else  ! level = model_level
      !--------------------------------------------------------------------
      ! Cannot average here to get rho_theta. Hence calculate by using
      ! the equation of state vis
      !         ___r                     P0
      !         rho   = -----------------------------------------
      !                 kappa * Cp * ____________________r * theta_v
      !                              [          kappa-1 ]
      !                              [ exner ** ------- ]
      !                              [          kappa   ]
      !-------------------------------------------------------------------
!$OMP  PARALLEL do DEFAULT(SHARED) private(i,j, rho1, rho2, DM,                &
!$OMP  exner_ave, thetav, rho_theta) SCHEDULE(STATIC)
  do j=1, rows
    do i = 1, row_length

      thetav = theta(i,j,level) * (1.0 +                                       &
                 ( q(i,j,level) * C_Virtual )                                  &
               - qcl(i,j,level)-qcf(i,j,level)  )

      exner_ave = (exner(i,j,level) ** ((kappa - 1.0)/kappa) +                 &
                   exner(i,j,level+1) ** ((kappa - 1.0)/kappa))/2.0
      rho_theta = Pref/( kappa * Cp * exner_ave * thetav )

      ! rho_theta is at the top theta level. We also need the value
      ! at the top rho level. This will be rho1
      rho1 = rho(i,j,model_levels)/(r_rho_levels(i,j,model_levels)*            &
                                    r_rho_levels(i,j,model_levels))

      ! rho2 will be the average of rho1 and rho_theta
      rho2 = ( rho1 + rho_theta ) * 0.5

      dm = rho2 * (r_theta_levels(i,j,level) -                                 &
                   r_rho_levels(i,j,level) )

      tracer(i, j) = tracer(i, j) + srce(i, j) * ts/dm
    end do
  end do
!$OMP end PARALLEL do

end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine trsrce
end module trsrce_mod
