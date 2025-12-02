! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
!   A module containing constants/parameters used in the dust scheme
!
module dust_parameters_mod

!
! Description:
!   This module contains declarations for constants and tunable
!   parameters used to diagnose the emission and deposition of
!   mineral dust
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Aerosols
!
! Code description:
!   Language: Fortran 90
!   This code is written to UM programming standards version 8.2.
!

use missing_data_mod, only: rmdi, imdi
use ereport_mod, only: ereport
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim
use errormessagelength_mod, only: errormessagelength
use constants_mod, only: r_def
use um_types, only: real_umphys

implicit none

!
! Parameters fundamental to the dust scheme
! (Details of dust properties/size divisions etc.)
! =======================================================================

! Prognostic mineral dust aerosol
logical :: l_dust
! Diagnostic mineral dust aerosol lifting
logical :: l_dust_diag
! Diagnostic mineral dust aerosol lifting without resistance calculations
! for dry deposition: calculate emission fluxes only
logical :: l_dust_flux_only

! Logical to contain whether using two-bin (true) or six-bin (false) dust
logical :: l_twobin_dust = .false.
!   (default to six-bin dust)

! Define which dust bins are active
logical :: l_dust_div1, l_dust_div2, l_dust_div3, l_dust_div4,                 &
           l_dust_div5, l_dust_div6

! Define which dust bins are active for LBCs
logical :: l_dust_div1_lbc, l_dust_div2_lbc, l_dust_div3_lbc,                  &
           l_dust_div4_lbc, l_dust_div5_lbc, l_dust_div6_lbc

! Define if dust emission potential is being read in to make a scaling
logical :: l_dust_emp_sc = .false.

! Number of discrete particle size divisions (bins)
integer :: ndiv              ! number of divisions that can be lifted
                             !  from the surface
integer, parameter :: ndivh=9! number of divisions that can be blown
                             !  horizontally along the surface and contribute
                             !  to the lifting of the 1st NDIV divisions
integer, parameter :: ndivl=6! number of divisions in dust soil ancillaries

integer :: i_dust = imdi                    ! dust scheme setting in namelist
integer, parameter :: i_dust_off = 0        ! No dust
integer, parameter :: i_dust_prognostic = 1 ! Prognostic dust
integer, parameter :: i_dust_diagnostic = 2 ! Diagnostic dust
integer, parameter :: i_dust_flux = 3       ! Dust emission fluxes only
                                            !  (no resitances or UM diagnostics)

! The size of particles included in each division
! Note that by using two arrays for the max/min diameter here we can set
! up overlapping divisions. We must take care to make them consistent

real(kind=real_umphys), allocatable ::  dmax(:)
        ! max diameter of particles in each div.
real(kind=real_umphys), allocatable ::  dmin(:)
        ! min diameter of particles in each div.
real(kind=real_umphys), allocatable ::  drep(:)
        ! representative particle diameter
!
! Physical properties of the dust
real(kind=real_umphys), parameter :: rhop = 2.65e+3
                                   ! density of a dust particle (quartz)
real(r_def), parameter :: rhop_def = real(rhop, r_def)
!
! Parameters used during dust emissions calculations
! =======================================================================
!
! Parameters based on observations/published research
real(kind=real_umphys), parameter :: ustd_bas(ndivh) =                         &
       [ 0.85, 0.72, 0.59, 0.46, 0.33, 0.16, 0.14, 0.18, 0.28 ]
                                      ! impact U*t derived from Bagnold (1941)
real(kind=real_umphys), parameter :: horiz_c = 2.61
                                      ! C in horizontal flux calc (White 1979)
real(kind=real_umphys), parameter :: vert_a = 13.4
                                      ! A in vertical flux calc(Gillette 1979)
real(kind=real_umphys), parameter :: vert_b = -6.0
                                       ! B in vertical flux calc(Gillette 1979)
real(kind=real_umphys), parameter :: vert_c = 0.01      ! cgs to si conversion

!
! Input variables needed when using the dust lifting scheme with 1 tile
!
real(kind=real_umphys)            :: z0_soil
                                      ! Roughness length over bare soil

!
! Tuning parameters - set by namelist
real(kind=real_umphys)            :: us_am = rmdi
                                      ! ustar correction (multiplic.)
real(kind=real_umphys)            :: sm_corr = rmdi
                                      ! soil moist. correction factor
real(kind=real_umphys)            :: horiz_d = rmdi
                                      ! Global tuning param for horizontal
                                      !  (and hence vertical) flux
! Tuning parameters - defined here
real(kind=real_umphys), parameter :: us_aa  = 0.0
                                      ! ustar correction (additional)
real(kind=real_umphys), parameter :: ust_aa = 0.0
                                      ! ustar_t correction (add)
real(kind=real_umphys), parameter :: ust_am = 1.0
                                      ! ustar_t correction (multi.)
!
! Limits used in dust flux calculations
real(kind=real_umphys)            :: u_s_min = 1.0e-5
                                       ! Minimum val of u_s_std,
                                      !  below which no dust is produced
                                      !  (avoids divide by zero problem)
real(kind=real_umphys), parameter :: clay_max = 0.1     ! Max clay fraction.
real(kind=real_umphys)            :: snowmin = 1.0e-6
                                       ! Min snow depth for dust
real(kind=real_umphys)            :: h_orog_limit = 150.0
                                       ! 1/2pk-trough height above which
                                      !  no dust is produced
real(kind=real_umphys), parameter :: fland_lim =0.99999
                                      ! No ems if fland<lim as windspeed
                                      !  too high at coastal points
!
! Switch to diagnose vertical flux using a fixed (user definable) size
!  distribution. If set to false, the vertical flux in each bin at a point
!  is poportional to the horizontal flux in that bin at that point.
logical         :: l_fix_size_dist = .false.

! does the clay_max term get used as the max, or a constant scaling
logical :: l_dust_clay_as_max = .false.

! Namelist components can't be allocatable, so set to six, but in two-bin
! case should only use the first two array positions.
real(kind=real_umphys) :: size_dist(6) = rmdi

!
! Switch to allow emission of dust on tiles other than the bare soil tile
! 0 allows emission only on bare soil, 1 uses a bare soil radiative frac
! from the LAI, as for surface albedo (scaling by dust_veg_sc_nojules),
! and other methods (2+ could be added later).
integer         :: dust_veg_emiss = imdi

! Parameters used during boundary layer mixing of dust
! =======================================================================
! Parameter for BL mixing of dust
real :: dust_bl_mixfac = 1.0

! Parameters used during the gravitational settling of dust
! =======================================================================
!

real(kind=real_umphys), parameter :: accf = 1.257
                                ! Cunningham correction factor term A
real(kind=real_umphys), parameter :: bccf = 0.4
                                ! Cunningham correction factor term B
real(kind=real_umphys), parameter :: cccf = -1.1
                                ! Cunningham correction factor term C

!
! The total dust deposition flux (required by ocean biogeochemistry)
! =======================================================================
!
real(kind=real_umphys), allocatable :: tot_dust_dep_flux(:,:)

!
! Parameters used during the scavenging of dust
! =======================================================================
!
real(kind=real_umphys), allocatable :: krain_dust(:)
real(kind=real_umphys), allocatable :: ksnow_dust(:)
!
! Parameters for dust diagnostics
real(kind=real_umphys) :: pwsdiag_sfc_em = rmdi
                              ! fraction of dust after emission in PWS diags
!
! RUN_Dust namelist via which non-parameter values herein can be set
! =======================================================================
!
namelist /RUN_Dust/                                                            &
     us_am, sm_corr, horiz_d, l_fix_size_dist, dust_veg_emiss,                 &
     i_dust, l_twobin_dust, l_dust_div1_lbc, l_dust_div2_lbc,                  &
     l_dust_div3_lbc, l_dust_div4_lbc, l_dust_div5_lbc, l_dust_div6_lbc,       &
     pwsdiag_sfc_em, l_dust_emp_sc, l_dust_clay_as_max, dust_bl_mixfac

!DrHook-related parameters
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1

character(len=*), parameter, private :: ModuleName='DUST_PARAMETERS_MOD'

contains
!
! Internal subroutine to check that entries to RUN_dust are consistent
! =======================================================================
!
subroutine dust_parameters_check( )
!
!   Check that a user-defined emissions size distribution is normalised
!
implicit none

if (l_fix_size_dist) then
  size_dist(:)=size_dist(:)/sum(size_dist)
end if

end subroutine dust_parameters_check


! Internal subroutine to set the default value of size_dist. This may be
! overwritten by setting size_dist in the namelist (hence it's in a separate
! subroutine) or by setting l_twobin_dust. If l_fix_size_dist is set it is
! possible to override size_dist in the namelist as that was the previous
! behaviour, and it's mandatory to use it in the case of l_twobin_dust. If
! l_fix_size_dist is not set the value here is irrelevant as it's not used

subroutine dust_size_dist_initialise( )

implicit none

size_dist(1:6) =                                                               &
      [ 0.0005, 0.0049, 0.0299, 0.2329, 0.4839, 0.2479 ]

end subroutine dust_size_dist_initialise

! Internal subroutine to load the correct parameters into the variables in this
! module depending on whether two- or six- bin dust is being used
subroutine dust_parameters_load( )

implicit none

integer           :: icode           ! Error code
character(len=errormessagelength) :: cmessage        ! Error message
real(kind=jprb) :: zhook_handle

character(len=*), parameter :: RoutineName='DUST_PARAMETERS_LOAD'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Convert i_dust into appropriate logical settings
select case(i_dust)
case (i_dust_off)
  l_dust      = .false.
  l_dust_diag = .false.
  l_dust_flux_only = .false.
case (i_dust_prognostic)
  l_dust      = .true.
  l_dust_diag = .false.
  l_dust_flux_only = .false.
case (i_dust_diagnostic)
  l_dust      = .false.
  l_dust_diag = .true.
  l_dust_flux_only = .false.
case (i_dust_flux)
  l_dust      = .false.
  l_dust_diag = .false.
  l_dust_flux_only = .true.
case DEFAULT
  icode = 2
  write(cmessage, '(A,I2,A)') 'i_dust = ', i_dust, ' is not a valid'           &
                            // ' setting'
  call ereport('DUST_PARAMETERS_LOAD', icode, cmessage)
end select

! Ascertain which dust bins are active
if (i_dust == i_dust_prognostic) then
  l_dust_div1 = .true.
  l_dust_div2 = .true.
  if (l_twobin_dust) then
    l_dust_div3 = .false.
    l_dust_div4 = .false.
    l_dust_div5 = .false.
    l_dust_div6 = .false.
    l_dust_div3_lbc = .false.
    l_dust_div4_lbc = .false.
    l_dust_div5_lbc = .false.
    l_dust_div6_lbc = .false.
  else
    l_dust_div3 = .true.
    l_dust_div4 = .true.
    l_dust_div5 = .true.
    l_dust_div6 = .true.
  end if
else
  l_dust_div1 = .false.
  l_dust_div2 = .false.
  l_dust_div3 = .false.
  l_dust_div4 = .false.
  l_dust_div5 = .false.
  l_dust_div6 = .false.
  l_dust_div1_lbc = .false.
  l_dust_div2_lbc = .false.
  l_dust_div3_lbc = .false.
  l_dust_div4_lbc = .false.
  l_dust_div5_lbc = .false.
  l_dust_div6_lbc = .false.
end if

! Mandatory use of fixed size distribution when using two-bin code
if (l_twobin_dust) l_fix_size_dist = .true.

! Define ndiv
if (l_twobin_dust) then
  ndiv = 2
else
  ndiv = 6
end if

! Allocate scavenging co-efficients
allocate(krain_dust(ndiv))
allocate(ksnow_dust(ndiv))
allocate(dmax(ndiv))
allocate(dmin(ndiv))
allocate(drep(ndiv))


! Set variables depending on how many bins
if (l_twobin_dust) then

  ! scav. coeff. for rain
  krain_dust(1:ndiv) = [ 4.0e-5, 3.0e-4 ]
  ! scav. coeff. for snow
  ksnow_dust(1:ndiv) = [ 4.0e-5, 3.0e-4 ]

  dmax(1:ndiv) =  [ 4.0e-6, 2.0e-5 ]
            ! max diameter of particles in each div.
  dmin(1:ndiv) =  [ 2.0e-7, 4.0e-6 ]
            ! min diameter of particles in each div.
  drep(1:ndiv) =  [ 2.0e-6, 8.0e-6 ]
            ! representative particle diameter

  ! Two-bin dust must have this size distribution
  size_dist(1:ndivl) =                                                         &
      [ 0.1800, 0.8200, 0.0000, 0.0000, 0.0000, 0.0000 ]

else ! Six-bin dust

  ! scav. coeff. for rain
  krain_dust(1:ndiv) =                                                         &
   [ 2.0e-5, 2.0e-5, 3.0e-5, 6.0e-5, 4.0e-4, 4.0e-4 ]
  ! scav. coeff. for snow
  ksnow_dust(1:ndiv) =                                                         &
   [ 2.0e-5, 2.0e-5, 3.0e-5, 6.0e-5, 4.0e-4, 4.0e-4 ]


  dmax(1:ndiv) =                                                               &
   [ 2.0e-7, 6.32456e-7, 2.0e-6, 6.32456e-6, 2.0e-5, 6.32456e-5 ]
            ! max diameter of particles in each div.
  dmin(1:ndiv) =                                                               &
   [ 6.32456e-8, 2.0e-7, 6.32456e-7, 2.0e-6, 6.32456e-6, 2.0e-5 ]
            ! min diameter of particles in each div.
  drep(1:ndiv) =                                                               &
                 [ 0.112468e-06, 0.355656e-06, 0.112468e-05,                   &
                    0.355656e-05, 0.112468e-04, 0.355656e-04 ]
            ! representative particle diameter


end if ! l_twobin_dust

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine dust_parameters_load


! internal subroutine to deallocate the dust arrays - provided for consistency,
! not actually called at present
subroutine dust_parameters_unload( )

implicit none
deallocate (drep)
deallocate (dmin)
deallocate (dmax)
deallocate (ksnow_dust)
deallocate (krain_dust)

end subroutine dust_parameters_unload

subroutine print_nlist_run_dust()
use umPrintMgr, only: umPrint
implicit none
character(len=50000) :: lineBuffer
real(kind=jprb) :: zhook_handle

character(len=*), parameter :: RoutineName='PRINT_NLIST_RUN_DUST'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call umPrint('Contents of namelist run_dust',                                  &
    src='dust_parameters_mod')

write(lineBuffer,'(A,F0.5)')' us_am = ',us_am
call umPrint(lineBuffer,src='dust_parameters_mod')
write(lineBuffer,'(A,F0.5)')' sm_corr = ',sm_corr
call umPrint(lineBuffer,src='dust_parameters_mod')
write(lineBuffer,'(A,F0.5)')' horiz_d = ',horiz_d
call umPrint(lineBuffer,src='dust_parameters_mod')
write(lineBuffer,'(A,L1)')' l_fix_size_dist = ',l_fix_size_dist
call umPrint(lineBuffer,src='dust_parameters_mod')
write(lineBuffer,'(A,I0)')' dust_veg_emiss = ',dust_veg_emiss
call umPrint(lineBuffer,src='dust_parameters_mod')
write(lineBuffer,'(A,L1)')' l_twobin_dust = ',l_twobin_dust
call umPrint(lineBuffer,src='dust_parameters_mod')
write(lineBuffer,'(A,L1)') ' l_dust_div1_lbc = ',l_dust_div1_lbc
call umPrint(lineBuffer,src='dust_parameters_mod')
write(lineBuffer,'(A,L1)') ' l_dust_div2_lbc = ',l_dust_div2_lbc
call umPrint(lineBuffer,src='dust_parameters_mod')
write(lineBuffer,'(A,L1)') ' l_dust_div3_lbc = ',l_dust_div3_lbc
call umPrint(lineBuffer,src='dust_parameters_mod')
write(lineBuffer,'(A,L1)') ' l_dust_div4_lbc = ',l_dust_div4_lbc
call umPrint(lineBuffer,src='dust_parameters_mod')
write(lineBuffer,'(A,L1)') ' l_dust_div5_lbc = ',l_dust_div5_lbc
call umPrint(lineBuffer,src='dust_parameters_mod')
write(lineBuffer,'(A,L1)') ' l_dust_div6_lbc = ',l_dust_div6_lbc
call umPrint(lineBuffer,src='dust_parameters_mod')
write(lineBuffer,'(A,L1)') ' l_dust_clay_as_max = ',l_dust_clay_as_max
call umPrint(lineBuffer,src='dust_parameters_mod')
write(lineBuffer,'(A,L1)') ' l_dust_emp_sc = ',l_dust_emp_sc
call umPrint(lineBuffer,src='dust_parameters_mod')
write(lineBuffer,'(A,F16.4)') ' pwsdiag_sfc_em = ',pwsdiag_sfc_em
call umPrint(lineBuffer,src='dust_parameters_mod')
write(lineBuffer,'(A,F16.4)') ' dust_bl_mixfac = ',dust_bl_mixfac
call umPrint(lineBuffer,src='dust_parameters_mod')

call umPrint('- - - - - - end of namelist - - - - - -',                        &
    src='dust_parameters_mod')

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine print_nlist_run_dust


end module dust_parameters_mod
