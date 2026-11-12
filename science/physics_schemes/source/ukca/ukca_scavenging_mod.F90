! *****************************COPYRIGHT*******************************
!
! (c) [University of Oxford] [2011]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Module for plume scavenging, contains definitions, and contained routines:
!    ukca_plume_scav   (called from convec2)
!    ukca_calc_aqueous (called from ukca_plume_scav)
!    ukca_set_conv_indices (called from addres)
!
!  Method:
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components provided by The University of Cambridge,
!  University of Leeds, University of Oxford, and The Met Office. See
!  www.ukca.ac.uk
!
!  Code Owner: Please refer to the UM file CodeOwners.txt
!  This file belongs in section: UKCA_UM
!
!  Code Description:
!    Language:  Fortran
!
! ######################################################################
module ukca_scavenging_mod

use ukca_mode_setup,           only: nmodes
use cv_run_mod,   only: i_convection_vn, i_convection_vn_6a,                   &
                        i_convection_vn_5a, i_cv_comorph
use parkind1, only: jprb, jpim
use yomhook, only: lhook, dr_hook
implicit none
private
save

real, parameter :: precip_lim = 1e-50   ! No calculation if ppn less than this
                                        ! also used to to test mass flux

! Information about how tracers are scavenged
type :: tr_scav_struct
  ! Index of first UKCA tracer in convection tracer array
  integer :: i_ukca_first
  ! Index of last UKCA  tracer in convection tracer array
  integer :: i_ukca_last
end type tr_scav_struct

! structure holding information for scavenging of all tracers
type (tr_scav_struct), public :: tracer_info

! Index of tracers for mode number, size is the number of modes.
! Holds the tracer number in the tracer_ukca array.
! Defaults to zero if the mode is not found.
integer, public, allocatable :: nmr_index_um(:)

! Index of tracers for mass mixing ratio of component, size is
! no. of modes * no. of components.
! Holds the appropriate tracer number in the tracer_ukca array.
! Defaults to -1 if the component is not found.
integer, public, allocatable :: mmr_index_um(:,:)

real    :: rscavn(nmodes)         ! Number scavenging ratios
real    :: rscavm(nmodes)         ! Mass scavenging ratios

public ukca_plume_scav, ukca_mode_scavcoeff, ukca_set_conv_indices

! Dr Hook parameters
integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1

character(len=*), parameter, private :: ModuleName='UKCA_SCAVENGING_MOD'

contains

! ######################################################################
subroutine ukca_set_conv_indices()
!
! set the entries in tracer_info to allow plume scavenging to
! correctly locate MODE tracers in the tot_tracer array in convection

! we need to set up mmr_index_um and nmr_index_um - the location
! in the tracer_ukca array of the number and mass mixing ratios
! Call from addres - outside of all OpenMP loops and only called once

use ereport_mod,               only: ereport
use errormessagelength_mod,    only: errormessagelength
use ukca_mode_tracer_maps_mod, only: arrayloc

use ukca_config_specification_mod,  only: glomap_variables

! Note that all values used from ukca_mode_setup are parameters so will be
! available before call to UKCA
use ukca_mode_setup,           only: nmodes, mode_names, cp_su, cp_cl,         &
                                     cp_bc, cp_oc, cp_du, cp_so, cp_nh4,       &
                                     cp_no3, cp_nn, cp_mp

use ukca_nmspec_mod,           only: nm_spec_active
use um_parcore,                only: mype
use umprintmgr,                only: umprint, ummessage, printstatus,          &
                                     prstatus_diag

implicit none

! Local variables

! Caution - pointers to type glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
integer, pointer :: ncp

integer :: icp                  ! component index
integer :: imode                ! mode index
integer :: icode                ! error code
integer :: i, j                 ! loop counter
#if defined(NAG_FORTRAN) && (NAG_FORTRAN == 7000000)
logical :: tracer_belongs
#endif
character(len=errormessagelength) :: cmessage   ! error message

real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='UKCA_SET_CONV_INDICES'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to type glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
ncp         => glomap_variables%ncp

if (.not. allocated(mmr_index_um)) then
  allocate(mmr_index_um(nmodes,ncp))
  mmr_index_um(:,:) = -1
end if

if (.not. allocated(nmr_index_um)) then
  allocate(nmr_index_um(nmodes))
  nmr_index_um(:) = 0
end if

! Construct the number and component mass index arrays
! by looping over all the active UKCA tracers
do i=1,size(nm_spec_active)
  ! (Slow) Workaround for NAG Fortran vn7.0 internal compiler error.
  ! Fixed in vn7.1
#if defined(NAG_FORTRAN) && (NAG_FORTRAN == 7000000)
  tracer_belongs = .false.
  do j=1, nmodes
    if ( mode_names(j) == nm_spec_active(i)(1:7) ) then
      tracer_belongs = .true.
    end if
  end do

  if ( tracer_belongs ) then
#else
  if ( any(mode_names(:) == nm_spec_active(i)(1:7)) ) then
#endif
      ! Tracer belongs to UKCA_MODE aerosol scheme

      ! first 7 characters of tracer name are the mode name
      ! imode is the position in the mode_names array
    imode = arrayloc(mode_names, nm_spec_active(i)(1:7))

    ! last 2 characters are the component name
    ! or 'N ' for number mixing ratio
    select case (nm_spec_active(i)(9:10))
    case ('N ')
      icp = 0
    case ('SU')
      icp = cp_su
    case ('SS')
      icp = cp_cl
    case ('BC')
      icp = cp_bc
    case ('OM')
      icp = cp_oc
    case ('DU')
      icp = cp_du
    case ('SO')
      icp = cp_so
    case ('NH')
      icp = cp_nh4
    case ('NT')
      icp = cp_no3
    case ('NN')
      icp = cp_nn
    case ('MP')
      icp = cp_mp
    case default
      icode = 1
      cmessage = ' Unknown component in case statement'
      write(ummessage,'(A40,A10)') cmessage,nm_spec_active(i)
      call umPrint(ummessage,src=RoutineName)
      call ereport(ModuleName//':'//RoutineName,icode,cmessage)
    end select

    ! now use this information to set up the indices for number and mass
    if (icp == 0) then
      ! This is a number density - set index in nmr_index
      nmr_index_um(imode) = i
    else if (icp > 0) then
      ! This is a component masss - set index in mmr_index_um
      mmr_index_um(imode,icp) = i
    end if
  end if
end do

! print out nmr and mmr indices
if (printstatus >= prstatus_diag .and. mype ==0) then
  write(ummessage,'(A)') 'number tracer maps'
  call umprint(ummessage,src=routinename)
  do i=1, nmodes
    write(ummessage,'(I6,1X,I6)') i, nmr_index_um(i)
    call umPrint(ummessage,src=routinename)
  end do
end if

if (printstatus >= prstatus_diag .and. mype ==0) then
  write(ummessage,'(A)') 'mass tracer maps'
  call umprint(ummessage,src=routinename)
  do j=1, ncp
    do i=1,nmodes
      write(ummessage,'(I6,1X,I6,1X,I6)') i, j, mmr_index_um(i, j)
      call umprint(ummessage,src=routinename)
    end do
  end do
end if

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine ukca_set_conv_indices

! ######################################################################
! Subroutine Interface:
subroutine ukca_plume_scav(ntra, npnts, trapkp1, xpkp1,                        &
                           prekp1, flxkp1)

! ----------------------------------------------------------------------
!  Calculate the change in parcel tracer content due to scavenging
!  by precipitation for UKCA tracers.
! ----------------------------------------------------------------------
use ukca_config_specification_mod, only: glomap_variables
use ukca_mode_setup,     only: mode_names
use ukca_nmspec_mod,     only: nm_spec_active
use planet_constants_mod, only: gg => g
use umprintmgr,          only: umprint, ummessage
use ereport_mod,         only: ereport

use errormessagelength_mod, only: errormessagelength

implicit none

integer, intent(in) :: ntra             ! Number of UKCA tracer variables
integer, intent(in) :: npnts            ! Number of points


real, intent(in)    ::                                                         &
  xpkp1(npnts)                ! parcel cloud condensate in layer k+1 (kg/kg)
                              ! liquid water for [6A] convection
                              ! liquid + frozen water for [4A/5A] convection

real, intent(in)    ::                                                         &
  prekp1(npnts)               ! precipitation from parcel as it rises from
                              ! layer k to k+1 (kg/m**2/s)
real, intent(in)    ::                                                         &
  flxkp1(npnts)               ! parcel massflux in layer k+1 (Pa/s)

real, intent(in out) ::                                                        &
  trapkp1(npnts,ntra)         ! parcel tracer content in layer k+1 (kg/kg)

! Local
integer :: ktra            ! counter
#if defined(NAG_FORTRAN) && (NAG_FORTRAN == 7000000)
integer :: j
logical :: is_a_mode_tracer
#endif
real :: traprekp1(npnts,ntra) ! Tracer content in precipitation produced
                              ! during ascent from layer K to K+1
                              ! [TR/kg(H2O), where tracer is TR/kg(air)]
real :: removal(npnts,ntra)   ! Removed tracer
real :: trapkp1_old(npnts,ntra) ! Temporary array to store original value of
                                ! trapkp1

! Caution - pointers to type glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
integer, pointer :: topmode

integer :: ifirst_ukca               ! index of first UKCA tracer in trapkp1
integer :: ilast_ukca                ! index of last UKCA tracer in trapkp1
integer :: iukca                     ! index of tracer in nm_spec_active
integer :: icode                     ! error code

character(len=errormessagelength) :: cmessage        ! error message

real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='UKCA_PLUME_SCAV'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to type glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
topmode     => glomap_variables%topmode

traprekp1(:,:) = 0.0
removal(:,:) = 0.0

! Calculate the aqueous phase mixing ratio for the UKCA tracers
! only if precipitiation is greater than a threshold
if (any(prekp1(:) > precip_lim)) then
  select case( i_convection_vn )
  case ( i_convection_vn_6a, i_cv_comorph )
    call ukca_calc_aqueous_6a(npnts, ntra, topmode,                            &
                              trapkp1, xpkp1(:),                               &
                              prekp1(:), flxkp1(:),                            &
                              traprekp1(:,:) )
  case ( i_convection_vn_5a )
    call ukca_calc_aqueous_5a(npnts, ntra, topmode,                            &
                              trapkp1,                                         &
                              xpkp1(:) + prekp1(:)*gg/flxkp1(:),               &
                              traprekp1(:,:) )
  case default
    cmessage = ' Plume scavenging is not supported for this convection version'
    icode = 1
    write(ummessage,'(A65,I4)') cmessage,i_convection_vn
    call umPrint(ummessage,src=RoutineName)
    call ereport(ModuleName//':'//RoutineName,icode,cmessage)
  end select       ! i_convection_vn
end if     ! prekp1(:) > precip_lim)


! Remove scavenged tracer along with water

! set tracer indices in to loop over
ifirst_ukca = tracer_info%i_ukca_first
ilast_ukca  = tracer_info%i_ukca_last

! Check whether there is precipitation - if not do nothing
if (any(prekp1(:) > precip_lim)) then

  ! loop over all UKCA tracers
  do ktra = ifirst_ukca,ilast_ukca

    ! for names in nm_spec_active need to index relative
    ! to the ukca_tracer array
    iukca = ktra - ifirst_ukca + 1

    ! check whether this is a MODE tracer by matching the
    ! first 7 characters of the tracer name against
    ! the list of mode names

    ! (Slow) Workaround for NAG Fortran vn7.0 internal compiler error.
    ! Fixed in vn7.1
#if defined(NAG_FORTRAN) && (NAG_FORTRAN == 7000000)
    is_a_mode_tracer = .false.
    do j=1, nmodes
      if (mode_names(j) == nm_spec_active(iukca)(1:7)) then
        is_a_mode_tracer = .true.
      end if
    end do
    if ( is_a_mode_tracer ) then
#else
    if ( any(mode_names(1:topmode) == nm_spec_active(iukca)(1:7)) ) then
#endif
        ! Only scavenge aerosols where there is a minimum
        ! amount of tracer in precip
      where (traprekp1(:,ktra) > 1.0e-30)

        ! calculate the loss of tracer
        removal(:,ktra) = traprekp1(:,ktra) * prekp1(:)*gg/flxkp1(:)
        ! store the previous value of the tracer
        trapkp1_old(:,ktra) = trapkp1(1:npnts,ktra)
        ! now reduce the tracer concentration by the removal amount
        trapkp1(1:npnts,ktra) = trapkp1(1:npnts,ktra) - removal(:,ktra)
      end where


      ! If negative tracer results, then limit the removal to give zero.
      where (trapkp1(1:npnts,ktra) < 0.0)
        trapkp1(1:npnts,ktra) = 0.0
        removal(:,ktra) = trapkp1_old(:,ktra)
      end where

    end if   ! Mode tracer
  end do     ! ktra
end if       ! any(prekp1 > 0)

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

return
end subroutine ukca_plume_scav

! ######################################################################
! Subroutine Interface:
subroutine ukca_calc_aqueous_5a(npnts, ntra, topmode, trapkp1, qcl,            &
                                aqueous)

!   Calculates the mixing ratio of UKCA tracers in the aqueous phase
!   for scavenging. For [4A/5A] convection scheme.

! UKCA_D1_DEFS etc won't have been initialised on the first call,
!  as convection happens before the first call to ukca_main1, so
!  we must only use constants here. See ukca_wetdp_addr

use ukca_config_specification_mod, only: glomap_variables

use umPrintMgr,      only: umPrint, ummessage
use ereport_mod,     only: ereport

use errormessagelength_mod, only: errormessagelength


implicit none

! Inputs
integer, intent(in) :: npnts                ! Number of grid boxes
integer, intent(in) :: ntra                 ! Number of UKCA MODE tracers
integer, intent(in) :: topmode              ! Last mode to scavenge over
real, intent(in)    :: trapkp1(npnts,ntra)  ! parcel tracer content in
                                            ! layer k+1 (kg/kg)
                                            ! n.b. all model tracers
real, intent(in)    :: qcl(npnts)           ! Liquid water mixing ratio
                                            ! [kg(liquid water)/kg(air)]

! Outputs
real, intent(out)   :: aqueous(npnts,ntra)  ! Tracer/water mixing ratio
                                            ! [TR(aq)/kg(water)]

! Local variables

! Caution - pointers to type glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
logical, pointer :: component(:,:)
logical, pointer :: mode (:)
integer, pointer :: ncp

integer :: j                    ! array index
integer :: icp                  ! component index
integer :: imode                ! mode index
! position of first UKCA tracer in all tracers array
integer :: ifirst_ukca
! position of last UKCA tracer in all tracers array
integer :: ilast_ukca
integer :: icode                ! error code
character (len=* ), parameter :: RoutineName='UKCA_CALC_AQUEOUS_5A'
character(len=errormessagelength) :: cmessage   ! error message

real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to type glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp

! Set this to zero for any other tracers (e.g. chemistry)
aqueous(:,:) = 0.0

ifirst_ukca = tracer_info%i_ukca_first
ilast_ukca = tracer_info%i_ukca_last

do imode=1,topmode
  if (mode(imode)) then
    if (nmr_index_um(imode) > 0) then
      ! For number tracers, [tracer] = [ptcl/molec(air)]
      !                              = [ptcl/kg(air)] * Boltzmann / R
      !                                i.e. [TR] = [ptcl] * Boltzmann / R

      ! Only calculate where non-negligible liquid water
      j = nmr_index_um(imode) + ifirst_ukca - 1
      if (j < ifirst_ukca .or. j > ilast_ukca) then
        cmessage = 'j is out of range for aerosol number'
        icode = abs(j)
        call umPrint(cmessage,src=RoutineName)
        write(ummessage,'(A,I5,1X,A,I5)') 'imode', imode, 'j= ', j
        call umPrint(ummessage,src=RoutineName)
        write(ummessage,'(A,I5,1X,I5)') 'should be between ',                  &
            ifirst_ukca, ilast_ukca
        call umPrint(ummessage,src=RoutineName)
        write(ummessage,'(A,I5)') 'icode: ',icode
        call umPrint(ummessage,src=RoutineName)
        call ereport(RoutineName,icode,cmessage)
      end if
      where (qcl(1:npnts) > 1e-10)
        aqueous(:,j) = rscavn(imode) *  trapkp1(1:npnts,j) / qcl(:)
      end where
      do icp=1,ncp
        if (component(imode,icp)) then
          if (mmr_index_um(imode,icp) > 0) then
            ! For mass tracers, [tracer] = [kg/kg(air)] i.e. [TR] = [kg]
            j = mmr_index_um(imode,icp) + ifirst_ukca - 1
            if (j < ifirst_ukca .or. j > ilast_ukca) then
              cmessage = 'j is out of range'
              icode = abs(j)
              call umPrint(cmessage,src=RoutineName)
              write(ummessage,'(A,I6,1X,I6,1x,A,I6)') 'imode, icp',            &
                  imode, icp, 'j= ', j
              call umPrint(ummessage,src=RoutineName)
              write(ummessage,'(A,I5,1X,I5)') 'should be between: ',           &
                  ifirst_ukca, ilast_ukca
              call umPrint(ummessage,src=RoutineName)
              write(ummessage,'(A,I5)') 'icode: ',icode
              call umPrint(ummessage,src=RoutineName)
              call ereport(RoutineName,icode,cmessage)
            end if
            where (qcl(:) > 1e-10)
              aqueous(:,j) = rscavm(imode) * trapkp1(1:npnts,j) / qcl(:)
            end where
          end if
        end if
      end do    ! icp
    end if     ! nmr_index_um > 0
  end if    ! mode(imode)
end do     ! imode

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ukca_calc_aqueous_5a

! ######################################################################
! Subroutine Interface:
subroutine ukca_calc_aqueous_6a(npnts, ntra, topmode, trapkp1,                 &
                                xpkp1, prekp1, flxkp1,                         &
                                aqueous)

!   Calculates the mixing ratio of UKCA tracers in the aqueous phase
!   for scavenging. For [6A] convection scheme.

! UKCA_D1_DEFS etc won't have been initialised on the first call,
!  as convection happens before the first call to ukca_main1, so
!  we must only use constants here. See ukca_wetdp_addr

use ukca_config_specification_mod, only: glomap_variables

use planet_constants_mod, only: gg => g
use umPrintMgr,      only: umPrint, ummessage
use ereport_mod,     only: ereport

use errormessagelength_mod, only: errormessagelength


implicit none

! Inputs
integer, intent(in) :: npnts                ! Number of grid boxes
integer, intent(in) :: ntra                 ! Number of UKCA MODE tracers
integer, intent(in) :: topmode              ! Last mode to scavenge over
real, intent(in)    :: trapkp1(npnts,ntra)  ! parcel tracer content in
                                            ! layer k+1 (kg/kg)
                                            ! n.b. all model tracers
real, intent(in)    ::                                                         &
  xpkp1(npnts)                ! parcel cloud water in layer k+1 (kg/kg)

real, intent(in)    ::                                                         &
  prekp1(npnts)               ! precipitation from parcel as it rises from
                              ! layer k to k+1 (kg/m**2/s)
real, intent(in)    ::                                                         &
  flxkp1(npnts)               ! parcel massflux in layer k+1 (Pa/s)


! Outputs
real, intent(out)   :: aqueous(npnts,ntra)  ! Tracer/water mixing ratio
                                            ! [TR(aq)/kg(water)]

! Local variables

! Caution - pointers to type glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
logical, pointer :: component(:,:)
logical, pointer :: mode (:)
integer, pointer :: ncp

real, parameter :: qcl_min = 1e-15      ! Minimum qcl for calculations
real, parameter :: trc_min = 1e-20      ! Minimum tracer for calculations

real          :: qcl(npnts)             ! Liquid water mixing ratio
                                        ! [TR/kg(air)]
integer :: j                    ! array index
integer :: icp                  ! component index
integer :: imode                ! mode index
integer :: ifirst_ukca          ! position of first UKCA tracer in
                                ! all tracers array
integer :: ilast_ukca           ! position of last UKCA tracer in
                                ! ukca tracer array
integer :: icode                ! error code
character (len=* ), parameter :: RoutineName='UKCA_CALC_AQUEOUS_6A'
character(len=errormessagelength) :: cmessage   ! error message
logical :: todo(npnts)          ! mask for calculation

real(kind=jprb)               :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! Caution - pointers to type glomap_variables%
!           have been included here to make the code easier to read
!           take care when making changes involving pointers
component   => glomap_variables%component
mode        => glomap_variables%mode
ncp         => glomap_variables%ncp

! Set this to zero for any other tracers (e.g. chemistry)
aqueous(:,:) = 0.0

!  Consistency check to avoid floating point error
where (prekp1 > precip_lim .and. flxkp1 > precip_lim)
  qcl(:) = xpkp1(:) + prekp1(:)*gg/flxkp1(:)
  todo(:) = .true.
else where
  todo(:) = .false.
  qcl(:) = 0.0
end where

ifirst_ukca = tracer_info%i_ukca_first
ilast_ukca = tracer_info%i_ukca_last

do imode=1,topmode
  if (mode(imode)) then
    if (nmr_index_um(imode) > 0) then
      ! For number tracers, [tracer] = [ptcl/molec(air)]
      !                              = [ptcl/kg(air)] * Boltzmann / R
      !                                i.e. [TR] = [ptcl] * Boltzmann / R

      j = nmr_index_um(imode) + ifirst_ukca - 1
      if (j < ifirst_ukca .or. j > ilast_ukca) then
        cmessage = 'j is out of range for aerosol number'
        icode = abs(j)
        call umPrint(cmessage,src=RoutineName)
        write(ummessage,'(A,I5,1X,A,I5)') 'imode', imode, 'j= ', j
        call umPrint(ummessage,src=RoutineName)
        write(ummessage,'(A,I5,1X,I5)') 'should be between ',                  &
            ifirst_ukca, ilast_ukca
        call umPrint(ummessage,src=RoutineName)
        write(ummessage,'(A,I5)') 'icode: ',icode
        call umPrint(ummessage,src=RoutineName)
        call ereport(RoutineName,icode,cmessage)
      end if
      where (todo .and. trapkp1(1:npnts,j) > trc_min .and. qcl(:) > qcl_min)
        aqueous(:,j) = rscavn(imode) * trapkp1(1:npnts,j) / qcl(:)
      end where
      do icp=1,ncp
        if (component(imode,icp)) then
          if (mmr_index_um(imode,icp) > 0) then
            ! For mass tracers, [tracer] = [kg/kg(air)] i.e. [TR] = [kg]
            j = mmr_index_um(imode,icp) + ifirst_ukca - 1
            if (j < ifirst_ukca .or. j > ilast_ukca) then
              cmessage = 'j is out of range'
              icode = abs(j)
              call umPrint(cmessage,src=RoutineName)
              write(ummessage,'(A,I6,1X,I6,1x,A,I6)') 'imode, icp',            &
                  imode, icp, 'j= ', j
              call umPrint(ummessage,src=RoutineName)
              write(ummessage,'(A,I5,1X,I5)') 'should be between: ',           &
                  ifirst_ukca, ilast_ukca
              call umPrint(ummessage,src=RoutineName)
              write(ummessage,'(A,I5)') 'icode: ',icode
              call umPrint(ummessage,src=RoutineName)
              call ereport(RoutineName,icode,cmessage)
            end if
            where (todo .and. trapkp1(1:npnts,j) > trc_min .and.               &
              qcl(:) > qcl_min)
              aqueous(:,j) = rscavm(imode) * trapkp1(1:npnts,j) / qcl(:)
            end where
          end if
        end if
      end do    ! icp
    end if     ! nmr_index > 0
  end if    ! mode(imode)
end do     ! imode

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ukca_calc_aqueous_6a

! ######################################################################
subroutine ukca_mode_scavcoeff()
! Set the scavenging co-efficients to be used in plume scavenging
! set depending on the switch i_mode_nucscav set in ukca_option_mod
! call after namelists read
! ######################################################################

use ereport_mod,               only: ereport
use umprintmgr,                only: umprint, ummessage
use errormessagelength_mod, only: errormessagelength
use ukca_option_mod, only: mode_aitsol_cvscav, i_mode_nucscav

implicit none

logical, parameter :: l_convective = .true. ! True if convective cloud

integer :: icode
real(kind=jprb)               :: zhook_handle

character(len=errormessagelength) :: cmessage   ! error message
character(len=*), parameter :: RoutineName='UKCA_MODE_SCAVCOEFF'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

select case (i_mode_nucscav)
case (1,3)
  if (l_convective) then
    ! Original GLOMAP-MODE scavenging ratios, plus scavenging
    ! of part of the soluble Aitken mode due to the stronger
    ! updraughts and higher supersaturations in convective cloud.
    !              1=NS  2=KS   3=AS  4=CS  5=KI  6=AI  7=CI  8=SI
    rscavn(:) = [ 0.00, mode_aitsol_cvscav, 1.00, 1.00, 0.00, 0.00,            &
                   0.00, 0.00 ]
    rscavm(:) = [ 0.00, mode_aitsol_cvscav, 1.00, 1.00, 0.00, 0.00,            &
                   0.00, 0.00 ]
  else
    ! Original GLOMAP-MODE scavenging ratios
    !              1=NS  2=KS  3=AS  4=CS  5=KI  6=AI  7=CI  8=SI
    rscavn(:) = [ 0.00, 0.00, 1.00, 1.00, 0.00, 0.00, 0.00, 0.00 ]
    rscavm(:) = rscavn(:)
  end if
case (2)
  ! ECHAM5-HAM scavenging ratios (Stier et al. 2005)
  if (l_convective) then
    ! for convective mixed cloud
    !              1=NS  2=KS  3=AS  4=CS  5=KI  6=AI  7=CI  8=SI
    rscavn(:) = [ 0.20, 0.60, 0.99, 0.99, 0.20, 0.40, 0.40, 0.40 ]
  else
    ! for stratiform liquid cloud
    !              1=NS  2=KS  3=AS  4=CS  5=KI  6=AI  7=CI  8=SI
    rscavn(:) = [ 0.10, 0.25, 0.85, 0.99, 0.20, 0.40, 0.40, 0.40 ]
  end if
  rscavm(:) = rscavn(:)
case default
  cmessage = 'i_mode_nucscav out of range'
  icode = 1
  write(ummessage,'(A40,I8)') cmessage,i_mode_nucscav
  call umPrint(ummessage,src=RoutineName)
  call ereport(RoutineName,icode,cmessage)
end select ! i_mode_nucscav

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine ukca_mode_scavcoeff

end module ukca_scavenging_mod
