! *****************************COPYRIGHT*******************************
!
! (c) [University of Cambridge] [2008]. All rights reserved.
! This routine has been licensed to the Met Office for use and
! distribution under the UKCA collaboration agreement, subject
! to the terms and conditions set out therein.
! [Met Office Ref SC138]
!
! *****************************COPYRIGHT*******************************
!
!  Description:
!   Module containing subroutine UKCA_VOLCANIC_SO2 that performs
!   interactive emissions of SO2 from explosive volcanic eruptions
!   into the stratosphere.
!
!   Requires an external csv file with eruption properties and
!   SO2 emissions from explosive eruptions.
!   Calls the Plumeria model to calculate plume height of the eruption
!   if the logical l_ukca_so2ems_plume is set as true.
!
!  UKCA is a community model supported by The Met Office and
!  NCAS, with components initially provided by The University of
!  Cambridge, University of Leeds and The Met Office. See
!  www.ukca.ac.uk
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: ukca_um
!
!  Code Description:
!    Language:  FORTRAN 90
!
! ######################################################################
!
module ukca_volcanic_so2_mod

implicit none

character(len=*), parameter, private :: ModuleName = 'UKCA_VOLCANIC_SO2_MOD'

contains

subroutine ukca_volcanic_so2                                                   &
          (so2_mmr, mass, row_length, rows, model_levels,                      &
           year, timestep, r_theta_levels, rel_humid_frac,                     &
           p_theta_levels, t_theta_levels,                                     &
           geopH_on_theta_mlevs, u_rho_levels, v_rho_levels,                   &
           plumeria_height)

use model_domain_mod,     only: model_type, mt_global, mt_cyclic_lam,          &
                                l_cartesian, l_regular
use cderived_mod,         only: delta_lambda, delta_phi
use planet_constants_mod, only: planet_radius
use conversions_mod,      only: pi_over_180, rsec_per_day, pi
use model_time_mod,       only: i_day_number, i_hour, i_minute, i_second
use ereport_mod,          only: ereport
use errormessagelength_mod, only: errormessagelength
use yomhook, only: lhook, dr_hook
use parkind1, only: jprb, jpim

use file_manager,        only: assign_file_unit, release_file_unit
use filenamelength_mod,  only: filenamelength
use ukca_option_mod,     only: nerupt, file_volc_so2, l_ukca_so2ems_plumeria
use plumeria_main_mod,   only: plumeria_main

implicit none

! input parameters:

integer, intent(in) :: row_length, rows, model_levels

integer, intent(in) :: year
real, intent(in) :: timestep

! input fields
real, intent(in) :: mass(row_length, rows, model_levels)
real, intent(in) :: r_theta_levels(row_length,rows,model_levels)
real, intent(in) :: rel_humid_frac(row_length,rows,model_levels)
real, intent(in) :: p_theta_levels(row_length,rows,model_levels)
real, intent(in) :: t_theta_levels(row_length,rows,model_levels)
real, intent(in) :: u_rho_levels(row_length,rows,model_levels)
real, intent(in) :: v_rho_levels(row_length,rows,model_levels)
real, intent(in) :: geopH_on_theta_mlevs(row_length,rows,model_levels)

! in/out field
! Mass mixing ratio of SO2 in kg(S)/kg(air)
real, intent(in out) :: so2_mmr(row_length, rows, model_levels)
real, intent(out) :: plumeria_height(row_length, rows)

! Parameters from explosive volcanic SO2 emission file
integer, parameter :: nparams = 9 ! number of eruption parameters
real, allocatable, save :: volc_emissions(:,:)
real, allocatable, save :: lat_volc(:)
real, allocatable, save :: lon_volc(:)
integer, allocatable, save :: year_volc(:)
integer, allocatable, save :: start_volc(:)
integer, allocatable, save :: end_volc(:)
real, allocatable, save :: magnitude(:)
real, allocatable, save :: height_center(:)
real, allocatable, save :: vent_height(:)
real, allocatable, save :: mer(:)

integer, allocatable, save :: lat_index(:)
integer, allocatable, save :: lon_index(:)

real, allocatable :: lat_volc_rad(:)
real, allocatable :: lon_volc_rad(:)

! Parameters for gaussian plume injection and plume height calculation
real, allocatable, save :: plume_thickness(:)
real, allocatable, save :: height_top_gaussian(:)
real, allocatable, save :: height_bot_gaussian(:)
real, allocatable :: weight_volc_vertdist(:) ! weights for gaussian
                                             ! injection profile

! Input parameters for Plumeria
real, allocatable, save :: wind_speed(:)
real, allocatable, save :: wind_dir(:)
real, allocatable, save :: u_theta_levels(:)
real, allocatable, save :: v_theta_levels(:)
real :: plume_height

! Variables relating to true geographic location
real, allocatable  ::  true_latitude (:,:)
real, allocatable  ::  true_longitude (:,:)

integer :: ukcavolc_unit
character(len=filenamelength) :: filename

real, allocatable, save :: emission(:) ! emission in kg(S)/timestep

logical :: trigger_calc       ! Triggers volcanic calculation on the first
                              ! timestep for the period of eruption, not on
                              ! the last timestep for the previous day
logical, save :: first = .true.

! local variables
integer :: i,j,k,l1,l2,n,ierr
character(len=256) :: format_input ! string for format descriptor for read

integer                             :: errcode=0     ! Error flag (0 = OK)
character(len=errormessagelength)   :: cmessage      ! Error return message

integer :: nrows = 0

integer(kind=jpim), parameter :: zhook_in  = 0
integer(kind=jpim), parameter :: zhook_out = 1
real(kind=jprb)               :: zhook_handle

character(len=*), parameter :: RoutineName='UKCA_VOLCANIC_SO2'

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

! On first call, determine locations on the UM grid where volcanos are located.
! The method assumes that the grid is a regular latitude-longitude grid with
! latitude constant on each row and equally-spaced between rows and
! longitude constant on each column and equally-spaced between columns.
if (first) then
  if ( .not. allocated(volc_emissions))                                        &
     allocate(volc_emissions(nparams,nerupt))
  if ( .not. allocated(lat_volc)) allocate(lat_volc(nerupt))
  if ( .not. allocated(lon_volc)) allocate(lon_volc(nerupt))
  if ( .not. allocated(year_volc)) allocate(year_volc(nerupt))
  if ( .not. allocated(start_volc)) allocate(start_volc(nerupt))
  if ( .not. allocated(end_volc)) allocate(end_volc(nerupt))
  if ( .not. allocated(magnitude)) allocate(magnitude(nerupt))
  if ( .not. allocated(height_center)) allocate(height_center(nerupt))
  if ( .not. allocated(vent_height)) allocate(vent_height(nerupt))
  if ( .not. allocated(height_top_gaussian))                                   &
     allocate(height_top_gaussian(nerupt))
  if ( .not. allocated(height_bot_gaussian))                                   &
     allocate(height_bot_gaussian(nerupt))
  if ( .not. allocated(plume_thickness)) allocate(plume_thickness(nerupt))
  if ( .not. allocated(mer)) allocate(mer(nerupt))
  if ( .not. allocated(lon_index)) allocate(lon_index(nerupt))
  if ( .not. allocated(lat_index)) allocate(lat_index(nerupt))
  if ( .not. allocated(lat_volc_rad)) allocate(lat_volc_rad(nerupt))
  if ( .not. allocated(lon_volc_rad)) allocate(lon_volc_rad(nerupt))
  if ( .not. allocated(emission)) allocate(emission(nerupt))

  filename = trim(adjustl(file_volc_so2))

  ! Open the file
  call assign_file_unit(filename, ukcavolc_unit, handler="fortran")
  open(ukcavolc_unit,file=filename,action='read')

  ! Format descriptors for explosive volcanic emissions input file
  ! The input data must follow the specified format descriptors
  ! See ctldata directory for sample input file
  format_input='(2(F7.3,1X),F4.0,1X,2(F3.0,1X),F14.1,1X,2(F7.1,1X),F12.1)'

  ! Check that nerupt equals to number of rows (eruptions) in input file
  read(unit=ukcavolc_unit,fmt='(A)') ! header row
  read_loop: do
    read(unit=ukcavolc_unit,fmt=format_input, iostat=ierr)
    if (ierr < 0) exit read_loop
    nrows = nrows + 1
  end do read_loop

  if (nerupt /= nrows) then
    cmessage = 'No of eruptions specified by nerupt /= no of eruptions in file'
    errcode = 1
    call ereport(RoutineName,errcode,cmessage)
  end if

  ! Read the input data into volc_emissions
  rewind ukcavolc_unit
  read(unit=ukcavolc_unit,fmt='(A)') ! header row
  do n = 1, nrows
    read(unit=ukcavolc_unit, fmt=format_input) volc_emissions(:,n)
  end do

  close(ukcavolc_unit)
  call release_file_unit(ukcavolc_unit, handler="fortran")

  ! Assign the columns of volc_emissions to each parameter
  lat_volc(:) = volc_emissions(1,:)
  lon_volc(:) = volc_emissions(2,:)
  year_volc(:) = int(volc_emissions(3,:))
  start_volc(:) = int(volc_emissions(4,:))
  end_volc(:) = int(volc_emissions(5,:))
  magnitude(:) = volc_emissions(6,:)
  vent_height(:) = volc_emissions(7,:)
  height_center(:)= volc_emissions(8,:)
  mer(:) = volc_emissions(9,:)

  ! each element in array is set to -1
  lon_index(:) = -1
  lat_index(:) = -1

  ! Check for compatible UM grid
  if (.not. (                                                                  &
       (model_type == mt_global .or.                                           &
         (model_type == mt_cyclic_lam .and. .not. l_cartesian)) .and.          &
       l_regular)) then
    cmessage = 'Cannot locate volcanos - model does not have a compatible grid'
    errcode  = 1
    call ereport (RoutineName, errcode, cmessage)
  end if

  ! revert from degrees to radians for latitude and longitude of volcanoes
  lat_volc_rad = lat_volc * pi_over_180
  lon_volc_rad = lon_volc * pi_over_180

  ! Find closest grid point on the UM grid.
  ! delta_phi and delta_lambda give the grid spacing for longitude
  ! and latitude respectively.
  do k=1,nerupt
    do i=1,row_length
      do j=1,rows
        if ((abs(true_latitude(i,j)-lat_volc_rad(k))<0.5*delta_phi) .and.      &
            (abs(true_longitude(i,j)-lon_volc_rad(k))<0.5*delta_lambda)) then
          lat_index(k) = j
          lon_index(k) = i
        end if
      end do
    end do
  end do

  emission = magnitude / (end_volc - start_volc + 1.0)                         &
               * timestep / rsec_per_day
  ! above is as applied to UKCA SO2 which has units of kgSO2/kgair

  first = .false.
end if !first

! Initialise the plume height from Plumeria with filled values -999
plumeria_height = -999.0

do k=1,nerupt
  ! check whether in processor domain of eruptive volcano
  if ((lat_index(k) > 0) .and. (year == year_volc(k))) then
    trigger_calc = .false.
    ! check whether in eruptive phase
    if ( i_day_number == start_volc(k) .and.                                   &
         i_hour + i_minute + i_second > 0 ) then
      trigger_calc = .true.
    end if
    if ( trigger_calc .or. (i_day_number > start_volc(k) .and.                 &
      i_day_number <= end_volc(k)) ) then
      i = lon_index(k)
      j = lat_index(k)

      ! If calling plumeria
      if (l_ukca_so2ems_plumeria) then
        if ( .not. allocated(u_theta_levels))                                  &
            allocate(u_theta_levels(model_levels))
        if ( .not. allocated(v_theta_levels))                                  &
            allocate(v_theta_levels(model_levels))
        if ( .not. allocated(wind_speed))allocate(wind_speed(model_levels))
        if ( .not. allocated(wind_dir))allocate(wind_dir(model_levels))

        ! Regrid u_rho_levels and v_rho_levels to theta levels
        do n=1,model_levels-1
          u_theta_levels(n) = (u_rho_levels(i,j,n) + u_rho_levels(i,j,n+1))/2.0
          v_theta_levels(n) = (v_rho_levels(i,j,n) + v_rho_levels(i,j,n+1))/2.0
        end do
        wind_speed = sqrt(u_theta_levels**2+v_theta_levels**2)
        wind_dir = atan2(u_theta_levels/wind_speed, v_theta_levels/wind_speed)
        wind_dir = 90-(wind_dir*180/pi + 180)

        ! Call Plumeria
        call plumeria_main(model_levels,                                       &
                        wind_speed,                                            &
                        wind_dir,                                              &
                        t_theta_levels(i,j,:),                                 &
                        p_theta_levels(i,j,:),                                 &
                        geopH_on_theta_mlevs(i,j,:),                           &
                        rel_humid_frac(i,j,:),                                 &
                        vent_height(k),                                        &
                        mer(k),                                                &
                        plume_height)

        height_center(k) = plume_height*1000+vent_height(k)
        plumeria_height(i,j) = height_center(k)

      end if

      ! Gaussian thickness of the plume, parameterized on the basis of
      ! 3D plume model simulations (Aubry et al. 2019)
      plume_thickness(k) = 0.108*(height_center(k)-vent_height(k))

      ! Plume top and bottom: limit injection to 4 gaussian widths
      height_top_gaussian(k) = height_center(k) + 4*plume_thickness(k)
      height_bot_gaussian(k) = height_center(k) - 4*plume_thickness(k)

      l1 = 1
      do while (r_theta_levels(i,j,l1) - planet_radius <=                      &
                height_bot_gaussian(k))
        l1 = l1 + 1
      end do
      l2 = model_levels
      do while (r_theta_levels(i,j,l2) - planet_radius >                       &
                height_top_gaussian(k))
        l2 = l2 - 1
      end do
      if (l2 < l1) l1 = l2
      ! express emission as increase in average increase in mass mixing ratio
      weight_volc_vertdist = (abs(r_theta_levels(i,j,l1+1:l2+1) -              &
                             r_theta_levels(i,j,l1-1:l2-1))/2) *               &
                             exp(-((r_theta_levels(i,j,l1:l2) -                &
                             planet_radius-height_center(k)) /                 &
                             plume_thickness(k))**2)
      so2_mmr(i,j,l1:l2) = so2_mmr(i,j,l1:l2) + emission(k) *                  &
                           (weight_volc_vertdist/sum(weight_volc_vertdist)) /  &
                           mass(i,j,l1:l2)
    end if
  end if
end do

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine ukca_volcanic_so2
end module ukca_volcanic_so2_mod
