!-----------------------------------------------------------------------------
! (c) Crown copyright 2018 Met Office. All rights reserved.
! The file LICENCE, distributed with this code, contains details of the terms
! under which the code may be used.
!-----------------------------------------------------------------------------
! @brief Module setting sizes of spectral arrays

module dimensions_spec_ucf

implicit none

! These are dimensions for arrays in the Socrates radiation code that have
! not yet been transitioned to dynamic allocation. It is intended that
! eventually these will be removed entirely. For now, they are set to the
! largest value necessary for spectral files that are expected to be used
! with LFRic.

integer, parameter :: npd_k_term = 14
! Number of esft terms
integer, parameter :: npd_type = 20
! Number of data types
integer, parameter :: npd_scale_variable = 4
! Number of scaling variables
integer, parameter :: npd_drop_type = 5
! Number of drop types
integer, parameter :: npd_ice_type = 16
! Number of ice crystal types
integer, parameter :: npd_cloud_parameter = 30
! Number of cloud parameters
integer, parameter :: npd_humidities = 21
! Number of humidities

end module dimensions_spec_ucf
