module tcs_cmt_params_sh

use um_types, only: real_umphys

implicit none

! Description:
!  Contains shallow CMT settings for turbulence CMT scheme
!
! (c) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! ------------------------------------------------------------------------------
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: Convection

! Parameters for calculation of cloud base stress

real(kind=real_umphys), parameter ::                                           &
  gamma_cmt_shall = 1.0/2.3       & ! Value from fit to CRM results
 ,delta_cmt_shall = (1.0-1.63/2.3)  ! value from fit to CRM results


end module tcs_cmt_params_sh
