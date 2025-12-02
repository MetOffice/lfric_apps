! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: top_level
!
! Purpose: Read in and hold information related to model domain. This
!          relates to the domain properties should not be confused
!          with grid properties, i.e. number of rows/columns etc.
!

module model_domain_mod

use missing_data_mod,       only: rmdi, imdi
use errormessagelength_mod, only: errormessagelength
use yomhook,  only: lhook, dr_hook
use parkind1, only: jprb, jpim

implicit none

!-----------------------------------------------------
! Parameters
!-----------------------------------------------------
! Possible values for namelist variable - model_type
! 0 => Small execs don't use this, so we have a placeholder value
! 1 => Global Model (i.e. cyclic EW, and NS overpole - a sphere)
! 2 => Limited Area Model (Classic style)
! 3 => Limited Area Model (Cyclic boundary conditions - EW only)
! 4 => Limited Area Model (Cyclic boundary conditions - EW and NS, (i.e. torus)
! 5 => Single Column Model (i.e. one horizontal point)

integer, parameter :: mt_smexe         = 0
integer, parameter :: mt_global        = 1
integer, parameter :: mt_lam           = 2
integer, parameter :: mt_cyclic_lam    = 3
integer, parameter :: mt_bi_cyclic_lam = 4
integer, parameter :: mt_single_column = 5
integer, parameter :: mt_lfric         = 10

! DrHook-related parameters
integer(kind=jpim), parameter, private :: zhook_in  = 0
integer(kind=jpim), parameter, private :: zhook_out = 1


! Item not in the model domain namelist
! It is read in and set up by the reconfiguration within the horizont namelist.
integer :: output_grid_stagger = 6

integer,parameter :: FH_GridStagger_C         = 3
integer,parameter :: FH_GridStagger_Endgame   = 6

!-----------------------------------------------------
! Items set in namelist
!-----------------------------------------------------

integer :: model_type = imdi
logical :: l_regular  = .false. ! true if not variable resolution
logical :: l_cartesian = .false.   ! True = Cartesian co-ordinates

namelist /model_domain/ model_type, l_regular, l_cartesian


character(len=*), parameter, private :: ModuleName='MODEL_DOMAIN_MOD'

contains

subroutine check_nml_model_domain()

! Description:
!   Subroutine to apply checks based on options read in model_domain namelist.

use chk_opts_mod, only: chk_var, def_src

implicit none

character (len=*), parameter       :: RoutineName = 'CHECK_NML_MODEL_DOMAIN'
real(kind=jprb)                    :: zhook_handle


if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)
def_src = RoutineName

! Check available model domain type is chosen
call chk_var( model_type, 'model_type'                                         &
            , [mt_global, mt_lam, mt_cyclic_lam, mt_bi_cyclic_lam,             &
               mt_single_column] )

def_src = ''
if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)
return
end subroutine check_nml_model_domain


subroutine print_nlist_model_domain()

use umPrintMgr, only: umPrint

implicit none

character(len=50000) :: lineBuffer
character(len=*), parameter    :: RoutineName = 'PRINT_NLIST_MODEL_DOMAIN'
real(kind=jprb)      :: zhook_handle

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_in,zhook_handle)

call umPrint('Contents of namelist model_domain', src='model_domain_mod')

write(lineBuffer,'(A,I2)') ' model_type = ', model_type
call umPrint(lineBuffer, src='model_domain_mod')
write(lineBuffer,'(A,L1)') ' l_regular  = ', l_regular
call umPrint(lineBuffer, src='model_domain_mod')
write(lineBuffer,'(A,L1)') ' l_cartesian  = ', l_cartesian
call umPrint(lineBuffer, src='model_domain_mod')
call umPrint('- - - - - - end of namelist - - - - - -', src='model_domain_mod')

if (lhook) call dr_hook(ModuleName//':'//RoutineName,zhook_out,zhook_handle)

end subroutine print_nlist_model_domain

end module model_domain_mod
