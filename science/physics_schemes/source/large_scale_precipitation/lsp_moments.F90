! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
!  Large-scale precipitation scheme. Conversion between moments of PSD

! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: large_scale_precipitation

!Note regarding variable precision:
!This module is used in various places beyond the lsp scheme which may require
!either 64 or 32 bit calculation. Therefore, we cannot use the real_lsprec
!parameter as seen generally in the lsp_* modules. Instead, we create an
!interface with both 32 and 64 bit versions available.

module lsp_moments_mod

use um_types, only: real_64, real_32

implicit none

private
public ::  lsp_moments

interface lsp_moments
  module procedure lsp_moments_64b, lsp_moments_32b
end interface

character(len=*), parameter, private :: ModuleName='LSP_MOMENTS_MOD'

contains

subroutine lsp_moments_64b( points, rho, t, qcf, cficei, n_out, moment_out )

  ! microphysics modules
use mphys_inputs_mod,    only: ai_mod => ai, bi_mod => bi

! General atmosphere modules
use conversions_mod,  only: zerodegc

! Dr Hook modules
use yomhook,          only: lhook, dr_hook
use parkind1,         only: jprb, jpim
use vectlib_mod,      only: powr_v    => powr_v_interface,                     &
                            oneover_v => oneover_v_interface,                  &
                            exp_v     => exp_v_interface

implicit none
character(len=*), parameter :: RoutineName='LSP_MOMENTS_64B'
integer, parameter :: prec = real_64
#include "lsp_moments.h"
return
end subroutine lsp_moments_64b

!==============================================================================

subroutine lsp_moments_32b( points, rho, t, qcf, cficei, n_out, moment_out )

  ! microphysics modules
use mphys_inputs_mod,    only: ai_mod => ai, bi_mod => bi

! General atmosphere modules
use conversions_mod,  only: zerodegc => zerodegc_32

! Dr Hook modules
use yomhook,          only: lhook, dr_hook
use parkind1,         only: jprb, jpim
use vectlib_mod,      only: powr_v    => powr_v_interface,                     &
                            oneover_v => oneover_v_interface,                  &
                            exp_v     => exp_v_interface

implicit none
character(len=*), parameter :: RoutineName='LSP_MOMENTS_32B'
integer, parameter :: prec = real_32
#include "lsp_moments.h"
return
end subroutine lsp_moments_32b

end module lsp_moments_mod
