! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file COPYRIGHT.txt
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
! A module to contain information about submodels

module Submodel_Mod

! Description:
!   Data module to contain information about submodels.
!   Largely unnecessary but needed for consistency across UM.
!
! Code Owner: Please refer to the UM file CodeOwners.txt
! This file belongs in section: top_level
!
! Code Description:
!   Language: FORTRAN 90
!   This code is written to UMDP 003 programming standards.

implicit none

!  1. Maximum internal model/submodel array sizes for this version.
! Max no. of internal models
integer, parameter  ::  N_Internal_Model_Max = 1

! Max no. of subm. dump parts
integer, parameter  ::  N_Submodel_Partition_Max = 1

! Max value of int. model id
integer, parameter  ::  Internal_Id_Max = N_Internal_Model_Max

! Max value of subm. dump id
integer, parameter  ::  Submodel_Id_Max = N_Submodel_Partition_Max


!  2. Atmos Internal Model identifiers
integer, parameter  ::  atmos_im = 1


! Atmos Submodel partition identifiers
integer, parameter  ::  atmos_sm = 1

!  3. Lists of internal models and their submodel dump partitions -
integer             :: N_Internal_Model  = 1    ! No. of internal models
                     ! For all forecast jobs this is 1 (atmos).
                     ! Not able to set as parameter due to use in utilities
                     ! and SCM.
integer, save  :: Internal_Model_List(N_Internal_Model_Max)

! Submodel identifier for each internal model in list
integer, save       :: Submodel_For_IM(N_Internal_Model_Max)

! Submodel number for each submodel id
integer, save       :: Submodel_For_SM(N_Internal_Model_Max)

!  4. Lists calculated in model from user interface supplied arrays -
!     - experiment specific.
! No of internal models in  each submodel partition indexed by sm identifier
integer, save  :: N_Internal_For_SM(Submodel_ID_Max)

! List of submodel partition identifiers
integer, save  :: Submodel_Partition_List(N_Submodel_Partition_Max)

! Submodel partition identifier indexed by internal model identifier
integer, save  :: Submodel_Partition_Index(Internal_ID_Max)

end module Submodel_Mod
