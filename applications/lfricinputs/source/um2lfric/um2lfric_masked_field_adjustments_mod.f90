! *****************************COPYRIGHT*******************************
! (C) Crown copyright Met Office. All rights reserved.
! For further details please refer to the file LICENCE
! which you should have received as part of this distribution.
! *****************************COPYRIGHT*******************************
MODULE um2lfric_masked_field_adjustments_mod

USE lfricinp_masked_field_adjust_type_mod, ONLY:                               &
                                              lfricinp_masked_field_adjust_type

IMPLICIT NONE

PRIVATE

TYPE(lfricinp_masked_field_adjust_type), PUBLIC  :: land_field_adjustments,    &
                                                    maritime_field_adjustments

END MODULE um2lfric_masked_field_adjustments_mod
