.. -----------------------------------------------------------------------------
     (c) Crown copyright 2025 Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   -----------------------------------------------------------------------------

*********
Transport
*********
This is an application that allows transport-only tests with any of the
transport methods to be run without having to use the full Gungho Dynamics
model.

A specified wind field is used to transport all fields. A number of fields are
transported:

#. A W3 field in a conservative (flux) form equation (equivalent to density in
   the full model)
#. A W3 field in a non-conservative (advective) form equation (equivalent to the
   wind components in the full model)
#. A Wtheta field in a non-conservative (advective) form equation (equivalent to
   the potential temperature in the full model)
#. A Wtheta field in a conservative (flux) form equation as a mixing ratio
   (equivalent to moisture in the full model)
