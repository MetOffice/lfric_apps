.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

.. _lfric_atm:

The LFRic atmosphere model - lfric_atm
======================================

The `lfric_atm` application is the LFRic atmosphere model used within the
|Momentum| system. It is built using the LFRic core infrastructure
code. It runs the Gungho dynamical core coupled to several physical
parametrizations. Some parametrizations are included in the LFRic apps code
base, others such as Jules, Socrates and UKCA are extracted from external
repositories.

For IO, the application uses the XIOS library.

The checkpoint/restart system is described :ref:`here <lfric_atm_checkpoint>`.
