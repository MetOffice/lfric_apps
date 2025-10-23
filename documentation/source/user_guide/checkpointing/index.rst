.. ------------------------------------------------------------------------------
     (c) Crown copyright Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _checkpoint_index:

Checkpoint and restart
======================

Checkpoint and restart describes the process of writing out a file containing
the model state during or at the end of one model run and then reading it back
into a new execution of the same model configuration to extend the model run
time or to recover from a technical failure.

Normally, it is a requirement that a long single run produces exactly the same
results as a short initial run plus a restarted run of the same total
length. Maintaining such equivalence is considered to be important by scientist
users.

Currently, support for checkpointing exists only for ``lfric_atm`` and
``gungho_model`` applications.

.. toctree::
    :maxdepth: 1

    lfric_atm_checkpoint
