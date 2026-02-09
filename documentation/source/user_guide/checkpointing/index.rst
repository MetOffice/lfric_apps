.. ------------------------------------------------------------------------------
     (c) Crown copyright Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------

.. _checkpoint_index:

Checkpoint and restart
======================

Checkpointing describes the process of writing out a checkpoint dump file
containing key application data during or at the end of an application run. It
should then be possible to restart the application with the same configuration,
but configured to start from the newly created checkpoint file. Checkpoint and
restart can be used to support very long runs such as those used for climate
projections, or to recover a run from an intermediate point following an
application failure.

Normally, it is a requirement that a long single run produces exactly the same
results as a short initial run plus a restarted run of the same total
length. Maintaining such equivalence is important to ensure reproducibility of
results regardless of choices made about when a run should write or read a
checkpoint dump.

Currently, support for checkpointing exists only for the ``lfric_atm`` and
``gungho_model`` applications.

.. toctree::
    :maxdepth: 1

    lfric_atm_checkpoint
