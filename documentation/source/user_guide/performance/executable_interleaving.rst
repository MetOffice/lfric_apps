.. ------------------------------------------------------------------------------
     (c) Crown copyright Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------
.. _executable_interleaving_user_index:

Interleaved Executables
+++++++++++++++++++++++

Interleaved configurations place multiple application executables on a single
compute node, allowing individual nodes to host a mixture of applications.

Interleaving is a more complicated job layout on the high performance computer,
it commonly requires targeted configuration using the job scheduler for the
target HPC.

Interleaving requirements
#########################

The interleaving functionaility of the launcher is activated and configured using
the following environment variables:

.. list-table::
   :widths: 8 8 8 25
   :header-rows: 1

   * - Key
     - Type
     - Default
     - Description
   * - INTERLEAVE_EXECS
     - Bool
     - False
     - When set to True, diverts launcher to use interleaving
   * - INTERLEAVING_METHOD
     - String
     - "uniform"
     - Use this to select from the :ref:`supported interleaving methods <Interleaving methods>`.
   * - EXEC_LABELS
     - String
     - NA
     - A space-separated list of executables to interleave in MPMD order

And for each executable defined in the ``EXEC_LABELS`` environment variable, the
following variables must be defined (all uppercase):

.. list-table::
   :widths: 8 8 8 25
   :header-rows: 1

   * - Key
     - Type
     - Default
     - Description
   * - <EXEC_NAME>_TOTAL_RANKS
     - Int
     - NA
     - The total number of ranks the executable will use.
   * - <EXEC_NAME>_RANKS_PER_NODE
     - Int
     - NA
     - The number of MPI ranks used by the executable on each node. Note that this is not the number of CPUs, it is the number of processors actually running the executable.
   * - <EXEC_NAME>_RANK_DEPTH
     - Int
     - NA
     - The number of CPUs to dedicate to each process running the executable.
   * - <EXEC_NAME>_OMP_NUM_THREADS
     - Int
     - NA
     - The number of OpenMP threads to run the executable using. This cannot exceed <EXEC_NAME>_RANK_DEPTH.
   * - <EXEC_NAME>_PATH
     - String
     - NA
     - The path to the executable.


Interleaving methods
####################

The interleaving infrastructure currently supports the following interleaving
methods:

1. Uniform interleaving
-----------------------

A single pattern of interleaving is applied to all nodes, uniformly.



Performance implications of interleaving
########################################

Interleaving executables can provide a further performance benefit.  This is due
to two different performance limitations:

* MPI Network communication
    * Messages between application processes placed on distinct nodes must use
      the network interconnect and networking implementation on the HPC.
      Interleaving *can* reduce the node-to-node messaging load when ranks with
      high levels of communication are located on the same node. For example,
      I/O server ranks relocated to their clients' nodes can
      significantly reduce messages on the network interconnect (compared to a
      set of segregated  I/O nodes, which must receive all of the
      data offloads via the network).

.. warning::
    Interleaving can increase node-to-node messaging, especially when interleaving
    executables that have low inter-application comms and high intra-application
    comms.

* Storage Throughput Parallelisation
    * Commonly individual HPC nodes will have limits on how much each node may
      parallelise I/O operations. Where more nodes are able to interact with the
      storage system directly, then more I/O throughput can be achieved.

.. seealso::
    :ref:`I/O Interleaving <IO_performance>` for good practices regarding
    interleaving XIOS and model executables.




