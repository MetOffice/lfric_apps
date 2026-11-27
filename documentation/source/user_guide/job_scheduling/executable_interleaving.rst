.. ------------------------------------------------------------------------------
     (c) Crown copyright Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------
.. _executable_interleaving_user_index:

Executable interleaving
=======================

Interleaved configurations place multiple application executables on a single
compute node, allowing individual nodes to host a mixture of applications.

Interleaving is a more complicated job layout on the high performance computer,
it commonly requires targeted configuration using the job scheduler for the
target HPC.

How it works
------------

The launcher leverages the ``PBS_NODEFILE``, which describes the nodes allocated
to the job, to construct an ``mpiexec`` hostfile containing the required rank-to-node
placement.

By default, executables are separated into distinct blocks of compute ranks. In
a typical hostfile (not usually seen in the standard launching method), you can
imagine it to look like a long list of consecutive node ids with as many entries
for each node as ``ppn`` allows.

When ``mpiexec`` launches the MPI processes, each process is assigned an MPI rank.
The hostfile determines which node each of those ranks is launched on. By
generating the hostfile in executable-sized blocks with the required number of
ranks-per-node, ranks belonging to different executables can be placed on the
same compute nodes.

Example hostfile
----------------

Consider two applications, A and B, sharing a set of compute nodes. Application
A uses 24 MPI ranks and application B uses 12 MPI ranks. The job is configured
to place six MPI ranks on each node.

Without interleaving, each application occupies its own set of nodes. The
hostfile can be thought of as containing the following rank-to-node placement.
The MPI rank and executable are shown here for illustration; the hostfile
itself contains only the node IDs.

.. code-block:: text

    MPI ranks    executable    node

    0-5          A             node01
    6-11         A             node02
    12-17        A             node03
    18-23        A             node04

    24-29        B             node05
    30-35        B             node06

The MPMD ``mpiexec`` command determines that ranks 0-23 run application A and
ranks 24-35 run application B. The hostfile determines the node on which each
of those ranks is launched.

With interleaving, the same ranks and executables are used, but the hostfile is
rearranged so that four ranks from application A and two ranks from application
B are placed on each node:

.. code-block:: text

    MPI ranks    executable    node

    0-3          A             node01
    4-7          A             node02
    8-11         A             node03
    12-15        A             node04
    16-19        A             node05
    20-23        A             node06

    24-25        B             node01
    26-27        B             node02
    28-29        B             node03
    30-31        B             node04
    32-33        B             node05
    34-35        B             node06

The ordering of the executables in the global MPI rank space has not changed:
application A still occupies ranks 0-23 and application B still occupies ranks
24-35. Only the node associated with each rank has changed.

The resulting placement on each node is therefore:

.. code-block:: text

    node01    A ranks 0-3      B ranks 24-25
    node02    A ranks 4-7      B ranks 26-27
    node03    A ranks 8-11     B ranks 28-29
    node04    A ranks 12-15    B ranks 30-31
    node05    A ranks 16-19    B ranks 32-33
    node06    A ranks 20-23    B ranks 34-35

The generated hostfile therefore controls the placement of MPI ranks across
nodes. It does not specify which CPUs within a node those ranks may use; this
is handled separately by per-rank CPU binding.

Rank binding
------------

The generated hostfile determines which node each MPI rank is placed on, but
does not control which CPUs within that node the rank may use. Additional
per-rank CPU binding is therefore applied after placement. Each rank uses its
global MPI rank, ranks-per-node, CPU starting offset and required CPU depth to
calculate its CPU range. The resulting affinity is applied with taskset before
the requested executable is launched. This ensures that co-located executables
use their intended, non-overlapping CPU regions.

.. _interleaving_requirements:

Interleaving requirements
=========================

.. _interleaving_requirements_env_vars:

Environment variables
---------------------

The interleaving functionality of the launcher is activated and configured using
the following environment variables, all of which are required:

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
     - Use this to select from the :ref:`supported interleaving methods <interleaving_methods>`.
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

Configuration validations
-------------------------

Interleaving configurations must adhere to the following rules:

1. The total requested processors per-node for all executables must not exceed the total number of available cores per-node.

.. _exec_interleaving_rose-stem:

Rose-stem
---------

For tasks launched via the rose-stem test suite:

1. The following task-specific settings must be defined. See
   :ref:`above <interleaving_requirements_env_vars>` for detialed descriptions
   of the environment variables set by these options.

.. list-table::
   :widths: 8 8 8 25
   :header-rows: 1

   * - Key
     - Type
     - Default
     - Description
   * - interleaved
     - Boolean
     - False
     - Activates interleaving launch method. Used to set ``INTERLEAVE_EXECS`` as True.
   * - mpi_parts
     - Int
     - NA
     - The number of ranks on which to run the application. Used to set <EXEC_NAME>_TOTAL_RANKS.
   * - task_ranks_per_node
     - Int
     - NA
     - Used to set the <EXEC_NAME>_RANKS_PER_NODE.
   * - task_ranks_depth_pad
     - Int
     - NA
     - The number of CPUs to dedicate to each model process running an executable. Used to set <EXEC_NAME>_RANK_DEPTH
   * - threads
     - Int
     - NA
     - The number of OpenMP threads to run the executable using. This cannot exceed task_ranks_depth_pad.
   * - xios_ppn
     - Int
     - 1 if interleaved is True, else NA
     - The number of XIOS processes to place on each node. Used to set XIOS_SERVER_RANKS_PER_NODE.
   * - xios_depth
     - Int
     - 1 if interleaved is True, else NA
     - The number of CPUs to dedicate to each process running the xios executable. Used to set XIOS_SERVER_RANK_DEPTH.
   * - xios_nodes
     - Int
     - 0 if interleaved is True, else see :ref:`rose-stem task options <rose-stem_task_options>`
     - The number of nodes to dedicate to the XIOS server.
   * - mpi_parts_xios
     - Int
     - NA if interleaved is True, else see :ref:`rose-stem task options <rose-stem_task_options>`
     - The number of ranks to dedicate to the XIOS server. Must be greater than 0 for interleaving. Used to set XIOS_SERVER_TOTAL_RANKS.
   * - xios_threads
     - Int
     - 1 if interleaved is True.
     - The number of OpenMP threads with which to run the XIOS server. This cannot exceed task_ranks_depth_pad. USed to set XIOS_SERVER_OMP_NUM_THREADS.


2. The test must be launched for one of the following supported :ref:`sites <new_site>`:

.. code-block::text

    meto

3. The application must not be part of the ``coupled_interface``. This is a
   temporary measure while support for setups other than model + XIOS is
   implemented into rose-stem.

.. _interleaving_methods:

Interleaving methods
====================

The interleaving infrastructure currently supports the following interleaving
methods:

1. Uniform interleaving
-----------------------

A single pattern of interleaving is applied to all nodes, uniformly and regardless
of whether the reserved resources are completely filled by the number of processes.


Take the following example using two executables, A and B:

Given a total of 12 available cores per node, A is given 8 cores and
B is allocated 4, totalling 12.

.. code-block::text

    rank id   1   2   3   4   5   6   7   8 | 9   10  11  12
    exec      A   A   A   A   A   A   A   A | B   B   B   B

The job uses 2 nodes, both using this pattern. However, while executable B
requests 8 processors, executable A only requests 12. The final job will look
like this:

.. code-blaock::text

    Node 0
    rank id   1   2   3   4   5   6   7   8 | 9   10  11  12
    exec      A   A   A   A   A   A   A   A | B   B   B   B

    Node 1
    rank id   1   2   3   4   5   6   7   8 | 9   10  11  12
    exec      A   A   A   A   -   -   -   - | B   B   B   B


Performance implications of interleaving
========================================

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
    :ref:`I/O Interleaving <interleaved_server_io>` for good practices regarding
    interleaving XIOS and model executables.
