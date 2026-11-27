.. ------------------------------------------------------------------------------
     (c) Crown copyright Met Office. All rights reserved.
     The file LICENCE, distributed with this code, contains details of the terms
     under which the code may be used.
   ------------------------------------------------------------------------------
.. _IO_performance_user_index:

I/O Performance
===============

Applications use XIOS, the XML I/O Server, to manage the majority of data input and
all data output. XIOS is maintained at
`<https://gitlab.in2p3.fr/ipsl/projets/xios-projects/xios>`_ with published
`documentation <https://ipsl.pages.in2p3.fr/projets/xios-projects/xios/>`_.

``lfric_apps`` applications use XIOS in multiple configurations.

Attached I/O
============

In ``Attached I/O`` mode, then the simulation or application ranks act as their own
I/O resource, using XIOS functionality to interact with the storage system.

This configuration works well for smaller work loads with a limited number of ranks,
perhaps running as a ``shared queue`` job.

In ``Attached I/O`` mode, then I/O operations are blocking, application code is limited
in how fast it can run and synchronise by the I/O block patterns.

This configuration often does not scale well with job size and complexity.

There is limited opportunity for I/O performance tuning.

Implementing Attached I/O
-------------------------

Launching applications in SPMD (Single Program Multiple Data) mode, with an
``mpiexec`` command (or similar) that only references the applications executable will
run any XIOS code in ``Attached I/O`` mode.

``lfric_apps`` provides site specific launch scripts, for example
`meto/launch-exe <https://github.com/MetOffice/lfric_apps/blob/main/rose-stem/site/meto/common/bin/launch-exe>`_
which configure the launch command based on the controlling environment variables.

Suites and tests may use this launch script and provide environment variables
to configure the ``mpiexec`` command.

Setting Environment variables for ``Attached I/O`` mode:

* ``XIOS_SERVER_MODE=False``

Server I/O
----------

In ``Server I/O`` mode, then the simulation or application ranks are run alongside
a set of I/O ranks, which are running ``xios_server.exe``, using ``MPI``  ``MPMD``
(Multiple programs, multiple data streams).

This configuration is significantly more scalable than ``Attached I/O`` mode.
It is commonly used for larger data volumes and for ``normal queue`` jobs running
on multiple HPC nodes.

These configurations are significantly more scalable than ``Attached I/O``.
In particular, application code is able to offload data to separate I/O processes,
allowing the application to continue activity with minimal impact from the storage
system's performance (as long as the I/O can continue fast enough to keep buffers
available to offload data).



Segregated Server I/O
=====================

``Segregated Server I/O`` configurations place application ranks on separate nodes from I/O server ranks.  There are nodes which are dedicated to only running ``xios_server.exe`` processes.

Implementing Segregated Server I/O
----------------------------------

Launching applications in MPMD (Multiple Programs, Multiple Data streams) mode, with an
``mpiexec`` command (or similar) that references the applications executable and the
``xios_server.exe`` application will run the application in
``Segregated Server I/O`` mode.

``lfric_apps`` provides site specific launch scripts, for example
`meto/launch-exe <https://github.com/MetOffice/lfric_apps/blob/main/rose-stem/site/meto/common/bin/launch-exe>`_
which configure the launch command based on the controlling environment variables.

Suites and tests may use this launch script and provide environment variables
to configure the ``mpiexec`` command.

Setting Environment variables for ``Segregated Server I/O`` mode:

* ``XIOS_SERVER_MODE=True``

.. _interleaved_server_io:

Interleaved Server I/O
======================

``Interleaved Server I/O`` configurations place application and I/O server
ranks across nodes, with individual nodes hosting a mixture of application and
I/O processes.

Interleaving is a more complicated job layout on the high performance computer,
it commonly requires targeted configuration using the job scheduler for the
target HPC.

Interleaving can provide a further I/O performance benefit.  This is due to two
different performance limitations:

* MPI Network communication
    * Messages between application ranks and I/O ranks use the interconnect and
      networking implementation on the HPC.  Interleaving reduces the node-by-node
      messaging load, with fewer data offload messages being sent to each node in the
      configuration (compared to a small number of segregated  I/O nodes, which must
      receive all of the data offloads).

* Storage Throughput Parallelisation
    * Commonly individual HPC nodes will have limits on how much each node may parallelise
      I/O.  Where more nodes are able to interact with the storage system directly, then
      more I/O throughput can be achieved.


Implementing Interleaved Server I/O
-----------------------------------

``lfric_apps`` provides site specific launch scripts, for example
`meto/launch-exe <https://github.com/MetOffice/lfric_apps/blob/main/rose-stem/site/meto/common/bin/launch-exe>`_
which configure the launch command based on the controlling environment variables.

Suites and tests may use this launch script and provide environment variables
to configure the ``mpiexec`` command.

It is recommended to dedicate an entire NUMA region to I/O server executables
(though not every CPU must be a server) to preserve memory locality on the
model-side.

Suites
------

For suites, the required environment variables for a given job can be found in
the :ref:`interleaving requirements <interleaving_requirements>` section of
the :ref:`executable interleaving <executable_interleaving_user_index>` page.

Rose-stem
---------

For tasks wishing to use interleaving, see :ref:`these required settings <exec_interleaving_rose-stem>`.
As an example, the task definition might add this for an ``lfric_atm`` task:

.. code-block::text

    "interleaved": true,
    "mpi_parts": 2352,
    "threads": 2,
    "task_ranks_depth_pad": 2,
    "task_ranks_per_node": 168,
    "use_xios": true,
    "mpi_parts_xios": 168,
    "xios_ppn": 6,
    "xios_depth": 4,


Storage System
==============

Performance characteristics of the storage system being used have large impacts
on I/O performance.  Implementation specific storage system interaction tuning
is a valuable optimisation surface.

Lustre Striping
---------------

Lustre is a commonly used storage provision in high performance computing.

Lustre provides a key performance optimisation tool:
`Lustre Striping <https://wiki.lustre.org/index.php/Configuring_Lustre_File_Striping>`_
.

Well tuned Lustre stripe configuration that is complementary to the I/O characteristics
of the application shows excellent return in performance measurements.
