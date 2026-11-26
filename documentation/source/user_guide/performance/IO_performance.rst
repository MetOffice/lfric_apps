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
------------

In ``Attached I/O`` mode, then the simulation or application ranks act as their own
I/O resource, using XIOS functionality to interact with the storage system.

This configuration works well for smaller work loads with a limited number of ranks,
perhaps running as a ``shared queue`` job.

In ``Attached I/O`` mode, then I/O operations are blocking, application code is limited
in how fast it can run and synchronise by the I/O block patterns.

This configuration often does not scale well with job size and complexity.

There is limited opportunity for I/O performance tuning.

Implementing Attached I/O
+++++++++++++++++++++++++

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
+++++++++++++++++++++

``Segregated Server I/O`` configurations place application ranks on separate nodes from I/O server ranks.  There are nodes which are dedicated to only running ``xios_server.exe`` processes.

Implementing Segregated Server I/O
##################################

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


Interleaved Server I/O
++++++++++++++++++++++

``Interleaved Server I/O`` configurations places application and I/O server ranks across nodes, with individual nodes hosting a mixture of application and I/O processes.

Interleaving is a more complicated job layout on the high performance computer, it commonly requires targeted configuration using the job scheduler for the target HPC.

Interleaving can provide a further I/O performance benefit.  This is due to two different
performance limitations:

* MPI Network communication
    * Messages between application ranks and I/O ranks use the interconnect and
      networking implementation on the HPC.  Interleaving reduces the node by node
      messaging load, with fewer data offload messages being sent to each node in the
      configuration (compared to a small number of segregated  I/O nodes, which must
      receive all of the data offloads).

* Storage Throughput Parallelisation
    * Commonly individual HPC nodes will have limits on how much each node may parallelise
      I/O.  Where more nodes are able to interact with the storage system directly, then
      more I/O throughput can be achieved.

Implementing Interleaved Server I/O
###################################

``lfric_apps`` provides site specific launch scripts, for example
`meto/launch-exe <https://github.com/MetOffice/lfric_apps/blob/main/rose-stem/site/meto/common/bin/launch-exe>`_
which configure the launch command based on the controlling environment variables.

Suites and tests may use this launch script and provide environment variables
to configure the ``mpiexec`` command.

Setting Environment variables for ``Segregated Server I/O`` mode:

* ``XIOS_SERVER_MODE=True``


Launch Commands and Job Schedulers
----------------------------------

A scheduler will manage work being run on HPC platforms, including:

* targeting queues,
* requesting cores and nodes,
* requesting node types.

A launch command is responsible for defining the job layout, including facets such as:

* ``--ppn`` Processes/ Ranks per compute node.
* ``-n``  Number of processes to start
* ``--depth`` CPUs per process / rank. e.g. Run <application> with two CPUs per process.

(see ``man mpiexec`` for more information)

The layout must be consistent with the job definition, particularly:

* the number of nodes being requested
* the number of cores per node being compatible with:
    * the processes / ranks per compute node
    * the number of ``threads`` for ``OMP``
    * the CPUs per process / rank


Storage System
--------------

Performance characteristics of the storage system being used have large impacts on I/O performance.  Implementation specific storage system interaction tuning is a valuable optimisation surface.

Lustre Striping
+++++++++++++++

Lustre is a commonly used storage provision in high performance computing.

Lustre provides a key performance optimisation tool:
`Lustre Striping <https://wiki.lustre.org/index.php/Configuring_Lustre_File_Striping>`_
.

Well tuned Lustre stripe configuration that is complementary to the I/O characteristics
of the application shows excellent return in performance measurements.
