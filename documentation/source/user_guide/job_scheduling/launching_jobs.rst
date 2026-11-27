
Launch Commands and Job Schedulers
++++++++++++++++++++++++++++++++++

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
