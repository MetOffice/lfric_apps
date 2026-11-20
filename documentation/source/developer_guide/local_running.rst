.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

Running LFRic Apps from the Command Line
========================================

Introduction
------------

When debugging application codes, it is often useful to be able to run a
particular configuration outside of the suite. The advantages are:

#. The job will often run much more quickly because it doesn't need to set up 
   the suite and create the environment in which to run it.
#. Although suites are very configurable, so can often be cut down to just the 
   run you are interested in, running from the command line allows you to run an
   even smaller part of the suite, with pre-built meshes and configuration     
   files.
#. As discussed above (in :ref:`command line builds <command_line_builds>`), it
   is possible to perform incremental builds. After the first full build, all
   subsequent builds are  very much faster.

This leads to a quicker turnaround and considerably shortens the debugging
cycle:

#. -> look at your results
#. -> make a change to the code to test something
#. -> recompile
#. -> re-run
#. -> look at your results
#. ... etc.

Run from the command line
-------------------------

To run a configuration from the command line, you will need to first run that
configuration from within a suite. Then you will need:

1. *An executable to run*. To get your own local build of an executable, see
   these :ref:`instructions <command_line_builds>`, or take a copy of the
   executable that the suite was run with from under the ``output/bin``
   directory within the ``$HOME/cylc-run`` directory produced by the suite.

.. note::

   Note that if if you are rebuilding the application  and you are trying to 
   replicate a run from the suite, you will have to use the same PSyclone 
   optimisation scripts.

2. *A configuration to run*. Take a copy of the directory that contains the
   configuration you want to run from under the ``work`` directory in the
   ``$HOME/cylc-run`` directory. Copy it into the root directory of the
   application in your working copy. You can give it any name you want to.

``cd`` into the new directory that contains your configuration and type

``../bin/<application> configuration.nml``

where ``<application>`` is the name of your application executable (e.g.
``lfric_atm`` or ``gungho_model``).

.. note::

 The directory you've copied probably contains links back into your cylc-run
 directory - so if housekeeping runs, or you ``cylc clean`` the suite,  those
 files will disappear. Also note that the iodef.xml file includes other xml
 files held in the cylc-run directory. The configuration.nml file may also
 contain references to files in the cylc-run directory. For example, the
 checkpoint stem name usually points to a location in the share directory under
 cylc-run. If you want to run locally over a longer period, consider taking
 copies of any files from cylc-run that are needed by your configuration, so
 they won't just disappear at some point.

Parallel, multi-core running
----------------------------

Many jobs in the suite run on multiple cores. Your desktop probably doesn't have
enough cores to be able to run them directly from the command line. If your job
needs more cores than your desktop provides, it will have to be run on a bigger
computer that may have a resource management system for managing batch queues,
such as Slurm or PBS. You will, therefore, need a script you can submit. For
example:

.. tab-set::

    .. tab-item:: Slurm

        .. code-block:: shell


            #!/bin/bash -l
            #SBATCH --job-name=lfric
            #SBATCH --output=job.out
            #SBATCH --error=job.err
            #SBATCH --time=5:00
            #SBATCH --export=NONE
            #SBATCH --mem-per-cpu=3G
            #SBATCH --ntasks=6
            #SBATCH --cpus-per-task=1

            module purge
            module use <path_to_modules>
            module load lfric

            export OMP_NUM_THREADS=1

            ulimit -s unlimited

            mpirun -n 6 ../bin/lfric_atm configuration.nml

        You can change the values of the variables in the slurm header. The
        values used in the suite you are replicating can be found at the top of
        the ``job`` file in the log directory under cylc-run.

        Useful slurm commands:

        #. Submit with:   ``sbatch <slurm_file>``
        #. Monitor with:  ``squeue -u <userid>``
        #. Kill with:     ``scancel <jobid>``

        .. note::

          Certainly within the Met Office (but it might also be true elsewhere):

          If your working copy is being held in a directory on a file system
          that is local to your desktop (as this is where builds will run
          fastest to keep the debug cycle short), you may not be able to run
          directly on Spice (Spice can't see your desktop's local filing
          system). Take a copy of the job onto somewhere Spice can see (like a
          data directory) and once you've rebuilt on your local filing system,
          copy the executable over to your data directory before submitting the
          job

        Once your job has run on Spice, All the output will be written to the
        directory you ran from (the directory that contains your configuration).
        You will find all log messages written to the PET files. The output and
        error messages will be written to ``job.out`` an ``job.err``
        respectively. Any XIOS errors will be written to  XIOS ``.out`` and
        ``.err`` files.

    .. tab-item:: PBS

        .. code-block:: shell

            #!/bin/bash
            #PBS -S /bin/bash
            #PBS -N lfric
            #PBS -o job.out
            #PBS -o job.err
            #PBS -q shared
            #PBS -l select=1:ncpus=6:mem=3gb
            #PBS -l walltime=5:00

            module use <path_to_modules>
            module switch PrgEnv-cray PrgEnv-cray/8.4.0
            module load cpe/23.05
            module switch cce cce/15.0.0
            module load lfric-cray/15.0.0/3.2

            export OMP_NUM_THREADS=1

            ulimit -s unlimited -c unlimited

            mpirun -n 6 ../bin/lfric_atm configuration.nml

        You can change the values of the variables in the PBS header. The
        values used in the suite you are replicating can be found at the top of
        the ``job`` file in the log directory under cylc-run.

        It is likely that this documentation will become out of date. Check the
        ``job`` file, described, above for the list of modules that are loaded
        by the job you are copying.

        Useful PBS commands:

        #. Submit with:   ``qsub <pbs_file>``
        #. Monitor with:  ``qstat -u <userid>``
        #. Kill with:     ``qdel <jobid>``


Tips for debugging
------------------

Once you are able to build and run from the command line, there are a number of
techniques and useful bits of information that can make your life debugging code
a lot easier.

Finding generated code
^^^^^^^^^^^^^^^^^^^^^^

Various parts of the code compiled for LFRic are generated. If you get a
backtrace that refers to line numbers in the generated code, you can find
the code it is referring to under the ``working`` directory within the
top-level directory of your application.

Shortening the run
^^^^^^^^^^^^^^^^^^

If the thing you are testing happens quickly, reduce the number of timesteps
being run by editing the ``configuration.nml``. The start and end timesteps are 
held in the ``&time`` namelist.

.. note::

  Because the entries in the configuration file are often listed in alphabetical
  order, ``timestep_end`` often appears before ``timestep_start``. Make sure
  the numbers are the right way around, or an error will occur.

Set the logging level to "debug"
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Again edit the ``configuration.nml`` and set ``run_log_level='debug'`` in the
``logging`` namelist. Not only does this print out more logging information
that may be useful in debugging, but it also adds a  flush to the logging output
buffer after each message is written. If your  code crashes while the important
log messages you are looking for are in the buffer, they are lost forever.
Flushing the log output buffer after every  message prevents any loss of
messages. It slows the run down - but this is a  debugging run, so it doesn't
really matter.

Set XIOS to write output files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To save disk space many of the suite jobs have the XIOS output turned off, so
you will not get ``xios_client_n.err`` or ``xios_client_n.out`` files. To turn
on these debugging files, edit the ``iodef.xml`` file and find where
``info_level`` and ``print_file`` are defined, then set them to ``100`` and
``.true.`` respectively:

.. code-block:: xml

        <variable id = "info_level" type = "int" >100</variable>
        <variable id = "print_file" type="bool">.true.</variable>

Adding "print" statements to the code
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You can check variables not currently outputted by adding print statements to
the code:

Make sure there is a use statement at the top of the code that includes the
logging modules

.. code-block:: fortran

 use log_mod, only: log_event, &
                    log_scratch_space, &
                    log_level_debug

Then add lines such as:

.. code-block:: fortran

 write(log_scratch_space,"('My variable = ',i0)")my_variable
 call log_event(log_scratch_space, log_level_debug)

Your "print" message will then be written to the PET file for each processor
that runs it.

If you want to check for bit-level changes in variables, it can be useful to
print out exactly how the variable is held in memory using the hexadecimal
format. For example, to write out a 64-bit real number in hexadecimal:

.. code-block:: fortran

 write(log_scratch_space,"('My real variable = ',z16)")my_real_variable
 call log_event(log_scratch_space, log_level_debug)

Editing preprocessed and generated files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The build system processes and generates source code using preprocessors such as
the PSyclone code generator. All modified and unmodified code is put into the
``working`` directory before the final compilation stage. If files are edited
and the build is rerun, then the build will process and recompile modified files
and files that depend on those modified files. Commonly, such a rebuild is much
quicker than the original build as only a few files need to be recompiled.

The original files can be modified in the working copy, or the generated files
can be modified from within the working directory. Pros and cons are as follows:

* If original files are modified, it is easy to scrap the working directory and
  start afresh without losing changes. Then again, to back out changes, it may
  be necessary to check out the code again.
* The working directory holds generated files such as the PSyclone-generated PSy
  layer files. Modifying these files can be very helpful when debugging
  problems. And changes to a generated file can quickly be discarded by
  rebuilding after touching the file that was used to generate it.

Note that if an original ``.x90`` file is modified, both PSyclone-generated
files (a ``.f90`` file and a ``_psy.f90`` file, both with the same prefix as the
``.x90`` file) in the working directory will be overwritten. Bear this in mind
if you have made any extensive changes to one of those files in the working
directory: it is worth occasionally copying your changes to a safe location to
avoid accidentally over-writing them.

The ability to edit PSyclone-generated code has been useful, in particular, when
debugging numerical errors that cause infeasible numbers in fields or that cause
unexpected divergences between two runs, as code can be added to the PSy layer
subroutines to query the contents of fields.

There are some things that do not work:

* Adding or removing files, or adding new use statements often cannot be done
  safely as doing so can upset the dependency analyser.
* The Configurator generates files for reading namelists which have a
  ``_config_mod.f90`` suffix. The above method does not work on these.

