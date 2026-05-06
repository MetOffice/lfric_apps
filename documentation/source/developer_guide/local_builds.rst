.. -----------------------------------------------------------------------------
    (c) Crown copyright Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------

LFRic Apps Command Line Builds
==============================

The LFRic makefile build system is reliant on having code from both lfric_core
and lfric_apps available when launching, preventing the makefile from exporting
lfric_core as part of the build process. This introduces a circular dependency
in the makefile build system, not present when everything was on a single repo.
As the makefile system is planned to be replaced by fab, it was decided to not
spend too much effort rewriting the makefiles and instead use a wrapper script
to export lfric_core and call the make command.

The wrapper script is located at ``/build/local_build.py``. It is intended that
all command line options available in the makefile system can be used here - if
any have been missed then as it's a python file it should be easily extendable.
The script can be launched from anywhere by calling it, eg:

``./local_build.py <PROJECT>``

where ``PROJECT`` is the name of the LFRic Apps project being built, eg.
``gungho_model``. This can be an application, science or interface section,
although some of these may only contain a subset of possible make targets. The
``PROJECT`` is always required. The script has been setup to behave as the old
make system did. It will create a ``working`` and a ``bin`` dir inside the
project directory by default. Its first step is to export the lfric_core repo
and then rsync it to the working dir so that incremental builds continue to
work.

At the Met Office you will need to load the lfric environments before starting
the script.

This table lists the command line arguments available:

+----------------------+-----------------------------------------------------------+
| *Option*             | *Default*                   | *Description*               |
+======================+=============================+=============================+
| ``-c --core_source`` | If not provided, the script | The lfric_core source       |
|                      | will parse the              | being used. Can be a        |
|                      | ``dependencies.yaml`` file  | remote github repository    |
|                      | and use the core source set | or a local clone.           |
|                      | there.                      |                             |
+----------------------+-----------------------------+-----------------------------+
| ``-w --working_dir`` | The project directory.      | The working directory       |
|                      |                             | for the build process.      |
+----------------------+-----------------------------+-----------------------------+
| ``-j --ncores``      | 4                           | Integer, the number of      |
|                      |                             | cores for the build task    |
+----------------------+-----------------------------+-----------------------------+
| ``-t --target``      | build                       | The makefile target, eg.    |
|                      |                             | build, unit-tests, clean    |
+----------------------+-----------------------------+-----------------------------+
| ``-o --optlevel``    | None, Uses the default set  | The optimisation level      |
|                      | in the makefile, usually    | of the build process.       |
|                      | ``fast-debug``.             |                             |
+----------------------+-----------------------------+-----------------------------+
| ``-p --psyclone``    | None, Uses the default set  | Value passed to the         |
|                      | in the makefile.            | ``PSYCLONE_TRANSFORMATION`` |
|                      |                             | variable in the makefile.   |
+----------------------+-----------------------------+-----------------------------+
| ``-v --verbose``     | None, or "off"              | Supplying this argument,    |
|                      |                             | which takes no argument,    |
|                      |                             | will request verbose output |
|                      |                             | from the makefile.          |
+----------------------+-----------------------------+-----------------------------+
| ``-m --mirrors``     | False                       | If True, this will attempt  |
| ``store_true``       |                             | to extract using local      |
|                      |                             | github mirrors              |
+----------------------+-----------------------------+-----------------------------+
| ``--mirror-loc``     | MetOffice Mirror Location   | The path to the github      |
|                      |                             | mirror location             |
+----------------------+-----------------------------+-----------------------------+

Incremental Builds
------------------

The local build script will attempt to build incrementally if a previous attempt
at the build exists. This should happen automatically if the working directory
is the same. If there are large changes then it may be sensible to start the
build afresh by cleaning the build ``-t clean`` (or deleting the working
directory).
