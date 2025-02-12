.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. _perftools_rose-stem:

Perftools in Rose Stem
======================

To make use of these changes, you just need to target the rose stem compiler option that includes perftools.

.. warning::

   Currently it only works on the XC40.

Create a group for Gungho Model or Lfric atm which follows the usual format, but with one exception, it points at:
* ``perftools-compiler``
* eg: ``app_app-conf_mesh_perftools-compiler_opt-bit``
* working eg: ``gungho_model_baroclinic-noio-C24_MG_xc40_perftools-intel_fast-debug-64bit``

The build, run and export will be handled by rose stem.

Two example groups have already been created to do a quick run with:

.. code-block::

   gungho_model_baroclinic-noio-C24_MG_xc40_perftools-intel_fast-debug-64bit
   lfric_atm_nwp_gal9-C48_MG_xc40_perftools-intel_fast-debug-64bit

You can extend any existing groups in the same way, to allow them to run with perftools.

The export task/application, pat_export
---------------------------------------

* Opens the pat_report (so that the binary can be safely thrown away)
* Dumps the text output from the pat_report to a named file
   * Name uses app, date, revision number and pat report id
* Export task can be modified to use:
   * An alternative output-path
   * The current revision number (this has to be manually added however)
   * A choice to move the Pat_report and exports to your datadir (or the alternative output-path)
