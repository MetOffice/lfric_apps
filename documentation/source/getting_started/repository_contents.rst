.. -----------------------------------------------------------------------------
    (c) Crown copyright 2025 Met Office. All rights reserved.
    The file LICENCE, distributed with this code, contains details of the terms
    under which the code may be used.
   -----------------------------------------------------------------------------
.. Some of the content of this file has been produced with the assistance of
   GitHub Copilot.

.. _repository_contents:

Contents of the Repository
--------------------------

The repository contains the following directories:

- The ``applications`` directory contains the different applications that have
  been developed using the LFRic Core and the science libraries. Applications
  include the main |Momentum| Atmosphere application as well as smaller
  applications that test subsections of different models (e.g. transport,
  solver, etc.).
- The ``build`` directory contains additional build scripts required by the
  modelling system.
- The ``documentation`` directory contains the source code for this
  documentation.
- The ``interfaces`` directory contains interface code to couple to sub models
  such as the JULES land surface model.
- The ``rose-meta`` dorectory contains links to the rose metadata for the
  different libraries and applications.
- The ``rose-stem`` directory contains the system test suites for development
  testing.
- The ``science`` directory contains several libraries of science code which can
  be used in the development of applications. For example, the **gungho**
  library provides the dynamical core used in the |Momentum| Atmosphere
  application.

Many of the directories contain directories of unit and integration
tests, directories of Rose metadata, and Makefiles for building
applications.
