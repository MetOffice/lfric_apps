# LFRic Apps Fab Build Scripts

Make sure you have Fab version 2.0.1 or later installed (in addition to all
LFRic apps requirements of course).

## Setting up Site- and Platform-specific Settings
The lfric_core repository must contain a setup for the site and platform,
see the documentation in ``lfric_core/lfric_build/README.md``.

## Building the Gravity Wave

In order to build Gravity Wave, change into the directory
```$LFRIC_APPS/applications/gravity_wave```,
and use the following command:

```
./fab_gravity_wave.py --core $LFRIC_CORE --nprocs 4 --site nci --platform gadi \
                      --suite gnu
```
Select an appropriate number of processes to run in parallel, and your site and platform.
If you don't have a default compiler suite in your site-specific setup (or
want to use a non-default suite), use the ``--suite`` option. Once the process is finished,
you should have a binary in the directory
```./fab-workspace/gravity_wave-full-COMPILER``` (where ```COMPILER``` is the compiler
used, e.g. ```mpif90-gfortran```).

Using ```./fab_gravity_wave.py -h``` will show a help message with all supported command line
options (and their default value). If a default value is listed using an environment
variables (```(default: $SITE or 'default')```), the corresponding environment variable
is used if no command line option has been specified.

A different compilation profile can be specified using ```--profile``` option. Note
that the available compilation profiles can vary from site to site (see
[Fab documentation](https://metoffice.github.io/fab/fab_base/config.html) for details).

If Fab has issues finding a compiler, you can use the Fab debug option
```--available-compilers```, which will list all compilers and linkers Fab has
identified as being available.

## Testing with pFUnit

Any application script can additionally inherit from the ``pfunit_mixin.py``
mixin, e.g.:
```
from lfric_apps_base import LFRicAppsBase  # noqa: E402
from pfunit_mixin import PfUnitMixin  # noqa: E402


class FabGravityWave(PfUnitMixin, LFRicAppsBase):

```
This will by default also build any tests in a ``unit-test`` directory.
Additionally, a new command line option ``--no-test`` is added if
building of the tests should be disabled.
