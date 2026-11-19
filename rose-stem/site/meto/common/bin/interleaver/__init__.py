"""XIOS colocation hostfile and binding generation.

This module provides placement strategies and hostfile/binding generation
for colocating XIOS servers with LFRIC atmosphere ranks.

Main entry points:
  - ColocationCLI.run() - Generate hostfile and binding files
  - BindingExecutor.run() - Execute command with per-rank CPU binding
"""

__version__ = "1.0.0"
