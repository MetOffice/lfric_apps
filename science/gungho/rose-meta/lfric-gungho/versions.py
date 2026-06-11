import sys

from metomi.rose.upgrade import MacroUpgrade  # noqa: F401

from .version30_31 import *


class UpgradeError(Exception):
    """Exception created when an upgrade fails."""

    def __init__(self, msg):
        self.msg = msg

    def __repr__(self):
        sys.tracebacklimit = 0
        return self.msg

    __str__ = __repr__


"""
Copy this template and complete to add your macro

class vnXX_txxx(MacroUpgrade):
    # Upgrade macro for <TICKET> by <Author>

    BEFORE_TAG = "vnX.X"
    AFTER_TAG = "vnX.X_txxx"

    def upgrade(self, config, meta_config=None):
        # Add settings
        return config, self.reports
"""

class vn31_t11(MacroUpgrade):
    # Upgrade macro for <TICKET> by <Author>

    BEFORE_TAG = "vn3.1"
    AFTER_TAG = "vn3.1_t11"

    def upgrade(self, config, meta_config=None):
        # Add settings
        partitioner = self.get_setting_value(config, ["namelist:partitioning","partitioner"])
        multigrid_chain_nitems = self.get_setting_value(config, ["namelist:multigrid","multigrid_chain_nitems"])
        if partitioner == "'planar'":
            self.add_setting(
                config,
                ["namelist:base_mesh","prime_mesh_name"],
                "'planar_l0'",
                forced=True
                )
            if multigrid_chain_nitems == "4":  
                self.add_setting(
                    config,
                    ["namelist:multigrid","chain_mesh_tags"],
                    "'planar_l0','planar_l1','planar_l2','planar_l3'",
                    forced=True
                    )
            elif multigrid_chain_nitems == "3":  
                self.add_setting(
                    config,
                    ["namelist:multigrid","chain_mesh_tags"],
                    "'planar_l0','planar_l1','planar_l2'",
                    forced=True
                    )
            elif multigrid_chain_nitems == "2":  
                self.add_setting(
                    config,
                    ["namelist:multigrid","chain_mesh_tags"],
                    "'planar_l0','planar_l1'",
                    forced=True
                    )
        else:
            self.add_setting(
                config,
                ["namelist:base_mesh","prime_mesh_name"],
                "'cubedsphere_l0'",
                forced=True
                )
            if multigrid_chain_nitems == "4":  
                self.add_setting(
                    config,
                    ["namelist:multigrid","chain_mesh_tags"],
                    "'cubedsphere_l0','cubedsphere_l1','cubedsphere_l2','cubedsphere_l3'",
                    forced=True
                    )
            elif multigrid_chain_nitems == "3":  
                self.add_setting(
                    config,
                    ["namelist:multigrid","chain_mesh_tags"],
                    "'cubedsphere_l0','cubedsphere_l1','cubedsphere_l2'",
                    forced=True
                    )
            elif multigrid_chain_nitems == "2":  
                self.add_setting(
                    config,
                    ["namelist:multigrid","chain_mesh_tags"],
                    "'cubedsphere_l0','cubedsphere_l1'",
                    forced=True
                    )
        partitioner_dest = self.get_setting_value(config, ["namelist:partitioning(destination)","partitioner"])
        destination_mesh_name = self.get_setting_value(config, ["namelist:lfric2lfric","destination_mesh_name"])
        if partitioner_dest == "'planar'":
            if destination_mesh_name == "'dynamics'":
                self.add_setting(
                    config,
                    ["namelist:lfric2lfric","destination_mesh_name"],
                    "'planar_l0",
                    forced=True
                    )
            elif destination_mesh_name == "'multigrid_l1'":
                self.add_setting(
                    config,
                    ["namelist:lfric2lfric","destination_mesh_name"],
                    "'planar_l1",
                    forced=True
                    )
        elif partitioner_dest == "'cubedsphere'":
            if destination_mesh_name == "'dynamics'":
                self.add_setting(
                    config,
                    ["namelist:lfric2lfric","destination_mesh_name"],
                    "'cubedsphere_l0",
                    forced=True
                    )
            elif destination_mesh_name == "'multigrid_l1'":
                self.add_setting(
                    config,
                    ["namelist:lfric2lfric","destination_mesh_name"],
                    "'cubedsphere_l1",
                    forced=True
                    )
        partitioner_source = self.get_setting_value(config, ["namelist:partitioning(source)","partitioner"])
        source_mesh_name = self.get_setting_value(config, ["namelist:lfric2lfric","source_mesh_name"])
        if partitioner_source == "'planar'":
            if source_mesh_name == "'dynamics'":
                self.add_setting(
                    config,
                    ["namelist:lfric2lfric","source_mesh_name"],
                    "'planar_l0",
                    forced=True
                    )
            elif source_mesh_name == "'multigrid_l1'":
                self.add_setting(
                    config,
                    ["namelist:lfric2lfric","source_mesh_name"],
                    "'planar_l1",
                    forced=True
                    )
        elif partitioner_source == "'cubedsphere'":
            if source_mesh_name == "'dynamics'":
                self.add_setting(
                    config,
                    ["namelist:lfric2lfric","source_mesh_name"],
                    "'cubedsphere_l0",
                    forced=True
                    )
            elif source_mesh_name == "'multigrid_l1'":
                self.add_setting(
                    config,
                    ["namelist:lfric2lfric","source_mesh_name"],
                    "'cubedsphere_l1",
                    forced=True
                    )
        return config, self.reports
