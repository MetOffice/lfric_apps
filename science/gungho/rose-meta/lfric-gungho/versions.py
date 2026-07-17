import sys

from metomi.rose.upgrade import MacroUpgrade  # noqa: F401

from .version31_32 import *


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

class vn32_t551(MacroUpgrade):
    """Upgrade macro for ticket #551 by Christine Johnson."""

    BEFORE_TAG = "vn3.2_t386"
    AFTER_TAG = "vn3.2_t551"

    def upgrade(self, config, meta_config=None):
        # Commands From: rose-meta/lfric-gungho
        partitioner = self.get_setting_value(
            config, ["namelist:partitioning", "partitioner"]
        )
        multigrid_chain_nitems = self.get_setting_value(
            config, ["namelist:multigrid", "multigrid_chain_nitems"]
        )
        if partitioner == "'planar'":
            self.change_setting_value(
                config,
                ["namelist:base_mesh", "prime_mesh_name"],
                "'l0_planar'",
            )
            if multigrid_chain_nitems == "4":
                self.add_setting(
                    config,
                    ["namelist:multigrid", "chain_mesh_tags"],
                    "'l0_planar','l1_planar','l2_planar','l3_planar'",
                    forced=True,
                )
            elif multigrid_chain_nitems == "3":
                self.add_setting(
                    config,
                    ["namelist:multigrid", "chain_mesh_tags"],
                    "'l0_planar','l1_planar','l2_planar'",
                    forced=True,
                )
            elif multigrid_chain_nitems == "2":
                self.add_setting(
                    config,
                    ["namelist:multigrid", "chain_mesh_tags"],
                    "'l0_planar','l1_planar'",
                    forced=True,
                )
        else:
            self.change_setting_value(
                config,
                ["namelist:base_mesh", "prime_mesh_name"],
                "'l0_cubedsphere'",
            )
            if multigrid_chain_nitems == "4":
                self.add_setting(
                    config,
                    ["namelist:multigrid", "chain_mesh_tags"],
                    "'l0_cubedsphere','l1_cubedsphere','l2_cubedsphere','l3_cubedsphere'",
                    forced=True,
                )
            elif multigrid_chain_nitems == "3":
                self.add_setting(
                    config,
                    ["namelist:multigrid", "chain_mesh_tags"],
                    "'l0_cubedsphere','l1_cubedsphere','l2_cubedsphere'",
                    forced=True,
                )
            elif multigrid_chain_nitems == "2":
                self.add_setting(
                    config,
                    ["namelist:multigrid", "chain_mesh_tags"],
                    "'l0_cubedsphere','l1_cubedsphere'",
                    forced=True,
                )
        partitioner_dest = self.get_setting_value(
            config, ["namelist:partitioning(destination)", "partitioner"]
        )
        destination_mesh_name = self.get_setting_value(
            config, ["namelist:lfric2lfric", "destination_mesh_name"]
        )
        if partitioner_dest == "'planar'":
            if destination_mesh_name == "'dynamics'":
                self.add_setting(
                    config,
                    ["namelist:lfric2lfric", "destination_mesh_name"],
                    "'l0_planar",
                    forced=True,
                )
            elif destination_mesh_name == "'multigrid_l1'":
                self.add_setting(
                    config,
                    ["namelist:lfric2lfric", "destination_mesh_name"],
                    "'l1_planar",
                    forced=True,
                )
        elif partitioner_dest == "'cubedsphere'":
            if destination_mesh_name == "'dynamics'":
                self.add_setting(
                    config,
                    ["namelist:lfric2lfric", "destination_mesh_name"],
                    "'l0_planar",
                    forced=True,
                )
            elif destination_mesh_name == "'multigrid_l1'":
                self.add_setting(
                    config,
                    ["namelist:lfric2lfric", "destination_mesh_name"],
                    "'l1_cubedsphere",
                    forced=True,
                )
        partitioner_source = self.get_setting_value(
            config, ["namelist:partitioning(source)", "partitioner"]
        )
        source_mesh_name = self.get_setting_value(
            config, ["namelist:lfric2lfric", "source_mesh_name"]
        )
        if partitioner_source == "'planar'":
            if source_mesh_name == "'dynamics'":
                self.add_setting(
                    config,
                    ["namelist:lfric2lfric", "source_mesh_name"],
                    "'l0_planar",
                    forced=True,
                )
            elif source_mesh_name == "'multigrid_l1'":
                self.add_setting(
                    config,
                    ["namelist:lfric2lfric", "source_mesh_name"],
                    "'l1_planar",
                    forced=True,
                )
        elif partitioner_source == "'cubedsphere'":
            if source_mesh_name == "'dynamics'":
                self.add_setting(
                    config,
                    ["namelist:lfric2lfric", "source_mesh_name"],
                    "'l0_cubedsphere",
                    forced=True,
                )
            elif source_mesh_name == "'multigrid_l1'":
                self.add_setting(
                    config,
                    ["namelist:lfric2lfric", "source_mesh_name"],
                    "'l1_cubedsphere",
                    forced=True,
                )
        if source_mesh_name == "'C12'":
                self.add_setting(
                    config,
                    ["namelist:lfric2lfric", "source_mesh_name"],
                    "'l1_cubedsphere",
                    forced=True,
                )

        aerosol_mesh_name = self.get_setting_value(
            config, ["namelist:multires_coupling", "aerosol_mesh_name"])        
        if aerosol_mesh_name == "'multigrid_l2'":
            self.add_setting(
                config,
                ["namelist:multires_coupling", "aerosol_mesh_name"],
                "'l2_cubedsphere",
                forced=True,
            )
        multires_coupling_mesh_tags = self.get_setting_value(
            config, ["namelist:multires_coupling", "multires_coupling_mesh_tags"])        
        if multires_coupling_mesh_tags == "'dynamics','multigrid_l2'":
            self.add_setting(
                config,
                ["namelist:multires_coupling", "multires_coupling_mesh_tags"],
                "'l0_cubedsphere','l2_cubedsphere'",
                forced=True,
            )
        orography_mesh_name = self.get_setting_value(
            config, ["namelist:multires_coupling", "orography_mesh_name"])        
        if orography_mesh_name == "'dynamics'":
            self.add_setting(
                config,
                ["namelist:multires_coupling", "orography_mesh_name"],
                "'l0_cubedsphere",
                forced=True,
            )
            
        physics_mesh_name = self.get_setting_value(
            config, ["namelist:multires_coupling", "physics_mesh_name"])
                        
        if physics_mesh_name == "'dynamics'":
            self.add_setting(
                config,
                ["namelist:multires_coupling", "physics_mesh_name"],
                "'l0_cubedsphere",
                forced=True,
            )

        blpert_mesh_name = self.get_setting_value(
            config, ["namelist:stochastic_physics", "blpert_mesh_name"])

        if blpert_mesh_name == "'multigrid_l3'":
            self.add_setting(
                config,
                ["namelist:stochastic_physics", "blpert_mesh_name"],
                "'l3_planar",
                forced=True,
            )

        return config, self.reports
