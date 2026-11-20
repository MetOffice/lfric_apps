##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
import os
import pytest

from conftest import get_nodelist
from interleaver.placements import StandardPlacement
from interleaver.config import InterleaverConfig


class TestStandardPlacement:

    def test_node_blocks(self, interleaver_env):
        '''
        Testing that the correct number of nodes are returned for each
        executable. Our test configuration uses 7 nodes. A standard placement
        strategy should split each executables ranks over 7 nodes.
        '''
        cfg = InterleaverConfig.from_env()
        nodelist = get_nodelist(cfg)

        placement = StandardPlacement(cfg, nodelist)

        node_blocks = placement.get_node_blocks()

        for _, item in node_blocks.items():
            assert len(item) == 7


    def test_nonuniform_blocks(self, interleaver_env):
        '''
        Testing that we get an expected number of nodes in the node block when
        reducing the number of TOTAL_RANKS by enough to effectively remove that
        executable from one of the nodes.
        In this case, when reducing the number of XIOS ranks to 36, we would
        populate 6 of the 7 nodes with XIOS executables and have no executables
        on the final node.
        We still expect the node to be included in the nodelist, however, since
        the total number of required nodes is the combination of all executables
        requirements.
        '''
        # using ppn=6 this would require 6 nodes
        os.environ['XIOS_TOTAL_RANKS'] = '36'

        cfg = InterleaverConfig.from_env()
        nodelist = get_nodelist(cfg)

        placement = StandardPlacement(cfg, nodelist)

        node_blocks = placement.get_node_blocks()

        assert len(node_blocks['lfric_atm']) == 7
        assert len(node_blocks['xios']) == 7

    def test_nonuniform_blocks_validation(self, interleaver_env):
        '''
        Testing that we get an expected number of nodes in the node block when
        reducing the number of TOTAL_RANKS by enough to effectively remove that
        executable from one of the nodes.
        In this case, when reducing the number of XIOS ranks to 36, we would
        populate 6 of the 7 nodes with XIOS executables and have no executables
        on the final node.
        We still expect the node to be included in the nodelist, however, since
        the total number of required nodes is the combination of all executables
        requirements.
        '''
        os.environ['XIOS_TOTAL_RANKS'] = '36'

        cfg = InterleaverConfig.from_env()
        nodelist = get_nodelist(cfg)

        placement = StandardPlacement(cfg, nodelist)

        # Same deal, as above - having excess nodes is fine
        placement.validate()

    def test_placement_hostfile_entries(self, interleaver_env):
        '''
        This tests that we get an expected hostfile entry for each executable.

        This is possibly the mst important step since it tells the launcher
        where to place each executable.

        '''
        cfg = InterleaverConfig.from_env()
        nodelist = get_nodelist(cfg)

        placement = StandardPlacement(cfg, nodelist)

        hostfile = placement.get_hostfile_entries()

                         # 6 ranks per node
        expected_xios = ['nid001000', 'nid001000', 'nid001000', 'nid001000',
                         'nid001000', 'nid001000',
                         'nid001001', 'nid001001', 'nid001001', 'nid001001',
                         'nid001001', 'nid001001',
                         'nid001002', 'nid001002', 'nid001002', 'nid001002',
                         'nid001002', 'nid001002',
                         'nid001003', 'nid001003', 'nid001003', 'nid001003',
                         'nid001003', 'nid001003',
                         'nid001004', 'nid001004', 'nid001004', 'nid001004',
                         'nid001004', 'nid001004',
                         'nid001005', 'nid001005', 'nid001005', 'nid001005',
                         'nid001005', 'nid001005',
                         'nid001006', 'nid001006', 'nid001006', 'nid001006',
                         'nid001006', 'nid001006']

        assert hostfile['xios'] == expected_xios

    def test_placement_total_hostfile(self, interleaver_env):
        '''
        This tests that we get an expected hostfile containing the final order
        of our node ids.

        This concatenates the hostfile entries in MPMD order.
        '''
        cfg = InterleaverConfig.from_env()
        nodelist = get_nodelist(cfg)

        placement = StandardPlacement(cfg, nodelist)

        total_host = placement.get_total_hostfile()

        assert 'xios' not in total_host

    def test_placement_binding_entries(self, interleaver_env):
        '''
        This tests that we get an expected binding list for each executable.

        This is where the cpu identifiers are assigned to each rank - the depth
        parameter.

        '''
        cfg = InterleaverConfig.from_env()
        nodelist = get_nodelist(cfg)

        placement = StandardPlacement(cfg, nodelist)

        bindings = placement.get_binding_entries()

                         # Using depth 4 for each rank, 6 ranks per node
        expected_xios = ['168-171', '172-175', '176-179', '180-183', '184-187', '188-191',
                         '168-171', '172-175', '176-179', '180-183', '184-187', '188-191',
                         '168-171', '172-175', '176-179', '180-183', '184-187', '188-191',
                         '168-171', '172-175', '176-179', '180-183', '184-187', '188-191',
                         '168-171', '172-175', '176-179', '180-183', '184-187', '188-191',
                         '168-171', '172-175', '176-179', '180-183', '184-187', '188-191',
                         '168-171', '172-175', '176-179', '180-183', '184-187', '188-191']

        assert bindings['xios'] == expected_xios

    def test_placement_total_binding(self, interleaver_env):

        '''
        This tests that we get an expected binding list for each executable.

        This is where the cpu identifiers are assigned to each rank - the depth
        parameter.

        '''
        cfg = InterleaverConfig.from_env()
        nodelist = get_nodelist(cfg)

        placement = StandardPlacement(cfg, nodelist)

        total_binding = placement.get_total_binding()

        assert 'xios' not in total_binding