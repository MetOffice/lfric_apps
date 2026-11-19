##############################################################################
# (c) Crown copyright 2026 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
import pytest

from interleaver.config import ExecutableConfig, ConfigError


class TestWorkingExecutableConfig:

    exe = ExecutableConfig(
            name='lfric_atm',
            total_ranks=240,
            ranks_per_node=120,
            rank_depth=1,
            omp_num_threads=1,
            binding_start_core=0
        )

    def test_env_prefix(self):
        assert self.exe.env_prefix == 'LFRIC_ATM'

    def test_cores_per_node(self):
        assert self.exe.cores_per_node == 120

    def test_binding_end_cores(self):
        assert self.exe.binding_end_core == 120

    def test_needed_nodes(self):
        assert self.exe.needed_nodes == 2

    def test_validate(self):
        assert self.exe.validate(node_cores=128) is None

    def test_binding_end_core_error(self):
        self.exe.binding_start_core = 10

        with pytest.raises(ConfigError) as err:
            self.exe.validate(node_cores=128)
            assert "lfric_atm requires cores 10-129, but CORES_PER_NODE is 128" in str(err)

    def test_cores_per_node_error(self):
        self.exe.ranks_per_node = 129
        self.exe.binding_start_core = 0

        with pytest.raises(ConfigError) as err:
            self.exe.validate(node_cores=128)
            assert "lfric_atm requires cores 0-128, but CORES_PER_NODE is 128" in str(err)

    def test_omp_gt_depth_error(self):
        self.exe.omp_num_threads = 8

        with pytest.raises(ConfigError) as err:
            self.exe.validate(node_cores=128)
            assert "LFRIC_ATM_OMP_NUM_THREADS must not exceed LFRIC_ATM_RANK_DEPTH" in str(err)