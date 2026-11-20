##############################################################################
# (c) Crown copyright 2025 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
import pytest
import os

from interleaver.config import InterleaverConfig, ConfigError, ExecutableConfig


class TestInterleaverConfig:

    def test_from_env(self, interleaver_env):
        '''
        '''
        expected_execs = {
        'lfric_atm': ExecutableConfig(
            name='lfric_atm',
            total_ranks=1176,
            ranks_per_node=168,
            rank_depth=1,
            omp_num_threads=1,
            binding_start_core=0),
        'xios': ExecutableConfig(
            name='xios',
            total_ranks=42,
            ranks_per_node=6,
            rank_depth=4,
            omp_num_threads=1,
            binding_start_core=168)
        }

        cfg = InterleaverConfig.from_env()

        assert cfg.placement_pattern == 'standard'
        assert len(cfg.separate_executables) == 0
        assert cfg.executables == expected_execs

    def test_get_int_return_default(self, interleaver_env):
        '''
        '''
        # was originally '1176'
        original = InterleaverConfig.get_int('LFRIC_ATM_TOTAL_RANKS',
                                        environ=os.environ, default=None)

        del os.environ['LFRIC_ATM_TOTAL_RANKS']

        modified = InterleaverConfig.get_int('LFRIC_ATM_TOTAL_RANKS',
                                        environ=os.environ, default='failure')

        assert original is not modified
        assert modified == 'failure'
        assert isinstance(original, int)

    def test_get_int_error(self, interleaver_env):
        '''
        '''
        os.environ['LFRIC_ATM_TOTAL_RANKS'] = 'not_an_int'

        with pytest.raises(ConfigError) as err:
            InterleaverConfig.get_int('LFRIC_ATM_TOTAL_RANKS',
                                      environ=os.environ)

            assert ('LFRIC_ATM_TOTAL_RANKS must be supplied as a positive '
                    'integer by the workflow' in str(err))