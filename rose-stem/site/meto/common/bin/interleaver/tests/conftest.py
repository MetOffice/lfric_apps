import pytest
import os
import tempfile


THIS_DIR = os.path.dirname(os.path.abspath(__file__))


def get_nodelist(interleaver_config):
    '''
    Helper function for generating the node list
    '''
    nodes = []
    seen = set()
    for line in interleaver_config.pbs_nodefile.read_text().splitlines():
        node = line.strip()
        if node and node not in seen:
            nodes.append(node)
            seen.add(node)

    return nodes


@pytest.fixture(scope='function')
def interleaver_env(monkeypatch):
    monkeypatch.setenv('EXECUTABLES', 'lfric_atm xios')
    monkeypatch.setenv('CORES_PER_NODE', '192')
    monkeypatch.setenv('PLACEMENT_PATTERN', 'standard')
    monkeypatch.setenv('SEPARATE_EXECUTABLES', '')
    monkeypatch.setenv('PBS_NODEFILE', os.path.join(THIS_DIR, 'test_files', 'PBS_NODEFILE'))
    monkeypatch.setenv('CYLC_TASK_WORK_DIR', str(tempfile.mkdtemp()))
    # exec specifics - LFRIC_ATM
    monkeypatch.setenv('LFRIC_ATM_TOTAL_RANKS', '1176')
    monkeypatch.setenv('LFRIC_ATM_RANKS_PER_NODE', '168')
    monkeypatch.setenv('LFRIC_ATM_RANK_DEPTH', '1')
    monkeypatch.setenv('LFRIC_ATM_OMP_NUM_THREADS', '1')
    monkeypatch.setenv('LFRIC_ATM_BINDING_START_CORE', '0')
    # exec specifics - XIOS
    monkeypatch.setenv('XIOS_TOTAL_RANKS', '42')
    monkeypatch.setenv('XIOS_RANKS_PER_NODE', '6')
    monkeypatch.setenv('XIOS_RANK_DEPTH', '4')
    monkeypatch.setenv('XIOS_OMP_NUM_THREADS', '1')
    monkeypatch.setenv('XIOS_BINDING_START_CORE', '0')


@pytest.fixture(scope='function')
def striped_interleaver_env(monkeypatch):
    monkeypatch.setenv('EXECUTABLES', 'lfric_atm xios nemo')
    monkeypatch.setenv('CORES_PER_NODE', '192')
    monkeypatch.setenv('PLACEMENT_PATTERN', 'standard')
    monkeypatch.setenv('SEPARATE_EXECUTABLES', 'nemo')
    monkeypatch.setenv('PBS_NODEFILE', os.path.join(THIS_DIR, 'test_files', 'ASYM_PBS_NODEFILE'))
    monkeypatch.setenv('CYLC_TASK_WORK_DIR', str(tempfile.mkdtemp()))
    # exec specifics - LFRIC_ATM
    monkeypatch.setenv('LFRIC_ATM_TOTAL_RANKS', '1176')
    monkeypatch.setenv('LFRIC_ATM_RANKS_PER_NODE', '168')
    monkeypatch.setenv('LFRIC_ATM_RANK_DEPTH', '1')
    monkeypatch.setenv('LFRIC_ATM_OMP_NUM_THREADS', '1')
    monkeypatch.setenv('LFRIC_ATM_BINDING_START_CORE', '0')
    # exec specifics - XIOS
    monkeypatch.setenv('XIOS_TOTAL_RANKS', '42')
    monkeypatch.setenv('XIOS_RANKS_PER_NODE', '6')
    monkeypatch.setenv('XIOS_RANK_DEPTH', '4')
    monkeypatch.setenv('XIOS_OMP_NUM_THREADS', '1')
    monkeypatch.setenv('XIOS_BINDING_START_CORE', '0')
    # exec specifics - NEMO
    monkeypatch.setenv('NEMO_TOTAL_RANKS', '192')
    monkeypatch.setenv('NEMO_RANKS_PER_NODE', '96')
    monkeypatch.setenv('NEMO_RANK_DEPTH', '2')
    monkeypatch.setenv('NEMO_OMP_NUM_THREADS', '2')
    monkeypatch.setenv('NEMO_BINDING_START_CORE', '0')

