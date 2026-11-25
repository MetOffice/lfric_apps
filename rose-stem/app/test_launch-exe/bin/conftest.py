import pytest


@pytest.fixture(scope='function')
def xios_server_false_env(monkeypatch):
    monkeypatch.setenv("XIOS_SERVER_MODE", "False")
    monkeypatch.setenv("XIOS_SERVER_RANKS", "0")
    monkeypatch.setenv("xios_nodes", "0")
    monkeypatch.setenv("mpi_parts_xios", "0")