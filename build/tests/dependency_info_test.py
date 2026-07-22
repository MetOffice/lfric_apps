##############################################################################
# (c) Crown copyright Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################

"""
This module tests dependency_info.
"""

from pathlib import Path

import pytest

from build.dependency_info import DependencyInfo, RepoInfo


def test_repo_info_append_and_iterate() -> None:
    """
    Check RepoInfo stores appended source/ref pairs and iterates over them.
    """
    repo_info = RepoInfo()

    repo_info.append("git@example.com:first.git", "first-ref")
    repo_info.append("git@example.com:second.git", "second-ref")

    assert repo_info.path == Path(".")
    assert list(repo_info) == [
        RepoInfo.SourceRef("git@example.com:first.git", "first-ref"),
        RepoInfo.SourceRef("git@example.com:second.git", "second-ref"),
    ]


@pytest.mark.parametrize("filename", [None, ""])
def test_dependency_info_empty_filename(filename) -> None:
    """
    Check that no filename creates an empty dependency set.
    """
    dependency_info = DependencyInfo(filename, [])

    assert dependency_info == {}
    assert dependency_info.get_repo_names() == []
    with pytest.raises(KeyError):
        dependency_info.get_repo_info("missing")


def test_dependency_info_reads_single_and_multiple_dependencies(
        tmp_path: Path) -> None:
    """
    Check that both single dependency definitions and lists are supported.
    """
    dependency_file = tmp_path / "dependencies.yaml"
    dependency_file.write_text(
        "lfric_core:\n"
        "    source:\n"
        "    ref:\n"
        "SimSys_Scripts:\n"
        "    - source: git@github.com:MetOffice/SimSys_Scripts.git\n"
        "      ref: cab3315147a3c7e8546dda559d3da0fccd702f29\n"
        "    - source: git@github.com:MetOffice/SimSys_Scripts-fork.git\n"
        "      ref: feature-branch\n",
        encoding="utf8"
    )

    dependency_info = DependencyInfo(dependency_file, [])

    assert dependency_info.get_repo_names() == ["lfric_core", "SimSys_Scripts"]
    assert list(dependency_info.get_repo_info("lfric_core")) == [
        RepoInfo.SourceRef(None, None)
    ]
    assert list(dependency_info.get_repo_info("SimSys_Scripts")) == [
        RepoInfo.SourceRef(
            "git@github.com:MetOffice/SimSys_Scripts.git",
            "cab3315147a3c7e8546dda559d3da0fccd702f29",
        ),
        RepoInfo.SourceRef(
            "git@github.com:MetOffice/SimSys_Scripts-fork.git",
            "feature-branch",
        ),
    ]


def test_dependency_info_only(tmp_path: Path) -> None:
    """
    Check that specifying 'only' works as expected
    """
    dependency_file = tmp_path / "dependencies.yaml"
    dependency_file.write_text(
        "repo1:\n"
        "    source: git@bgithub.com/repo1\n"
        "    ref: 1\n"
        "repo2:\n"
        "    source: git@bgithub.com/repo2\n"
        "    ref: 2\n"
        "repo3:\n"
        "    source: git@bgithub.com/repo3\n"
        "    ref: 3\n",
        encoding="utf8"
    )

    dependency_info = DependencyInfo(dependency_file, [])
    assert dependency_info.get_repo_names() == ["repo1", "repo2", "repo3"]

    dependency_info = DependencyInfo(dependency_file, ["repo1"])
    assert dependency_info.get_repo_names() == ["repo1"]

    dependency_info = DependencyInfo(dependency_file, ["repo1", "repo2"])
    assert dependency_info.get_repo_names() == ["repo1", "repo2"]

    dependency_info = DependencyInfo(dependency_file, ["repo1", "repo2",
                                                       "repo3"])
    assert dependency_info.get_repo_names() == ["repo1", "repo2", "repo3"]


@pytest.mark.parametrize(
    "yaml_text, expected_message",
    [
        (
            "test_repo:\n"
            "    ref: test-ref\n",
            "does not contain a 'source' definition for repo 'test_repo'",
        ),
        (
            "test_repo:\n"
            "    source: git@example.com:test.git\n",
            "does not contain a 'ref' definition for repo 'test_repo'",
        ),
    ],
)
def test_dependency_info_rejects_missing_required_keys(
        tmp_path: Path,
        yaml_text: str,
        expected_message: str) -> None:
    """
    Check that malformed dependency entries raise a helpful RuntimeError.
    """
    dependency_file = tmp_path / "dependencies.yaml"
    dependency_file.write_text(yaml_text, encoding="utf8")

    with pytest.raises(RuntimeError, match=expected_message):
        DependencyInfo(dependency_file, [])
