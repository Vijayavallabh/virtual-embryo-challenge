"""The test helpers' git never acts on a repository the environment points at (docs/PLAN.md decision log
2026-11-19 (t): a suite run with GIT_DIR set to the real checkout put 22 empty commits on its main)."""

from conftest import GIT_LOCATION_VARS, git, git_env, make_repo


def test_a_test_git_ignores_a_repository_named_by_the_environment(tmp_path, monkeypatch):
    other = tmp_path / "other"
    make_repo(other)
    before = git(other, "rev-list", "--count", "HEAD")
    mine = tmp_path / "mine"
    make_repo(mine)
    monkeypatch.setenv("GIT_DIR", str(other / ".git"))
    monkeypatch.setenv("GIT_WORK_TREE", str(other))
    assert not set(GIT_LOCATION_VARS) & set(git_env())
    git(mine, "commit", "-q", "--allow-empty", "-m", "probe")
    assert git(mine, "log", "-1", "--format=%s") == "probe"
    monkeypatch.delenv("GIT_DIR")
    monkeypatch.delenv("GIT_WORK_TREE")
    assert git(other, "rev-list", "--count", "HEAD") == before  # nothing landed in the other repository
