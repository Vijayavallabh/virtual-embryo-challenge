import json

import pytest

from vec.data import AuthError, DataFile, load_token


def test_token_file_must_be_private(tmp_path, monkeypatch):
    monkeypatch.delenv("VEC_TOKEN", raising=False)
    path = tmp_path / "token"
    path.write_text("abc\n")
    path.chmod(0o644)
    with pytest.raises(AuthError, match="readable by others"):
        load_token(path)
    path.chmod(0o600)
    assert load_token(path) == "abc"


def test_token_accepts_pasted_session_json(tmp_path, monkeypatch):
    monkeypatch.delenv("VEC_TOKEN", raising=False)
    path = tmp_path / "token"
    path.write_text(json.dumps({"token": "xyz", "account": {}, "expiresAt": 0}))
    path.chmod(0o600)
    assert load_token(path) == "xyz"


@pytest.mark.parametrize("key", ["../etc/passwd", "a/../../b", "x y.h5ad", "a;rm -rf.h5ad"])
def test_unsafe_keys_rejected(key):
    with pytest.raises(ValueError):
        _ = DataFile(key=key, task="T1", role="train", stage="E8.5", modality="scRNA", bytes=1).relpath


def test_tilde_never_reaches_a_quoted_remote_path():
    from vec.data import resolve_dest

    assert resolve_dest("~/v/proj/data/raw", remote=True) == "v/proj/data/raw"
    assert resolve_dest("/abs/path", remote=True) == "/abs/path"
    assert not resolve_dest("~/x", remote=False).startswith("~")


def test_data_module_never_posts():
    import inspect

    import vec.data

    src = inspect.getsource(vec.data)
    assert "method=\"POST\"" not in src and "data=" not in src.replace("metadata", "")


def test_safe_key_kept():
    assert DataFile(key="T1/E8.5_RNA.h5ad", task="T1", role="train", stage="E8.5", modality="scRNA",
                    bytes=1).relpath == "T1/E8.5_RNA.h5ad"
