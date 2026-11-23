"""experiments/t1_ext/geo.py: the two faults of the first real plan run (290da3f), on synthetic input only (no network):

(a) the repository answered one metadata request with HTTP 403 and the same address with 200 a moment later: the fetch
    retries a bounded number of times, with a pause, on 403, 429, 5xx and timeouts, stays on the one host over https and
    follows no redirect;
(b) the SOFT parser refused a whole series at the repository's standard field name "Sample_contact_zip/postal_code": the
    standard slash names are accepted and every other odd key still refuses.
"""

from __future__ import annotations

import io
import sys
import urllib.error
import urllib.request
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "experiments" / "t1_ext"
if str(EXT) not in sys.path:
    sys.path.insert(0, str(EXT))
import geo  # noqa: E402
import t1x_common as cm  # noqa: E402

URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE999nnn/GSE999011/soft/GSE999011_family.soft.gz"


class _Resp(io.BytesIO):
    def __init__(self, body: bytes, url: str = URL):
        super().__init__(body)
        self._url = url

    def geturl(self):
        return self._url


def fake_opener(script):
    """Each call takes the next item of `script`: bytes (a 200 with that body), an int (an HTTP error with that code),
    'timeout', or ('moved', url) (a response whose final address is another)."""
    calls = []

    def open_url(req, timeout):
        calls.append((req.full_url, timeout, dict(req.header_items())))
        item = script[len(calls) - 1]
        if isinstance(item, bytes):
            return _Resp(item)
        if isinstance(item, int):
            raise urllib.error.HTTPError(req.full_url, item, "synthetic", {}, io.BytesIO(b""))
        if item == "timeout":
            raise urllib.error.URLError(TimeoutError("timed out"))
        if item == "read-timeout":
            raise TimeoutError("read timed out")
        return _Resp(b"x", item[1])
    return open_url, calls


@pytest.mark.parametrize("first", [403, 429, 500, 503, "timeout", "read-timeout"])
def test_a_retryable_answer_is_retried_after_a_pause(tmp_path, first):
    open_url, calls = fake_opener([first, b"body"])
    pauses = []
    n, h = geo.get(URL, tmp_path / "a.soft.gz", 100, open_url=open_url, sleep=pauses.append)
    assert (n, (tmp_path / "a.soft.gz").read_bytes()) == (4, b"body") and len(h) == 64
    assert len(calls) == 2 and pauses == [geo.PAUSE_S]
    assert all(c[0] == URL for c in calls)
    assert not (tmp_path / "a.soft.gz.part").exists()


def test_the_retry_is_bounded(tmp_path):
    open_url, calls = fake_opener([403] * 10)
    pauses = []
    with pytest.raises(cm.Refused, match=f"after {geo.TRIES} tries"):
        geo.get(URL, tmp_path / "a", 100, open_url=open_url, sleep=pauses.append)
    assert len(calls) == geo.TRIES and len(pauses) == geo.TRIES - 1
    assert 2 <= geo.TRIES <= 5 and geo.PAUSE_S > 0
    assert not (tmp_path / "a").exists() and not (tmp_path / "a.part").exists()


@pytest.mark.parametrize("code", [400, 401, 404, 410])
def test_other_errors_are_not_retried(tmp_path, code):
    open_url, calls = fake_opener([code, b"body"])
    with pytest.raises(cm.Refused, match="not retried"):
        geo.get(URL, tmp_path / "a", 100, open_url=open_url, sleep=lambda s: None)
    assert len(calls) == 1 and not (tmp_path / "a").exists()


def test_the_byte_limit_is_not_retried(tmp_path):
    open_url, calls = fake_opener([b"x" * 200, b"x"])
    with pytest.raises(cm.Refused, match="byte limit"):
        geo.get(URL, tmp_path / "a", 100, open_url=open_url, sleep=lambda s: None)
    assert len(calls) == 1 and not (tmp_path / "a.part").exists()


def test_a_response_from_another_address_refuses(tmp_path):
    open_url, _ = fake_opener([("moved", "https://example.org/x")])
    with pytest.raises(cm.Refused):
        geo.get(URL, tmp_path / "a", 100, open_url=open_url, sleep=lambda s: None)
    open_url, _ = fake_opener([("moved", URL.replace("GSE999011_family", "GSE999012_family"))])
    with pytest.raises(cm.Refused, match="not the requested address"):
        geo.get(URL, tmp_path / "a", 100, open_url=open_url, sleep=lambda s: None)
    assert not (tmp_path / "a").exists()


def test_the_real_opener_follows_no_redirect():
    handler = [h for h in geo.opener().handlers if isinstance(h, urllib.request.HTTPRedirectHandler)]
    assert len(handler) == 1 and isinstance(handler[0], geo._NoRedirect)
    req = urllib.request.Request(URL)
    for target in ("https://ftp.ncbi.nlm.nih.gov/other", "https://example.org/x", "http://ftp.ncbi.nlm.nih.gov/x"):
        with pytest.raises(cm.Refused, match="redirect"):
            handler[0].redirect_request(req, None, 302, "Found", {}, target)


def test_only_the_one_host_over_https(tmp_path):
    open_url, calls = fake_opener([b"x"])
    for bad in ("https://example.org/x", "http://ftp.ncbi.nlm.nih.gov/x", "https://ftp.ncbi.nlm.nih.gov.evil.org/x"):
        with pytest.raises(cm.Refused):
            geo.get(bad, tmp_path / "a", 100, open_url=open_url, sleep=lambda s: None)
    assert calls == []
    assert geo.HOSTS == ("ftp.ncbi.nlm.nih.gov",)


def test_the_request_names_nobody(tmp_path):
    open_url, calls = fake_opener([b"x"])
    geo.get(URL, tmp_path / "a", 100, open_url=open_url, sleep=lambda s: None)
    ua = {k.lower(): v for k, v in calls[0][2].items()}["user-agent"]
    assert ua == geo.USER_AGENT and "vec" not in ua.lower() and "embryo" not in ua.lower()


# ------------------------------------------------------------------------------------------------ (b) slash field names


def _soft(extra_sample_line):
    return "\n".join(["^SERIES = GSE999011", "!Series_title = synthetic",
                      "!Series_contact_zip/postal_code = 00000",
                      "^SAMPLE = GSM9991101", "!Sample_title = E8.5 heart",
                      "!Sample_contact_zip/postal_code = 00000", extra_sample_line]) + "\n"


def test_the_standard_slash_field_names_are_accepted():
    d = geo.parse_soft(_soft("!Sample_characteristics_ch1 = developmental stage: E8.5"))
    assert d["samples"]["GSM9991101"]["contact_zip/postal_code"] == ["00000"]
    assert d["series"]["contact_zip/postal_code"] == ["00000"]
    # contact fields are never read for a stage
    assert "contact_zip/postal_code" not in cm.stage_fields(d["samples"]["GSM9991101"])
    assert "contact_zip/postal_code" in cm.record_text("SAMPLE", "GSM9991101", d["samples"]["GSM9991101"])


@pytest.mark.parametrize("line", ["!Sample_characteristics/stage = E10.5", "!Sample_contact_zip/postal = 0",
                                  "!Sample_zip/postal_code = 0", "!Sample_contact_zip/postal_code/x = 0",
                                  "!Sample_contact zip = 0", "!Sample_contact_zip/postal_codé = 0",
                                  "!Sample_contact_zip\\postal_code = 0"])
def test_every_other_odd_key_still_refuses(line):
    with pytest.raises(cm.Refused, match="not a plain key"):
        geo.parse_soft(_soft(line))


def test_the_slash_list_is_literal():
    assert geo.SLASH_KEYS == frozenset({"contact_zip/postal_code"})
