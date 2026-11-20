"""The committed gene-id map of the bulk Ensembl route (experiments/t3_ext/panel_ids.tsv) and its pin."""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
EXT = REPO / "experiments" / "t3_ext"
if str(EXT) not in sys.path:
    sys.path.insert(0, str(EXT))
import t3x_common as cm  # noqa: E402

ENS = re.compile(r"ENSMUSG\d{11}")


def test_the_committed_map_hashes_to_its_pin_and_loads():
    assert cm.PANEL_IDS_SHA256 == "5fc8e6676b29faadf430d2569886490bd32a6d44a3ab2b8e679ec8babdd36a07"
    assert cm.sha256_file(cm.PANEL_IDS) == cm.PANEL_IDS_SHA256
    m = cm.load_panel_ids()
    assert len(m["symbols"]) == 500 and len(set(m["symbols"])) == 500
    assert m["header"]["source_file"] == "GSE277839/GSM8532456/GSM8532456_U_E9_5_WT_E_features.tsv.gz"
    assert m["header"]["reference_genes"] == "500"
    # every panel gene has one id, no id serves two genes, and no chromosome is stated (a 10x features file has none)
    assert all(i is not None and ENS.fullmatch(i) for i in m["ids"]) and len(set(m["ids"])) == 500
    assert all(c is None for c in m["chromosome"])


def test_the_source_sample_is_a_pinned_sample_of_a_cleared_single_cell_dataset():
    ds = cm.dataset("GSE277839")
    assert ds["status"] == "cleared" and not cm.is_bulk(ds)
    pin = [p for p in ds["pinned"] if p["accession"] == "GSM8532456"]
    assert len(pin) == 1 and pin[0]["arm"] == "control"
    assert "GSM8532456_U_E9_5_WT_E_features.tsv.gz" in pin[0]["files"]


def test_the_extras_hold_the_sex_genes_and_every_pinned_line_of_the_many_line_set():
    m = cm.load_panel_ids()
    for g in cm.SEX_GENES:
        assert ENS.fullmatch(m["extra"][g])
    lines = sorted({p["gene"] for p in cm.dataset("E-MTAB-7322")["pinned"]})
    assert len(lines) == 25 and all(ln in m["extra"] for ln in lines)
    # two lines carry no id in the source (a microRNA and a symbol the source's annotation no longer uses): the
    # genotype guard then refuses their units, which is the fail-closed outcome
    assert sorted(ln for ln in lines if m["extra"][ln] is None) == ["Mir96", "Zfyve20"]
    ids = [m["extra"][ln] for ln in lines if m["extra"][ln] is not None]
    assert all(ENS.fullmatch(i) for i in ids) and len(set(ids)) == 23


def test_a_changed_map_or_an_unset_pin_is_refused(tmp_path, monkeypatch):
    bad = tmp_path / "panel_ids.tsv"
    bad.write_bytes(cm.PANEL_IDS.read_bytes() + b"\n")
    monkeypatch.setattr(cm, "PANEL_IDS", bad)
    with pytest.raises(cm.Refused, match="not the pinned PANEL_IDS_SHA256"):
        cm.load_panel_ids()
    monkeypatch.setattr(cm, "PANEL_IDS_SHA256", None)
    with pytest.raises(cm.Refused, match="PANEL_IDS_SHA256 is None"):
        cm.load_panel_ids()
