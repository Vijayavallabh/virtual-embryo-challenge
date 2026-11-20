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
# The ids the resource's own sample table states for its pinned lines (metadata, column "Term Accession Number";
# the reviewer's point L1): an independent source for the extras of the committed map.
RESOURCE_IDS = {"Adamts3": "ENSMUSG00000043635", "Bap1": "ENSMUSG00000021901", "Brd2": "ENSMUSG00000024335",
                "Camsap3": "ENSMUSG00000044433", "Ccdc6": "ENSMUSG00000048701", "Cnot4": "ENSMUSG00000038784",
                "Col4a3bp": "ENSMUSG00000021669", "Cyp11a1": "ENSMUSG00000032323", "Dennd4c": "ENSMUSG00000038024",
                "Edc4": "ENSMUSG00000036270", "Ehmt1": "ENSMUSG00000036893", "Fcho2": "ENSMUSG00000041685",
                "Fryl": "ENSMUSG00000070733", "H13": "ENSMUSG00000019188", "Morc2a": "ENSMUSG00000034543",
                "Nadk2": "ENSMUSG00000022253", "Otud7b": "ENSMUSG00000038495", "Pdzk1": "ENSMUSG00000038298",
                "Psph": "ENSMUSG00000029446", "Sh3pxd2a": "ENSMUSG00000053617", "Smg9": "ENSMUSG00000002210",
                "Ssr2": "ENSMUSG00000041355", "Tial1": "ENSMUSG00000030846"}
SEX_IDS = {"Xist": "ENSMUSG00000086503", "Ddx3y": "ENSMUSG00000069045", "Eif2s3y": "ENSMUSG00000069049",
           "Kdm5d": "ENSMUSG00000056673", "Uty": "ENSMUSG00000068457"}


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


def test_the_extras_equal_the_ids_the_resource_states_for_its_lines():
    m = cm.load_panel_ids()
    assert {k: m["extra"][k] for k in RESOURCE_IDS} == RESOURCE_IDS
    assert {k: m["extra"][k] for k in SEX_IDS} == SEX_IDS and set(SEX_IDS) == set(cm.SEX_GENES)
    # no extra's id is a panel id
    assert not (set(RESOURCE_IDS.values()) | set(SEX_IDS.values())) & set(m["ids"])
