"""configs/t2_chamber_classes.yaml (the organ cut's label config, file D; row R21) and its loader in vec.t2_build.

The config ports the committed rule of configs/t1_chamber_classes.yaml (fixed 2026-10-19) to the T2 labels: the
categories, scopes and group bases are T1's, every lineage-map label is listed, each label sits at its group's base
unless an exception names it (rule 2's name exceptions, rule 3's default for null-base groups), and the T2 labels the
committed per-label lists name are all covered. The loader refuses an unmapped label. Committed files only; no data
file is read."""

from __future__ import annotations

import csv
import hashlib
import sys
from pathlib import Path

import numpy as np
import pytest
import yaml

from vec import t1_scope, t2_build
from vec.prep import lineage_map

REPO = Path(__file__).resolve().parents[1]
CONFIG = REPO / "configs" / "t2_chamber_classes.yaml"
CHAMBER = ("cm", "endocardium", "epicardium", "vascular_endothelium", "blood")


def _labels_of(path: Path) -> set[str]:
    with open(path, newline="") as f:
        return {row["label"] for row in csv.DictReader(f)}


def _committed_t2_labels() -> dict[str, set[str]]:
    """The T2 heart labels the committed per-label files name, by stage (E8.25_late has no such file)."""
    sys.path.insert(0, str(Path(__file__).resolve().parent))
    import test_prep

    return {"E10.5": _labels_of(REPO / "experiments/labels/results/E10.5__vs__E9.5.csv") | test_prep.VALIDATION_LABELS[
                "T2/E10.5"],
            "E9.5": _labels_of(REPO / "experiments/harmonise/results/ref_perclass_E9.5.csv"),
            "E8.75": _labels_of(REPO / "experiments/harmonise/results/ref_perclass_E8.75.csv"),
            "E8.5": set(test_prep.VALIDATION_LABELS["T2/E8.5"])}


@pytest.fixture(scope="module")
def spec():
    return t2_build.load_t2_chamber_spec()


def test_the_kept_categories_are_t1s_chamber_scope(spec):
    t1 = t1_scope.load_spec()
    assert spec.kept == CHAMBER == tuple(t1.in_scope("chamber"))
    assert spec.categories == t1.categories
    assert spec.group_base == t1.group_base
    assert t2_build.ORGAN_KEPT_SCOPE == "chamber"


def test_every_lineage_map_label_is_listed_and_nothing_else(spec):
    assert set(spec.labels) == set(lineage_map()["labels"])


def test_each_label_follows_the_rule(spec):
    """A label sits at its group's base unless `exceptions` names it; every label of a null-base group is named there;
    no exception is stale."""
    lm = lineage_map()["labels"]
    for lab, cat in spec.labels.items():
        base = spec.group_base[lm[lab]]
        if base is None:
            assert lab in spec.exceptions, lab
        elif cat != base:
            assert lab in spec.exceptions, (lab, cat, base)
    for lab in spec.exceptions:
        assert lab in spec.labels
        assert spec.group_base[lm[lab]] != spec.labels[lab] or spec.group_base[lm[lab]] is None, lab


def test_the_rule_2_and_rule_3_placements(spec):
    """Rule 2's name exceptions and rule 3's default, as T1's committed file places the same labels."""
    t1 = t1_scope.load_spec()
    assert spec.labels["ST"] == "surrounding" == t1.labels["ST"]
    assert spec.labels["Early VIC"] == "endocardium" == t1.labels["Early VIC"]
    for lab in ("LPM", "Unknown", "Artificial"):
        assert spec.labels[lab] == "surrounding"
    for lab in ("Cush-EndoMT-Endo", "Valve-EndoMT-Endo"):
        assert spec.labels[lab] == "endocardium"
    for lab in ("Epi", "Proepi", "Epicardium", "Proepicardium"):
        assert spec.labels[lab] == "epicardium"
    for lab in ("Blood", "Blood Progenitor", "Macrophage"):
        assert spec.labels[lab] == "blood"
    shared = set(spec.labels) & set(t1.labels)
    assert {lab for lab in shared if spec.labels[lab] != t1.labels[lab]} == set()  # T1's own labels agree


def test_the_core_lies_inside_the_kept_categories(spec):
    """Every core label of the lineage map (CM plus endocardium) is kept: the organ cut's core is always kept."""
    lm = lineage_map()
    core = {lab for lab, g in lm["labels"].items() if lm["groups"][g]["core"]}
    assert core and all(spec.labels[lab] in spec.kept for lab in core)


def test_every_committed_t2_heart_label_is_mapped(spec):
    by_stage = _committed_t2_labels()
    assert all(by_stage.values())
    for stage, labels in by_stage.items():
        assert labels <= set(spec.labels), (stage, sorted(labels - set(spec.labels)))
    # what the organ cut keeps of E10.5's labels, by the rule (no share of any stage)
    kept = {lab for lab in by_stage["E10.5"] if spec.labels[lab] in spec.kept}
    assert kept == {"V-CM", "A-CM", "SV-CM", "OFT-CM", "Chamber-Endo", "Cush-EndoMT-Endo", "Valve-EndoMT-Endo",
                    "Epi", "BEC", "early-VEC", "Great Artery Endoth", "HEM-Endoth"}


def test_the_spec_records_the_files_sha256(spec):
    assert spec.sha256 == hashlib.sha256(CONFIG.read_bytes()).hexdigest()
    assert Path(spec.path) == CONFIG


def test_the_loader_refuses_an_unmapped_label(spec):
    with pytest.raises(t2_build.UnmappedT2LabelError) as e:
        t2_build.t2_chamber_categories(["V-CM", "Mystery", "BEC", "Mystery", "Other"], spec)
    assert e.value.labels == {"Mystery": 2, "Other": 1}
    assert "Mystery" in str(e.value) and "refused" in str(e.value)


def test_the_loader_maps_labels_in_order(spec):
    got = t2_build.t2_chamber_categories(np.array(["Dorsal-PAM", "V-CM", "ST", "Epi", "V-CM"]), spec)
    assert list(got) == ["surrounding", "cm", "surrounding", "epicardium", "cm"]


def test_the_loader_refuses_undefined_categories(tmp_path):
    raw = yaml.safe_load(CONFIG.read_text())
    raw["labels"]["V-CM"] = "heart_muscle"
    bad = tmp_path / "bad.yaml"
    bad.write_text(yaml.safe_dump(raw))
    with pytest.raises(ValueError, match="undefined category"):
        t2_build.load_t2_chamber_spec(bad)
    raw = yaml.safe_load(CONFIG.read_text())
    for c in raw["categories"].values():
        c["scope"] = "adjacent"
    bad.write_text(yaml.safe_dump(raw))
    with pytest.raises(ValueError, match="no category has scope"):
        t2_build.load_t2_chamber_spec(bad)
