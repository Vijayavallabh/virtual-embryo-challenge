"""Stage 8 of experiments/t3_ext (campaign/t3-ext-4): the team's held list in the code (t3x_common.HELD_EXTRA), joined
into the derived held set, and GSE205117 held (the first reviewer of ec86a1b, lesser point 2). Synthetic YAML through
the function argument only: nothing here touches the network or real data."""

from __future__ import annotations

import hashlib
import importlib.util
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
EXT = ROOT / "experiments" / "t3_ext"
if str(EXT) not in sys.path:
    sys.path.insert(0, str(EXT))
import t3x_common as cm  # noqa: E402


def _load(name):
    spec = importlib.util.spec_from_file_location(f"t3_ext_held_{name}", EXT / f"{name}.py")
    m = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = m
    spec.loader.exec_module(m)
    return m


fe = _load("fetch")

# the sha256 of the team's held-accession file (its accessions sorted, one per line); fixed 2026-11-21
HELD_LIST_SHA256 = "5df13a82927231c28227232e061452651099050ad02d6a6b6ca1adcc21852048"
HELD = sorted(cm.HELD_EXTRA)

HEAD = """version: 1
qc_floors:
  min_confident_share_per_sample: 0.5
datasets:
"""


def _entry(acc, screen, extra=""):
    return f"""  - accession: {acc}
    status: cleared
    tier: unrelated
    study: {acc}
    design: germline
    control_relation: same_series
    repository: geo
    screen_series: [{", ".join(screen)}]
    stages_used: [E8.5]
    sample_file_patterns: ['matrix\\.mtx(\\.gz)?$']
{extra}    pinned: []
"""


def _sister(acc, screen):
    return (f"    control_series:\n      accession: {acc}\n      screen_series: [{', '.join(screen)}]\n"
            '      relation: {same_lab: "a", same_protocol: "b", same_publication: "c"}\n      reader_note: "d"\n')


def _write(tmp_path, body) -> Path:
    p = tmp_path / "datasets.yaml"
    p.write_text(HEAD + body)
    return p


def test_held_extra_is_exactly_the_team_list():
    assert len(cm.HELD_EXTRA) == 82
    assert hashlib.sha256(("\n".join(HELD) + "\n").encode()).hexdigest() == HELD_LIST_SHA256
    assert cm.NEVER_FETCH <= cm.HELD_EXTRA
    assert "GSE205117" in cm.HELD_EXTRA


def test_the_held_set_holds_every_listed_accession(tmp_path):
    y = _write(tmp_path, _entry("GSE930001", ["GSE930001"]))
    allds = cm.load_datasets(y)
    assert cm.HELD_EXTRA <= cm.held_series(allds)
    assert cm.HELD_EXTRA <= fe.held_set(allds["GSE930001"])
    assert cm.HELD_EXTRA <= cm.held_series(cm.load_datasets())


@pytest.mark.parametrize("acc", HELD)
def test_every_listed_accession_is_refused_as_a_dataset_a_screen_series_and_a_control_series(tmp_path, acc):
    # as the dataset's own series: dataset() refuses it whatever the YAML says, and load_datasets refuses it cleared
    y = _write(tmp_path, _entry(acc, [acc]))
    with pytest.raises(cm.Refused, match="held or excluded series|NEVER_FETCH"):
        cm.load_datasets(y)
    with pytest.raises(cm.Refused, match="HELD_EXTRA|NEVER_FETCH"):
        cm.dataset(acc, y)
    # as a further screened series of a cleared dataset
    y = _write(tmp_path, _entry("GSE930001", ["GSE930001", acc]))
    with pytest.raises(cm.Refused, match="held or excluded series"):
        cm.load_datasets(y)
    # as a control series (its accession, or a series of its screen_series)
    for sis, screen in ((acc, [acc]), ("GSE930002", ["GSE930002", acc])):
        y = _write(tmp_path, _entry("GSE930001", ["GSE930001"], _sister(sis, screen)))
        with pytest.raises(cm.Refused):
            cm.load_datasets(y)
    # a sample whose record lists it is not pinnable
    ds = cm.load_datasets(_write(tmp_path, _entry("GSE930001", ["GSE930001"])))["GSE930001"]
    why = fe.series_membership(ds, {"series_id": ["GSE930001", acc]})
    assert why and "held or excluded" in why


def test_a_control_series_on_the_held_list_names_the_held_set_in_its_refusal(tmp_path):
    y = _write(tmp_path, _entry("GSE930001", ["GSE930001"], _sister("GSE205117", ["GSE205117"])))
    with pytest.raises(cm.Refused, match=r"held or excluded series \['GSE205117'\]"):
        cm.load_datasets(y)


def test_gse205117_is_held_in_the_committed_file():
    ds = cm.load_datasets()
    assert ds["GSE205117"]["status"] == "held_phenocopy"
    assert ds["GSE205117"]["status_reason"] == ("the two readers did not clear its review term (decision log "
                                                "2026-11-20)")
    with pytest.raises(cm.Refused, match="HELD_EXTRA"):
        cm.dataset("GSE205117")
    assert "GSE205117" not in (EXT / "run_dgx.sh").read_text().split("ACCS=")[1].splitlines()[0]
    # no cleared dataset screens a listed accession
    for a, d in ds.items():
        if d["status"] == "cleared":
            assert not set(cm.all_screened(d)) & cm.HELD_EXTRA, a
