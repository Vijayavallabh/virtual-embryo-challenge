"""Protected stage windows (docs/DATA_POLICY.md) and stage-label parsing, shared by the external-data filters.

Windows are in embryonic days (E-days). A row of external data may be written only if its stage lies outside the
window of every task the data will serve (a model shared across tasks respects the union of their windows):

  T1, T2-heart (extrapolation)   excluded after E9.5 up to and including E13.5      (9.5, 13.5]
  T2-heart (interpolation)       excluded strictly between E8.25 and E8.75          (8.25, 8.75)
  T2-embryo (interpolation)      excluded strictly between E7.25 and E8.0           (7.25, 8.0)

T3 is protected by genotype (no Ctnnb1 or Gata4 perturbation data, other alleles or phenocopies), not by a stage window,
so asking for T3 here raises instead of silently passing. Postnatal labels map to E-days as P<n> -> 19.5 + n (birth
taken as E19.5); every postnatal day is far outside every window, so the exact offset never changes a decision.
Unknown or malformed labels raise (fail closed).
"""

from __future__ import annotations

import re
from collections.abc import Iterable

# (low, high, low_inclusive, high_inclusive): a stage s is inside when low <(=) s <(=) high
WINDOWS: dict[str, tuple[tuple[float, float, bool, bool], ...]] = {
    "T1": ((9.5, 13.5, False, True),),
    "T2-heart": ((9.5, 13.5, False, True), (8.25, 8.75, False, False)),
    "T2-embryo": ((7.25, 8.0, False, False),),
}
TASKS = tuple(WINDOWS)
BIRTH_EDAY = 19.5
_STAGE = re.compile(r"^\s*([EePp])\s*(\d+(?:\.\d+)?)\s*$")


def parse_stage(label: object) -> float:
    """E-day of a stage label: 'E9.5' -> 9.5, 'E14.333' -> 14.333, 'E18' -> 18.0, 'P0' -> 19.5, 'P9' -> 28.5."""
    m = _STAGE.match(str(label))
    if not m:
        raise ValueError(f"unparseable stage label {label!r}")
    value = float(m.group(2))
    return value if m.group(1) in "Ee" else BIRTH_EDAY + value


def check_tasks(tasks: Iterable[str]) -> tuple[str, ...]:
    tasks = tuple(tasks)
    if not tasks:
        raise ValueError("no task given: name every task the data will serve")
    for t in tasks:
        if t == "T3":
            raise NotImplementedError("T3 is protected by genotype, not by a stage window: screen genotypes separately")
        if t not in WINDOWS:
            raise ValueError(f"unknown task {t!r}; expected one of {TASKS}")
    return tasks


def in_window(eday: float, task: str) -> bool:
    """True if a stage (in E-days) lies inside one of the task's protected windows."""
    for low, high, low_inc, high_inc in WINDOWS[check_tasks([task])[0]]:
        above = eday >= low if low_inc else eday > low
        below = eday <= high if high_inc else eday < high
        if above and below:
            return True
    return False


def allowed(eday: float, tasks: Iterable[str]) -> bool:
    """True if a stage lies outside the windows of every task in `tasks` (their union)."""
    return not any(in_window(eday, t) for t in check_tasks(tasks))


def allowed_label(label: object, tasks: Iterable[str]) -> bool:
    return allowed(parse_stage(label), tasks)
