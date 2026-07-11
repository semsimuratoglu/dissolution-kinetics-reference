"""Canonical labels for the three distinct Polli dissolution comparison sets.

Do not conflate these in outputs or golden checks:

  (a) polli_slow_similar — slow metoprolol 6-unit mean vs slow Unit 6; f2≈91, SIMILAR;
      DDSolver Table-2 golden validation set.
  (b) polli_fast_hc_demo — fast Lopressor-like reference only; HC unconstrained demo.
  (c) polli_fast_vs_slow_dissimilar — fast ref vs slow test; f2≈22, DISSIMILAR;
      dissimilarity demonstration pair.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class PolliSet:
    label: str
    ref_file: str
    test_file: str | None
    publication: str
    profile_character: str
    role: str
    legacy_ref_files: tuple[str, ...] = ()
    legacy_test_files: tuple[str, ...] = ()
    expected_f2_primary: float | None = None
    similarity_verdict: str | None = None


POLLI_SLOW_SIMILAR = PolliSet(
    label="polli_slow_similar",
    ref_file="polli_slow_similar_ref.csv",
    test_file="polli_slow_similar_test.csv",
    legacy_ref_files=(
        "polli_ref.csv",
        "01_polli_metoprolol_WEIBULL.csv",
    ),
    legacy_test_files=(
        "polli_test.csv",
        "02b_polli_TEST.csv",
    ),
    publication="Polli 1997; Zhang 2010 DDSolver supplementary (Table 2)",
    profile_character="Slow metoprolol — 6-unit mean ref vs slow Unit 6 test",
    role="Deterministic-equivalence validation vs golden_polli_slow_similar.json",
    expected_f2_primary=91.02,
    similarity_verdict="SIMILAR",
)

POLLI_FAST_HC_DEMO = PolliSet(
    label="polli_fast_hc_demo",
    ref_file="polli_fast_hc_demo_ref.csv",
    test_file=None,
    legacy_ref_files=("02a_polli_REFERENCE.csv",),
    publication="Polli benchmark / Lopressor-like fast formulation",
    profile_character="Fast Lopressor-like (near-complete by 15 min)",
    role="Hixson-Crowell unconstrained-fit demonstration (kHC≈0.033, R²≈0.83)",
)

POLLI_FAST_VS_SLOW_DISSIMILAR = PolliSet(
    label="polli_fast_vs_slow_dissimilar",
    ref_file="polli_fast_vs_slow_dissimilar_ref.csv",
    test_file="polli_fast_vs_slow_dissimilar_test.csv",
    legacy_ref_files=("02a_polli_REFERENCE.csv",),
    legacy_test_files=("02b_polli_TEST.csv",),
    publication="Polli 1997 benchmark (fast vs slow formulation contrast)",
    profile_character="Fast ref (02a-like) vs slow test (02b-like)",
    role="Dissimilarity demonstration pair (f2≈22)",
    expected_f2_primary=22.49,
    similarity_verdict="DISSIMILAR",
)

ALL_POLLI_SETS: tuple[PolliSet, ...] = (
    POLLI_SLOW_SIMILAR,
    POLLI_FAST_HC_DEMO,
    POLLI_FAST_VS_SLOW_DISSIMILAR,
)

_GOLDEN_PAIR = POLLI_SLOW_SIMILAR


def _matches_ref(ps: PolliSet, ref_file: str) -> bool:
    return ref_file == ps.ref_file or ref_file in ps.legacy_ref_files


def _matches_test(ps: PolliSet, test_file: str) -> bool:
    return test_file == ps.test_file or test_file in ps.legacy_test_files


def resolve_polli_label(ref_file: str, test_file: str | None = None) -> str | None:
    """Return canonical polli set label for a ref[/test] file pair, or None."""
    ref_file = ref_file.strip()
    if test_file is not None:
        test_file = test_file.strip()
        for ps in ALL_POLLI_SETS:
            if ps.test_file is None:
                continue
            if _matches_ref(ps, ref_file) and _matches_test(ps, test_file):
                return ps.label
        return None
    for ps in ALL_POLLI_SETS:
        if ps.test_file is None and _matches_ref(ps, ref_file):
            return ps.label
    return None


def pair_key_label(pair_key: str) -> str:
    """Map a ref:test pair key to an unambiguous polli label when applicable."""
    if ":" not in pair_key:
        return pair_key
    ref_file, test_file = pair_key.split(":", 1)
    label = resolve_polli_label(ref_file, test_file)
    return label if label is not None else pair_key


def is_golden_polli_pair(ref_file: str, test_file: str) -> bool:
    return resolve_polli_label(ref_file, test_file) == POLLI_SLOW_SIMILAR.label


def golden_section_key(golden: dict) -> str | None:
    if "polli_slow_similar" in golden:
        return "polli_slow_similar"
    if "polli_metoprolol" in golden:
        return "polli_metoprolol"
    return None
