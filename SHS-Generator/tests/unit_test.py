"""In-process unit tests for shs_generator.

Unlike integration_test.py (which drives the CLI via subprocess), these import
the module directly and call functions/methods in-process. They are fast (no
per-test interpreter startup) and can reach internals the CLI cannot exercise
-- most notably mutate_multiplet, which needs a base triple that a dot-bracket
structure can never produce.

These lock the CURRENT intended behaviour so a refactor can't silently change
it. Run just these with:  pytest tests/unit_test.py
"""

import sys
from pathlib import Path

import pytest

# Make the package dir (parent of tests/) importable, then import the module.
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import shs_generator as shs  # noqa: E402


# A pair (as a 2-char string) is structurally valid if it is Watson-Crick or a
# G-U wobble. mutate_pair / _partner_for must only ever emit these.
VALID_PAIRS = {"AU", "UA", "GC", "CG", "GU", "UG"}


def make_args(**overrides):
    """Build an argparse.Namespace with the real shs_generator defaults, then
    apply overrides. Uses the actual parser so the defaults stay in sync with
    the source rather than being hand-copied here."""
    saved_argv = sys.argv
    sys.argv = ["shs_generator.py"]
    try:
        args = shs.parse_args()
    finally:
        sys.argv = saved_argv
    for key, value in overrides.items():
        setattr(args, key, value)
    return args


@pytest.fixture
def make_gen(tmp_path):
    """Factory for MsaGenerator instances. Defaults to a deterministic seed and
    a tmp output dir (MsaGenerator.__init__ mkdir's it)."""
    def _make(**overrides):
        overrides.setdefault("seed", 0)
        overrides.setdefault("output_json_dir", str(tmp_path / "out"))
        return shs.MsaGenerator(make_args(**overrides))
    return _make


############################### pair_indices ###############################

def test_pair_indices_nested():
    g = shs.MsaGenerator(make_args(seed=0, output_json_dir="/tmp"))
    assert g.pair_indices("(())") == {0: 3, 3: 0, 1: 2, 2: 1}


def test_pair_indices_pseudoknot_multiple_families():
    # () and [] families are tracked on independent stacks, so crossing
    # (pseudoknotted) pairs are all recovered.
    g = shs.MsaGenerator(make_args(seed=0, output_json_dir="/tmp"))
    pairs = g.pair_indices("((..[[..))..]]")
    assert pairs == {0: 9, 9: 0, 1: 8, 8: 1, 4: 13, 13: 4, 5: 12, 12: 5}


def test_pair_indices_unbalanced_closer_is_skipped():
    # An extra ')' with an empty stack is skipped (logged, not fatal); the
    # balanced pairs are still returned.
    g = shs.MsaGenerator(make_args(seed=0, output_json_dir="/tmp"))
    assert g.pair_indices("(()))") == {0: 3, 3: 0, 1: 2, 2: 1}


############################### pair normalisation ###############################

def test_normalize_pairs_from_list_dedupes_and_drops_self_pairs():
    # list of [i, j, ...] rows -> sorted unique (i<j) tuples, self-pairs dropped.
    assert shs._normalize_pairs([[5, 9, 0], [1, 3, 0], [3, 3], [9, 5]]) == [(1, 3), (5, 9)]


def test_normalize_pairs_from_dict():
    assert shs._normalize_pairs({0: 3, 3: 0, 1: 2}) == [(0, 3), (1, 2)]


def test_convert_pred_pairs_list_to_bidirectional_dict():
    assert shs.convert_pred_pairs([[1, 4], [2, 3]]) == {1: 4, 4: 1, 2: 3, 3: 2}


def test_convert_pred_pairs_dict_is_rebuilt_and_int_keyed():
    # A dict input is now rebuilt via the shared parser (no longer passed
    # through), so a well-formed bidirectional dict is preserved by value and
    # non-int keys/values are int-converted.
    assert shs.convert_pred_pairs({1: 4, 4: 1}) == {1: 4, 4: 1}
    assert shs.convert_pred_pairs({"2": "3"}) == {2: 3, 3: 2}


############################### mutate_pair / _partner_for ###############################

def test_mutate_pair_only_emits_valid_pairs(make_gen):
    g = make_gen(wobble_prob=0.1)
    for a in "AUGC":
        for b in "AUGC":
            for _ in range(50):
                assert g.mutate_pair(a, b) in VALID_PAIRS


def test_mutate_pair_wobble_always_forces_gu(make_gen):
    g = make_gen(wobble_prob=1.0)
    for _ in range(50):
        assert g.mutate_pair("C", "G") in {"GU", "UG"}


def test_partner_for_watson_crick_without_wobble(make_gen):
    g = make_gen(wobble_prob=0.0)
    assert (g._partner_for("A"), g._partner_for("U"),
            g._partner_for("G"), g._partner_for("C")) == ("U", "A", "C", "G")


def test_partner_for_wobble_shifts_gu(make_gen):
    g = make_gen(wobble_prob=1.0)
    assert g._partner_for("G") == "U"
    assert g._partner_for("U") == "G"


############################### mutate_multiplet (CLI-unreachable) ###############################

def test_mutate_multiplet_chain_keeps_every_edge_valid(make_gen):
    # Base triple as a chain: position 1 is the hub shared by (0,1) and (1,2).
    g = make_gen(wobble_prob=0.1)
    positions, edges = (0, 1, 2), [(0, 1), (1, 2)]
    for _ in range(200):
        assigned = g.mutate_multiplet("GCGC", positions, edges)
        assert set(assigned) == set(positions)
        for a, b in edges:
            assert assigned[a] + assigned[b] in VALID_PAIRS


def test_mutate_multiplet_hub_keeps_every_edge_valid(make_gen):
    # Hub: position 1 pairs with three partners (0, 2, 3).
    g = make_gen(wobble_prob=0.1)
    positions, edges = (0, 1, 2, 3), [(0, 1), (1, 2), (1, 3)]
    for _ in range(200):
        assigned = g.mutate_multiplet("GCGC", positions, edges)
        assert set(assigned) == set(positions)
        for a, b in edges:
            assert assigned[a] + assigned[b] in VALID_PAIRS


def test_mutate_multiplet_covers_all_positions_even_with_cycle(make_gen):
    # A closing edge of a cycle is left as-is (not re-mutated), but every
    # position must still receive an assignment.
    g = make_gen(wobble_prob=0.1)
    positions, edges = (0, 1, 2), [(0, 1), (1, 2), (0, 2)]
    assigned = g.mutate_multiplet("GCG", positions, edges)
    assert set(assigned) == set(positions)


############################### derive_multiplets ###############################

def test_derive_multiplets_disabled_returns_empty(make_gen):
    g = make_gen(triplet_prob=0.0)
    assert g.derive_multiplets([(0, 1), (1, 2)], 3) == []


def test_derive_multiplets_promotes_connected_component(make_gen):
    # Shared node 1 -> connected component {0,1,2} of size 3 -> one multiplet.
    g = make_gen(triplet_prob=1.0)
    assert g.derive_multiplets([(0, 1), (1, 2)], 3) == [
        {"positions": (0, 1, 2), "edges": [(0, 1), (1, 2)]}
    ]


def test_derive_multiplets_ignores_ordinary_pairs(make_gen):
    # Two disjoint size-2 pairs are ordinary pairs, not multiplets.
    g = make_gen(triplet_prob=1.0)
    assert g.derive_multiplets([(0, 5), (1, 4)], 6) == []


def test_derive_multiplets_is_deterministic_under_seed(make_gen):
    pair_list = [(0, 1), (1, 2), (5, 6), (6, 7)]
    a = make_gen(seed=123, triplet_prob=0.5).derive_multiplets(pair_list, 8)
    b = make_gen(seed=123, triplet_prob=0.5).derive_multiplets(pair_list, 8)
    assert a == b


############################### generate_msa structural invariants ###############################

def test_generate_msa_row_count_and_query_first(make_gen):
    g = make_gen(N=5)
    msa = g.generate_msa("GGCC", {0: 3, 3: 0, 1: 2, 2: 1})
    assert len(msa) == 5           # query row + N-1 generated rows
    assert msa[0] == "GGCC"        # first row is always the unmutated query


def test_generate_msa_preserves_fully_paired_stem_when_keep_is_one(make_gen):
    # stem_keep_prob=1.0 on a fully paired structure -> no position ever mutates.
    g = make_gen(N=6, stem_keep_prob=1.0)
    seq = "GGCC"
    pairs = g.pair_indices("(())")
    msa = g.generate_msa(seq, pairs)
    assert msa == [seq] * 6


def test_generate_msa_deletes_all_loop_positions_when_forced(make_gen):
    # All-unpaired structure with deletion_prob_loop=1.0 and every other
    # insertion/deletion knob at 0 -> each non-query row is fully deleted.
    g = make_gen(
        N=4,
        insertion_prob_loop=0.0,
        long_insertion_prob=0.0,
        long_deletion_prob=0.0,
        deletion_prob_loop=1.0,
    )
    seq = "ACGU"
    pairs = g.pair_indices("....")  # {} -> all positions unpaired
    msa = g.generate_msa(seq, pairs)
    assert msa[0] == seq
    assert all(row == "----" for row in msa[1:])


############################### build_output_name ###############################

def test_build_output_name_always_includes_stem_keep(make_gen):
    # stem_keep_prob must appear regardless of branch or triplet setting.
    g = make_gen(pdb_id="X", N=3, seed=7, stem_keep_prob=0.99,
                 triplet_prob=0.0, structure_predictor=None)
    g.max_insertion_length = 0
    g.max_deletion_length = 0
    name = g.build_output_name()
    assert "_stemkeep_0.99_" in name
    assert "_triplet_" not in name  # triplet disabled -> no triplet suffix


def test_build_output_name_appends_triplet_suffix_when_enabled(make_gen):
    g = make_gen(pdb_id="X", triplet_prob=0.5, triplet_keep_prob=0.99)
    g.max_insertion_length = 1
    g.max_deletion_length = 1
    name = g.build_output_name()
    assert "_stemkeep_" in name
    assert "_triplet_0.5_keep_0.99" in name


############################### plot_final_features ###############################

def test_plot_final_features_writes_valid_png(make_gen, tmp_path, monkeypatch):
    # plot_final_features saves to ./msa_final_plots/msa_final_<pdb_id>.png,
    # relative to the CWD -> chdir into tmp_path so it lands there.
    import msa_plotting
    monkeypatch.chdir(tmp_path)
    g = make_gen(N=5)
    seq = "GGGAACCCGGGAACCC"
    pairs = g.pair_indices("(((..[[[)))..]]]")
    msa = g.generate_msa(seq, pairs)
    msa_plotting.plot_final_features(msa, seq, pairs, "PLOTTEST")
    png = tmp_path / "msa_final_plots" / "msa_final_PLOTTEST.png"
    assert png.exists(), f"Expected plot at {png}"
    assert png.stat().st_size > 0
    assert png.read_bytes()[:8] == b"\x89PNG\r\n\x1a\n", "Output is not a valid PNG"


############################### msa_plotting helpers ###############################
# Imported locally so the rest of the unit suite stays matplotlib-free.

def test_pad_lowercase_drops_lowercase_and_pads():
    import msa_plotting
    assert msa_plotting.pad_lowercase("AbCdE", 8) == "ACE-----"


def test_bool_left_deletion_marks_gaps():
    import msa_plotting
    assert list(msa_plotting.bool_left_deletion("A-G-")) == [False, True, False, True]
