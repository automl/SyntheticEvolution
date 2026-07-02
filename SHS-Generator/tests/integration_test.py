
import json
import subprocess
import sys
from pathlib import Path


################################## Helper functions for integration tests ##################################

def sanity_check(json_file, N, dialect="alphafold3", filename=None):
    data = json.loads(json_file.read_text(encoding="utf-8"))
    # check name if given
    if filename is not None:
        assert json_file.name == filename, f"Expected output filename {filename}, got {json_file.name}"
    assert data["dialect"] == dialect, f"Expected dialect {dialect}, got {data['dialect']}"

    for chain in data.get("sequences", []):
        if not "rna" in chain:
            continue
        assert "unpairedMsa" in chain["rna"], "No unpairedMsa found in RNA chain"
        msa_lines = chain["rna"]["unpairedMsa"].split(">")[1:]
        assert len(msa_lines) == N, (
            f"Expected {N} generated unpaired MSA sequences, got {len(msa_lines)}"
        )

def equivalent_json_files(new_file, expected_file):
    data1 = json.loads(new_file.read_text(encoding="utf-8"))
    data2 = json.loads(expected_file.read_text(encoding="utf-8"))
    assert data1 == data2, (f"Generated JSON did not match expected fixture {expected_file}")


def get_output_json(output_dir):
    json_files = list(output_dir.glob("*.json"))
    assert len(json_files) == 1, f"Expected exactly one JSON output, found {len(json_files)}"
    return json_files[0]


def run_generator(output_dir, *cli_args):
    """Invoke shs_generator.py with the given CLI args plus --output_json_dir.
    Returns the completed process. main() logs-and-swallows exceptions, so the
    process exits 0 even on failure; callers inspect stderr / the output dir."""
    script = Path(__file__).resolve().parents[1] / "shs_generator.py"
    return subprocess.run(
        [sys.executable, str(script), *cli_args, "--output_json_dir", str(output_dir)],
        check=True,
        capture_output=True,
        text=True,
    )


def fixture(name):
    """Path to a fixture file living next to this test module."""
    return Path(__file__).with_name(name)


def check_output_matches(output_dir, N, expected_name, expected_file):
    """Assert the single output JSON has the expected name, N MSA rows, and is
    byte-for-byte equivalent to the recorded fixture."""
    output_json = get_output_json(output_dir)
    sanity_check(output_json, N, filename=expected_name)
    equivalent_json_files(output_json, expected_file)

################################### Integration tests ##################################

def test_basic(tmp_path):
    output_dir = tmp_path / "out"
    run_generator(
        output_dir,
        "--rna-seq", "ACGU", "--protein-seq", "MK", "--structure", "(())",
        "--pdb_id", "TEST", "-N", "2", "--seed", "1",
    )
    output_json = get_output_json(output_dir)
    sanity_check(output_json, 2)
    assert output_json.exists(), f"Expected output JSON to exist at {output_json}"


def test_basic_matches_current_state(tmp_path):
    output_dir = tmp_path / "out"
    run_generator(
        output_dir,
        "--rna-seq", "ACGU", "--protein-seq", "MK", "--structure", "(())",
        "--pdb_id", "TEST", "-N", "2", "--seed", "1",
    )
    expected_name = "TEST_custom_rnamsa_N2_seed1_insl_0.2_dell_0.8_inss_0.01_dels_0.01_lins_0.05_ldels_0.05_maxinslen_0_maxdellen_0_wobble_0.1_None.json"
    check_output_matches(output_dir, 2, expected_name,
                         fixture("test_basic_matches_current_state_out.json"))


def test_advanced_matches_current_state(tmp_path):
    output_dir = tmp_path / "out"
    run_generator(
        output_dir,
        "--rna-seq", "CCUGGUAUUGCAGUACCUCCAGGU",
        "--protein-seq", "VKLTAELIEQAARRKVKLKERQEAEKMFK",
        "--pdb_id", "TEST", "-N", "10", "--seed", "42",
        "--structure_predictor", "rnafold",
    )
    expected_name = "TEST_custom_rnamsa_N10_seed42_insl_0.2_dell_0.8_inss_0.01_dels_0.01_lins_0.05_ldels_0.05_maxinslen_2_maxdellen_2_wobble_0.1_rnafold.json"
    check_output_matches(output_dir, 10, expected_name,
                         fixture("test_advanced_matches_current_state_out.json"))


def test_json_input_matches_current_state(tmp_path):
    output_dir = tmp_path / "out"
    run_generator(
        output_dir,
        "--input_json_path", str(fixture("test_json_input_matches_current_state_in.json")),
        "--pdb_id", "TEST", "-N", "10", "--seed", "42",
        "--structure_predictor", "rnafold",
    )
    expected_name = "TEST_custom_rnamsa_N10_seed42_insl_0.2_dell_0.8_inss_0.01_dels_0.01_lins_0.05_ldels_0.05_maxinslen_2_maxdellen_2_wobble_0.1_stemkeep_0.99_rnafold.json"
    check_output_matches(output_dir, 10, expected_name,
                         fixture("test_json_input_matches_current_state_out.json"))


def test_triplet_pseudoknot_matches_current_state(tmp_path):
    # Combined coverage for two paths the other tests miss:
    #   * pseudoknots: a multi-family dot-bracket structure "(((..[[[)))..]]]"
    #     exercises pair_indices' handling of the () and [] bracket families
    #     (6 base pairs across two crossing stems).
    #   * triplet mode: --triplet-prob > 0 enables the triplet co-mutation
    #     plumbing (derive_multiplets + the "_triplet_..._keep_..." name suffix).
    #     A dot-bracket structure can never produce a real multiplet (every
    #     position pairs with exactly one partner), so this locks the
    #     triplet-enabled code path and naming, not mutate_multiplet itself;
    #     exercising that needs a predictor emitting overlapping pairs.
    output_dir = tmp_path / "out"
    run_generator(
        output_dir,
        "--rna-seq", "GGGAACCCGGGAACCC", "--protein-seq", "MK",
        "--structure", "(((..[[[)))..]]]",
        "--pdb_id", "TEST", "-N", "5", "--seed", "1", "--triplet-prob", "0.5",
    )
    expected_name = "TEST_custom_rnamsa_N5_seed1_insl_0.2_dell_0.8_inss_0.01_dels_0.01_lins_0.05_ldels_0.05_maxinslen_1_maxdellen_1_wobble_0.1_triplet_0.5_keep_0.99_None.json"
    check_output_matches(output_dir, 5, expected_name,
                         fixture("test_triplet_pseudoknot_matches_current_state_out.json"))


def test_max_chains_skips_large_input(tmp_path):
    # The input fixture has 6 chains; --max_chains 3 must trigger the early-return
    # skip path in the JSON-input branch so that no output JSON is written.
    output_dir = tmp_path / "out"
    result = run_generator(
        output_dir,
        "--input_json_path", str(fixture("test_json_input_matches_current_state_in.json")),
        "--pdb_id", "TEST", "-N", "2", "--seed", "1", "--max_chains", "3",
    )
    assert not list(output_dir.glob("*.json")), (
        "Expected no output JSON when chain count exceeds --max_chains"
    )
    assert "max_chains=3" in result.stderr


########################## Error-path regression tests (run only these: pytest -k error) ##########################
# These lock the CURRENT intended behaviour of the failure paths so a refactor
# cannot silently change them. main() logs and swallows exceptions (process
# exits 0), so a failed run is observed as "no output JSON written" plus a
# specific log line on stderr.

def test_error_missing_structure_and_predictor(tmp_path):
    # Sequences given, but neither --structure nor --structure_predictor.
    output_dir = tmp_path / "out"
    result = run_generator(
        output_dir,
        "--rna-seq", "ACGU", "--protein-seq", "MK",
        "--pdb_id", "TEST", "-N", "2", "--seed", "1",
    )
    assert not list(output_dir.glob("*.json"))
    assert "Either a structure predictor or a structure must be provided." in result.stderr


def test_error_missing_sequences(tmp_path):
    # No --input_json_path and an incomplete sequence pair (--protein-seq absent).
    output_dir = tmp_path / "out"
    result = run_generator(
        output_dir,
        "--rna-seq", "ACGU", "--structure", "(())",
        "--pdb_id", "TEST", "-N", "2", "--seed", "1",
    )
    assert not list(output_dir.glob("*.json"))
    assert "Either --rna-seq and --protein-seq or --input_json_path must be provided." in result.stderr


def test_error_unknown_predictor(tmp_path):
    # An unrecognised --structure_predictor must be rejected with a clear error.
    output_dir = tmp_path / "out"
    result = run_generator(
        output_dir,
        "--rna-seq", "ACGU", "--protein-seq", "MK",
        "--pdb_id", "TEST", "-N", "2", "--seed", "1",
        "--structure_predictor", "foobar",
    )
    assert not list(output_dir.glob("*.json"))
    assert "Unknown structure predictor: foobar" in result.stderr


def test_error_invalid_input_json_path(tmp_path):
    # A missing --input_json_path file surfaces as a load failure and no output.
    output_dir = tmp_path / "out"
    missing = tmp_path / "does_not_exist.json"
    result = run_generator(
        output_dir,
        "--input_json_path", str(missing),
        "--pdb_id", "TEST", "-N", "2", "--seed", "1",
    )
    assert not list(output_dir.glob("*.json"))
    assert "Failed to load JSON" in result.stderr


def test_error_unbalanced_structure_recovers(tmp_path):
    # An extra closing bracket is logged as an error but is NOT fatal: pair_indices
    # skips the unmatched closer and generation still produces output. This locks
    # the current lenient recovery so a refactor doesn't turn it into a hard fail.
    output_dir = tmp_path / "out"
    result = run_generator(
        output_dir,
        "--rna-seq", "ACGUA", "--protein-seq", "MK",
        "--pdb_id", "TEST", "-N", "2", "--seed", "1",
        "--structure", "(()))",
    )
    assert "Unbalanced structure at position 4" in result.stderr
    output_json = get_output_json(output_dir)
    sanity_check(output_json, 2)

