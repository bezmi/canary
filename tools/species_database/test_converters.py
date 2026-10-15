"""Converts the test excerpts of NASA's thermo.inp, Cantera's airNASA9.yaml
and a GDTk-layout file, and checks that each gives the same molar masses,
charges, temperature ranges, blend widths and coefficients as the shipped
five-species air file (with its transport data for the GDTk layout).
Converts the GDTk-layout reaction-mechanism, energy-exchange and
collision-integral samples and checks the converted values.

Usage: test_converters.py <canary source directory>
"""

import copy
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import yaml
import jsonschema

import species_json

root = Path(sys.argv[1])
tools = root / "tools" / "species_database"
data = root / "tests" / "data" / "species"
shipped = json.loads((root / "database" / "data" / "air5.json").read_text())["species"]
blend = ["--blend", "1000=400", "6000=1000"]

runs = {
    "thermo.inp": [tools / "thermo_inp_to_json.py", data / "thermo_inp_excerpt.inp",
                   "--species", "N2", "O2", "NO", "--roles", "N2=heavy", "O2=heavy", "NO=heavy", *blend],
    "Cantera": [tools / "cantera_to_json.py", data / "cantera_excerpt.yaml", "--species", "N2", "O2", "NO",
                "--roles", "N2=heavy", "O2=heavy", "NO=heavy",
                "--atomic-weights", "N=14.0067", "O=15.9994", *blend],
    "GDTk layout": [tools / "gdtk_to_json.py", "--database", data / "gdtk", "--species", "N2", "--transport",
                    "--electronic", "--roles", "N2=heavy"],
}

expected_species = {
    "thermo.inp": {"N2", "O2", "NO"},
    "Cantera": {"N2", "O2", "NO"},
    "GDTk layout": {"N2"},
}

failures = 0
with tempfile.TemporaryDirectory() as tmp:
    for label, command in runs.items():
        out = Path(tmp) / "out.json"
        subprocess.run([sys.executable, *map(str, command), "--source", "test", "-o", str(out)], check=True)
        converted = json.loads(out.read_text())
        assert converted["format"] == "canary-species-2"
        assert set(converted["species"]) == expected_species[label], f"{label}: unexpected species set"
        if label == "GDTk layout":
            assert converted["species"]["N2"]["electronic_levels"] == shipped["N2"]["electronic_levels"]
            assert converted["species"]["N2"]["transport"] == shipped["N2"]["transport"]
        for name, record in converted["species"].items():
            expected = shipped[name]
            for key in ("molar_mass", "charge", "role", "elements"):
                if record[key] != expected[key]:
                    failures += 1
                    print(f"{label} {name}: {key} {record[key]} != {expected[key]}")
            for key in ("model", "temperature_ranges", "blend_widths", "coefficients"):
                if record["thermo"][key] != expected["thermo"][key]:
                    failures += 1
                    print(f"{label} {name}: thermo {key} differs")
            if "electronic_levels" in record and record["electronic_levels"] != expected["electronic_levels"]:
                failures += 1
                print(f"{label} {name}: electronic levels differ")
            if "transport" in record and record["transport"] != expected["transport"]:
                failures += 1
                print(f"{label} {name}: transport differs")
            print(f"{label} {name}: checked")

# Duplicate atomic weights and species names must be rejected before output.
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "duplicate.json"
    command = list(runs["Cantera"])
    command.insert(command.index("--blend"), "N=14.0067")
    result = subprocess.run([sys.executable, *map(str, command), "--source", "test", "-o", str(out)],
                            capture_output=True, text=True)
    assert result.returncode != 0 and "atomic weight for N: given twice" in result.stderr
    assert not out.exists()

    with (data / "cantera_excerpt.yaml").open(encoding="utf-8") as stream:
        records = yaml.load(stream, Loader=yaml.BaseLoader)
    records["species"].append(records["species"][0])
    duplicate = Path(tmp) / "duplicate.yaml"
    duplicate.write_text(yaml.safe_dump(records), encoding="utf-8")
    command = list(runs["Cantera"])
    command[1] = duplicate
    result = subprocess.run([sys.executable, *map(str, command), "--source", "test", "-o", str(out)],
                            capture_output=True, text=True)
    assert result.returncode != 0 and "species N2 is given twice" in result.stderr
    assert str(duplicate) in result.stderr
    assert not out.exists()

def reject_conversion(command, out, diagnostic):
    result = subprocess.run([sys.executable, *map(str, command), "--source", "test", "-o", str(out)],
                            capture_output=True, text=True)
    assert result.returncode != 0 and diagnostic in result.stderr, result.stderr
    assert "Traceback" not in result.stderr, result.stderr
    assert not out.exists()


# Every importer requires one unambiguous declaration for each selected species.
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "roles.json"
    for label, original in runs.items():
        command = list(original)
        start = command.index("--roles")
        end = start + 1
        while end < len(command) and not str(command[end]).startswith("--"):
            end += 1
        del command[start:end]
        reject_conversion(command, out, "--roles")
        valid = [name + "=heavy" for name in sorted(expected_species[label])]
        for declarations, diagnostic in (
            (valid + [valid[0]], "given twice"),
            (valid + ["absent=heavy"], "species is not selected"),
            (["N2=unspecified"], "expected NAME=heavy or NAME=electron"),
            (["N2=unknown"], "expected NAME=heavy or NAME=electron"),
            (["N2"], "expected NAME=heavy or NAME=electron"),
            (["=heavy"], "expected NAME=heavy or NAME=electron"),
            (["N2=heavy=electron"], "expected NAME=heavy or NAME=electron"),
        ):
            reject_conversion(command + ["--roles", *declarations], out, diagnostic)
        if len(valid) > 1:
            reject_conversion(command + ["--roles", *valid[1:]], out, "missing declaration")
        invalid_charge = ["N2=electron" if role == "N2=heavy" else role for role in valid]
        reject_conversion(command + ["--roles", *invalid_charge], out, "electron role requires charge -1")
print("Explicit converter roles: checked")


# Explicit roles survive conversion independently of labels and source ordering.
with tempfile.TemporaryDirectory() as tmp:
    source = Path(tmp) / "roles.yaml"
    out = Path(tmp) / "roles.json"
    thermo = {"model": "NASA9", "temperature-ranges": [200, 20000],
              "data": [[0, 0, 2.5, 0, 0, 0, 0, 0, 0]]}
    source.write_text(yaml.safe_dump({"species": [
        {"name": "carrier", "composition": {"E": 1}, "thermo": thermo},
        {"name": "e-", "composition": {"N": 1}, "thermo": thermo},
    ]}), encoding="utf-8")
    command = [sys.executable, str(tools / "cantera_to_json.py"), str(source),
               "--species", "e-", "carrier", "--roles", "carrier=electron", "e-=heavy",
               "--atomic-weights", "N=14.0067", "E=0.000548579909", *blend,
               "--source", "test", "-o", str(out)]
    subprocess.run(command, check=True)
    records = json.loads(out.read_text())["species"]
    assert records["e-"]["role"] == "heavy" and records["e-"]["charge"] == 0
    assert records["carrier"]["role"] == "electron" and records["carrier"]["charge"] == -1
    assert records["carrier"]["elements"] == {}
    assert records["e-"]["elements"] == {"N": 1}

    # A GDTk-layout electron with a user-chosen label and a constant cp/R fit.
    (Path(tmp) / "carrier.lua").write_text(
        "db.carrier.charge = -1\n"
        "db.carrier.M = {value = 5.48579909e-7}\n"
        "db.carrier.ceaThermoCoeffs = {nsegments = 1, T_break_points = {200, 20000}, "
        "segment0 = {0, 0, 2.5, 0, 0, 0, 0, 0, 0}}\n", encoding="utf-8")
    subprocess.run([sys.executable, str(tools / "gdtk_to_json.py"), "--database", tmp,
                    "--species", "carrier", "--roles", "carrier=electron",
                    "--source", "test", "-o", str(out)], check=True)
    record = json.loads(out.read_text())["species"]["carrier"]
    assert record["role"] == "electron" and record["charge"] == -1

    # A complete fixed-width electron record in thermo.inp layout.
    metadata = " 1" + " " * 8 + "E " + "  1.00" + " " * 34 + "0.0005485799".rjust(13)
    powers = "".join(f"{p:5.1f}" for p in [-2, -1, 0, 1, 2, 3, 4])
    head = f"{200:11.3f}{20000:11.3f}7" + powers
    first = "".join(f"{v:16.8E}" for v in [0, 0, 2.5, 0, 0])
    second = "".join(f"{v:16.8E}" for v in [0, 0, 0, 0, 0])
    source = Path(tmp) / "electron.inp"
    source.write_text("\n".join(["carrier".ljust(18), metadata, head, first, second]) + "\n")
    subprocess.run([sys.executable, str(tools / "thermo_inp_to_json.py"), str(source),
                    "--species", "carrier", "--roles", "carrier=electron", *blend,
                    "--source", "test", "-o", str(out)], check=True)
    record = json.loads(out.read_text())["species"]["carrier"]
    assert record["role"] == "electron" and record["charge"] == -1
print("Renamed electron conversion: checked")


# Reject invalid masses and dimensions before creating an output file.
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "invalid.json"
    for weight in ("0", "-1", "NaN", "sNaN", "Infinity", "-Infinity", "invalid"):
        command = ["N=" + weight if str(arg) == "N=14.0067" else arg for arg in runs["Cantera"]]
        reject_conversion(command, out, "atomic weight for N: must be finite and positive")
    records = yaml.load((data / "cantera_excerpt.yaml").read_text(), Loader=yaml.BaseLoader)
    for count in ("0", "-1", "NaN", "sNaN", "Infinity", "1e400", "1e999999", "1e-400", "invalid"):
        changed = copy.deepcopy(records)
        changed["species"][0]["composition"] = {"N": count}
        source = Path(tmp) / "invalid.yaml"
        source.write_text(yaml.safe_dump(changed))
        command = list(runs["Cantera"])
        command[1] = source
        reject_conversion(command, out, "N2: molar mass must be finite and positive")
    for field, value, diagnostic in (
        ("temperature-ranges", [], "requires at least two temperature boundaries"),
        ("temperature-ranges", [200], "requires at least two temperature boundaries"),
        ("data", [], "coefficient rows must match"),
        ("data", [[0] * 9] * 2, "coefficient rows must match"),
        ("data", [[0] * 9] * 4, "coefficient rows must match"),
        ("data", [[0] * 8] * 3, "requires nine coefficients"),
        ("data", [[0] * 10] * 3, "requires nine coefficients"),
    ):
        changed = copy.deepcopy(records)
        changed["species"][0]["thermo"][field] = value
        source = Path(tmp) / "invalid.yaml"
        source.write_text(yaml.safe_dump(changed))
        command = list(runs["Cantera"])
        command[1] = source
        reject_conversion(command, out, "N2: NASA9 " + diagnostic)
    # Signed electron counts in ion compositions remain supported.
    changed = copy.deepcopy(records)
    changed["species"] = [changed["species"][0]]
    changed["species"][0]["composition"] = {"N": "2", "E": "-1"}
    source = Path(tmp) / "ion.yaml"
    source.write_text(yaml.safe_dump(changed))
    command = [tools / "cantera_to_json.py", source, "--species", "N2",
               "--roles", "N2=heavy", "--atomic-weights", "N=14.0067", "E=0.000548579909", *blend]
    subprocess.run([sys.executable, *map(str, command), "--source", "test", "-o", str(out)], check=True)
    record = json.loads(out.read_text())["species"]["N2"]
    assert record["charge"] == 1.0
    assert record["molar_mass"] == species_json.kg_per_mol("28.012851420091")
print("Cantera mass and NASA9 dimensions: checked")

for widths in ([], [1, 2]):
    try:
        species_json.nasa9_thermo([200, 1000, 6000], [[0] * 9] * 2, widths)
    except ValueError as error:
        assert "blend widths must match" in str(error)
    else:
        raise AssertionError("invalid blend-width count accepted")

# Check every record truncation and each fixed-width polynomial exponent.
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "invalid.json"
    original = (data / "thermo_inp_excerpt.inp").read_text().splitlines()
    first = next(i for i, line in enumerate(original) if line[:18].split()[:1] == ["N2"])
    lines = original[first:first + 11]
    source = Path(tmp) / "invalid.inp"
    command = [tools / "thermo_inp_to_json.py", source, "--species", "N2", "--roles", "N2=heavy", *blend]
    for length in range(1, len(lines)):
        source.write_text("\n".join(lines[:length]) + "\n")
        reject_conversion(command, out, "N2: incomplete")
    for count in (" 0", "-1", "xx"):
        changed = list(lines)
        changed[1] = count + changed[1][2:]
        source.write_text("\n".join(changed) + "\n")
        reject_conversion(command, out, "temperature interval count")
    for record, length in ((1, 64), (2, 57), (3, 79), (4, 79)):
        changed = list(lines)
        changed[record] = changed[record][:length]
        source.write_text("\n".join(changed) + "\n")
        reject_conversion(command, out, "incomplete")
    for start, end, replacement in [(22, 23, "6"), *[(23 + 5*j, 28 + 5*j, "  5.0") for j in range(7)]]:
        changed = list(lines)
        changed[2] = changed[2][:start] + replacement + changed[2][end:]
        source.write_text("\n".join(changed) + "\n")
        reject_conversion(command, out, "N2: interval 1: NASA9 requires seven powers")
    changed = list(lines)
    changed[3] = "invalid".ljust(16) + changed[3][16:]
    source.write_text("\n".join(changed) + "\n")
    reject_conversion(command, out, "N2: interval 1: invalid NASA9 numeric field")
    # D/d exponents are accepted without an intermediate binary rounding of mass.
    for exponent in ("D", "d"):
        changed = [line.replace("D", exponent) for line in lines]
        changed[1] = changed[1][:52] + ("2.80134" + exponent + "+01").rjust(13) + changed[1][65:]
        changed[2] = ("2.000" + exponent + "+02").rjust(11) + changed[2][11:]
        source.write_text("\n".join(changed) + "\n")
        subprocess.run([sys.executable, *map(str, command), "--source", "test", "-o", str(out)], check=True)
        record = json.loads(out.read_text())["species"]["N2"]
        assert record["molar_mass"] == shipped["N2"]["molar_mass"]
        for key in ("temperature_ranges", "blend_widths", "coefficients"):
            assert record["thermo"][key] == shipped["N2"]["thermo"][key]
        out.unlink()
print("NASA9 complete records, powers and Fortran numbers: checked")

# GDTk layout with the bracket form of a species name (N+ in N_plus.lua),
# with its electronic levels converted to temperatures: equal to the shipped
# eleven-species air file.
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "ion.json"
    subprocess.run([sys.executable, str(tools / "gdtk_to_json.py"), "--database", str(data / "gdtk"),
                    "--species", "N+", "--roles", "N+=heavy", "--electronic", "--source", "test", "-o", str(out)], check=True)
    converted = json.loads(out.read_text())["species"]["N+"]
    expected = json.loads((root / "database" / "data" / "air11.json").read_text())["species"]["N+"]
    for key in ("molar_mass", "charge", "role", "electronic_levels"):
        if converted[key] != expected[key]:
            failures += 1
            print(f"GDTk bracket form N+: {key} {converted[key]} != {expected[key]}")
    for key in ("temperature_ranges", "blend_widths", "coefficients"):
        if converted["thermo"][key] != expected["thermo"][key]:
            failures += 1
            print(f"GDTk bracket form N+: thermo {key} differs")
    print("GDTk bracket form N+: checked")

# Reaction mechanism converter: CGS to SI by (1e-6)^(order - 1), third-body
# efficiencies, the default backward rate from the equilibrium constant, and
# Park rates as Arrhenius rates with a translational exponent.
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "reactions.json"
    subprocess.run([sys.executable, str(tools / "gdtk_reactions_to_json.py"), str(data / "gdtk_mechanism_sample.inp"),
                    "--source", "test", "-o", str(out)], check=True)
    mechanism = json.loads(out.read_text())
    checks = [
        (mechanism["format"], "canary-reactions-1"),
        (mechanism["temperature_limits"], {"lower": 300.0, "upper": 30000.0}),
        (len(mechanism["reactions"]), 5),
        (mechanism["reactions"][0]["reactants"], {"N2": 2.0}),
        (mechanism["reactions"][0]["products"], {"N": 2.0, "N2": 1.0}),
        (mechanism["reactions"][0]["forward"], {"model": "arrhenius", "A": 1.216e14, "n": -1.214,
                                                "activation_temperature": 113200.0}),
        (mechanism["reactions"][0]["backward"], {"model": "equilibrium"}),
        (mechanism["reactions"][1]["forward"]["A"], 1.45e9),
        (mechanism["reactions"][1]["third_body"]["efficiencies"]["NO"], 0.664827586),
        (mechanism["reactions"][2]["backward"], {"model": "arrhenius", "A": 1.6e7, "n": 0.1,
                                                 "activation_temperature": 0.0}),
        (mechanism["reactions"][3]["forward"], {"model": "arrhenius", "A": 2.0e15, "n": -1.5,
                                                "activation_temperature": 59500.0,
                                                "translational_exponent": 0.5}),
        ("translational_exponent" in mechanism["reactions"][0]["forward"], False),
        (mechanism["reactions"][4]["reactants"], {"N": 1.0, "e-": 1.0}),
        (mechanism["reactions"][4]["products"], {"N+": 1.0, "e-": 2.0}),
        (mechanism["reactions"][4]["forward"], {"model": "arrhenius", "A": 2.5e28, "n": -3.82,
                                                "activation_temperature": 168600.0, "translational_exponent": 0.0,
                                                "other_temperature": "electron"}),
        (mechanism["reactions"][4]["backward"], {"model": "equilibrium", "other_temperature": "electron"}),
    ]
    for got, expected in checks:
        if got != expected:
            failures += 1
            print(f"reaction converter: {got} != {expected}")
    print("reaction converter: checked")

# Energy-exchange converter: V-T entries without and with the Park (1993),
# Kim and Jo and Park (1985) corrections, an impartial-dissociation C-V
# entry (and one with negative dissociation energy, dropped), elastic
# electron exchange with a neutral and an ion, and Candler's V-E.
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "exchange.json"
    subprocess.run([sys.executable, str(tools / "gdtk_exchange_to_json.py"), str(data / "gdtk_exchange_sample.lua"),
                    "--source", "test", "-o", str(out)], check=True)
    exchange = json.loads(out.read_text())
    checks = [
        (exchange["format"], "canary-energy-exchange-1"),
        (exchange["units"], {"energy": "J/mol", "temperature": "K", "area": "m^2"}),
        (exchange["vibration_translation"], [
            {"species": "N2", "partner": "N", "model": "landau_teller",
             "relaxation_time": {"model": "millikan_white", "a": 180.0, "b": 0.0262,
                                 "high_temperature_correction": {"model": "park_1993", "cross_section": 1.0e-21,
                                                                 "reference_temperature": 50000.0}}},
            {"species": "O2", "partner": "O2", "model": "landau_teller",
             "relaxation_time": {"model": "millikan_white", "a": 138.0, "b": 0.03}},
            {"species": "N2", "partner": "N2", "model": "landau_teller",
             "relaxation_time": {"model": "millikan_white", "a": 221.53, "b": 0.029,
                                 "high_temperature_correction": {"model": "kim_2021", "cross_section": 1.7e-23,
                                                                 "exponent": 0.55}}},
            {"species": "NO", "partner": "N", "model": "landau_teller",
             "relaxation_time": {"model": "millikan_white", "a": 49.5, "b": 0.042,
                                 "high_temperature_correction": {"model": "park_1985", "cross_section": 5.8e-20}}},
        ]),
        (exchange["chemistry_vibration"], [
            {"reaction": 0, "species": "N2", "model": "impartial_dissociation",
             "dissociation_energy": 945360.0, "vibrational_temperature": 3393.44},
        ]),
        (exchange["electron_translation"], [
            {"partner": "N2", "model": "elastic",
             "cross_section": {"model": "gnoffo_neutral", "a": 7.5e-20, "b": 5.5e-24, "c": -1e-28,
                               "maximum_temperature": 30000.0}},
            {"partner": "N+", "model": "elastic", "cross_section": {"model": "coulomb"}},
        ]),
        (exchange["vibration_electron"], [
            {"species": "N2", "model": "landau_teller",
             "relaxation_time": {"model": "candler", "a_low": 7.5, "b_low": -57.0, "c_low": 98.7,
                                 "a_high": 2.36, "b_high": -17.9, "c_high": 24.35,
                                 "blend_lower": 6950.0, "blend_upper": 7050.0}},
        ]),
        (exchange["electron_impact_ionisation"], []),
    ]
    for got, expected in checks:
        if got != expected:
            failures += 1
            print(f"exchange converter: {got} != {expected}")
    print("exchange converter: checked")

# Collision-integral converter: D shifted to m^2 by ln(1e-20), pairs in file
# order, equal to the shipped five-species file for the same pairs.
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "collision_integrals.json"
    subprocess.run([sys.executable, str(tools / "gdtk_collision_integrals_to_json.py"),
                    str(data / "gdtk_collision_integrals_sample.lua"), "--species", "N2", "O2",
                    "--source", "test", "-o", str(out)], check=True)
    converted = json.loads(out.read_text())
    shipped_pairs = json.loads((root / "database" / "data" / "air5_collision_integrals.json").read_text())["pairs"]
    checks = [
        (converted["format"], "canary-collision-integrals-1"),
        (converted["units"], {"temperature": "K", "area": "m^2"}),
        ([p["species"] for p in converted["pairs"]], [["N2", "N2"], ["O2", "N2"], ["O2", "O2"]]),
        (converted["pairs"], shipped_pairs[:3]),
    ]
    for got, expected in checks:
        if got != expected:
            failures += 1
            print(f"collision-integral converter: {got} != {expected}")
    print("collision-integral converter: checked")



# Species header generation preserves all three composition states.
with tempfile.TemporaryDirectory() as tmp:
    file = json.loads((root / "database/data/air11.json").read_text())
    del file["species"]["N"]["elements"]
    file["species"]["O"]["elements"] = {}
    source = Path(tmp) / "species.json"
    header = Path(tmp) / "species.hpp"
    source.write_text(json.dumps(file))
    command = [sys.executable, str(tools / "species_to_header.py"), str(source),
               "--species", "N2", "N", "O", "e-", "--name", "example", "--thermo", "-o", str(header)]
    subprocess.run(command, check=True)
    text = header.read_text()
    assert 'example_elements_0' in text and '{"N", 2}' in text
    assert 'example_elements_1' not in text
    assert 'std::array<ElementCount, 0> example_elements_2' in text
    assert 'std::array<ElementCount, 0> example_elements_3' in text
    source.write_text(json.dumps(file).replace('"elements": {"N": 2}', '"elements": {"N": 2, "\\u004e": 3}', 1))
    result = subprocess.run(command, text=True, capture_output=True)
    assert result.returncode and 'duplicate JSON key' in result.stderr
    result = subprocess.run([sys.executable, str(tools / "check_species_files.py"),
                             str(root / "database/species.schema.json"), str(source)],
                            text=True, capture_output=True)
    assert result.returncode and 'duplicate JSON key' in result.stderr
print("Species elemental composition and generated metadata: checked")

schema = json.loads((root / "database/species.schema.json").read_text())
validator = jsonschema.Draft202012Validator(schema)
file = json.loads((root / "database/data/air5.json").read_text())
for invalid in ({"N": 0}, {"N": -1}, {"N": 1.5}, {"N": True}, {"N": "2"},
                {"N": 4294967296}, {"": 1}, []):
    file["species"]["N2"]["elements"] = invalid
    assert list(validator.iter_errors(file)), invalid
    try:
        species_json.species_file(file["species"])
    except SystemExit:
        pass
    else:
        raise AssertionError(f"invalid element composition accepted: {invalid}")
for valid in ({}, {"N": 2}, {"N": 4294967295}, None):
    file["species"]["N2"]["elements"] = valid
    assert not list(validator.iter_errors(file)), valid
print("Element schema and import validation boundaries: checked")

# Every JSON converter shares the regular-file output contract.
output_runs = [*runs.values(),
               [tools / "gdtk_reactions_to_json.py", data / "gdtk_mechanism_sample.inp"],
               [tools / "gdtk_exchange_to_json.py", data / "gdtk_exchange_sample.lua"],
               [tools / "gdtk_collision_integrals_to_json.py", data / "gdtk_collision_integrals_sample.lua",
                "--species", "N2", "O2"]]
with tempfile.TemporaryDirectory() as tmp:
    directory = Path(tmp)
    original = directory / "original.json"
    original.write_text('{"original": true}\n', encoding="utf-8")
    out = directory / "output.json"
    out.symlink_to(original)
    for command in output_runs:
        result = subprocess.run([sys.executable, *map(str, command), "--source", "test", "-o", str(out)],
                                capture_output=True, text=True)
        assert result.returncode and "symlinks are not supported" in result.stderr, result.stderr
        assert str(out) in result.stderr
        assert out.is_symlink()
        assert original.read_text() == '{"original": true}\n'
        assert not list(directory.glob(".canary-*.tmp"))
print("JSON converter destination checks: checked")

sys.exit(1 if failures else 0)
