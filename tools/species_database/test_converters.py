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

import species_json

root = Path(sys.argv[1])
tools = root / "tools" / "species_database"
data = root / "tests" / "data" / "species"
shipped = json.loads((root / "database" / "data" / "air5.json").read_text())["species"]
blend = ["--blend", "1000=400", "6000=1000"]

runs = {
    "thermo.inp": [tools / "thermo_inp_to_json.py", data / "thermo_inp_excerpt.inp",
                   "--species", "N2", "O2", "NO", *blend],
    "Cantera": [tools / "cantera_to_json.py", data / "cantera_excerpt.yaml", "--species", "N2", "O2", "NO",
                "--atomic-weights", "N=14.0067", "O=15.9994", *blend],
    "GDTk layout": [tools / "gdtk_to_json.py", "--database", data / "gdtk", "--species", "N2", "--transport",
                    "--electronic"],
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
        assert converted["format"] == "canary-species-1"
        assert set(converted["species"]) == expected_species[label], f"{label}: unexpected species set"
        if label == "GDTk layout":
            assert converted["species"]["N2"]["electronic_levels"] == shipped["N2"]["electronic_levels"]
            assert converted["species"]["N2"]["transport"] == shipped["N2"]["transport"]
        for name, record in converted["species"].items():
            expected = shipped[name]
            for key in ("molar_mass", "charge"):
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
               "--atomic-weights", "N=14.0067", "E=0.000548579909", *blend]
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
    command = [tools / "thermo_inp_to_json.py", source, "--species", "N2", *blend]
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
                    "--species", "N+", "--electronic", "--source", "test", "-o", str(out)], check=True)
    converted = json.loads(out.read_text())["species"]["N+"]
    expected = json.loads((root / "database" / "data" / "air11.json").read_text())["species"]["N+"]
    for key in ("molar_mass", "charge", "electronic_levels"):
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

sys.exit(1 if failures else 0)
