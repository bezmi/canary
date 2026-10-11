"""Converts the test excerpts of NASA's thermo.inp, Cantera's airNASA9.yaml
and a GDTk-layout file, and checks that each gives the same molar masses,
charges, temperature ranges, blend widths and coefficients as the shipped
five-species air file (with its transport data for the GDTk layout).
Converts the GDTk-layout reaction-mechanism, energy-exchange and
collision-integral samples and checks the converted values.

Usage: test_converters.py <gaslib source directory>
"""

import json
import subprocess
import sys
import tempfile
from pathlib import Path

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
    "GDTk layout": [tools / "gdtk_to_json.py", "--database", data / "gdtk", "--species", "N2", "--transport"],
}

failures = 0
with tempfile.TemporaryDirectory() as tmp:
    for label, command in runs.items():
        out = Path(tmp) / "out.json"
        subprocess.run([sys.executable, *map(str, command), "--source", "test", "-o", str(out)], check=True)
        converted = json.loads(out.read_text())
        assert converted["format"] == "gaslib-species-1"
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
            if "transport" in record and record["transport"] != expected["transport"]:
                failures += 1
                print(f"{label} {name}: transport differs")
            print(f"{label} {name}: checked")

# Reaction mechanism converter: CGS to SI by (1e-6)^(order - 1), third-body
# efficiencies, the default backward rate from the equilibrium constant, and
# Park rates as Arrhenius rates with a translational exponent.
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "reactions.json"
    subprocess.run([sys.executable, str(tools / "gdtk_reactions_to_json.py"), str(data / "gdtk_mechanism_sample.inp"),
                    "--source", "test", "-o", str(out)], check=True)
    mechanism = json.loads(out.read_text())
    checks = [
        (mechanism["format"], "gaslib-reactions-1"),
        (mechanism["temperature_limits"], {"lower": 300.0, "upper": 30000.0}),
        (len(mechanism["reactions"]), 4),
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
    ]
    for got, expected in checks:
        if got != expected:
            failures += 1
            print(f"reaction converter: {got} != {expected}")
    print("reaction converter: checked")

# Energy-exchange converter: V-T entries with and without Park's
# correction, and an impartial-dissociation C-V entry.
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "exchange.json"
    subprocess.run([sys.executable, str(tools / "gdtk_exchange_to_json.py"), str(data / "gdtk_exchange_sample.lua"),
                    "--source", "test", "-o", str(out)], check=True)
    exchange = json.loads(out.read_text())
    checks = [
        (exchange["format"], "gaslib-energy-exchange-1"),
        (exchange["units"], {"energy": "J/mol", "temperature": "K", "area": "m^2"}),
        (exchange["vibration_translation"], [
            {"species": "N2", "partner": "N", "model": "landau_teller",
             "relaxation_time": {"model": "millikan_white", "a": 180.0, "b": 0.0262,
                                 "high_temperature_correction": {"model": "park_1993", "cross_section": 1.0e-21,
                                                                 "reference_temperature": 50000.0}}},
            {"species": "O2", "partner": "O2", "model": "landau_teller",
             "relaxation_time": {"model": "millikan_white", "a": 138.0, "b": 0.03}},
        ]),
        (exchange["chemistry_vibration"], [
            {"reaction": 0, "species": "N2", "model": "impartial_dissociation",
             "dissociation_energy": 945360.0, "vibrational_temperature": 3393.44},
        ]),
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
        (converted["format"], "gaslib-collision-integrals-1"),
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
