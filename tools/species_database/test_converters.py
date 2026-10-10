"""Converts the test excerpts of NASA's thermo.inp, Cantera's airNASA9.yaml
and a GDTk-layout file, and checks that each gives the same molar masses,
charges, temperature ranges, blend widths and coefficients as the shipped
five-species air file.

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
    "GDTk layout": [tools / "gdtk_to_json.py", "--database", data / "gdtk", "--species", "N2"],
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
            print(f"{label} {name}: checked")

# Reaction mechanism converter: CGS to SI by (1e-6)^(order - 1), third-body
# efficiencies, and the default backward rate from the equilibrium constant.
with tempfile.TemporaryDirectory() as tmp:
    out = Path(tmp) / "reactions.json"
    subprocess.run([sys.executable, str(tools / "gdtk_reactions_to_json.py"), str(data / "gdtk_mechanism_sample.inp"),
                    "--source", "test", "-o", str(out)], check=True)
    mechanism = json.loads(out.read_text())
    checks = [
        (mechanism["format"], "gaslib-reactions-1"),
        (mechanism["temperature_limits"], {"lower": 300.0, "upper": 30000.0}),
        (len(mechanism["reactions"]), 3),
        (mechanism["reactions"][0]["reactants"], {"N2": 2.0}),
        (mechanism["reactions"][0]["products"], {"N": 2.0, "N2": 1.0}),
        (mechanism["reactions"][0]["forward"], {"model": "arrhenius", "A": 1.216e14, "n": -1.214,
                                                "activation_temperature": 113200.0}),
        (mechanism["reactions"][0]["backward"], {"model": "equilibrium"}),
        (mechanism["reactions"][1]["forward"]["A"], 1.45e9),
        (mechanism["reactions"][1]["third_body"]["efficiencies"]["NO"], 0.664827586),
        (mechanism["reactions"][2]["backward"], {"model": "arrhenius", "A": 1.6e7, "n": 0.1,
                                                 "activation_temperature": 0.0}),
    ]
    for got, expected in checks:
        if got != expected:
            failures += 1
            print(f"reaction converter: {got} != {expected}")
    print("reaction converter: checked")

sys.exit(1 if failures else 0)
