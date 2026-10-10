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

sys.exit(1 if failures else 0)
