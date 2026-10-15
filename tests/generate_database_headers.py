"""Generate source-valid database headers for the C++ construction tests."""

import copy
import json
from pathlib import Path
import subprocess
import sys

root, output = map(Path, sys.argv[1:])
output.mkdir(parents=True, exist_ok=True)
tools = root / "tools/species_database"
data = root / "database/data"
names = ["N2", "O2", "N", "O", "NO"]


def generate(tool, *arguments):
    subprocess.run([sys.executable, str(tools / tool), *map(str, arguments)], check=True, timeout=30)


def species(source, name):
    generate("species_to_header.py", source, "--species", *names, "--name", name,
             "--thermo", "--electronic", "-o", output / f"{name}.hpp")


species(data / "air5.json", "generated_species")
generate("transport_to_header.py", data / "air5.json", data / "air5_collision_integrals.json",
         "--species", *names, "--name", "generated", "--segments", "3", "-o", output / "generated_transport.hpp")
generate("reactions_to_header.py", data / "air5_park.json", "--name", "generated_reactions",
         "--alias", "GeneratedReactions", "--species-header", "generated_species.hpp",
         "--species-array", "generated_species_species", "-o", output / "generated_reactions.hpp")
generate("exchange_to_header.py", data / "air5_park_exchange.json", "--name", "generated_exchange",
         "--alias", "GeneratedExchange", "--species-header", "generated_species.hpp",
         "--species-array", "generated_species_species", "--mechanism-header", "generated_reactions.hpp",
         "--mechanism-alias", "GeneratedReactions", "-o", output / "generated_exchange.hpp")

original = json.loads((data / "air5.json").read_text())
for name in ("negative_capacity", "small_mass"):
    fixture = copy.deepcopy(original)
    record = fixture["species"]["N2"]
    if name == "negative_capacity":
        # Finite coefficients and valid dimensions do not establish positive cp.
        record["thermo"].update(temperature_ranges=[200.0, 20000.0], blend_widths=[],
                               coefficients=[[0.0, 0.0, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]])
    else:
        # The positive SI mass becomes unrepresentable in a float SI model.
        record["molar_mass"] = 1e-30
    source = output / f"{name}.json"
    source.write_text(json.dumps(fixture))
    species(source, name)
