"""Converts NASA9 species from a Cantera YAML file to a gaslib species file.

Cantera files give each species' composition but not its molar mass, so the
atomic weights (g/mol) are required: --atomic-weights N=14.0067 O=15.9994 E=5.48579909e-4
An electron count E n gives charge -n. The blend width of each interior
boundary is required as for thermo_inp_to_json.py.

Usage: cantera_to_json.py air.yaml --species N2 O2 --atomic-weights ... --blend ... --source TEXT -o out.json
"""

import argparse
from decimal import Decimal

import yaml

import species_json


def convert(entry, weights, widths):
    name = entry["name"]
    thermo = entry["thermo"]
    if thermo.get("model") != "NASA9":
        raise SystemExit(f"{name}: thermo model is {thermo.get('model')}, not NASA9")
    composition = entry["composition"]
    missing = [element for element in composition if element not in weights]
    if missing:
        raise SystemExit(f"{name}: no atomic weight for {', '.join(missing)}")
    grams = sum(Decimal(str(count)) * weights[element] for element, count in composition.items())
    ranges = [float(t) for t in thermo["temperature-ranges"]]
    coefficients = [[float(a) for a in row] for row in thermo["data"]]
    return {
        "molar_mass": float(grams / Decimal(1000)),
        "charge": -float(composition["E"]) if "E" in composition else 0.0,
        "thermo": species_json.nasa9_thermo(
            ranges, coefficients, species_json.blend_widths_for(name, ranges, widths), thermo.get("note")
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("yaml_file")
    parser.add_argument("--species", required=True, nargs="+")
    parser.add_argument("--atomic-weights", required=True, nargs="+", help="element=g/mol")
    parser.add_argument("--blend", required=True, nargs="+", help="boundary=width in K")
    parser.add_argument("--source", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    weights = {k: Decimal(v) for k, v in (item.split("=") for item in args.atomic_weights)}
    widths = species_json.parse_blend_widths(args.blend)
    with open(args.yaml_file, encoding="utf-8") as stream:
        # BaseLoader keeps every scalar a string: YAML 1.1 would read the
        # species name NO as false.
        data = yaml.load(stream, Loader=yaml.BaseLoader)
    entries = {entry["name"]: entry for entry in data["species"]}
    species = {}
    for name in args.species:
        if name not in entries:
            raise SystemExit(f"{name}: not in {args.yaml_file}")
        species[name] = convert(entries[name], weights, widths)
    species_json.dump(species_json.species_file(species, reference=args.source), args.output)


if __name__ == "__main__":
    main()
