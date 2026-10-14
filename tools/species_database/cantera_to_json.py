"""Converts NASA9 species from a Cantera YAML file to a canary species file.

Cantera files give each species' composition but not its molar mass, so the
atomic weights (g/mol) are required: --atomic-weights N=14.0067 O=15.9994 E=5.48579909e-4
An electron count E n gives charge -n. The blend width of each interior
boundary is required as for thermo_inp_to_json.py.

Usage: cantera_to_json.py air.yaml --species N2 O2 --atomic-weights ... --blend ... --source TEXT -o out.json
"""

import argparse
import math
from decimal import Decimal, DecimalException, InvalidOperation

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
    try:
        grams = sum(Decimal(str(count)) * weights[element] for element, count in composition.items())
        mass = float(grams / Decimal(1000))
    except (DecimalException, OverflowError):
        raise SystemExit(f"{name}: molar mass must be finite and positive") from None
    if not math.isfinite(mass) or mass <= 0:
        raise SystemExit(f"{name}: molar mass must be finite and positive")
    ranges = [float(t) for t in thermo["temperature-ranges"]]
    coefficients = [[float(a) for a in row] for row in thermo["data"]]
    return {
        "molar_mass": mass,
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
    weights = {}
    for item in args.atomic_weights:
        element, value = item.split("=")
        if element in weights:
            raise SystemExit(f"atomic weight for {element}: given twice")
        try:
            weight = Decimal(value)
        except InvalidOperation:
            raise SystemExit(f"atomic weight for {element}: must be finite and positive") from None
        if not weight.is_finite() or weight <= 0:
            raise SystemExit(f"atomic weight for {element}: must be finite and positive")
        weights[element] = weight
    widths = species_json.parse_blend_widths(args.blend)
    with open(args.yaml_file, encoding="utf-8") as stream:
        # BaseLoader keeps every scalar a string: YAML 1.1 would read the
        # species name NO as false.
        data = yaml.load(stream, Loader=yaml.BaseLoader)
    entries = {}
    for entry in data["species"]:
        name = entry["name"]
        if name in entries:
            raise SystemExit(f"{args.yaml_file}: species {name} is given twice")
        entries[name] = entry
    species = {}
    for name in args.species:
        if name not in entries:
            raise SystemExit(f"{name}: not in {args.yaml_file}")
        try:
            species[name] = convert(entries[name], weights, widths)
        except ValueError as error:
            raise SystemExit(f"{name}: {error}") from None
    species_json.dump(species_json.species_file(species, reference=args.source), args.output)


if __name__ == "__main__":
    main()
