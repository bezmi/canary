"""Converts species from GDTk species-database Lua files to a gaslib species file.

Reads db.<name>.M, db.<name>.charge and db.<name>.ceaThermoCoeffs (NASA-9
segments, break points and blend ranges) from <database>/<name>.lua.

Usage: gdtk_to_json.py --database DIR --species N2 O2 N O NO --source TEXT -o out.json
"""

import argparse
import re
from pathlib import Path

import species_json


def table(text, name, key):
    """Body of the Lua table db.<name>.<key> = { ... } (balanced braces)."""
    match = re.search(rf"db\.{re.escape(name)}\.{key}\s*=\s*\{{", text)
    if not match:
        raise SystemExit(f"{name}: no {key} table")
    depth, start = 1, match.end()
    for i in range(start, len(text)):
        depth += {"{": 1, "}": -1}.get(text[i], 0)
        if depth == 0:
            return text[start:i]
    raise SystemExit(f"{name}: unterminated {key} table")


def numbers(body):
    return [species_json.number(x) for x in body.split(",") if x.strip()]


def field(body, key):
    match = re.search(rf"\b{key}\s*=\s*\{{([^}}]*)\}}", body)
    if not match:
        raise SystemExit(f"no {key} in table")
    return match.group(1)


def quoted(body, key):
    match = re.search(rf"\b{key}\s*=\s*['\"]([^'\"]*)['\"]", body)
    return match.group(1) if match else None


def convert(path, name):
    text = path.read_text(encoding="utf-8")
    mass = table(text, name, "M")
    value = re.search(r"value\s*=\s*([0-9.eE+-]+)", mass).group(1)
    charge = re.search(rf"db\.{re.escape(name)}\.charge\s*=\s*([0-9.eE+-]+)", text)
    thermo = table(text, name, "ceaThermoCoeffs")
    segments = int(re.search(r"nsegments\s*=\s*(\d+)", thermo).group(1))
    ranges = numbers(field(thermo, "T_break_points"))
    blends = numbers(field(thermo, "T_blend_ranges")) if segments > 1 else []
    coefficients = [numbers(field(thermo, f"segment{k}")) for k in range(segments)]
    if len(ranges) != segments + 1 or any(len(c) != 9 for c in coefficients):
        raise SystemExit(f"{name}: inconsistent ceaThermoCoeffs")
    record = {
        "molar_mass": species_json.number(value),
        "charge": species_json.number(charge.group(1)) if charge else 0.0,
        "thermo": species_json.nasa9_thermo(ranges, coefficients, blends, quoted(thermo, "ref") or quoted(thermo, "reference")),
    }
    mass_reference = quoted(mass, "reference")
    if mass_reference:
        record["notes"] = f"molar mass: {mass_reference}"
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--database", required=True, type=Path)
    parser.add_argument("--species", required=True, nargs="+")
    parser.add_argument("--source", required=True, help="where the database files came from (version, path)")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    species = {name: convert(args.database / f"{name}.lua", name) for name in args.species}
    species_json.dump(species_json.species_file(species, reference=args.source), args.output)


if __name__ == "__main__":
    main()
