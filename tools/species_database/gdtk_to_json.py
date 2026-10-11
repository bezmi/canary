"""Converts species from GDTk species-database Lua files to a gaslib species file.

Reads db.<name>.M, db.<name>.charge and db.<name>.ceaThermoCoeffs (NASA-9
segments, break points and blend ranges) from <database>/<name>.lua. With
--transport, also db.<name>.ceaViscosity and db.<name>.ceaThermCond (CEA
transport fits ln(x) = A ln T + B / T + C / T^2 + D with the viscosity in
micropoise and the conductivity in microwatt / (cm K); D is converted to SI by
adding ln(1e-7) and ln(1e-4)) and the Lennard-Jones db.<name>.sigma
(Angstrom, converted to m) and db.<name>.epsilon (K).

Usage: gdtk_to_json.py --database DIR --species N2 O2 N O NO [--transport]
           --source TEXT -o out.json
"""

import argparse
import math
import re
from decimal import Decimal
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


def cea_curve(text, name, key, unit_scale):
    """db.<name>.<key> CEA transport fit as a gaslib transport curve in SI."""
    body = table(text, name, key)
    segments = int(re.search(r"nsegments\s*=\s*(\d+)", body).group(1))
    ranges, coefficients = [], []
    for k in range(segments):
        segment = field(body, f"segment{k}")
        values = dict(re.findall(r"\b(T_lower|T_upper|A|B|C|D)\s*=\s*([0-9.eE+-]+)", segment))
        lower, upper = species_json.number(values["T_lower"]), species_json.number(values["T_upper"])
        if ranges and ranges[-1] != lower:
            raise SystemExit(f"{name}: {key} segments are not contiguous")
        if not ranges:
            ranges.append(lower)
        ranges.append(upper)
        coefficients.append([species_json.number(values["A"]), species_json.number(values["B"]),
                             species_json.number(values["C"]),
                             species_json.number(values["D"]) + math.log(unit_scale)])
    return {"model": "cea", "temperature_ranges": ranges, "coefficients": coefficients}


def value_and_reference(text, name, key):
    body = table(text, name, key)
    value = re.search(r"value\s*=\s*([0-9.eE+-]+)", body).group(1)
    return value, quoted(body, "reference")


def transport(text, name):
    sigma, sigma_reference = value_and_reference(text, name, "sigma")
    epsilon, _ = value_and_reference(text, name, "epsilon")
    lennard_jones = {"sigma": float(Decimal(sigma) * Decimal("1e-10")), "epsilon": species_json.number(epsilon)}
    if sigma_reference:
        lennard_jones["reference"] = sigma_reference
    return {
        "viscosity": cea_curve(text, name, "ceaViscosity", 1.0e-7),
        "thermal_conductivity": cea_curve(text, name, "ceaThermCond", 1.0e-4),
        "lennard_jones": lennard_jones,
        "notes": "CEA transport fits with D converted to SI (Pa s, W/(m K)) by adding ln(1e-7) and ln(1e-4)",
    }


def convert(path, name, with_transport=False):
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
    if with_transport:
        record["transport"] = transport(text, name)
    mass_reference = quoted(mass, "reference")
    if mass_reference:
        record["notes"] = f"molar mass: {mass_reference}"
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--database", required=True, type=Path)
    parser.add_argument("--species", required=True, nargs="+")
    parser.add_argument("--transport", action="store_true", help="also convert the transport data")
    parser.add_argument("--source", required=True, help="where the database files came from (version, path)")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    species = {name: convert(args.database / f"{name}.lua", name, args.transport) for name in args.species}
    species_json.dump(species_json.species_file(species, reference=args.source), args.output)


if __name__ == "__main__":
    main()
