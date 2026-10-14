"""Shared output structure of canary species files (format canary-species-1).

Numbers are written by Python's json module in the shortest form that reads
back to the same double. Decimal strings from source files are converted with
float(), and molar masses given in g/mol are divided by 1000 in decimal before
rounding, so that "28.0134" g/mol and "28.0134e-3" kg/mol give the same double.
"""

import json
from decimal import Decimal

FORMAT = "canary-species-1"
UNITS = {"molar_mass": "kg/mol", "energy": "J/mol", "temperature": "K", "length": "m"}


def kg_per_mol(grams_per_mol: str) -> float:
    return float(Decimal(grams_per_mol.replace("D", "E").replace("d", "e")) / Decimal(1000))


def number(text: str) -> float:
    return float(text.strip().replace("D", "E").replace("d", "e"))


def nasa9_thermo(ranges, coefficients, blend_widths, reference=None, notes=None):
    if len(ranges) < 2:
        raise ValueError("NASA9 requires at least two temperature boundaries")
    if len(coefficients) != len(ranges) - 1:
        raise ValueError("NASA9 coefficient rows must match the temperature intervals")
    if any(len(row) != 9 for row in coefficients):
        raise ValueError("NASA9 requires nine coefficients per temperature interval")
    if len(blend_widths) != len(ranges) - 2:
        raise ValueError("NASA9 blend widths must match the interior temperature boundaries")
    thermo = {
        "model": "nasa9",
        "temperature_ranges": ranges,
        "blend_widths": blend_widths,
        "coefficients": coefficients,
    }
    if reference:
        thermo["reference"] = reference
    if notes:
        thermo["notes"] = notes
    return thermo


def blend_widths_for(name, ranges, widths_by_boundary):
    """Blend widths for the interior boundaries of ranges, looked up by
    boundary temperature in widths_by_boundary ({temperature: width})."""
    widths = []
    for boundary in ranges[1:-1]:
        if boundary not in widths_by_boundary:
            raise SystemExit(f"{name}: no blend width given for the boundary at {boundary} K")
        widths.append(widths_by_boundary[boundary])
    return widths


def parse_blend_widths(items):
    """'1000=400' style arguments to {1000.0: 400.0}."""
    widths = {}
    for item in items:
        boundary, width = item.split("=")
        widths[float(boundary)] = float(width)
    return widths


def species_file(species, reference=None, notes=None):
    out = {"format": FORMAT, "units": UNITS, "species": species}
    if reference:
        out["reference"] = reference
    if notes:
        out["notes"] = notes
    return out


def dump(data, path):
    with open(path, "w", encoding="utf-8") as stream:
        json.dump(data, stream, indent=3)
        stream.write("\n")
