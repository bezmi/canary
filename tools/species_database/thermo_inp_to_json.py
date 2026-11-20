# SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
# Copyright (C) 2026 Shahzeb Imran and the canary contributors
"""Converts species from a NASA Glenn thermo.inp file (CEA format, McBride,
Zehe and Gordon, NASA/TP-2002-211556, appendix C) to a canary species file.

The file gives no blending between temperature intervals, so the blend width
of each interior boundary is required: --blend 1000=400 6000=1000 ...

Every selected species requires an explicit --roles NAME=heavy|electron declaration.
Roles are independent of species names, ordering and charge; electron roles require charge -1.

Usage: thermo_inp_to_json.py thermo.inp --species N2 O2 NO --roles N2=heavy O2=heavy NO=heavy --blend 1000=400 6000=1000 --source TEXT -o out.json
"""

import argparse
from decimal import InvalidOperation

import species_json


def find_species(lines, name):
    for i, line in enumerate(lines):
        if line.startswith("!") or line[:18].split()[:1] != [name]:
            continue
        if i + 1 >= len(lines) or len(lines[i + 1]) < 65:
            raise SystemExit(f"{name}: incomplete species metadata record")
        try:
            intervals = int(lines[i + 1][0:2])
        except ValueError:
            raise SystemExit(f"{name}: invalid temperature interval count") from None
        if intervals <= 0:
            raise SystemExit(f"{name}: temperature interval count must be positive")
        if len(lines) < i + 2 + 3 * intervals:
            raise SystemExit(f"{name}: incomplete NASA9 records for {intervals} temperature intervals")
        return i, intervals
    raise SystemExit(f"{name}: not found in thermo.inp")


def charge(record):
    """Charge from the formula field: 5 pairs of (element A2, count F6.2); an
    electron count E n gives charge -n."""
    total = 0.0
    for k in range(5):
        element = record[10 + 8 * k: 12 + 8 * k].strip()
        count = record[12 + 8 * k: 18 + 8 * k].strip()
        if element.upper() == "E" and count:
            total -= species_json.number(count)
    return total


def elements(record, name):
    """Nuclear entries from the formula field; E is charge metadata."""
    entries = []
    for k in range(5):
        element = record[10 + 8 * k:12 + 8 * k].strip()
        count = record[12 + 8 * k:18 + 8 * k].strip()
        if element and element.upper() != "E":
            entries.append((element, species_json.number(count)))
    return species_json.element_counts(entries, name)


def convert(lines, name, widths):
    i, intervals = find_species(lines, name)
    record = lines[i + 1]
    ranges, coefficients = [], []
    for k in range(intervals):
        head = lines[i + 2 + 3 * k]
        first = lines[i + 3 + 3 * k]
        second = lines[i + 4 + 3 * k]
        if len(head) < 58 or len(first) < 80 or len(second) < 80:
            raise SystemExit(f"{name}: interval {k + 1}: incomplete NASA9 record")
        try:
            count = int(head[22:23])
            powers = [species_json.number(head[23 + 5 * j:28 + 5 * j]) for j in range(7)]
            low, high = species_json.number(head[0:11]), species_json.number(head[11:22])
            a = [species_json.number(first[16 * j:16 * (j + 1)]) for j in range(5)]
            a += [species_json.number(second[0:16]), species_json.number(second[16:32])]
            a += [species_json.number(second[48:64]), species_json.number(second[64:80])]
        except ValueError:
            raise SystemExit(f"{name}: interval {k + 1}: invalid NASA9 numeric field") from None
        if count != 7 or powers != [-2, -1, 0, 1, 2, 3, 4]:
            raise SystemExit(f"{name}: interval {k + 1}: NASA9 requires seven powers -2, -1, 0, 1, 2, 3, 4")
        if not ranges:
            ranges.append(low)
        elif ranges[-1] != low:
            raise SystemExit(f"{name}: intervals are not contiguous at {low} K")
        ranges.append(high)
        coefficients.append(a)
    return {
        "molar_mass": species_json.kg_per_mol(record[52:65].strip()),
        "charge": charge(record),
        "elements": elements(record, name),
        "thermo": species_json.nasa9_thermo(
            ranges, coefficients, species_json.blend_widths_for(name, ranges, widths), lines[i][18:].strip() or None
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("thermo_inp")
    parser.add_argument("--species", required=True, nargs="+")
    parser.add_argument("--blend", required=True, nargs="+", help="boundary=width in K")
    parser.add_argument("--source", required=True)
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--roles", required=True, nargs="+", metavar="NAME=ROLE",
                        help="explicit heavy or electron role for every selected species")
    args = parser.parse_args()
    roles = species_json.parse_roles(args.roles, args.species)
    with open(args.thermo_inp, encoding="utf-8") as stream:
        lines = [line.rstrip("\n") for line in stream]
    widths = species_json.parse_blend_widths(args.blend)
    species = {}
    for name in args.species:
        try:
            species[name] = convert(lines, name, widths)
        except (ValueError, InvalidOperation) as error:
            raise SystemExit(f"{name}: invalid NASA9 record: {error}") from None
    for name, record in species.items():
        record["role"] = roles[name]
    species_json.dump(species_json.species_file(species, reference=args.source), args.output)


if __name__ == "__main__":
    main()
