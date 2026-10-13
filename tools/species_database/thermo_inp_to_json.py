"""Converts species from a NASA Glenn thermo.inp file (CEA format, McBride,
Zehe and Gordon, NASA/TP-2002-211556, appendix A) to a canary species file.

The file gives no blending between temperature intervals, so the blend width
of each interior boundary is required: --blend 1000=400 6000=1000 ...

Usage: thermo_inp_to_json.py thermo.inp --species N2 O2 NO --blend 1000=400 6000=1000 --source TEXT -o out.json
"""

import argparse

import species_json


def find_species(lines, name):
    for i, line in enumerate(lines[:-1]):
        if line.startswith("!") or line[:18].split()[:1] != [name]:
            continue
        nxt = lines[i + 1]
        try:
            intervals = int(nxt[0:2])
            float(nxt[52:65])
        except ValueError:
            continue
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
            total -= float(count)
    return total


def convert(lines, name, widths):
    i, intervals = find_species(lines, name)
    record = lines[i + 1]
    ranges, coefficients = [], []
    for k in range(intervals):
        head = lines[i + 2 + 3 * k]
        first = lines[i + 3 + 3 * k]
        second = lines[i + 4 + 3 * k]
        low, high = species_json.number(head[0:11]), species_json.number(head[11:22])
        if not ranges:
            ranges.append(low)
        elif ranges[-1] != low:
            raise SystemExit(f"{name}: intervals are not contiguous at {low} K")
        ranges.append(high)
        a = [species_json.number(first[16 * j: 16 * (j + 1)]) for j in range(5)]
        a += [species_json.number(second[0:16]), species_json.number(second[16:32])]
        a += [species_json.number(second[48:64]), species_json.number(second[64:80])]
        coefficients.append(a)
    return {
        "molar_mass": species_json.kg_per_mol(record[52:65].strip()),
        "charge": charge(record),
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
    args = parser.parse_args()
    with open(args.thermo_inp, encoding="utf-8") as stream:
        lines = [line.rstrip("\n") for line in stream]
    widths = species_json.parse_blend_widths(args.blend)
    species = {name: convert(lines, name, widths) for name in args.species}
    species_json.dump(species_json.species_file(species, reference=args.source), args.output)


if __name__ == "__main__":
    main()
