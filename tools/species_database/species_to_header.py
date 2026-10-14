"""Writes species of a canary species file (format canary-species-2) as a C++
header of host input records: with --thermo, the NASA-9 species
(canary::Nasa9Species, see include/canary/models/thermally_perfect.hpp); with
--electronic, the electronic levels (canary::ElectronicLevelsInput, see
include/canary/models/multi_temperature.hpp) of every species, zero levels for
a species without any.

Numbers are written with Python's repr, which gives the shortest decimal that
reads back to the same double. Species are written in the order given by
--species. The output is meant to be formatted with clang-format afterwards.

Usage: species_to_header.py species.json --species N2 O2 ... --name air11
           [--thermo --segments 4] [--electronic --levels 8] -o header.hpp
"""

import argparse
from pathlib import Path

import species_json

from reactions_to_header import comment, number, string

FORMAT = "canary-species-2"


def array(items):
    """A std::array initialiser with trailing commas: {{a, b,},}."""
    items = list(items)
    return "{{" + ", ".join(items) + ",},}" if items else "{}"


def padded(values, size, zero):
    values = list(values)
    if len(values) > size:
        raise SystemExit(f"{len(values)} entries do not fit {size}")
    return values + [zero] * (size - len(values))


def thermo_entry(name, record, segments, elements_name=None):
    role = record.get("role")
    if role not in ("heavy", "electron"):
        raise SystemExit(f"{name}: role must be heavy or electron")
    species_json.validate_elements(record, name)
    metadata = f".role = SpeciesRole::{role}"
    if elements_name is not None:
        metadata += f", .elements = std::span<const ElementCount>{{detail::{elements_name}}}"
    thermo = record["thermo"]
    if thermo["model"] != "nasa9":
        raise SystemExit(f"{name}: thermo model must be nasa9")
    count = len(thermo["coefficients"])
    breaks = padded(thermo["temperature_ranges"], segments + 1, 0.0)
    blends = padded(thermo["blend_widths"], segments, 0.0)
    rows = padded(thermo["coefficients"], segments, [0.0] * 9)
    return "\n".join([
        "{",
        f".name = {string(name)},",
        f".molar_mass = {number(record['molar_mass'])},",
        f".charge = {number(record['charge'])},",
        f".metadata = {{{metadata}}},",
        f".fit = {{.segments = {count}, .breaks = {array(number(x) for x in breaks)}, "
        f".blends = {array(number(x) for x in blends)}, "
        f".a = {{{{{', '.join(array(number(x) for x in row) for row in rows)},}},}},}},",
        "},",
    ])


def levels_entry(name, record, levels):
    data = record.get("electronic_levels") or {"degeneracies": [], "temperatures": []}
    g, theta = data["degeneracies"], data["temperatures"]
    return "\n".join([
        "{",
        f".name = {string(name)},",
        f".levels = {len(g)},",
        f".degeneracies = {array(str(int(x)) for x in padded(g, levels, 0))},",
        f".temperatures = {array(number(x) for x in padded(theta, levels, 0.0))},",
        "},",
    ])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("species_file")
    parser.add_argument("--species", nargs="+", required=True)
    parser.add_argument("--name", required=True, help="prefix of the generated arrays")
    parser.add_argument("--thermo", action="store_true")
    parser.add_argument("--segments", type=int, default=4)
    parser.add_argument("--electronic", action="store_true")
    parser.add_argument("--levels", type=int, default=8)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    data = species_json.load(args.species_file)
    if data["format"] != FORMAT:
        raise SystemExit(f"{args.species_file}: not a {FORMAT} file")
    species = data["species"]
    source = Path(args.species_file).name
    guard = f"CANARY_DATA_{Path(args.output).stem.upper()}_HPP"
    includes = ["#include <array>"]
    blocks = []
    count = len(args.species)
    if args.thermo:
        includes.append("#include <canary/models/thermally_perfect.hpp>")
        includes.append("#include <span>")
        element_names = {}
        element_arrays = []
        for i, name in enumerate(args.species):
            record = species[name]
            species_json.validate_elements(record, name)
            if record.get("elements") is not None:
                identifier = f"{args.name}_elements_{i}"
                element_names[name] = identifier
                elements = species_json.element_counts(record["elements"].items(), name)
                values = [f"{{{string(e)}, {n}}}" for e, n in elements.items()]
                element_arrays.append(f"inline constexpr std::array<ElementCount, {len(values)}> {identifier} {array(values)};")
        if element_arrays:
            blocks.append("namespace detail\n{\n" + "\n".join(element_arrays) + "\n}  // namespace detail")
        blocks.append(f"""// NASA-9 data of {", ".join(args.species)}, in this order: molar masses in
// kg/mol, charges in elementary charges, fits in SI (kelvin); entries past a
// species' segment count are zero.
inline constexpr std::array<Nasa9Species<{args.segments}>, {count}> {args.name}_species {{{{
{chr(10).join(thermo_entry(n, species[n], args.segments, element_names.get(n)) for n in args.species)}
}},}};""")
    if args.electronic:
        includes.append("#include <canary/models/multi_temperature.hpp>")
        blocks.append(f"""// Electronic levels of {", ".join(args.species)}, in this order: degeneracies
// and level temperatures E_k / k_B in K; entries past a species' level count
// are zero.
inline constexpr std::array<ElectronicLevelsInput<{args.levels}>, {count}> {args.name}_electronic_levels {{{{
{chr(10).join(levels_entry(n, species[n], args.levels) for n in args.species)}
}},}};""")
    reference = data.get("reference", "")
    header = f"""#ifndef {guard}
#define {guard}

// Generated from database/data/{source} by
// tools/species_database/species_to_header.py; do not edit.
{comment(f"Reference: {reference}") if reference else ""}

{chr(10).join(includes)}

namespace canary
{{

{(chr(10) * 2).join(blocks)}

}}  // namespace canary

#endif  // {guard}
"""
    Path(args.output).write_text(header, encoding="utf-8")


if __name__ == "__main__":
    main()
