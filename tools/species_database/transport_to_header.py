"""Writes the transport data of a canary species file (format canary-species-1,
the "transport" record of each species) and a collision-integral file (format
canary-collision-integrals-1) as a C++ header of host input records
(canary::SpeciesTransportInput and canary::CollisionPairInput, see
include/canary/transport/mixture.hpp).

Numbers are written with Python's repr, which gives the shortest decimal
that reads back to the same double. Species are written in the order given
by --species; collision pairs in file order. The output is meant to be
formatted with clang-format afterwards.

Usage: transport_to_header.py species.json collision_integrals.json
           --species N2 O2 N O NO --name air5 --segments 3
           -o include/canary/data/air5_transport.hpp
"""

import argparse
import json
from pathlib import Path

from reactions_to_header import comment, number, string

SPECIES_FORMAT = "canary-species-1"
COLLISION_FORMAT = "canary-collision-integrals-1"


def array(items):
    """A std::array initialiser with trailing commas: {{a, b,},}."""
    return "{{" + ", ".join(items) + ",},}"


def curve(name, record, segments):
    if record["model"] != "cea":
        raise SystemExit(f"{name}: unsupported transport curve model {record['model']}")
    ranges, coefficients = record["temperature_ranges"], record["coefficients"]
    if len(coefficients) > segments or len(ranges) != len(coefficients) + 1:
        raise SystemExit(f"{name}: {len(coefficients)} segments do not fit {segments}")
    breaks = ranges + [0.0] * (segments + 1 - len(ranges))
    rows = coefficients + [[0.0] * 4] * (segments - len(coefficients))
    row_text = ", ".join(array(number(x) for x in row) for row in rows)
    return (f"{{.segments = {len(coefficients)}, "
            f".breaks = {array(number(x) for x in breaks)}, "
            f".coefficients = {{{{{row_text},}},}},}}")


def species_entry(name, record, segments):
    transport = record.get("transport")
    if not transport or not {"viscosity", "thermal_conductivity", "lennard_jones"} <= transport.keys():
        raise SystemExit(f"{name}: needs viscosity, thermal_conductivity and lennard_jones transport data")
    lj = transport["lennard_jones"]
    return "\n".join([
        "{",
        f".name = {string(name)},",
        f".viscosity = {curve(name, transport['viscosity'], segments)},",
        f".conductivity = {curve(name, transport['thermal_conductivity'], segments)},",
        f".lennard_jones = {{.sigma = {number(lj['sigma'])}, .epsilon = {number(lj['epsilon'])},}},",
        "},",
    ])


def fit(record):
    return "{" + ", ".join(f".{k} = {number(record[k])}" for k in "ABCD") + ",}"


def pair_entry(record):
    a, b = record["species"]
    return (f"{{.first = {string(a)}, .second = {string(b)}, "
            f".pi_omega_11 = {fit(record['pi_omega_11'])}, .pi_omega_22 = {fit(record['pi_omega_22'])},}},")


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("species_file")
    parser.add_argument("collision_file")
    parser.add_argument("--species", nargs="+", required=True)
    parser.add_argument("--name", required=True, help="prefix of the generated arrays")
    parser.add_argument("--segments", type=int, required=True, help="segment capacity of the curves")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    species = json.loads(Path(args.species_file).read_text(encoding="utf-8"))
    collisions = json.loads(Path(args.collision_file).read_text(encoding="utf-8"))
    if species["format"] != SPECIES_FORMAT or collisions["format"] != COLLISION_FORMAT:
        raise SystemExit("not a species file and a collision-integral file")
    entries = [species_entry(n, species["species"][n], args.segments) for n in args.species]
    wanted = set(args.species)
    pairs = [p for p in collisions["pairs"] if set(p["species"]) <= wanted]
    guard = f"CANARY_DATA_{args.name.upper()}_TRANSPORT_HPP"
    species_source = Path(args.species_file).name
    collision_source = Path(args.collision_file).name
    references = "; ".join(r for r in (species.get("reference"), collisions.get("reference")) if r)

    header = f"""#ifndef {guard}
#define {guard}

// Generated from database/data/{species_source} and
// database/data/{collision_source} by
// tools/species_database/transport_to_header.py; do not edit.
{comment(f"References: {references}") if references else ""}

#include <array>
#include <canary/transport/mixture.hpp>

namespace canary
{{

// Transport data of {", ".join(args.species)}, in this order, in SI units:
// CEA viscosity and conductivity fits (ln(x / SI unit) = A ln T + B / T +
// C / T^2 + D) and Lennard-Jones sigma in m and epsilon / k_B in K.
inline constexpr std::array<SpeciesTransportInput<{args.segments}>, {len(entries)}> {args.name}_species_transport {{{{
{chr(10).join(entries)}
}},}};

// Collision integrals of the {len(pairs)} pairs of these species, in the order of
// {collision_source}: ln(pi Omega / m^2) = A (ln T)^3 + B (ln T)^2 + C ln T + D.
inline constexpr std::array<CollisionPairInput, {len(pairs)}> {args.name}_collision_integrals {{{{
{chr(10).join(pair_entry(p) for p in pairs)}
}},}};

}}  // namespace canary

#endif  // {guard}
"""
    Path(args.output).write_text(header, encoding="utf-8")


if __name__ == "__main__":
    main()
