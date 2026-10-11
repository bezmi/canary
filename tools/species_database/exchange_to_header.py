"""Writes a gaslib energy-exchange file (format gaslib-energy-exchange-1) as a
C++ header of host input records (gaslib::EnergyExchangeInput, see
include/gaslib/energy_exchange.hpp) and a builder of the energy exchange for
one thermo model's species and reaction mechanism.

Numbers are written with Python's repr, which gives the shortest decimal
that reads back to the same double. Entries keep the order of the file. The
output is meant to be formatted with clang-format afterwards.

Usage: exchange_to_header.py exchange.json --name air5_park_exchange
           --alias Air5ParkExchange --species-header gaslib/data/air5.hpp
           --species-array air5_species --mechanism-alias Air5Park
           --mechanism-header gaslib/data/air5_park.hpp
           -o include/gaslib/data/air5_park_exchange.hpp
"""

import argparse
import json
from pathlib import Path

from reactions_to_header import comment, number, string

FORMAT = "gaslib-energy-exchange-1"
UNITS = {"energy": "J/mol", "temperature": "K", "area": "m^2"}


def vibration_translation(entry):
    if entry["model"] != "landau_teller" or entry["relaxation_time"]["model"] != "millikan_white":
        raise SystemExit(f"{entry['species']}-{entry['partner']}: unsupported V-T model")
    time = entry["relaxation_time"]
    correction = time.get("high_temperature_correction")
    if correction is not None and correction["model"] != "park_1993":
        raise SystemExit(f"{entry['species']}-{entry['partner']}: unsupported correction")
    park = (f"{{.cross_section = {number(correction['cross_section'])}, "
            f".reference_temperature = {number(correction['reference_temperature'])},}}"
            if correction else "{.cross_section = 0.0, .reference_temperature = 0.0,}")
    return "\n".join([
        "{",
        f".species = {string(entry['species'])},",
        f".partner = {string(entry['partner'])},",
        f".millikan_white = {{.a = {number(time['a'])}, .b = {number(time['b'])},}},",
        f".park_correction = {'true' if correction else 'false'},",
        f".park = {park},",
        "},",
    ])


def chemistry_vibration(entry):
    if entry["model"] != "impartial_dissociation":
        raise SystemExit(f"reaction {entry['reaction']}: unsupported C-V model")
    return "\n".join([
        "{",
        f".reaction = {int(entry['reaction'])},",
        f".species = {string(entry['species'])},",
        f".dissociation_energy = {number(entry['dissociation_energy'])},",
        f".vibrational_temperature = {number(entry['vibrational_temperature'])},",
        "},",
    ])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("exchange")
    parser.add_argument("--name", required=True, help="name of the input records and of make_<name>")
    parser.add_argument("--alias", required=True, help="name of the energy-exchange type alias")
    parser.add_argument("--species-header", required=True, help="header of the thermo model's species")
    parser.add_argument("--species-array", required=True, help="array of gaslib::Nasa9Species in that header")
    parser.add_argument("--mechanism-header", required=True, help="header of the reaction mechanism")
    parser.add_argument("--mechanism-alias", required=True, help="mechanism type alias in that header")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    data = json.loads(Path(args.exchange).read_text(encoding="utf-8"))
    if data["format"] != FORMAT or data["units"] != UNITS:
        raise SystemExit(f"{args.exchange}: not a {FORMAT} file in SI units")
    vt = data["vibration_translation"]
    cv = data["chemistry_vibration"]
    guard = f"GASLIB_DATA_{args.name.upper()}_HPP"
    reference = data.get("reference", "")
    source = Path(args.exchange).name

    header = f"""#ifndef {guard}
#define {guard}

// Generated from database/data/{source} by
// tools/species_database/exchange_to_header.py; do not edit.
{comment(f"Reference: {reference}") if reference else ""}

#include <array>
#include <concepts>
#include <cstddef>
#include <gaslib/energy_exchange.hpp>
#include <gaslib/species.hpp>
#include <gaslib/units.hpp>
#include <string_view>
#include <{args.species_header}>
#include <{args.mechanism_header}>

namespace gaslib
{{

// The {len(vt)} vibration-translation pairs and {len(cv)} chemistry-vibration
// couplings of {source} in file order, in SI units per mole.
inline constexpr EnergyExchangeInput<{len(vt)}, {len(cv)}> {args.name} {{
.vibration_translation = {{
{{
{chr(10).join(vibration_translation(e) for e in vt)}
}},
}},
.chemistry_vibration = {{
{{
{chr(10).join(chemistry_vibration(e) for e in cv)}
}},
}},
}};

// Energy exchange of {args.name} over the species of {args.species_array}.
template <std::floating_point P, UnitSystem U>
using {args.alias} = EnergyExchange<P, U, {args.species_array}.size(), {len(vt)}, {len(cv)}>;

namespace detail
{{
inline constexpr auto {args.name}_species_names = [] {{
    std::array<std::string_view, {args.species_array}.size()> names {{}};
    for (std::size_t i = 0; i < names.size(); ++i) {{ names[i] = {args.species_array}[i].name; }}
    return names;
}}();
}}  // namespace detail

// {args.alias} model of {args.name} for the species of {args.species_array}, in
// their order, with their species data, and the reactions of mechanism, in the
// unit system of the scales. See make_energy_exchange.
template <std::floating_point P, UnitSystem U>
{args.alias}<P, U> make_{args.name}(const SpeciesData<P, U, {args.species_array}.size()>& species,
                                const {args.mechanism_alias}<P, U>& mechanism,
                                const UnitScales<P, U>& scales)
{{
    return make_energy_exchange({args.name}, detail::{args.name}_species_names, species, mechanism, scales);
}}

}}  // namespace gaslib

#endif  // {guard}
"""
    Path(args.output).write_text(header, encoding="utf-8")


if __name__ == "__main__":
    main()
