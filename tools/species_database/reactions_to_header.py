"""Writes a canary reactions file (format canary-reactions-1) as a C++ header
of host input records (canary::MechanismInput, see include/canary/reactions.hpp)
and a builder of the mechanism for one thermo model's species.

Numbers are written with Python's repr, which gives the shortest decimal
that reads back to the same double, so the header holds exactly the values
of the file. Species lists keep the order of the file. The output is meant
to be formatted with clang-format afterwards.

Usage: reactions_to_header.py reactions.json --name air5_kim_jo
           --alias Air5KimJo --species-header canary/data/air5.hpp
           --species-array air5_species -o include/canary/data/air5_kim_jo.hpp
"""

import argparse
import json
from pathlib import Path

import reaction_json


def number(value):
    text = repr(float(value))
    if text in ("inf", "-inf", "nan"):
        raise SystemExit(f"number {text} has no C++ literal")
    return text


def string(text):
    # JSON string escapes are valid C++ escapes for ASCII text.
    if not text.isascii():
        raise SystemExit(f"non-ASCII text: {text!r}")
    return json.dumps(text)


def comment(text, width=76):
    """text as // comment lines of at most width characters."""
    lines, line = [], ""
    for word in text.split():
        if line and len(line) + 1 + len(word) > width:
            lines.append(line)
            line = word
        else:
            line = f"{line} {word}" if line else word
    if line:
        lines.append(line)
    return "\n".join(f"// {item}" for item in lines)


def coefficients(values):
    entries = "".join(f"{{.name = {string(name)}, .value = {number(v)},}}," for name, v in values.items())
    return f"{{.count = {len(values)}, .entries = {{{{{entries}}},}},}}"


def arrhenius(rate):
    """An ArrheniusInput; an equilibrium rate carries only its other_temperature."""
    other = rate.get("other_temperature") if rate else None
    other_text = f" .other_temperature = {string(other)}," if other is not None else ""
    if rate is None or rate["model"] == "equilibrium":
        if other is not None:
            return "{.A = 0.0, .n = 0.0, .activation_temperature = 0.0, .translational_exponent = 1.0," + other_text + "}"
        return "{.A = 0.0, .n = 0.0, .activation_temperature = 0.0,}"
    exponent = rate.get("translational_exponent")
    if exponent is None and other is not None:
        exponent = 1.0
    return (f"{{.A = {number(rate['A'])}, .n = {number(rate['n'])}, "
            f".activation_temperature = {number(rate['activation_temperature'])},"
            + (f" .translational_exponent = {number(exponent)}," if exponent is not None else "") + other_text + "}")


def reaction(record):
    backward = record["backward"]
    equilibrium = backward["model"] == "equilibrium"
    third_body = record.get("third_body")
    return "\n".join([
        "{",
        f".equation = {string(record['equation'])},",
        f".reactants = {coefficients(record['reactants'])},",
        f".products = {coefficients(record['products'])},",
        f".forward = {arrhenius(record['forward'])},",
        f".backward_model = BackwardRateModel::{'equilibrium' if equilibrium else 'arrhenius'},",
        f".backward = {arrhenius(backward)},",
        f".third_body = {'true' if third_body is not None else 'false'},",
        f".efficiencies = {coefficients(third_body['efficiencies'] if third_body else {})},",
        "},",
    ])


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("reactions")
    parser.add_argument("--name", required=True, help="name of the input records and of make_<name>")
    parser.add_argument("--alias", required=True, help="name of the mechanism type alias")
    parser.add_argument("--species-header", required=True, help="header of the thermo model's species")
    parser.add_argument("--species-array", required=True, help="array of canary::Nasa9Species in that header")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    data = reaction_json.load(args.reactions)
    reactions = data["reactions"]
    capacity = max(max(len(r["reactants"]), len(r["products"]),
                       len((r.get("third_body") or {}).get("efficiencies", {}))) for r in reactions)
    count = len(reactions)
    guard = f"CANARY_DATA_{args.name.upper()}_HPP"
    reference = data.get("reference", "")
    source = Path(args.reactions).name

    header = f"""#ifndef {guard}
#define {guard}

// Generated from database/data/{source} by
// tools/species_database/reactions_to_header.py; do not edit.
{comment(f"Reference: {reference}") if reference else ""}

#include <array>
#include <concepts>
#include <cstddef>
#include <canary/quantity.hpp>
#include <canary/reactions.hpp>
#include <canary/species.hpp>
#include <canary/units.hpp>
#include <string_view>
#include <{args.species_header}>

namespace canary
{{

// The {count} reactions of {source} in file order, in SI units per mole,
// with the temperature limits of rate evaluation in K.
inline constexpr MechanismInput<{count}, {capacity}> {args.name} {{
.temperature_lower = {number(data["temperature_limits"]["lower"])},
.temperature_upper = {number(data["temperature_limits"]["upper"])},
.reactions = {{
{{
{chr(10).join(reaction(r) for r in reactions)}
}},
}},
}};

// Mechanism of {args.name} over the species of {args.species_array}.
template <std::floating_point P, UnitSystem U>
using {args.alias} = ReactionMechanism<P, U, {args.species_array}.size(), {count}>;

namespace detail
{{
inline constexpr auto {args.name}_species_names = [] {{
    std::array<std::string_view, {args.species_array}.size()> names {{}};
    for (const SpeciesIndex i : species_indices(0, names.size())) {{ names[i.value()] = {args.species_array}[i.value()].name; }}
    return names;
}}();
}}  // namespace detail

// {args.alias} model of {args.name} for the species of {args.species_array}, in
// their order, with their species data, in the unit system of the scales,
// with standard pressure p0. See make_reaction_mechanism.
template <std::floating_point P, UnitSystem U>
{args.alias}<P, U> make_{args.name}(const SpeciesData<P, U, {args.species_array}.size()>& species,
                                const UnitScales<P, U>& scales,
                                Pressure<P, U> standard_pressure)
{{
    return make_reaction_mechanism({args.name}, detail::{args.name}_species_names, species, scales, standard_pressure);
}}

// As above, with the standard pressure 1 bar (1e5 Pa).
template <std::floating_point P, UnitSystem U>
{args.alias}<P, U> make_{args.name}(const SpeciesData<P, U, {args.species_array}.size()>& species,
                                const UnitScales<P, U>& scales)
{{
    return make_reaction_mechanism({args.name}, detail::{args.name}_species_names, species, scales);
}}

// As above for the thermo model gas with several energy modes: rates'
// other temperatures resolve to its modes. See make_reaction_mechanism.
template <std::floating_point P, UnitSystem U, typename G>
{args.alias}<P, U> make_{args.name}(const SpeciesData<P, U, {args.species_array}.size()>& species,
                                const UnitScales<P, U>& scales,
                                Pressure<P, U> standard_pressure,
                                const G& gas)
{{
    return make_reaction_mechanism({args.name}, detail::{args.name}_species_names, species, scales, standard_pressure, gas);
}}

}}  // namespace canary

#endif  // {guard}
"""
    Path(args.output).write_text(header, encoding="utf-8")


if __name__ == "__main__":
    main()
