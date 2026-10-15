"""Writes a canary energy-exchange file (format canary-energy-exchange-1) as a
C++ header of host input records (canary::EnergyExchangeInput, see
include/canary/energy_exchange.hpp) and a builder of the energy exchange for
one thermo model's species and reaction mechanism.

Numbers are written with Python's repr, which gives the shortest decimal
that reads back to the same double. Entries keep the order of the file. The
output is meant to be formatted with clang-format afterwards.

Usage: exchange_to_header.py exchange.json --name air5_park_exchange
           --alias Air5ParkExchange --species-header canary/data/air5.hpp
           --species-array air5_species --mechanism-alias Air5Park
           --mechanism-header canary/data/air5_park.hpp
           -o include/canary/data/air5_park_exchange.hpp
"""

import argparse
from pathlib import Path

import exchange_json
import species_json
from reactions_to_header import comment, number, string


def _optional_number(record, field):
    value = record.get(field)
    return number(value if value is not None else 0.0)


def vibration_translation(entry):
    time = entry["relaxation_time"]
    correction = time.get("high_temperature_correction") or {"model": "none", "cross_section": 0.0}
    correction_text = (f"{{.model = RelaxationCorrection::{correction['model']}, "
                       f".cross_section = {number(correction['cross_section'])}, "
                       f".reference_temperature = {_optional_number(correction, 'reference_temperature')}, "
                       f".exponent = {_optional_number(correction, 'exponent')},}}")
    return "\n".join([
        "{",
        f".species = {string(entry['species'])},",
        f".partner = {string(entry['partner'])},",
        f".millikan_white = {{.a = {number(time['a'])}, .b = {number(time['b'])},}},",
        f".correction = {correction_text},",
        "},",
    ])


def chemistry_vibration(entry):
    return "\n".join([
        "{",
        f".reaction = {int(entry['reaction'])},",
        f".species = {string(entry['species'])},",
        f".dissociation_energy = {number(entry['dissociation_energy'])},",
        f".vibrational_temperature = {number(entry['vibrational_temperature'])},",
        "},",
    ])


def electron_translation(entry):
    section = entry["cross_section"]
    return (f"{{.partner = {string(entry['partner'])}, "
            f".cross_section = ElectronCrossSection::{section['model']}, "
            + "".join(f".{k} = {_optional_number(section, k)}, " for k in ("a", "b", "c", "maximum_temperature"))
            + "},")


def vibration_electron(entry):
    time = entry["relaxation_time"]
    keys = ("a_low", "b_low", "c_low", "a_high", "b_high", "c_high", "blend_lower", "blend_upper")
    fit = "{" + "".join(f".{k} = {number(time[k])}, " for k in keys) + "}"
    return f"{{.species = {string(entry['species'])}, .relaxation_time = {fit},}},"


def ionisation(entry, species):
    record = species.get(entry["species"], {})
    if "ionisation_energy" not in record:
        raise SystemExit(f"{entry['species']}: no ionisation_energy in the species file")
    return (f"{{.reaction = {int(entry['reaction'])}, .species = {string(entry['species'])}, "
            f".ionisation_energy = {number(record['ionisation_energy'])},}},")


def listed(items):
    return "{{\n" + "\n".join(items) + "\n},}" if items else "{}"


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("exchange")
    parser.add_argument("--name", required=True, help="name of the input records and of make_<name>")
    parser.add_argument("--alias", required=True, help="name of the energy-exchange type alias")
    parser.add_argument("--species-header", required=True, help="header of the thermo model's species")
    parser.add_argument("--species-array", required=True, help="array of canary::Nasa9Species in that header")
    parser.add_argument("--mechanism-header", required=True, help="header of the reaction mechanism")
    parser.add_argument("--mechanism-alias", required=True, help="mechanism type alias in that header")
    parser.add_argument("--species-file", help="species file with ionisation energies (needed for C-E entries)")
    parser.add_argument("--model-alias", help="multi-temperature model type alias, for exchanges with electrons")
    parser.add_argument("--model-header", help="header of that model")
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    data = exchange_json.load(args.exchange)
    vt = data["vibration_translation"]
    cv = data["chemistry_vibration"]
    et = data.get("electron_translation", [])
    ve = data.get("vibration_electron", [])
    ce = data.get("electron_impact_ionisation", [])
    electrons = bool(et or ve or ce)
    if electrons and not (args.model_alias and args.model_header):
        raise SystemExit("an exchange with electrons needs --model-alias and --model-header")
    species = {}
    if args.species_file:
        species_data = species_json.load_validated(args.species_file)
        species = species_data["species"]
    guard = f"CANARY_DATA_{args.name.upper()}_HPP"
    reference = data.get("reference", "")
    source = Path(args.exchange).name

    includes = ["#include <array>", "#include <concepts>", "#include <cstddef>",
                "#include <canary/energy_exchange.hpp>", "#include <canary/species.hpp>",
                "#include <canary/units.hpp>", "#include <string_view>",
                f"#include <{args.species_header}>", f"#include <{args.mechanism_header}>"]
    if electrons:
        includes.append(f"#include <{args.model_header}>")
    names = f"""namespace detail
{{
inline constexpr auto {args.name}_species_names = [] {{
    std::array<std::string_view, {args.species_array}.size()> names {{}};
    for (const SpeciesIndex i : species_indices(0, names.size())) {{ names[i.value()] = {args.species_array}[i.value()].name; }}
    return names;
}}();
}}  // namespace detail"""
    if electrons:
        tail = f"""// The electron exchanges of {source}: {len(et)} elastic
// electron-heavy pairs, {len(ve)} vibration-electron exchanges and {len(ce)}
// electron-impact ionisations (ionisation energies in J/mol from the
// species file).
inline constexpr ElectronExchangeInput<{len(et)}, {len(ve)}, {len(ce)}> {args.name}_electrons {{
.electron_translation = {listed([electron_translation(e) for e in et])},
.vibration_electron = {listed([vibration_electron(e) for e in ve])},
.electron_impact_ionisation = {listed([ionisation(e, species) for e in ce])},
}};

// Energy exchange of {args.name} over the species and modes of {args.model_alias}.
template <std::floating_point P, UnitSystem U>
using {args.alias} = MultiTemperatureExchange<P, U, {args.species_array}.size(), {args.model_alias}<P, U>::mode_count,
  {len(vt)}, {len(cv)}, {len(et)}, {len(ve)}, {len(ce)}>;

{names}

// {args.alias} model of {args.name} for the species of {args.species_array}, in
// their order, with their species data, the reactions of mechanism and the
// modes of gas, in the unit system of the scales. See
// make_multi_temperature_exchange.
template <std::floating_point P, UnitSystem U>
{args.alias}<P, U> make_{args.name}(const SpeciesData<P, U, {args.species_array}.size()>& species,
                                const {args.mechanism_alias}<P, U>& mechanism,
                                const {args.model_alias}<P, U>& gas,
                                const UnitScales<P, U>& scales)
{{
    return make_multi_temperature_exchange({args.name}, {args.name}_electrons, detail::{args.name}_species_names,
                                           species, mechanism, gas, scales);
}}"""
    else:
        tail = f"""// Energy exchange of {args.name} over the species of {args.species_array}.
template <std::floating_point P, UnitSystem U>
using {args.alias} = EnergyExchange<P, U, {args.species_array}.size(), {len(vt)}, {len(cv)}>;

{names}

// {args.alias} model of {args.name} for the species of {args.species_array}, in
// their order, with their species data, and the reactions of mechanism, in the
// unit system of the scales. See make_energy_exchange.
template <std::floating_point P, UnitSystem U>
{args.alias}<P, U> make_{args.name}(const SpeciesData<P, U, {args.species_array}.size()>& species,
                                const {args.mechanism_alias}<P, U>& mechanism,
                                const UnitScales<P, U>& scales)
{{
    return make_energy_exchange({args.name}, detail::{args.name}_species_names, species, mechanism, scales);
}}"""
    header = f"""#ifndef {guard}
#define {guard}

// Generated from database/data/{source} by
// tools/species_database/exchange_to_header.py; do not edit.
{comment(f"Reference: {reference}") if reference else ""}

{chr(10).join(includes)}

namespace canary
{{

// The {len(vt)} vibration-translation pairs and {len(cv)} chemistry-vibration
// couplings of {source} in file order, in SI units per mole.
inline constexpr EnergyExchangeInput<{len(vt)}, {len(cv)}> {args.name} {{
.vibration_translation = {listed([vibration_translation(e) for e in vt])},
.chemistry_vibration = {listed([chemistry_vibration(e) for e in cv])},
}};

{tail}

}}  // namespace canary

#endif  // {guard}
"""
    Path(args.output).write_text(header, encoding="utf-8")


if __name__ == "__main__":
    main()
