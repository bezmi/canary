"""Converts species from GDTk species-database Lua files to a canary species file.

Reads db.<name>.M, db.<name>.charge and db.<name>.ceaThermoCoeffs (NASA-9
segments, break points and blend ranges) from <database>/<name>.lua. With
--transport, also db.<name>.ceaViscosity and db.<name>.ceaThermCond (CEA
transport fits ln(x) = A ln T + B / T + C / T^2 + D with the viscosity in
micropoise and the conductivity in microwatt / (cm K); D is converted to SI by
adding ln(1e-7) and ln(1e-4)) and the Lennard-Jones db.<name>.sigma
(Angstrom, converted to m) and db.<name>.epsilon (K). With --electronic,
also db.<name>.electronic_levels (degeneracies g and level energies Te in
cm^-1, converted to temperatures in K with h c / k_B). Entries are written
db.<name> or db['<name>']; the file of species N2+ is N2_plus.lua and that
of e- is e_minus.lua.

With --ionisation-energy N=14.53413 ..., the given ionisation energies
(in eV, converted to J/mol with e N_A) are added; GDTk's database has none.

Every selected species requires an explicit --roles NAME=heavy|electron declaration.
Roles are independent of species names, ordering and charge; electron roles require charge -1.

Usage: gdtk_to_json.py --database DIR --species N2 O2 N O NO [--transport]
           [--electronic] [--ionisation-energy N=14.53413 O=13.61805]
           --roles N2=heavy O2=heavy N=heavy O=heavy NO=heavy --source TEXT -o out.json
"""

import argparse
import math
import re
from decimal import Decimal
from pathlib import Path

import species_json

# h c / k_B in cm K (exact in SI since 2019): converts a level energy in cm^-1
# to a temperature in K.
SECOND_RADIATION_CONSTANT = Decimal("6.62607015e-34") * Decimal("299792458") * 100 / Decimal("1.380649e-23")


# e N_A in J/(mol eV) (exact in SI since 2019).
MOLAR_ELECTRONVOLT = Decimal("1.602176634e-19") * Decimal("6.02214076e23")


def entry(name):
    """Regex for the start of db.<name> or db['<name>']."""
    return rf"db(?:\.{re.escape(name)}|\[\s*['\"]{re.escape(name)}['\"]\s*\])"


def file_name(name):
    return name.replace("+", "_plus").replace("e-", "e_minus") + ".lua"


def table(text, name, key):
    """Body of the Lua table db.<name>.<key> = { ... } (balanced braces)."""
    match = re.search(rf"{entry(name)}\.{key}\s*=\s*\{{", text)
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
    """db.<name>.<key> CEA transport fit as a canary transport curve in SI."""
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


def electronic_levels(text, name):
    """db.<name>.electronic_levels as degeneracies and temperatures in K."""
    body = table(text, name, "electronic_levels")
    g = [int(species_json.number(x)) for x in field(field_body(body, "g"), "value").split(",") if x.strip()]
    te = [x.strip() for x in field(field_body(body, "Te"), "value").split(",") if x.strip()]
    levels = {"reference": "NIST, as in the GDTk species database"}
    if len(g) != len(te):
        # GDTk uses the first min(len(g), len(Te)) levels.
        count = min(len(g), len(te))
        levels["notes"] = (f"the source lists {len(g)} degeneracies and {len(te)} level energies; "
                           f"the first {count} of each are used, as GDTk does")
        g, te = g[:count], te[:count]
    levels["degeneracies"] = g
    levels["temperatures"] = [float(Decimal(x) * SECOND_RADIATION_CONSTANT) for x in te]
    return {k: levels[k] for k in ("degeneracies", "temperatures", "reference", "notes") if k in levels}


def field_body(body, key):
    """Body of the nested table key = { ... } (balanced braces)."""
    match = re.search(rf"\b{key}\s*=\s*\{{", body)
    if not match:
        raise SystemExit(f"no {key} in table")
    depth, start = 1, match.end()
    for i in range(start, len(body)):
        depth += {"{": 1, "}": -1}.get(body[i], 0)
        if depth == 0:
            return body[start:i]
    raise SystemExit(f"unterminated {key} table")


def convert(path, name, with_transport=False, with_electronic=False):
    text = path.read_text(encoding="utf-8")
    mass = table(text, name, "M")
    value = re.search(r"value\s*=\s*([0-9.eE+-]+)", mass).group(1)
    charge = re.search(rf"{entry(name)}\.charge\s*=\s*([0-9.eE+-]+)", text)
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
    if with_electronic and re.search(rf"{entry(name)}\.electronic_levels\s*=", text):
        record["electronic_levels"] = electronic_levels(text, name)
    mass_reference = quoted(mass, "reference")
    if mass_reference:
        record["notes"] = f"molar mass: {mass_reference}"
    return record


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--database", required=True, type=Path)
    parser.add_argument("--species", required=True, nargs="+")
    parser.add_argument("--transport", action="store_true", help="also convert the transport data")
    parser.add_argument("--electronic", action="store_true", help="also convert the electronic levels")
    parser.add_argument("--ionisation-energy", nargs="*", default=[], metavar="SPECIES=EV",
                        help="ionisation energies in eV to add")
    parser.add_argument("--source", required=True, help="where the database files came from (version, path)")
    parser.add_argument("-o", "--output", required=True)
    parser.add_argument("--roles", required=True, nargs="+", metavar="NAME=ROLE",
                        help="explicit heavy or electron role for every selected species")
    args = parser.parse_args()
    roles = species_json.parse_roles(args.roles, args.species)
    species = {name: convert(args.database / file_name(name), name, args.transport, args.electronic)
               for name in args.species}
    for item in args.ionisation_energy:
        name, value = item.split("=")
        if name not in species:
            raise SystemExit(f"ionisation energy for {name}, which is not converted")
        species[name]["ionisation_energy"] = float(Decimal(value) * MOLAR_ELECTRONVOLT)
    for name, record in species.items():
        record["role"] = roles[name]
    species_json.dump(species_json.species_file(species, reference=args.source), args.output)


if __name__ == "__main__":
    main()
