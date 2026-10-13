"""Converts a GDTk energy-exchange file as written by prep-kinetics (a Lua
table mechanism[k] = {...}) to a canary energy-exchange file (format
canary-energy-exchange-1, SI units).

Converted entries:
- type 'V-T', rate 'Landau-Teller', with relaxation_time
  {model='Millikan-White', a=, b=} or a high-temperature correction of it,
  {model='ParkHTC' | 'ParkHTC2' | 'KimHTC', sigma=, [exponent=,]
  submodel={model='Millikan-White', a=, b=}}: a vibration_translation entry
  of species p with collision partner q. ParkHTC is Park's correction with a
  constant cross-section sigma (m^2), ParkHTC2 Park's 1993 form with sigma at
  the reference temperature 50000 K, KimHTC Kim and Jo's (2021) sigma
  T^exponent.
- type 'C-V', rate 'Marrone-Treanor', with coupling_model
  {model='ImpartialDissociation', D=, Thetav=}: a chemistry_vibration entry
  of species p in reaction reaction_index (0-based, in mechanism order),
  with dissociation energy D (J/mol) and vibrational temperature Thetav (K).
- type 'E-T', rate 'ElectronExchange', with exchange_cross_section
  {type='GnoffoNeutral', a=, b=, c=} (sigma = a + b T_e + c T_e^2, T_e limited
  to 30000 K) or {type='Coulomb'}: an electron_translation entry of partner
  q.
- type 'E-V', rate 'Landau-Teller', with relaxation_time {model='Candler',
  a_low=, ..., c_high=}: a vibration_electron entry of species p, with the
  blending range 6950-7050 K.
C-V entries with a non-positive dissociation energy (GDTk writes them for
ions) are dropped, and so, with --reactions, are C-V entries of a species
that has no net change in its reaction; both have no effect. Any other
entry is an error. With --reactions and --ionisation-label, the
reactions of that canary reactions file with that label become
electron_impact_ionisation entries (the neutral reactant that becomes an
ion); GDTk's files do not contain them.

Usage: gdtk_exchange_to_json.py exchange.lua [--reactions reactions.json
           --ionisation-label EII] --source TEXT -o out.json
"""

import argparse
import json
import re

FORMAT = "canary-energy-exchange-1"
UNITS = {"energy": "J/mol", "temperature": "K", "area": "m^2"}
# Reference temperature of the cross-section of GDTk's ParkHTC2 correction.
PARK_REFERENCE_TEMPERATURE = 50000.0
# Upper temperature limit of GDTk's GnoffoNeutral cross-section fits.
GNOFFO_MAXIMUM_TEMPERATURE = 30000.0
# Blending range of GDTk's Candler relaxation-time fits, K.
CANDLER_BLEND = (6950.0, 7050.0)


def blocks(text):
    """Bodies of the mechanism[k] = {...} entries, in file order."""
    out = []
    for match in re.finditer(r"mechanism\[(\d+)\]\s*=\s*\{", text):
        depth, start = 1, match.end()
        for i in range(start, len(text)):
            depth += {"{": 1, "}": -1}.get(text[i], 0)
            if depth == 0:
                out.append(text[start:i])
                break
    return out


def string_field(body, key):
    match = re.search(rf"\b{key}\s*=\s*'([^']*)'", body)
    if not match:
        raise SystemExit(f"missing {key} in: {body.strip()[:80]}")
    return match.group(1)


def number_field(body, key):
    match = re.search(rf"\b{key}\s*=\s*([0-9.eE+-]+)", body)
    if not match:
        raise SystemExit(f"missing {key} in: {body.strip()[:80]}")
    return float(match.group(1))


def millikan_white(body):
    if string_field(body, "model") != "Millikan-White":
        raise SystemExit(f"unsupported relaxation time: {body.strip()[:80]}")
    return {"model": "millikan_white", "a": number_field(body, "a"), "b": number_field(body, "b")}


def vibration_translation(body):
    if string_field(body, "rate") != "Landau-Teller":
        raise SystemExit(f"unsupported V-T rate: {body.strip()[:80]}")
    relaxation = re.search(r"relaxation_time\s*=\s*(\{.*\})", body, re.S).group(1)
    model = re.match(r"\{\s*model\s*=\s*'([^']*)'", relaxation).group(1)
    if model in ("ParkHTC", "ParkHTC2", "KimHTC"):
        sub = re.search(r"submodel\s*=\s*(\{[^}]*\})", relaxation).group(1)
        time = millikan_white(sub)
        outer = re.sub(r"submodel\s*=\s*\{[^}]*\}", "", relaxation)
        correction = {"cross_section": number_field(outer, "sigma")}
        if model == "ParkHTC":
            correction = {"model": "park_1985", **correction}
        elif model == "ParkHTC2":
            correction = {"model": "park_1993", **correction, "reference_temperature": PARK_REFERENCE_TEMPERATURE}
        else:
            correction = {"model": "kim_2021", **correction, "exponent": number_field(outer, "exponent")}
        time["high_temperature_correction"] = correction
    else:
        time = millikan_white(relaxation)
    return {"species": string_field(body, "p"), "partner": string_field(body, "q"),
            "model": "landau_teller", "relaxation_time": time}


def electron_translation(body):
    if string_field(body, "rate") != "ElectronExchange" or string_field(body, "p") != "e-":
        raise SystemExit(f"unsupported E-T entry: {body.strip()[:80]}")
    section = re.search(r"exchange_cross_section\s*=\s*(\{[^}]*\})", body).group(1)
    kind = string_field(section, "type")
    if kind == "GnoffoNeutral":
        cross_section = {"model": "gnoffo_neutral", "a": number_field(section, "a"), "b": number_field(section, "b"),
                         "c": number_field(section, "c"), "maximum_temperature": GNOFFO_MAXIMUM_TEMPERATURE}
    elif kind == "Coulomb":
        cross_section = {"model": "coulomb"}
    else:
        raise SystemExit(f"unsupported exchange cross-section {kind}")
    return {"partner": string_field(body, "q"), "model": "elastic", "cross_section": cross_section}


def vibration_electron(body):
    relaxation = re.search(r"relaxation_time\s*=\s*(\{[^}]*\})", body).group(1)
    if string_field(body, "rate") != "Landau-Teller" or string_field(relaxation, "model") != "Candler":
        raise SystemExit(f"unsupported E-V entry: {body.strip()[:80]}")
    time = {"model": "candler"}
    for key in ("a_low", "b_low", "c_low", "a_high", "b_high", "c_high"):
        time[key] = number_field(relaxation, key)
    time["blend_lower"], time["blend_upper"] = CANDLER_BLEND
    return {"species": string_field(body, "p"), "model": "landau_teller", "relaxation_time": time}


def ionisations(reactions_file, label):
    """electron_impact_ionisation entries for the reactions with label."""
    reactions = json.loads(open(reactions_file, encoding="utf-8").read())["reactions"]
    out = []
    for index, reaction in enumerate(reactions):
        if reaction.get("label") != label:
            continue
        ions = [s for s in reaction["products"] if s.endswith("+")]
        neutral = [s for s in reaction["reactants"] if s + "+" in ions]
        if len(neutral) != 1:
            raise SystemExit(f"{reaction['equation']}: no single neutral that becomes an ion")
        out.append({"reaction": index, "species": neutral[0]})
    if not out:
        raise SystemExit(f"no reactions with label {label}")
    return out


def chemistry_vibration(body):
    if string_field(body, "rate") != "Marrone-Treanor":
        raise SystemExit(f"unsupported C-V rate: {body.strip()[:80]}")
    coupling = re.search(r"coupling_model\s*=\s*(\{[^}]*\})", body).group(1)
    if string_field(coupling, "model") != "ImpartialDissociation":
        raise SystemExit(f"unsupported coupling model: {coupling}")
    return {"reaction": int(number_field(body, "reaction_index")), "species": string_field(body, "p"),
            "model": "impartial_dissociation", "dissociation_energy": number_field(coupling, "D"),
            "vibrational_temperature": number_field(coupling, "Thetav")}


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("exchange")
    parser.add_argument("--reactions", help="canary reactions file for electron_impact_ionisation entries")
    parser.add_argument("--ionisation-label", default="EII")
    parser.add_argument("--source", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    text = "\n".join(line.split("--", 1)[0] for line in open(args.exchange, encoding="utf-8").read().splitlines())
    out = {"format": FORMAT, "units": UNITS, "vibration_translation": [], "chemistry_vibration": [],
           "electron_translation": [], "vibration_electron": [], "electron_impact_ionisation": []}
    for body in blocks(text):
        kind = string_field(body, "type")
        if kind == "V-T":
            out["vibration_translation"].append(vibration_translation(body))
        elif kind == "C-V":
            entry = chemistry_vibration(body)
            if entry["dissociation_energy"] > 0.0:
                out["chemistry_vibration"].append(entry)
        elif kind == "E-T":
            out["electron_translation"].append(electron_translation(body))
        elif kind == "E-V":
            out["vibration_electron"].append(vibration_electron(body))
        else:
            raise SystemExit(f"unsupported mechanism type {kind}")
    if args.reactions:
        reactions = json.loads(open(args.reactions, encoding="utf-8").read())["reactions"]

        def changes(entry):
            reaction = reactions[entry["reaction"]]
            species = entry["species"]
            return reaction["products"].get(species, 0.0) != reaction["reactants"].get(species, 0.0)

        out["chemistry_vibration"] = [e for e in out["chemistry_vibration"] if changes(e)]
        out["electron_impact_ionisation"] = ionisations(args.reactions, args.ionisation_label)
    out["reference"] = args.source
    with open(args.output, "w", encoding="utf-8") as stream:
        json.dump(out, stream, indent=3)
        stream.write("\n")


if __name__ == "__main__":
    main()
