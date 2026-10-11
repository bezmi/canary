"""Converts a GDTk energy-exchange file as written by prep-kinetics (a Lua
table mechanism[k] = {...}) to a gaslib energy-exchange file (format
gaslib-energy-exchange-1, SI units).

Converted entries:
- type 'V-T', rate 'Landau-Teller', with relaxation_time
  {model='Millikan-White', a=, b=} or {model='ParkHTC2', sigma=,
  submodel={model='Millikan-White', a=, b=}}: a vibration_translation entry
  of species p with collision partner q. ParkHTC2 is Park's 1993
  high-temperature correction with cross-section sigma (m^2) at the
  reference temperature 50000 K.
- type 'C-V', rate 'Marrone-Treanor', with coupling_model
  {model='ImpartialDissociation', D=, Thetav=}: a chemistry_vibration entry
  of species p in reaction reaction_index (0-based, in mechanism order),
  with dissociation energy D (J/mol) and vibrational temperature Thetav (K).
Any other entry is an error.

Usage: gdtk_exchange_to_json.py exchange.lua --source TEXT -o out.json
"""

import argparse
import json
import re

FORMAT = "gaslib-energy-exchange-1"
UNITS = {"energy": "J/mol", "temperature": "K", "area": "m^2"}
# Reference temperature of the cross-section of GDTk's ParkHTC2 correction.
PARK_REFERENCE_TEMPERATURE = 50000.0


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
    if model == "ParkHTC2":
        sub = re.search(r"submodel\s*=\s*(\{[^}]*\})", relaxation).group(1)
        time = millikan_white(sub)
        time["high_temperature_correction"] = {
            "model": "park_1993",
            "cross_section": number_field(relaxation, "sigma"),
            "reference_temperature": PARK_REFERENCE_TEMPERATURE,
        }
    else:
        time = millikan_white(relaxation)
    return {"species": string_field(body, "p"), "partner": string_field(body, "q"),
            "model": "landau_teller", "relaxation_time": time}


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
    parser.add_argument("--source", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    text = "\n".join(line.split("--", 1)[0] for line in open(args.exchange, encoding="utf-8").read().splitlines())
    out = {"format": FORMAT, "units": UNITS, "vibration_translation": [], "chemistry_vibration": []}
    for body in blocks(text):
        kind = string_field(body, "type")
        if kind == "V-T":
            out["vibration_translation"].append(vibration_translation(body))
        elif kind == "C-V":
            out["chemistry_vibration"].append(chemistry_vibration(body))
        else:
            raise SystemExit(f"unsupported mechanism type {kind}")
    out["reference"] = args.source
    with open(args.output, "w", encoding="utf-8") as stream:
        json.dump(out, stream, indent=3)
        stream.write("\n")


if __name__ == "__main__":
    main()
