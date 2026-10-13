"""Converts a GDTk collision-integral data file (cis['A:B'] = {pi_Omega_11 =
{A=, B=, C=, D=}, pi_Omega_22 = {...}}) to a canary collision-integral file
(format canary-collision-integrals-1, SI units).

GDTk's fits give pi Omega in Angstrom^2 as exp(D) T^(A (ln T)^2 + B ln T + C),
that is ln(pi Omega) = A (ln T)^3 + B (ln T)^2 + C ln T + D; D is converted to
m^2 by adding ln(1e-20). With --species, only pairs of the listed species are
written, and every such pair must be in the file.

Usage: gdtk_collision_integrals_to_json.py data.lua [--species N2 O2 ...]
           --source TEXT -o out.json
"""

import argparse
import json
import math
import re

FORMAT = "canary-collision-integrals-1"
UNITS = {"temperature": "K", "area": "m^2"}


def fit(body, key):
    match = re.search(rf"\b{key}\s*=\s*\{{([^}}]*)\}}", body)
    if not match:
        raise SystemExit(f"no {key} in: {body.strip()[:80]}")
    values = dict(re.findall(r"\b([ABCD])\s*=\s*([0-9.eE+-]+)", match.group(1)))
    return {"A": float(values["A"]), "B": float(values["B"]), "C": float(values["C"]),
            "D": float(values["D"]) + math.log(1.0e-20)}


def pairs(text):
    out = []
    for match in re.finditer(r"cis\[\s*['\"]([^'\":]+):([^'\"]+)['\"]\s*\]\s*=\s*\{", text):
        depth, start = 1, match.end()
        for i in range(start, len(text)):
            depth += {"{": 1, "}": -1}.get(text[i], 0)
            if depth == 0:
                body = text[start:i]
                break
        out.append({"species": [match.group(1), match.group(2)],
                    "pi_omega_11": fit(body, "pi_Omega_11"), "pi_omega_22": fit(body, "pi_Omega_22")})
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("data")
    parser.add_argument("--species", nargs="+")
    parser.add_argument("--source", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    text = "\n".join(line.split("--", 1)[0] for line in open(args.data, encoding="utf-8").read().splitlines())
    found = pairs(text)
    if args.species:
        wanted = set(args.species)
        found = [p for p in found if set(p["species"]) <= wanted]
        have = {frozenset(p["species"]) for p in found}
        for i, a in enumerate(args.species):
            for b in args.species[: i + 1]:
                if frozenset((a, b)) not in have:
                    raise SystemExit(f"no collision integrals for {a}:{b}")
    out = {"format": FORMAT, "units": UNITS, "pairs": found, "reference": args.source}
    with open(args.output, "w", encoding="utf-8") as stream:
        json.dump(out, stream, indent=3)
        stream.write("\n")


if __name__ == "__main__":
    main()
