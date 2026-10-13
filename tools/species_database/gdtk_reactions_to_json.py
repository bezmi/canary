"""Converts a GDTk reaction-mechanism input file (Reaction{...} entries,
pre-exponential factors in CGS units) to a canary reactions file
(format canary-reactions-1, SI units).

Each Reaction{...} entry gives an equation such as 'N2 + N2 <=> N + N + N2'
or 'NO + M <=> N + O + M', a forward rate fr={'Arrhenius', A=, n=, C=}
(k = A T^n exp(-C / T), A in (cm^3/mol)^(order - 1) / s, C in K) or
fr={'Park', A=, n=, C=, s=} (the same at the rate-controlling temperature
T^s T_ve^(1 - s); written as an Arrhenius rate with translational_exponent
s; mode="X" names the species whose vibrational temperature is T_ve, and
rateControllingTemperature="electron-electronic" the electron temperature,
written as other_temperature "X" or "electron"), optionally
a backward rate br={...} of the same form (otherwise the backward rate comes
from the equilibrium constant; br={'fromEqConst',
rateControllingTemperature="electron-electronic"} evaluates it at the
electron temperature, written as other_temperature "electron"), optionally collider efficiencies
efficiencies={['N2']=1.0, ...} for reactions with M, and a label. The
reaction order counts every reactant, including M. Temperature limits come
from Config{tempLimits={lower=, upper=}}.

Usage: gdtk_reactions_to_json.py mechanism.inp --source TEXT -o out.json
"""

import argparse
import json
import re
from decimal import Decimal

FORMAT = "canary-reactions-1"
UNITS = {"quantity": "mol", "length": "m", "time": "s", "temperature": "K"}


def entries(text, keyword):
    """Bodies of keyword{...} blocks with balanced braces."""
    out, pos = [], 0
    while True:
        match = re.compile(rf"\b{keyword}\s*\{{").search(text, pos)
        if not match:
            return out
        depth, start = 1, match.end()
        for i in range(start, len(text)):
            depth += {"{": 1, "}": -1}.get(text[i], 0)
            if depth == 0:
                out.append(text[start:i])
                pos = i + 1
                break


def strip_comments(text):
    return "\n".join(line.split("--", 1)[0] for line in text.splitlines())


def side(expression):
    """'N + N + N2' -> ({'N': 2, 'N2': 1}, has_M)."""
    counts, third_body = {}, False
    for term in re.split(r"\s+\+\s+", expression.strip()):
        term = term.strip()
        match = re.fullmatch(r"(\d*)\s*([A-Za-z][A-Za-z0-9_+-]*?)", term)
        if not match:
            raise SystemExit(f"cannot read term '{term}'")
        coefficient, name = int(match.group(1) or 1), match.group(2)
        if name == "M":
            third_body = True
        else:
            counts[name] = counts.get(name, 0) + coefficient
    return counts, third_body


# Energy-mode names of GDTk mechanism files and the canary other_temperature
# they mean.
OTHER_TEMPERATURES = {"electron-electronic": "electron"}


def other_temperature(body):
    """other_temperature from mode="X" (a species) or
    rateControllingTemperature="..." in a rate table, or None."""
    mode = re.search(r"\bmode\s*=\s*['\"]([^'\"]+)['\"]", body)
    rct = re.search(r"\brateControllingTemperature\s*=\s*['\"]([^'\"]+)['\"]", body)
    if mode and rct:
        raise SystemExit(f"both mode and rateControllingTemperature in: {body.strip()[:80]}")
    if rct:
        if rct.group(1) not in OTHER_TEMPERATURES:
            raise SystemExit(f"unknown rateControllingTemperature {rct.group(1)}")
        return OTHER_TEMPERATURES[rct.group(1)]
    if mode:
        return OTHER_TEMPERATURES.get(mode.group(1), mode.group(1))
    return None


def arrhenius(body, order):
    match = re.search(r"\{\s*['\"](Arrhenius|Park)['\"]\s*,([^}]*)\}", body)
    if not match:
        raise SystemExit(f"unsupported rate form in: {body.strip()[:60]}")
    fields = dict(re.findall(r"(\w+)\s*=\s*([0-9.eE+-]+)", match.group(2)))
    a_cgs = Decimal(fields["A"])
    # (cm^3/mol)^(order - 1) / s -> (m^3/mol)^(order - 1) / s
    a_si = a_cgs * (Decimal("1e-6") ** (order - 1))
    rate = {
        "model": "arrhenius",
        "A": float(a_si),
        "n": float(fields.get("n", "0")),
        "activation_temperature": float(fields.get("C", "0")),
    }
    other = other_temperature(match.group(2))
    if match.group(1) == "Park":
        rate["translational_exponent"] = float(fields["s"])
    elif other is not None:
        # An Arrhenius rate at the named temperature alone.
        rate["translational_exponent"] = 0.0
    if other is not None:
        rate["other_temperature"] = other
    return rate


def convert_reaction(body):
    equation = re.search(r"['\"]([^'\"]*<=>[^'\"]*|[^'\"]*=>[^'\"]*)['\"]", body).group(1)
    if "<=>" not in equation:
        raise SystemExit(f"only reversible reactions (<=>) are converted: {equation}")
    left, right = (s.strip() for s in equation.split("<=>"))
    reactants, m_left = side(left)
    products, m_right = side(right)
    if m_left != m_right:
        raise SystemExit(f"M on one side only: {equation}")
    order = sum(reactants.values()) + (1 if m_left else 0)
    forward_body = re.search(r"\bfr\s*=\s*(\{[^}]*\})", body).group(1)
    reaction = {
        "equation": equation,
        "reactants": {k: float(v) for k, v in reactants.items()},
        "products": {k: float(v) for k, v in products.items()},
        "forward": arrhenius(forward_body, order),
    }
    backward = re.search(r"\bbr\s*=\s*(\{[^}]*\})", body)
    if backward and re.search(r"['\"]fromEqConst['\"]", backward.group(1)):
        reaction["backward"] = {"model": "equilibrium"}
        other = other_temperature(backward.group(1))
        if other is not None:
            reaction["backward"]["other_temperature"] = other
    elif backward:
        backward_order = sum(products.values()) + (1 if m_right else 0)
        reaction["backward"] = arrhenius(backward.group(1), backward_order)
    else:
        reaction["backward"] = {"model": "equilibrium"}
    if m_left:
        efficiencies = re.search(r"\befficiencies\s*=\s*\{([^}]*)\}", body)
        pairs = re.findall(r"\[\s*['\"]([^'\"]+)['\"]\s*\]\s*=\s*([0-9.eE+-]+)", efficiencies.group(1)) if efficiencies else []
        reaction["third_body"] = {"efficiencies": {k: float(v) for k, v in pairs}}
    label = re.search(r"\blabel\s*=\s*['\"]([^'\"]*)['\"]", body)
    if label:
        reaction["label"] = label.group(1)
    return reaction


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("mechanism")
    parser.add_argument("--source", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()
    text = strip_comments(open(args.mechanism, encoding="utf-8").read())
    limits = {"lower": 300.0, "upper": 50000.0}
    for config in entries(text, "Config"):
        match = re.search(r"tempLimits\s*=\s*\{([^}]*)\}", config)
        if match:
            limits.update({k: float(v) for k, v in re.findall(r"(lower|upper)\s*=\s*([0-9.eE+-]+)", match.group(1))})
    out = {
        "format": FORMAT,
        "units": UNITS,
        "temperature_limits": limits,
        "reactions": [convert_reaction(body) for body in entries(text, "Reaction")],
        "reference": args.source,
    }
    with open(args.output, "w", encoding="utf-8") as stream:
        json.dump(out, stream, indent=3)
        stream.write("\n")


if __name__ == "__main__":
    main()
