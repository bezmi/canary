"""Shared structural and rate-model validation for reaction JSON tools.

Checks cover file structure, model fields and SI input domains. Integer storage,
converted-parameter representability and species-dependent constraints remain
the responsibility of make_reaction_mechanism.
"""

import math

from _schema_validation import validate_structure

import species_json


def _validate_rate(rate, role, context):
    model = rate["model"]
    if role == "backward" and model == "equilibrium":
        for field in ("A", "n", "activation_temperature", "translational_exponent"):
            if rate.get(field) is not None:
                raise SystemExit(f"{context}: equilibrium backward rate forbids {field}")
    elif model == "arrhenius":
        for field in ("A", "n", "activation_temperature"):
            if rate.get(field) is None:
                raise SystemExit(f"{context}: {role} Arrhenius rate requires {field}")
        location = f"{context}: {role} rate"
        if not math.isfinite(rate["A"]) or rate["A"] <= 0:
            raise SystemExit(f"{location}: A must be finite and positive")
        if not math.isfinite(rate["n"]):
            raise SystemExit(f"{location}: n must be finite")
        theta = rate["activation_temperature"]
        if not math.isfinite(theta) or theta < 0:
            raise SystemExit(f"{location}: activation_temperature must be finite and non-negative")
        exponent = rate.get("translational_exponent")
        if exponent is not None and (not math.isfinite(exponent) or not 0 <= exponent <= 1):
            raise SystemExit(f"{location}: translational_exponent must be finite and between 0 and 1")
    else:
        supported = "arrhenius" if role == "forward" else "arrhenius or equilibrium"
        raise SystemExit(f"{context}: {role} rate model must be {supported}, got {model!r}")


def _validate_stoichiometry(coefficients, role, context):
    if not coefficients:
        raise SystemExit(f"{context}: {role} must contain at least one species")
    for name, value in coefficients.items():
        if not math.isfinite(value) or value <= 0 or value != int(value):
            raise SystemExit(f"{context}: {role} coefficient of {name!r} must be a positive integer")


def validate(data, source):
    """Check a decoded file, raising SystemExit with source/record context."""
    validate_structure(data, "reactions.schema.json", source)
    if data["format"] != "canary-reactions-1":
        raise SystemExit(f"{source}: format must be canary-reactions-1")
    if data["units"] != {"quantity": "mol", "length": "m", "time": "s", "temperature": "K"}:
        raise SystemExit(f"{source}: units must be quantity mol, length m, time s and temperature K")
    lower = float(data["temperature_limits"]["lower"])
    upper = float(data["temperature_limits"]["upper"])
    if not math.isfinite(lower) or not math.isfinite(upper) or not 0 < lower < upper:
        raise SystemExit(f"{source}: temperature limits must be finite and satisfy 0 < lower < upper")
    if not data["reactions"]:
        raise SystemExit(f"{source}: at least one reaction is required")
    for index, record in enumerate(data["reactions"]):
        context = f"{source}: reaction {index} ({record['equation']})"
        _validate_stoichiometry(record["reactants"], "reactants", context)
        _validate_stoichiometry(record["products"], "products", context)
        third_body = record.get("third_body")
        if third_body is not None:
            for name, efficiency in third_body["efficiencies"].items():
                if not math.isfinite(efficiency) or efficiency < 0:
                    raise SystemExit(f"{context}: third-body efficiency of {name!r} must be finite and non-negative")
        _validate_rate(record["forward"], "forward", context)
        _validate_rate(record["backward"], "backward", context)


def load(path):
    """Decode without duplicate keys and validate before returning the file."""
    data = species_json.load(path)
    validate(data, path)
    return data
