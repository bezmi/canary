"""Shared structural and rate-model validation for reaction JSON tools.

These checks match the database reader and mechanism_input conversion.
Numerical domains, converted parameters and species-dependent constraints
remain the responsibility of make_reaction_mechanism.
"""

import json
from pathlib import Path

import jsonschema

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
    else:
        supported = "arrhenius" if role == "forward" else "arrhenius or equilibrium"
        raise SystemExit(f"{context}: {role} rate model must be {supported}, got {model!r}")


def validate(data, source):
    """Check a decoded file, raising SystemExit with source/record context."""
    schema_path = Path(__file__).resolve().parents[2] / "database" / "reactions.schema.json"
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    for error in jsonschema.Draft202012Validator(schema).iter_errors(data):
        location = "/".join(map(str, error.path))
        raise SystemExit(f"{source}: {location}: {error.message}")
    if data["format"] != "canary-reactions-1":
        raise SystemExit(f"{source}: format must be canary-reactions-1")
    if data["units"] != {"quantity": "mol", "length": "m", "time": "s", "temperature": "K"}:
        raise SystemExit(f"{source}: units must be quantity mol, length m, time s and temperature K")
    for index, record in enumerate(data["reactions"]):
        context = f"{source}: reaction {index} ({record['equation']})"
        _validate_rate(record["forward"], "forward", context)
        _validate_rate(record["backward"], "backward", context)


def load(path):
    """Decode without duplicate keys and validate before returning the file."""
    data = species_json.load(path)
    validate(data, path)
    return data
