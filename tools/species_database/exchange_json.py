"""Shared structural, model-field and SI input-domain validation for energy-exchange tools.

These checks match the file converters and source parameter domains in the
C++ builders. Species-dependent, interval and converted-parameter checks
belong to model construction.
"""

import math

from _schema_validation import validate_structure
import species_json


def _expect_model(record, expected, context):
    if record["model"] not in expected:
        raise SystemExit(f"{context}: model must be {' or '.join(expected)}, got {record['model']!r}")


def _require_field(record, field, required, context):
    present = record.get(field) is not None
    if present != required:
        action = "requires" if required else "forbids"
        raise SystemExit(f"{context}: {record['model']} {action} {field}")


def _finite_field(record, field, context):
    value = float(record[field])
    if not math.isfinite(value):
        raise SystemExit(f"{context}: {field} must be finite")
    return value


def _positive_field(record, field, context):
    value = _finite_field(record, field, context)
    if value <= 0:
        raise SystemExit(f"{context}: {field} must be positive")
    return value


def validate(data, source):
    """Check a decoded file, raising SystemExit with source/record context."""
    validate_structure(data, "energy_exchange.schema.json", source)
    if data["format"] != "canary-energy-exchange-1":
        raise SystemExit(f"{source}: format must be canary-energy-exchange-1")
    if data["units"] != {"energy": "J/mol", "temperature": "K", "area": "m^2"}:
        raise SystemExit(f"{source}: units must be energy J/mol, temperature K and area m^2")
    for name in ("chemistry_vibration", "electron_impact_ionisation"):
        for index, record in enumerate(data[name]):
            if record["reaction"] < 0:
                raise SystemExit(f"{source}: {name} {index} ({record['species']}): reaction index must be non-negative")
    for index, record in enumerate(data["vibration_translation"]):
        context = f"{source}: vibration_translation {index} ({record['species']}-{record['partner']})"
        _expect_model(record, ("landau_teller",), context)
        time = record["relaxation_time"]
        _expect_model(time, ("millikan_white",), context)
        for field in ("a", "b"):
            _finite_field(time, field, f"{context}: relaxation_time")
        correction = time.get("high_temperature_correction")
        if correction is not None:
            _expect_model(correction, ("park_1985", "park_1993", "kim_2021"), context)
            _require_field(correction, "reference_temperature", correction["model"] == "park_1993", context)
            _require_field(correction, "exponent", correction["model"] == "kim_2021", context)
            _positive_field(correction, "cross_section", context)
            if correction["model"] == "park_1993":
                _positive_field(correction, "reference_temperature", context)
            elif correction["model"] == "kim_2021":
                _finite_field(correction, "exponent", context)
    for index, record in enumerate(data["chemistry_vibration"]):
        context = f"{source}: chemistry_vibration {index} ({record['species']})"
        _expect_model(record, ("impartial_dissociation",), context)
        for field in ("dissociation_energy", "vibrational_temperature"):
            _positive_field(record, field, context)
    for index, record in enumerate(data["electron_translation"]):
        context = f"{source}: electron_translation {index} ({record['partner']})"
        _expect_model(record, ("elastic",), context)
        section = record["cross_section"]
        _expect_model(section, ("gnoffo_neutral", "coulomb"), context)
        if section["model"] == "gnoffo_neutral":
            for field in ("a", "b", "c", "maximum_temperature"):
                _require_field(section, field, True, context)
            for field in ("a", "b", "c"):
                _finite_field(section, field, f"{context}: cross_section")
            _positive_field(section, "maximum_temperature", context)
    for index, record in enumerate(data["vibration_electron"]):
        context = f"{source}: vibration_electron {index} ({record['species']})"
        _expect_model(record, ("landau_teller",), context)
        time = record["relaxation_time"]
        _expect_model(time, ("power_exponential",), context)
        for field in ("constant", "exponent"):
            _finite_field(time, field, f"{context}: relaxation_time")
        _positive_field(time, "reference_temperature", f"{context}: relaxation_time")
        coefficients = time["inverse_coefficients"]
        if len(coefficients) != 3 or not all(math.isfinite(value) for value in coefficients):
            raise SystemExit(f"{context}: relaxation_time inverse_coefficients must be 3 finite numbers")


def load(path):
    """Decode without duplicate keys and validate before returning the file."""
    data = species_json.load(path)
    validate(data, path)
    return data
