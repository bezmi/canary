"""Structural and conversion checks for collision-integral JSON tools."""

from _schema_validation import validate_structure
import species_json


def validate(data, source):
    """Validate file metadata and exact pair dimensions before conversion."""
    validate_structure(data, "collision_integrals.schema.json", source)
    if data["format"] != "canary-collision-integrals-1":
        raise SystemExit(f"{source}: format must be canary-collision-integrals-1")
    if data["units"] != {"temperature": "K", "area": "m^2"}:
        raise SystemExit(f"{source}: collision units must be temperature K and area m^2")
    for index, pair in enumerate(data["pairs"]):
        if len(pair["species"]) != 2:
            raise SystemExit(f"{source}: pair {index}: requires exactly two species names")


def load(path):
    """Decode without duplicate keys and validate before returning the file."""
    data = species_json.load(path)
    validate(data, path)
    return data
