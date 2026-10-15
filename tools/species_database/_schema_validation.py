"""Structural validation against the repository's generated JSON schemas."""

import json
from pathlib import Path

import jsonschema


def validate_structure(data, schema_name, source):
    schema_path = Path(__file__).resolve().parents[2] / "database" / schema_name
    schema = json.loads(schema_path.read_text(encoding="utf-8"))
    for error in jsonschema.Draft202012Validator(schema).iter_errors(data):
        location = "/".join(map(str, error.path))
        raise SystemExit(f"{source}: {location}: {error.message}")
