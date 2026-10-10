"""Validates species files against the JSON Schema of gaslib species files.

Usage: check_species_files.py <schema.json> <species file>...
"""

import json
import sys

import jsonschema

schema = json.loads(open(sys.argv[1], encoding="utf-8").read())
jsonschema.Draft202012Validator.check_schema(schema)
failures = 0
for path in sys.argv[2:]:
    data = json.loads(open(path, encoding="utf-8").read())
    errors = list(jsonschema.Draft202012Validator(schema).iter_errors(data))
    for error in errors:
        print(f"{path}: {'/'.join(map(str, error.path))}: {error.message}")
    failures += len(errors)
    print(f"{path}: {'valid' if not errors else 'INVALID'}")
sys.exit(1 if failures else 0)
