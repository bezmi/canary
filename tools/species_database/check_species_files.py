# SPDX-License-Identifier: AGPL-3.0-or-later WITH AdditionRef-Canary-Cerberus-Exception
# Copyright (C) 2026 Shahzeb Imran and the canary contributors
"""Validates species, reactions, energy-exchange and collision-integral files
against the structural JSON Schemas of canary database files.

Success does not establish physical validity or that a model can be constructed.
Apply C++ input conversion and model construction checks after reading a file; see
database/README.md for the validation responsibilities.

Usage: check_species_files.py <schema.json> <file>... [-- <schema.json> <file>...]...

Each group of arguments, separated by --, is a schema followed by the files
that must conform to it.
"""

import json
import sys

import jsonschema

import species_json


def groups(arguments):
    group = []
    for argument in arguments:
        if argument == "--":
            yield group
            group = []
        else:
            group.append(argument)
    yield group


argument_groups = list(groups(sys.argv[1:]))
if any(len(group) < 2 for group in argument_groups):
    raise SystemExit(
        "usage: check_species_files.py <schema.json> <file>... "
        "[-- <schema.json> <file>...]...\n"
        "each group requires a schema and at least one input file"
    )

failures = 0
for schema_path, *paths in argument_groups:
    schema = json.loads(open(schema_path, encoding="utf-8").read())
    jsonschema.Draft202012Validator.check_schema(schema)
    for path in paths:
        data = species_json.load(path)
        errors = list(jsonschema.Draft202012Validator(schema).iter_errors(data))
        for error in errors:
            print(f"{path}: {'/'.join(map(str, error.path))}: {error.message}")
        failures += len(errors)
        print(f"{path}: {'valid' if not errors else 'INVALID'} ({schema_path})")
sys.exit(1 if failures else 0)
